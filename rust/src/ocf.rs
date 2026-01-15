//! Orientation-aware cfDNA Fragmentation (OCF) calculation
//!
//! Measures fragment end orientation patterns around open chromatin regions (OCRs).
//! OCF captures tissue-of-origin signals based on strand asymmetry at regulatory elements.
//!
//! ## Key Features
//! - **Strand asymmetry**: OCF = (Upstream - Downstream) / (Upstream + Downstream)
//! - **Region types**: Calculates per-OCR and synchronized (pooled) scores
//! - **GC correction**: Per-fragment weight adjustment via correction factors
//! - **Panel mode**: On-target and off-target split outputs
//!
//! ## Output
//! - {sample}.OCF.tsv: Per-region scores and fragment counts
//! - {sample}.OCF.sync.tsv: Synchronized regions for cross-sample comparison

use pyo3::prelude::*;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, Write};
use std::collections::HashMap;


/// Calculate Orientation-aware cfDNA Fragmentation (OCF)
/// 
/// # Arguments
/// * `bed_path` - Path to the input .bed.gz file
/// * `ocr_path` - Path to the Open Chromatin Regions (OCR) BED file
/// * `output_dir` - Directory to save output files
#[pyfunction]
pub fn calculate_ocf(
    _py: Python,
    bed_path: PathBuf,
    ocr_path: PathBuf,
    output_dir: PathBuf,
) -> PyResult<()> {
    let mut chrom_map = ChromosomeMap::new();
    let consumer = OcfConsumer::new(&ocr_path, &mut chrom_map, None, None)  // No target_regions for legacy API
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        
    let analyzer = FragmentAnalyzer::new(consumer, 100_000);
    let final_consumer = analyzer.process_file(&bed_path, &mut chrom_map, false)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        
    final_consumer.write_output(&output_dir)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    Ok(())
}

use crate::bed::{Fragment, ChromosomeMap};
use crate::engine::{FragmentConsumer, FragmentAnalyzer};
use std::sync::Arc;
use std::path::Path;
use anyhow::{Result, Context};
use coitrees::{COITree, IntervalNode, IntervalTree};
use crate::gc_correction::CorrectionFactors;

#[derive(Clone, Default)]
struct LabelStats {
    left_pos: Vec<f64>,  // Size 2000
    right_pos: Vec<f64>, // Size 2000
    total_starts: f64,
    total_ends: f64,
}

impl LabelStats {
    fn new() -> Self {
        Self {
            left_pos: vec![0.0; 2000],
            right_pos: vec![0.0; 2000],
            total_starts: 0.0,
            total_ends: 0.0,
        }
    }
    
    fn merge(&mut self, other: &Self) {
        for (a, b) in self.left_pos.iter_mut().zip(other.left_pos.iter()) { *a += *b; }
        for (a, b) in self.right_pos.iter_mut().zip(other.right_pos.iter()) { *a += *b; }
        self.total_starts += other.total_starts;
        self.total_ends += other.total_ends;
    }
}

// Region with pre-mapped label ID
struct OcrRegionInfo {
    // chrom stored in tree map key
    start: u64,
    end: u64,
    label_id: usize,
}

/// Parse a single OCR line
fn parse_ocr_line(
    line: &str,
    chrom_map: &mut ChromosomeMap,
    label_to_id: &mut HashMap<String, usize>,
    labels: &mut Vec<String>,
    nodes_by_chrom: &mut HashMap<u32, Vec<IntervalNode<usize, u32>>>,
    region_infos: &mut Vec<OcrRegionInfo>,
) {
    let line = line.trim();
    if line.is_empty() || line.starts_with('#') {
        return;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 4 {
        return;
    }
    
    let chrom = fields[0];
    let start: u64 = fields[1].parse().unwrap_or(0);
    let end: u64 = fields[2].parse().unwrap_or(0);
    let label = fields[3].to_string();
    
    // Map Label
    let label_id = if let Some(&id) = label_to_id.get(&label) {
        id
    } else {
        let id = labels.len();
        label_to_id.insert(label.clone(), id);
        labels.push(label);
        id
    };
    
    // Map Chrom
    let chrom_norm = chrom.trim_start_matches("chr");
    let chrom_id = chrom_map.get_id(chrom_norm);
    
    // Store Region Info
    let info_idx = region_infos.len();
    region_infos.push(OcrRegionInfo {
        start,
        end,
        label_id,
    });
    
    // Add to Tree nodes
    // Standard overlap query.
    // COITree (closed): [start, end-1]
    let s = start as u32;
    let e = end as u32;
    let e_closed = if e > s { e - 1 } else { s };
    
    nodes_by_chrom.entry(chrom_id).or_default().push(
        IntervalNode::new(s as i32, e_closed as i32, info_idx)
    );
}

#[derive(Clone)]
pub struct OcfConsumer {
    // Shared state
    trees: Arc<HashMap<u32, COITree<usize, u32>>>,
    region_infos: Arc<Vec<OcrRegionInfo>>, // Index -> Info
    labels: Arc<Vec<String>>, // LabelID -> Name
    factors: Option<Arc<CorrectionFactors>>,
    
    // Target regions for on/off-target split (panel data)
    target_tree: Option<Arc<HashMap<u32, COITree<(), u32>>>>,
    
    // Thread-local state
    // Indexed by LabelID
    // Off-target stats (primary - unbiased)
    stats: Vec<LabelStats>,
    // On-target stats (for comparison - PCR biased)
    stats_on: Vec<LabelStats>,
}

impl OcfConsumer {
    pub fn new(
        ocr_path: &Path, 
        chrom_map: &mut ChromosomeMap, 
        factors: Option<Arc<CorrectionFactors>>,
        target_regions: Option<Arc<HashMap<u32, COITree<(), u32>>>>,
    ) -> Result<Self> {
        // 1. Parse OCR File
        // Format: chrom start end label
        let reader = crate::bed::get_reader(ocr_path)?;
        
        let mut label_to_id: HashMap<String, usize> = HashMap::new();
        let mut labels: Vec<String> = Vec::new();
        
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        let mut region_infos = Vec::new(); // Implicitly indexed by order of insertion
        
        for line in reader.lines() {
            let line = line?;
            parse_ocr_line(&line, chrom_map, &mut label_to_id, &mut labels, &mut nodes_by_chrom, &mut region_infos);
        }
        
        // Build Trees
        let mut trees = HashMap::new();
        for (cid, nodes) in nodes_by_chrom {
            trees.insert(cid, COITree::new(&nodes));
        }
        
        // Init stats (off-target and on-target)
        let stats = vec![LabelStats::new(); labels.len()];
        let stats_on = vec![LabelStats::new(); labels.len()];
        
        if target_regions.is_some() {
            log::info!("OCF: Panel mode enabled, on/off-target split active");
        }
        
        Ok(Self {
            trees: Arc::new(trees),
            region_infos: Arc::new(region_infos),
            labels: Arc::new(labels),
            factors,
            target_tree: target_regions,
            stats,
            stats_on,
        })
    }
    
    /// Check if a fragment overlaps any target region
    fn is_on_target(&self, fragment: &Fragment) -> bool {
        if let Some(ref target_tree) = self.target_tree {
            if let Some(tree) = target_tree.get(&fragment.chrom_id) {
                let start = fragment.start as i32;
                let end = if fragment.end > fragment.start { (fragment.end - 1) as i32 } else { start };
                let mut found = false;
                tree.query(start, end, |_| { found = true; });
                return found;
            }
        }
        false
    }
    
    pub fn write_output(&self, output_dir: &Path) -> Result<()> {
        // Write off-target stats (primary - unbiased for biomarkers)
        self.write_stats(&self.stats, output_dir, "all.ocf.tsv", "all.sync.tsv")?;
        
        // Write on-target stats if we have any data
        let has_on_target_data = self.stats_on.iter().any(|s| s.total_starts > 0.0 || s.total_ends > 0.0);
        if has_on_target_data {
            self.write_stats(&self.stats_on, output_dir, "all.ocf.ontarget.tsv", "all.sync.ontarget.tsv")?;
            log::info!("OCF: Wrote on-target output files");
        }
        
        Ok(())
    }
    
    fn write_stats(&self, stats: &[LabelStats], output_dir: &Path, summary_name: &str, sync_name: &str) -> Result<()> {
        // Output 1: OCF summary
        let summary_path = output_dir.join(summary_name);
        let mut summary_file = File::create(&summary_path)
            .with_context(|| format!("Failed to create summary file: {:?}", summary_path))?;
        writeln!(summary_file, "tissue\tOCF")?;

        // Output 2: sync positions
        let sync_path = output_dir.join(sync_name);
        let mut sync_file = File::create(&sync_path)
             .with_context(|| format!("Failed to create sync file: {:?}", sync_path))?;
        writeln!(sync_file, "tissue\tposition\tleft_count\tleft_norm\tright_count\tright_norm")?;

        // Sort labels
        let mut label_indices: Vec<usize> = (0..self.labels.len()).collect();
        label_indices.sort_by_key(|&i| &self.labels[i]);
        
        for &id in &label_indices {
            let label = &self.labels[id];
            let s = &stats[id];
            
            let ts = if s.total_starts > 0.0 { s.total_starts / 10000.0 } else { 1.0 };
            let te = if s.total_ends > 0.0 { s.total_ends / 10000.0 } else { 1.0 };
            
            let peak = 60;
            let bin_width = 10;
            let mut trueends = 0.0;
            let mut background = 0.0;
            
            for k in 0..2000 {
                let l_count = s.left_pos[k];
                let r_count = s.right_pos[k];
                let l_norm = l_count / ts;
                let r_norm = r_count / te;
                let loc = k as i64 - 1000;
                
                writeln!(sync_file, "{}\t{}\t{:.2}\t{:.6}\t{:.2}\t{:.6}", 
                    label, loc, l_count, l_norm, r_count, r_norm)?;
                
                if loc >= -peak - bin_width && loc <= -peak + bin_width {
                    trueends += r_norm;
                    background += l_norm;
                } else if loc >= peak - bin_width && loc <= peak + bin_width {
                    trueends += l_norm;
                    background += r_norm;
                }
            }
            
            let ocf_score = trueends - background;
            writeln!(summary_file, "{}\t{:.6}", label, ocf_score)?;
        }
        
        Ok(())
    }
}

impl FragmentConsumer for OcfConsumer {
    fn name(&self) -> &str {
        "OCF"
    }

    fn consume(&mut self, fragment: &Fragment) {
        // Check if fragment is on-target (for routing to correct stats vector)
        let on_target = self.is_on_target(fragment);
        
        if let Some(tree) = self.trees.get(&fragment.chrom_id) {
            let start = fragment.start as i32;
            let end_closed = if fragment.end > fragment.start { (fragment.end - 1) as i32 } else { start };
            
            tree.query(start, end_closed, |node| {
                let region_idx = node.metadata.to_owned();
                let region = &self.region_infos[region_idx];
                
                // Route to correct stats vector based on on-target status
                let label_stats = if on_target {
                    &mut self.stats_on[region.label_id]
                } else {
                    &mut self.stats[region.label_id]
                };
                
                // Logic per Legacy:
                let r_start = fragment.start; // u64
                let r_end = fragment.end;     // u64
                let reg_start = region.start; // u64
                let reg_end = region.end;     // u64
                
                let gc_pct = (fragment.gc * 100.0).round() as u8;
                let weight = if let Some(ref factors) = self.factors {
                        factors.get_factor(fragment.length, gc_pct)
                } else { 1.0 };
                
                // Starts
                if r_start >= reg_start {
                    let s = (r_start - reg_start) as usize;
                    if s < 2000 {
                        label_stats.left_pos[s] += weight;
                        label_stats.total_starts += weight;
                    } else {
                        label_stats.total_starts += weight;
                    }
                }
                
                // Ends
                if r_end <= reg_end {
                    let diff = r_end.wrapping_sub(reg_start).wrapping_add(1);
                    let e = diff as usize;
                    
                    if e < 2000 {
                         label_stats.right_pos[e] += weight;
                         label_stats.total_ends += weight;
                    } else {
                         label_stats.total_ends += weight;
                    }
                }
            });
        }
    }

    fn merge(&mut self, other: Self) {
        // Merge off-target stats
        for (i, other_s) in other.stats.iter().enumerate() {
            self.stats[i].merge(other_s);
        }
        // Merge on-target stats
        for (i, other_s) in other.stats_on.iter().enumerate() {
            self.stats_on[i].merge(other_s);
        }
    }
}
