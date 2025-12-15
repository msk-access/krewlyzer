use pyo3::prelude::*;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
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
    let consumer = OcfConsumer::new(&ocr_path, &mut chrom_map)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        
    let analyzer = FragmentAnalyzer::new(consumer, 100_000);
    let final_consumer = analyzer.process_file(&bed_path, &mut chrom_map)
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

#[derive(Clone, Default)]
struct LabelStats {
    left_pos: Vec<u64>,  // Size 2000
    right_pos: Vec<u64>, // Size 2000
    total_starts: u64,
    total_ends: u64,
}

impl LabelStats {
    fn new() -> Self {
        Self {
            left_pos: vec![0; 2000],
            right_pos: vec![0; 2000],
            total_starts: 0,
            total_ends: 0,
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

#[derive(Clone)]
pub struct OcfConsumer {
    // Shared state
    trees: Arc<HashMap<u32, COITree<usize, u32>>>,
    region_infos: Arc<Vec<OcrRegionInfo>>, // Index -> Info
    labels: Arc<Vec<String>>, // LabelID -> Name
    
    // Thread-local state
    // Indexed by LabelID
    stats: Vec<LabelStats>,
}

impl OcfConsumer {
    pub fn new(ocr_path: &Path, chrom_map: &mut ChromosomeMap) -> Result<Self> {
        // 1. Parse OCR File
        // Format: chrom start end label
        let file = File::open(ocr_path)
            .with_context(|| format!("Failed to open OCR file: {:?}", ocr_path))?;
        let reader = BufReader::new(file);
        
        let mut label_to_id: HashMap<String, usize> = HashMap::new();
        let mut labels: Vec<String> = Vec::new();
        
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        let mut region_infos = Vec::new(); // Implicitly indexed by order of insertion
        
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 4 { continue; }
            
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
        
        // Build Trees
        let mut trees = HashMap::new();
        for (cid, nodes) in nodes_by_chrom {
            trees.insert(cid, COITree::new(&nodes));
        }
        
        // Init stats
        let stats = vec![LabelStats::new(); labels.len()];
        
        Ok(Self {
            trees: Arc::new(trees),
            region_infos: Arc::new(region_infos),
            labels: Arc::new(labels),
            stats,
        })
    }
    
    pub fn write_output(&self, output_dir: &Path) -> Result<()> {
        // Output 1: all.ocf.csv
        let summary_path = output_dir.join("all.ocf.tsv");
        let mut summary_file = File::create(&summary_path)
            .with_context(|| format!("Failed to create summary file: {:?}", summary_path))?;
        writeln!(summary_file, "tissue\tOCF")?;

        // Output 2: all.sync.tsv
        let sync_path = output_dir.join("all.sync.tsv");
        let mut sync_file = File::create(&sync_path)
             .with_context(|| format!("Failed to create sync file: {:?}", sync_path))?;
        writeln!(sync_file, "tissue\tposition\tleft_count\tleft_norm\tright_count\tright_norm")?;

        // Sort labels
        // We have self.labels (Vec<String>) indexed by ID.
        // We want output sorted by Label Name.
        let mut label_indices: Vec<usize> = (0..self.labels.len()).collect();
        label_indices.sort_by_key(|&i| &self.labels[i]);
        
        for &id in &label_indices {
            let label = &self.labels[id];
            let s = &self.stats[id];
            
            let ts = if s.total_starts > 0 { s.total_starts as f64 / 10000.0 } else { 1.0 };
            let te = if s.total_ends > 0 { s.total_ends as f64 / 10000.0 } else { 1.0 };
            
            let peak = 60;
            let bin_width = 10;
            let mut trueends = 0.0;
            let mut background = 0.0;
            
            for k in 0..2000 {
                let l_count = s.left_pos[k] as f64;
                let r_count = s.right_pos[k] as f64;
                let l_norm = l_count / ts;
                let r_norm = r_count / te;
                let loc = k as i64 - 1000;
                
                writeln!(sync_file, "{}\t{}\t{}\t{:.6}\t{}\t{:.6}", 
                    label, loc, l_count as u64, l_norm, r_count as u64, r_norm)?;
                
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
        if let Some(tree) = self.trees.get(&fragment.chrom_id) {
            let start = fragment.start as i32;
            let end_closed = if fragment.end > fragment.start { (fragment.end - 1) as i32 } else { start };
            
            tree.query(start, end_closed, |node| {
                let region_idx = node.metadata.to_owned();
                let region = &self.region_infos[region_idx];
                let label_stats = &mut self.stats[region.label_id];
                
                // Logic per Legacy:
                let r_start = fragment.start; // u64
                let r_end = fragment.end;     // u64
                let reg_start = region.start; // u64
                let reg_end = region.end;     // u64
                
                // Starts
                if r_start >= reg_start {
                    let s = (r_start - reg_start) as usize;
                    if s < 2000 {
                        label_stats.left_pos[s] += 1;
                        label_stats.total_starts += 1;
                    } else {
                        // "else { total_starts++ }" implies we count even if out of window,
                        // IF the read overlaps the region.
                        // (Legacy logic confirmed this branch exists)
                        label_stats.total_starts += 1;
                    }
                }
                
                // Ends
                if r_end <= reg_end {
                    // Legacy: e = (r_end - region.start + 1)
                    // r_end is u64.
                    let diff = r_end.wrapping_sub(reg_start).wrapping_add(1);
                    // Check logic: if r_end < reg_start-1, diff wraps?
                    // But we are inside `if r_end <= reg_end`.
                    // Is it possible `r_end < reg_start`?
                    // Fragment overlaps Region.
                    // Overlap: `r_start < reg_end` and `r_end > reg_start`.
                    // So `r_end > reg_start`. Thus `r_end - reg_start` >= 1.
                    // So `diff` >= 2? (since +1).
                    // Or if r_end == reg_start + 1? diff = 2.
                    // Wait, if r_end > reg_start. Minimally r_end = reg_start + 1.
                    // So diff >= 2.
                    // Legacy check: `if e >= 0 && e < 2000`. `e` was likely i64.
                    // Here we used u64 for starts.
                    // Let's safe cast.
                    let e = diff as usize;
                    
                    if e < 2000 {
                         label_stats.right_pos[e] += 1;
                         label_stats.total_ends += 1;
                    } else {
                         label_stats.total_ends += 1;
                    }
                }
            });
        }
    }

    fn merge(&mut self, other: Self) {
        for (i, other_s) in other.stats.iter().enumerate() {
            self.stats[i].merge(other_s);
        }
    }
}
