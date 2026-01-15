//! Fragment Size Distribution (FSD) calculation
//!
//! Computes fragment size histograms per chromosome arm for cancer biomarker analysis.
//! Key features:
//! - **Size bins**: 65bp to 400bp (1bp resolution)
//! - **Arm-level aggregation**: Separate histograms for p/q arms
//! - **GC correction**: Per-fragment weight adjustment via correction factors
//! - **Panel mode**: On-target and off-target split outputs
//!
//! Output: TSV with rows=arms (e.g., chr1p, chr1q), columns=fragment counts per size.

use pyo3::prelude::*;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, Write};
use std::collections::HashMap;
use crate::bed;
use log::info;


/// Calculate Fragment Size Distribution (FSD)
/// 
/// # Arguments
/// * `bed_path` - Path to the input .bed.gz file (from motif)
/// * `arms_path` - Path to the chromosome arms BED file
/// * `output_path` - Path to the output TSV file
#[pyfunction]
pub fn calculate_fsd(
    bed_path: PathBuf,
    arms_path: PathBuf,
    output_path: PathBuf,
) -> PyResult<()> {
    // 1. Parse Arms
    let regions = parse_regions_file(&arms_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    // 2. Setup Engine
    let mut chrom_map = ChromosomeMap::new();
    let consumer = FsdConsumer::new(regions, &mut chrom_map, None, None);
    
    // 3. Process
    let analyzer = FragmentAnalyzer::new(consumer, 100_000);
    let final_consumer = analyzer.process_file(&bed_path, &mut chrom_map, false)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        
    // 4. Write Output
    final_consumer.write_output(&output_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    Ok(())
}

use crate::bed::{Fragment, ChromosomeMap, Region};
use crate::engine::{FragmentConsumer, FragmentAnalyzer};
use std::sync::Arc;
use std::path::Path;
use anyhow::{Result, Context};
use coitrees::{COITree, IntervalNode, IntervalTree};

use crate::gc_correction::CorrectionFactors;

/// FSD Consumer with on-target/off-target split support
/// 
/// When target_regions are provided (panel data like MSK-ACCESS),
/// fragments are split into two separate histograms:
/// - Off-target: Unbiased global signal (written to main output)
/// - On-target: Capture-biased local signal (written to .ontarget.tsv)
#[derive(Clone)]
pub struct FsdConsumer {
    // Shared state
    trees: Arc<HashMap<u32, COITree<usize, u32>>>,
    regions: Arc<Vec<Region>>,
    factors: Option<Arc<CorrectionFactors>>,
    
    // Target regions for on/off-target split (panel data)
    target_tree: Option<Arc<HashMap<u32, COITree<(), u32>>>>,
    
    // Thread-local state: 67 bins (65-400, step 5)
    // Off-target histograms (default output)
    histograms_off: Vec<Vec<f64>>, 
    totals_off: Vec<f64>,
    
    // On-target histograms (only populated if target_tree is Some)
    histograms_on: Vec<Vec<f64>>,
    totals_on: Vec<f64>,
}

impl FsdConsumer {
    pub fn new(
        regions: Vec<Region>, 
        chrom_map: &mut ChromosomeMap, 
        factors: Option<Arc<CorrectionFactors>>,
        target_regions: Option<Arc<HashMap<u32, COITree<(), u32>>>>
    ) -> Self {
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        let n_regions = regions.len();
        
        for (i, region) in regions.iter().enumerate() {
            // Map chromosome
            let chrom_norm = region.chrom.trim_start_matches("chr");
            let chrom_id = chrom_map.get_id(chrom_norm);
            
            // Add to tree
            let start = region.start as u32;
            let end = region.end as u32; 
            let end_closed = if end > start { end - 1 } else { start };
            
            nodes_by_chrom.entry(chrom_id).or_default().push(
                IntervalNode::new(start as i32, end_closed as i32, i)
            );
        }
        
        let mut trees = HashMap::new();
        for (chrom_id, nodes) in nodes_by_chrom {
            trees.insert(chrom_id, COITree::new(&nodes));
        }
        
        Self {
            trees: Arc::new(trees),
            regions: Arc::new(regions),
            factors,
            target_tree: target_regions,
            histograms_off: vec![vec![0.0; 67]; n_regions],
            totals_off: vec![0.0; n_regions],
            histograms_on: vec![vec![0.0; 67]; n_regions],
            totals_on: vec![0.0; n_regions],
        }
    }
    
    /// Check if fragment overlaps any target region
    fn is_on_target(&self, fragment: &Fragment) -> bool {
        if let Some(ref target_tree) = self.target_tree {
            if let Some(tree) = target_tree.get(&fragment.chrom_id) {
                let start = fragment.start as i32;
                let end_closed = if fragment.end > fragment.start { 
                    (fragment.end - 1) as i32 
                } else { 
                    start 
                };
                
                let mut found = false;
                tree.query(start, end_closed, |_| { found = true; });
                return found;
            }
        }
        false
    }
    
    /// Write histogram to file (raw GC-weighted counts)
    fn write_histogram(&self, histograms: &[Vec<f64>], totals: &[f64], output_path: &Path) -> Result<()> {
        let mut file = File::create(output_path)
            .with_context(|| format!("Failed to create output file: {:?}", output_path))?;
            
        // Header: region + 67 bin columns + total
        let mut header_cols = vec!["region".to_string()];
        for s in (65..400).step_by(5) {
            header_cols.push(format!("{}-{}", s, s + 4));
        }
        header_cols.push("total".to_string());
        writeln!(file, "{}", header_cols.join("\t"))?;
        
        for (i, region) in self.regions.iter().enumerate() {
            let region_str = format!("{}:{}-{}", region.chrom, region.start, region.end);
            let mut row = vec![region_str];
            
            // Output raw GC-weighted counts (not frequencies)
            for count in &histograms[i] {
                row.push(format!("{:.4}", *count));
            }
            row.push(format!("{:.4}", totals[i]));
            
            writeln!(file, "{}", row.join("\t"))?;
        }
        
        Ok(())
    }
    
    pub fn write_output(&self, output_path: &Path) -> Result<()> {
        // Log summary
        let total_off: f64 = self.totals_off.iter().sum();
        let total_on: f64 = self.totals_on.iter().sum();
        
        info!("FSD: {} arms, {:.0} off-target frags, {:.0} on-target frags", 
            self.regions.len(), total_off, total_on);
        
        // Write off-target (main output)
        self.write_histogram(&self.histograms_off, &self.totals_off, output_path)?;
        
        // Write on-target if targets were provided
        if self.target_tree.is_some() && total_on > 0.0 {
            let stem = output_path.file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("output");
            let on_path = output_path.with_file_name(format!("{}.ontarget.tsv", stem));
            self.write_histogram(&self.histograms_on, &self.totals_on, &on_path)?;
            info!("FSD on-target: {:?}", on_path);
        }
        
        Ok(())
    }
}

impl FragmentConsumer for FsdConsumer {
    fn name(&self) -> &str {
        "FSD"
    }

    fn consume(&mut self, fragment: &Fragment) {
        // Length check [65, 400)
        let len = fragment.length;
        if len >= 65 && len < 400 {
            let bin_idx = ((len - 65) / 5) as usize;
            if bin_idx >= 67 { return; }
            
            if let Some(tree) = self.trees.get(&fragment.chrom_id) {
                let start = fragment.start as i32;
                let end_closed = if fragment.end > fragment.start { 
                    (fragment.end - 1) as i32 
                } else { 
                    start 
                };
                
                // Compute GC weight once
                let gc_pct = (fragment.gc * 100.0).round() as u8;
                let weight = if let Some(ref factors) = self.factors {
                    factors.get_factor(len, gc_pct)
                } else { 1.0 };
                
                // Check if on-target
                let on_target = self.is_on_target(fragment);
                
                tree.query(start, end_closed, |node| {
                    let idx = node.metadata.to_owned();
                    
                    if on_target {
                        if let Some(hist) = self.histograms_on.get_mut(idx) {
                            hist[bin_idx] += weight;
                            self.totals_on[idx] += weight;
                        }
                    } else {
                        if let Some(hist) = self.histograms_off.get_mut(idx) {
                            hist[bin_idx] += weight;
                            self.totals_off[idx] += weight;
                        }
                    }
                });
            }
        }
    }

    fn merge(&mut self, other: Self) {
        for i in 0..self.histograms_off.len() {
            for (my_bin, other_bin) in self.histograms_off[i].iter_mut().zip(other.histograms_off[i].iter()) {
                *my_bin += *other_bin;
            }
            self.totals_off[i] += other.totals_off[i];
            
            for (my_bin, other_bin) in self.histograms_on[i].iter_mut().zip(other.histograms_on[i].iter()) {
                *my_bin += *other_bin;
            }
            self.totals_on[i] += other.totals_on[i];
        }
    }
}

/// Parse Arms/Regions file (Chrom Start End) - supports both plain and BGZF-compressed
pub fn parse_regions_file(path: &Path) -> Result<Vec<Region>> {
    let reader = bed::get_reader(path)?;
    
    let mut regions = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(region) = parse_region_line(&line) {
            regions.push(region);
        }
    }
    Ok(regions)
}

/// Parse a single region line
fn parse_region_line(line: &str) -> Option<Region> {
    let line = line.trim();
    if line.is_empty() || line.starts_with('#') {
        return None;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 3 {
        return None;
    }
    let chrom = fields[0].to_string();
    let start: u64 = fields[1].parse().ok()?;
    let end: u64 = fields[2].parse().ok()?;
    Some(Region::new(chrom, start, end))
}

/// Build target regions tree from BED file
/// Used for on/off-target split in panel data
pub fn load_target_regions(path: &Path, chrom_map: &mut ChromosomeMap) -> Result<HashMap<u32, COITree<(), u32>>> {
    let reader = bed::get_reader(path)?;
    let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<(), u32>>> = HashMap::new();
    
    for line in reader.lines() {
        let line = line?;
        if let Some(region) = parse_region_line(&line) {
            let chrom_norm = region.chrom.trim_start_matches("chr");
            let chrom_id = chrom_map.get_id(chrom_norm);
            
            let start = region.start as i32;
            let end_closed = if region.end > region.start { 
                (region.end - 1) as i32 
            } else { 
                start 
            };
            
            nodes_by_chrom.entry(chrom_id).or_default().push(
                IntervalNode::new(start, end_closed, ())
            );
        }
    }
    
    let mut trees = HashMap::new();
    for (chrom_id, nodes) in nodes_by_chrom {
        trees.insert(chrom_id, COITree::new(&nodes));
    }
    
    info!("Loaded {} target region chromosomes", trees.len());
    Ok(trees)
}
