use pyo3::prelude::*;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::collections::HashMap;


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
    let consumer = FsdConsumer::new(regions, &mut chrom_map);
    
    // 3. Process
    let analyzer = FragmentAnalyzer::new(consumer, 100_000);
    let final_consumer = analyzer.process_file(&bed_path, &mut chrom_map)
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

#[derive(Clone)]
pub struct FsdConsumer {
    // Shared state
    trees: Arc<HashMap<u32, COITree<usize, u32>>>,
    regions: Arc<Vec<Region>>, // To retrieve original region info (name/chrom) for output
    
    // Thread-local state
    // Index matches regions index
    // Histogram: 67 bins (65-400, step 5)
    histograms: Vec<Vec<u32>>, 
    totals: Vec<u32>,
}

impl FsdConsumer {
    pub fn new(regions: Vec<Region>, chrom_map: &mut ChromosomeMap) -> Self {
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        let mut histograms = Vec::with_capacity(regions.len());
        let mut totals = Vec::with_capacity(regions.len());
        
        for (i, region) in regions.iter().enumerate() {
            // Init counters
            histograms.push(vec![0; 67]);
            totals.push(0);
            
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
            histograms,
            totals,
        }
    }
    
    pub fn write_output(&self, output_path: &Path) -> Result<()> {
        let mut file = File::create(output_path)
            .with_context(|| format!("Failed to create output file: {:?}", output_path))?;
            
        // Header
        let mut header_cols = Vec::new();
        header_cols.push("region".to_string());
        for s in (65..400).step_by(5) {
            header_cols.push(format!("{}-{}", s, s + 4));
        }
        writeln!(file, "{}", header_cols.join("\t"))?;
        
        for (i, region) in self.regions.iter().enumerate() {
             let region_str = format!("{}:{}-{}", region.chrom, region.start, region.end);
             let mut row = Vec::new();
             row.push(region_str);
             
             let total = self.totals[i] as f64;
             for count in &self.histograms[i] {
                 if total > 0.0 {
                     row.push(format!("{:.6}", *count as f64 / total));
                 } else {
                     row.push("0.0".to_string());
                 }
             }
             writeln!(file, "{}", row.join("\t"))?;
        }
        
        Ok(())
    }
}

impl FragmentConsumer for FsdConsumer {
    fn name(&self) -> &str {
        "FSD"
    }

    fn consume(&mut self, fragment: &Fragment) {
        // Binning Logic
        // Length check [65, 400)
        let len = fragment.length;
        if len >= 65 && len < 400 {
             let bin_idx = ((len - 65) / 5) as usize;
             if bin_idx >= 67 { return; } // Should be covered by check but safety first
             
             if let Some(tree) = self.trees.get(&fragment.chrom_id) {
                 let start = fragment.start as i32;
                 let end_closed = if fragment.end > fragment.start { (fragment.end - 1) as i32 } else { start };
                 
                 tree.query(start, end_closed, |node| {
                     let idx = node.metadata.to_owned();
                     // Safety: idx comes from initialization which matches histograms size
                     if let Some(hist) = self.histograms.get_mut(idx) {
                         hist[bin_idx] += 1;
                         self.totals[idx] += 1;
                     }
                 });
             }
        }
    }

    fn merge(&mut self, other: Self) {
        for (i, (my_hist, other_hist)) in self.histograms.iter_mut().zip(other.histograms.iter()).enumerate() {
            for (my_bin, other_bin) in my_hist.iter_mut().zip(other_hist.iter()) {
                *my_bin += *other_bin;
            }
            self.totals[i] += other.totals[i];
        }
    }
}

/// Parse Arms/Regions file (Chrom Start End)
pub fn parse_regions_file(path: &Path) -> Result<Vec<Region>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut regions = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 { continue; }
        let chrom = fields[0].to_string();
        let start: u64 = fields[1].parse().unwrap_or(0);
        let end: u64 = fields[2].parse().unwrap_or(0);
        regions.push(Region::new(chrom, start, end));
    }
    Ok(regions)
}
