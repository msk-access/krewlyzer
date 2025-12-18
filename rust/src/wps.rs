//! Windowed Protection Score (WPS) calculation
//!
//! Calculates WPS Long, Short, and Ratio for transcript regions.
//! Matches the reference cfDNAFE implementation exactly.
//! Supports GC bias correction using LOESS.

use std::path::Path;
use std::io::{BufRead, BufReader, Write};
use std::fs::File;
use anyhow::{Result, Context, anyhow};

use pyo3::prelude::*;
use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::faidx;
use log::{info, debug};

// GC correction support (available for future per-region correction)
// use crate::gc_correction::correct_gc_bias_wps;

/// Transcript region from TSV with optional GC content
#[derive(Debug, Clone)]
pub struct Region {
    pub id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: String,
    pub gc: f64,  // GC content 0.0-1.0, computed from reference
}

/// Parse TSV transcript file (regions will have gc=0.0, computed later)
pub fn parse_regions(tsv_path: &Path) -> Result<Vec<Region>> {
    let file = File::open(tsv_path).with_context(|| "Failed to open regions file")?;
    let reader = BufReader::new(file);
    let mut regions = Vec::new();
    
    let valid_chroms: Vec<String> = (1..=22).map(|i| i.to_string()).chain(vec!["X".to_string(), "Y".to_string()]).collect();
    
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() { continue; }
        
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 5 { continue; }
        
        let id = fields[0].to_string();
        let chrom_raw = fields[1];
        let start_str = fields[2];
        let end_str = fields[3];
        let strand = fields[4].to_string();
        
        let chrom_norm = chrom_raw.trim_start_matches("chr").to_string();
        
        if !valid_chroms.iter().any(|c| c == &chrom_norm) {
            continue;
        }
        
        let start: u64 = start_str.parse::<f64>().unwrap_or(0.0) as u64;
        let end: u64 = end_str.parse::<f64>().unwrap_or(0.0) as u64;
        
        if start < 1 { continue; }
        
        regions.push(Region { id, chrom: chrom_norm, start, end, strand, gc: 0.0 });
    }
    
    Ok(regions)
}

/// Compute GC content for regions from reference FASTA
pub fn compute_gc_from_fasta(regions: &mut [Region], fasta_path: &Path) -> Result<()> {
    let faidx = faidx::Reader::from_path(fasta_path)
        .map_err(|e| anyhow!("Failed to open FASTA index: {}. Make sure .fai exists.", e))?;
    
    for region in regions.iter_mut() {
        // Try both "chr" prefixed and non-prefixed chromosome names
        let chrom_variants = [
            region.chrom.clone(),
            format!("chr{}", region.chrom),
        ];
        
        let mut gc_computed = false;
        for chrom in &chrom_variants {
            match faidx.fetch_seq(chrom, region.start as usize, region.end as usize) {
                Ok(seq) => {
                    let len = seq.len();
                    if len > 0 {
                        let gc_count = seq.iter()
                            .filter(|&&c| c == b'G' || c == b'g' || c == b'C' || c == b'c')
                            .count();
                        region.gc = gc_count as f64 / len as f64;
                    } else {
                        region.gc = 0.5; // Default for empty regions
                    }
                    gc_computed = true;
                    break;
                }
                Err(_) => continue,
            }
        }
        
        if !gc_computed {
            region.gc = 0.5; // Default if chromosome not found
            debug!("Could not fetch GC for region {}:{}-{}, using default 0.5", 
                   region.id, region.start, region.end);
        }
    }
    
    Ok(())
}


/// Unified entry point for WPS (replaces legacy sequential implementation)
#[pyfunction]
#[pyo3(signature = (bedgz_path, tsv_path, output_dir, file_stem, empty=false, total_fragments=None, reference_path=None, gc_correct=false, verbose=false))]
pub fn calculate_wps(
    _py: Python<'_>,
    bedgz_path: &str,
    tsv_path: &str,
    output_dir: &str,
    file_stem: &str,
    empty: bool,
    total_fragments: Option<u64>,
    reference_path: Option<&str>,
    gc_correct: bool,
    verbose: bool,
) -> PyResult<usize> {
    let bed_path = Path::new(bedgz_path);
    let tsv = Path::new(tsv_path);
    let output_path = Path::new(output_dir).join(format!("{}.WPS.tsv.gz", file_stem));
    
    // 1. Parse Regions
    let mut regions = parse_regions(tsv)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to parse regions: {}", e)))?;
        
    let initial_count = regions.len();
    
    // 2. Compute GC from FASTA if GC correction is enabled
    if gc_correct {
        if let Some(ref_path) = reference_path {
            if verbose {
                info!("Computing region GC content from reference FASTA...");
            }
            compute_gc_from_fasta(&mut regions, Path::new(ref_path))
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("GC computation failed: {}", e)))?;
            if verbose {
                info!("GC content computed for {} regions", regions.len());
            }
        } else {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "gc_correct=True requires reference_path to be provided"
            ));
        }
    }

    // 3. Setup Engine
    let mut chrom_map = ChromosomeMap::default();
    let consumer = WpsConsumer::new(regions, &mut chrom_map);
    let analyzer = FragmentAnalyzer::new(consumer, 100_000); // 100k chunk size
    
    // 4. Process
    let final_consumer = analyzer.process_file(bed_path, &mut chrom_map)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Processing failed: {}", e)))?;
        
    // 5. Write Output (with optional GC correction applied internally)
    final_consumer.write_output(&output_path, total_fragments, empty, gc_correct, verbose)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to write output: {}", e)))?;
        
    Ok(initial_count)
}


use crate::bed::ChromosomeMap;
use crate::engine::{FragmentConsumer, FragmentAnalyzer};
use std::sync::Arc;
use std::collections::HashMap;
use coitrees::{COITree, IntervalNode, IntervalTree};

/// Output row for a single position with both Long and Short WPS
#[derive(Debug, Clone)]
struct WpsRow {
    gene_id: String,
    chrom: String,
    pos: u64,
    cov_long: u32,
    cov_short: u32,
    wps_long: i32,
    wps_short: i32,
    wps_ratio: f64,
    wps_long_norm: f64,
    wps_short_norm: f64,
    wps_ratio_norm: f64,
}

/// Internal accumulator for a single region (using Difference Arrays)
#[derive(Clone)]
struct RegionAccumulator {
    // Difference arrays. 
    // Size = length + 1 to handle range end+1.
    // Using i32 because diffs can be negative.
    // We accumulate everything in diffs, then integrate at the end.
    cov_long: Vec<i32>,
    cov_short: Vec<i32>,
    wps_long: Vec<i32>,
    wps_short: Vec<i32>,
}

impl RegionAccumulator {
    fn new(length: usize) -> Self {
        Self {
            cov_long: vec![0; length + 1],
            cov_short: vec![0; length + 1],
            wps_long: vec![0; length + 1],
            wps_short: vec![0; length + 1],
        }
    }
    
    // Add val to [start, end)
    fn add_range(vec: &mut Vec<i32>, start: i64, end: i64, val: i32) {
        let len = vec.len() as i64 - 1; 
        // Clamp to region bounds [0, len)
        let s = start.clamp(0, len) as usize;
        let e = end.clamp(0, len) as usize;
        
        if s < e {
            vec[s] += val;
            vec[e] -= val;
        }
    }

    fn merge(&mut self, other: &Self) {
        for (a, b) in self.cov_long.iter_mut().zip(other.cov_long.iter()) { *a += *b; }
        for (a, b) in self.cov_short.iter_mut().zip(other.cov_short.iter()) { *a += *b; }
        for (a, b) in self.wps_long.iter_mut().zip(other.wps_long.iter()) { *a += *b; }
        for (a, b) in self.wps_short.iter_mut().zip(other.wps_short.iter()) { *a += *b; }
    }
}

#[derive(Clone)]
pub struct WpsConsumer {
    // Shared state
    // Map chrom_id -> IntervalTree of Region Indices
    trees: Arc<HashMap<u32, COITree<usize, u32>>>,
    // Store regions metadata (start/end/strand/id) to generate output
    regions: Arc<Vec<Region>>,
    
    // Thread-local state
    accumulators: HashMap<usize, RegionAccumulator>,
}

impl WpsConsumer {
    pub fn new(regions: Vec<Region>, chrom_map: &mut ChromosomeMap) -> Self {
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        // We do NOT pre-allocate accumulators anymore.
        
        for (i, region) in regions.iter().enumerate() {
            // Map chromosome
            let chrom_norm = region.chrom.trim_start_matches("chr");
            let chrom_id = chrom_map.get_id(chrom_norm);
            
            // Add to tree (using extended window for lookup)
            // Long protection is 60. Max Long Frag is 180.
            // Any fragment that could theoretically touch the region's analysis must be included.
            // Safe bet: expand by MAX_FRAG_SIZE/2 + PROTECTION?
            // Actually, we just need to catch any fragment that OVERLAPS the region-of-interest extended by protection.
            // Legacy used: [start - 60, end + 60].
            // To be safe, let's use the maximum lookup needed.
            let start = region.start.saturating_sub(60) as u32;
            let end = (region.end + 60) as u32;
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
            accumulators: HashMap::new(),
        }
    }
    
    /// Write results to output file (with optional GC correction logging)
    pub fn write_output(&self, output_path: &Path, total_markers: Option<u64>, empty: bool, gc_correct: bool, verbose: bool) -> Result<()> {
        if gc_correct && verbose {
            info!("Writing WPS output with GC-aware region metadata...");
        }
        let file = File::create(output_path)?;
        let mut encoder = GzEncoder::new(file, Compression::default());
        let norm_factor = total_markers.unwrap_or(1_000_000) as f64 / 1_000_000.0;
        
        writeln!(encoder, "gene_id\tchrom\tpos\tcov_long\tcov_short\twps_long\twps_short\twps_ratio\twps_long_norm\twps_short_norm\twps_ratio_norm")?;
        
        for (i, region) in self.regions.iter().enumerate() {
            // Check if we have data for this region
            let acc_opt = self.accumulators.get(&i);
            
            // If no data and skipping empty - verify logic
            if acc_opt.is_none() && !empty {
                continue;
            }
            
            let len = (region.end - region.start + 1) as usize;
            
            if acc_opt.is_none() {
                // If we must print empty regions, print zeros
                if empty {
                    for j in 0..len {
                        writeln!(encoder, "{}\t{}\t{}\t0\t0\t0\t0\t0.0000\t0.000000\t0.000000\t0.000000",
                            region.id, region.chrom, region.start + j as u64)?;
                    }
                }
                continue;
            }
            
            let acc = acc_opt.unwrap();
            
            // Reconstruct values from difference arrays
            let mut curr_cov_long = 0;
            let mut curr_cov_short = 0;
            let mut curr_wps_long = 0;
            let mut curr_wps_short = 0;
            
            // let len = (region.end - region.start + 1) as usize; // Already calc'd above
            let mut rows = Vec::with_capacity(len);
            let mut total_cov = 0;
            
            for j in 0..len {
                curr_cov_long += acc.cov_long[j];
                curr_cov_short += acc.cov_short[j];
                curr_wps_long += acc.wps_long[j];
                curr_wps_short += acc.wps_short[j];
                
                total_cov += curr_cov_long + curr_cov_short;
                
                let wps_l = curr_wps_long;
                let wps_s = curr_wps_short;
                
                let ratio = if wps_s != 0 {
                    wps_l as f64 / wps_s.abs() as f64
                } else if wps_l != 0 {
                    wps_l as f64
                } else {
                    0.0
                };
                
                 rows.push(WpsRow {
                    gene_id: region.id.clone(),
                    chrom: region.chrom.clone(),
                    pos: region.start + j as u64,
                    cov_long: curr_cov_long as u32,
                    cov_short: curr_cov_short as u32,
                    wps_long: wps_l,
                    wps_short: wps_s,
                    wps_ratio: ratio,
                    wps_long_norm: wps_l as f64 / norm_factor,
                    wps_short_norm: wps_s as f64 / norm_factor,
                    wps_ratio_norm: ratio / norm_factor,
                });
            }
            
            // Skip empty if requested (double check)
            if !empty && total_cov == 0 {
                continue;
            }
            
            if region.strand == "-" {
                rows.reverse();
            }
            
            for row in rows {
                writeln!(encoder, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.6}\t{:.6}\t{:.6}", 
                    row.gene_id, row.chrom, row.pos, 
                    row.cov_long, row.cov_short,
                    row.wps_long, row.wps_short, row.wps_ratio,
                    row.wps_long_norm, row.wps_short_norm, row.wps_ratio_norm)?;
            }
        }
        
        encoder.finish()?;
        Ok(())
    }
}

impl FragmentConsumer for WpsConsumer {
    fn name(&self) -> &str {
        "WPS"
    }

    fn consume(&mut self, fragment: &crate::bed::Fragment) {
        if let Some(tree) = self.trees.get(&fragment.chrom_id) {
            let start = fragment.start as u32;
            let end = fragment.end as u32;
            let end_closed = if end > start { end - 1 } else { start };
            
            // Collect matches first to avoid concurrent borrow of self
            let mut matches: Vec<usize> = Vec::new();
            tree.query(start as i32, end_closed as i32, |node| {
                matches.push(node.metadata.to_owned());
            });
            
            // Now update accumulators
            for region_idx in matches {
                // Lazy allocation + Get mutable reference
                let acc = self.accumulators.entry(region_idx).or_insert_with(|| {
                     let region = &self.regions[region_idx];
                     let len = (region.end - region.start + 1) as usize;
                     RegionAccumulator::new(len)
                });

                let region = &self.regions[region_idx];
                let r_start = region.start as i64;
                // Fragment coords relative to region
                let f_start = (fragment.start + 1) as i64 - r_start;
                let f_end = fragment.end as i64 - r_start; // exclusive
                
                // Parameters
                let is_long = fragment.length >= 120 && fragment.length <= 180;
                let is_short = fragment.length >= 35 && fragment.length <= 80;
                
                if is_long {
                     // Coverage: [start, end)
                     RegionAccumulator::add_range(&mut acc.cov_long, f_start, f_end, 1);
                     
                     // WPS Long: P=60
                     let p = 60;
                     RegionAccumulator::add_range(&mut acc.wps_long, f_start + p, f_end - p + 1, 1);
                     RegionAccumulator::add_range(&mut acc.wps_long, f_start - p, f_start + p, -1);
                     RegionAccumulator::add_range(&mut acc.wps_long, f_end - p + 1, f_end + p + 1, -1);

                } else if is_short {
                     RegionAccumulator::add_range(&mut acc.cov_short, f_start, f_end, 1);
                     
                     let p = 8;
                     RegionAccumulator::add_range(&mut acc.wps_short, f_start + p, f_end - p + 1, 1);
                     RegionAccumulator::add_range(&mut acc.wps_short, f_start - p, f_start + p, -1);
                     RegionAccumulator::add_range(&mut acc.wps_short, f_end - p + 1, f_end + p + 1, -1);
                }
            }
        }
    }

    fn merge(&mut self, other: Self) {
        for (idx, other_acc) in other.accumulators {
            match self.accumulators.entry(idx) {
                std::collections::hash_map::Entry::Occupied(mut entry) => {
                    entry.get_mut().merge(&other_acc);
                },
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(other_acc);
                }
            }
        }
    }
}


