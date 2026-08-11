//! Windowed Protection Score (WPS) calculation
//!
//! Calculates WPS for transcript/CTCF anchor regions and Alu background hierarchies.
//!
//! ## Features
//! - **Dual-stream analysis**: WPS-Nuc (nucleosome, 120bp window) and WPS-TF (TF, 16bp window)
//! - **Savitzky-Golay smoothing**: Applied via sci-rs (foreground: window=11, order=3; background: window=7, order=3)
//! - **FFT periodicity extraction**: Nucleosome Repeat Length (NRL) and quality score via realfft
//! - **NRL deviation scoring**: Deviation from expected 190bp with penalty (nrl_deviation_bp, adjusted_score)
//! - **GC bias correction**: LOESS-based correction using correction_factors.csv
//!
//! ## Background WPS Parquet Schema
//! - group_id: Hierarchical group (e.g., Global_All, Chr1_H, AluJb)
//! - stacked_wps_nuc/wps_tf: Aggregated WPS profiles (200 bins)
//! - nrl_bp: Nucleosome Repeat Length in bp (expected ~190bp)
//! - nrl_deviation_bp: |nrl_bp - 190|
//! - periodicity_score: Raw SNR-based quality (0-1)
//! - adjusted_score: periodicity_score × deviation_penalty


use std::path::{Path, PathBuf};
use std::io::{BufRead, Write};
use std::fs::File;
use anyhow::{Result, Context, anyhow};

use pyo3::prelude::*;
use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::faidx;
use log::{info, debug, warn, error};
use rayon::prelude::*;

// GC correction support 
use crate::gc_correction::CorrectionFactors;

/// Configuration for WPS dual-stream analysis
/// 
/// Implements weighted fragment classification for:
/// - WPS-Nuc (Nucleosome): 120bp window, weighted fragments
/// - WPS-TF (Transcription Factor): 16bp window, short fragments
#[derive(Debug, Clone)]
pub struct WpsConfig {
    /// Protection window for nucleosome signal (default: 60bp, total window 120bp)
    pub nuc_window: i64,
    /// Protection window for TF signal (default: 8bp, total window 16bp)
    pub tf_window: i64,
    /// Primary nucleosome fragment range (weight 1.0)
    pub nuc_primary_min: u64,
    pub nuc_primary_max: u64,
    /// Secondary nucleosome fragment range (weight 0.5)
    pub nuc_secondary_min: u64,
    pub nuc_secondary_max: u64,
    /// TF fragment range (weight 1.0)
    pub tf_min: u64,
    pub tf_max: u64,
}

impl Default for WpsConfig {
    fn default() -> Self {
        Self {
            nuc_window: 60,
            tf_window: 8,
            nuc_primary_min: 160,
            nuc_primary_max: 175,
            nuc_secondary_min: 120,
            nuc_secondary_max: 180,
            tf_min: 35,
            tf_max: 80,
        }
    }
}

impl WpsConfig {
    /// Get nucleosome weight for a fragment length
    /// - 1.0 for primary range [160, 175]
    /// - 0.5 for secondary range [120, 159] or [176, 180]
    /// - 0.0 otherwise
    pub fn nuc_weight(&self, length: u64) -> f64 {
        if length >= self.nuc_primary_min && length <= self.nuc_primary_max {
            1.0
        } else if length >= self.nuc_secondary_min && length <= self.nuc_secondary_max {
            0.5
        } else {
            0.0
        }
    }
    
    /// Get TF weight for a fragment length
    /// - 1.0 for [35, 80]
    /// - 0.0 otherwise
    pub fn tf_weight(&self, length: u64) -> f64 {
        if length >= self.tf_min && length <= self.tf_max {
            1.0
        } else {
            0.0
        }
    }
}

/// Anchor region from BED6 file for WPS computation
/// 
/// BED6 format: chrom, start, end, name, score, strand
/// TSS entries have +/- strand (requires vector flipping for ML)
/// CTCF entries have "." strand (symmetric, no flipping)
#[derive(Debug, Clone)]
pub struct Region {
    pub id: String,      // Name from column 4 (e.g., "TSS|TP53|ENST00000xxx" or "CTCF|1:12345")
    pub chrom: String,
    pub start: u64,      // Window start (expanded from center)
    pub end: u64,        // Window end (expanded from center)
    pub strand: String,  // +, -, or . (dot for unstranded)
    pub gc: f64,         // GC content 0.0-1.0, computed from reference
    pub center: u64,     // Original center point from BED
    pub is_minus_strand: bool, // For ML vector flipping
}

impl Region {
    /// Check if this region needs vector reversal for ML
    /// TSS on minus strand should have WPS vector reversed
    pub fn should_reverse_for_ml(&self) -> bool {
        self.strand == "-"
    }
}

/// Parse BED6 anchor file with window expansion
/// 
/// BED6: chrom, start, end, name, score, strand
/// Input is 1bp centered regions, expanded to ±window_bp
pub fn parse_regions(tsv_path: &Path) -> Result<Vec<Region>> {
    parse_regions_with_window(tsv_path, 1000) // Default ±1000bp window
}

/// Parse BED6 with configurable window size
/// Supports both plain text and BGZF-compressed (.gz) BED files
pub fn parse_regions_with_window(bed_path: &Path, window_bp: u64) -> Result<Vec<Region>> {
    // Use shared BGZF-aware reader
    let reader = crate::bed::get_reader(bed_path)?;
    
    let mut regions = Vec::new();
    
    let valid_chroms: Vec<String> = (1..=22).map(|i| i.to_string()).chain(vec!["X".to_string(), "Y".to_string()]).collect();
    
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') { continue; }
        
        let fields: Vec<&str> = line.split('\t').collect();
        
        // Support both BED6 and legacy TSV formats
        if fields.len() >= 6 {
            // BED6: chrom, start, end, name, score, strand
            let chrom_raw = fields[0];
            let start: u64 = fields[1].parse().unwrap_or(0);
            let end: u64 = fields[2].parse().unwrap_or(0);
            let id = fields[3].to_string();
            // fields[4] is score (ignored)
            let strand = fields[5].to_string();
            
            let chrom_norm = chrom_raw.trim_start_matches("chr").to_string();
            
            if !valid_chroms.iter().any(|c| c == &chrom_norm) {
                continue;
            }
            
            // Calculate center (midpoint of 1bp window from BED)
            let center = (start + end) / 2;
            
            // Expand window around center
            let window_start = center.saturating_sub(window_bp);
            let window_end = center + window_bp;
            
            let is_minus = strand == "-";
            
            regions.push(Region { 
                id, 
                chrom: chrom_norm, 
                start: window_start, 
                end: window_end, 
                strand, 
                gc: 0.0,
                center,
                is_minus_strand: is_minus,
            });
        } else if fields.len() >= 5 {
            // Legacy TSV: id, chrom, start, end, strand
            let id = fields[0].to_string();
            let chrom_raw = fields[1];
            let start: u64 = fields[2].parse::<f64>().unwrap_or(0.0) as u64;
            let end: u64 = fields[3].parse::<f64>().unwrap_or(0.0) as u64;
            let strand = fields[4].to_string();
            
            let chrom_norm = chrom_raw.trim_start_matches("chr").to_string();
            
            if !valid_chroms.iter().any(|c| c == &chrom_norm) {
                continue;
            }
            
            if start < 1 { continue; }
            
            let center = (start + end) / 2;
            let is_minus = strand == "-";
            
            regions.push(Region { 
                id, 
                chrom: chrom_norm, 
                start, 
                end, 
                strand, 
                gc: 0.0,
                center,
                is_minus_strand: is_minus,
            });
        }
    }
    
    info!("Loaded {} regions from {:?}", regions.len(), bed_path);
    let tss_count = regions.iter().filter(|r| r.id.starts_with("TSS|")).count();
    let ctcf_count = regions.iter().filter(|r| r.id.starts_with("CTCF|")).count();
    if tss_count > 0 || ctcf_count > 0 {
        info!("  TSS: {}, CTCF: {}", tss_count, ctcf_count);
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
    let consumer = WpsConsumer::new(regions, &mut chrom_map, None); // No pre-loaded factors for legacy API
    let analyzer = FragmentAnalyzer::new(consumer, 100_000); // 100k chunk size
    
    // 4. Process
    let final_consumer = analyzer.process_file(bed_path, &mut chrom_map, false)
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

/// Bait mask for panel data - marks positions near capture bait edges as unreliable
/// 
/// When --target-regions is provided with panel data:
/// - Positions inside baits but away from edges: capture_mask = 1 (reliable)
/// - Positions inside baits but near edges (within trim_bp): capture_mask = 0 (unreliable)
/// - Positions outside baits: capture_mask = 0 (unreliable)
#[derive(Clone)]
pub struct BaitMask {
    /// COITree for efficient bait overlap queries
    tree: HashMap<u32, COITree<(u64, u64), u32>>, // Store (start, end) as metadata
    /// Number of bp to trim from each bait edge (default 50)
    trim_bp: u64,
}

impl BaitMask {
    /// Create BaitMask from target regions BED file
    pub fn from_bed(bed_path: &Path, chrom_map: &mut ChromosomeMap, trim_bp: u64) -> Result<Self> {
        use std::io::BufRead;
        
        // Use bed::get_reader to handle both plain and gzipped BED files
        let reader = crate::bed::get_reader(bed_path)
            .with_context(|| format!("Failed to open target regions: {:?}", bed_path))?;
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<(u64, u64), u32>>> = HashMap::new();
        
        let valid_chroms: Vec<String> = (1..=22).map(|i| i.to_string()).chain(vec!["X".to_string(), "Y".to_string()]).collect();
        
        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') { continue; }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 { continue; }
            
            let chrom_raw = fields[0];
            let start: u64 = fields[1].parse().unwrap_or(0);
            let end: u64 = fields[2].parse().unwrap_or(0);
            
            let chrom_norm = chrom_raw.trim_start_matches("chr");
            if !valid_chroms.iter().any(|c| c == chrom_norm) { continue; }
            
            let chrom_id = chrom_map.get_id(chrom_norm);
            let end_closed = if end > start { end - 1 } else { start };
            
            nodes_by_chrom.entry(chrom_id).or_default().push(
                IntervalNode::new(start as i32, end_closed as i32, (start, end))
            );
        }
        
        let mut tree = HashMap::new();
        let total_baits: usize = nodes_by_chrom.values().map(|v| v.len()).sum();
        for (chrom_id, nodes) in nodes_by_chrom {
            tree.insert(chrom_id, COITree::new(&nodes));
        }
        
        info!("BaitMask: Loaded {} target regions, trim_bp={}", total_baits, trim_bp);
        
        Ok(Self { tree, trim_bp })
    }
    
    /// Check if a genomic position is in a reliable (non-edge) region of a bait
    /// 
    /// Returns: (is_on_target, is_reliable)
    /// - is_on_target: true if position overlaps any bait
    /// - is_reliable: true if position is inside a bait AND >= effective_trim away from edges
    /// 
    /// Adaptive Safety: effective_trim = min(user_trim, bait_length / 4)
    /// This ensures we never mask more than 50% of a small exon (25% per side)
    pub fn check_position(&self, chrom_id: u32, pos: u64) -> (bool, bool) {
        if let Some(tree) = self.tree.get(&chrom_id) {
            let mut is_on_target = false;
            let mut is_reliable = false;
            
            tree.query(pos as i32, pos as i32, |node| {
                let (bait_start, bait_end) = (node.metadata.0, node.metadata.1);
                is_on_target = true;
                
                // Adaptive trim: min(user_trim, length/4) to handle small exons
                let bait_length = bait_end.saturating_sub(bait_start);
                let max_safe_trim = bait_length / 4;  // Never trim more than 25% per side
                let effective_trim = std::cmp::min(self.trim_bp, max_safe_trim);
                
                // Check if position is away from edges
                let dist_from_start = pos.saturating_sub(bait_start);
                let dist_from_end = bait_end.saturating_sub(pos);
                
                if dist_from_start >= effective_trim && dist_from_end >= effective_trim {
                    is_reliable = true;
                }
            });
            
            (is_on_target, is_reliable)
        } else {
            (false, false)
        }
    }
}

/// Output row for a single position with both Long and Short WPS
#[derive(Debug, Clone)]
struct WpsRow {
    gene_id: String,
    chrom: String,
    pos: u64,
    strand: String,      // +, -, or .
    region_type: String, // TSS or CTCF
    capture_mask: u8,    // 1 = reliable, 0 = edge/off-target
    cov_long: f64,
    cov_short: f64,
    wps_long: f64,
    wps_short: f64,
    wps_ratio: f64,
    wps_long_norm: f64,
    wps_short_norm: f64,
    wps_ratio_norm: f64,
    prot_frac_nuc: f64,  // Protection fraction for nucleosome: spanning/(spanning+endpoints)
    prot_frac_tf: f64,   // Protection fraction for TF
}

/// Internal accumulator for a single region (using Difference Arrays)
/// 
/// Protection Fraction = spanning / (spanning + endpoints)
/// - spanning: fragments that fully span a position (contribute +1 to WPS)
/// - endpoints: fragments with start/end at a position (contribute -1 to WPS)
#[derive(Clone)]
struct RegionAccumulator {
    // Difference arrays. 
    // Size = length + 1 to handle range end+1.
    // Using f64 because weights are float
    // We accumulate everything in diffs, then integrate at the end.
    cov_long: Vec<f64>,
    cov_short: Vec<f64>,
    wps_long: Vec<f64>,
    wps_short: Vec<f64>,
    
    // For protection fraction calculation
    // Track spanning vs endpoints separately
    spanning_long: Vec<f64>,   // Central protection window
    endpoints_long: Vec<f64>,  // Fragment start/end flanks
    spanning_short: Vec<f64>,
    endpoints_short: Vec<f64>,
}

impl RegionAccumulator {
    fn new(length: usize) -> Self {
        Self {
            cov_long: vec![0.0; length + 1],
            cov_short: vec![0.0; length + 1],
            wps_long: vec![0.0; length + 1],
            wps_short: vec![0.0; length + 1],
            spanning_long: vec![0.0; length + 1],
            endpoints_long: vec![0.0; length + 1],
            spanning_short: vec![0.0; length + 1],
            endpoints_short: vec![0.0; length + 1],
        }
    }
    
    // Add val to [start, end)
    fn add_range(vec: &mut [f64], start: i64, end: i64, val: f64) {
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
        // Protection fraction arrays
        for (a, b) in self.spanning_long.iter_mut().zip(other.spanning_long.iter()) { *a += *b; }
        for (a, b) in self.endpoints_long.iter_mut().zip(other.endpoints_long.iter()) { *a += *b; }
        for (a, b) in self.spanning_short.iter_mut().zip(other.spanning_short.iter()) { *a += *b; }
        for (a, b) in self.endpoints_short.iter_mut().zip(other.endpoints_short.iter()) { *a += *b; }
    }
}

#[derive(Clone)]
pub struct WpsConsumer {
    // Shared state
    // Map chrom_id -> IntervalTree of Region Indices
    trees: Arc<HashMap<u32, COITree<usize, u32>>>,
    // Store regions metadata (start/end/strand/id) to generate output
    regions: Arc<Vec<Region>>,
    // GC Correction factors
    factors: Option<Arc<CorrectionFactors>>,
    // WPS configuration for dual-stream weighted analysis
    config: Arc<WpsConfig>,
    // Bait mask for panel edge detection (optional)
    bait_mask: Option<Arc<BaitMask>>,
    
    // Thread-local state
    accumulators: HashMap<usize, RegionAccumulator>,
    // Fragment counters for logging
    nuc_fragments: u64,
    tf_fragments: u64,
}

impl WpsConsumer {
    pub fn new(regions: Vec<Region>, chrom_map: &mut ChromosomeMap, factors: Option<Arc<CorrectionFactors>>) -> Self {
        Self::with_config(regions, chrom_map, factors, WpsConfig::default())
    }
    
    /// Create WpsConsumer with custom configuration
    pub fn with_config(regions: Vec<Region>, chrom_map: &mut ChromosomeMap, factors: Option<Arc<CorrectionFactors>>, config: WpsConfig) -> Self {
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        
        info!("WPS: Initializing with {} regions", regions.len());
        debug!("WPS config: nuc_window={}, tf_window={}, primary_range=[{},{}]", 
            config.nuc_window, config.tf_window, config.nuc_primary_min, config.nuc_primary_max);
        
        for (i, region) in regions.iter().enumerate() {
            // Map chromosome (handle both chr-prefixed and non-prefixed)
            let chrom_norm = region.chrom.trim_start_matches("chr");
            let chrom_id = chrom_map.get_id(chrom_norm);
            
            // Add to tree (using extended window for lookup)
            // Expand by protection window to catch all overlapping fragments
            let start = region.start.saturating_sub(config.nuc_window as u64) as u32;
            let end = (region.end + config.nuc_window as u64) as u32;
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
            factors,
            config: Arc::new(config),
            bait_mask: None,
            nuc_fragments: 0,
            tf_fragments: 0,
        }
    }
    
    /// Set bait mask for panel edge detection
    pub fn with_bait_mask(mut self, bait_mask: BaitMask) -> Self {
        self.bait_mask = Some(Arc::new(bait_mask));
        self
    }
    
    /// Write results to output file (with optional GC correction)
    /// 
    /// When gc_correct=true:
    /// 1. Compute mean WPS long/short per region  
    /// 2. Fit LOESS: expected_wps = f(region_gc)
    /// 3. Divide each position's WPS by region's expected
    pub fn write_output(&self, output_path: &Path, total_markers: Option<u64>, empty: bool, gc_correct: bool, verbose: bool) -> Result<()> {
        // Log summary before writing
        let regions_with_data = self.accumulators.len();
        let total_regions = self.regions.len();
        info!("WPS: {} of {} regions with coverage ({:.1}%)", 
            regions_with_data, total_regions, 
            (regions_with_data as f64 / total_regions as f64) * 100.0);
        info!("WPS: {} nucleosome fragments, {} TF fragments processed", 
            self.nuc_fragments, self.tf_fragments);
        
        let file = File::create(output_path)?;
        let mut encoder = GzEncoder::new(file, Compression::default());
        let norm_factor = total_markers.unwrap_or(1_000_000) as f64 / 1_000_000.0;
        
        // If GC correction requested, compute expected values per region
        if gc_correct && verbose {
            info!("Legacy WPS GC correction ignored (using upstream fragment weighting).");
        }
        let _gc_correction_factors: Option<Vec<(f64, f64)>> = None;
        
        writeln!(encoder, "gene_id\tchrom\tpos\tstrand\tregion_type\tcapture_mask\tcov_long\tcov_short\twps_long\twps_short\twps_ratio\twps_long_norm\twps_short_norm\twps_ratio_norm\tprot_frac_nuc\tprot_frac_tf")?;
        
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
            let mut curr_cov_long = 0.0;
            let mut curr_cov_short = 0.0;
            let mut curr_wps_long = 0.0;
            let mut curr_wps_short = 0.0;
            // Protection fraction accumulators
            let mut curr_spanning_long = 0.0;
            let mut curr_endpoints_long = 0.0;
            let mut curr_spanning_short = 0.0;
            let mut curr_endpoints_short = 0.0;
            
            let mut rows = Vec::with_capacity(len);
            let mut total_cov = 0.0;
            
            for j in 0..len {
                curr_cov_long += acc.cov_long[j];
                curr_cov_short += acc.cov_short[j];
                curr_wps_long += acc.wps_long[j];
                curr_wps_short += acc.wps_short[j];
                curr_spanning_long += acc.spanning_long[j];
                curr_endpoints_long += acc.endpoints_long[j];
                curr_spanning_short += acc.spanning_short[j];
                curr_endpoints_short += acc.endpoints_short[j];
                
                total_cov += curr_cov_long + curr_cov_short;
                
                // Legacy correction handled upstream
                let wps_l_corrected = curr_wps_long;
                let wps_s_corrected = curr_wps_short;
                
                let ratio = if wps_s_corrected.abs() > 1e-9 {
                    wps_l_corrected / wps_s_corrected.abs()
                } else if wps_l_corrected.abs() > 1e-9 {
                    wps_l_corrected
                } else {
                    0.0
                };
                
                // Determine region type from name prefix
                let region_type = if region.id.starts_with("TSS|") {
                    "TSS".to_string()
                } else if region.id.starts_with("CTCF|") {
                    "CTCF".to_string()
                } else {
                    "OTHER".to_string()
                };
                
                // Compute capture_mask using bait_mask if available
                let capture_mask = if let Some(ref mask) = self.bait_mask {
                    let chrom_id = *self.trees.keys().find(|_| true).unwrap_or(&0); // Get chrom_id
                    let pos = region.start + j as u64;
                    let (_on_target, is_reliable) = mask.check_position(chrom_id, pos);
                    if is_reliable { 1u8 } else { 0u8 }
                } else {
                    1u8 // No panel data = assume all reliable
                };
                
                // Compute protection fraction: spanning / (spanning + endpoints)
                let prot_frac_nuc = {
                    let total = curr_spanning_long + curr_endpoints_long;
                    if total > 1e-9 { curr_spanning_long / total } else { 0.0 }
                };
                let prot_frac_tf = {
                    let total = curr_spanning_short + curr_endpoints_short;
                    if total > 1e-9 { curr_spanning_short / total } else { 0.0 }
                };
                
                rows.push(WpsRow {
                    gene_id: region.id.clone(),
                    chrom: region.chrom.clone(),
                    pos: region.start + j as u64,
                    strand: region.strand.clone(),
                    region_type,
                    capture_mask,
                    cov_long: curr_cov_long,
                    cov_short: curr_cov_short,
                    wps_long: wps_l_corrected,
                    wps_short: wps_s_corrected,
                    wps_ratio: ratio,
                    wps_long_norm: wps_l_corrected / norm_factor,
                    wps_short_norm: wps_s_corrected / norm_factor,
                    wps_ratio_norm: ratio / norm_factor,
                    prot_frac_nuc,
                    prot_frac_tf,
                });
            }
            
            // Skip empty if requested (double check)
            if !empty && total_cov == 0.0 {
                continue;
            }
            
            if region.strand == "-" {
                rows.reverse();
            }
            
            for row in rows {
                writeln!(encoder, "{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{:.4}	{:.6}	{:.6}	{:.6}	{:.4}	{:.4}", 
                    row.gene_id, row.chrom, row.pos, row.strand, row.region_type, row.capture_mask,
                    row.cov_long, row.cov_short,
                    row.wps_long, row.wps_short, row.wps_ratio,
                    row.wps_long_norm, row.wps_short_norm, row.wps_ratio_norm,
                    row.prot_frac_nuc, row.prot_frac_tf)?;
            }
        }
        
        encoder.finish()?;
        Ok(())
    }
    
    /// Write results as ML-ready Parquet with vector columns
    /// 
    /// Instead of per-position rows, aggregates into:
    /// - 200 bins (10bp each from 2000bp window)
    /// - Vector columns: wps_nuc[200], wps_tf[200], capture_mask[200], etc.
    pub fn write_parquet(&self, output_path: &Path, _total_markers: Option<u64>) -> Result<()> {
        use arrow::array::{ArrayRef, Float32Builder, Int32Builder, Int8Builder, ListBuilder, StringBuilder};
        use arrow::datatypes::{DataType, Field, Schema};
        use arrow::record_batch::RecordBatch;
        use parquet::arrow::ArrowWriter;
        use parquet::file::properties::WriterProperties;
        use std::sync::Arc as StdArc;
        
        const NUM_BINS: usize = 200;
        const WINDOW_SIZE: usize = 2000;
        const BIN_SIZE: usize = WINDOW_SIZE / NUM_BINS; // 10bp per bin
        
        info!("WPS: Writing Parquet with {} bins per region", NUM_BINS);
        
        // Define schema - inner List items must be nullable:true to match ListBuilder defaults
        let schema = Schema::new(vec![
            Field::new("region_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("center", DataType::Int32, false),
            Field::new("strand", DataType::Utf8, false),
            Field::new("region_type", DataType::Utf8, false),
            Field::new("wps_nuc", DataType::List(StdArc::new(Field::new("item", DataType::Float32, true))), false),
            Field::new("wps_tf", DataType::List(StdArc::new(Field::new("item", DataType::Float32, true))), false),
            Field::new("capture_mask", DataType::List(StdArc::new(Field::new("item", DataType::Int8, true))), false),
            Field::new("prot_frac_nuc", DataType::List(StdArc::new(Field::new("item", DataType::Float32, true))), false),
            Field::new("prot_frac_tf", DataType::List(StdArc::new(Field::new("item", DataType::Float32, true))), false),
            Field::new("local_depth", DataType::Float32, false),
        ]);
        
        // Builders
        let mut region_id_builder = StringBuilder::new();
        let mut chrom_builder = StringBuilder::new();
        let mut center_builder = Int32Builder::new();
        let mut strand_builder = StringBuilder::new();
        let mut region_type_builder = StringBuilder::new();
        let mut wps_nuc_builder = ListBuilder::new(Float32Builder::new());
        let mut wps_tf_builder = ListBuilder::new(Float32Builder::new());
        let mut mask_builder = ListBuilder::new(Int8Builder::new());
        let mut prot_nuc_builder = ListBuilder::new(Float32Builder::new());
        let mut prot_tf_builder = ListBuilder::new(Float32Builder::new());
        let mut depth_builder = Float32Builder::new();
        
        // =========================================================================
        // PARALLEL REGION PROCESSING
        // =========================================================================
        // Region binning is CPU-bound and embarrassingly parallel.
        // Each region independently: integrates diff arrays → bins → smooths.
        // We collect results, then sequentially write to Parquet builders.
        // =========================================================================
        
        use std::time::Instant;
        let parallel_start = Instant::now();
        
        /// Result from processing a single region (collected for sequential write)
        #[derive(Debug)]
        struct RegionResult {
            region_idx: usize,
            region_id: String,
            chrom: String,
            center: i32,
            strand: String,
            region_type: String,
            binned_wps_nuc: Vec<f32>,
            binned_wps_tf: Vec<f32>,
            binned_mask: Vec<i8>,
            binned_prot_nuc: Vec<f32>,
            binned_prot_tf: Vec<f32>,
            local_depth: f32,
        }
        
        // Clone shared references for parallel access
        let regions = &self.regions;
        let accumulators = &self.accumulators;
        let bait_mask = &self.bait_mask;
        
        // Parallel computation of all region binned vectors
        let results: Vec<Option<RegionResult>> = regions
            .par_iter()
            .enumerate()
            .map(|(i, region)| {
                // Skip regions without accumulator data
                let acc = accumulators.get(&i)?;
                let len = (region.end - region.start + 1) as usize;
                
                // -------------------------------------------------------------
                // Step 1: Integrate difference arrays to get raw WPS values
                // This is O(len) per region, ~2000 iterations for standard windows
                // -------------------------------------------------------------
                let mut wps_long_raw = Vec::with_capacity(len);
                let mut wps_short_raw = Vec::with_capacity(len);
                let mut spanning_long_raw = Vec::with_capacity(len);
                let mut endpoints_long_raw = Vec::with_capacity(len);
                let mut spanning_short_raw = Vec::with_capacity(len);
                let mut endpoints_short_raw = Vec::with_capacity(len);
                
                let mut curr_wps_long = 0.0;
                let mut curr_wps_short = 0.0;
                let mut curr_span_long = 0.0;
                let mut curr_end_long = 0.0;
                let mut curr_span_short = 0.0;
                let mut curr_end_short = 0.0;
                let mut curr_cov = 0.0;
                
                for j in 0..len {
                    curr_wps_long += acc.wps_long[j];
                    curr_wps_short += acc.wps_short[j];
                    curr_span_long += acc.spanning_long[j];
                    curr_end_long += acc.endpoints_long[j];
                    curr_span_short += acc.spanning_short[j];
                    curr_end_short += acc.endpoints_short[j];
                    curr_cov += acc.cov_long[j] + acc.cov_short[j];
                    
                    wps_long_raw.push(curr_wps_long);
                    wps_short_raw.push(curr_wps_short);
                    spanning_long_raw.push(curr_span_long);
                    endpoints_long_raw.push(curr_end_long);
                    spanning_short_raw.push(curr_span_short);
                    endpoints_short_raw.push(curr_end_short);
                }
                
                // Skip empty regions (no coverage)
                if curr_cov == 0.0 {
                    return None;
                }
                
                // -------------------------------------------------------------
                // Step 2: Bin aggregation (2000bp → 200 bins of 10bp each)
                // -------------------------------------------------------------
                let mut binned_wps_nuc = vec![0.0f32; NUM_BINS];
                let mut binned_wps_tf = vec![0.0f32; NUM_BINS];
                let mut binned_mask = vec![1i8; NUM_BINS];
                let mut binned_prot_nuc = vec![0.0f32; NUM_BINS];
                let mut binned_prot_tf = vec![0.0f32; NUM_BINS];
                
                for bin_idx in 0..NUM_BINS {
                    let start_pos = bin_idx * BIN_SIZE;
                    let end_pos = (start_pos + BIN_SIZE).min(len);
                    
                    if start_pos >= len {
                        break;
                    }
                    
                    let mut sum_wps_nuc = 0.0;
                    let mut sum_wps_tf = 0.0;
                    let mut sum_span_nuc = 0.0;
                    let mut sum_end_nuc = 0.0;
                    let mut sum_span_tf = 0.0;
                    let mut sum_end_tf = 0.0;
                    let mut mask_ok = true;
                    let count = (end_pos - start_pos) as f32;
                    
                    for j in start_pos..end_pos {
                        sum_wps_nuc += wps_long_raw[j];
                        sum_wps_tf += wps_short_raw[j];
                        sum_span_nuc += spanning_long_raw[j];
                        sum_end_nuc += endpoints_long_raw[j];
                        sum_span_tf += spanning_short_raw[j];
                        sum_end_tf += endpoints_short_raw[j];
                        
                        // Check capture mask for panel edge detection
                        if let Some(ref mask) = bait_mask {
                            let chrom_id = 0u32; // Simplified - need proper chrom_id lookup
                            let pos = region.start + j as u64;
                            let (_, reliable) = mask.check_position(chrom_id, pos);
                            if !reliable {
                                mask_ok = false;
                            }
                        }
                    }
                    
                    binned_wps_nuc[bin_idx] = (sum_wps_nuc / count as f64) as f32;
                    binned_wps_tf[bin_idx] = (sum_wps_tf / count as f64) as f32;
                    binned_mask[bin_idx] = if mask_ok { 1 } else { 0 };
                    
                    let total_nuc = sum_span_nuc + sum_end_nuc;
                    let total_tf = sum_span_tf + sum_end_tf;
                    binned_prot_nuc[bin_idx] = if total_nuc > 0.0 {
                        (sum_span_nuc / total_nuc) as f32
                    } else {
                        0.0
                    };
                    binned_prot_tf[bin_idx] = if total_tf > 0.0 {
                        (sum_span_tf / total_tf) as f32
                    } else {
                        0.0
                    };
                }
                
                // -------------------------------------------------------------
                // Step 3: Strand reversal for minus strand regions
                // -------------------------------------------------------------
                if region.strand == "-" {
                    binned_wps_nuc.reverse();
                    binned_wps_tf.reverse();
                    binned_mask.reverse();
                    binned_prot_nuc.reverse();
                    binned_prot_tf.reverse();
                }
                
                // -------------------------------------------------------------
                // Step 4: Savitzky-Golay smoothing (window=11, order=3)
                // Matches scipy.signal.savgol_filter defaults
                // -------------------------------------------------------------
                use sci_rs::signal::filter::savgol_filter_dyn;
                const SAVGOL_WINDOW: usize = 11;
                const SAVGOL_POLYORDER: usize = 3;
                
                if binned_wps_nuc.len() >= SAVGOL_WINDOW {
                    binned_wps_nuc = savgol_filter_dyn(
                        binned_wps_nuc.iter().map(|x| *x as f64),
                        SAVGOL_WINDOW,
                        SAVGOL_POLYORDER,
                        None,
                        None,
                    )
                    .into_iter()
                    .map(|x| x as f32)
                    .collect();
                    
                    binned_wps_tf = savgol_filter_dyn(
                        binned_wps_tf.iter().map(|x| *x as f64),
                        SAVGOL_WINDOW,
                        SAVGOL_POLYORDER,
                        None,
                        None,
                    )
                    .into_iter()
                    .map(|x| x as f32)
                    .collect();
                }
                
                // Determine region type from ID prefix
                let region_type = if region.id.starts_with("TSS|") {
                    "TSS"
                } else if region.id.starts_with("CTCF|") {
                    "CTCF"
                } else {
                    "OTHER"
                };
                
                Some(RegionResult {
                    region_idx: i,
                    region_id: region.id.clone(),
                    chrom: region.chrom.clone(),
                    center: region.center as i32,
                    strand: region.strand.clone(),
                    region_type: region_type.to_string(),
                    binned_wps_nuc,
                    binned_wps_tf,
                    binned_mask,
                    binned_prot_nuc,
                    binned_prot_tf,
                    local_depth: (curr_cov / len as f64) as f32,
                })
            })
            .collect();
        
        let parallel_elapsed = parallel_start.elapsed();
        
        // =========================================================================
        // PARALLEL DIAGNOSTICS & ORDER VALIDATION
        // =========================================================================
        let valid_results: Vec<&RegionResult> = results.iter().filter_map(|r| r.as_ref()).collect();
        let n_computed = valid_results.len();
        
        // Validate output ordering (detect any par_iter ordering issues)
        let mut order_valid = true;
        let mut last_idx: Option<usize> = None;
        for result in &valid_results {
            if let Some(last) = last_idx {
                if result.region_idx <= last {
                    warn!("⚠️ ORDER VIOLATION: Region {} appears after {} (parallel ordering issue)", 
                        result.region_idx, last);
                    order_valid = false;
                }
            }
            last_idx = Some(result.region_idx);
        }
        
        if order_valid {
            debug!("Parallel region processing: {} regions computed in {:.2}s, order validated", 
                n_computed, parallel_elapsed.as_secs_f64());
        }
        
        // Log parallel summary
        let n_threads = rayon::current_num_threads();
        info!("WPS parallel binning: {} regions in {:.2}s ({} threads, {:.0} regions/sec)",
            n_computed, parallel_elapsed.as_secs_f64(), n_threads,
            n_computed as f64 / parallel_elapsed.as_secs_f64().max(0.001));
        
        // =========================================================================
        // SEQUENTIAL PARQUET WRITE
        // =========================================================================
        // Arrow builders are not thread-safe, so we write sequentially.
        // Order is preserved by par_iter().enumerate().collect() pattern.
        // =========================================================================
        let mut processed_count = 0usize;
        
        for result in valid_results {
            region_id_builder.append_value(&result.region_id);
            chrom_builder.append_value(&result.chrom);
            center_builder.append_value(result.center);
            strand_builder.append_value(&result.strand);
            region_type_builder.append_value(&result.region_type);
            
            // Append list values
            for v in &result.binned_wps_nuc {
                wps_nuc_builder.values().append_value(*v);
            }
            wps_nuc_builder.append(true);
            
            for v in &result.binned_wps_tf {
                wps_tf_builder.values().append_value(*v);
            }
            wps_tf_builder.append(true);
            
            for v in &result.binned_mask {
                mask_builder.values().append_value(*v);
            }
            mask_builder.append(true);
            
            for v in &result.binned_prot_nuc {
                prot_nuc_builder.values().append_value(*v);
            }
            prot_nuc_builder.append(true);
            
            for v in &result.binned_prot_tf {
                prot_tf_builder.values().append_value(*v);
            }
            prot_tf_builder.append(true);
            
            depth_builder.append_value(result.local_depth);
            processed_count += 1;
        }

        
        debug!("Loop complete. Processed {} regions with data", processed_count);
        
        // Build arrays
        let arrays: Vec<ArrayRef> = vec![
            StdArc::new(region_id_builder.finish()),
            StdArc::new(chrom_builder.finish()),
            StdArc::new(center_builder.finish()),
            StdArc::new(strand_builder.finish()),
            StdArc::new(region_type_builder.finish()),
            StdArc::new(wps_nuc_builder.finish()),
            StdArc::new(wps_tf_builder.finish()),
            StdArc::new(mask_builder.finish()),
            StdArc::new(prot_nuc_builder.finish()),
            StdArc::new(prot_tf_builder.finish()),
            StdArc::new(depth_builder.finish()),
        ];
        
        // Log array lengths for debugging
        debug!("Array lengths: region_id={}, chrom={}, center={}, strand={}, type={}, wps_nuc={}, wps_tf={}, mask={}, prot_nuc={}, prot_tf={}, depth={}",
            arrays[0].len(), arrays[1].len(), arrays[2].len(), arrays[3].len(), arrays[4].len(),
            arrays[5].len(), arrays[6].len(), arrays[7].len(), arrays[8].len(), arrays[9].len(), arrays[10].len());
        
        let batch = match RecordBatch::try_new(StdArc::new(schema), arrays) {
            Ok(b) => b,
            Err(e) => {
                error!("RecordBatch::try_new failed: {:?}", e);
                return Err(anyhow::anyhow!("RecordBatch creation failed: {}", e));
            }
        };
        
        debug!("Created batch with {} rows, writing to {:?}", batch.num_rows(), output_path);
        
        // Write Parquet
        let file = File::create(output_path)
            .context(format!("Creating output file: {:?}", output_path))?;
        let props = WriterProperties::builder().build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))
            .context("Creating ArrowWriter")?;
        writer.write(&batch)
            .context("Writing batch to Parquet")?;
        writer.close()
            .context("Closing Parquet writer")?;
        
        info!("WPS: Wrote {} regions to {:?}", batch.num_rows(), output_path);
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
            
            // Get weighted classification for this fragment
            // Nucleosome: 1.0 for [160,175], 0.5 for [120,159] or [176,180]
            // TF: 1.0 for [35,80]
            let nuc_len_weight = self.config.nuc_weight(fragment.length);
            let tf_len_weight = self.config.tf_weight(fragment.length);
            
            // Skip fragments that don't contribute to either stream
            if nuc_len_weight == 0.0 && tf_len_weight == 0.0 {
                return;
            }
            
            // GC correction weight
            let gc_weight = if let Some(ref factors) = self.factors {
                let gc_pct = (fragment.gc * 100.0).round() as u8;
                factors.get_factor(fragment.length, gc_pct)
            } else {
                1.0
            };
            
            // Track fragment counts
            if nuc_len_weight > 0.0 {
                self.nuc_fragments += 1;
            }
            if tf_len_weight > 0.0 {
                self.tf_fragments += 1;
            }
            
            // Now update accumulators
            for region_idx in matches {
                // Lazy allocation
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
                
                // Nucleosome stream (WPS-Nuc)
                if nuc_len_weight > 0.0 {
                    let weight = nuc_len_weight * gc_weight;
                    
                    // Coverage: [start, end)
                    RegionAccumulator::add_range(&mut acc.cov_long, f_start, f_end, weight);
                    
                    // WPS: protection window from config
                    let p = self.config.nuc_window;
                    RegionAccumulator::add_range(&mut acc.wps_long, f_start + p, f_end - p + 1, weight);
                    RegionAccumulator::add_range(&mut acc.wps_long, f_start - p, f_start + p, -weight);
                    RegionAccumulator::add_range(&mut acc.wps_long, f_end - p + 1, f_end + p + 1, -weight);
                    
                    // Protection fraction tracking
                    // Spanning: central window [start+p, end-p]
                    RegionAccumulator::add_range(&mut acc.spanning_long, f_start + p, f_end - p + 1, weight);
                    // Endpoints: flanking windows [start-p, start+p) and [end-p, end+p)
                    RegionAccumulator::add_range(&mut acc.endpoints_long, f_start - p, f_start + p, weight);
                    RegionAccumulator::add_range(&mut acc.endpoints_long, f_end - p + 1, f_end + p + 1, weight);
                }
                
                // TF stream (WPS-TF)
                if tf_len_weight > 0.0 {
                    let weight = tf_len_weight * gc_weight;
                    
                    RegionAccumulator::add_range(&mut acc.cov_short, f_start, f_end, weight);
                    
                    let p = self.config.tf_window;
                    RegionAccumulator::add_range(&mut acc.wps_short, f_start + p, f_end - p + 1, weight);
                    RegionAccumulator::add_range(&mut acc.wps_short, f_start - p, f_start + p, -weight);
                    RegionAccumulator::add_range(&mut acc.wps_short, f_end - p + 1, f_end + p + 1, -weight);
                    
                    // Protection fraction tracking for TF
                    RegionAccumulator::add_range(&mut acc.spanning_short, f_start + p, f_end - p + 1, weight);
                    RegionAccumulator::add_range(&mut acc.endpoints_short, f_start - p, f_start + p, weight);
                    RegionAccumulator::add_range(&mut acc.endpoints_short, f_end - p + 1, f_end + p + 1, weight);
                }
            }
        }
    }

    fn merge(&mut self, other: Self) {
        // Merge fragment counters
        self.nuc_fragments += other.nuc_fragments;
        self.tf_fragments += other.tf_fragments;
        
        // Merge accumulators
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


/// Alu region for background stacking with subfamily info
#[derive(Debug, Clone)]
pub struct AluRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: String,
    pub subfamily: String,  // AluY, AluS, AluJ, AluOther
}

/// Parse Alu BED7 file with subfamily column
/// Format: chrom, start, end, name, score, strand, subfamily
pub fn parse_alu_regions(bed_path: &Path) -> Result<Vec<AluRegion>> {
    // Use shared get_reader for consistent BGZF/gzip handling
    let reader = crate::bed::get_reader(bed_path)?;
    
    let mut regions = Vec::new();
    let valid_chroms: Vec<String> = (1..=22).map(|i| i.to_string())
        .chain(vec!["X".to_string(), "Y".to_string()]).collect();
    
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') { continue; }
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 7 { continue; }
        
        let chrom = fields[0].trim_start_matches("chr").to_string();
        if !valid_chroms.contains(&chrom) { continue; }
        
        let start: u64 = fields[1].parse().unwrap_or(0);
        let end: u64 = fields[2].parse().unwrap_or(0);
        let strand = fields[5].to_string();
        let subfamily = fields[6].to_string();
        
        regions.push(AluRegion { chrom, start, end, strand, subfamily });
    }
    
    info!("Parsed {} Alu regions from {:?}", regions.len(), bed_path);
    Ok(regions)
}

/// Background WPS Consumer with Hierarchical Stacking
/// 
/// Stacks Alu elements into multiple groups:
/// - Global_All: All Alus stacked (~770K elements)
/// - Family_AluY/AluS/AluJ: By evolutionary age
/// - Chr{N}_All: Per-chromosome stacking (24 groups)
#[derive(Clone)]
pub struct WpsBackgroundConsumer {
    trees: Arc<HashMap<u32, COITree<usize, u32>>>,  // region_idx only (Copy-friendly)
    regions: Arc<Vec<AluRegion>>,
    config: Arc<WpsConfig>,
    
    // Hierarchical accumulators: group_name -> (wps_nuc, wps_tf, count, nuc_frags, tf_frags)
    accumulators: HashMap<String, (Vec<f64>, Vec<f64>, u64, u64, u64)>,
}

impl WpsBackgroundConsumer {
    /// Length of the Alu body itself (~300bp).
    const ALU_BODY_BP: usize = 300;
    /// Flank stacked on each side of the Alu body.
    ///
    /// The nucleosome repeat length we are trying to measure is ~190bp. A
    /// 300bp profile spans only ~1.5 repeats, which is not enough to resolve
    /// that period at all: with 10bp bins the DFT period grid is
    /// `profile_len / i`, and for a 300bp profile the ONLY index landing in a
    /// nucleosomal search band is i = 2 (150bp). That made `nrl_bp`,
    /// `periodicity_score` and `adjusted_score` constants independent of the
    /// data. Stacking +/- 850bp of flank gives a 2000bp profile spanning
    /// ~10 nucleosome repeats.
    const FLANK_BP: usize = 850;
    /// Total stacked profile: 850 + 300 + 850 = 2000bp, Alu body centred.
    const PROFILE_LENGTH: usize = Self::FLANK_BP * 2 + Self::ALU_BODY_BP;
    
    pub fn new(regions: Vec<AluRegion>, chrom_map: &mut crate::bed::ChromosomeMap) -> Self {
        use coitrees::{COITree, IntervalNode};
        
        // Initialize accumulators for all groups
        let mut accumulators: HashMap<String, (Vec<f64>, Vec<f64>, u64, u64, u64)> = HashMap::new();
        
        // Global_All
        accumulators.insert("Global_All".to_string(), (
            vec![0.0; Self::PROFILE_LENGTH],
            vec![0.0; Self::PROFILE_LENGTH],
            0, 0, 0
        ));
        
        // Family groups
        for family in &["Family_AluY", "Family_AluS", "Family_AluJ", "Family_AluOther"] {
            accumulators.insert(family.to_string(), (
                vec![0.0; Self::PROFILE_LENGTH],
                vec![0.0; Self::PROFILE_LENGTH],
                0, 0, 0
            ));
        }
        
        // Chromosome groups
        for chr in (1..=22).map(|i| i.to_string()).chain(vec!["X".to_string(), "Y".to_string()]) {
            accumulators.insert(format!("Chr{}_All", chr), (
                vec![0.0; Self::PROFILE_LENGTH],
                vec![0.0; Self::PROFILE_LENGTH],
                0, 0, 0
            ));
        }
        
        // Build interval tree with just region index (Copy-friendly)
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        
        for (i, region) in regions.iter().enumerate() {
            let chrom_id = chrom_map.get_id(&region.chrom);
            // Widen by FLANK_BP so fragments sitting in the flanking window
            // (not just inside the Alu body) are retrieved and stacked.
            let lo = (region.start as i64 - Self::FLANK_BP as i64).max(0) as i32;
            let hi = (region.end as i64 + Self::FLANK_BP as i64) as i32;
            let node = IntervalNode::new(lo, hi, i);
            nodes_by_chrom.entry(chrom_id).or_default().push(node);
        }
        
        let mut trees = HashMap::new();
        for (chrom_id, nodes) in nodes_by_chrom {
            trees.insert(chrom_id, COITree::new(&nodes));
        }
        
        info!("WPS Background: Loaded {} Alu regions for hierarchical stacking ({} groups)", 
              regions.len(), accumulators.len());
        
        Self {
            trees: Arc::new(trees),
            regions: Arc::new(regions),
            config: Arc::new(WpsConfig::default()),
            accumulators,
        }
    }
    
    /// Helper to add signal to a group accumulator
    fn add_to_group(&mut self, group_name: &str, pos: usize, nuc_weight: f64, tf_weight: f64, is_valid_frag: bool) {
        if let Some((wps_nuc, wps_tf, count, nuc_frags, tf_frags)) = self.accumulators.get_mut(group_name) {
            if pos < Self::PROFILE_LENGTH {
                wps_nuc[pos] += nuc_weight;
                wps_tf[pos] += tf_weight;
            }
            if is_valid_frag {
                *count += 1;
                if nuc_weight > 0.0 { *nuc_frags += 1; }
                if tf_weight > 0.0 { *tf_frags += 1; }
            }
        }
    }
    
    /// Write hierarchical background metrics to Parquet (~30 rows)
    pub fn write_parquet(&self, output_path: &Path) -> Result<()> {
        use arrow::array::{
            ArrayRef, BooleanBuilder, Float32Builder, Int64Builder, ListBuilder, StringBuilder,
        };
        use arrow::datatypes::{DataType, Field, Schema};
        use arrow::record_batch::RecordBatch;
        use parquet::arrow::ArrowWriter;
        use parquet::file::properties::WriterProperties;
        use std::sync::Arc as StdArc;
        
        const BIN_SIZE: usize = 10;
        // 2000bp profile / 10bp bins = 200 bins.
        const NUM_BINS: usize = WpsBackgroundConsumer::PROFILE_LENGTH / BIN_SIZE;
        
        info!("WPS Background: Writing hierarchical Parquet ({} groups)", self.accumulators.len());
        
        // Schema - inner List items must be nullable:true to match ListBuilder defaults
        let schema = Schema::new(vec![
            Field::new("group_id", DataType::Utf8, false),
            Field::new("stacked_wps_nuc", DataType::List(StdArc::new(Field::new("item", DataType::Float32, true))), false),
            Field::new("stacked_wps_tf", DataType::List(StdArc::new(Field::new("item", DataType::Float32, true))), false),
            Field::new("alu_count", DataType::Int64, false),
            Field::new("mean_wps_nuc", DataType::Float32, false),
            Field::new("mean_wps_tf", DataType::Float32, false),
            Field::new("nrl_bp", DataType::Float32, false),               // Nucleosome Repeat Length in bp
            Field::new("nrl_deviation_bp", DataType::Float32, false),     // Deviation from expected (190bp)
            Field::new("periodicity_score", DataType::Float32, false),    // Raw quality score 0-1 (SNR-based)
            Field::new("adjusted_score", DataType::Float32, false),       // Deviation-penalized score
            // True when nrl_bp is the search-band edge rather than a peak:
            // a right-censored estimate. Consumers must be able to tell
            // "at least this" from "exactly this".
            Field::new("nrl_at_band_limit", DataType::Boolean, false),
            Field::new("fragment_ratio", DataType::Float32, false),
        ]);
        
        // Builders
        let mut group_builder = StringBuilder::new();
        let mut wps_nuc_builder = ListBuilder::new(Float32Builder::new());
        let mut wps_tf_builder = ListBuilder::new(Float32Builder::new());
        let mut count_builder = Int64Builder::new();
        let mut mean_nuc_builder = Float32Builder::new();
        let mut mean_tf_builder = Float32Builder::new();
        let mut nrl_builder = Float32Builder::new();             // NRL in bp
        let mut nrl_deviation_builder = Float32Builder::new();    // Deviation from 190bp
        let mut periodicity_builder = Float32Builder::new();      // Raw quality score 0-1
        let mut adjusted_score_builder = Float32Builder::new();   // Deviation-penalized score
        let mut at_band_limit_builder = BooleanBuilder::new();    // nrl_bp is a censoring bound
        let mut frag_ratio_builder = Float32Builder::new();
        
        // Sort groups for consistent output order
        let mut groups: Vec<_> = self.accumulators.keys().cloned().collect();
        groups.sort();
        
        for group_name in &groups {
            let (wps_nuc, wps_tf, count, nuc_frags, tf_frags) = &self.accumulators[group_name];
            
            // Skip empty groups
            if *count == 0 { continue; }
            
            // Bin the profile
            let mut binned_nuc = vec![0.0f32; NUM_BINS];
            let mut binned_tf = vec![0.0f32; NUM_BINS];
            
            for bin_idx in 0..NUM_BINS {
                let start = bin_idx * BIN_SIZE;
                let end = (start + BIN_SIZE).min(Self::PROFILE_LENGTH);
                
                let mut sum_nuc = 0.0;
                let mut sum_tf = 0.0;
                for j in start..end {
                    sum_nuc += wps_nuc[j];
                    sum_tf += wps_tf[j];
                }
                
                let cnt = (end - start) as f64;
                binned_nuc[bin_idx] = (sum_nuc / cnt) as f32;
                binned_tf[bin_idx] = (sum_tf / cnt) as f32;
            }
            
            // Apply Savitzky-Golay smoothing to stacked profiles
            // Use smaller window for 30-bin background (window=7, order=3)
            use sci_rs::signal::filter::savgol_filter_dyn;
            const BG_SAVGOL_WINDOW: usize = 7;
            const BG_SAVGOL_POLYORDER: usize = 3;
            
            if binned_nuc.len() >= BG_SAVGOL_WINDOW {
                binned_nuc = savgol_filter_dyn(
                    binned_nuc.iter().map(|x| *x as f64),
                    BG_SAVGOL_WINDOW, BG_SAVGOL_POLYORDER, None, None
                ).into_iter().map(|x| x as f32).collect();
                
                binned_tf = savgol_filter_dyn(
                    binned_tf.iter().map(|x| *x as f64),
                    BG_SAVGOL_WINDOW, BG_SAVGOL_POLYORDER, None, None
                ).into_iter().map(|x| x as f32).collect();
            }
            
            // Metrics
            let mean_nuc = wps_nuc.iter().sum::<f64>() / Self::PROFILE_LENGTH as f64;
            let mean_tf = wps_tf.iter().sum::<f64>() / Self::PROFILE_LENGTH as f64;
            let frag_ratio = if *nuc_frags > 0 { *tf_frags as f32 / *nuc_frags as f32 } else { 0.0 };
            let (nrl_bp, periodicity_score, nrl_at_band_limit) =
                Self::estimate_periodicity(&binned_nuc);
            
            // Calculate deviation from expected NRL (190bp) and penalized score
            const EXPECTED_NRL_BP: f32 = 190.0;
            // Documented as 50bp in docs/features/core/wps.md; the code used
            // 20bp, which combined with the pinned 150bp NRL forced
            // adjusted_score to 0.0 for every sample.
            const TOLERANCE_BP: f32 = 50.0;
            let nrl_deviation = (nrl_bp - EXPECTED_NRL_BP).abs();
            let deviation_penalty = (1.0 - nrl_deviation / TOLERANCE_BP).max(0.0);
            let adjusted_score = periodicity_score * deviation_penalty;
            
            // Append row
            group_builder.append_value(group_name);
            for v in &binned_nuc { wps_nuc_builder.values().append_value(*v); }
            wps_nuc_builder.append(true);
            for v in &binned_tf { wps_tf_builder.values().append_value(*v); }
            wps_tf_builder.append(true);
            count_builder.append_value(*count as i64);
            mean_nuc_builder.append_value(mean_nuc as f32);
            mean_tf_builder.append_value(mean_tf as f32);
            nrl_builder.append_value(nrl_bp);
            nrl_deviation_builder.append_value(nrl_deviation);
            periodicity_builder.append_value(periodicity_score);
            adjusted_score_builder.append_value(adjusted_score);
            at_band_limit_builder.append_value(nrl_at_band_limit);
            frag_ratio_builder.append_value(frag_ratio);
        }
        
        let arrays: Vec<ArrayRef> = vec![
            StdArc::new(group_builder.finish()),
            StdArc::new(wps_nuc_builder.finish()),
            StdArc::new(wps_tf_builder.finish()),
            StdArc::new(count_builder.finish()),
            StdArc::new(mean_nuc_builder.finish()),
            StdArc::new(mean_tf_builder.finish()),
            StdArc::new(nrl_builder.finish()),
            StdArc::new(nrl_deviation_builder.finish()),
            StdArc::new(periodicity_builder.finish()),
            StdArc::new(adjusted_score_builder.finish()),
            StdArc::new(at_band_limit_builder.finish()),
            StdArc::new(frag_ratio_builder.finish()),
        ];
        
        let batch = RecordBatch::try_new(StdArc::new(schema), arrays)?;
        
        let file = File::create(output_path)?;
        let props = WriterProperties::builder().build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))?;
        writer.write(&batch)?;
        writer.close()?;
        
        if batch.num_rows() == 0 {
            log::warn!("WPS Background: No hierarchical groups written - fragments may not overlap Alu regions");
        }
        info!("WPS Background: Wrote {} hierarchical groups to {:?}", batch.num_rows(), output_path);
        Ok(())
    }
    
    /// Extract periodicity from WPS profile using FFT
    /// 
    /// Computes the Nucleosome Repeat Length (NRL) and quality score from a WPS profile.
    /// 
    /// ## Algorithm
    /// 1. **Detrend**: Subtract mean to remove DC component
    /// 2. **Z-score normalize**: Standardize to mean=0, std=1
    /// 3. **Hann window**: Apply to reduce spectral leakage
    /// 4. **FFT**: Compute real-to-complex FFT via realfft
    /// 5. **Peak finding**: Find dominant peak in 140-250bp period range
    /// 6. **Quality score**: SNR-based score capped at 1.0 (SNR > 3 = clear periodicity)
    /// 
    /// ## Returns
    ///
    /// `(nrl_bp, quality_score, at_band_limit)`:
    ///
    /// - `nrl_bp`: Nucleosome Repeat Length in bp (expected ~190bp for healthy cfDNA)
    /// - `quality_score`: SNR-based quality 0-1 (higher = clearer periodicity signal)
    /// - `at_band_limit`: the peak sat on the edge of the search band, so
    ///   `nrl_bp` is that bound rather than a measurement -- a right-censored
    ///   estimate. Reported separately because the value alone cannot
    ///   distinguish "at least this" from "exactly this", and on real plasma
    ///   this is 14-43% of Alu groups depending on assay.
    fn estimate_periodicity(profile: &[f32]) -> (f32, f32, bool) {
        use realfft::RealFftPlanner;
        use std::f64::consts::PI;
        
        let n = profile.len();
        if n < 10 { return (0.0, 0.0, false); }
        
        // Convert to f64 for precision
        let mut data: Vec<f64> = profile.iter().map(|x| *x as f64).collect();
        
        // 1. Detrend (remove linear trend)
        let n_f = n as f64;
        let x_mean = (n_f - 1.0) / 2.0;
        let y_mean = data.iter().sum::<f64>() / n_f;
        
        let mut slope_numer = 0.0;
        let mut slope_denom = 0.0;
        for (i, _) in data.iter().enumerate().take(n) {
            let xi = i as f64 - x_mean;
            let yi = data[i] - y_mean;
            slope_numer += xi * yi;
            slope_denom += xi * xi;
        }
        let slope = if slope_denom.abs() > 1e-9 { slope_numer / slope_denom } else { 0.0 };
        let intercept = y_mean - slope * x_mean;
        
        for (i, val) in data.iter_mut().enumerate().take(n) {
            *val -= slope * i as f64 + intercept;
        }
        
        // 2. Z-score normalize
        let mean: f64 = data.iter().sum::<f64>() / n_f;
        let variance: f64 = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n_f;
        let std = variance.sqrt();
        
        if std < 1e-9 {
            return (0.0, 0.0, false); // Flat profile
        }
        
        for val in data.iter_mut() {
            *val /= std;
        }
        
        // 3. Apply Hann window
        for (i, val) in data.iter_mut().enumerate().take(n) {
            let window = 0.5 * (1.0 - (2.0 * PI * i as f64 / (n_f - 1.0)).cos());
            *val *= window;
        }
        
        // 4. Compute FFT using realfft, zero-padded for frequency resolution
        const PAD_FACTOR: usize = 8;
        let n_fft = n * PAD_FACTOR;

        let mut planner = RealFftPlanner::<f64>::new();
        let fft = planner.plan_fft_forward(n_fft);

        // Hann-windowed signal in the first n samples, zeros thereafter.
        let mut indata = vec![0.0f64; n_fft];
        indata[..n].copy_from_slice(&data);
        let mut spectrum = fft.make_output_vec();

        if fft.process(&mut indata, &mut spectrum).is_err() {
            return (0.0, 0.0, false);
        }
        
        // 5. Find peak in the nucleosomal period band.
        //
        // freq[i] = i / (n_fft * bin_size)  =>  period[i] = n_fft * bin_size / i
        //
        // The native grid is coarse, so we zero-pad by PAD_FACTOR before the
        // transform. Zero-padding interpolates the spectrum: it does not add
        // information, but it does let the peak be located between native
        // bins, and it gives the SNR denominator more than one sample to
        // average over. Without it a 2000bp profile offers only 5 candidate
        // periods in the band (200.0, 181.8, 166.7, 153.8, 142.9bp).
        let bin_size_bp = 10.0;
        let min_period_bp = 140.0;
        let max_period_bp = 250.0;
        
        let amplitudes: Vec<f64> = spectrum.iter().map(|c| c.norm()).collect();
        
        let mut peak_amplitude = 0.0;
        let mut peak_idx = 0;
        let mut sum_amplitude = 0.0;
        let mut count_in_range = 0;
        // Track the eligible index range so a peak sitting on the band edge can
        // be recognised. Comparing the resulting period against 250.0 would be
        // a float comparison against a value the grid only approximates; the
        // index is exact.
        let mut first_eligible = usize::MAX;
        let mut last_eligible = 0usize;

        for (i, &amp) in amplitudes.iter().enumerate().skip(1) {
            let period_bp = (n_fft as f64 * bin_size_bp) / i as f64;

            if period_bp >= min_period_bp && period_bp <= max_period_bp {
                sum_amplitude += amp;
                count_in_range += 1;
                first_eligible = first_eligible.min(i);
                last_eligible = last_eligible.max(i);

                if amp > peak_amplitude {
                    peak_amplitude = amp;
                    peak_idx = i;
                }
            }
        }

        if count_in_range == 0 || peak_idx == 0 {
            return (0.0, 0.0, false);
        }

        // The argmax landing on either end of the band means the true peak is
        // outside it, or there is no peak and the spectrum's slope decided the
        // winner. Either way the reported period is the edge of the search
        // window rather than a measurement -- right-censored, not missing.
        // Measured on real plasma this is 14-43% of Alu groups depending on
        // assay, and those groups are indistinguishable from interior ones by
        // periodicity score or fragment support, so consumers cannot infer it.
        let at_band_limit = peak_idx == first_eligible || peak_idx == last_eligible;
        
        // 6. Calculate SNR and quality score
        let background = sum_amplitude / count_in_range as f64;
        let snr = if background > 0.0 { peak_amplitude / background } else { 0.0 };
        
        // Calculate NRL (nucleosome repeat length) from peak index
        let nrl_bp = (n_fft as f64 * bin_size_bp) / peak_idx as f64;
        
        // Quality score: SNR-based, capped at 1.0 (SNR > 3 indicates clear periodicity)
        let quality_score = (snr / 3.0).min(1.0) as f32;
        
        (nrl_bp as f32, quality_score, at_band_limit)
    }
}

impl crate::engine::FragmentConsumer for WpsBackgroundConsumer {
    fn name(&self) -> &str {
        "WPS_Background"
    }

    fn consume(&mut self, fragment: &crate::bed::Fragment) {
        if let Some(tree) = self.trees.get(&fragment.chrom_id) {
            let start = fragment.start as u32;
            let end = fragment.end as u32;
            let end_closed = if end > start { end - 1 } else { start };
            
            let frag_len = fragment.end - fragment.start;
            let nuc_weight = self.config.nuc_weight(frag_len);
            let tf_weight = self.config.tf_weight(frag_len);
            
            if nuc_weight == 0.0 && tf_weight == 0.0 { return; }
            
            // Collect overlapping region indices first (avoids borrow conflict)
            let mut overlaps: Vec<usize> = Vec::new();
            // coitrees 0.4.0: node.metadata is &usize on macOS but usize on CI Linux.
            // .clone() handles both; suppress clippy::clone_on_copy for the usize case.
            #[allow(clippy::clone_on_copy)]
            tree.query(start as i32, end_closed as i32, |node| {
                overlaps.push(node.metadata.clone());
            });
            
            // Process each overlapping region
            for region_idx in overlaps {
                let region = &self.regions[region_idx];
                let is_minus = region.strand == "-";
                let subfamily = region.subfamily.clone();
                let chrom = region.chrom.clone();
                let r_start = region.start;
                
                // Fragment position relative to the Alu, shifted into the
                // flanked profile frame: index FLANK_BP == Alu start.
                let f_start = (fragment.start + 1) as i64 - r_start as i64
                    + Self::FLANK_BP as i64;
                let f_end = fragment.end as i64 - r_start as i64
                    + Self::FLANK_BP as i64;
                
                let p_nuc = self.config.nuc_window;
                
                // Spanning positions
                for pos in (f_start + p_nuc).max(0)..=(f_end - p_nuc).min(Self::PROFILE_LENGTH as i64 - 1) {
                    let idx = if is_minus {
                        Self::PROFILE_LENGTH as i64 - 1 - pos
                    } else {
                        pos
                    } as usize;
                    
                    // Add to Global_All
                    self.add_to_group("Global_All", idx, nuc_weight, tf_weight, false);
                    
                    // Add to Family group
                    let family_key = format!("Family_{}", subfamily);
                    self.add_to_group(&family_key, idx, nuc_weight, tf_weight, false);
                    
                    // Add to Chromosome group
                    let chr_key = format!("Chr{}_All", chrom);
                    self.add_to_group(&chr_key, idx, nuc_weight, tf_weight, false);
                }
                
                // Count once per region overlap
                if let Some((_, _, count, nuc_frags, tf_frags)) = self.accumulators.get_mut("Global_All") {
                    *count += 1;
                    if nuc_weight > 0.0 { *nuc_frags += 1; }
                    if tf_weight > 0.0 { *tf_frags += 1; }
                }
                let family_key = format!("Family_{}", subfamily);
                if let Some((_, _, count, nuc_frags, tf_frags)) = self.accumulators.get_mut(&family_key) {
                    *count += 1;
                    if nuc_weight > 0.0 { *nuc_frags += 1; }
                    if tf_weight > 0.0 { *tf_frags += 1; }
                }
                let chr_key = format!("Chr{}_All", chrom);
                if let Some((_, _, count, nuc_frags, tf_frags)) = self.accumulators.get_mut(&chr_key) {
                    *count += 1;
                    if nuc_weight > 0.0 { *nuc_frags += 1; }
                    if tf_weight > 0.0 { *tf_frags += 1; }
                }
            }
        }
    }

    fn merge(&mut self, other: Self) {
        for (group_name, (other_nuc, other_tf, other_count, other_nuc_frags, other_tf_frags)) in other.accumulators {
            if let Some((self_nuc, self_tf, self_count, self_nuc_frags, self_tf_frags)) = self.accumulators.get_mut(&group_name) {
                for i in 0..Self::PROFILE_LENGTH {
                    self_nuc[i] += other_nuc[i];
                    self_tf[i] += other_tf[i];
                }
                *self_count += other_count;
                *self_nuc_frags += other_nuc_frags;
                *self_tf_frags += other_tf_frags;
            }
        }
    }
}

// NOTE: `apply_pon_zscore` was removed here, along with the parquet
// read/write imports that existed only to serve it.
//
// It emitted a single scalar z per anchor from `wps_long_std`/`wps_short_std`
// -- the v1.0 scalar baseline field names, while every shipped PON is v2.0
// vector format -- and pushed 0.0 wherever the sigma was not positive. It was
// also unreachable: `wps_processor.py` gated the call on `pon._source_path`,
// an attribute nothing in the codebase ever set.
//
// So it was dead, schema-obsolete, and fabricating. Scoring now lives in
// `core/wps_pon.py`, which emits per-position z vectors plus three derived
// shape quantities -- because measurement showed a scalar reduction of the
// profile is the wrong statistic: adjacent positions have lag-1
// autocorrelation 0.986.

// =============================================================================
// TESTS
// =============================================================================

#[cfg(test)]
mod background_periodicity_tests {
    use super::*;

    const BIN_SIZE_BP: f64 = 10.0;

    /// Build a binned background profile carrying a known period.
    fn synth_profile(period_bp: f64, n_bins: usize, amplitude: f64, noise: f64) -> Vec<f32> {
        (0..n_bins)
            .map(|i| {
                let x = i as f64 * BIN_SIZE_BP;
                // deterministic pseudo-noise so the test is reproducible
                let jitter = ((i * 2654435761) % 1000) as f64 / 1000.0 - 0.5;
                (amplitude * (2.0 * std::f64::consts::PI * x / period_bp).cos()
                    + noise * jitter) as f32
            })
            .collect()
    }

    #[test]
    fn profile_length_spans_multiple_nucleosome_repeats() {
        // A 300bp profile cannot resolve a ~190bp period; 2000bp spans ~10.
        assert_eq!(WpsBackgroundConsumer::PROFILE_LENGTH, 2000);
        let repeats = WpsBackgroundConsumer::PROFILE_LENGTH as f64 / 190.0;
        assert!(repeats >= 5.0, "profile spans only {repeats:.1} repeats");
    }

    #[test]
    fn recovers_a_190bp_nucleosome_repeat_length() {
        let n_bins = WpsBackgroundConsumer::PROFILE_LENGTH / 10;
        let profile = synth_profile(190.0, n_bins, 1.0, 0.05);

        let (nrl_bp, score, _) = WpsBackgroundConsumer::estimate_periodicity(&profile);

        assert!(
            (nrl_bp - 190.0).abs() < 10.0,
            "expected NRL near 190bp, got {nrl_bp}"
        );
        assert!(score > 0.0, "expected non-zero periodicity score, got {score}");
    }

    #[test]
    fn distinguishes_different_periods() {
        // Regression: NRL used to be pinned at 150.0 for every input, so two
        // genuinely different profiles produced identical output.
        let n_bins = WpsBackgroundConsumer::PROFILE_LENGTH / 10;

        let (nrl_a, _, _) = WpsBackgroundConsumer::estimate_periodicity(
            &synth_profile(190.0, n_bins, 1.0, 0.05),
        );
        let (nrl_b, _, _) = WpsBackgroundConsumer::estimate_periodicity(
            &synth_profile(160.0, n_bins, 1.0, 0.05),
        );

        assert!(
            (nrl_a - nrl_b).abs() > 10.0,
            "NRL must track the input period, got {nrl_a} and {nrl_b}"
        );
        assert!(nrl_a > nrl_b, "190bp input should yield a longer NRL than 160bp");
    }

    #[test]
    fn a_period_inside_the_band_is_not_flagged_as_censored() {
        let n_bins = WpsBackgroundConsumer::PROFILE_LENGTH / 10;
        let (nrl_bp, _, at_limit) =
            WpsBackgroundConsumer::estimate_periodicity(&synth_profile(190.0, n_bins, 1.0, 0.05));

        assert!(
            !at_limit,
            "a clean 190bp signal sits well inside 140-250 and must not be \
             reported as censored (got nrl_bp={nrl_bp})"
        );
    }

    #[test]
    fn a_profile_with_no_periodicity_is_flagged_as_censored() {
        // The point of the flag, and the case that actually occurs. A profile
        // that is a monotonic trend has no nucleosomal peak at all, so the
        // spectrum's low-frequency slope wins the argmax and it lands on the
        // band ceiling. nrl_bp then reads 250.0, which is the edge of the
        // search window rather than a measurement -- indistinguishable, by
        // value alone, from a genuine 250bp repeat.
        //
        // Note this is NOT triggered by a long period: 400bp and 2000bp
        // sinusoids both resolve to 225-235bp inside the band. Only the
        // absence of periodicity pins the edge, which is why the flag means
        // "no peak found" rather than "period too long".
        let n_bins = WpsBackgroundConsumer::PROFILE_LENGTH / 10;
        let ramp: Vec<f32> = (0..n_bins)
            .map(|i| i as f32 / n_bins as f32)
            .collect();

        let (nrl_bp, _, at_limit) = WpsBackgroundConsumer::estimate_periodicity(&ramp);

        assert!(
            at_limit,
            "a trend with no periodicity must be flagged censored, got \
             nrl_bp={nrl_bp} with at_band_limit=false"
        );
    }

    #[test]
    fn a_long_period_still_resolves_inside_the_band() {
        // Guards the boundary of the previous test: the flag must not fire
        // merely because the input period is long, or it would discard real
        // measurements.
        let n_bins = WpsBackgroundConsumer::PROFILE_LENGTH / 10;
        let (nrl_bp, _, at_limit) =
            WpsBackgroundConsumer::estimate_periodicity(&synth_profile(400.0, n_bins, 1.0, 0.05));

        assert!(
            !at_limit,
            "a 400bp input resolves to a peak inside the band and must not be \
             flagged censored (got nrl_bp={nrl_bp})"
        );
    }

    #[test]
    fn adjusted_score_is_not_forced_to_zero_at_expected_nrl() {
        // With NRL == 190 the deviation penalty must be 1.0 (no penalty).
        const EXPECTED_NRL_BP: f32 = 190.0;
        const TOLERANCE_BP: f32 = 50.0;

        let n_bins = WpsBackgroundConsumer::PROFILE_LENGTH / 10;
        let (nrl_bp, score, _) = WpsBackgroundConsumer::estimate_periodicity(
            &synth_profile(190.0, n_bins, 1.0, 0.05),
        );

        let penalty = (1.0 - (nrl_bp - EXPECTED_NRL_BP).abs() / TOLERANCE_BP).max(0.0);
        let adjusted = score * penalty;

        assert!(penalty > 0.5, "penalty collapsed at the expected NRL: {penalty}");
        assert!(adjusted > 0.0, "adjusted_score must be non-zero for healthy-like input");
    }

    #[test]
    fn flat_profile_yields_no_periodicity() {
        let profile = vec![0.0f32; WpsBackgroundConsumer::PROFILE_LENGTH / 10];
        let (nrl_bp, score, _) = WpsBackgroundConsumer::estimate_periodicity(&profile);
        assert_eq!((nrl_bp, score), (0.0, 0.0));
    }
}

// ===========================================================================
// PON scoring for the WPS foreground
// ===========================================================================
//
// Ported from `core/wps_pon.py`, which held the only PON z-score still
// computed in Python. `architecture.md` puts "PON z-score / log-ratio",
// "loops over >1000 rows" and "row-level computation" on the Rust side, and
// this was all three: 89,034 anchors, a +/-30 lag search per anchor, about
// 5.4M correlations, ~13 minutes per sample.
//
// It was written in Python first for good reason -- three of its decisions
// changed under measurement, and would not have survived being written here
// first. `log1p` amplitude, because raw range correlates +0.512 with depth.
// Fisher-z on the correlation, because bounded r left 302 of 400 anchors
// unable to reach +2. And no z-score at all on the phase shift, because its
// intraclass correlation came out at 0.479. The Python version stays on as
// the equivalence oracle in tests/unit/test_rust_python_equivalence.py.
//
// This is a bug-for-bug port. Where the Python is arguably wrong -- ties in
// the lag search resolving to the most negative lag -- the behaviour is
// reproduced and commented, not corrected, because an oracle you have
// deliberately diverged from cannot tell a porting error from an intended
// change. Measured incidence of those ties: 0 in 317 anchors.

/// How far the phase search looks, in positions. Mirrors `PHASE_MAX_LAG`.
const PHASE_MAX_LAG: i32 = 30;

/// Peak-to-trough range on a log scale, or NaN from fewer than two finite points.
fn wps_log_amplitude(profile: &[f64]) -> f64 {
    let finite: Vec<f64> = profile.iter().copied().filter(|v| v.is_finite()).collect();
    if finite.len() < 2 {
        return f64::NAN;
    }
    let max = finite.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min = finite.iter().cloned().fold(f64::INFINITY, f64::min);
    (max - min).ln_1p()
}

/// Pearson r over positions finite in both, or NaN.
///
/// The `1e-12` floor matches `np.std(a) < 1e-12` in the Python: a constant
/// window has no correlation to report, and dividing by its spread would
/// manufacture one.
fn wps_shape_correlation(profile: &[f64], baseline: &[f64]) -> f64 {
    let pairs: Vec<(f64, f64)> = profile
        .iter()
        .zip(baseline.iter())
        .filter(|(a, b)| a.is_finite() && b.is_finite())
        .map(|(a, b)| (*a, *b))
        .collect();
    if pairs.len() < 3 {
        return f64::NAN;
    }
    let n = pairs.len() as f64;
    let mean_a = pairs.iter().map(|(a, _)| a).sum::<f64>() / n;
    let mean_b = pairs.iter().map(|(_, b)| b).sum::<f64>() / n;
    // Population std, as numpy's default -- ddof=0. The correlation is
    // scale-free so the choice cancels, but it is matched for exactness.
    let var_a = pairs.iter().map(|(a, _)| (a - mean_a).powi(2)).sum::<f64>() / n;
    let var_b = pairs.iter().map(|(_, b)| (b - mean_b).powi(2)).sum::<f64>() / n;
    if var_a.sqrt() < 1e-12 || var_b.sqrt() < 1e-12 {
        return f64::NAN;
    }
    let cov = pairs
        .iter()
        .map(|(a, b)| (a - mean_a) * (b - mean_b))
        .sum::<f64>()
        / n;
    cov / (var_a.sqrt() * var_b.sqrt())
}

/// `arctanh(r)`, clipped as the Python clips, so a perfect correlation is finite.
fn wps_fisher_z(r: f64) -> f64 {
    if !r.is_finite() {
        return f64::NAN;
    }
    r.clamp(-0.999_999, 0.999_999).atanh()
}

/// Displacement of the sample against the baseline, and whether the search
/// ended on its own window edge.
///
/// The edge flag is `nrl_at_band_limit` one level down: a search that stops at
/// its boundary has found the boundary, not a measurement.
fn wps_phase_shift(profile: &[f64], baseline: &[f64], max_lag: i32) -> (f64, bool) {
    let n = profile.len().min(baseline.len());
    if (n as i32) < 4 * max_lag {
        return (f64::NAN, false);
    }
    let lo = max_lag as usize;
    let core = &baseline[lo..n - lo];
    let mut best_r = f64::NEG_INFINITY;
    let mut best_lag = 0i32;
    for k in -max_lag..=max_lag {
        let start = (lo as i32 + k) as usize;
        let window = &profile[start..start + core.len()];
        let r = wps_shape_correlation(window, core);
        // Strictly greater, so a tie keeps the earliest lag scanned -- the most
        // negative. Reproduced from the Python deliberately; see the header.
        if r.is_finite() && r > best_r {
            best_r = r;
            best_lag = k;
        }
    }
    if !best_r.is_finite() {
        return (f64::NAN, false);
    }
    (best_lag as f64, best_lag.abs() == max_lag)
}

/// `(x - mean) / std` per position, NaN where sigma is unusable.
///
/// Never a substituted 1.0: the builder writes NaN for positions whose spread
/// it could not measure, and a default here would undo that at the read side.
fn wps_z_vector(profile: &[f64], mean: &[f64], std: &[f64]) -> Vec<f64> {
    let n = profile.len().min(mean.len()).min(std.len());
    (0..n)
        .map(|i| {
            // Sigma alone, as `compute_z_vector` gates in model.py. A non-finite
            // profile or mean needs no test of its own: the subtraction carries
            // the NaN through, which is what the Python relies on too.
            if std[i].is_finite() && std[i] > 0.0 {
                (profile[i] - mean[i]) / std[i]
            } else {
                f64::NAN
            }
        })
        .collect()
}

/// `zscore_or_nan` for the derived scalar statistics.
fn wps_scalar_z(observed: f64, mean: f64, std: f64) -> f64 {
    if !std.is_finite() || std <= 0.0 || !mean.is_finite() || !observed.is_finite() {
        return f64::NAN;
    }
    (observed - mean) / std
}

#[cfg(test)]
mod pon_scoring_tests {
    use super::*;

    #[test]
    fn amplitude_needs_two_finite_points() {
        assert!(wps_log_amplitude(&[]).is_nan());
        assert!(wps_log_amplitude(&[1.0]).is_nan());
        assert!(wps_log_amplitude(&[f64::NAN, 2.0]).is_nan());
        // ln(1 + (5 - 1))
        assert!((wps_log_amplitude(&[1.0, 5.0]) - 4.0f64.ln_1p()).abs() < 1e-12);
    }

    #[test]
    fn correlation_refuses_a_constant_window() {
        let flat = [2.0, 2.0, 2.0, 2.0];
        let ramp = [1.0, 2.0, 3.0, 4.0];
        assert!(wps_shape_correlation(&flat, &ramp).is_nan(), "no spread to correlate");
        assert!((wps_shape_correlation(&ramp, &ramp) - 1.0).abs() < 1e-12);
        let down = [4.0, 3.0, 2.0, 1.0];
        assert!((wps_shape_correlation(&ramp, &down) + 1.0).abs() < 1e-12);
    }

    #[test]
    fn correlation_needs_three_shared_finite_points() {
        let a = [1.0, 2.0, f64::NAN, f64::NAN];
        let b = [1.0, 2.0, 3.0, 4.0];
        assert!(wps_shape_correlation(&a, &b).is_nan());
    }

    #[test]
    fn fisher_z_is_finite_at_a_perfect_correlation() {
        assert!(wps_fisher_z(1.0).is_finite(), "the clip is what makes this finite");
        assert!(wps_fisher_z(f64::NAN).is_nan());
        assert!((wps_fisher_z(0.0)).abs() < 1e-12);
    }

    #[test]
    fn a_zero_sigma_gives_no_z() {
        assert!(wps_scalar_z(1.0, 0.5, 0.0).is_nan());
        assert!(wps_scalar_z(1.0, 0.5, f64::NAN).is_nan());
        assert!(wps_scalar_z(f64::NAN, 0.5, 1.0).is_nan());
        assert!((wps_scalar_z(1.5, 0.5, 0.5) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn the_z_vector_is_absent_where_sigma_is() {
        let z = wps_z_vector(&[1.0, 2.0, 3.0], &[0.0, 0.0, 0.0], &[1.0, 0.0, f64::NAN]);
        assert!((z[0] - 1.0).abs() < 1e-12);
        assert!(z[1].is_nan(), "a zero sigma is not a divisor");
        assert!(z[2].is_nan());
    }

    #[test]
    fn a_profile_too_short_for_the_search_gets_no_shift() {
        let short = vec![1.0; 4 * PHASE_MAX_LAG as usize - 1];
        let (lag, hit) = wps_phase_shift(&short, &short, PHASE_MAX_LAG);
        assert!(lag.is_nan());
        assert!(!hit);
    }

    #[test]
    fn a_shifted_profile_reports_its_shift() {
        // Sign convention, taken from the Python oracle rather than intuition.
        //
        // For lag k the search compares `profile[30+k ..]` against
        // `baseline[30 .. n-30]`. A profile that leads the baseline by 5 --
        // `profile[i] == baseline[i+5]` -- therefore matches at k = -5, since
        // `profile[30-5 ..] == baseline[30 ..]`. The lag is where the *window*
        // sits, not how far the signal moved.
        //
        // This test first asserted +5 on my assumption and failed. Confirmed
        // against `phase_shift` in wps_pon.py on the identical input: -5.
        let n = 200usize;
        let base: Vec<f64> = (0..n).map(|i| (i as f64 / 12.0).sin()).collect();
        let shifted: Vec<f64> = (0..n).map(|i| ((i as f64 + 5.0) / 12.0).sin()).collect();
        let (lag, hit) = wps_phase_shift(&shifted, &base, PHASE_MAX_LAG);
        assert_eq!(lag, -5.0, "sign convention changed; check against the oracle");
        assert!(!hit, "5 is well inside the +/-30 window");
    }

    /// A helper matching `_wave` in tests/unit/test_wps_pon.py.
    fn wave(n: usize, period: f64, offset: f64, scale: f64) -> Vec<f64> {
        (0..n)
            .map(|i| scale * ((i as f64 - offset) / period).sin())
            .collect()
    }

    #[test]
    fn amplitude_is_not_a_coverage_measurement() {
        // Raw peak-to-trough range measures how deeply the sample was
        // sequenced. Measured on the real cohort: raw amplitude correlates
        // +0.512 with local_depth and is skewed 11.6; the log form drops that
        // to -0.036 and 1.6. A z on the raw scale would rank samples by depth.
        let shallow = wps_log_amplitude(&wave(200, 12.0, 0.0, 1.0));
        let deep = wps_log_amplitude(&wave(200, 12.0, 0.0, 100.0));
        assert!(deep > shallow, "a deeper profile must still register larger");
        assert!(
            deep / shallow < 10.0,
            "100x the depth must not be ~100x the statistic ({:.1}x)",
            deep / shallow
        );
    }

    #[test]
    fn the_fisher_transform_removes_a_binding_ceiling() {
        // A correlation is bounded at 1.0. On the real cohort the shape
        // correlation sits at mean 0.844 with sigma 0.099, so the largest
        // attainable positive z is ~1.5 and 302 of 400 anchors could not reach
        // +2 however tumour-like the sample. On the Fisher scale there is no
        // bound to run into.
        assert!(
            (1.0 - 0.844) / 0.0985 < 2.0,
            "the raw ceiling should be binding"
        );
        assert!(wps_fisher_z(0.999) > wps_fisher_z(0.99));
        assert!(wps_fisher_z(0.99) > wps_fisher_z(0.9));
    }

    #[test]
    fn correlation_sees_a_wrong_shape() {
        // The point of the statistic: same amplitude, different structure.
        let r = wps_shape_correlation(&wave(200, 12.0, 0.0, 1.0), &wave(200, 40.0, 0.0, 1.0));
        assert!(r < 0.9, "a different period should not correlate at {r}");
    }

    #[test]
    fn the_search_reports_ending_on_its_own_edge() {
        // `nrl_at_band_limit` one level down: a search that stops at its
        // boundary has found the boundary, not a measurement.
        let ramp: Vec<f64> = (0..200).map(|i| i as f64).collect();
        let inverted: Vec<f64> = ramp.iter().map(|v| -v).collect();
        let (_, hit) = wps_phase_shift(&inverted, &ramp, PHASE_MAX_LAG);
        assert!(hit, "a search ending on the edge must say so");
    }

    #[test]
    fn an_identical_profile_reports_no_shift() {
        let n = 200usize;
        let base: Vec<f64> = (0..n).map(|i| (i as f64 / 12.0).sin()).collect();
        let (lag, hit) = wps_phase_shift(&base, &base, PHASE_MAX_LAG);
        assert_eq!(lag, 0.0);
        assert!(!hit);
    }
}

// ---------------------------------------------------------------------------
// Reading the PON blocks
// ---------------------------------------------------------------------------

use arrow::record_batch::RecordBatch;

/// One anchor's baseline profile: the mean and sigma vectors.
struct WpsVectorEntry {
    mean: Vec<f64>,
    std: Vec<f64>,
}

/// One anchor's baselines for the derived scalars, as `(mean, std)` pairs.
///
/// Mirrors `WPS_SHAPE_STATS` in `pon/model.py`. Kept as a struct rather than a
/// map so that adding a statistic there and forgetting it here is a compile
/// error instead of a column of NaN.
#[derive(Clone, Copy)]
struct WpsShapeEntry {
    log_amplitude: (f64, f64),
    shape_corr_fisher: (f64, f64),
}

/// Everything `apply_pon_zscore` needs out of the PON parquet.
struct WpsPonBaselines {
    /// region_id -> the 200-element mean and sigma profiles
    vectors: HashMap<String, WpsVectorEntry>,
    /// region_id -> the derived-scalar baselines; empty for a pre-0.9.0 PON
    shapes: HashMap<String, WpsShapeEntry>,
}

/// Materialise a Utf8 column as owned strings, tolerating either offset width.
///
/// `table` and `region_id` arrive as `large_string` from the Python builder and
/// as plain `string` from anything written by a bare `pandas.to_parquet` --
/// logically identical, differing only in whether offsets are i32 or i64. A
/// reader that downcasts to one of them reads **nothing** from the other, and
/// reads it silently, because an empty baseline is a legitimate state that the
/// caller degrades on rather than raises.
///
/// That is not hypothetical: every PON this project ships is `large_string`,
/// and `pon_model.rs::PonModel::load` downcasts to `StringArray` only. It is
/// unreachable today, which is the only reason the blindness has never shown
/// up in an output.
fn batch_utf8_column(batch: &RecordBatch, name: &str) -> Option<Vec<Option<String>>> {
    use arrow::array::{Array, LargeStringArray, StringArray};

    let col = batch.column_by_name(name)?;
    if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
        return Some(
            (0..a.len())
                .map(|i| (!a.is_null(i)).then(|| a.value(i).to_string()))
                .collect(),
        );
    }
    if let Some(a) = col.as_any().downcast_ref::<LargeStringArray>() {
        return Some(
            (0..a.len())
                .map(|i| (!a.is_null(i)).then(|| a.value(i).to_string()))
                .collect(),
        );
    }
    warn!(
        "WPS PON: column '{name}' is {:?}, which is not a string type; treating \
         it as absent",
        col.data_type()
    );
    None
}

/// Materialise a numeric column as f64, NaN where null.
///
/// Accepts f64 and f32 because the PON's own writers disagree: the Python
/// builder emits doubles, the Rust builder floats.
fn batch_f64_column(batch: &RecordBatch, name: &str) -> Option<Vec<f64>> {
    use arrow::array::{Array, Float32Array, Float64Array};

    let col = batch.column_by_name(name)?;
    if let Some(a) = col.as_any().downcast_ref::<Float64Array>() {
        return Some(
            (0..a.len())
                .map(|i| if a.is_null(i) { f64::NAN } else { a.value(i) })
                .collect(),
        );
    }
    if let Some(a) = col.as_any().downcast_ref::<Float32Array>() {
        return Some(
            (0..a.len())
                .map(|i| if a.is_null(i) { f64::NAN } else { a.value(i) as f64 })
                .collect(),
        );
    }
    warn!(
        "WPS PON: column '{name}' is {:?}, which is not a float type; treating \
         it as absent",
        col.data_type()
    );
    None
}

/// Materialise a list-of-float column as one `Vec<f64>` per row.
///
/// The shipped PONs store `list<double>`; the Rust builder's own fixtures store
/// `list<float>`. Both are read, and anything else is reported rather than
/// quietly yielding no anchors.
fn batch_vector_column(batch: &RecordBatch, name: &str) -> Option<Vec<Option<Vec<f64>>>> {
    use arrow::array::{Array, Float32Array, Float64Array, ListArray};

    let col = batch.column_by_name(name)?;
    let list = match col.as_any().downcast_ref::<ListArray>() {
        Some(list) => list,
        None => {
            warn!(
                "WPS PON: column '{name}' is {:?}, not a list; treating it as absent",
                col.data_type()
            );
            return None;
        }
    };
    let mut out = Vec::with_capacity(list.len());
    for row in 0..list.len() {
        if list.is_null(row) {
            out.push(None);
            continue;
        }
        let inner = list.value(row);
        if let Some(a) = inner.as_any().downcast_ref::<Float64Array>() {
            out.push(Some(
                (0..a.len())
                    .map(|i| if a.is_null(i) { f64::NAN } else { a.value(i) })
                    .collect(),
            ));
        } else if let Some(a) = inner.as_any().downcast_ref::<Float32Array>() {
            out.push(Some(
                (0..a.len())
                    .map(|i| if a.is_null(i) { f64::NAN } else { a.value(i) as f64 })
                    .collect(),
            ));
        } else {
            warn!(
                "WPS PON: '{name}' holds {:?}, not floats; treating it as absent",
                inner.data_type()
            );
            return None;
        }
    }
    Some(out)
}

/// Read the vector and shape baselines out of a long-format PON parquet.
///
/// Projected down to the eight columns actually needed. A PON is a single wide
/// table of ~120 columns in which each block populates its own few, so reading
/// it whole to reach `wps_nuc_mean` would materialise every other block's
/// columns for all ~130k rows.
fn load_wps_pon_baselines(
    pon_path: &Path,
    vector_table: &str,
    column: &str,
) -> Result<WpsPonBaselines> {
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use parquet::arrow::ProjectionMask;

    let mean_col = format!("{column}_mean");
    let std_col = format!("{column}_std");
    let wanted = [
        "table",
        "region_id",
        mean_col.as_str(),
        std_col.as_str(),
        "log_amplitude_mean",
        "log_amplitude_std",
        "shape_corr_fisher_mean",
        "shape_corr_fisher_std",
    ];

    let file = File::open(pon_path)
        .with_context(|| format!("Failed to open PON file: {pon_path:?}"))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .context("Creating the PON parquet reader")?;

    // Project by *root* index, not by leaf name: a list column's leaf is named
    // "element", so a name-based mask would miss the vectors entirely.
    let roots = builder.parquet_schema().root_schema().get_fields();
    let indices: Vec<usize> = roots
        .iter()
        .enumerate()
        .filter(|(_, f)| wanted.contains(&f.name()))
        .map(|(i, _)| i)
        .collect();
    let found: Vec<&str> = roots
        .iter()
        .map(|f| f.name())
        .filter(|n| wanted.contains(n))
        .collect();
    for name in wanted.iter() {
        if !found.contains(name) {
            debug!("WPS PON: '{name}' is not a column of this PON");
        }
    }
    let mask = ProjectionMask::roots(builder.parquet_schema(), indices);
    let reader = builder
        .with_projection(mask)
        .build()
        .context("Building the projected PON reader")?;

    let mut vectors: HashMap<String, WpsVectorEntry> = HashMap::new();
    let mut shapes: HashMap<String, WpsShapeEntry> = HashMap::new();
    let mut saw_vector_table = false;
    let mut saw_shape_table = false;

    for batch in reader {
        let batch = batch.context("Reading a PON batch")?;
        let tables = match batch_utf8_column(&batch, "table") {
            Some(tables) => tables,
            // Not a long-format PON at all. Loud, because every block would
            // otherwise come back empty and read as "PON has no WPS baseline".
            None => return Err(anyhow!("PON parquet has no readable 'table' column")),
        };
        let ids = batch_utf8_column(&batch, "region_id");
        let means = batch_vector_column(&batch, &mean_col);
        let stds = batch_vector_column(&batch, &std_col);
        let amp_mean = batch_f64_column(&batch, "log_amplitude_mean");
        let amp_std = batch_f64_column(&batch, "log_amplitude_std");
        let corr_mean = batch_f64_column(&batch, "shape_corr_fisher_mean");
        let corr_std = batch_f64_column(&batch, "shape_corr_fisher_std");

        for row in 0..batch.num_rows() {
            let table = match tables[row].as_deref() {
                Some(table) => table,
                None => continue,
            };
            let region_id = match ids.as_ref().and_then(|c| c[row].clone()) {
                Some(id) => id,
                None => continue,
            };
            if table == vector_table {
                saw_vector_table = true;
                if let (Some(mean), Some(std)) = (
                    means.as_ref().and_then(|c| c[row].clone()),
                    stds.as_ref().and_then(|c| c[row].clone()),
                ) {
                    vectors.insert(region_id, WpsVectorEntry { mean, std });
                }
            } else if table == "wps_shape_baseline" {
                saw_shape_table = true;
                let at = |c: &Option<Vec<f64>>| c.as_ref().map_or(f64::NAN, |v| v[row]);
                shapes.insert(
                    region_id,
                    WpsShapeEntry {
                        log_amplitude: (at(&amp_mean), at(&amp_std)),
                        shape_corr_fisher: (at(&corr_mean), at(&corr_std)),
                    },
                );
            }
        }
    }

    // Distinguish "the block is not in this PON" from "the block is there but
    // its columns did not parse". Both give zero anchors; only one is a bug.
    if saw_vector_table && vectors.is_empty() {
        warn!(
            "WPS PON: the '{vector_table}' block is present but no anchor \
             yielded both '{mean_col}' and '{std_col}'. The PON may predate \
             the vector format, or the columns may have an unexpected type."
        );
    }
    if saw_shape_table && shapes.is_empty() {
        warn!("WPS PON: 'wps_shape_baseline' is present but parsed to no anchors.");
    }
    info!(
        "WPS PON: loaded {} anchors from '{vector_table}' and {} from \
         'wps_shape_baseline'",
        vectors.len(),
        shapes.len()
    );
    Ok(WpsPonBaselines { vectors, shapes })
}

// ---------------------------------------------------------------------------
// The entry point
// ---------------------------------------------------------------------------

/// Score a WPS table against the PON, writing an extended copy.
///
/// Adds, per anchor:
///
/// | column | |
/// |---|---|
/// | `{column}_z`                | 200-element z vector, null where the anchor is absent |
/// | `wps_log_amplitude`         | + `_z` |
/// | `wps_shape_corr`            | + `_z`, computed on the Fisher scale |
/// | `wps_phase_shift_bp`        | raw displacement, deliberately not z-scored |
/// | `wps_phase_at_search_limit` | the boundary flag |
///
/// Every input column is carried through untouched: the schema is taken from
/// the file rather than restated here, so a new column added by the WPS writer
/// survives this step without anyone remembering to update it.
///
/// Returns the number of anchors scored. Writes nothing and returns 0 when the
/// PON has no matching vector baseline, matching `region_entropy::apply_pon_zscore`
/// so the Python caller can degrade to a raw output rather than assert on a
/// file that was never created.
#[pyfunction]
#[pyo3(signature = (wps_path, pon_parquet_path, output_path, baseline_table="wps_baseline", column="wps_nuc"))]
pub fn apply_pon_zscore(
    wps_path: PathBuf,
    pon_parquet_path: PathBuf,
    output_path: PathBuf,
    baseline_table: &str,
    column: &str,
) -> PyResult<usize> {
    apply_pon_zscore_inner(
        &wps_path,
        &pon_parquet_path,
        &output_path,
        baseline_table,
        column,
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("WPS PON scoring: {e:#}")))
}

fn apply_pon_zscore_inner(
    wps_path: &Path,
    pon_parquet_path: &Path,
    output_path: &Path,
    baseline_table: &str,
    column: &str,
) -> Result<usize> {
    use arrow::array::{
        ArrayRef, BooleanBuilder, Float64Builder, ListBuilder,
    };
    use arrow::datatypes::{DataType, Field, Schema};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use parquet::arrow::ArrowWriter;
    use parquet::basic::Compression;
    use parquet::file::properties::WriterProperties;

    let baselines = load_wps_pon_baselines(pon_parquet_path, baseline_table, column)?;
    if baselines.vectors.is_empty() {
        warn!(
            "WPS PON: no '{baseline_table}' anchors; {} keeps its raw profile \
             and gets no z-scores. Rebuild the PON with build-pon.",
            wps_path.display()
        );
        return Ok(0);
    }
    if baselines.shapes.is_empty() {
        info!(
            "WPS PON: no 'wps_shape_baseline' (PON built before 0.9.0); the \
             per-position z vectors are still written, the derived shape \
             z-scores are not."
        );
    }

    let file = File::open(wps_path)
        .with_context(|| format!("Failed to open WPS table: {wps_path:?}"))?;
    let wps_builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .context("Creating the WPS parquet reader")?;
    let input_schema = wps_builder.schema().clone();
    let reader = wps_builder.build().context("Building the WPS parquet reader")?;

    // Same order as the Python oracle wrote them, so a reader that indexes by
    // position rather than by name sees no change across the port.
    let z_field = format!("{column}_z");
    let mut fields: Vec<Field> = input_schema.fields().iter().map(|f| f.as_ref().clone()).collect();
    fields.push(Field::new(
        &z_field,
        DataType::List(Arc::new(Field::new("element", DataType::Float64, true))),
        true,
    ));
    for name in [
        "wps_log_amplitude",
        "wps_log_amplitude_z",
        "wps_shape_corr",
        "wps_shape_corr_z",
        "wps_phase_shift_bp",
    ] {
        fields.push(Field::new(name, DataType::Float64, true));
    }
    fields.push(Field::new("wps_phase_at_search_limit", DataType::Boolean, true));
    let out_schema = Arc::new(Schema::new(fields));

    // Snappy, matching what the Python oracle wrote via pandas. The default is
    // uncompressed, and these tables are 200 doubles an anchor across ~89k
    // anchors -- the one place in this pipeline where the codec is load-bearing.
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();
    let out_file = File::create(output_path)
        .with_context(|| format!("Creating {output_path:?}"))?;
    let mut writer = ArrowWriter::try_new(out_file, out_schema.clone(), Some(props))
        .context("Creating the WPS output writer")?;

    let mut n_rows = 0usize;
    let mut n_scored = 0usize;
    let mut n_absent = 0usize;
    let mut n_mismatched = 0usize;
    let mut n_at_limit = 0usize;

    for batch in reader {
        let batch = batch.context("Reading a WPS batch")?;
        let ids = batch_utf8_column(&batch, "region_id")
            .ok_or_else(|| anyhow!("WPS table has no readable 'region_id' column"))?;
        let profiles = batch_vector_column(&batch, column)
            .ok_or_else(|| anyhow!("WPS table has no readable '{column}' column"))?;

        let rows = batch.num_rows();
        // `.with_field`, not the default: ListBuilder names its inner field
        // "item" while the schema above declares "element" -- and arrow
        // compares list fields by name, so the two disagree at assembly. The
        // name is "element" because that is what the Python oracle wrote, and
        // a downstream reader that checks the type would see a changed schema
        // otherwise.
        let mut z_builder = ListBuilder::new(Float64Builder::new())
            .with_field(Arc::new(Field::new("element", DataType::Float64, true)));
        let mut amp = Float64Builder::with_capacity(rows);
        let mut amp_z = Float64Builder::with_capacity(rows);
        let mut corr = Float64Builder::with_capacity(rows);
        let mut corr_z = Float64Builder::with_capacity(rows);
        let mut shift = Float64Builder::with_capacity(rows);
        let mut limit = BooleanBuilder::with_capacity(rows);

        for row in 0..rows {
            n_rows += 1;
            let profile = profiles[row].clone().unwrap_or_default();

            // Amplitude needs no baseline, so it is measured for every anchor
            // -- including those absent from the PON, which then carry a raw
            // reading with no z beside it.
            let amplitude = wps_log_amplitude(&profile);
            amp.append_value(amplitude);

            let entry = ids[row]
                .as_deref()
                .and_then(|id| baselines.vectors.get(id));
            // A baseline profile of a different length than the sample is a
            // corrupt pairing, not a measurement. Truncating to the shorter of
            // the two -- which is what zipping them would do -- yields a
            // full-looking z vector computed from misaligned positions, and
            // nothing downstream could tell. The Python raises a broadcast
            // error here and takes the whole sample down with it; this reports
            // the count and scores the rest.
            let entry = entry.filter(|e| {
                let ok = e.mean.len() == profile.len() && e.std.len() == profile.len();
                if !ok {
                    n_mismatched += 1;
                }
                ok
            });
            let entry = match entry {
                Some(entry) => entry,
                None => {
                    n_absent += 1;
                    z_builder.append_null();
                    corr.append_value(f64::NAN);
                    shift.append_value(f64::NAN);
                    limit.append_value(false);
                    amp_z.append_value(f64::NAN);
                    corr_z.append_value(f64::NAN);
                    continue;
                }
            };

            for z in wps_z_vector(&profile, &entry.mean, &entry.std) {
                z_builder.values().append_value(z);
            }
            z_builder.append(true);

            let r = wps_shape_correlation(&profile, &entry.mean);
            corr.append_value(r);
            let (lag, hit) = wps_phase_shift(&profile, &entry.mean, PHASE_MAX_LAG);
            shift.append_value(lag);
            limit.append_value(hit);
            n_at_limit += usize::from(hit);
            n_scored += 1;

            match ids[row].as_deref().and_then(|id| baselines.shapes.get(id)) {
                Some(shape) => {
                    amp_z.append_value(wps_scalar_z(
                        amplitude,
                        shape.log_amplitude.0,
                        shape.log_amplitude.1,
                    ));
                    corr_z.append_value(wps_scalar_z(
                        wps_fisher_z(r),
                        shape.shape_corr_fisher.0,
                        shape.shape_corr_fisher.1,
                    ));
                }
                // Absent from the shape baseline reads the same as unmeasured.
                None => {
                    amp_z.append_value(f64::NAN);
                    corr_z.append_value(f64::NAN);
                }
            }
        }

        let mut columns: Vec<ArrayRef> = batch.columns().to_vec();
        columns.push(Arc::new(z_builder.finish()));
        columns.push(Arc::new(amp.finish()));
        columns.push(Arc::new(amp_z.finish()));
        columns.push(Arc::new(corr.finish()));
        columns.push(Arc::new(corr_z.finish()));
        columns.push(Arc::new(shift.finish()));
        columns.push(Arc::new(limit.finish()));
        let out_batch = RecordBatch::try_new(out_schema.clone(), columns)
            .context("Assembling the scored WPS batch")?;
        writer.write(&out_batch).context("Writing a scored batch")?;
    }
    writer.close().context("Closing the WPS output writer")?;

    info!(
        "WPS PON: {n_scored}/{n_rows} anchors scored ({})",
        wps_path.file_name().unwrap_or_default().to_string_lossy()
    );
    if n_absent > 0 {
        warn!(
            "WPS PON: {n_absent}/{n_rows} anchors are absent from the baseline \
             and get no z-score. Anchors backed by fewer than 3 samples are \
             dropped at build time, so this is expected to be large for duplex \
             PONs."
        );
    }
    if n_mismatched > 0 {
        warn!(
            "WPS PON: {n_mismatched}/{n_rows} anchors have a baseline profile \
             of a different length than the sample and were left unscored. \
             This is a PON/sample mismatch -- check that the PON was built \
             with the same anchor set."
        );
    }
    if n_at_limit > 0 {
        warn!(
            "WPS PON: {n_at_limit}/{n_rows} anchors hit the +/-{PHASE_MAX_LAG} \
             phase-search window. wps_phase_shift_bp is the edge of the search \
             there, not a measurement -- see wps_phase_at_search_limit."
        );
    }
    Ok(n_scored)
}
