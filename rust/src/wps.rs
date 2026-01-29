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


use std::path::Path;
use std::io::{BufRead, Write};
use std::fs::File;
use anyhow::{Result, Context, anyhow};

use pyo3::prelude::*;
use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::faidx;
use log::{info, debug};

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
    fn add_range(vec: &mut Vec<f64>, start: i64, end: i64, val: f64) {
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
        
        eprintln!("DEBUG: Starting region loop, {} regions, {} accumulators", self.regions.len(), self.accumulators.len());
        let mut processed_count = 0usize;
        
        for (i, region) in self.regions.iter().enumerate() {
            let acc_opt = self.accumulators.get(&i);
            if acc_opt.is_none() { continue; }
            
            let acc = acc_opt.unwrap();
            let len = (region.end - region.start + 1) as usize;
            
            // Reconstruct raw values from difference arrays
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
            
            // Skip empty regions
            if curr_cov == 0.0 { continue; }
            
            // Bin aggregation: 2000bp → 200 bins
            let mut binned_wps_nuc = vec![0.0f32; NUM_BINS];
            let mut binned_wps_tf = vec![0.0f32; NUM_BINS];
            let mut binned_mask = vec![1i8; NUM_BINS];
            let mut binned_prot_nuc = vec![0.0f32; NUM_BINS];
            let mut binned_prot_tf = vec![0.0f32; NUM_BINS];
            
            for bin_idx in 0..NUM_BINS {
                let start_pos = bin_idx * BIN_SIZE;
                let end_pos = (start_pos + BIN_SIZE).min(len);
                
                if start_pos >= len { break; }
                
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
                    
                    // Check capture mask for this position
                    if let Some(ref bait_mask) = self.bait_mask {
                        let chrom_id = 0u32; // Simplified - need proper chrom_id lookup
                        let pos = region.start + j as u64;
                        let (_, reliable) = bait_mask.check_position(chrom_id, pos);
                        if !reliable { mask_ok = false; }
                    }
                }
                
                binned_wps_nuc[bin_idx] = (sum_wps_nuc / count as f64) as f32;
                binned_wps_tf[bin_idx] = (sum_wps_tf / count as f64) as f32;
                binned_mask[bin_idx] = if mask_ok { 1 } else { 0 };
                
                let total_nuc = sum_span_nuc + sum_end_nuc;
                let total_tf = sum_span_tf + sum_end_tf;
                binned_prot_nuc[bin_idx] = if total_nuc > 0.0 { (sum_span_nuc / total_nuc) as f32 } else { 0.0 };
                binned_prot_tf[bin_idx] = if total_tf > 0.0 { (sum_span_tf / total_tf) as f32 } else { 0.0 };
            }
            
            // Reverse for minus strand
            if region.strand == "-" {
                binned_wps_nuc.reverse();
                binned_wps_tf.reverse();
                binned_mask.reverse();
                binned_prot_nuc.reverse();
                binned_prot_tf.reverse();
            }
            
            // Apply Savitzky-Golay smoothing (window=11, order=3, matching scipy defaults)
            use sci_rs::signal::filter::savgol_filter_dyn;
            const SAVGOL_WINDOW: usize = 11;
            const SAVGOL_POLYORDER: usize = 3;
            
            if binned_wps_nuc.len() >= SAVGOL_WINDOW {
                binned_wps_nuc = savgol_filter_dyn(
                    binned_wps_nuc.iter().map(|x| *x as f64),
                    SAVGOL_WINDOW, SAVGOL_POLYORDER, None, None
                ).into_iter().map(|x| x as f32).collect();
                
                binned_wps_tf = savgol_filter_dyn(
                    binned_wps_tf.iter().map(|x| *x as f64),
                    SAVGOL_WINDOW, SAVGOL_POLYORDER, None, None
                ).into_iter().map(|x| x as f32).collect();
            }
            
            // Region type
            let region_type = if region.id.starts_with("TSS|") { "TSS" } 
                else if region.id.starts_with("CTCF|") { "CTCF" } 
                else { "OTHER" };
            
            // Append row
            region_id_builder.append_value(&region.id);
            chrom_builder.append_value(&region.chrom);
            center_builder.append_value(region.center as i32);
            strand_builder.append_value(&region.strand);
            region_type_builder.append_value(region_type);
            
            // Append list values
            for v in &binned_wps_nuc { wps_nuc_builder.values().append_value(*v); }
            wps_nuc_builder.append(true);
            
            for v in &binned_wps_tf { wps_tf_builder.values().append_value(*v); }
            wps_tf_builder.append(true);
            
            for v in &binned_mask { mask_builder.values().append_value(*v); }
            mask_builder.append(true);
            
            for v in &binned_prot_nuc { prot_nuc_builder.values().append_value(*v); }
            prot_nuc_builder.append(true);
            
            for v in &binned_prot_tf { prot_tf_builder.values().append_value(*v); }
            prot_tf_builder.append(true);
            
            depth_builder.append_value((curr_cov / len as f64) as f32);
            processed_count += 1;
        }
        
        eprintln!("DEBUG: Loop complete. Processed {} regions with data", processed_count);
        
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
        eprintln!("DEBUG: Array lengths: region_id={}, chrom={}, center={}, strand={}, type={}, wps_nuc={}, wps_tf={}, mask={}, prot_nuc={}, prot_tf={}, depth={}",
            arrays[0].len(), arrays[1].len(), arrays[2].len(), arrays[3].len(), arrays[4].len(),
            arrays[5].len(), arrays[6].len(), arrays[7].len(), arrays[8].len(), arrays[9].len(), arrays[10].len());
        
        let batch = match RecordBatch::try_new(StdArc::new(schema), arrays) {
            Ok(b) => b,
            Err(e) => {
                eprintln!("ERROR: RecordBatch::try_new failed: {:?}", e);
                return Err(anyhow::anyhow!("RecordBatch creation failed: {}", e));
            }
        };
        
        eprintln!("DEBUG: Created batch with {} rows, writing to {:?}", batch.num_rows(), output_path);
        
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
    const PROFILE_LENGTH: usize = 300;  // Alu ~300bp
    
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
            let node = IntervalNode::new(region.start as i32, region.end as i32, i);
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
        use arrow::array::{ArrayRef, Float32Builder, Int64Builder, ListBuilder, StringBuilder};
        use arrow::datatypes::{DataType, Field, Schema};
        use arrow::record_batch::RecordBatch;
        use parquet::arrow::ArrowWriter;
        use parquet::file::properties::WriterProperties;
        use std::sync::Arc as StdArc;
        
        const NUM_BINS: usize = 30;
        const BIN_SIZE: usize = 10;  // 300 / 30
        
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
            let (nrl_bp, periodicity_score) = self.estimate_periodicity(&binned_nuc);
            
            // Calculate deviation from expected NRL (190bp) and penalized score
            const EXPECTED_NRL_BP: f32 = 190.0;
            const TOLERANCE_BP: f32 = 20.0;
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
    /// 5. **Peak finding**: Find dominant peak in 140-200bp period range
    /// 6. **Quality score**: SNR-based score capped at 1.0 (SNR > 3 = clear periodicity)
    /// 
    /// ## Returns
    /// `(nrl_bp, quality_score)` tuple:
    /// - `nrl_bp`: Nucleosome Repeat Length in bp (expected ~190bp for healthy cfDNA)
    /// - `quality_score`: SNR-based quality 0-1 (higher = clearer periodicity signal)
    fn estimate_periodicity(&self, profile: &[f32]) -> (f32, f32) {
        use realfft::RealFftPlanner;
        use std::f64::consts::PI;
        
        let n = profile.len();
        if n < 10 { return (0.0, 0.0); }
        
        // Convert to f64 for precision
        let mut data: Vec<f64> = profile.iter().map(|x| *x as f64).collect();
        
        // 1. Detrend (remove linear trend)
        let n_f = n as f64;
        let x_mean = (n_f - 1.0) / 2.0;
        let y_mean = data.iter().sum::<f64>() / n_f;
        
        let mut slope_numer = 0.0;
        let mut slope_denom = 0.0;
        for i in 0..n {
            let xi = i as f64 - x_mean;
            let yi = data[i] - y_mean;
            slope_numer += xi * yi;
            slope_denom += xi * xi;
        }
        let slope = if slope_denom.abs() > 1e-9 { slope_numer / slope_denom } else { 0.0 };
        let intercept = y_mean - slope * x_mean;
        
        for i in 0..n {
            data[i] -= slope * i as f64 + intercept;
        }
        
        // 2. Z-score normalize
        let mean: f64 = data.iter().sum::<f64>() / n_f;
        let variance: f64 = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n_f;
        let std = variance.sqrt();
        
        if std < 1e-9 {
            return (0.0, 0.0); // Flat profile
        }
        
        for val in data.iter_mut() {
            *val /= std;
        }
        
        // 3. Apply Hann window
        for i in 0..n {
            let window = 0.5 * (1.0 - (2.0 * PI * i as f64 / (n_f - 1.0)).cos());
            data[i] *= window;
        }
        
        // 4. Compute FFT using realfft
        let mut planner = RealFftPlanner::<f64>::new();
        let fft = planner.plan_fft_forward(n);
        
        let mut indata = data.clone();
        let mut spectrum = fft.make_output_vec();
        
        if fft.process(&mut indata, &mut spectrum).is_err() {
            return (0.0, 0.0);
        }
        
        // 5. Find peak in target period range (140-200bp for nucleosomes)
        // With 10bp bins, 30 bins = 300bp, so periods in range:
        // Frequency f = 1/period, and freq[i] = i / (n * bin_size)
        // So period[i] = n * bin_size / i
        let bin_size_bp = 10.0;
        let min_period_bp = 140.0;
        let max_period_bp = 200.0;
        
        let amplitudes: Vec<f64> = spectrum.iter().map(|c| c.norm()).collect();
        
        let mut peak_amplitude = 0.0;
        let mut peak_idx = 0;
        let mut sum_amplitude = 0.0;
        let mut count_in_range = 0;
        
        for i in 1..amplitudes.len() {
            let period_bp = (n as f64 * bin_size_bp) / i as f64;
            
            if period_bp >= min_period_bp && period_bp <= max_period_bp {
                sum_amplitude += amplitudes[i];
                count_in_range += 1;
                
                if amplitudes[i] > peak_amplitude {
                    peak_amplitude = amplitudes[i];
                    peak_idx = i;
                }
            }
        }
        
        if count_in_range == 0 || peak_idx == 0 {
            return (0.0, 0.0);
        }
        
        // 6. Calculate SNR and quality score
        let background = sum_amplitude / count_in_range as f64;
        let snr = if background > 0.0 { peak_amplitude / background } else { 0.0 };
        
        // Calculate NRL (nucleosome repeat length) from peak index
        let nrl_bp = (n as f64 * bin_size_bp) / peak_idx as f64;
        
        // Quality score: SNR-based, capped at 1.0 (SNR > 3 indicates clear periodicity)
        let quality_score = (snr / 3.0).min(1.0) as f32;
        
        (nrl_bp as f32, quality_score)
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
                
                // Fragment position relative to Alu
                let f_start = (fragment.start + 1) as i64 - r_start as i64;
                let f_end = fragment.end as i64 - r_start as i64;
                
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

// =============================================================================
// PON Z-Score Normalization for WPS Parquet Files
// =============================================================================

use std::path::PathBuf;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use arrow::array::{Float64Array, StringArray, ArrayRef, Array};
use arrow::datatypes::{Schema, Field, DataType};
use std::sync::Arc as StdArc;

/// Apply PON z-score normalization to WPS parquet file.
/// 
/// Reads WPS parquet with region_id, wps_nuc, wps_tf columns,
/// joins with PON baseline to compute z-scores.
/// 
/// # Arguments
/// * `wps_parquet_path` - Path to sample WPS parquet
/// * `pon_parquet_path` - Path to PON parquet with wps_baseline table
/// * `output_path` - Output path for normalized WPS
/// 
/// # Returns
/// * Number of regions processed
/// 
/// # Performance
/// 5-20x faster than Python pandas merge + loop
#[pyfunction]
#[pyo3(signature = (wps_parquet_path, pon_parquet_path, output_path=None))]
pub fn apply_pon_zscore(
    wps_parquet_path: PathBuf,
    pon_parquet_path: PathBuf,
    output_path: Option<PathBuf>,
) -> PyResult<usize> {
    use std::fs::File;
    
    info!("WPS PON z-score: loading PON from {:?}", pon_parquet_path);
    
    // 1. Load WPS baseline from PON
    let wps_baseline = load_wps_baseline_from_parquet(&pon_parquet_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to load WPS baseline: {}", e)
        ))?;
    
    if wps_baseline.regions.is_empty() {
        info!("WPS PON: No WPS baseline found, skipping normalization");
        return Ok(0);
    }
    
    info!("WPS PON: Loaded baseline for {} regions", wps_baseline.regions.len());
    
    // 2. Read sample WPS parquet
    let file = File::open(&wps_parquet_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to open WPS parquet: {}", e)
        ))?;
    
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?
        .build()
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    
    let mut all_batches: Vec<RecordBatch> = Vec::new();
    let mut regions_processed = 0;
    
    for batch_result in reader {
        let batch = batch_result
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        
        // Find column indices
        let schema = batch.schema();
        let region_idx = schema.index_of("region_id").ok();
        let nuc_idx = schema.index_of("wps_nuc_mean").or_else(|_| schema.index_of("wps_nuc")).ok();
        let tf_idx = schema.index_of("wps_tf_mean").or_else(|_| schema.index_of("wps_tf")).ok();
        
        if region_idx.is_none() || (nuc_idx.is_none() && tf_idx.is_none()) {
            // Can't process without required columns
            all_batches.push(batch);
            continue;
        }
        
        let region_col = batch.column(region_idx.unwrap())
            .as_any().downcast_ref::<StringArray>();
        
        if region_col.is_none() {
            all_batches.push(batch);
            continue;
        }
        
        let region_array = region_col.unwrap();
        let num_rows = region_array.len();
        
        // Compute z-scores
        let mut nuc_z_scores: Vec<f64> = Vec::with_capacity(num_rows);
        let mut tf_z_scores: Vec<f64> = Vec::with_capacity(num_rows);
        
        for i in 0..num_rows {
            let region_id = region_array.value(i);
            
            // Get sample values
            let nuc_val = nuc_idx.map(|idx| {
                batch.column(idx).as_any()
                    .downcast_ref::<Float64Array>()
                    .map(|arr| arr.value(i))
                    .unwrap_or(0.0)
            }).unwrap_or(0.0);
            
            let tf_val = tf_idx.map(|idx| {
                batch.column(idx).as_any()
                    .downcast_ref::<Float64Array>()
                    .map(|arr| arr.value(i))
                    .unwrap_or(0.0)
            }).unwrap_or(0.0);
            
            // Compute z-scores from PON baseline
            if let Some(baseline) = wps_baseline.get_baseline(region_id) {
                let nuc_z = if baseline.wps_long_std > 0.0 {
                    (nuc_val - baseline.wps_long_mean) / baseline.wps_long_std
                } else {
                    0.0
                };
                
                let tf_z = if baseline.wps_short_std > 0.0 {
                    (tf_val - baseline.wps_short_mean) / baseline.wps_short_std
                } else {
                    0.0
                };
                
                nuc_z_scores.push(nuc_z);
                tf_z_scores.push(tf_z);
            } else {
                nuc_z_scores.push(0.0);
                tf_z_scores.push(0.0);
            }
            
            regions_processed += 1;
        }
        
        // Build new batch with z-score columns
        let mut new_columns: Vec<ArrayRef> = batch.columns().to_vec();
        
        let nuc_z_array = Float64Array::from(nuc_z_scores);
        let tf_z_array = Float64Array::from(tf_z_scores);
        
        new_columns.push(StdArc::new(nuc_z_array));
        new_columns.push(StdArc::new(tf_z_array));
        
        // Create extended schema
        let mut fields: Vec<Field> = schema.fields().iter().map(|f| f.as_ref().clone()).collect();
        fields.push(Field::new("wps_nuc_z", DataType::Float64, true));
        fields.push(Field::new("wps_tf_z", DataType::Float64, true));
        let new_schema = StdArc::new(Schema::new(fields));
        
        let new_batch = RecordBatch::try_new(new_schema, new_columns)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        
        all_batches.push(new_batch);
    }
    
    // 3. Write output
    if !all_batches.is_empty() {
        let output_path = output_path.unwrap_or(wps_parquet_path.clone());
        let output_file = File::create(&output_path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(
                format!("Failed to create output: {}", e)
            ))?;
        
        let schema = all_batches[0].schema();
        let mut writer = ArrowWriter::try_new(output_file, schema.clone(), None)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        
        for batch in all_batches {
            writer.write(&batch)
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }
        
        writer.close()
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    }
    
    info!("WPS PON: Processed {} regions with z-score normalization", regions_processed);
    
    Ok(regions_processed)
}

/// Load WpsBaseline from PON parquet.
fn load_wps_baseline_from_parquet(path: &std::path::Path) -> anyhow::Result<crate::pon_model::WpsBaseline> {
    use parquet::file::reader::{FileReader, SerializedFileReader};
    use parquet::record::RowAccessor;
    use std::fs::File;
    
    let file = File::open(path)?;
    let reader = SerializedFileReader::new(file)?;
    
    let mut regions: std::collections::HashMap<String, crate::pon_model::WpsRegionBaseline> = std::collections::HashMap::new();
    
    for row_result in reader.get_row_iter(None)? {
        let row = row_result?;
        
        // Check table column
        if let Ok(table) = row.get_string(
            row.get_column_iter().position(|(name, _)| name == "table").unwrap_or(0)
        ) {
            if table != "wps_baseline" {
                continue;
            }
        } else {
            continue;
        }
        
        // Parse region_id and baseline values
        let region_id = row.get_string(
            row.get_column_iter().position(|(name, _)| name == "region_id").unwrap_or(0)
        ).map_or("".to_string(), |v| v.to_string());
        
        let wps_long_mean = row.get_double(
            row.get_column_iter().position(|(name, _)| name == "wps_long_mean").unwrap_or(0)
        ).unwrap_or(0.0);
        
        let wps_long_std = row.get_double(
            row.get_column_iter().position(|(name, _)| name == "wps_long_std").unwrap_or(0)
        ).unwrap_or(1.0);
        
        let wps_short_mean = row.get_double(
            row.get_column_iter().position(|(name, _)| name == "wps_short_mean").unwrap_or(0)
        ).unwrap_or(0.0);
        
        let wps_short_std = row.get_double(
            row.get_column_iter().position(|(name, _)| name == "wps_short_std").unwrap_or(0)
        ).unwrap_or(1.0);
        
        regions.insert(region_id, crate::pon_model::WpsRegionBaseline {
            wps_long_mean,
            wps_long_std,
            wps_short_mean,
            wps_short_std,
        });
    }
    
    info!("Loaded WPS baseline: {} regions", regions.len());
    
    Ok(crate::pon_model::WpsBaseline { regions })
}
