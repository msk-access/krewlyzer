//! GC Reference Asset Generation
//!
//! Generates pre-computed GC correction assets:
//! 1. valid_regions.bed - Curated 100kb bins excluding blacklist/gaps
//! 2. ref_genome_GC.parquet - Expected fragment counts per (length, GC)
//!
//! These assets are generated once per reference genome and shipped with krewlyzer.

use std::path::Path;
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::fs::File;
use std::collections::HashMap;

use anyhow::{Result, anyhow, Context};
use log::{info, debug, warn};
use pyo3::prelude::*;
use rust_htslib::faidx;

/// Length bins for GC correction (20bp bins from 60-400bp = 17 bins)
pub const LENGTH_BINS: [(u32, u32); 17] = [
    (60, 80), (80, 100), (100, 120), (120, 140), (140, 160),
    (160, 180), (180, 200), (200, 220), (220, 240), (240, 260),
    (260, 280), (280, 300), (300, 320), (320, 340), (340, 360),
    (360, 380), (380, 400),
];

/// Get the length bin index for a given fragment length
/// Returns None if the length is outside the range [60, 400]
pub fn get_length_bin_index(length: u32) -> Option<usize> {
    if length < 60 || length >= 400 {
        return None;
    }
    Some(((length - 60) / 20) as usize)
}

/// GC correction factors lookup table
/// Indexed by [length_bin][gc_percent]
#[derive(Debug, Clone)]
pub struct GcCorrectionFactors {
    /// Correction factors: [length_bin][gc_percent] -> factor
    /// gc_percent is 0-100 (101 values)
    pub factors: Vec<Vec<f64>>,
    
    /// Number of fragments observed per (length_bin, gc_percent)
    pub counts: Vec<Vec<u64>>,
}

impl GcCorrectionFactors {
    /// Create a new empty correction factors table
    pub fn new() -> Self {
        Self {
            factors: vec![vec![1.0; 101]; 17],
            counts: vec![vec![0; 101]; 17],
        }
    }
    
    /// Get the correction factor for a given length and GC content
    /// 
    /// # Arguments
    /// * `length` - Fragment length in bp
    /// * `gc` - GC content as fraction 0.0-1.0
    /// 
    /// # Returns
    /// * Correction factor (1.0 if no data available)
    pub fn get_factor(&self, length: u32, gc: f64) -> f64 {
        let bin_idx = match get_length_bin_index(length) {
            Some(idx) => idx,
            None => return 1.0,
        };
        
        let gc_percent = (gc * 100.0).round() as usize;
        let gc_idx = gc_percent.min(100);
        
        self.factors[bin_idx][gc_idx]
    }
}

/// Expected fragment counts from reference genome
/// Used for rate = observed / expected calculation
#[derive(Debug, Clone)]
pub struct RefGenomeGc {
    /// Expected counts: [length_bin][gc_percent] -> expected_count
    pub expected: Vec<Vec<f64>>,
    
    /// Reference genome name (e.g., "hg38")
    pub genome_name: String,
}

impl RefGenomeGc {
    /// Load reference GC expected values from a Parquet file
    /// 
    /// # Arguments
    /// * `path` - Path to the .parquet file
    /// 
    /// # Returns
    /// * Loaded RefGenomeGc or error
    pub fn load(path: &Path) -> Result<Self> {
        // TODO: Implement loading from Parquet file
        // Schema: length_bin (u32), gc_percent (u8), expected_count (f64)
        info!("Loading ref_genome_GC from {:?}", path);
        
        Err(anyhow!("ref_genome_GC.parquet loading not yet implemented"))
    }
    
    /// Get expected count for a given length and GC content
    pub fn get_expected(&self, length: u32, gc: f64) -> f64 {
        let bin_idx = match get_length_bin_index(length) {
            Some(idx) => idx,
            None => return 1.0,
        };
        
        let gc_percent = (gc * 100.0).round() as usize;
        let gc_idx = gc_percent.min(100);
        
        self.expected[bin_idx][gc_idx]
    }
}

/// Valid genomic regions for GC correction estimation
/// Excludes blacklist regions, assembly gaps, and low-mappability areas
#[derive(Debug, Clone)]
pub struct ValidRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub gc_content: f64,
}

/// Load valid regions from a BED file
/// 
/// # Arguments
/// * `bed_path` - Path to the valid_regions.bed file
/// 
/// # Returns
/// * Vector of ValidRegion structs
pub fn load_valid_regions(bed_path: &Path) -> Result<Vec<ValidRegion>> {
    info!("Loading valid regions from {:?}", bed_path);
    
    let file = File::open(bed_path)
        .with_context(|| format!("Failed to open valid regions file: {:?}", bed_path))?;
    let reader = BufReader::new(file);
    
    let mut regions = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        
        let chrom = fields[0].trim_start_matches("chr").to_string();
        let start: u64 = fields[1].parse().unwrap_or(0);
        let end: u64 = fields[2].parse().unwrap_or(0);
        
        // GC content may be in 4th column if pre-computed
        let gc_content = if fields.len() >= 4 {
            fields[3].parse().unwrap_or(0.0)
        } else {
            0.0
        };
        
        regions.push(ValidRegion {
            chrom,
            start,
            end,
            gc_content,
        });
    }
    
    info!("Loaded {} valid regions", regions.len());
    Ok(regions)
}

/// Generate valid regions BED file by excluding blacklist and gap regions
/// 
/// # Arguments
/// * `reference_path` - Path to reference FASTA
/// * `blacklist_path` - Path to ENCODE blacklist BED
/// * `output_path` - Path to write valid_regions.bed
/// * `bin_size` - Size of each bin (default: 100000)
/// 
/// # Returns
/// * Number of valid regions generated
#[pyfunction]
#[pyo3(signature = (reference_path, blacklist_path, output_path, bin_size=100000))]
pub fn generate_valid_regions(
    reference_path: &str,
    blacklist_path: &str,
    output_path: &str,
    bin_size: u64,
) -> PyResult<usize> {
    info!("Generating valid regions from reference: {}", reference_path);
    info!("Blacklist: {}", blacklist_path);
    info!("Bin size: {} bp", bin_size);
    
    // Load reference index to get chromosome lengths
    let faidx = faidx::Reader::from_path(reference_path)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to open reference FASTA: {}. Make sure .fai exists.", e)
        ))?;
    
    // Load blacklist regions
    let blacklist = load_blacklist(Path::new(blacklist_path))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to load blacklist: {}", e)
        ))?;
    
    // Valid chromosomes (1-22, X, Y)
    let valid_chroms: Vec<String> = (1..=22)
        .map(|i| i.to_string())
        .chain(vec!["X".to_string(), "Y".to_string()])
        .collect();
    
    let mut regions: Vec<ValidRegion> = Vec::new();
    
    // TODO: Get chromosome lengths from faidx and generate bins
    // For now, placeholder implementation
    
    // Write output BED
    let output_file = File::create(output_path)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to create output file: {}", e)
        ))?;
    let mut writer = BufWriter::new(output_file);
    
    // Write header
    writeln!(writer, "# Valid regions for GC correction estimation")
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    writeln!(writer, "# Excludes: ENCODE blacklist, assembly gaps, low mappability")
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    
    // TODO: Implement actual region generation
    let region_count = regions.len();
    
    info!("Generated {} valid regions", region_count);
    Ok(region_count)
}

/// Load blacklist regions from BED file
fn load_blacklist(path: &Path) -> Result<HashMap<String, Vec<(u64, u64)>>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open blacklist file: {:?}", path))?;
    let reader = BufReader::new(file);
    
    let mut blacklist: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        
        let chrom = fields[0].trim_start_matches("chr").to_string();
        let start: u64 = fields[1].parse().unwrap_or(0);
        let end: u64 = fields[2].parse().unwrap_or(0);
        
        blacklist.entry(chrom).or_insert_with(Vec::new).push((start, end));
    }
    
    // Sort intervals by start position
    for intervals in blacklist.values_mut() {
        intervals.sort_by_key(|i| i.0);
    }
    
    info!("Loaded blacklist with {} chromosomes", blacklist.len());
    Ok(blacklist)
}

/// Generate reference genome GC expected counts
/// 
/// # Algorithm
/// For each length bin (60-80, 80-100, ..., 380-400):
///   For each GC percent (0-100):
///     Count how many theoretical fragments of that length have that GC
///     
/// This represents the "expected" distribution under uniform sampling.
/// 
/// # Arguments
/// * `reference_path` - Path to reference FASTA
/// * `valid_regions_path` - Path to valid_regions.bed
/// * `output_path` - Path to write ref_genome_GC.parquet
/// 
/// # Returns
/// * Total fragments counted
#[pyfunction]
#[pyo3(signature = (reference_path, valid_regions_path, output_path))]
pub fn generate_ref_genome_gc(
    reference_path: &str,
    valid_regions_path: &str,
    output_path: &str,
) -> PyResult<u64> {
    info!("Generating ref_genome_GC from reference: {}", reference_path);
    info!("Using valid regions: {}", valid_regions_path);
    
    // Load reference
    let faidx = faidx::Reader::from_path(reference_path)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to open reference FASTA: {}", e)
        ))?;
    
    // Load valid regions
    let regions = load_valid_regions(Path::new(valid_regions_path))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to load valid regions: {}", e)
        ))?;
    
    // Initialize counts matrix [length_bin][gc_percent]
    let mut counts: Vec<Vec<u64>> = vec![vec![0; 101]; 17];
    let mut total_fragments: u64 = 0;
    
    info!("Processing {} regions...", regions.len());
    
    // TODO: Implement actual GC counting from reference
    // For each region:
    //   For each position in region:
    //     For each length bin:
    //       Extract sequence, compute GC, increment count
    
    // TODO: Save as Parquet file (consistent with PON format)
    // Schema: length_bin (u32), gc_percent (u8), expected_count (u64)
    
    info!("Generated ref_genome_GC with {} total fragments", total_fragments);
    Ok(total_fragments)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_length_bin_index() {
        assert_eq!(get_length_bin_index(60), Some(0));
        assert_eq!(get_length_bin_index(79), Some(0));
        assert_eq!(get_length_bin_index(80), Some(1));
        assert_eq!(get_length_bin_index(399), Some(16));
        assert_eq!(get_length_bin_index(400), None);
        assert_eq!(get_length_bin_index(59), None);
    }
    
    #[test]
    fn test_gc_correction_factors() {
        let factors = GcCorrectionFactors::new();
        
        // Default factor should be 1.0
        assert_eq!(factors.get_factor(100, 0.5), 1.0);
        assert_eq!(factors.get_factor(200, 0.3), 1.0);
    }
}
