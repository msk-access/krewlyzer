//! GC Reference Asset Generation
//!
//! Generates pre-computed GC correction assets:
//! 1. valid_regions.bed - Curated 100kb bins excluding problematic regions
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

/// Valid genomic regions for GC correction estimation
/// Excludes problematic regions, assembly gaps, and low-mappability areas
#[derive(Debug, Clone)]
pub struct ValidRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub gc_content: f64,
}

/// Load valid regions from a BED file (plain or BGZF compressed)
/// 
/// # Arguments
/// * `bed_path` - Path to the valid_regions.bed or .bed.gz file
/// 
/// # Returns
/// * Vector of ValidRegion structs
pub fn load_valid_regions(bed_path: &Path) -> Result<Vec<ValidRegion>> {
    use noodles::bgzf;
    use std::io::BufRead;
    
    info!("Loading valid regions from {:?}", bed_path);
    
    let file = File::open(bed_path)
        .with_context(|| format!("Failed to open valid regions file: {:?}", bed_path))?;
    
    // Check if file is BGZF compressed based on extension
    let is_bgzf = bed_path.extension()
        .map(|ext| ext == "gz")
        .unwrap_or(false);
    
    let mut regions = Vec::new();
    
    // Create appropriate reader based on compression
    if is_bgzf {
        let mut reader = bgzf::io::Reader::new(file);
        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            if let Some(region) = parse_bed_line(&line) {
                regions.push(region);
            }
            line.clear();
        }
    } else {
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            if let Some(region) = parse_bed_line(&line) {
                regions.push(region);
            }
        }
    }
    
    info!("Loaded {} valid regions", regions.len());
    Ok(regions)
}

/// Parse a single BED line into a ValidRegion
fn parse_bed_line(line: &str) -> Option<ValidRegion> {
    let line = line.trim();
    if line.starts_with('#') || line.is_empty() {
        return None;
    }
    
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 3 {
        return None;
    }
    
    let chrom = fields[0].trim_start_matches("chr").to_string();
    let start: u64 = fields[1].parse().ok()?;
    let end: u64 = fields[2].parse().ok()?;
    
    // GC content may be in 4th column if pre-computed
    let gc_content = if fields.len() >= 4 {
        fields[3].parse().unwrap_or(0.0)
    } else {
        0.0
    };
    
    Some(ValidRegion {
        chrom,
        start,
        end,
        gc_content,
    })
}

/// Generate valid regions BED file by excluding problematic and gap regions
/// 
/// # Arguments
/// * `reference_path` - Path to reference FASTA
/// * `exclude_regions_path` - Path to exclude regions BED (e.g., ENCODE list)
/// * `output_path` - Path to write valid_regions.bed
/// * `bin_size` - Size of each bin (default: 100000)
/// 
/// # Returns
/// * Number of valid regions generated
#[pyfunction]
#[pyo3(signature = (reference_path, exclude_regions_path, output_path, bin_size=100000))]
pub fn generate_valid_regions(
    reference_path: &str,
    exclude_regions_path: &str,
    output_path: &str,
    bin_size: u64,
) -> PyResult<usize> {
    info!("Generating valid regions from reference: {}", reference_path);
    info!("Exclude regions: {}", exclude_regions_path);
    info!("Bin size: {} bp", bin_size);
    
    // Load reference index to get chromosome lengths
    let faidx = faidx::Reader::from_path(reference_path)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to open reference FASTA: {}. Make sure .fai exists.", e)
        ))?;
    
    // Load exclude regions
    let exclude_regions = load_exclude_regions(Path::new(exclude_regions_path))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to load exclude regions: {}", e)
        ))?;
    
    // Valid chromosomes (1-22, X, Y) - try with and without chr prefix
    let valid_chroms: Vec<String> = (1..=22)
        .map(|i| i.to_string())
        .chain(vec!["X".to_string(), "Y".to_string()])
        .collect();
    
    let mut regions: Vec<ValidRegion> = Vec::new();
    
    // Get number of sequences in the reference
    let n_seqs = faidx.n_seqs();
    info!("Reference has {} sequences", n_seqs);
    
    // Build a map of chromosome name -> length from reference
    let mut chrom_lengths: HashMap<String, u64> = HashMap::new();
    for i in 0..n_seqs {
        if let Some(name) = faidx.seq_name(i as i32).ok() {
            // Normalize chromosome name (remove "chr" prefix if present)
            let normalized = name.trim_start_matches("chr").to_string();
            
            // Only include valid chromosomes
            if valid_chroms.contains(&normalized) {
                // Get sequence length
                let len = faidx.fetch_seq_len(&name);
                chrom_lengths.insert(normalized, len);
            }
        }
    }
    
    info!("Found {} valid chromosomes in reference", chrom_lengths.len());
    
    // Generate bins for each chromosome, excluding problematic regions
    for chrom in &valid_chroms {
        if let Some(&chrom_len) = chrom_lengths.get(chrom) {
            let excluded = exclude_regions.get(chrom).cloned().unwrap_or_default();
            
            // Generate bins of bin_size, excluding those that overlap with excluded regions
            let mut pos: u64 = 0;
            while pos + bin_size <= chrom_len {
                let bin_start = pos;
                let bin_end = pos + bin_size;
                
                // Check if this bin overlaps with any excluded region
                let overlaps_excluded = excluded.iter().any(|(excl_start, excl_end)| {
                    // Overlap if: bin_start < excl_end AND bin_end > excl_start
                    bin_start < *excl_end && bin_end > *excl_start
                });
                
                if !overlaps_excluded {
                    regions.push(ValidRegion {
                        chrom: chrom.clone(),
                        start: bin_start,
                        end: bin_end,
                        gc_content: 0.0,  // Will be computed later
                    });
                }
                
                pos += bin_size;
            }
        }
    }
    
    info!("Generated {} valid regions (excluded {} problematic bins)", 
          regions.len(), 
          chrom_lengths.values().map(|l| l / bin_size).sum::<u64>() as usize - regions.len());
    
    // Write output BED.gz (BGZF compressed)
    use noodles::bgzf;
    use std::io::Write;
    
    let output_file = File::create(output_path)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to create output file: {}", e)
        ))?;
    let mut writer = bgzf::io::Writer::new(output_file);
    
    // Write header
    writeln!(writer, "# Valid regions for GC correction estimation")
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    writeln!(writer, "# Excludes: problematic regions, assembly gaps, low mappability")
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    writeln!(writer, "#chrom\tstart\tend")
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    
    // Write regions
    for region in &regions {
        writeln!(writer, "{}\t{}\t{}", region.chrom, region.start, region.end)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    }
    
    // Finish BGZF stream
    writer.finish()
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to finish BGZF stream: {}", e)
        ))?;
    
    let region_count = regions.len();
    
    info!("Generated {} valid regions", region_count);
    Ok(region_count)
}

/// Load exclude regions from BED file
fn load_exclude_regions(path: &Path) -> Result<HashMap<String, Vec<(u64, u64)>>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open exclude regions file: {:?}", path))?;
    let reader = BufReader::new(file);
    
    let mut exclude_regions: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    
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
        
        exclude_regions.entry(chrom).or_insert_with(Vec::new).push((start, end));
    }
    
    // Sort intervals by start position
    for intervals in exclude_regions.values_mut() {
        intervals.sort_by_key(|i| i.0);
    }
    
    info!("Loaded exclude regions with {} chromosomes", exclude_regions.len());
    Ok(exclude_regions)
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
    
    // Progress tracking
    let total_regions = regions.len();
    let mut processed_regions = 0;
    
    // For each valid region, scan for GC content at each position for each length
    for region in &regions {
        // Convert chromosome name to the format expected by the FASTA
        // Try both with and without "chr" prefix
        let chrom_name = region.chrom.clone();
        let chrom_with_chr = format!("chr{}", region.chrom);
        
        // Try to fetch the sequence for this region
        let seq_result = faidx.fetch_seq(&chrom_name, region.start as usize, region.end as usize - 1);
        let seq = match seq_result {
            Ok(s) => s,
            Err(_) => {
                // Try with chr prefix
                match faidx.fetch_seq(&chrom_with_chr, region.start as usize, region.end as usize - 1) {
                    Ok(s) => s,
                    Err(_) => {
                        // Skip this region if we can't fetch it
                        continue;
                    }
                }
            }
        };
        
        let seq_len = seq.len();
        
        // For each length bin, sample positions and count GC
        for (bin_idx, &(min_len, max_len)) in LENGTH_BINS.iter().enumerate() {
            // Use the midpoint of each length bin for sampling
            let frag_len = ((min_len + max_len) / 2) as usize;
            
            // Sample every 100bp to speed up (still statistically representative)
            let step = 100;
            
            let mut pos = 0;
            while pos + frag_len <= seq_len {
                // Get the subsequence for this theoretical fragment
                let frag_seq = &seq[pos..pos + frag_len];
                
                // Count GC content
                let gc_count = frag_seq.iter()
                    .filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c')
                    .count();
                
                // Skip if too many Ns (uncertain bases)
                let n_count = frag_seq.iter()
                    .filter(|&&b| b == b'N' || b == b'n')
                    .count();
                
                if n_count * 10 > frag_len {
                    // More than 10% Ns, skip this position
                    pos += step;
                    continue;
                }
                
                // Calculate GC percent (0-100)
                let valid_bases = frag_len - n_count;
                let gc_percent = if valid_bases > 0 {
                    ((gc_count as f64 / valid_bases as f64) * 100.0).round() as usize
                } else {
                    50  // Default to 50% if no valid bases
                };
                
                let gc_idx = gc_percent.min(100);
                
                // Increment count
                counts[bin_idx][gc_idx] += 1;
                total_fragments += 1;
                
                pos += step;
            }
        }
        
        processed_regions += 1;
        
        // Log progress every 1000 regions
        if processed_regions % 1000 == 0 {
            info!("Processed {}/{} regions ({:.1}%), {} fragments so far", 
                  processed_regions, total_regions,
                  (processed_regions as f64 / total_regions as f64) * 100.0,
                  total_fragments);
        }
    }
    
    info!("Completed GC counting: {} total fragments across {} regions", 
          total_fragments, processed_regions);
    
    // Save as Parquet file (consistent with PON format)
    // Schema: length_bin_min (u32), length_bin_max (u32), gc_percent (u8), expected_count (u64)
    use arrow::array::{UInt32Array, UInt8Array, UInt64Array};
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::ArrowWriter;
    use std::sync::Arc;
    
    // Build arrays from counts
    let mut length_bin_min_vec: Vec<u32> = Vec::new();
    let mut length_bin_max_vec: Vec<u32> = Vec::new();
    let mut gc_percent_vec: Vec<u8> = Vec::new();
    let mut expected_count_vec: Vec<u64> = Vec::new();
    
    for (bin_idx, &(min_len, max_len)) in LENGTH_BINS.iter().enumerate() {
        for gc_percent in 0..=100u8 {
            let count = counts[bin_idx][gc_percent as usize];
            if count > 0 {  // Only write non-zero counts to save space
                length_bin_min_vec.push(min_len);
                length_bin_max_vec.push(max_len);
                gc_percent_vec.push(gc_percent);
                expected_count_vec.push(count);
            }
        }
    }
    
    // Create arrow arrays
    let length_bin_min_array = UInt32Array::from(length_bin_min_vec);
    let length_bin_max_array = UInt32Array::from(length_bin_max_vec);
    let gc_percent_array = UInt8Array::from(gc_percent_vec);
    let expected_count_array = UInt64Array::from(expected_count_vec);
    
    // Create schema
    let schema = Arc::new(Schema::new(vec![
        Field::new("length_bin_min", DataType::UInt32, false),
        Field::new("length_bin_max", DataType::UInt32, false),
        Field::new("gc_percent", DataType::UInt8, false),
        Field::new("expected_count", DataType::UInt64, false),
    ]));
    
    // Create record batch
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(length_bin_min_array),
            Arc::new(length_bin_max_array),
            Arc::new(gc_percent_array),
            Arc::new(expected_count_array),
        ],
    ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
        format!("Failed to create record batch: {}", e)
    ))?;
    
    // Write to Parquet file
    let output_file = std::fs::File::create(output_path)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to create output file: {}", e)
        ))?;
    
    let mut writer = ArrowWriter::try_new(output_file, schema, None)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to create Parquet writer: {}", e)
        ))?;
    
    writer.write(&batch).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
        format!("Failed to write Parquet batch: {}", e)
    ))?;
    
    writer.close().map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
        format!("Failed to close Parquet file: {}", e)
    ))?;
    
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
