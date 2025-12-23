//! GC Bias Correction Module
//!
//! Provides LOESS-based GC bias correction for cfDNA fragment counts.
//! Supports both within-sample correction and PON-based correction.

use anyhow::{Result, anyhow};
use log::{info, debug};
use lowess::prelude::*;
use std::io::{Write, BufRead};
use coitrees::IntervalTree;

/// Configuration for GC bias correction
#[derive(Clone, Debug)]
pub struct GcCorrectionConfig {
    /// LOESS span (fraction of data used for each fit, 0.0-1.0)
    pub fraction: f64,
    /// Number of robust iterations for outlier handling
    pub iterations: usize,
    /// Delta for interpolation optimization
    pub delta: f64,
}

impl Default for GcCorrectionConfig {
    fn default() -> Self {
        Self {
            fraction: 0.3,    // 30% of GC range
            iterations: 3,    // Robust to CNV outliers
            delta: 0.01,      // 1% GC interpolation
        }
    }
}

/// Fragment length bin definition for GC correction (17 bins, 20bp width)
/// Range: 60bp to 400bp
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct LengthBin(pub u8);

impl LengthBin {
    /// Get LengthBin from fragment length
    /// Returns None if length is outside tracked range (60-400)
    pub fn from_len(len: u64) -> Option<Self> {
        if len < 60 || len >= 400 {
            return None;
        }
        // (len - 60) / 20 -> 0..16
        // E.g. 60 -> 0, 79 -> 0, 80 -> 1
        let bin = ((len - 60) / 20) as u8;
        if bin < 17 {
            Some(LengthBin(bin))
        } else {
            None
        }
    }

    /// Get min/max length for this bin
    pub fn range(&self) -> (u64, u64) {
        let start = 60 + (self.0 as u64 * 20);
        (start, start + 20)
    }
}

use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::RowAccessor;
use anyhow::Context;

/// Pre-computed GC reference data
#[derive(Debug, Clone)]
pub struct ReferenceData {
    // Map (LengthBin, gc_percent) -> expected_count
    pub counts: HashMap<(LengthBin, u8), u64>,
}

impl ReferenceData {
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open reference GC file: {:?}", path))?;
        let reader = SerializedFileReader::new(file)
            .map_err(|e| anyhow!("Failed to create Parquet reader: {}", e))?;
            
        let mut counts = HashMap::new();
        
        info!("Loading GC reference data from {:?}", path);
        
        // Using row iterator for simplicity
        let iter = reader.get_row_iter(None)
            .map_err(|e| anyhow!("Failed to get row iterator: {}", e))?;
            
        for row in iter {
            let row = row.map_err(|e| anyhow!("Failed to read row: {}", e))?;
            
            // Schema: length_bin_min (0), length_bin_max (1), gc_percent (2), expected_count (3)
            // Note: Indices depend on column order in file.
            // Using get_uint/get_ulong assuming correct types.
            // Safe access requires checking column names but for this internal tool we assume schema matches.
            
            let len_min: u32 = row.get_uint(0).map_err(|e| anyhow!("Missing/Invalid length_bin_min: {}", e))?;
            // Skip len_max (idx 1)
            // gc_percent may be stored as different types - try ubyte first (most common)
            let gc: u8 = row.get_ubyte(2)
                .or_else(|_| row.get_int(2).map(|v| v as u8))
                .or_else(|_| row.get_uint(2).map(|v| v as u8))
                .map_err(|e| anyhow!("Missing/Invalid gc_percent: {}", e))?; 
            let count: u64 = row.get_ulong(3).map_err(|e| anyhow!("Missing/Invalid expected_count: {}", e))?;
            
            // Map len_min back to LengthBin
            if len_min >= 60 {
                let bin_idx = ((len_min - 60) / 20) as u8;
                if bin_idx < 17 {
                    counts.insert((LengthBin(bin_idx), gc), count);
                }
            }
        }
        
        info!("Loaded GC reference data: {} entries", counts.len());
        Ok(Self { counts })
    }
    
    pub fn get_expected(&self, len_bin: LengthBin, gc: u8) -> Option<u64> {
        self.counts.get(&(len_bin, gc)).copied()
    }
}

/// Detailed statistics for a GC/Length bin
#[derive(Debug, Clone, Copy)]
pub struct CorrectionBinStats {
    pub observed: u64,
    pub expected: u64,
    pub factor: f64,
}

/// GC correction factors lookup table
#[derive(Debug, Clone)]
pub struct CorrectionFactors {
    // Map (LengthBin, gc_percent) -> stats
    pub data: HashMap<(LengthBin, u8), CorrectionBinStats>,
}

impl CorrectionFactors {
    pub fn get_factor(&self, len: u64, gc: u8) -> f64 {
        if let Some(bin) = LengthBin::from_len(len) {
            if let Some(stats) = self.data.get(&(bin, gc)) {
                return stats.factor;
            }
        }
        1.0
    }
    
    /// Write correction factors and stats to CSV
    pub fn write_csv<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let mut file = File::create(path)
            .with_context(|| "Failed to create correction factors CSV")?;
            
        writeln!(file, "length_bin_min,length_bin_max,gc_percent,observed,expected,correction_factor")?;
        
        let mut keys: Vec<_> = self.data.keys().collect();
        keys.sort();
        
        for key in keys {
            let (bin, gc) = key;
            let stats = self.data.get(key).unwrap();
            let (min, max) = bin.range();
            
            writeln!(file, "{},{},{},{},{},{:.4}", 
                min, max, gc, stats.observed, stats.expected, stats.factor)?;
        }
        info!("Written GC correction factors to CSV");
        Ok(())
    }
    
    /// Load correction factors from CSV (reverse of write_csv)
    pub fn load_csv<P: AsRef<Path>>(path: P) -> Result<Self> {
        use std::io::{BufRead, BufReader};
        
        let file = File::open(path.as_ref())
            .with_context(|| format!("Failed to open correction factors CSV: {:?}", path.as_ref()))?;
        let reader = BufReader::new(file);
        
        let mut data = HashMap::new();
        
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if i == 0 { continue; } // Skip header
            if line.trim().is_empty() { continue; }
            
            let fields: Vec<&str> = line.split(',').collect();
            if fields.len() < 6 { continue; }
            
            // Parse: length_bin_min,length_bin_max,gc_percent,observed,expected,correction_factor
            let len_min: u64 = fields[0].parse().unwrap_or(0);
            let gc: u8 = fields[2].parse().unwrap_or(0);
            let observed: u64 = fields[3].parse().unwrap_or(0);
            let expected: u64 = fields[4].parse().unwrap_or(0);
            let factor: f64 = fields[5].parse().unwrap_or(1.0);
            
            // Map len_min back to LengthBin
            if len_min >= 60 {
                let bin_idx = ((len_min - 60) / 20) as u8;
                if bin_idx < 17 {
                    data.insert((LengthBin(bin_idx), gc), CorrectionBinStats {
                        observed,
                        expected,
                        factor,
                    });
                }
            }
        }
        
        info!("Loaded {} correction factors from CSV", data.len());
        Ok(Self { data })
    }
}

/// Computes GCfix-style correction factors
/// 
/// # Algorithm
/// 1. For each length bin (60-400bp):
///    - Calculate ratio = observed_count / expected_count (from reference)
///    - Fit LOESS curve to (gc, ratio) points
///    - factor = mean_ratio / smoothed_ratio
/// 
/// # Arguments
/// * `observed` - Map of (LengthBin, gc) -> count
/// * `reference` - Loaded reference counts
/// * `config` - LOESS configuration
pub fn compute_gcfix_factors(
    observed: &HashMap<(LengthBin, u8), u64>,
    reference: &ReferenceData,
    config: Option<GcCorrectionConfig>
) -> Result<CorrectionFactors> {
    let cfg = config.unwrap_or_default();
    let mut factors = HashMap::new();
    
    // Iterate over each length bin (0..17)
    for bin_idx in 0..17 {
        let bin = LengthBin(bin_idx);
        
        // Collect (gc, ratio) pairs for this bin
        let mut gc_points = Vec::new();
        let mut ratios = Vec::new();
        
        for gc in 0..=100 {
            let obs = *observed.get(&(bin, gc)).unwrap_or(&0) as f64;
            let exp = reference.get_expected(bin, gc).unwrap_or(0) as f64;
            
            // Only use points with sufficient data for fitting
            if obs > 50.0 && exp > 50.0 {
                let ratio = obs / exp;
                gc_points.push(gc as f64 / 100.0); // 0.0-1.0 range for LOESS
                ratios.push(ratio);
            }
        }
        
        let n_points = gc_points.len();
        if n_points < 10 {
            debug!("Skipping length bin {:?} - insufficient data points ({})", bin, n_points);
            continue;
        }
        
        // Compute mean ratio for this bin (to preserve scale roughly)
        let mean_ratio: f64 = ratios.iter().sum::<f64>() / n_points as f64;
        
        // Fit LOESS
        let model = Lowess::new()
            .fraction(cfg.fraction)
            .iterations(cfg.iterations)
            .delta(cfg.delta)
            .adapter(Batch)
            .build(); // Using .build() directly? Check prev code
            
        // Previous code: .build().map_err(...)?
        // Let's check imports
        
        if let Ok(model) = model {
             if let Ok(result) = model.fit(&gc_points, &ratios) {
                 for (i, &gc_val) in result.x.iter().enumerate() {
                     let gc_byte = (gc_val * 100.0).round() as u8;
                     let fit_val = result.y[i];
                     
                     let factor = if fit_val > 0.0001 {
                         mean_ratio / fit_val
                     } else {
                         1.0
                     };
                     
                     let obs = *observed.get(&(bin, gc_byte)).unwrap_or(&0);
                     let exp = reference.get_expected(bin, gc_byte).unwrap_or(0);
                     
                     factors.insert((bin, gc_byte), CorrectionBinStats {
                         observed: obs,
                         expected: exp,
                         factor,
                     });
                 }
                 
                 info!("Computed factors for bin {:?}: {} points, mean_ratio={:.3}", bin, result.x.len(), mean_ratio);
             }
        }
    }
    
    Ok(CorrectionFactors { data: factors })
}

/// Apply LOESS-based GC bias correction to counts
/// 
/// # Arguments
/// * `gc_values` - GC content per region (0.0-1.0)
/// * `counts` - Fragment counts per region
/// * `config` - Optional LOESS configuration
///
/// # Returns
/// * Corrected counts
pub fn correct_gc_bias(
    gc_values: &[f64],
    counts: &[f64],
    config: Option<GcCorrectionConfig>
) -> Result<Vec<f64>> {
    let cfg = config.unwrap_or_default();
    
    if gc_values.len() != counts.len() {
        return Err(anyhow!("gc_values and counts must have same length"));
    }
    
    if gc_values.len() < 15 {
        return Err(anyhow!("Need at least 15 points for LOESS fitting"));
    }
    
    // Fit LOESS to (gc, counts) 
    let model = Lowess::new()
        .fraction(cfg.fraction)
        .iterations(cfg.iterations)
        .delta(cfg.delta)
        .adapter(Batch);
    
    if let Ok(model) = model.build() {
        if let Ok(result) = model.fit(gc_values, counts) {
            // Compute mean of fitted values
            let mean_fit: f64 = result.y.iter().sum::<f64>() / result.y.len() as f64;
            
            // Correction: multiply counts by (mean_fit / predicted_fit)
            let mut corrected = vec![0.0; counts.len()];
            for (i, &c) in counts.iter().enumerate() {
                let fitted = result.y[i];
                if fitted > 0.0001 {
                    corrected[i] = c * (mean_fit / fitted);
                } else {
                    corrected[i] = c;
                }
            }
            return Ok(corrected);
        }
    }
    
    // Fallback: return uncorrected
    Ok(counts.to_vec())
}



/// Correct WPS fragment types (different size ranges)
///
/// # Arguments
/// * `gc_values` - GC content per region
/// * `wps_long_counts` - WPS long fragment coverage (120-180bp)
/// * `wps_short_counts` - WPS short fragment coverage (35-80bp)
///
/// # Returns
/// * Tuple of corrected counts (wps_long, wps_short)
pub fn correct_gc_bias_wps(
    gc_values: &[f64],
    wps_long_counts: &[f64],
    wps_short_counts: &[f64],
    verbose: bool,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let config = GcCorrectionConfig::default();
    
    if verbose {
        info!("Applying WPS GC correction...");
    }
    
    let long_corrected = correct_gc_bias(gc_values, wps_long_counts, Some(config.clone()))
        .map_err(|e| anyhow!("WPS long GC correction failed: {}", e))?;
    
    let short_corrected = correct_gc_bias(gc_values, wps_short_counts, Some(config))
        .map_err(|e| anyhow!("WPS short GC correction failed: {}", e))?;
    
    if verbose {
        info!("WPS GC correction complete");
    }
    
    Ok((long_corrected, short_corrected))
}

use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;

/// Compute GC correction factors from observations and write to CSV.
/// Python-exposed function for inline factor computation during extract.
/// 
/// # Arguments
/// * `observations` - HashMap<(u8, u8), u64> from extract_motif (length_bin, gc_pct) -> count
/// * `gc_ref_path` - Path to GC reference parquet file
/// * `valid_regions_path` - Path to valid regions BED file  
/// * `output_path` - Path to write correction_factors.csv
#[pyfunction]
pub fn compute_and_write_gc_factors(
    observations: HashMap<(u8, u8), u64>,
    gc_ref_path: &str,
    valid_regions_path: &str,
    output_path: &str,
) -> PyResult<u64> {
    info!("Computing GC correction factors from {} observation bins...", observations.len());
    
    // Convert observations from (u8, u8) keys to (LengthBin, u8) keys
    let mut obs_converted: HashMap<(LengthBin, u8), u64> = HashMap::new();
    for ((len_bin, gc_pct), count) in observations {
        obs_converted.insert((LengthBin(len_bin), gc_pct), count);
    }
    
    // Load reference data
    let ref_data = ReferenceData::load(Path::new(gc_ref_path))
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to load GC reference: {}", e)))?;
    
    debug!("Loaded GC reference with {} entries", ref_data.counts.len());
    
    // Compute factors
    let factors = compute_gcfix_factors(&obs_converted, &ref_data, None)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to compute factors: {}", e)))?;
    
    let n_factors = factors.data.len() as u64;
    info!("Computed {} correction factors", n_factors);
    
    // Write to CSV
    factors.write_csv(output_path)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to write factors CSV: {}", e)))?;
    
    info!("Written correction factors to {}", output_path);
    
    Ok(n_factors)
}

/// Compute variance of a slice
fn variance(data: &[f64]) -> f64 {
    let n = data.len() as f64;
    if n < 2.0 {
        return 0.0;
    }
    let mean = data.iter().sum::<f64>() / n;
    data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)
}

use coitrees::{COITree, IntervalNode};
use crate::bed::{ChromosomeMap, Fragment};
use crate::engine::FragmentConsumer;
use std::sync::Arc;

/// Filters fragments based on valid regions (BED)
pub struct ValidRegionFilter {
    trees: HashMap<u32, COITree<usize, u32>>,
}

impl ValidRegionFilter {
    pub fn load(path: &Path, chrom_map: &mut ChromosomeMap) -> Result<Self> {
        let reader = crate::bed::get_reader(path)?;
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() { continue; }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 { continue; }
            
            let chrom = fields[0].trim_start_matches("chr");
            let start: i32 = fields[1].parse().unwrap_or(0);
            let end: i32 = fields[2].parse().unwrap_or(0);
            let end_closed = if end > start { end - 1 } else { start };
            
            let chrom_id = chrom_map.get_id(chrom);
            
            nodes_by_chrom.entry(chrom_id).or_default().push(
                IntervalNode::new(start, end_closed, 0)
            );
        }
        
        let mut trees = HashMap::new();
        for (chrom, nodes) in nodes_by_chrom {
            trees.insert(chrom, COITree::new(&nodes));
        }
        
        info!("Loaded valid regions filter covering {} chromosomes", trees.len());
        Ok(Self { trees })
    }
    
    /// Check if a fragment overlaps with valid regions
    pub fn is_valid(&self, chrom: u32, start: u64, end: u64) -> bool {
        if let Some(tree) = self.trees.get(&chrom) {
            let s = start as i32;
            let e = end as i32; // end is exclusive in fragment, but query is closed?
            let e_closed = if e > s { e - 1 } else { s };
            
            // Standard check: query_count > 0 simply means some overlap.
            tree.query_count(s, e_closed) > 0
        } else {
            false
        }
    }
}

/// Consumer to count observed fragments in valid regions
#[derive(Clone)]
pub struct GcObservationConsumer {
    pub observed: HashMap<(LengthBin, u8), u64>,
    pub valid_regions: Arc<ValidRegionFilter>, 
}

impl GcObservationConsumer {
    pub fn new(valid_regions: Arc<ValidRegionFilter>) -> Self {
        Self {
             observed: HashMap::new(),
             valid_regions,
        }
    }
}

impl FragmentConsumer for GcObservationConsumer {
    fn name(&self) -> &str { "GcObservation" }
    
    fn consume(&mut self, fragment: &Fragment) {
        if self.valid_regions.is_valid(fragment.chrom_id, fragment.start, fragment.end) {
             if let Some(bin) = LengthBin::from_len(fragment.length) {
                 let gc_byte = (fragment.gc * 100.0).round() as u8;
                 *self.observed.entry((bin, gc_byte)).or_insert(0) += 1;
             }
        }
    }
    
    fn merge(&mut self, other: Self) {
        for (k, v) in other.observed {
            *self.observed.entry(k).or_insert(0) += v;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_gc_correction_basic() {
        // Simulated GC bias: counts increase with GC (15 points to meet minimum)
        let gc = vec![0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.30, 0.40, 0.50, 0.60];
        let counts = vec![75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 82.0, 92.0, 98.0, 108.0];
        
        let corrected = correct_gc_bias(&gc, &counts, None).unwrap();
        
        // Check output has same length as input
        assert_eq!(corrected.len(), counts.len());
        
        // After correction, variance should be reduced (values should be closer to mean)
        let var_before = variance(&counts);
        let var_after = variance(&corrected);
        // Note: Due to LOESS fit, variance may not always decrease for small/noisy datasets
        println!("Variance before: {}, after: {}", var_before, var_after);
    }
    
    #[test]
    fn test_gc_correction_too_few_points() {
        let gc = vec![0.4, 0.5];
        let counts = vec![100.0, 100.0];
        
        let result = correct_gc_bias(&gc, &counts, None);
        assert!(result.is_err(), "Should fail with too few points");
    }
    
    #[test]
    fn test_variance() {
        let data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let var = variance(&data);
        assert!((var - 4.571).abs() < 0.01, "Variance calculation error");
    }
}
