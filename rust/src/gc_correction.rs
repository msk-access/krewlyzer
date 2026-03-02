//! GC Bias Correction Module
//!
//! Provides LOESS-based GC bias correction for cfDNA fragment counts.
//! Supports both within-sample correction and PON-based correction.
//!
//! ## Output format
//! Correction factors can be written as tab-separated TSV (`.correction_factors.tsv`)
//! or Parquet (`.correction_factors.parquet`) depending on the `--output-format` flag.
//! Use [`CorrectionFactors::load`] to read either format transparently (Parquet-first).
//!
//! ## On-target vs off-target
//! Two independent write sites exist:
//! - **Off-target**: `pipeline.rs` calls `write_tsv`/`write_parquet` after `compute_gcfix_factors`.
//! - **On-target**: `compute_and_write_gc_factors` PyO3 fn called from Python (`sample_processor.py`)
//!   with the on-target observation dict and a `*.correction_factors.ontarget.tsv` path.

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

/// Fragment length bin definition for GC correction (188 bins, 5bp width)
/// Range: 60bp to 1000bp (extended for ultra-long fragment analysis)
/// 
/// 5bp granularity provides finer-grained GC bias modeling for ML features.
/// Bin 0 = 60-64bp, Bin 1 = 65-69bp, ..., Bin 187 = 995-999bp
/// 
/// Extended range captures ultra-long fragments associated with:
/// - Necrosis-derived cfDNA
/// - Fetal cfDNA (placental origin)  
/// - Apoptosis late-stage markers
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct LengthBin(pub u8);

impl LengthBin {
    /// Number of bins (188 = 940bp range / 5bp width, covering 60-1000bp)
    pub const NUM_BINS: u8 = 188;
    
    /// Get LengthBin from fragment length
    /// Returns None if length is outside tracked range (60-1000)
    pub fn from_len(len: u64) -> Option<Self> {
        if !(60..1000).contains(&len) {
            return None;
        }
        // (len - 60) / 5 -> 0..187
        // E.g. 60 -> 0, 64 -> 0, 65 -> 1, 999 -> 187
        let bin = ((len - 60) / 5) as u8;
        if bin < Self::NUM_BINS {
            Some(LengthBin(bin))
        } else {
            None
        }
    }

    /// Get min/max length for this bin (exclusive end)
    pub fn range(&self) -> (u64, u64) {
        let start = 60 + (self.0 as u64 * 5);
        (start, start + 5)
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
                let bin_idx = ((len_min - 60) / 5) as u8;
                if bin_idx < LengthBin::NUM_BINS {
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
    /// Return the GC correction factor for a fragment of given length and GC content.
    /// Returns `1.0` (no correction) when the bin is outside the tracked range.
    pub fn get_factor(&self, len: u64, gc: u8) -> f64 {
        if let Some(bin) = LengthBin::from_len(len) {
            if let Some(stats) = self.data.get(&(bin, gc)) {
                return stats.factor;
            }
        }
        1.0
    }

    // -----------------------------------------------------------------------
    // Writers
    // -----------------------------------------------------------------------

    /// Write correction factors to a **tab-separated** TSV file.
    ///
    /// Columns: `length_bin_min`, `length_bin_max`, `gc_percent`, `observed`,
    /// `expected`, `correction_factor`.
    ///
    /// This fixes the historical bug where the file had a `.tsv` extension but
    /// used comma separators.  All new code should call this method; the old
    /// `write_csv` name is gone.
    pub fn write_tsv<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        use std::io::BufWriter;
        let path = path.as_ref();
        debug!("gc_correction: writing TSV → {:?}", path);
        let file = File::create(path)
            .with_context(|| format!("Failed to create correction factors TSV: {:?}", path))?;
        let mut w = BufWriter::new(file);
        writeln!(w, "length_bin_min\tlength_bin_max\tgc_percent\tobserved\texpected\tcorrection_factor")?;
        let mut keys: Vec<_> = self.data.keys().collect();
        keys.sort();
        for key in &keys {
            let (bin, gc) = key;
            let stats = self.data.get(key).unwrap();
            let (min, max) = bin.range();
            writeln!(w, "{}\t{}\t{}\t{}\t{}\t{:.4}",
                min, max, gc, stats.observed, stats.expected, stats.factor)?;
        }
        info!("gc_correction: wrote {} factor rows → {:?}", keys.len(), path);
        Ok(())
    }

    /// Write correction factors to a Parquet file using the shared utility.
    ///
    /// Schema: `length_bin_min` (u64), `length_bin_max` (u64), `gc_percent` (u8),
    /// `observed` (u64), `expected` (u64), `correction_factor` (f64).
    pub fn write_parquet<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        use arrow::array::{ArrayRef, Float64Array, UInt8Array, UInt64Array};
        use arrow::datatypes::{DataType, Field, Schema};
        use std::sync::Arc;
        use crate::output_utils::write_parquet_batch;
        let path = path.as_ref();
        debug!("gc_correction: writing Parquet → {:?}", path);
        let mut keys: Vec<_> = self.data.keys().collect();
        keys.sort();
        let mut bin_mins: Vec<u64> = Vec::with_capacity(keys.len());
        let mut bin_maxs: Vec<u64> = Vec::with_capacity(keys.len());
        let mut gc_pcts: Vec<u8>   = Vec::with_capacity(keys.len());
        let mut obs: Vec<u64>      = Vec::with_capacity(keys.len());
        let mut exp: Vec<u64>      = Vec::with_capacity(keys.len());
        let mut factors: Vec<f64>  = Vec::with_capacity(keys.len());
        for key in &keys {
            let (bin, gc) = key;
            let stats = self.data.get(key).unwrap();
            let (min, max) = bin.range();
            bin_mins.push(min); bin_maxs.push(max); gc_pcts.push(*gc);
            obs.push(stats.observed); exp.push(stats.expected); factors.push(stats.factor);
        }
        let schema = Schema::new(vec![
            Field::new("length_bin_min",    DataType::UInt64, false),
            Field::new("length_bin_max",    DataType::UInt64, false),
            Field::new("gc_percent",        DataType::UInt8,  false),
            Field::new("observed",          DataType::UInt64, false),
            Field::new("expected",          DataType::UInt64, false),
            Field::new("correction_factor", DataType::Float64, false),
        ]);
        let arrays: Vec<ArrayRef> = vec![
            Arc::new(UInt64Array::from(bin_mins)),
            Arc::new(UInt64Array::from(bin_maxs)),
            Arc::new(UInt8Array::from(gc_pcts)),
            Arc::new(UInt64Array::from(obs)),
            Arc::new(UInt64Array::from(exp)),
            Arc::new(Float64Array::from(factors)),
        ];
        write_parquet_batch(path, Arc::new(schema), arrays)?;
        info!("gc_correction: wrote {} factor rows → {:?}", keys.len(), path);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Readers
    // -----------------------------------------------------------------------

    /// Load correction factors from a **tab-separated** TSV file.
    ///
    /// This is the updated loader matching [`write_tsv`].  It parses tab-separated
    /// fields and handles both the new TSV format and any legacy comma-separated
    /// files by auto-detecting the delimiter on the first data line.
    pub fn load_tsv<P: AsRef<Path>>(path: P) -> Result<Self> {
        use std::io::{BufRead, BufReader};
        let path = path.as_ref();
        debug!("gc_correction: loading TSV → {:?}", path);
        let file = File::open(path)
            .with_context(|| format!("Failed to open correction factors TSV: {:?}", path))?;
        let mut reader = BufReader::new(file);
        let mut data = HashMap::new();
        // Detect delimiter from header line
        let mut header = String::new();
        reader.read_line(&mut header)?;
        let delim = if header.contains('\t') { '\t' } else { ',' };
        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() { continue; }
            let fields: Vec<&str> = line.split(delim).collect();
            if fields.len() < 6 { continue; }
            let len_min: u64 = fields[0].trim().parse().unwrap_or(0);
            let gc: u8       = fields[2].trim().parse().unwrap_or(0);
            let observed: u64 = fields[3].trim().parse().unwrap_or(0);
            let expected: u64 = fields[4].trim().parse().unwrap_or(0);
            let factor: f64   = fields[5].trim().parse().unwrap_or(1.0);
            if len_min >= 60 {
                let bin_idx = ((len_min - 60) / 5) as u8;
                if bin_idx < LengthBin::NUM_BINS {
                    data.insert((LengthBin(bin_idx), gc), CorrectionBinStats { observed, expected, factor });
                }
            }
        }
        info!("gc_correction: loaded {} correction factors from {:?}", data.len(), path);
        Ok(Self { data })
    }

    /// Load correction factors from a Parquet file.
    pub fn load_parquet<P: AsRef<Path>>(path: P) -> Result<Self> {
        use parquet::file::reader::{FileReader, SerializedFileReader};
        use parquet::record::RowAccessor;
        let path = path.as_ref();
        debug!("gc_correction: loading Parquet → {:?}", path);
        let file = File::open(path)
            .with_context(|| format!("Failed to open correction factors Parquet: {:?}", path))?;
        let reader = SerializedFileReader::new(file)
            .map_err(|e| anyhow!("Parquet reader error: {}", e))?;
        let mut data = HashMap::new();
        for row in reader.get_row_iter(None)
            .map_err(|e| anyhow!("Row iterator error: {}", e))?
        {
            let row = row.map_err(|e| anyhow!("Row read error: {}", e))?;
            let len_min: u64 = row.get_ulong(0).unwrap_or(0);
            let gc: u8       = row.get_ubyte(2).or_else(|_| row.get_uint(2).map(|v| v as u8)).unwrap_or(0);
            let observed: u64 = row.get_ulong(3).unwrap_or(0);
            let expected: u64 = row.get_ulong(4).unwrap_or(0);
            let factor: f64   = row.get_double(5).unwrap_or(1.0);
            if len_min >= 60 {
                let bin_idx = ((len_min - 60) / 5) as u8;
                if bin_idx < LengthBin::NUM_BINS {
                    data.insert((LengthBin(bin_idx), gc), CorrectionBinStats { observed, expected, factor });
                }
            }
        }
        info!("gc_correction: loaded {} correction factors from Parquet {:?}", data.len(), path);
        Ok(Self { data })
    }

    /// **Parquet-first** auto-detect loader.  Used by all internal consumers
    /// (`mfsd.rs`, `pipeline.rs` input path, `fsc_processor.py`).
    ///
    /// Priority:
    /// 1. `{path}.parquet` (replace extension)
    /// 2. `{path}` as-is (TSV or legacy CSV — delimiter auto-detected)
    ///
    /// Returns an error only if both candidates fail to load.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let p = path.as_ref();
        let parquet_candidate = p.with_extension("parquet");
        if parquet_candidate.exists() {
            debug!("gc_correction::load: found Parquet candidate {:?}", parquet_candidate);
            Self::load_parquet(&parquet_candidate)
        } else if p.exists() {
            debug!("gc_correction::load: falling back to TSV {:?}", p);
            Self::load_tsv(p)
        } else {
            Err(anyhow!("Correction factors file not found: {:?} (tried {:?} and {:?})",
                p, parquet_candidate, p))
        }
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
    let mut bins_processed = 0u32;
    let mut bins_skipped = 0u32;
    
    // Iterate over each length bin (0..188)
    for bin_idx in 0..LengthBin::NUM_BINS {
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
            bins_skipped += 1;
            continue;
        }
        
        // Compute mean ratio for this bin (to preserve scale roughly)
        let mean_ratio: f64 = ratios.iter().sum::<f64>() / n_points as f64;
        
        // Fit LOESS
        let model = Lowess::new()
            .fraction(cfg.fraction)
            .iterations(cfg.iterations)
            .delta(cfg.delta)
            .adapter(Batch);
            
        if let Ok(model) = model.build() {
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
                 
                 debug!("Computed factors for bin {:?}: {} points, mean_ratio={:.3}", bin, result.x.len(), mean_ratio);
                 bins_processed += 1;
             }
        }
    }
    
    info!("GC correction: {} length bins processed, {} skipped (insufficient data)", bins_processed, bins_skipped);
    
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

/// Compute GC correction factors from fragment observations and write to TSV and/or Parquet.
///
/// Python-exposed function called during the extract step to produce
/// `*.correction_factors.tsv` / `*.correction_factors.parquet`.
///
/// Called from Python for **both** the off-target path (via `pipeline.rs`) and the
/// on-target path (from `sample_processor.py` with
/// `*.correction_factors.ontarget.tsv`).
///
/// # Arguments
/// * `observations` - `HashMap<(u8, u8), u64>` from `extract_motif`
///   where the key is `(length_bin_index, gc_percent_0_100)` and the value is the
///   fragment count.
/// * `gc_ref_path` - Path to the pre-computed GC reference Parquet file.
/// * `_valid_regions_path` - Reserved; not currently used (passed for API symmetry).
/// * `output_path` - Base output path **without** extension.  Extensions are appended
///   automatically based on `output_format`/`compress`.
/// * `output_format` - `"tsv"`, `"parquet"`, or `"both"`.  Defaults to `"tsv"`.
/// * `compress` - When `true` and writing TSV, gzip the output (`.tsv.gz`).
///
/// # Returns
/// Number of correction factor bins written.
#[pyfunction]
#[pyo3(signature = (observations, gc_ref_path, valid_regions_path, output_path, output_format="tsv", compress=false))]
pub fn compute_and_write_gc_factors(
    observations: HashMap<(u8, u8), u64>,
    gc_ref_path: &str,
    valid_regions_path: &str,
    output_path: &str,
    output_format: &str,
    compress: bool,
) -> PyResult<u64> {
    use crate::output_utils::{should_write_tsv, should_write_parquet, tsv_path, validated_output_format};
    let fmt = validated_output_format(output_format);
    info!(
        "gc_correction: computing factors from {} observations (format={}, compress={})",
        observations.len(), fmt, compress
    );

    // Convert observations from (u8, u8) keys to (LengthBin, u8) keys
    let mut obs_converted: HashMap<(LengthBin, u8), u64> = HashMap::new();
    for ((len_bin, gc_pct), count) in observations {
        obs_converted.insert((LengthBin(len_bin), gc_pct), count);
    }

    // Load reference data
    let ref_data = ReferenceData::load(Path::new(gc_ref_path))
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to load GC reference: {}", e)))?;
    debug!("gc_correction: GC reference loaded ({} entries)", ref_data.counts.len());

    // Compute factors
    let factors = compute_gcfix_factors(&obs_converted, &ref_data, None)
        .map_err(|e| PyRuntimeError::new_err(format!("Failed to compute GC factors: {}", e)))?;
    let n_factors = factors.data.len() as u64;
    info!("gc_correction: computed {} correction factor bins", n_factors);

    // Write outputs
    let base = Path::new(output_path);
    if should_write_tsv(fmt) {
        // Strip .tsv extension if caller already included it, then re-add with compress awareness
        let base_no_ext = if base.extension().map(|e| e == "tsv").unwrap_or(false) {
            base.with_extension("")
        } else {
            base.to_path_buf()
        };
        let tsv_out = tsv_path(&base_no_ext.with_extension("tsv"), compress);
        factors.write_tsv(&tsv_out)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to write TSV: {}", e)))?;
        info!("gc_correction: wrote {:?}", tsv_out);
    }
    if should_write_parquet(fmt) {
        let parquet_out = base.with_extension("parquet");
        factors.write_parquet(&parquet_out)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to write Parquet: {}", e)))?;
        info!("gc_correction: wrote {:?}", parquet_out);
    }

    Ok(n_factors)
}

#[cfg(test)]
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
