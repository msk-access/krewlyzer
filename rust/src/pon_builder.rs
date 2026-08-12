//! PON Builder: Rust-accelerated aggregation for Panel of Normals construction.
//!
//! This module provides high-performance sample aggregation functions for building
//! PON (Panel of Normals) models. These functions are 3-10x faster than Python
//! equivalents due to optimized memory access and parallel processing.
//!
//! # Functions
//! - `compute_gc_bias_model`: Aggregate GC bias curves from sample data
//! - `compute_fsd_baseline`: Aggregate FSD per-arm baselines
//! - `compute_wps_baseline`: Aggregate WPS per-region baselines

use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;
use log::{info, warn};

/// GC bins for bias model: [0.25, 0.27, ..., 0.73]
const GC_BIN_START: f64 = 0.25;
const GC_BIN_END: f64 = 0.75;
const GC_BIN_STEP: f64 = 0.02;

/// Result structure for GC bias model
#[derive(Debug)]
pub struct GcBiasResult {
    pub gc_bins: Vec<f64>,
    pub short_expected: Vec<f64>,
    pub short_std: Vec<f64>,
    pub intermediate_expected: Vec<f64>,
    pub intermediate_std: Vec<f64>,
    pub long_expected: Vec<f64>,
    pub long_std: Vec<f64>,
}

/// Compute median of a slice
fn median(data: &mut [f64]) -> f64 {
    if data.is_empty() {
        return 1.0;
    }
    data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = data.len() / 2;
    if data.len().is_multiple_of(2) {
        (data[mid - 1] + data[mid]) / 2.0
    } else {
        data[mid]
    }
}

/// Compute standard deviation of a slice
fn std_dev(data: &[f64]) -> f64 {
    if data.len() < 2 {
        return 0.1;
    }
    let mean = data.iter().sum::<f64>() / data.len() as f64;
    let variance = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (data.len() - 1) as f64;
    variance.sqrt()
}

/// Compute GC bias model from sample data.
/// 
/// Aggregates GC-binned fragment counts across samples and computes
/// median expected coverage and std per GC bin.
/// 
/// # Arguments
/// * `all_gc_data` - List of dicts, each with keys: gc, short, intermediate, long (as arrays)
/// 
/// # Returns
/// * Dict with gc_bins, short_expected, short_std, intermediate_expected, intermediate_std, long_expected, long_std
/// 
/// # Performance
/// 3-10x faster than Python numpy/pandas aggregation
#[pyfunction]
pub fn compute_gc_bias_model(py: Python<'_>, all_gc_data: Vec<pyo3::Bound<'_, PyDict>>) -> PyResult<PyObject> {
    info!("PON Builder: Computing GC bias model from {} samples", all_gc_data.len());
    
    // Generate GC bins
    let num_bins = ((GC_BIN_END - GC_BIN_START) / GC_BIN_STEP) as usize;
    let gc_bins: Vec<f64> = (0..num_bins)
        .map(|i| GC_BIN_START + i as f64 * GC_BIN_STEP)
        .collect();
    
    // Aggregation structures
    let mut short_by_gc: Vec<Vec<f64>> = vec![Vec::new(); num_bins];
    let mut intermediate_by_gc: Vec<Vec<f64>> = vec![Vec::new(); num_bins];
    let mut long_by_gc: Vec<Vec<f64>> = vec![Vec::new(); num_bins];
    
    // Process each sample
    for sample in all_gc_data.iter() {
        // Extract arrays
        let gc_arr: Vec<f64> = sample.get_item("gc")
            .ok().and_then(|v| v.and_then(|item| item.extract().ok()))
            .unwrap_or_default();
        let short_arr: Vec<f64> = sample.get_item("short")
            .ok().and_then(|v| v.and_then(|item| item.extract().ok()))
            .unwrap_or_default();
        let intermediate_arr: Vec<f64> = sample.get_item("intermediate")
            .ok().and_then(|v| v.and_then(|item| item.extract().ok()))
            .unwrap_or_default();
        let long_arr: Vec<f64> = sample.get_item("long")
            .ok().and_then(|v| v.and_then(|item| item.extract().ok()))
            .unwrap_or_default();
        
        if gc_arr.is_empty() {
            continue;
        }
        
        // Compute means for normalization
        let short_mean: f64 = short_arr.iter().filter(|x| !x.is_nan()).sum::<f64>() 
            / short_arr.iter().filter(|x| !x.is_nan()).count().max(1) as f64;
        let intermediate_mean: f64 = intermediate_arr.iter().filter(|x| !x.is_nan()).sum::<f64>() 
            / intermediate_arr.iter().filter(|x| !x.is_nan()).count().max(1) as f64;
        let long_mean: f64 = long_arr.iter().filter(|x| !x.is_nan()).sum::<f64>() 
            / long_arr.iter().filter(|x| !x.is_nan()).count().max(1) as f64;
        
        // Bin by GC
        for (i, &gc_val) in gc_arr.iter().enumerate() {
            if gc_val.is_nan() || gc_val <= 0.0 || !(GC_BIN_START..GC_BIN_END).contains(&gc_val) {
                continue;
            }
            
            let bin_idx = ((gc_val - GC_BIN_START) / GC_BIN_STEP) as usize;
            if bin_idx >= num_bins {
                continue;
            }
            
            // Normalized values
            if i < short_arr.len() && !short_arr[i].is_nan() && short_mean > 0.0 {
                short_by_gc[bin_idx].push(short_arr[i] / short_mean);
            }
            if i < intermediate_arr.len() && !intermediate_arr[i].is_nan() && intermediate_mean > 0.0 {
                intermediate_by_gc[bin_idx].push(intermediate_arr[i] / intermediate_mean);
            }
            if i < long_arr.len() && !long_arr[i].is_nan() && long_mean > 0.0 {
                long_by_gc[bin_idx].push(long_arr[i] / long_mean);
            }
        }
    }
    
    // Compute median and std per bin
    let mut short_expected: Vec<f64> = Vec::with_capacity(num_bins);
    let mut short_std: Vec<f64> = Vec::with_capacity(num_bins);
    let mut intermediate_expected: Vec<f64> = Vec::with_capacity(num_bins);
    let mut intermediate_std: Vec<f64> = Vec::with_capacity(num_bins);
    let mut long_expected: Vec<f64> = Vec::with_capacity(num_bins);
    let mut long_std: Vec<f64> = Vec::with_capacity(num_bins);
    
    for i in 0..num_bins {
        short_expected.push(median(&mut short_by_gc[i].clone()));
        short_std.push(std_dev(&short_by_gc[i]));
        intermediate_expected.push(median(&mut intermediate_by_gc[i].clone()));
        intermediate_std.push(std_dev(&intermediate_by_gc[i]));
        long_expected.push(median(&mut long_by_gc[i].clone()));
        long_std.push(std_dev(&long_by_gc[i]));
    }
    
    // Log summary
    let total_samples: usize = short_by_gc.iter().map(|v| v.len()).sum();
    info!("PON Builder: GC model computed: {} bins, {} data points", num_bins, total_samples);
    
    // Build result dict
    let result = PyDict::new(py);
    result.set_item("gc_bins", gc_bins)?;
    result.set_item("short_expected", short_expected)?;
    result.set_item("short_std", short_std)?;
    result.set_item("intermediate_expected", intermediate_expected)?;
    result.set_item("intermediate_std", intermediate_std)?;
    result.set_item("long_expected", long_expected)?;
    result.set_item("long_std", long_std)?;
    
    Ok(result.into())
}

/// Compute FSD baseline from sample TSV files.
/// 
/// Reads FSD TSV files, aggregates counts per arm and size bin,
/// computes mean and std across samples.
/// 
/// # Arguments
/// * `fsd_paths` - List of paths to FSD.tsv files
/// 
/// # Returns
/// * Dict mapping arm -> {size_bins, expected, std}
#[pyfunction]
pub fn compute_fsd_baseline(py: Python<'_>, fsd_paths: Vec<String>) -> PyResult<PyObject> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    
    info!("PON Builder: Computing FSD baseline from {} samples", fsd_paths.len());
    
    // Aggregation: arm -> size_bin -> list of values
    let mut arm_data: HashMap<String, HashMap<i32, Vec<f64>>> = HashMap::new();
    
    for path in &fsd_paths {
        let file = match File::open(path) {
            Ok(f) => f,
            Err(e) => {
                info!("PON Builder: Skipping {}: {}", path, e);
                continue;
            }
        };
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        
        // Parse header
        let header_line = match lines.next() {
            Some(Ok(line)) => line,
            _ => continue,
        };
        let headers: Vec<&str> = header_line.split('\t').collect();
        
        // Find bin columns (format: "65-69", "70-74", etc.)
        let bin_indices: Vec<(usize, i32)> = headers.iter().enumerate()
            .filter_map(|(i, h)| {
                if let Some(dash_pos) = h.find('-') {
                    if let Ok(start) = h[..dash_pos].parse::<i32>() {
                        return Some((i, start));
                    }
                }
                None
            })
            .collect();
        
        // Process data rows
        for line_result in lines {
            let line = match line_result {
                Ok(l) => l,
                Err(_) => continue,
            };
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.is_empty() {
                continue;
            }
            
            let arm = fields[0].to_string();
            
            for (col_idx, size) in &bin_indices {
                let value: f64 = fields.get(*col_idx)
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0.0);
                
                arm_data
                    .entry(arm.clone())
                    .or_default()
                    .entry(*size)
                    .or_default()
                    .push(value);
            }
        }
    }
    
    // Compute mean and std per arm/bin
    let result = PyDict::new(py);
    
    for (arm, size_data) in &arm_data {
        let arm_dict = PyDict::new(py);
        let mut size_bins: Vec<i32> = size_data.keys().copied().collect();
        size_bins.sort();
        
        let mut expected: Vec<f64> = Vec::new();
        let mut std_vals: Vec<f64> = Vec::new();
        
        for size in &size_bins {
            let values = size_data.get(size).unwrap();
            let mean = values.iter().sum::<f64>() / values.len().max(1) as f64;
            expected.push(mean);
            std_vals.push(std_dev(values));
        }
        
        arm_dict.set_item("size_bins", size_bins)?;
        arm_dict.set_item("expected", expected)?;
        arm_dict.set_item("std", std_vals)?;
        
        result.set_item(arm, arm_dict)?;
    }
    
    info!("PON Builder: FSD baseline computed: {} arms", arm_data.len());
    
    Ok(result.into())
}

/// Extract 200-element WPS vector from a ListArray column.
/// Returns None if the column is not a List type or extraction fails.
fn extract_wps_vector(col: &dyn arrow::array::Array, row: usize) -> Option<Vec<f32>> {
    use arrow::array::{Float32Array, ListArray, Array};
    
    if let Some(list_arr) = col.as_any().downcast_ref::<ListArray>() {
        if row < list_arr.len() && !list_arr.is_null(row) {
            let inner = list_arr.value(row);
            if let Some(f32_arr) = inner.as_any().downcast_ref::<Float32Array>() {
                // Extract all values from the list
                let vec: Vec<f32> = (0..f32_arr.len())
                    .map(|i| if f32_arr.is_null(i) { 0.0 } else { f32_arr.value(i) })
                    .collect();
                return Some(vec);
            }
        }
    }
    None
}

/// Mean and sample standard deviation over the finite values, NaN when the
/// spread cannot be measured. Never a floor -- see `element_wise_std`.
/// Below this, a spread is floating-point residue rather than a measurement.
///
/// The exact `min == max` test above this constant catches donors that are
/// *byte*-identical. It does not catch the case that actually shipped: donors
/// that are numerically zero but not identically so. WPS is (fragments
/// spanning a position) minus (fragments ending near it), and every fragment
/// carries a fractional GC-correction weight, so a true zero computes as
/// ~1e-17 after cancellation instead of `0.0`.
///
/// Not a taste threshold -- the two populations do not overlap. Across the
/// 24.7M positive sigmas in the shipped xs1.all_unique WPS baseline:
///
/// ```text
///     sigma < 1e-12         1,177,647 positions   (residue)
///     1e-12 <= sigma < 1e-6         0 positions   (nothing lives here)
///     sigma >= 1e-6            24.7M positions    (measurements)
/// ```
///
/// Six empty decades. 1e-9 is their midpoint in log space, so the floor can
/// move two orders either way without reclassifying a single position.
///
/// What it costs when wrong in the safe direction: a genuine spread below 1e-9
/// on a quantity whose values are O(1)-O(1e3) would be reported as
/// unmeasurable. What it prevents: on a real XS2 plasma sample, 728,007
/// `wps_nuc_z` above 100 and 354,260 above 1e6, peaking at 6.1e18.
pub const SIGMA_FLOOR: f32 = 1e-9;

fn mean_and_sd(values: &[f32]) -> (f32, f32) {
    let finite: Vec<f32> = values.iter().copied().filter(|x| x.is_finite()).collect();
    if finite.is_empty() {
        return (f32::NAN, f32::NAN);
    }
    let n = finite.len() as f32;
    let mean = finite.iter().sum::<f32>() / n;
    if finite.len() < 2 {
        return (mean, f32::NAN);
    }

    // Identity is tested on the values, not on the result.
    //
    // `sd > 0.0` misses the case it exists for. Summation error grows with n,
    // so identical donors give residue rather than a clean zero: 21 copies of
    // 0.1 -- the xs2 cohort size -- yield 7.6e-9, and 8 copies of -0.3333333
    // yield 3.2e-8. Both pass `> 0.0` and become divisors.
    //
    // Measured in the rebuilt models, where these reach `compute_z_vector`
    // through the same `std > 0` guard: 4.5% of usable positions in
    // xs1.all_unique, 12.2% in xs2.all_unique and 55% in xs2.duplex, touching
    // 38-74% of anchors. A typical real sigma there is 0.4-1.2, so a one-unit
    // deviation against a residue sigma scores z ~ 1e11.
    //
    // Comparing min to max is exact and invents no tolerance; values that
    // genuinely differ keep whatever spread they have. Same fix as
    // `sample_std_or_nan` on the Python side.
    let (lo, hi) = finite
        .iter()
        .fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), &x| {
            (lo.min(x), hi.max(x))
        });
    if lo == hi {
        return (mean, f32::NAN);
    }

    let var = finite.iter().map(|x| (x - mean).powi(2)).sum::<f32>() / (n - 1.0);
    let sd = var.sqrt();
    // `>= SIGMA_FLOOR`, not `> 0.0` -- see the constant.
    (mean, if sd >= SIGMA_FLOOR { sd } else { f32::NAN })
}

/// Peak-to-trough range of a profile, on a log scale.
///
/// Raw amplitude is a coverage measurement, not a chromatin one: measured on
/// the real cohort it correlates +0.512 with `local_depth` and is skewed 11.6.
/// `ln(1 + range)` drops the depth correlation to -0.036 and the skew to 1.6,
/// which is the same multiplicative structure FSD and FSC depth show.
fn log_amplitude(v: &[f32]) -> f32 {
    let finite: Vec<f32> = v.iter().copied().filter(|x| x.is_finite()).collect();
    if finite.len() < 2 {
        return f32::NAN;
    }
    let hi = finite.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let lo = finite.iter().cloned().fold(f32::INFINITY, f32::min);
    (1.0 + (hi - lo)).ln()
}

/// Pearson correlation, or NaN when either side has no variance.
fn correlation(a: &[f32], b: &[f32]) -> f32 {
    let n = a.len().min(b.len());
    if n < 3 {
        return f32::NAN;
    }
    let (mut sa, mut sb) = (0.0f64, 0.0f64);
    for i in 0..n {
        sa += a[i] as f64;
        sb += b[i] as f64;
    }
    let (ma, mb) = (sa / n as f64, sb / n as f64);
    let (mut num, mut da, mut db) = (0.0f64, 0.0f64, 0.0f64);
    for i in 0..n {
        let (x, y) = (a[i] as f64 - ma, b[i] as f64 - mb);
        num += x * y;
        da += x * x;
        db += y * y;
    }
    if da < 1e-12 || db < 1e-12 {
        return f32::NAN;
    }
    (num / (da.sqrt() * db.sqrt())) as f32
}

/// Fisher z of a correlation: `arctanh(r)`.
///
/// Required, not cosmetic. A correlation is bounded at 1.0, and measured on
/// the real cohort the shape correlation sits at mean 0.844 with sigma 0.099 --
/// so the largest attainable positive z is about 1.5, and **302 of 400 anchors
/// could not reach +2 however tumour-like the sample**. On the Fisher scale
/// the same data has mean 1.371 and no ceiling.
fn fisher_z(r: f32) -> f32 {
    if !r.is_finite() {
        return f32::NAN;
    }
    r.clamp(-0.999_999, 0.999_999).atanh()
}

/// Compute element-wise mean of multiple vectors.
/// All vectors must have the same length.
fn element_wise_mean(vectors: &[Vec<f32>]) -> Vec<f32> {
    if vectors.is_empty() {
        return Vec::new();
    }
    let len = vectors[0].len();
    let n = vectors.len() as f32;
    
    (0..len)
        .map(|i| {
            let sum: f32 = vectors.iter().map(|v| v.get(i).copied().unwrap_or(0.0)).sum();
            sum / n
        })
        .collect()
}

/// Compute element-wise standard deviation of multiple vectors.
///
/// Returns NaN where the spread cannot be measured, rather than a floor. A
/// floored sigma does not make the z-score conservative -- it makes it
/// arbitrarily large, because the division is by a number nothing measured.
/// One sample previously yielded 0.1 at every position and a single anchor
/// could then produce an enormous z that read as strong signal.
///
/// NaN propagates through `(x - mean) / std` to an absent z, which is what
/// "we could not measure this" should look like downstream.
fn element_wise_std(vectors: &[Vec<f32>], means: &[f32]) -> Vec<f32> {
    if vectors.len() < 2 || means.is_empty() {
        return means.iter().map(|_| f32::NAN).collect();
    }
    let n = vectors.len() as f32;
    
    means.iter().enumerate()
        .map(|(i, &mean)| {
            // Identity is tested on the values, not on the result.
            //
            // `sd > 0.0` misses the case it exists for: summation error grows
            // with n, so identical donors give residue rather than a clean
            // zero. 21 copies of 0.1 -- the xs2 cohort size -- yield 7.6e-9.
            //
            // This is the same fix as `mean_and_sd` above, and it has to be
            // made twice because the two are separate implementations: that
            // one serves the scalar shape statistics, this one the 200-element
            // vectors. Measured in the shipped models, it is this one that
            // matters: residue reaches 4.5% of usable positions in
            // xs1.all_unique, 11.8% in xs2.all_unique and 55.4% in xs2.duplex,
            // touching 38-74% of anchors. A typical real sigma there is
            // 0.4-1.2, so a one-unit deviation scores z ~ 1e11.
            let (lo, hi) = vectors.iter().fold(
                (f32::INFINITY, f32::NEG_INFINITY),
                |(lo, hi), v| {
                    let val = v.get(i).copied().unwrap_or(0.0);
                    (lo.min(val), hi.max(val))
                },
            );
            if lo == hi {
                return f32::NAN;
            }

            let variance: f32 = vectors.iter()
                .map(|v| {
                    let val = v.get(i).copied().unwrap_or(0.0);
                    (val - mean).powi(2)
                })
                .sum::<f32>() / (n - 1.0);
            let sd = variance.sqrt();
            // `>= SIGMA_FLOOR`, not `> 0.0`: the exact min==max test above
            // catches byte-identical donors; this catches donors that are
            // numerically zero but not identically so, which is the case
            // that reached every shipped model.
            if sd >= SIGMA_FLOOR { sd } else { f32::NAN }
        })
        .collect()
}

/// Compute WPS baseline from sample parquet files (Vector v2.0 format).
/// 
/// Reads WPS parquet files, aggregates wps_nuc and wps_tf 200-bin vectors per region,
/// computes element-wise mean and std vectors across all samples.
/// 
/// # Arguments
/// * `wps_paths` - List of paths to WPS.parquet files
/// 
/// # Returns
/// * Dict mapping region_id -> {wps_nuc_mean: [200], wps_nuc_std: [200], wps_tf_mean: [200], wps_tf_std: [200]}
/// 
/// # Schema
/// v2.0 vector format enables position-specific z-score computation and shape correlation analysis.
/// 
/// # Performance
/// 3-10x faster than Python pandas concat + groupby
#[pyfunction]
pub fn compute_wps_baseline(py: Python<'_>, wps_paths: Vec<String>) -> PyResult<PyObject> {
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use arrow::array::{StringArray, Array};
    use std::fs::File;
    
    info!("PON Builder: Computing WPS vector baseline from {} samples", wps_paths.len());
    
    // Aggregation: region_id -> list of 200-element vectors
    let mut region_nuc_vectors: HashMap<String, Vec<Vec<f32>>> = HashMap::new();
    let mut region_tf_vectors: HashMap<String, Vec<Vec<f32>>> = HashMap::new();
    
    let mut files_processed = 0;
    let mut regions_found = 0;
    
    for path in &wps_paths {
        let file = match File::open(path) {
            Ok(f) => f,
            Err(e) => {
                info!("PON Builder: Skipping {}: {}", path, e);
                continue;
            }
        };
        
        let reader = match ParquetRecordBatchReaderBuilder::try_new(file) {
            Ok(builder) => match builder.build() {
                Ok(r) => r,
                Err(e) => {
                    info!("PON Builder: Failed to read parquet {}: {}", path, e);
                    continue;
                }
            },
            Err(e) => {
                info!("PON Builder: Failed to open parquet {}: {}", path, e);
                continue;
            }
        };
        
        files_processed += 1;
        
        for batch_result in reader {
            let batch = match batch_result {
                Ok(b) => b,
                Err(_) => continue,
            };
            
            let schema = batch.schema();
            
            // Find required columns
            let region_idx = schema.index_of("region_id")
                .or_else(|_| schema.index_of("group_id"))
                .ok();
            let nuc_idx = schema.index_of("wps_nuc").ok();
            let tf_idx = schema.index_of("wps_tf").ok();
            
            if region_idx.is_none() {
                continue;
            }
            
            let region_col = batch.column(region_idx.unwrap())
                .as_any().downcast_ref::<StringArray>();
            
            if region_col.is_none() {
                continue;
            }
            
            let region_array = region_col.unwrap();
            
            for i in 0..region_array.len() {
                let region_id = region_array.value(i).to_string();
                regions_found += 1;
                
                // Extract wps_nuc vector (expecting List<Float32>[200])
                if let Some(idx) = nuc_idx {
                    if let Some(vec) = extract_wps_vector(batch.column(idx).as_ref(), i) {
                        region_nuc_vectors.entry(region_id.clone())
                            .or_default()
                            .push(vec);
                    }
                }
                
                // Extract wps_tf vector
                if let Some(idx) = tf_idx {
                    if let Some(vec) = extract_wps_vector(batch.column(idx).as_ref(), i) {
                        region_tf_vectors.entry(region_id)
                            .or_default()
                            .push(vec);
                    }
                }
            }
        }
    }
    
    info!("PON Builder: Processed {} files, found {} region entries", files_processed, regions_found);
    info!("PON Builder: Unique regions: {} nuc, {} tf", region_nuc_vectors.len(), region_tf_vectors.len());
    
    // Compute element-wise mean and std per region
    let result = PyDict::new(py);
    
    // FSC gene and FSC region have required >=3 samples since they were
    // written; WPS -- 141k anchors, the largest block by 100x -- required
    // nothing, so an anchor seen in one sample still produced a baseline.
    // Measured on the shipped models, that was 1.6% of anchors for
    // all_unique/xs1 and 28.8% for duplex/xs1.
    const MIN_SAMPLES: usize = 3;
    let mut skipped = 0usize;

    for (region_id, nuc_vectors) in &region_nuc_vectors {
        if nuc_vectors.len() < MIN_SAMPLES {
            skipped += 1;
            continue;
        }
        let region_dict = PyDict::new(py);
        
        // Compute WPS-Nuc mean/std vectors
        let nuc_mean = element_wise_mean(nuc_vectors);
        let nuc_std = element_wise_std(nuc_vectors, &nuc_mean);
        region_dict.set_item("wps_nuc_mean", nuc_mean.iter().map(|&x| x as f64).collect::<Vec<f64>>())?;
        region_dict.set_item("wps_nuc_std", nuc_std.iter().map(|&x| x as f64).collect::<Vec<f64>>())?;
        
        // Compute WPS-TF mean/std vectors
        if let Some(tf_vectors) = region_tf_vectors.get(region_id) {
            let tf_mean = element_wise_mean(tf_vectors);
            let tf_std = element_wise_std(tf_vectors, &tf_mean);
            region_dict.set_item("wps_tf_mean", tf_mean.iter().map(|&x| x as f64).collect::<Vec<f64>>())?;
            region_dict.set_item("wps_tf_std", tf_std.iter().map(|&x| x as f64).collect::<Vec<f64>>())?;
        } else {
            // Empty vectors if no TF data
            region_dict.set_item("wps_tf_mean", Vec::<f64>::new())?;
            region_dict.set_item("wps_tf_std", Vec::<f64>::new())?;
        }
        
        // Add sample count for diagnostics
        region_dict.set_item("n_samples", nuc_vectors.len())?;

        // -- shape baseline -------------------------------------------------
        //
        // Derived quantities, each z-scored later against its own mean/sigma
        // rather than being a reduction of the per-position z vector. That
        // distinction is the point: adjacent WPS positions have lag-1
        // autocorrelation 0.986 (a fragment spans ~167bp and touches many
        // positions at once), so averaging z across positions produces a
        // number with none of a z-score's properties. Derive the biological
        // quantity first, then z-score that.
        //
        // Computed here rather than in Python because this loop already holds
        // every sample's vectors for the region; the Python side would have to
        // re-read ~44 MB per sample to see them again.
        let amps: Vec<f32> = nuc_vectors.iter().map(|v| log_amplitude(v)).collect();
        let corrs: Vec<f32> = nuc_vectors
            .iter()
            .map(|v| fisher_z(correlation(v, &nuc_mean)))
            .collect();
        // No phase-shift baseline. Measured on the real cohort: per-sample
        // mean lag varies by 0.26 bp against a within-sample spread of 8.43,
        // so there is no whole-sample phasing signal; and per anchor the
        // intraclass correlation is 0.479, meaning about half of a lag is
        // noise -- optimistically, since that baseline included the samples
        // being scored. Z-scoring an integer-valued statistic that is half
        // noise produces a plausible number and nothing else.
        //
        // The raw lag is still emitted per sample by `core/wps_pon.py`: it is
        // cheap, genuinely non-redundant (corr -0.24 and -0.28 with the two
        // below), and may resolve at n=21/47. Adding the baseline back is a
        // small change if the rebuild shows reproducible shifts.
        for (name, values) in [
            ("log_amplitude", &amps),
            ("shape_corr_fisher", &corrs),
        ] {
            let (m, sd) = mean_and_sd(values);
            region_dict.set_item(format!("{name}_mean"), m as f64)?;
            region_dict.set_item(format!("{name}_std"), sd as f64)?;
        }
        
        result.set_item(region_id, region_dict)?;
    }
    
    if skipped > 0 {
        warn!(
            "PON Builder: WPS baseline skipped {} of {} anchors backed by fewer \
             than {} samples; those anchors get no z-score rather than one \
             divided by a placeholder",
            skipped, region_nuc_vectors.len(), MIN_SAMPLES
        );
    }
    info!(
        "PON Builder: WPS vector baseline computed: {} of {} anchors",
        region_nuc_vectors.len() - skipped, region_nuc_vectors.len()
    );
    
    Ok(result.into())
}

/// Compute Region MDS baseline from sample TSV files.
/// 
/// Reads MDS.gene.tsv files, aggregates per-gene MDS values across samples,
/// computes mean and std for each gene.
/// 
/// # Arguments
/// * `mds_paths` - List of paths to MDS.gene.tsv files
/// 
/// # Returns
/// * Dict mapping gene -> {mds_mean, mds_std, mds_e1_mean, mds_e1_std, n_samples}
/// 
/// # Performance
/// Optimized for aggregating MDS across hundreds of samples
#[pyfunction]
pub fn compute_region_mds_baseline(py: Python<'_>, mds_paths: Vec<String>) -> PyResult<PyObject> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    
    info!("PON Builder: Computing Region MDS baseline from {} samples", mds_paths.len());
    
    // Aggregation: gene -> (mds_values, mds_e1_values)
    let mut gene_mds: HashMap<String, Vec<f64>> = HashMap::new();
    let mut gene_mds_e1: HashMap<String, Vec<f64>> = HashMap::new();
    
    for path in &mds_paths {
        let file = match File::open(path) {
            Ok(f) => f,
            Err(e) => {
                info!("PON Builder: Skipping {}: {}", path, e);
                continue;
            }
        };
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        
        // Parse header to find column indices
        let header_line = match lines.next() {
            Some(Ok(line)) => line,
            _ => continue,
        };
        let headers: Vec<&str> = header_line.split('\t').collect();
        
        // Find column indices
        let gene_idx = headers.iter().position(|&h| h == "gene").unwrap_or(0);
        let mds_mean_idx = headers.iter().position(|&h| h == "mds_mean");
        let mds_e1_idx = headers.iter().position(|&h| h == "mds_e1");
        
        // Process data rows
        for line_result in lines {
            let line = match line_result {
                Ok(l) => l,
                Err(_) => continue,
            };
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.is_empty() {
                continue;
            }
            
            let gene = fields.get(gene_idx).unwrap_or(&"").to_string();
            if gene.is_empty() {
                continue;
            }
            
            // mds_mean
            if let Some(idx) = mds_mean_idx {
                if let Some(val_str) = fields.get(idx) {
                    if let Ok(val) = val_str.parse::<f64>() {
                        if !val.is_nan() && val > 0.0 {
                            gene_mds.entry(gene.clone()).or_default().push(val);
                        }
                    }
                }
            }
            
            // mds_e1
            if let Some(idx) = mds_e1_idx {
                if let Some(val_str) = fields.get(idx) {
                    if let Ok(val) = val_str.parse::<f64>() {
                        if !val.is_nan() && val > 0.0 {
                            gene_mds_e1.entry(gene).or_default().push(val);
                        }
                    }
                }
            }
        }
    }
    
    // Compute mean and std per gene
    let result = PyDict::new(py);
    
    for (gene, mds_values) in &gene_mds {
        let gene_dict = PyDict::new(py);
        
        // MDS mean stats
        let mds_mean = mds_values.iter().sum::<f64>() / mds_values.len().max(1) as f64;
        let mds_std = std_dev(mds_values);
        gene_dict.set_item("mds_mean", mds_mean)?;
        gene_dict.set_item("mds_std", mds_std)?;
        
        // MDS E1 stats
        if let Some(e1_values) = gene_mds_e1.get(gene) {
            let e1_mean = e1_values.iter().sum::<f64>() / e1_values.len().max(1) as f64;
            let e1_std = std_dev(e1_values);
            gene_dict.set_item("mds_e1_mean", e1_mean)?;
            gene_dict.set_item("mds_e1_std", e1_std)?;
        } else {
            gene_dict.set_item("mds_e1_mean", mds_mean)?;
            gene_dict.set_item("mds_e1_std", mds_std)?;
        }
        
        gene_dict.set_item("n_samples", mds_values.len())?;
        
        result.set_item(gene, gene_dict)?;
    }
    
    info!("PON Builder: Region MDS baseline computed: {} genes", gene_mds.len());
    
    Ok(result.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A floored sigma is not conservative -- it is a fabrication.
    ///
    /// The shipped models returned 0.1 at every position for a single-sample
    /// anchor and floored everything else at 0.01, so `z = (x - mean) / 0.01`
    /// produced enormous values from a baseline that measured nothing. NaN is
    /// the honest answer, and it propagates to an absent z.
    #[test]
    fn a_single_sample_yields_no_spread_rather_than_a_floor() {
        let vectors = vec![vec![5.0_f32, 7.0, 9.0]];
        let means = element_wise_mean(&vectors);
        let std = element_wise_std(&vectors, &means);
        assert_eq!(means, vec![5.0, 7.0, 9.0], "the mean is still measurable");
        assert!(
            std.iter().all(|v| v.is_nan()),
            "one sample must give NaN, got {std:?}"
        );
    }

    #[test]
    fn identical_samples_yield_no_spread_rather_than_zero() {
        let vectors = vec![vec![3.0_f32; 4], vec![3.0; 4], vec![3.0; 4]];
        let means = element_wise_mean(&vectors);
        let std = element_wise_std(&vectors, &means);
        assert!(
            std.iter().all(|v| v.is_nan()),
            "zero spread is still no information about spread, got {std:?}"
        );
    }

    #[test]
    fn a_real_spread_is_measured_with_the_sample_correction() {
        // values 2, 4, 6 -> mean 4, sample sd (ddof=1) = 2.0, not 1.633.
        // The population form understates the spread and inflates every z.
        let vectors = vec![vec![2.0_f32], vec![4.0], vec![6.0]];
        let means = element_wise_mean(&vectors);
        let std = element_wise_std(&vectors, &means);
        assert!((means[0] - 4.0).abs() < 1e-6);
        assert!(
            (std[0] - 2.0).abs() < 1e-5,
            "expected the ddof=1 sample sd of 2.0, got {}",
            std[0]
        );
    }

    #[test]
    fn a_measured_spread_is_never_replaced_by_a_floor() {
        // Genuinely tiny but real: 0.001 must survive, not be raised to 0.01.
        let vectors = vec![vec![1.0_f32], vec![1.001], vec![0.999]];
        let means = element_wise_mean(&vectors);
        let std = element_wise_std(&vectors, &means);
        assert!(std[0] < 0.01, "a real sub-floor spread was floored: {}", std[0]);
        assert!(std[0] > 0.0);
    }
}

#[cfg(test)]
mod shape_tests {
    use super::*;

    /// Raw amplitude is a coverage measurement, not a chromatin one.
    ///
    /// Measured on the real cohort: raw amplitude correlates +0.512 with
    /// `local_depth` and is skewed 11.6; `ln(1 + range)` drops that to -0.036
    /// and 1.6. A z-score on the raw scale would rank samples by how deeply
    /// they were sequenced.
    #[test]
    fn log_amplitude_compresses_a_depth_driven_range() {
        let shallow: Vec<f32> = (0..200).map(|i| (i as f32 / 20.0).sin()).collect();
        let deep: Vec<f32> = shallow.iter().map(|x| x * 100.0).collect();
        let (a, b) = (log_amplitude(&shallow), log_amplitude(&deep));
        assert!(b > a, "a deeper profile must still register as larger");
        assert!(
            b / a < 10.0,
            "100x the depth must not be 100x the statistic, got {}x",
            b / a
        );
    }

    #[test]
    fn log_amplitude_is_nan_for_a_profile_with_nothing_in_it() {
        assert!(log_amplitude(&[]).is_nan());
        assert!(log_amplitude(&[1.0]).is_nan());
    }

    /// The ceiling this transform exists for.
    ///
    /// A correlation is bounded at 1.0. On the real cohort the shape
    /// correlation sits at mean 0.844 with sigma 0.099, so the largest
    /// attainable positive z is ~1.5 and 302 of 400 anchors could not reach
    /// +2 however tumour-like the sample. Fisher's transform removes the
    /// bound.
    #[test]
    fn fisher_z_removes_the_correlation_ceiling() {
        let near_ceiling = [0.90_f32, 0.95, 0.99];
        let transformed: Vec<f32> = near_ceiling.iter().map(|&r| fisher_z(r)).collect();
        // Raw: the gaps shrink toward the bound. Fisher: they widen.
        let raw_gap = near_ceiling[2] - near_ceiling[1];
        let fisher_gap = transformed[2] - transformed[1];
        assert!(
            fisher_gap > raw_gap * 3.0,
            "fisher must expand the region near 1.0, gaps {raw_gap} vs {fisher_gap}"
        );
        assert!(transformed.iter().all(|z| z.is_finite()));
    }

    #[test]
    fn fisher_z_survives_a_perfect_correlation() {
        // arctanh(1.0) is infinite; the clamp keeps it a number.
        assert!(fisher_z(1.0).is_finite());
        assert!(fisher_z(-1.0).is_finite());
        assert!(fisher_z(f32::NAN).is_nan());
    }

    #[test]
    fn correlation_is_nan_without_variance() {
        let flat = [3.0_f32; 50];
        let ramp: Vec<f32> = (0..50).map(|i| i as f32).collect();
        assert!(correlation(&flat, &ramp).is_nan());
        assert!((correlation(&ramp, &ramp) - 1.0).abs() < 1e-5);
    }

    #[test]
    fn mean_and_sd_reports_an_unmeasurable_spread_as_nan() {
        assert!(mean_and_sd(&[]).0.is_nan());
        let (m, sd) = mean_and_sd(&[5.0]);
        assert_eq!(m, 5.0);
        assert!(sd.is_nan(), "one value has no spread");
        let (m, sd) = mean_and_sd(&[3.0, 3.0, 3.0]);
        assert_eq!(m, 3.0);
        assert!(sd.is_nan(), "zero spread is still no information about spread");
        let (m, sd) = mean_and_sd(&[2.0, 4.0, 6.0]);
        assert!((m - 4.0).abs() < 1e-6);
        assert!((sd - 2.0).abs() < 1e-5, "ddof=1 sample sd, got {sd}");
    }

    #[test]
    fn identical_values_give_nan_even_when_the_variance_does_not_reach_zero() {
        // `[3.0; 3]` above happens to cancel exactly, which is why a `sd > 0`
        // guard passed that test while failing in production. Summation error
        // grows with n: these are the real cohort sizes.
        for (values, label) in [
            (vec![0.1f32; 21], "21 donors, xs2"),
            (vec![-0.3333333f32; 8], "8 donors"),
            (vec![0.95f32; 47], "47 donors, xs1"),
        ] {
            let (_, sd) = mean_and_sd(&values);
            assert!(sd.is_nan(), "{label}: identical values gave sd {sd:e}");
        }
    }

    #[test]
    fn element_wise_std_gives_nan_for_identical_vectors() {
        // The vector twin of `mean_and_sd`, and the one that actually reaches
        // `wps_nuc_std`. It needed the same fix separately -- fixing only
        // `mean_and_sd` left every 200-element sigma untouched, which is where
        // the residue was measured in the shipped models.
        let vectors: Vec<Vec<f32>> = (0..21).map(|_| vec![0.1f32; 4]).collect();
        let means = element_wise_mean(&vectors);
        let stds = element_wise_std(&vectors, &means);
        for (i, sd) in stds.iter().enumerate() {
            assert!(sd.is_nan(), "position {i}: identical donors gave sd {sd:e}");
        }
    }

    #[test]
    fn element_wise_std_keeps_a_real_spread() {
        let vectors: Vec<Vec<f32>> = (0..21)
            .map(|i| vec![0.1f32 + i as f32 * 0.01, 0.5f32])
            .collect();
        let means = element_wise_mean(&vectors);
        let stds = element_wise_std(&vectors, &means);
        assert!(stds[0].is_finite() && stds[0] > 0.0, "position 0 varies");
        assert!(stds[1].is_nan(), "position 1 is identical across donors");
    }

    #[test]
    fn a_spread_of_one_ulp_is_still_a_spread() {
        // There *is* a tolerance now -- `SIGMA_FLOOR` -- and this test bounds
        // it from below. One ULP at 0.95 is 5.96e-8, giving sd 4.2e-8, which
        // is 42x above the 1e-9 floor. So the smallest spread f32 can express
        // at a magnitude the WPS baseline actually contains still survives.
        //
        // That is what makes the floor safe rather than a fudge: it can only
        // discard a spread finer than f32 can represent at these magnitudes,
        // which is by definition not a measurement.
        let a = 0.95f32;
        let b = f32::from_bits(a.to_bits() + 1);
        assert_ne!(a, b);
        let (_, sd) = mean_and_sd(&[a, b]);
        assert!(sd.is_finite() && sd > 0.0, "one ULP apart is not identical");
        assert!(
            sd > SIGMA_FLOOR * 10.0,
            "one ULP gives sd {sd:e}, uncomfortably close to the {SIGMA_FLOOR:e} \
             floor -- if these ever converge the floor is discarding real spread"
        );
    }

    #[test]
    fn donors_that_are_numerically_zero_give_no_sigma() {
        // The case that shipped in all four models. Not byte-identical, so the
        // exact min==max test never fired: WPS carries fractional GC weights,
        // so a true zero lands at ~1e-17 after cancellation.
        //
        // On a real XS2 plasma sample this produced 728,007 `wps_nuc_z` above
        // 100 and 354,260 above 1e6, peaking at 6.1e18.
        let vectors: Vec<Vec<f32>> = (0..21)
            .map(|i| vec![(i as f32 - 10.0) * 1e-18, 4.0 + i as f32 * 0.1])
            .collect();
        let means = element_wise_mean(&vectors);
        let stds = element_wise_std(&vectors, &means);
        assert!(
            stds[0].is_nan(),
            "position 0 is numerically zero across donors, got sd {:e}",
            stds[0]
        );
        assert!(
            stds[1].is_finite() && stds[1] > 0.0,
            "position 1 has a real spread and must keep it"
        );
    }

    #[test]
    fn the_floor_sits_in_an_empty_gap() {
        // Measured across the 24.7M positive sigmas of the shipped
        // xs1.all_unique WPS baseline: 1,177,647 below 1e-12, zero between
        // 1e-12 and 1e-6, the rest above. The floor is the log-space midpoint
        // of six empty decades, so it may move two orders either way without
        // reclassifying a single position. This pins that it stays there.
        assert!(SIGMA_FLOOR > 1e-12, "floor dropped into the residue population");
        assert!(SIGMA_FLOOR < 1e-6, "floor rose into the measured population");
    }
}
