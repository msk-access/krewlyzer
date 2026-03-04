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
use log::{debug, info};


/// Calculate Fragment Size Distribution (FSD)
/// 
/// # Arguments
/// * `bed_path` - Path to the input .bed.gz file (from motif)
/// * `arms_path` - Path to the chromosome arms BED file
/// * `output_path` - Path to the output TSV file
/// * `target_regions` - Optional path to target regions BED for panel mode
#[pyfunction]
#[pyo3(signature = (bed_path, arms_path, output_path, target_regions=None))]
pub fn calculate_fsd(
    bed_path: PathBuf,
    arms_path: PathBuf,
    output_path: PathBuf,
    target_regions: Option<PathBuf>,
) -> PyResult<()> {
    // 1. Parse Arms
    let regions = parse_regions_file(&arms_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    // 2. Setup Engine
    let mut chrom_map = ChromosomeMap::new();
    
    // Load target regions if provided (panel mode)
    let target_tree = if let Some(ref target_path) = target_regions {
        info!("Panel mode: loading target regions from {:?}", target_path);
        let tree = load_target_regions(target_path.as_path(), &mut chrom_map)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        Some(Arc::new(tree))
    } else {
        None
    };
    
    let consumer = FsdConsumer::new(regions, &mut chrom_map, None, target_tree);
    
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
    
    /// Write FSD output in TSV format (delegates to write_output_format).
    ///
    /// # Deprecated
    /// Use `write_output_format(output_path, "tsv", false)` for explicit format control.
    pub fn write_output(&self, output_path: &Path) -> Result<()> {
        self.write_output_format(output_path, "tsv", false)
    }

    /// Write FSD histogram to Parquet (output_utils::write_parquet_batch).
    ///
    /// Schema: region (Utf8) | 65-69 .. 395-399 (Float64 × 67) | total (Float64)
    ///
    /// # Errors
    /// Returns an error if Parquet serialization or I/O fails.
    fn write_parquet_output(&self, histograms: &[Vec<f64>], totals: &[f64], output_path: &Path) -> Result<()> {
        use arrow::array::{ArrayRef, Float64Array, StringArray};
        use arrow::datatypes::{DataType, Field, Schema};
        use crate::output_utils::write_parquet_batch;
        use std::sync::Arc;

        // Build column vectors
        let mut region_col: Vec<String> = Vec::with_capacity(self.regions.len());
        let mut bin_cols: Vec<Vec<f64>> = (0..67).map(|_| Vec::with_capacity(self.regions.len())).collect();
        let mut total_col: Vec<f64> = Vec::with_capacity(self.regions.len());

        for (i, region) in self.regions.iter().enumerate() {
            region_col.push(format!("{}:{}-{}", region.chrom, region.start, region.end));
            for b in 0..67 {
                bin_cols[b].push(histograms[i][b]);
            }
            total_col.push(totals[i]);
        }

        // Build Arrow schema
        let mut fields = vec![Field::new("region", DataType::Utf8, false)];
        for s in (65..400usize).step_by(5) {
            fields.push(Field::new(format!("{}-{}", s, s + 4), DataType::Float64, false));
        }
        fields.push(Field::new("total", DataType::Float64, false));
        let schema = Arc::new(Schema::new(fields));

        // Build Arrow arrays
        let mut arrays: Vec<ArrayRef> = Vec::with_capacity(69);
        arrays.push(Arc::new(StringArray::from(region_col)) as ArrayRef);
        for col in bin_cols {
            arrays.push(Arc::new(Float64Array::from(col)) as ArrayRef);
        }
        arrays.push(Arc::new(Float64Array::from(total_col)) as ArrayRef);

        write_parquet_batch(output_path, Arc::new(Schema::new(schema.fields().to_vec())), arrays)
    }

    /// Write FSD output in the requested format (TSV, Parquet, or both).
    ///
    /// Dispatches between TSV (`write_histogram`) and Parquet (`write_parquet_output`).
    pub fn write_output_format(&self, output_path: &Path, output_format: &str, compress: bool) -> Result<()> {
        use crate::output_utils::{should_write_tsv, should_write_parquet, tsv_path, validated_output_format};
        let fmt = validated_output_format(output_format);

        debug!("FSD: write_output_format(fmt={}, compress={}) → {:?}", fmt, compress, output_path);

        let total_off: f64 = self.totals_off.iter().sum();
        let total_on: f64 = self.totals_on.iter().sum();

        info!("FSD: {} arms, {:.0} off-target frags, {:.0} on-target frags",
            self.regions.len(), total_off, total_on);

        // Determine output paths.
        // IMPORTANT: Do not use with_extension("") then re-add .tsv — Rust's with_extension()
        // treats anything after the last dot as the extension, so 'test_sample.FSD.tsv'
        // → with_extension("") → 'test_sample.FSD' → with_extension("tsv") → 'test_sample.tsv'
        // (replaces 'FSD' as the extension). Instead, derive Parquet path by replacing .tsv,
        // and use tsv_path(output_path) directly for the TSV write.
        let parquet_path = if output_path.extension().is_some_and(|e| e == "tsv") {
            output_path.with_extension("parquet")
        } else {
            let mut p = output_path.to_path_buf();
            p.set_extension("parquet");
            p
        };
        // TSV path is output_path itself (optionally gzipped)
        let tsv_out = tsv_path(output_path, compress);

        // Off-target
        if should_write_tsv(fmt) {
            self.write_histogram(&self.histograms_off, &self.totals_off, &tsv_out)?;
        }
        if should_write_parquet(fmt) {
            self.write_parquet_output(&self.histograms_off, &self.totals_off, &parquet_path)?;
        }

        // On-target (panel mode)
        if self.target_tree.is_some() && total_on > 0.0 {
            // Build on-target paths by inserting '.ontarget' before the final extension
            let stem = output_path
                .file_stem().and_then(|s| s.to_str()).unwrap_or("output");
            let parent = output_path.parent().unwrap_or_else(|| std::path::Path::new("."));
            let on_tsv = parent.join(format!("{}.ontarget.tsv", stem));
            let on_parquet = parent.join(format!("{}.ontarget.parquet", stem));
            if should_write_tsv(fmt) {
                self.write_histogram(&self.histograms_on, &self.totals_on, &tsv_path(&on_tsv, compress))?;
            }
            if should_write_parquet(fmt) {
                self.write_parquet_output(&self.histograms_on, &self.totals_on, &on_parquet)?;
            }
            info!("FSD on-target: {:?}", on_tsv);
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
        if (65..400).contains(&len) {
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
                    } else if let Some(hist) = self.histograms_off.get_mut(idx) {
                        hist[bin_idx] += weight;
                        self.totals_off[idx] += weight;
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

// =============================================================================
// PON Log-Ratio Normalization (High-Performance Rust Implementation)
// =============================================================================

use crate::pon_model::FsdBaseline;

/// Apply PON log-ratio normalization to FSD TSV file.
/// 
/// Reads raw FSD counts, computes log2(sample / PON_expected) with pseudocount,
/// and writes enriched output with log-ratio columns.
/// 
/// # Arguments
/// * `fsd_input_path` - Path to raw FSD.tsv from calculate_fsd
/// * `pon_parquet_path` - Path to PON .parquet file containing fsd_baseline
/// * `output_path` - Output path for normalized FSD (optional, defaults to input)
/// 
/// # Returns
/// * Number of arms processed
/// 
/// # Performance
/// 10-50x faster than Python iterrows() implementation
#[pyfunction]
#[pyo3(signature = (fsd_input_path, pon_parquet_path, output_path=None))]
pub fn apply_pon_logratio(
    fsd_input_path: PathBuf,
    pon_parquet_path: PathBuf,
    output_path: Option<PathBuf>,
) -> PyResult<usize> {
    info!("FSD PON log-ratio: loading PON from {:?}", pon_parquet_path);
    
    // 1. Load FSD baseline from PON parquet
    let fsd_baseline = load_fsd_baseline_from_parquet(&pon_parquet_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to load FSD baseline: {}", e)
        ))?;
    
    if fsd_baseline.arms.is_empty() {
        info!("FSD PON: No FSD baseline found in PON, skipping normalization");
        return Ok(0);
    }
    
    info!("FSD PON: Loaded baseline for {} arms", fsd_baseline.arms.len());
    
    // 2. Read input FSD TSV
    let input_file = File::open(&fsd_input_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to open FSD input: {}", e)
        ))?;
    let reader = std::io::BufReader::new(input_file);
    
    // Parse header and data
    let mut lines: Vec<String> = Vec::new();
    for line in reader.lines() {
        lines.push(line.map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?);
    }
    
    if lines.is_empty() {
        return Ok(0);
    }
    
    // Parse header
    let header_line = &lines[0];
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
    
    info!("FSD PON: Found {} size bin columns", bin_indices.len());
    
    // 3. Process each arm row
    let output_path = output_path.unwrap_or(fsd_input_path.clone());
    let out_file = File::create(&output_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to create output file: {}", e)
        ))?;
    let mut writer = std::io::BufWriter::new(out_file);
    
    // Write extended header
    let mut new_headers = headers.iter().map(|s| s.to_string()).collect::<Vec<_>>();
    for (_, size) in &bin_indices {
        new_headers.push(format!("{}-{}_logR", size, size + 4));
    }
    new_headers.push("pon_stability".to_string());
    writeln!(writer, "{}", new_headers.join("\t"))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    
    let mut arms_processed = 0;
    let pseudocount = 1.0_f64;
    
    for line in lines.iter().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }
        
        let arm = fields[0];
        let mut new_row = fields.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        
        // Compute log-ratios for each bin
        let mut stds: Vec<f64> = Vec::new();
        
        for (col_idx, size) in &bin_indices {
            let sample_val: f64 = fields.get(*col_idx)
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.0);
            
            let log_ratio = if let Some((expected, std)) = fsd_baseline.get_stats(arm, *size) {
                if std > 0.0 {
                    stds.push(std);
                }
                if expected > 0.0 {
                    ((sample_val + pseudocount) / (expected + pseudocount)).log2()
                } else {
                    0.0
                }
            } else {
                0.0
            };
            
            new_row.push(format!("{:.6}", log_ratio));
        }
        
        // Compute PON stability score (inverse variance)
        let stability = if !stds.is_empty() {
            let avg_var = stds.iter().map(|s| s * s).sum::<f64>() / stds.len() as f64;
            1.0 / (avg_var + 0.01)
        } else {
            1.0
        };
        new_row.push(format!("{:.6}", stability));
        
        writeln!(writer, "{}", new_row.join("\t"))
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        
        arms_processed += 1;
    }
    
    info!("FSD PON: Processed {} arms with log-ratio normalization", arms_processed);
    
    Ok(arms_processed)
}

/// Load FsdBaseline from PON parquet file.
/// 
/// Parses the fsd_baseline table from the PON parquet.
fn load_fsd_baseline_from_parquet(path: &Path) -> Result<FsdBaseline> {
    use parquet::file::reader::{FileReader, SerializedFileReader};
    use parquet::record::RowAccessor;
    
    let file = File::open(path)?;
    let reader = SerializedFileReader::new(file)?;
    
    let mut arms: HashMap<String, crate::pon_model::FsdArmBaseline> = HashMap::new();
    
    for row_result in reader.get_row_iter(None)? {
        let row = row_result?;
        
        // Check if this row is from fsd_baseline table
        if let Ok(table) = row.get_string(
            row.get_column_iter().position(|(name, _)| name == "table").unwrap_or(0)
        ) {
            if table != "fsd_baseline" {
                continue;
            }
        } else {
            continue;
        }
        
        // Parse arm, size_bin, expected, std
        let arm = row.get_string(
            row.get_column_iter().position(|(name, _)| name == "arm").unwrap_or(0)
        ).map_or("".to_string(), |v| v.to_string());
        
        let size_bin = row.get_int(
            row.get_column_iter().position(|(name, _)| name == "size_bin").unwrap_or(0)
        ).unwrap_or(0);
        
        let expected = row.get_double(
            row.get_column_iter().position(|(name, _)| name == "expected").unwrap_or(0)
        ).unwrap_or(0.0);
        
        let std = row.get_double(
            row.get_column_iter().position(|(name, _)| name == "std").unwrap_or(0)
        ).unwrap_or(1.0);
        
        // Add to arm baseline
        let baseline = arms.entry(arm).or_insert_with(|| crate::pon_model::FsdArmBaseline {
            size_bins: Vec::new(),
            expected: Vec::new(),
            std: Vec::new(),
        });
        
        baseline.size_bins.push(size_bin);
        baseline.expected.push(expected);
        baseline.std.push(std);
    }
    
    // Sort each arm's bins
    for baseline in arms.values_mut() {
        let mut indices: Vec<usize> = (0..baseline.size_bins.len()).collect();
        indices.sort_by_key(|&i| baseline.size_bins[i]);
        
        let sorted_bins: Vec<i32> = indices.iter().map(|&i| baseline.size_bins[i]).collect();
        let sorted_exp: Vec<f64> = indices.iter().map(|&i| baseline.expected[i]).collect();
        let sorted_std: Vec<f64> = indices.iter().map(|&i| baseline.std[i]).collect();
        
        baseline.size_bins = sorted_bins;
        baseline.expected = sorted_exp;
        baseline.std = sorted_std;
    }
    
    info!("Loaded FSD baseline: {} arms", arms.len());
    
    Ok(FsdBaseline { arms })
}

