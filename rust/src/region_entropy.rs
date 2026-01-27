//! Region Entropy Module for krewlyzer
//!
//! Calculates Shannon entropy of fragment size distributions aggregated by genomic regions.
//! Used for TFBS (transcription factor binding site) and ATAC (cancer peak) analysis.
//!
//! # Output Format
//! - `{sample}.{type}.tsv`: Entropy scores per region label
//!   Columns: label, count, mean_size, entropy

use pyo3::prelude::*;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, Write};
use std::collections::HashMap;
use std::sync::Arc;
use anyhow::{Result, Context};
use coitrees::{COITree, IntervalNode, IntervalTree};

use crate::bed::{Fragment, ChromosomeMap};
use crate::engine::{FragmentConsumer, FragmentAnalyzer};
use crate::gc_correction::CorrectionFactors;

// Maximum fragment length to track for entropy (covers cfDNA signal 65-400bp, extended to 1000bp)
const MAX_LEN: usize = 1000;

/// Per-label entropy statistics
#[derive(Clone)]
pub struct EntropyStats {
    /// Histogram: index = fragment length (0..1000), value = weighted count
    len_counts: Vec<f64>,
    /// Total weighted fragment count
    total_count: f64,
}

impl Default for EntropyStats {
    fn default() -> Self {
        Self::new()
    }
}

impl EntropyStats {
    pub fn new() -> Self {
        Self {
            len_counts: vec![0.0; MAX_LEN + 1],
            total_count: 0.0,
        }
    }

    /// Merge another EntropyStats into this one
    pub fn merge(&mut self, other: &Self) {
        for (i, count) in other.len_counts.iter().enumerate() {
            self.len_counts[i] += count;
        }
        self.total_count += other.total_count;
    }

    /// Calculate Shannon entropy of the fragment size distribution
    /// 
    /// H = -Σ p(x) * log2(p(x))
    /// 
    /// Higher entropy = more diverse size distribution
    /// Lower entropy = more uniform size distribution
    pub fn calculate_entropy(&self) -> f64 {
        if self.total_count <= 0.0 {
            return 0.0;
        }

        let mut entropy = 0.0;
        for &count in &self.len_counts {
            if count > 0.0 {
                let p = count / self.total_count;
                entropy -= p * p.log2();
            }
        }
        entropy
    }

    /// Calculate mean fragment size
    pub fn mean_size(&self) -> f64 {
        if self.total_count <= 0.0 {
            return 0.0;
        }
        let sum_len: f64 = self
            .len_counts
            .iter()
            .enumerate()
            .map(|(len, count)| len as f64 * count)
            .sum();
        sum_len / self.total_count
    }
}

/// Region Entropy Consumer
/// 
/// Processes fragments and computes size entropy per region label.
/// Implements FragmentConsumer trait for parallel processing.
#[derive(Clone)]
pub struct RegionEntropyConsumer {
    /// Interval trees: ChromID -> IntervalTree -> label_id
    trees: Arc<HashMap<u32, COITree<usize, u32>>>,
    /// Label names indexed by label_id
    labels: Arc<Vec<String>>,
    /// GC correction factors (optional)
    factors: Option<Arc<CorrectionFactors>>,
    /// Per-label statistics
    stats: Vec<EntropyStats>,
}

impl RegionEntropyConsumer {
    /// Create a new RegionEntropyConsumer
    ///
    /// # Arguments
    /// * `region_path` - Path to BED file with regions (chrom, start, end, label)
    /// * `chrom_map` - Chromosome ID mapping
    /// * `factors` - Optional GC correction factors
    pub fn new(
        region_path: &Path,
        chrom_map: &mut ChromosomeMap,
        factors: Option<Arc<CorrectionFactors>>,
    ) -> Result<Self> {
        let reader = crate::bed::get_reader(region_path)?;
        let mut label_map: HashMap<String, usize> = HashMap::new();
        let mut labels: Vec<String> = Vec::new();
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 4 {
                continue;
            }

            let chrom = cols[0].trim_start_matches("chr");
            let start: i32 = cols[1].parse().unwrap_or(0);
            let end: i32 = cols[2].parse().unwrap_or(0);
            let label = cols[3].to_string();

            // Get or create label ID
            let label_id = *label_map.entry(label.clone()).or_insert_with(|| {
                labels.push(label);
                labels.len() - 1
            });

            let chrom_id = chrom_map.get_id(chrom);
            nodes_by_chrom
                .entry(chrom_id)
                .or_default()
                .push(IntervalNode::new(start, end - 1, label_id)); // closed interval
        }

        // Build interval trees
        let mut trees = HashMap::new();
        for (k, v) in nodes_by_chrom {
            trees.insert(k, COITree::new(&v));
        }

        log::info!(
            "RegionEntropy: Loaded {} labels from {} regions",
            labels.len(),
            trees.values().map(|t| t.len()).sum::<usize>()
        );

        Ok(Self {
            trees: Arc::new(trees),
            stats: vec![EntropyStats::new(); labels.len()],
            labels: Arc::new(labels),
            factors,
        })
    }

    /// Write entropy output to TSV file
    pub fn write_output(&self, output_path: &Path) -> Result<()> {
        let mut f = File::create(output_path)
            .with_context(|| format!("Failed to create output file: {:?}", output_path))?;
        writeln!(f, "label\tcount\tmean_size\tentropy")?;

        // Sort output by label name for consistency
        let mut indices: Vec<usize> = (0..self.labels.len()).collect();
        indices.sort_by_key(|&i| &self.labels[i]);

        let mut n_written = 0;
        for i in indices {
            let s = &self.stats[i];
            if s.total_count > 0.0 {
                writeln!(
                    f,
                    "{}\t{:.1}\t{:.2}\t{:.4}",
                    self.labels[i],
                    s.total_count,
                    s.mean_size(),
                    s.calculate_entropy()
                )?;
                n_written += 1;
            }
        }

        log::info!("RegionEntropy: Wrote {} labels to {:?}", n_written, output_path);
        Ok(())
    }
}

impl FragmentConsumer for RegionEntropyConsumer {
    fn name(&self) -> &str {
        "RegionEntropy"
    }

    fn consume(&mut self, frag: &Fragment) {
        if let Some(tree) = self.trees.get(&frag.chrom_id) {
            let start = frag.start as i32;
            let end = (frag.end - 1) as i32; // closed interval

            // Query for overlapping regions
            tree.query(start, end, |node| {
                let label_id = node.metadata.to_owned();

                // Calculate GC weight
                let gc_pct = (frag.gc * 100.0).round() as u8;
                let weight = self
                    .factors
                    .as_ref()
                    .map(|f| f.get_factor(frag.length, gc_pct))
                    .unwrap_or(1.0);

                // Update histogram
                let len_idx = (frag.length as usize).min(MAX_LEN);
                self.stats[label_id].len_counts[len_idx] += weight;
                self.stats[label_id].total_count += weight;
            });
        }
    }

    fn merge(&mut self, other: Self) {
        for (i, s) in other.stats.iter().enumerate() {
            self.stats[i].merge(s);
        }
    }
}

// ============================================================================
// Python Entry Points
// ============================================================================

/// Run region entropy analysis on a BED.gz fragment file
///
/// # Arguments
/// * `bed_path` - Path to input fragment BED.gz file
/// * `region_path` - Path to region BED file (chrom, start, end, label)
/// * `output_path` - Path to output TSV file (off-target/WGS-style)
/// * `gc_correction_path` - Optional path to GC correction factors TSV
/// * `target_regions_path` - Optional path to target regions BED (panel mode)
/// * `silent` - Suppress progress bar
///
/// # Panel Mode Behavior
/// When `target_regions_path` is provided:
/// - Writes off-target entropy to `output_path` (WGS-like, uses all off-target reads)
/// - Writes on-target entropy to `output_path.ontarget` suffix
/// - This matches the dual-output pattern used by FSD, FSC, OCF
#[pyfunction]
#[pyo3(signature = (bed_path, region_path, output_path, gc_correction_path=None, target_regions_path=None, silent=false))]
pub fn run_region_entropy(
    _py: Python,
    bed_path: String,
    region_path: String,
    output_path: String,
    gc_correction_path: Option<String>,
    target_regions_path: Option<String>,
    silent: bool,
) -> PyResult<(usize, usize)> {
    let mut chrom_map = ChromosomeMap::new();

    // Load GC factors if provided
    let factors = if let Some(p) = gc_correction_path {
        match CorrectionFactors::load_csv(&p) {
            Ok(f) => {
                log::info!("RegionEntropy: Loaded GC correction factors from {}", p);
                Some(Arc::new(f))
            }
            Err(e) => {
                log::warn!("RegionEntropy: Failed to load GC factors: {}", e);
                None
            }
        }
    } else {
        None
    };

    // Load target regions for panel mode (on/off-target split)
    let is_panel_mode = target_regions_path.is_some();
    let target_tree: Option<HashMap<u32, COITree<(), u32>>> = if let Some(ref p) = target_regions_path {
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<(), u32>>> = HashMap::new();
        
        if let Ok(reader) = crate::bed::get_reader(std::path::Path::new(p)) {
            for line in reader.lines().flatten() {
                if line.starts_with('#') { continue; }
                let cols: Vec<&str> = line.split('\t').collect();
                if cols.len() >= 3 {
                    let chrom = cols[0].trim_start_matches("chr");
                    let start: i32 = cols[1].parse().unwrap_or(0);
                    let end: i32 = cols[2].parse().unwrap_or(0);
                    let chrom_id = chrom_map.get_id(chrom);
                    nodes_by_chrom.entry(chrom_id).or_default()
                        .push(IntervalNode::new(start, end - 1, ()));
                }
            }
        }
        
        let mut trees = HashMap::new();
        for (k, v) in nodes_by_chrom {
            trees.insert(k, COITree::new(&v));
        }
        
        let n_regions: usize = trees.values().map(|t| t.len()).sum();
        log::info!("RegionEntropy: Panel mode enabled - {} target regions loaded", n_regions);
        Some(trees)
    } else {
        None
    };

    // Create consumers (one for off-target, one for on-target in panel mode)
    let consumer_off = RegionEntropyConsumer::new(
        std::path::Path::new(&region_path),
        &mut chrom_map,
        factors.clone(),
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    let mut consumer_on = if is_panel_mode {
        Some(RegionEntropyConsumer::new(
            std::path::Path::new(&region_path),
            &mut chrom_map,
            factors.clone(),
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?)
    } else {
        None
    };

    // Create a dual consumer wrapper for processing
    let dual_consumer = DualRegionEntropyConsumer {
        off_target: consumer_off,
        on_target: consumer_on.take(),
        target_tree,
    };

    // Process fragments
    let analyzer = FragmentAnalyzer::new(dual_consumer, 100_000);
    let final_consumer = analyzer
        .process_file(std::path::Path::new(&bed_path), &mut chrom_map, silent)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    // Write off-target output (primary/WGS-style)
    let n_off = final_consumer.off_target.stats.iter().filter(|s| s.total_count > 0.0).count();
    final_consumer.off_target
        .write_output(std::path::Path::new(&output_path))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    // Write on-target output (panel mode only)
    let n_on = if let Some(ref on_target) = final_consumer.on_target {
        // Construct proper ontarget path: foo.TFBS.raw.tsv -> foo.TFBS.ontarget.raw.tsv
        let ontarget_path = if output_path.contains(".raw.tsv") {
            output_path.replace(".raw.tsv", ".ontarget.raw.tsv")
        } else {
            output_path.replace(".tsv", ".ontarget.tsv")
        };
        
        let n = on_target.stats.iter().filter(|s| s.total_count > 0.0).count();
        on_target
            .write_output(std::path::Path::new(&ontarget_path))
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        
        log::info!("RegionEntropy: Panel mode - {} off-target, {} on-target labels", n_off, n);
        n
    } else {
        0
    };

    Ok((n_off, n_on))
}

/// Dual consumer wrapper for on/off-target split processing
#[derive(Clone)]
struct DualRegionEntropyConsumer {
    off_target: RegionEntropyConsumer,
    on_target: Option<RegionEntropyConsumer>,
    target_tree: Option<HashMap<u32, COITree<(), u32>>>,
}

impl FragmentConsumer for DualRegionEntropyConsumer {
    fn name(&self) -> &str {
        "DualRegionEntropy"
    }

    fn consume(&mut self, frag: &Fragment) {
        // Determine if fragment overlaps target regions (panel mode)
        let is_on_target = if let Some(ref trees) = self.target_tree {
            if let Some(tree) = trees.get(&frag.chrom_id) {
                let start = frag.start as i32;
                let end = (frag.end - 1) as i32;
                let mut overlaps = false;
                tree.query(start, end, |_| { overlaps = true; });
                overlaps
            } else {
                false
            }
        } else {
            false
        };

        // Route fragment to appropriate consumer
        if is_on_target {
            if let Some(ref mut on_target) = self.on_target {
                on_target.consume(frag);
            }
        } else {
            self.off_target.consume(frag);
        }
    }

    fn merge(&mut self, other: Self) {
        self.off_target.merge(other.off_target);
        if let (Some(ref mut on), Some(other_on)) = (&mut self.on_target, other.on_target) {
            on.merge(other_on);
        }
    }
}

// =============================================================================
// PON Z-SCORE NORMALIZATION
// =============================================================================

use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::RowAccessor;
use std::io::BufWriter;
use std::path::PathBuf;

/// Baseline statistics for a single label.
#[derive(Debug, Clone)]
struct EntropyBaselineStats {
    entropy_mean: f64,
    entropy_std: f64,
}

/// Load TFBS or ATAC baseline from PON Parquet file.
fn load_entropy_baseline_from_parquet(
    pon_path: &Path,
    table_name: &str,
) -> Result<HashMap<String, EntropyBaselineStats>> {
    let file = File::open(pon_path)
        .with_context(|| format!("Failed to open PON file: {:?}", pon_path))?;
    
    let reader = SerializedFileReader::new(file)
        .with_context(|| "Failed to create Parquet reader")?;
    
    let mut baseline = HashMap::new();
    
    for row_result in reader.get_row_iter(None)? {
        let row = row_result.with_context(|| "Failed to read row")?;
        
        // Check if this row belongs to the specified table
        if let Ok(table) = row.get_string(
            row.get_column_iter().position(|(name, _)| name == "table").unwrap_or(0)
        ) {
            if table != table_name {
                continue;
            }
        } else {
            continue;
        }
        
        // Extract label, entropy_mean, entropy_std
        let label = row.get_string(
            row.get_column_iter().position(|(name, _)| name == "label").unwrap_or(0)
        ).map_or("".to_string(), |v| v.to_string());
        
        let entropy_mean = row.get_double(
            row.get_column_iter().position(|(name, _)| name == "entropy_mean").unwrap_or(0)
        ).unwrap_or(0.0);
        
        let entropy_std = row.get_double(
            row.get_column_iter().position(|(name, _)| name == "entropy_std").unwrap_or(0)
        ).unwrap_or(1.0);
        
        if !label.is_empty() {
            baseline.insert(label, EntropyBaselineStats { 
                entropy_mean, 
                entropy_std 
            });
        }
    }
    
    Ok(baseline)
}

/// Apply PON z-score normalization to region entropy TSV file.
/// 
/// Reads entropy TSV with label, count, mean_size, entropy columns,
/// joins with PON baseline to compute z-scores.
/// 
/// # Arguments
/// * `entropy_path` - Path to sample entropy TSV
/// * `pon_parquet_path` - Path to PON Parquet
/// * `output_path` - Output path (None = overwrite input)
/// * `baseline_table` - Table name in PON: "tfbs_baseline" or "atac_baseline"
/// 
/// # Returns
/// Number of labels with z-scores computed
#[pyfunction]
#[pyo3(signature = (entropy_path, pon_parquet_path, output_path=None, baseline_table="tfbs_baseline"))]
pub fn apply_pon_zscore(
    entropy_path: PathBuf,
    pon_parquet_path: PathBuf,
    output_path: Option<PathBuf>,
    baseline_table: &str,
) -> PyResult<usize> {
    log::info!("RegionEntropy PON z-score: loading {} from {:?}", baseline_table, pon_parquet_path);
    
    // 1. Load baseline from PON
    let baseline = load_entropy_baseline_from_parquet(&pon_parquet_path, baseline_table)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to load {} baseline: {}", baseline_table, e)
        ))?;
    
    if baseline.is_empty() {
        log::info!("RegionEntropy PON: No {} baseline found, skipping normalization", baseline_table);
        return Ok(0);
    }
    
    log::info!("RegionEntropy PON: Loaded baseline for {} labels", baseline.len());
    
    // 2. Read entropy TSV
    let file = File::open(&entropy_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to open entropy file: {}", e)
        ))?;
    
    let reader = std::io::BufReader::new(file);
    let mut lines: Vec<String> = Vec::new();
    let mut header: Option<String> = None;
    
    for line in reader.lines() {
        let line = line.map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        if header.is_none() {
            header = Some(line);
        } else {
            lines.push(line);
        }
    }
    
    // 3. Parse header
    let header_str = header.unwrap_or_default();
    let header_cols: Vec<&str> = header_str.split('\t').collect();
    
    let label_idx = header_cols.iter().position(|&c| c == "label").unwrap_or(0);
    let entropy_idx = header_cols.iter().position(|&c| c == "entropy").unwrap_or(3);
    
    // 4. Compute z-scores
    let mut output_lines: Vec<String> = Vec::new();
    let mut n_matched = 0;
    let mut n_total = 0;
    
    for line in &lines {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() <= entropy_idx {
            output_lines.push(format!("{}\tNaN", line));
            continue;
        }
        
        n_total += 1;
        let label = cols[label_idx].trim();
        let entropy_value: f64 = cols[entropy_idx].trim().parse().unwrap_or(0.0);
        
        // Compute z-score
        let z_score = if let Some(stats) = baseline.get(label) {
            n_matched += 1;
            if stats.entropy_std > 1e-9 {
                (entropy_value - stats.entropy_mean) / stats.entropy_std
            } else {
                0.0
            }
        } else {
            f64::NAN
        };
        
        output_lines.push(format!("{}\t{:.6}", line, z_score));
    }
    
    // 5. Write output
    let out_path = output_path.unwrap_or_else(|| entropy_path.clone());
    let out_file = File::create(&out_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to create output file: {}", e)
        ))?;
    
    let mut writer = BufWriter::new(out_file);
    
    // Write header with new column
    writeln!(writer, "{}\tz_score", header_str)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    
    // Write data lines
    for line in output_lines {
        writeln!(writer, "{}", line)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    }
    
    log::info!(
        "RegionEntropy PON z-score: {}/{} labels matched ({:?})",
        n_matched, n_total, out_path.file_name().unwrap_or_default()
    );
    
    Ok(n_matched)
}
