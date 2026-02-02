//! Orientation-aware cfDNA Fragmentation (OCF) calculation
//!
//! Measures fragment end orientation patterns around open chromatin regions (OCRs).
//! OCF captures tissue-of-origin signals based on strand asymmetry at regulatory elements.
//!
//! ## Key Features
//! - **Strand asymmetry**: OCF = (Upstream - Downstream) / (Upstream + Downstream)
//! - **Region types**: Calculates per-OCR and synchronized (pooled) scores
//! - **GC correction**: Per-fragment weight adjustment via correction factors
//! - **Triple-output (panel mode)**: Inclusive routing produces:
//!   - ALL fragments (primary, WGS-comparable)
//!   - On-target fragments only
//!   - Off-target fragments only
//!
//! ## Output Files (6 in panel mode, 2 in WGS mode)
//! - all.ocf.tsv: Summary scores (all fragments)
//! - all.sync.tsv: Detailed sync data (all fragments)
//! - all.ocf.ontarget.tsv: Summary scores (on-target only)
//! - all.sync.ontarget.tsv: Detailed sync data (on-target only)
//! - all.ocf.offtarget.tsv: Summary scores (off-target only)
//! - all.sync.offtarget.tsv: Detailed sync data (off-target only)

use pyo3::prelude::*;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, Write};
use std::collections::HashMap;


/// Calculate Orientation-aware cfDNA Fragmentation (OCF)
/// 
/// # Arguments
/// * `bed_path` - Path to the input .bed.gz file
/// * `ocr_path` - Path to the Open Chromatin Regions (OCR) BED file
/// * `output_dir` - Directory to save output files
#[pyfunction]
pub fn calculate_ocf(
    _py: Python,
    bed_path: PathBuf,
    ocr_path: PathBuf,
    output_dir: PathBuf,
) -> PyResult<()> {
    let mut chrom_map = ChromosomeMap::new();
    let consumer = OcfConsumer::new(&ocr_path, &mut chrom_map, None, None)  // No target_regions for legacy API
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        
    let analyzer = FragmentAnalyzer::new(consumer, 100_000);
    let final_consumer = analyzer.process_file(&bed_path, &mut chrom_map, false)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        
    final_consumer.write_output(&output_dir)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    Ok(())
}

use crate::bed::{Fragment, ChromosomeMap};
use crate::engine::{FragmentConsumer, FragmentAnalyzer};
use std::sync::Arc;
use std::path::Path;
use anyhow::{Result, Context};
use coitrees::{COITree, IntervalNode, IntervalTree};
use crate::gc_correction::CorrectionFactors;

#[derive(Clone, Default)]
struct LabelStats {
    left_pos: Vec<f64>,  // Size 2000
    right_pos: Vec<f64>, // Size 2000
    total_starts: f64,
    total_ends: f64,
}

impl LabelStats {
    fn new() -> Self {
        Self {
            left_pos: vec![0.0; 2000],
            right_pos: vec![0.0; 2000],
            total_starts: 0.0,
            total_ends: 0.0,
        }
    }
    
    fn merge(&mut self, other: &Self) {
        for (a, b) in self.left_pos.iter_mut().zip(other.left_pos.iter()) { *a += *b; }
        for (a, b) in self.right_pos.iter_mut().zip(other.right_pos.iter()) { *a += *b; }
        self.total_starts += other.total_starts;
        self.total_ends += other.total_ends;
    }
}

// Region with pre-mapped label ID
struct OcrRegionInfo {
    // chrom stored in tree map key
    start: u64,
    end: u64,
    label_id: usize,
}

/// Parse a single OCR line
fn parse_ocr_line(
    line: &str,
    chrom_map: &mut ChromosomeMap,
    label_to_id: &mut HashMap<String, usize>,
    labels: &mut Vec<String>,
    nodes_by_chrom: &mut HashMap<u32, Vec<IntervalNode<usize, u32>>>,
    region_infos: &mut Vec<OcrRegionInfo>,
) {
    let line = line.trim();
    if line.is_empty() || line.starts_with('#') {
        return;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 4 {
        return;
    }
    
    let chrom = fields[0];
    let start: u64 = fields[1].parse().unwrap_or(0);
    let end: u64 = fields[2].parse().unwrap_or(0);
    let label = fields[3].to_string();
    
    // Map Label
    let label_id = if let Some(&id) = label_to_id.get(&label) {
        id
    } else {
        let id = labels.len();
        label_to_id.insert(label.clone(), id);
        labels.push(label);
        id
    };
    
    // Map Chrom
    let chrom_norm = chrom.trim_start_matches("chr");
    let chrom_id = chrom_map.get_id(chrom_norm);
    
    // Store Region Info
    let info_idx = region_infos.len();
    region_infos.push(OcrRegionInfo {
        start,
        end,
        label_id,
    });
    
    // Add to Tree nodes
    // Standard overlap query.
    // COITree (closed): [start, end-1]
    let s = start as u32;
    let e = end as u32;
    let e_closed = if e > s { e - 1 } else { s };
    
    nodes_by_chrom.entry(chrom_id).or_default().push(
        IntervalNode::new(s as i32, e_closed as i32, info_idx)
    );
}

#[derive(Clone)]
pub struct OcfConsumer {
    // Shared state
    trees: Arc<HashMap<u32, COITree<usize, u32>>>,
    region_infos: Arc<Vec<OcrRegionInfo>>, // Index -> Info
    labels: Arc<Vec<String>>, // LabelID -> Name
    factors: Option<Arc<CorrectionFactors>>,
    
    // Target regions for on/off-target split (panel data)
    target_tree: Option<Arc<HashMap<u32, COITree<(), u32>>>>,
    
    // Thread-local state - Indexed by LabelID
    // ALL fragments (primary output, WGS-comparable)
    stats_all: Vec<LabelStats>,
    // Off-target fragments only (panel mode)
    stats_off: Vec<LabelStats>,
    // On-target fragments only (panel mode)
    stats_on: Vec<LabelStats>,
}

impl OcfConsumer {
    pub fn new(
        ocr_path: &Path, 
        chrom_map: &mut ChromosomeMap, 
        factors: Option<Arc<CorrectionFactors>>,
        target_regions: Option<Arc<HashMap<u32, COITree<(), u32>>>>,
    ) -> Result<Self> {
        // 1. Parse OCR File
        // Format: chrom start end label
        let reader = crate::bed::get_reader(ocr_path)?;
        
        let mut label_to_id: HashMap<String, usize> = HashMap::new();
        let mut labels: Vec<String> = Vec::new();
        
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        let mut region_infos = Vec::new(); // Implicitly indexed by order of insertion
        
        for line in reader.lines() {
            let line = line?;
            parse_ocr_line(&line, chrom_map, &mut label_to_id, &mut labels, &mut nodes_by_chrom, &mut region_infos);
        }
        
        // Build Trees
        let mut trees = HashMap::new();
        for (cid, nodes) in nodes_by_chrom {
            trees.insert(cid, COITree::new(&nodes));
        }
        
        // Init stats (all, off-target, on-target)
        let stats_all = vec![LabelStats::new(); labels.len()];
        let stats_off = vec![LabelStats::new(); labels.len()];
        let stats_on = vec![LabelStats::new(); labels.len()];
        
        if target_regions.is_some() {
            log::info!("OCF: Panel mode enabled, triple output active (all/on/off)");
        }
        
        Ok(Self {
            trees: Arc::new(trees),
            region_infos: Arc::new(region_infos),
            labels: Arc::new(labels),
            factors,
            target_tree: target_regions,
            stats_all,
            stats_off,
            stats_on,
        })
    }
    
    /// Check if a fragment overlaps any target region
    fn is_on_target(&self, fragment: &Fragment) -> bool {
        if let Some(ref target_tree) = self.target_tree {
            if let Some(tree) = target_tree.get(&fragment.chrom_id) {
                let start = fragment.start as i32;
                let end = if fragment.end > fragment.start { (fragment.end - 1) as i32 } else { start };
                let mut found = false;
                tree.query(start, end, |_| { found = true; });
                return found;
            }
        }
        false
    }
    
    pub fn write_output(&self, output_dir: &Path) -> Result<()> {
        // Write primary stats (ALL fragments - WGS-comparable)
        self.write_stats(&self.stats_all, output_dir, "all.ocf.tsv", "all.sync.tsv")?;
        log::info!("OCF: Wrote primary output files (all fragments)");
        
        // Write on-target stats if we have any data (panel mode)
        let has_on_target_data = self.stats_on.iter().any(|s| s.total_starts > 0.0 || s.total_ends > 0.0);
        if has_on_target_data {
            self.write_stats(&self.stats_on, output_dir, "all.ocf.ontarget.tsv", "all.sync.ontarget.tsv")?;
            log::info!("OCF: Wrote on-target output files");
        }
        
        // Write off-target stats if we have any data (panel mode)
        let has_off_target_data = self.stats_off.iter().any(|s| s.total_starts > 0.0 || s.total_ends > 0.0);
        if has_off_target_data {
            self.write_stats(&self.stats_off, output_dir, "all.ocf.offtarget.tsv", "all.sync.offtarget.tsv")?;
            log::info!("OCF: Wrote off-target output files");
        }
        
        Ok(())
    }
    
    fn write_stats(&self, stats: &[LabelStats], output_dir: &Path, summary_name: &str, sync_name: &str) -> Result<()> {
        // Output 1: OCF summary
        let summary_path = output_dir.join(summary_name);
        let mut summary_file = File::create(&summary_path)
            .with_context(|| format!("Failed to create summary file: {:?}", summary_path))?;
        writeln!(summary_file, "tissue\tOCF")?;

        // Output 2: sync positions
        let sync_path = output_dir.join(sync_name);
        let mut sync_file = File::create(&sync_path)
             .with_context(|| format!("Failed to create sync file: {:?}", sync_path))?;
        writeln!(sync_file, "tissue\tposition\tleft_count\tleft_norm\tright_count\tright_norm")?;

        // Sort labels
        let mut label_indices: Vec<usize> = (0..self.labels.len()).collect();
        label_indices.sort_by_key(|&i| &self.labels[i]);
        
        for &id in &label_indices {
            let label = &self.labels[id];
            let s = &stats[id];
            
            let ts = if s.total_starts > 0.0 { s.total_starts / 10000.0 } else { 1.0 };
            let te = if s.total_ends > 0.0 { s.total_ends / 10000.0 } else { 1.0 };
            
            let peak = 60;
            let bin_width = 10;
            let mut trueends = 0.0;
            let mut background = 0.0;
            
            for k in 0..2000 {
                let l_count = s.left_pos[k];
                let r_count = s.right_pos[k];
                let l_norm = l_count / ts;
                let r_norm = r_count / te;
                let loc = k as i64 - 1000;
                
                writeln!(sync_file, "{}\t{}\t{:.2}\t{:.6}\t{:.2}\t{:.6}", 
                    label, loc, l_count, l_norm, r_count, r_norm)?;
                
                if loc >= -peak - bin_width && loc <= -peak + bin_width {
                    trueends += r_norm;
                    background += l_norm;
                } else if loc >= peak - bin_width && loc <= peak + bin_width {
                    trueends += l_norm;
                    background += r_norm;
                }
            }
            
            let ocf_score = trueends - background;
            writeln!(summary_file, "{}\t{:.6}", label, ocf_score)?;
        }
        
        Ok(())
    }
}

impl FragmentConsumer for OcfConsumer {
    fn name(&self) -> &str {
        "OCF"
    }

    fn consume(&mut self, fragment: &Fragment) {
        // Check if fragment overlaps target regions (panel mode routing)
        let is_panel_mode = self.target_tree.is_some();
        let on_target = if is_panel_mode { self.is_on_target(fragment) } else { false };
        
        if let Some(tree) = self.trees.get(&fragment.chrom_id) {
            let start = fragment.start as i32;
            let end_closed = if fragment.end > fragment.start { (fragment.end - 1) as i32 } else { start };
            
            tree.query(start, end_closed, |node| {
                let region_idx = node.metadata.to_owned();
                let region = &self.region_infos[region_idx];
                
                // Compute weight once
                let gc_pct = (fragment.gc * 100.0).round() as u8;
                let weight = if let Some(ref factors) = self.factors {
                    factors.get_factor(fragment.length, gc_pct)
                } else { 1.0 };
                
                // Helper closure to update stats
                let update_stats = |label_stats: &mut LabelStats| {
                    let r_start = fragment.start;
                    let r_end = fragment.end;
                    let reg_start = region.start;
                    let reg_end = region.end;
                    
                    // Starts
                    if r_start >= reg_start {
                        let s = (r_start - reg_start) as usize;
                        if s < 2000 {
                            label_stats.left_pos[s] += weight;
                        }
                        label_stats.total_starts += weight;
                    }
                    
                    // Ends
                    if r_end <= reg_end {
                        let diff = r_end.wrapping_sub(reg_start).wrapping_add(1);
                        let e = diff as usize;
                        if e < 2000 {
                            label_stats.right_pos[e] += weight;
                        }
                        label_stats.total_ends += weight;
                    }
                };
                
                // INCLUSIVE routing: ALL fragments go to stats_all
                update_stats(&mut self.stats_all[region.label_id]);
                
                // Panel mode: additionally route to on-target or off-target
                if is_panel_mode {
                    if on_target {
                        update_stats(&mut self.stats_on[region.label_id]);
                    } else {
                        update_stats(&mut self.stats_off[region.label_id]);
                    }
                }
            });
        }
    }

    fn merge(&mut self, other: Self) {
        // Merge all stats (primary)
        for (i, other_s) in other.stats_all.iter().enumerate() {
            self.stats_all[i].merge(other_s);
        }
        // Merge off-target stats
        for (i, other_s) in other.stats_off.iter().enumerate() {
            self.stats_off[i].merge(other_s);
        }
        // Merge on-target stats
        for (i, other_s) in other.stats_on.iter().enumerate() {
            self.stats_on[i].merge(other_s);
        }
    }
}

// =============================================================================
// PON Z-SCORE NORMALIZATION
// =============================================================================

use parquet::file::reader::FileReader;
use parquet::file::reader::SerializedFileReader;
use parquet::record::RowAccessor;

/// OCF baseline statistics for a single tissue/region.
#[derive(Debug, Clone)]
struct OcfBaselineStats {
    ocf_mean: f64,
    ocf_std: f64,
}

/// Load OCF baseline from PON Parquet file.
/// 
/// Reads the `ocf_baseline` table from the PON file and returns a HashMap
/// mapping tissue/region names to their baseline statistics.
fn load_ocf_baseline_from_parquet(pon_path: &Path) -> Result<HashMap<String, OcfBaselineStats>> {
    let file = File::open(pon_path)
        .with_context(|| format!("Failed to open PON file: {:?}", pon_path))?;
    
    let reader = SerializedFileReader::new(file)
        .with_context(|| "Failed to create Parquet reader")?;
    
    let mut baseline = HashMap::new();
    
    for row_result in reader.get_row_iter(None)? {
        let row = row_result.with_context(|| "Failed to read row")?;
        
        // Check if this row belongs to ocf_baseline table
        if let Ok(table) = row.get_string(
            row.get_column_iter().position(|(name, _)| name == "table").unwrap_or(0)
        ) {
            if table != "ocf_baseline" {
                continue;
            }
        } else {
            continue;
        }
        
        // Extract region_id, ocf_mean, ocf_std
        let region_id = row.get_string(
            row.get_column_iter().position(|(name, _)| name == "region_id").unwrap_or(0)
        ).map_or("".to_string(), |v| v.to_string());
        
        let ocf_mean = row.get_double(
            row.get_column_iter().position(|(name, _)| name == "ocf_mean").unwrap_or(0)
        ).unwrap_or(0.0);
        
        let ocf_std = row.get_double(
            row.get_column_iter().position(|(name, _)| name == "ocf_std").unwrap_or(0)
        ).unwrap_or(1.0);
        
        if !region_id.is_empty() {
            baseline.insert(region_id, OcfBaselineStats { 
                ocf_mean, 
                ocf_std 
            });
        }
    }
    
    Ok(baseline)
}

/// Apply PON z-score normalization to OCF TSV file.
/// 
/// Reads OCF TSV with tissue/region and OCF columns, joins with PON baseline
/// to compute z-scores: z = (observed - mean) / std
/// 
/// # Arguments
/// * `ocf_path` - Path to sample OCF TSV (tissue\tOCF format)
/// * `pon_parquet_path` - Path to PON Parquet with ocf_baseline table
/// * `output_path` - Output path (None = overwrite input)
/// 
/// # Returns
/// Number of regions with z-scores computed
/// 
/// # Performance
/// 10-50x faster than Python pandas iterrows implementation
#[pyfunction]
#[pyo3(signature = (ocf_path, pon_parquet_path, output_path=None))]
pub fn apply_pon_zscore(
    ocf_path: PathBuf,
    pon_parquet_path: PathBuf,
    output_path: Option<PathBuf>,
) -> PyResult<usize> {
    use std::io::BufWriter;
    
    log::info!("OCF PON z-score: loading PON from {:?}", pon_parquet_path);
    
    // 1. Load OCF baseline from PON
    let baseline = load_ocf_baseline_from_parquet(&pon_parquet_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to load OCF baseline: {}", e)
        ))?;
    
    if baseline.is_empty() {
        log::info!("OCF PON: No OCF baseline found, skipping normalization");
        return Ok(0);
    }
    
    log::info!("OCF PON: Loaded baseline for {} tissues", baseline.len());
    
    // 2. Read OCF TSV
    let file = File::open(&ocf_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to open OCF file: {}", e)
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
    
    // 3. Parse header to find column indices
    let header_str = header.unwrap_or_default();
    let header_cols: Vec<&str> = header_str.split('\t').collect();
    
    // Find tissue/region column (first column) and OCF column
    let tissue_idx = 0;  // First column is typically tissue/region
    let ocf_idx = header_cols.iter().position(|&c| 
        c.eq_ignore_ascii_case("ocf") || c.eq_ignore_ascii_case("ocf_score") || c.eq_ignore_ascii_case("score")
    ).unwrap_or(1);  // Default to second column
    
    // 4. Compute z-scores
    let mut output_lines: Vec<String> = Vec::new();
    let mut n_matched = 0;
    let mut n_total = 0;
    
    for line in &lines {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() <= ocf_idx {
            output_lines.push(format!("{}\tNaN", line));
            continue;
        }
        
        n_total += 1;
        let tissue = cols[tissue_idx].trim();
        let ocf_value: f64 = cols[ocf_idx].trim().parse().unwrap_or(0.0);
        
        // Compute z-score
        let z_score = if let Some(stats) = baseline.get(tissue) {
            n_matched += 1;
            if stats.ocf_std > 1e-9 {
                (ocf_value - stats.ocf_mean) / stats.ocf_std
            } else {
                0.0
            }
        } else {
            f64::NAN
        };
        
        output_lines.push(format!("{}\t{:.6}", line, z_score));
    }
    
    // 5. Write output
    let out_path = output_path.unwrap_or_else(|| ocf_path.clone());
    let out_file = File::create(&out_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to create output file: {}", e)
        ))?;
    
    let mut writer = BufWriter::new(out_file);
    
    // Write header with new column
    writeln!(writer, "{}\tocf_z", header_str)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    
    // Write data lines
    for line in output_lines {
        writeln!(writer, "{}", line)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    }
    
    log::info!(
        "OCF PON z-score: {}/{} regions matched ({:?})",
        n_matched, n_total, out_path.file_name().unwrap_or_default()
    );
    
    Ok(n_matched)
}
