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
/// * `output_path` - Path to output TSV file
/// * `gc_correction_path` - Optional path to GC correction factors TSV
#[pyfunction]
#[pyo3(signature = (bed_path, region_path, output_path, gc_correction_path=None, silent=false))]
pub fn run_region_entropy(
    _py: Python,
    bed_path: String,
    region_path: String,
    output_path: String,
    gc_correction_path: Option<String>,
    silent: bool,
) -> PyResult<()> {
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

    // Create consumer
    let consumer = RegionEntropyConsumer::new(
        std::path::Path::new(&region_path),
        &mut chrom_map,
        factors,
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    // Process fragments
    let analyzer = FragmentAnalyzer::new(consumer, 100_000);
    let final_consumer = analyzer
        .process_file(std::path::Path::new(&bed_path), &mut chrom_map, silent)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    // Write output
    final_consumer
        .write_output(std::path::Path::new(&output_path))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    Ok(())
}
