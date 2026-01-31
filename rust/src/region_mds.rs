//! Region MDS (Motif Diversity Score) Module for krewlyzer
//!
//! Calculates per-region Motif Diversity Score from BAM files.
//! Each region (gene exon, target) gets an MDS value based on the 4-mer
//! distribution of fragment ends overlapping that region.
//!
//! # Output Format
//! - `{sample}.MDS.exon.tsv`: Per-exon/target MDS scores
//! - `{sample}.MDS.gene.tsv`: Gene-level aggregated MDS scores
//!
//! # Key Features
//! - Direct BAM processing (requires sequence data for motif extraction)
//! - E1 (first exon) detection by genomic position
//! - Flexible gene BED format handling (panel vs WGS)
//! - GC-weighted fragment counting (optional)

use pyo3::prelude::*;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, Write};
use std::collections::HashMap;
use std::sync::Arc;
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use rayon::prelude::*;
use coitrees::{COITree, IntervalNode, IntervalTree};
use indicatif::{ProgressBar, ProgressStyle};
use std::time::Duration;
use log::{info, debug, warn};

use crate::motif_utils::{reverse_complement, kmer4_to_index, calculate_mds};

// =============================================================================
// DATA STRUCTURES
// =============================================================================

/// Information about a single genomic region (exon/target)
#[derive(Debug, Clone)]
pub struct RegionInfo {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub gene: String,
    pub name: String,
    pub strand: char,
    pub chrom_id: u32,
    pub is_e1: bool,  // True if this is the first exon (E1) for its gene
}

/// Per-region MDS statistics
#[derive(Clone)]
pub struct RegionMdsStats {
    /// Histogram: 256 possible 4-mers
    motif_counts: [u64; 256],
    /// Total fragment count
    total_count: u64,
}

impl Default for RegionMdsStats {
    fn default() -> Self {
        Self::new()
    }
}

impl RegionMdsStats {
    pub fn new() -> Self {
        Self {
            motif_counts: [0u64; 256],
            total_count: 0,
        }
    }

    /// Add a 4-mer observation (optionally weighted by GC correction)
    #[inline]
    pub fn add_motif(&mut self, kmer_idx: usize, weight: f64) {
        self.motif_counts[kmer_idx] += weight as u64;
        self.total_count += 1;
    }

    /// Calculate MDS for this region
    pub fn mds(&self) -> f64 {
        calculate_mds(&self.motif_counts)
    }
}

// =============================================================================
// GENE BED PARSING
// =============================================================================

/// Detect gene BED format by column count
#[derive(Debug, Clone, Copy)]
pub enum GeneBedFormat {
    /// Panel format: chrom, start, end, gene, name (5 columns)
    Panel,
    /// WGS format: chrom, start, end, ens_id, refseq_id, gene, exon_num, strand (8 columns)
    Wgs,
    /// Unknown format
    Unknown,
}

/// Parse gene BED file and return regions with format-aware key generation
fn parse_gene_bed(
    path: &Path,
    chrom_map: &mut HashMap<String, u32>,
) -> Result<(Vec<RegionInfo>, GeneBedFormat)> {
    let reader = crate::bed::get_reader(path)?;
    let mut regions = Vec::new();
    let mut chrom_counter = 0u32;
    let mut detected_format = GeneBedFormat::Unknown;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();

        // Detect format on first data line
        if matches!(detected_format, GeneBedFormat::Unknown) {
            detected_format = match cols.len() {
                5 => GeneBedFormat::Panel,
                8 => GeneBedFormat::Wgs,
                _ => {
                    warn!("Gene BED has {} columns, treating as Panel format", cols.len());
                    GeneBedFormat::Panel
                }
            };
            info!("Gene BED format detected: {:?}", detected_format);
        }

        let chrom = cols[0].trim_start_matches("chr").to_string();
        let start: u64 = cols[1].parse().unwrap_or(0);
        let end: u64 = cols[2].parse().unwrap_or(0);

        let (gene, name, strand) = match detected_format {
            GeneBedFormat::Panel => {
                // Panel: gene=col3, name=col4
                (
                    cols.get(3).unwrap_or(&"").to_string(),
                    cols.get(4).unwrap_or(&"").to_string(),
                    '+',
                )
            }
            GeneBedFormat::Wgs => {
                // WGS: gene=col5, name=gene:exonN, strand=col7
                let gene = cols.get(5).unwrap_or(&"").to_string();
                let exon_num = cols.get(6).unwrap_or(&"0");
                let strand = cols.get(7).and_then(|s| s.chars().next()).unwrap_or('+');
                let name = format!("{}:exon{}", gene, exon_num);
                (gene, name, strand)
            }
            GeneBedFormat::Unknown => {
                continue;
            }
        };

        // Get or create chrom ID
        let chrom_id = *chrom_map.entry(chrom.clone()).or_insert_with(|| {
            let id = chrom_counter;
            chrom_counter += 1;
            id
        });

        regions.push(RegionInfo {
            chrom,
            start,
            end,
            gene,
            name,
            strand,
            chrom_id,
            is_e1: false,  // Set later by identify_e1_regions()
        });
    }

    info!("Loaded {} regions from gene BED", regions.len());
    Ok((regions, detected_format))
}

/// Identify and mark E1 (first exon) for each gene by genomic position.
/// 
/// E1 (first exon) serves as a proxy for promoter regions which are
/// Nucleosome Depleted Regions (NDRs). Per Helzer et al. (2025), E1
/// has stronger cancer signal than whole-gene averages.
/// 
/// Returns the count of E1 regions identified.
fn identify_e1_regions(regions: &mut [RegionInfo]) -> usize {
    // Group by gene, find first by position (lowest start coordinate)
    let mut gene_first: HashMap<String, (usize, u64)> = HashMap::new();

    for (idx, region) in regions.iter().enumerate() {
        gene_first
            .entry(region.gene.clone())
            .and_modify(|(best_idx, best_start)| {
                if region.start < *best_start {
                    *best_idx = idx;
                    *best_start = region.start;
                }
            })
            .or_insert((idx, region.start));
    }

    // Mark E1 regions
    for (best_idx, _) in gene_first.values() {
        regions[*best_idx].is_e1 = true;
    }

    let e1_count = gene_first.len();
    debug!("Identified {} E1 regions from {} genes", e1_count, gene_first.len());
    e1_count
}

// =============================================================================
// BAM PROCESSING
// =============================================================================

/// Build interval trees from regions for fast overlap lookup
fn build_interval_trees(regions: &[RegionInfo]) -> HashMap<u32, COITree<usize, u32>> {
    let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();

    for (idx, region) in regions.iter().enumerate() {
        nodes_by_chrom
            .entry(region.chrom_id)
            .or_default()
            .push(IntervalNode::new(region.start as i32, (region.end - 1) as i32, idx));
    }

    let mut trees = HashMap::new();
    for (chrom_id, nodes) in nodes_by_chrom {
        trees.insert(chrom_id, COITree::new(&nodes));
    }
    trees
}

/// Extract 4-mer from fragment end (strand-aware)
fn extract_end_motif(seq: &[u8], is_reverse: bool) -> Option<usize> {
    if seq.len() < 4 {
        return None;
    }

    let kmer = if is_reverse {
        // Reverse strand: take last 4 bases and reverse complement
        let end_seq = &seq[seq.len() - 4..];
        reverse_complement(end_seq)
    } else {
        // Forward strand: take first 4 bases
        seq[0..4].to_vec()
    };

    kmer4_to_index(&kmer)
}

/// Chunk definition for parallel BAM processing
#[allow(dead_code)]
struct BamChunk {
    tid: u32,
    chrom: String,
    chrom_id: u32,
    start: u64,
    end: u64,
}

/// Result from processing one chunk
struct ChunkResult {
    stats: Vec<RegionMdsStats>,
}

impl ChunkResult {
    fn new(n_regions: usize) -> Self {
        Self {
            stats: vec![RegionMdsStats::new(); n_regions],
        }
    }
}

// =============================================================================
// MAIN ENTRY POINT
// =============================================================================

/// Run region MDS analysis on a BAM file
///
/// # Arguments
/// * `bam_path` - Path to indexed BAM file
/// * `fasta_path` - Path to reference FASTA (unused, kept for API consistency)
/// * `gene_bed_path` - Path to gene BED file (panel or WGS format)
/// * `output_exon_path` - Path for per-exon output TSV
/// * `output_gene_path` - Path for gene-level output TSV
/// * `e1_only` - If true, only output E1 (first exon) results
/// * `mapq` - Minimum mapping quality
/// * `min_len` - Minimum fragment length
/// * `max_len` - Maximum fragment length
/// * `silent` - Suppress progress bar
///
/// # Returns
/// Tuple of (n_regions, n_genes) processed
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (bam_path, _fasta_path, gene_bed_path, output_exon_path, output_gene_path, e1_only=false, mapq=20, min_len=65, max_len=400, silent=false))]
pub fn run_region_mds(
    _py: Python,
    bam_path: String,
    _fasta_path: String,  // Reserved for future GC-aware processing
    gene_bed_path: String,
    output_exon_path: String,
    output_gene_path: String,
    e1_only: bool,
    mapq: u8,
    min_len: u32,
    max_len: u32,
    silent: bool,
) -> PyResult<(usize, usize)> {
    info!("Starting region-MDS analysis");
    info!("  BAM: {}", bam_path);
    info!("  Gene BED: {}", gene_bed_path);
    info!("  E1 only: {}", e1_only);

    // Parse gene BED
    let mut chrom_map: HashMap<String, u32> = HashMap::new();
    let (mut regions, _format) = parse_gene_bed(Path::new(&gene_bed_path), &mut chrom_map)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to parse gene BED: {}", e)))?;

    if regions.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("No regions found in gene BED"));
    }

    // Identify E1 regions (mark is_e1 = true)
    let n_e1 = identify_e1_regions(&mut regions);
    info!("  E1 regions: {}", n_e1);

    // Apply E1 filtering if requested
    // This reduces processing time and output size for promoter-focused analysis
    if e1_only {
        let original_count = regions.len();
        regions.retain(|r| r.is_e1);
        info!("  E1 filtering: {} -> {} regions (E1 only)", original_count, regions.len());
    }

    // Build interval trees for fast overlap lookup
    let trees = Arc::new(build_interval_trees(&regions));
    let regions_arc = Arc::new(regions);

    // Open BAM and get header info
    let bam_reader = bam::IndexedReader::from_path(&bam_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open BAM: {}", e)))?;

    let header = bam_reader.header();

    // Build BAM TID -> our chrom_id mapping
    let mut tid_to_chrom_id: HashMap<u32, u32> = HashMap::new();
    for (tid, name_bytes) in header.target_names().iter().enumerate() {
        let name = String::from_utf8_lossy(name_bytes);
        let name_stripped = name.trim_start_matches("chr");
        if let Some(&chrom_id) = chrom_map.get(name_stripped) {
            tid_to_chrom_id.insert(tid as u32, chrom_id);
        }
    }

    // Create chunks for parallel processing
    let chunk_size = 10_000_000u64; // 10Mb chunks
    let mut chunks = Vec::new();

    for (tid, name_bytes) in header.target_names().iter().enumerate() {
        let tid = tid as u32;
        let name = String::from_utf8_lossy(name_bytes).to_string();
        let name_stripped = name.trim_start_matches("chr").to_string();

        if let Some(&chrom_id) = chrom_map.get(&name_stripped) {
            let len = header.target_len(tid).unwrap_or(0);
            let mut start = 0;
            while start < len {
                let end = (start + chunk_size).min(len);
                chunks.push(BamChunk {
                    tid,
                    chrom: name.clone(),
                    chrom_id,
                    start,
                    end,
                });
                start = end;
            }
        }
    }

    info!("Split into {} chunks for parallel processing", chunks.len());

    // Progress bar
    let pb = if silent {
        ProgressBar::hidden()
    } else {
        let pb = ProgressBar::new(chunks.len() as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"));
        pb.enable_steady_tick(Duration::from_millis(100));
        pb
    };

    // Parallel processing
    let n_regions = regions_arc.len();
    let results: Vec<ChunkResult> = chunks.par_iter().map(|chunk| {
        let mut result = ChunkResult::new(n_regions);

        // Thread-local BAM reader
        let mut local_bam = match bam::IndexedReader::from_path(&bam_path) {
            Ok(b) => b,
            Err(_) => return result,
        };

        // Fetch region
        if local_bam.fetch((chunk.tid, chunk.start, chunk.end)).is_err() {
            pb.inc(1);
            return result;
        }

        // Get tree for this chromosome
        let tree = match trees.get(&chunk.chrom_id) {
            Some(t) => t,
            None => {
                pb.inc(1);
                return result;
            }
        };

        for record_res in local_bam.records() {
            let record = match record_res {
                Ok(r) => r,
                Err(_) => continue,
            };

            // Standard filters
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary()
                || record.is_quality_check_failed() || record.is_duplicate() {
                continue;
            }
            if record.mapq() < mapq {
                continue;
            }

            // Fragment length filter
            let tlen = record.insert_size().abs() as u64;
            if tlen < min_len as u64 || tlen > max_len as u64 {
                continue;
            }

            // Only process R1 to avoid double counting
            if !record.is_first_in_template() {
                continue;
            }

            // Fragment coordinates
            let frag_start = record.pos() as i32;
            let frag_end = frag_start + tlen as i32;

            // Deduplicate by chunk boundary
            if (frag_start as u64) < chunk.start || (frag_start as u64) >= chunk.end {
                continue;
            }

            // Get sequence and extract motif
            let seq = record.seq().as_bytes();
            let motif_idx = match extract_end_motif(&seq, record.is_reverse()) {
                Some(idx) => idx,
                None => continue,
            };

            // Find overlapping regions
            tree.query(frag_start, frag_end, |node| {
                let region_idx = node.metadata.to_owned();
                result.stats[region_idx].add_motif(motif_idx, 1.0);
            });
        }

        pb.inc(1);
        result
    }).collect();

    pb.finish_with_message("Done processing BAM");

    // Merge results from all chunks
    info!("Merging results from {} chunks...", results.len());
    let mut final_stats: Vec<RegionMdsStats> = vec![RegionMdsStats::new(); n_regions];

    for chunk_result in results {
        for (idx, chunk_stats) in chunk_result.stats.into_iter().enumerate() {
            for (motif_idx, count) in chunk_stats.motif_counts.iter().enumerate() {
                final_stats[idx].motif_counts[motif_idx] += count;
            }
            final_stats[idx].total_count += chunk_stats.total_count;
        }
    }

    // Write exon-level output
    let exon_path = Path::new(&output_exon_path);
    write_exon_output(exon_path, &regions_arc, &final_stats)?;

    // Write gene-level output
    let gene_path = Path::new(&output_gene_path);
    let n_genes = write_gene_output(gene_path, &regions_arc, &final_stats)?;

    info!("Region-MDS complete: {} regions, {} genes", n_regions, n_genes);
    Ok((n_regions, n_genes))
}

/// Write per-exon/target MDS output
fn write_exon_output(
    path: &Path,
    regions: &[RegionInfo],
    stats: &[RegionMdsStats],
) -> PyResult<()> {
    let mut file = File::create(path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to create exon output: {}", e)))?;

    writeln!(file, "gene\tname\tchrom\tstart\tend\tstrand\tn_fragments\tmds")?;

    for (region, stat) in regions.iter().zip(stats.iter()) {
        let mds = if stat.total_count > 0 { stat.mds() } else { 0.0 };
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
            region.gene, region.name, region.chrom, region.start, region.end,
            region.strand, stat.total_count, mds
        )?;
    }

    info!("Wrote exon output to {}", path.display());
    Ok(())
}

/// Write gene-level aggregated MDS output
fn write_gene_output(
    path: &Path,
    regions: &[RegionInfo],
    stats: &[RegionMdsStats],
) -> PyResult<usize> {
    // Aggregate by gene
    let mut gene_data: HashMap<String, Vec<(usize, u64, f64, u64)>> = HashMap::new(); // (idx, start, mds, count)

    for (idx, (region, stat)) in regions.iter().zip(stats.iter()).enumerate() {
        let mds = if stat.total_count > 0 { stat.mds() } else { 0.0 };
        gene_data
            .entry(region.gene.clone())
            .or_default()
            .push((idx, region.start, mds, stat.total_count));
    }

    let mut file = File::create(path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to create gene output: {}", e)))?;

    writeln!(file, "gene\tn_exons\tn_fragments\tmds_mean\tmds_e1\tmds_std")?;

    for (gene, regions_data) in &gene_data {
        // Sort by position to find E1
        let mut sorted_data = regions_data.clone();
        sorted_data.sort_by_key(|(_, start, _, _)| *start);

        let mds_values: Vec<f64> = sorted_data.iter()
            .filter(|(_, _, _, count)| *count > 0)
            .map(|(_, _, mds, _)| *mds)
            .collect();

        if mds_values.is_empty() {
            continue;
        }

        let n_exons = sorted_data.len();
        let n_fragments: u64 = sorted_data.iter().map(|(_, _, _, c)| c).sum();
        let mds_mean = mds_values.iter().sum::<f64>() / mds_values.len() as f64;
        let mds_e1 = sorted_data.iter()
            .find(|(_, _, _, count)| *count > 0)
            .map(|(_, _, mds, _)| *mds)
            .unwrap_or(0.0);

        let mds_std = if mds_values.len() > 1 {
            let variance = mds_values.iter()
                .map(|x| (x - mds_mean).powi(2))
                .sum::<f64>() / mds_values.len() as f64;
            variance.sqrt()
        } else {
            0.0
        };

        writeln!(
            file,
            "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}",
            gene, n_exons, n_fragments, mds_mean, mds_e1, mds_std
        )?;
    }

    let n_genes = gene_data.len();
    info!("Wrote gene output to {} ({} genes)", path.display(), n_genes);
    Ok(n_genes)
}
