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
    /// Legacy panel format: chrom, start, end, gene, name (5 columns).
    ///
    /// Carries no strand, so E1 cannot be resolved correctly from it -- see
    /// the note in `parse_gene_bed`.
    Panel,
    /// Legacy WGS format: chrom, start, end, ens_id, refseq_id, gene,
    /// exon_num, strand (8 columns). `exon_num` here is *coordinate*-ordered
    /// and must not be used to identify the first exon.
    Wgs,
    /// Current format from `scripts/build_gene_bed.py` (11 columns):
    /// chrom, start, end, gene, name, transcript_id, exon_number, strand,
    /// is_e1, is_alt_e1, is_first_captured.
    ///
    /// `exon_number` is transcription-ordered and `is_e1` is precomputed, so
    /// nothing here has to be re-derived at runtime.
    Annotated,
    /// Unknown format
    Unknown,
}

/// Column count of [`GeneBedFormat::Annotated`], i.e. of `HEADER` in
/// `scripts/build_gene_bed.py`. Detection is by width, so the two must agree.
const ANNOTATED_COLUMNS: usize = 11;

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
                ANNOTATED_COLUMNS => GeneBedFormat::Annotated,
                5 => {
                    // No strand column exists in this format, so E1 falls back
                    // to lowest-coordinate and is wrong for every minus-strand
                    // gene. Say so rather than producing a plausible number.
                    warn!(
                        "Gene BED has 5 columns (legacy panel format) and carries \
                         no strand: E1 will be the lowest-coordinate region for \
                         every gene, which is the last exon on the minus strand. \
                         Regenerate with scripts/build_gene_bed.py."
                    );
                    GeneBedFormat::Panel
                }
                8 => GeneBedFormat::Wgs,
                n => {
                    warn!(
                        "Gene BED has {} columns, which matches no known format; \
                         treating as legacy panel (gene=col4, name=col5, no \
                         strand). E1 will not be strand-aware.",
                        n
                    );
                    GeneBedFormat::Panel
                }
            };
            info!("Gene BED format detected: {:?}", detected_format);
        }

        let chrom = cols[0].trim_start_matches("chr").to_string();
        let start: u64 = cols[1].parse().unwrap_or(0);
        let end: u64 = cols[2].parse().unwrap_or(0);

        // `prebuilt_e1` is Some only for the annotated format, where the flag
        // was resolved at build time against a GENCODE transcript. For the two
        // legacy formats it stays None and E1 is derived from coordinates.
        let (gene, name, strand, prebuilt_e1) = match detected_format {
            GeneBedFormat::Annotated => {
                // gene=col3, name=col4, strand=col7, is_e1=col8
                let gene = cols.get(3).unwrap_or(&"").to_string();
                let name = cols.get(4).unwrap_or(&"").to_string();
                let strand = cols.get(7).and_then(|s| s.chars().next()).unwrap_or('+');
                (gene, name, strand, Some(cols.get(8) == Some(&"1")))
            }
            GeneBedFormat::Panel => {
                // Legacy panel: gene=col3, name=col4. No strand column exists,
                // so '+' is a placeholder, not a claim -- warned about above.
                (
                    cols.get(3).unwrap_or(&"").to_string(),
                    cols.get(4).unwrap_or(&"").to_string(),
                    '+',
                    None,
                )
            }
            GeneBedFormat::Wgs => {
                // Legacy WGS: gene=col5, name=gene:exonN, strand=col7.
                // col6 is coordinate-ordered and deliberately not used for E1.
                let gene = cols.get(5).unwrap_or(&"").to_string();
                let exon_num = cols.get(6).unwrap_or(&"0");
                let strand = cols.get(7).and_then(|s| s.chars().next()).unwrap_or('+');
                let name = format!("{}:exon{}", gene, exon_num);
                (gene, name, strand, None)
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
            // For the annotated format this is final; identify_e1_regions()
            // leaves it alone. For the legacy formats it is a placeholder that
            // identify_e1_regions() fills in from coordinates and strand.
            is_e1: prebuilt_e1.unwrap_or(false),
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
    // If the asset already carries is_e1, trust it and derive nothing.
    //
    // The build-time flag comes from a GENCODE transcript's exon_number, which
    // is the actual answer. The coordinate heuristic below is only a proxy: it
    // returns the most 5' *captured* region, which on a targeted panel is
    // usually an internal exon rather than the first one.
    let prebuilt = regions.iter().filter(|r| r.is_e1).count();
    if prebuilt > 0 {
        debug!(
            "Using {} precomputed E1 region(s) from the gene BED; not deriving \
             from coordinates",
            prebuilt
        );
        return prebuilt;
    }

    // Group by gene and find the TRANSCRIPTIONALLY first exon.
    //
    // E1 is defined by transcription order, not genomic coordinate order:
    //   '+' strand -> lowest  start coordinate
    //   '-' strand -> highest start coordinate
    //
    // Ignoring strand (the historical behaviour) silently reported the LAST
    // exon as E1 for every minus-strand gene, i.e. roughly half of any panel.
    // Note this can only be as good as the strand it is given: the legacy
    // 5-column panel format has none, so every gene there looks plus-stranded.
    let mut gene_first: HashMap<String, (usize, u64)> = HashMap::new();

    for (idx, region) in regions.iter().enumerate() {
        let is_minus = region.strand == '-';
        gene_first
            .entry(region.gene.clone())
            .and_modify(|(best_idx, best_start)| {
                let better = if is_minus {
                    region.start > *best_start
                } else {
                    region.start < *best_start
                };
                if better {
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
    py: Python,
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

    // Release GIL before par_iter to prevent pyo3-log deadlock
    let n_regions = regions_arc.len();
    let results: Vec<ChunkResult> = py.allow_threads(|| {
        chunks.par_iter().map(|chunk| {
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
            let tlen_signed = record.insert_size();
            let tlen = tlen_signed.unsigned_abs();
            if tlen < min_len as u64 || tlen > max_len as u64 {
                continue;
            }

            // Only process R1 to avoid double counting
            if !record.is_first_in_template() {
                continue;
            }

            // Deduplicate by chunk boundary BEFORE deriving fragment
            // coordinates. Chunks overlap, so each read must be claimed by
            // exactly one of them; the read's own alignment start is the only
            // key that guarantees that. Keying on the fragment start instead
            // would move with orientation and could claim a read twice or not
            // at all.
            if (record.pos() as u64) < chunk.start || (record.pos() as u64) >= chunk.end {
                continue;
            }

            // Fragment coordinates span the outermost bases of the pair, so
            // pos() + |tlen| is wrong when R1 is the rightmost mate -- see the
            // same correction in extract_motif.rs. These feed the interval
            // query below, so a shifted fragment is attributed to the wrong
            // exon.
            let (frag_start, frag_end) = if tlen_signed >= 0 {
                let s = record.pos() as i32;
                (s, s + tlen as i32)
            } else {
                let e = record.cigar().end_pos() as i32;
                (e - tlen as i32, e)
            };

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
    }).collect()
    }); // end py.allow_threads

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
        // Sort by genomic position for stable, reproducible output ordering.
        // NOTE: E1 is NOT re-derived here -- it is read from the `is_e1` flag
        // set by identify_e1_regions(), which is strand-aware. Re-deriving it
        // by coordinate (the historical behaviour) both ignored strand and
        // fell through to the next exon whenever E1 had zero fragments.
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
        // Strict E1: the region flagged by identify_e1_regions() for this gene.
        // 0.0 means "E1 had no fragments", NOT "the next exon's value".
        let mds_e1 = sorted_data.iter()
            .find(|(idx, _, _, _)| regions[*idx].is_e1)
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

// =============================================================================
// TESTS
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Write a gene BED to a temp file and parse it, so the format-detection
    /// path is actually exercised.
    ///
    /// The existing E1 tests build `RegionInfo` by hand and never reach
    /// `parse_gene_bed`. That gap let a real regression through: adding
    /// columns to the bundled assets took them from 8 to 11 fields, which
    /// matched no arm of the detector and silently fell back to the legacy
    /// panel branch -- where strand is a hardcoded '+', so strand-aware E1
    /// stopped working for WGS without a single test failing.
    fn parse_str(body: &str) -> (Vec<RegionInfo>, GeneBedFormat) {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("krewlyzer_gene_bed_{}.bed", body.len()));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        drop(f);
        let mut map = HashMap::new();
        let out = parse_gene_bed(&path, &mut map).unwrap();
        let _ = std::fs::remove_file(&path);
        out
    }

    #[test]
    fn the_annotated_format_is_detected_by_its_width() {
        let body = "#chrom\tstart\tend\tgene\tname\ttranscript_id\texon_number\t\
strand\tis_e1\tis_alt_e1\tis_first_captured\n\
chr1\t100\t200\tTP53\ttile_a\tENST1\t1\t-\t1\t0\t1\n";
        let (regions, fmt) = parse_str(body);
        assert!(
            matches!(fmt, GeneBedFormat::Annotated),
            "11-column asset detected as {:?}, not Annotated -- it would fall \
             back to the legacy panel branch and lose strand",
            fmt
        );
        assert_eq!(regions[0].gene, "TP53");
        assert_eq!(regions[0].strand, '-', "strand must come from column 8");
        assert!(regions[0].is_e1, "is_e1 must be read from column 9");
    }

    #[test]
    fn the_annotated_format_carries_strand_for_both_assay_families() {
        // The regression in full: a minus-strand gene whose E1 is at its
        // highest coordinate. Under the legacy panel branch every gene reads
        // as '+' and the lowest coordinate wins.
        let body = "chr1\t100\t200\tBRCA1\ta\tENST1\t3\t-\t0\t0\t0\n\
chr1\t900\t1000\tBRCA1\tb\tENST1\t1\t-\t1\t0\t1\n";
        let (mut regions, fmt) = parse_str(body);
        assert!(matches!(fmt, GeneBedFormat::Annotated));
        assert!(regions.iter().all(|r| r.strand == '-'));

        let n = identify_e1_regions(&mut regions);
        assert_eq!(n, 1);
        let e1: Vec<_> = regions.iter().filter(|r| r.is_e1).collect();
        assert_eq!(e1.len(), 1);
        assert_eq!(
            e1[0].start, 900,
            "E1 should be the highest-coordinate region on a minus-strand gene"
        );
    }

    #[test]
    fn a_precomputed_e1_flag_is_not_overwritten() {
        // The build-time flag comes from a GENCODE exon_number; the coordinate
        // heuristic is only a proxy. Re-deriving would discard the better
        // answer -- and on a panel would move E1 onto an internal exon.
        let body = "chr1\t100\t200\tMYC\ta\tENST1\t5\t+\t0\t0\t1\n\
chr1\t900\t1000\tMYC\tb\tENST1\t1\t+\t1\t0\t0\n";
        let (mut regions, _) = parse_str(body);
        identify_e1_regions(&mut regions);
        let e1: Vec<_> = regions.iter().filter(|r| r.is_e1).collect();
        assert_eq!(e1.len(), 1);
        assert_eq!(
            e1[0].start, 900,
            "the precomputed flag was discarded and E1 re-derived from \
             coordinates, which would pick the lowest start on a + strand"
        );
    }

    #[test]
    fn legacy_formats_still_parse() {
        let (panel, fmt) = parse_str("chr1\t100\t200\tTP53\texon_01\n");
        assert!(matches!(fmt, GeneBedFormat::Panel));
        assert_eq!(panel[0].gene, "TP53");

        let (wgs, fmt) =
            parse_str("chr1\t100\t200\tENSG1\tNM_1\tTP53\t2\t-\n");
        assert!(matches!(fmt, GeneBedFormat::Wgs));
        assert_eq!(wgs[0].gene, "TP53");
        assert_eq!(wgs[0].strand, '-');
    }

    #[test]
    fn annotated_column_count_matches_the_builder() {
        // Detection is by width, so this constant and the HEADER in
        // scripts/build_gene_bed.py must stay in step.
        let header = "chrom start end gene name transcript_id exon_number \
strand is_e1 is_alt_e1 is_first_captured";
        assert_eq!(header.split_whitespace().count(), ANNOTATED_COLUMNS);
    }

    fn region(gene: &str, name: &str, start: u64, strand: char) -> RegionInfo {
        RegionInfo {
            chrom: "1".to_string(),
            start,
            end: start + 100,
            gene: gene.to_string(),
            name: name.to_string(),
            strand,
            chrom_id: 0,
            is_e1: false,
        }
    }

    #[test]
    fn e1_is_lowest_coordinate_on_plus_strand() {
        let mut regions = vec![
            region("TP53", "exon_01", 1_000, '+'),
            region("TP53", "exon_02", 2_000, '+'),
            region("TP53", "exon_03", 3_000, '+'),
        ];
        let n = identify_e1_regions(&mut regions);
        assert_eq!(n, 1);
        assert!(regions[0].is_e1, "plus-strand E1 must be the lowest start");
        assert!(!regions[1].is_e1);
        assert!(!regions[2].is_e1);
    }

    #[test]
    fn e1_is_highest_coordinate_on_minus_strand() {
        // Regression: strand was previously ignored, so the LAST exon of every
        // minus-strand gene was reported as E1.
        let mut regions = vec![
            region("BRCA1", "exon_01", 1_000, '-'),
            region("BRCA1", "exon_02", 2_000, '-'),
            region("BRCA1", "exon_03", 3_000, '-'),
        ];
        let n = identify_e1_regions(&mut regions);
        assert_eq!(n, 1);
        assert!(
            regions[2].is_e1,
            "minus-strand E1 must be the HIGHEST start (transcription order)"
        );
        assert!(!regions[0].is_e1);
        assert!(!regions[1].is_e1);
    }

    #[test]
    fn e1_is_per_gene_and_strand_aware_when_mixed() {
        let mut regions = vec![
            region("PLUS", "e1", 100, '+'),
            region("PLUS", "e2", 900, '+'),
            region("MINUS", "e1", 200, '-'),
            region("MINUS", "e2", 800, '-'),
        ];
        let n = identify_e1_regions(&mut regions);
        assert_eq!(n, 2, "one E1 per gene");
        assert!(regions[0].is_e1, "PLUS E1 = lowest start");
        assert!(!regions[1].is_e1);
        assert!(!regions[2].is_e1);
        assert!(regions[3].is_e1, "MINUS E1 = highest start");
    }
}
