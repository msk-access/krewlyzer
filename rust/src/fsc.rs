//! Fragment Size Coverage (FSC) and Ratio (FSR) calculation
//!
//! Counts fragments in genomic bins by 6 biologically-meaningful size categories:
//! - Ultra-short: 65-100bp (di-nucleosomal debris, early apoptosis markers)
//! - Core-short: 101-149bp (sub-nucleosomal, associated with specific chromatin states)
//! - Mono-nucleosomal: 150-220bp (classic cfDNA peak, nucleosome-protected)
//! - Di-nucleosomal: 221-260bp (two nucleosomes, transitional)
//! - Long: 261-400bp (multi-nucleosomal, associated with necrosis)
//! - Ultra-long: 401-1000bp (necrosis, fetal cfDNA, apoptosis late-stage)
//!
//! These 6 channels are non-overlapping for ML feature generation.
//!
//! Genome-bin and gene/region aggregation share one classifier
//! ([`SizeChannel::of`]). They used to carry independent copies of the
//! boundaries and had drifted: the gene path put 150-259bp in `mono_nucl` and
//! 260-399bp in `di_nucl`, so a column named `di_nucl` held the genome path's
//! `long` range, and gene `long` held `ultra_long`. Same names, different
//! meaning, in tables a consumer reads side by side.

use std::path::Path;
use std::io::BufRead;
use std::fs::File;
use anyhow::Result;

use pyo3::prelude::*;
use numpy::{PyArray1, IntoPyArray};
use coitrees::{COITree, IntervalNode, IntervalTree};
use std::sync::Arc;
use std::collections::HashMap;
use std::io::Write;

use crate::bed::{Region, Fragment, ChromosomeMap};
use crate::engine::{FragmentConsumer, FragmentAnalyzer};
use crate::gc_correction::CorrectionFactors;

/// Smallest and largest fragment length carried by any channel.
///
/// Both the genome-bin and the gene/region path filter to this range before
/// classifying, so the six channels partition every counted fragment.
pub const MIN_FRAGMENT_LEN: u64 = 65;
pub const MAX_FRAGMENT_LEN: u64 = 1000;

/// The six non-overlapping size channels.
///
/// Bounds live here and nowhere else. Duplicating them is what let the
/// gene-level bands drift out of step with the genome-level ones.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SizeChannel {
    /// 65-100bp: di-nucleosomal debris, early apoptosis markers
    UltraShort,
    /// 101-149bp: sub-nucleosomal, associated with specific chromatin states
    CoreShort,
    /// 150-220bp: classic cfDNA peak, nucleosome-protected
    MonoNucl,
    /// 221-260bp: two nucleosomes, transitional
    DiNucl,
    /// 261-400bp: multi-nucleosomal, associated with necrosis
    Long,
    /// 401-1000bp: necrosis, fetal cfDNA, late-stage apoptosis
    UltraLong,
}

impl SizeChannel {
    /// Classify a fragment length.
    ///
    /// Returns `None` outside `MIN_FRAGMENT_LEN..=MAX_FRAGMENT_LEN`; callers
    /// filter first, so `None` means the caller's filter and these bounds have
    /// diverged rather than that the fragment is merely uninteresting.
    #[inline]
    pub fn of(len: u64) -> Option<Self> {
        match len {
            0..=64 => None,
            65..=100 => Some(Self::UltraShort),
            101..=149 => Some(Self::CoreShort),
            150..=220 => Some(Self::MonoNucl),
            221..=260 => Some(Self::DiNucl),
            261..=400 => Some(Self::Long),
            401..=1000 => Some(Self::UltraLong),
            _ => None,
        }
    }
}

/// Result of FSC/FSR calculation for a single bin
///
/// Contains 6 non-overlapping fragment size channels optimized for ML features.
#[derive(Debug, Clone, Default)]
pub struct BinResult {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    // 6 biologically-meaningful channels (non-overlapping, weighted counts)
    pub ultra_short_count: f64,  // 65-100bp: di-nucleosomal debris
    pub core_short_count: f64,   // 101-149bp: sub-nucleosomal
    pub mono_nucl_count: f64,    // 150-220bp: mono-nucleosomal (classic cfDNA peak)
    pub di_nucl_count: f64,      // 221-260bp: di-nucleosomal
    pub long_count: f64,         // 261-400bp: multi-nucleosomal
    pub ultra_long_count: f64,   // 401-1000bp: necrosis, fetal, apoptosis late-stage
    pub total_count: f64,        // All fragments 65-1000bp
    /// Counted into `total_count` but into no channel.
    ///
    /// Unreachable while the length filter and [`SizeChannel::of`] agree.
    /// Accumulated rather than logged at the point of detection: this runs
    /// per fragment inside the rayon loop, and logging there routes through
    /// pyo3-log into Python, which deadlocks under the GIL when it fires at
    /// volume. Reported once, with a count, at write time.
    pub unclassified_count: f64,
    pub mean_gc: f64,
    // Internal use for mean calc
    pub gc_sum: f64,
    pub gc_count: f64,
}

/// Parse a BED file (plain or BGZF-compressed) to get regions (bins)
pub fn parse_bin_file(bin_path: &Path) -> Result<Vec<Region>> {
    use crate::bed;
    let reader = bed::get_reader(bin_path)?;
    
    let mut regions = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(region) = parse_bed_line_to_region(&line) {
            regions.push(region);
        }
    }
    
    Ok(regions)
}

/// Parse a single BED line to a Region
fn parse_bed_line_to_region(line: &str) -> Option<Region> {
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

/// FscConsumer for the Unified Engine
/// 
/// When target_regions are provided (panel data like MSK-ACCESS),
/// counts are split into off-target (primary) and on-target (secondary) outputs.
/// Off-target data is used for biomarker calculation (unbiased by PCR).
/// On-target data is only useful for CNV calling.
#[derive(Clone)]
pub struct FscConsumer {
    // Read-only shared state
    trees: Arc<HashMap<u32, COITree<usize, u32>>>, // ChromID -> (IntervalTree<Data=usize, Coords=u32>)
    factors: Option<Arc<CorrectionFactors>>,
    
    // Target regions for on/off-target split (panel data)
    target_tree: Option<Arc<HashMap<u32, COITree<(), u32>>>>,
    
    // Thread-local state
    // Off-target counts (default output - unbiased)
    counts: Vec<BinResult>,
    // On-target counts (only populated if target_tree is Some)
    on_target_counts: Vec<BinResult>,
}

impl FscConsumer {
    pub fn new(
        regions: &[Region], 
        chrom_map: &mut ChromosomeMap, 
        factors: Option<Arc<CorrectionFactors>>,
        target_regions: Option<Arc<HashMap<u32, COITree<(), u32>>>>
    ) -> Self {
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        let mut counts = Vec::with_capacity(regions.len());
        let mut on_target_counts = Vec::with_capacity(regions.len());
        
        for (i, region) in regions.iter().enumerate() {
            // Init count for this region (off-target)
            counts.push(BinResult {
                chrom: region.chrom.clone(),
                start: region.start,
                end: region.end,
                ..Default::default()
            });
            // Init count for on-target
            on_target_counts.push(BinResult {
                chrom: region.chrom.clone(),
                start: region.start,
                end: region.end,
                ..Default::default()
            });
            
            // Normalize chromosome
            let chrom_norm = region.chrom.trim_start_matches("chr");
            let chrom_id = chrom_map.get_id(chrom_norm);
            
            // COITree uses closed intervals [start, end]. BED uses semi-open [start, end).
            // Convert to closed by subtracting 1 from end.
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
            factors,
            target_tree: target_regions,
            counts,
            on_target_counts,
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
    
    /// Write counts to file
    fn write_counts(&self, counts: &[BinResult], path: &Path) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = std::io::BufWriter::new(file);
        
        writeln!(writer, "chrom\tstart\tend\tultra_short\tcore_short\tmono_nucl\tdi_nucl\tlong\tultra_long\ttotal\tmean_gc")?;
        
        for bin in counts {
            let mean_gc = if bin.gc_count > 0.0 { bin.gc_sum / bin.gc_count } else { 0.0 };
            writeln!(writer, "{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.4}",
                bin.chrom, bin.start, bin.end,
                bin.ultra_short_count, bin.core_short_count, bin.mono_nucl_count, 
                bin.di_nucl_count, bin.long_count, bin.ultra_long_count, bin.total_count, mean_gc
            )?;
        }
        warn_if_unclassified(
            &format!("FSC ({})", path.display()),
            counts.iter().map(|b| b.unclassified_count).sum(),
            counts.iter().map(|b| b.total_count).sum(),
        );
        Ok(())
    }
    
    pub fn write_output(&self, path: &Path) -> Result<()> {
        use log::info;
        
        // Log summary
        let total_off: f64 = self.counts.iter().map(|b| b.total_count).sum();
        let total_on: f64 = self.on_target_counts.iter().map(|b| b.total_count).sum();
        
        info!("FSC: {} bins, {:.0} off-target frags, {:.0} on-target frags", 
            self.counts.len(), total_off, total_on);
        
        // Write off-target (main output - unbiased for biomarkers)
        self.write_counts(&self.counts, path)?;
        
        // Write on-target if targets were provided and there's data
        if self.target_tree.is_some() && total_on > 0.0 {
            let stem = path.file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("output");
            let on_path = path.with_file_name(format!("{}.ontarget.tsv", stem));
            self.write_counts(&self.on_target_counts, &on_path)?;
            info!("FSC on-target: {:?}", on_path);
        }
        
        Ok(())
    }
}

impl FragmentConsumer for FscConsumer {
    fn name(&self) -> &str {
        "FSC"
    }

    fn consume(&mut self, fragment: &Fragment) {
        // Check if this fragment is on-target (for routing to correct count vector)
        let on_target = self.is_on_target(fragment);
        
        if let Some(tree) = self.trees.get(&fragment.chrom_id) {
            let start = fragment.start as u32;
            let end = fragment.end as u32;
            let end_closed = if end > start { end - 1 } else { start };
            
            // Query for overlapping bins
            tree.query(start as i32, end_closed as i32, |node| {
                let bin_idx = node.metadata.to_owned();
                
                // Accept fragments in 65-1000bp range (extended for ultra-long analysis)
                if fragment.length >= MIN_FRAGMENT_LEN && fragment.length <= MAX_FRAGMENT_LEN {
                    let gc_pct = (fragment.gc * 100.0).round() as u8;
                    let weight = if let Some(ref factors) = self.factors {
                         factors.get_factor(fragment.length, gc_pct)
                    } else {
                         1.0
                    };

                    // Route to correct vector (on-target or off-target)
                    // SAFETY: bin_idx comes from initialization, guaranteed valid.
                    unsafe {
                        let res = if on_target {
                            self.on_target_counts.get_unchecked_mut(bin_idx)
                        } else {
                            self.counts.get_unchecked_mut(bin_idx)
                        };
                        
                        res.total_count += weight;
                        res.gc_sum += fragment.gc * weight;
                        res.gc_count += weight;
                        
                        // 6 non-overlapping channels for ML features
                        match SizeChannel::of(fragment.length) {
                            Some(SizeChannel::UltraShort) => res.ultra_short_count += weight,
                            Some(SizeChannel::CoreShort) => res.core_short_count += weight,
                            Some(SizeChannel::MonoNucl) => res.mono_nucl_count += weight,
                            Some(SizeChannel::DiNucl) => res.di_nucl_count += weight,
                            Some(SizeChannel::Long) => res.long_count += weight,
                            Some(SizeChannel::UltraLong) => res.ultra_long_count += weight,
                            // Unreachable behind the length filter above. Landing
                            // in total but in no channel would silently break the
                            // partition every ratio depends on; see the field docs
                            // for why this counts instead of logging here.
                            None => res.unclassified_count += weight,
                        }
                    }
                }
            });
        }
    }

    fn merge(&mut self, other: Self) {
        // Merge off-target counts
        for (i, other_bin) in other.counts.into_iter().enumerate() {
            let my_bin = &mut self.counts[i];
            my_bin.ultra_short_count += other_bin.ultra_short_count;
            my_bin.core_short_count += other_bin.core_short_count;
            my_bin.mono_nucl_count += other_bin.mono_nucl_count;
            my_bin.di_nucl_count += other_bin.di_nucl_count;
            my_bin.long_count += other_bin.long_count;
            my_bin.ultra_long_count += other_bin.ultra_long_count;
            my_bin.unclassified_count += other_bin.unclassified_count;
            my_bin.total_count += other_bin.total_count;
            my_bin.gc_sum += other_bin.gc_sum;
            my_bin.gc_count += other_bin.gc_count;
        }
        
        // Merge on-target counts
        for (i, other_bin) in other.on_target_counts.into_iter().enumerate() {
            let my_bin = &mut self.on_target_counts[i];
            my_bin.ultra_short_count += other_bin.ultra_short_count;
            my_bin.core_short_count += other_bin.core_short_count;
            my_bin.mono_nucl_count += other_bin.mono_nucl_count;
            my_bin.di_nucl_count += other_bin.di_nucl_count;
            my_bin.long_count += other_bin.long_count;
            my_bin.ultra_long_count += other_bin.ultra_long_count;
            my_bin.unclassified_count += other_bin.unclassified_count;
            my_bin.total_count += other_bin.total_count;
            my_bin.gc_sum += other_bin.gc_sum;
            my_bin.gc_count += other_bin.gc_count;
        }
    }
}

// Re-implement the original function logic using the legacy approach for now 
// (or delete if replaced, but let's keep it for compatibility until switched).
// ACTUALLY, I will add `count_fragments_unified` and call IT if a flag is set?
// Or just replace the implementation of `count_fragments_by_bins`.
// I will replace the implementation of `count_fragments_sequential` to use the legacy logic 
// (it is already there).
// I won't touch the original logic in this file, I just appended the new logic.

// ... (Original logic omitted for brevity in prompt, but I keep it in file)








/// Python-exposed function to calculate FSC/FSR with 5 ML-ready channels
/// Returns: (ultra_shorts, core_shorts, mono_nucls, di_nucls, longs, totals, gcs)
#[pyfunction]
#[pyo3(signature = (bedgz_path, bin_path))]
pub fn count_fragments_by_bins(
    py: Python<'_>,
    bedgz_path: &str,
    bin_path: &str,
) -> PyResult<(
    Py<PyArray1<u32>>,  // ultra-short counts (65-100bp)
    Py<PyArray1<u32>>,  // core-short counts (101-149bp)
    Py<PyArray1<u32>>,  // mono-nucleosomal counts (150-220bp)
    Py<PyArray1<u32>>,  // di-nucleosomal counts (221-260bp)
    Py<PyArray1<u32>>,  // long counts (261-400bp)
    Py<PyArray1<u32>>,  // total counts (65-400bp)
    Py<PyArray1<f64>>,  // mean GC
)> {
    let bed_path = Path::new(bedgz_path);
    let bin_path_p = Path::new(bin_path);
    
    // 1. Prepare
    let regions = parse_bin_file(bin_path_p)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        
    let mut chrom_map = ChromosomeMap::new();
    let consumer = FscConsumer::new(&regions, &mut chrom_map, None, None);  // No target_regions for legacy API
    
    // 2. Run Engine (Using logic similar to Unified Pipeline but for single consumer)
    let engine = FragmentAnalyzer::new(consumer, 100_000); // 100k chunks
    let mut final_consumer = engine.process_file(bed_path, &mut chrom_map, false)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        
    // 3. Convert Results
    // Calculate mean GC (Logic ported from legacy struct handling)
    for bin in &mut final_consumer.counts {
        bin.mean_gc = if bin.gc_count > 0.0 {
            bin.gc_sum / bin.gc_count
        } else {
            f64::NAN
        };
    }
    
    let results = final_consumer.counts;
    let ultra_shorts: Vec<u32> = results.iter().map(|r| r.ultra_short_count as u32).collect();
    let core_shorts: Vec<u32> = results.iter().map(|r| r.core_short_count as u32).collect();
    let mono_nucls: Vec<u32> = results.iter().map(|r| r.mono_nucl_count as u32).collect();
    let di_nucls: Vec<u32> = results.iter().map(|r| r.di_nucl_count as u32).collect();
    let longs: Vec<u32> = results.iter().map(|r| r.long_count as u32).collect();
    let totals: Vec<u32> = results.iter().map(|r| r.total_count as u32).collect();
    let gcs: Vec<f64> = results.iter().map(|r| r.mean_gc).collect();
    
    Ok((
        ultra_shorts.into_pyarray(py).into(),
        core_shorts.into_pyarray(py).into(),
        mono_nucls.into_pyarray(py).into(),
        di_nucls.into_pyarray(py).into(),
        longs.into_pyarray(py).into(),
        totals.into_pyarray(py).into(),
        gcs.into_pyarray(py).into(),
    ))
}

// =============================================================================
// GENE-LEVEL FSC AGGREGATION
// =============================================================================

use std::io::BufWriter;
use std::path::PathBuf;

/// Per-gene or per-region aggregated counts
#[derive(Debug, Clone, Default)]
struct GeneResult {
    gene: String,
    n_regions: usize,
    total_bp: u64,
    ultra_short: f64,
    core_short: f64,
    mono_nucl: f64,
    di_nucl: f64,
    long: f64,
    ultra_long: f64,
    total: f64,
    /// See [`BinResult::unclassified_count`] -- counted, not logged in the loop.
    unclassified: f64,
}

impl GeneResult {
    fn add_fragment(&mut self, len: u64, weight: f64) {
        self.total += weight;
        match SizeChannel::of(len) {
            Some(SizeChannel::UltraShort) => self.ultra_short += weight,
            Some(SizeChannel::CoreShort) => self.core_short += weight,
            Some(SizeChannel::MonoNucl) => self.mono_nucl += weight,
            Some(SizeChannel::DiNucl) => self.di_nucl += weight,
            Some(SizeChannel::Long) => self.long += weight,
            Some(SizeChannel::UltraLong) => self.ultra_long += weight,
            None => self.unclassified += weight,
        }
    }

    fn merge(&mut self, other: &Self) {
        self.ultra_short += other.ultra_short;
        self.core_short += other.core_short;
        self.mono_nucl += other.mono_nucl;
        self.di_nucl += other.di_nucl;
        self.long += other.long;
        self.ultra_long += other.ultra_long;
        self.total += other.total;
        self.unclassified += other.unclassified;
    }
}

/// Warn once if any fragment was counted into `total` but into no channel.
///
/// The six ratios are computed as `channel / total`, so an unclassified
/// fragment makes them sum to less than 1 without any single column looking
/// wrong. Reported here rather than per fragment: this is one call per output
/// file instead of one per fragment inside the rayon loop, where routing
/// through pyo3-log into Python deadlocks under the GIL.
fn warn_if_unclassified(scope: &str, unclassified: f64, total: f64) {
    if unclassified > 0.0 {
        log::warn!(
            "{}: {:.0} of {:.0} fragments matched no size channel despite passing \
             the {}-{}bp filter; the six ratios will not sum to 1. The length \
             filter and SizeChannel::of have diverged.",
            scope,
            unclassified,
            total,
            MIN_FRAGMENT_LEN,
            MAX_FRAGMENT_LEN
        );
    }
}

/// Gene region from BED file
#[derive(Debug, Clone)]
struct GeneRegion {
    chrom: String,
    start: u64,
    end: u64,
    gene: String,
    name: String,
    /// `+`, `-`, or `.` when the source carried no strand.
    strand: String,
    /// Overlaps the canonical transcript's exon 1. `None` means the source
    /// could not say -- a legacy gene BED, or a target BED parsed on the fly.
    /// Deliberately not conflated with `Some(false)`: absent is not "no".
    is_e1: Option<bool>,
    /// Overlaps another basic protein-coding transcript's exon 1.
    is_alt_e1: Option<bool>,
}

/// Parse the gene BED written by `fsc_processor.aggregate_by_gene`:
/// chrom, start, end, gene, [name], [strand], [is_e1], [is_alt_e1].
///
/// Columns six to eight are what let E1 selection be strand-aware. A 5-column
/// file still parses; its flags are `None`, and `filter_fsc_to_e1` falls back
/// to the coordinate heuristic with a warning rather than silently claiming an
/// E1 it cannot know.
fn parse_gene_bed(path: &Path) -> Result<Vec<GeneRegion>> {
    let reader = crate::bed::get_reader(path)?;
    let mut regions = Vec::new();

    // `.` is unknown, `1`/`0` are known. Anything else is malformed and is
    // treated as unknown rather than guessed at.
    fn flag(field: Option<&&str>) -> Option<bool> {
        match field {
            Some(&"1") => Some(true),
            Some(&"0") => Some(false),
            _ => None,
        }
    }

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }
        
        let chrom = fields[0].to_string();
        let start: u64 = fields[1].parse().unwrap_or(0);
        let end: u64 = fields[2].parse().unwrap_or(0);
        let gene = fields[3].to_string();
        let name = if fields.len() > 4 { fields[4].to_string() } else { gene.clone() };
        let strand = fields.get(5).map(|s| s.to_string()).unwrap_or_else(|| ".".to_string());

        regions.push(GeneRegion {
            chrom,
            start,
            end,
            gene,
            name,
            strand,
            is_e1: flag(fields.get(6)),
            is_alt_e1: flag(fields.get(7)),
        });
    }
    
    log::info!("Parsed {} gene regions from {:?}", regions.len(), path);
    Ok(regions)
}

/// Gene FSC Consumer for fragment aggregation
#[derive(Clone)]
pub struct GeneFscConsumer {
    // Interval tree: ChromID -> (start, end) -> region_idx
    trees: Arc<HashMap<u32, COITree<usize, u32>>>,
    // All regions
    regions: Arc<Vec<GeneRegion>>,
    // GC correction factors
    factors: Option<Arc<CorrectionFactors>>,
    // Per-region counts (thread-local, merged after)
    region_counts: Vec<GeneResult>,
    // Total fragments for normalization
    total_fragments: u64,
}

impl GeneFscConsumer {
    #[allow(private_interfaces)]
    pub fn new(
        regions: Vec<GeneRegion>,
        chrom_map: &mut ChromosomeMap,
        factors: Option<Arc<CorrectionFactors>>,
    ) -> Self {
        let mut nodes_by_chrom: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        let mut region_counts = Vec::with_capacity(regions.len());
        
        for (i, region) in regions.iter().enumerate() {
            region_counts.push(GeneResult {
                gene: region.gene.clone(),
                n_regions: 1,
                total_bp: region.end - region.start,
                ..Default::default()
            });
            // Normalize chromosome (strip chr prefix) consistent with engine.rs
            let chrom_norm = region.chrom.trim_start_matches("chr");
            let chrom_id = chrom_map.get_id(chrom_norm);
            let start = region.start as i32;
            let end = (region.end.saturating_sub(1)) as i32;
            
            nodes_by_chrom.entry(chrom_id).or_default().push(
                IntervalNode::new(start, end, i)
            );
        }
        
        let mut trees = HashMap::new();
        for (chrom_id, nodes) in nodes_by_chrom {
            trees.insert(chrom_id, COITree::new(&nodes));
        }
        
        log::info!("GeneFSC: Built trees for {} regions", regions.len());
        
        Self {
            trees: Arc::new(trees),
            regions: Arc::new(regions),
            factors,
            region_counts,
            total_fragments: 0,
        }
    }
    
    /// Aggregate per-region counts to per-gene
    fn aggregate_to_genes(&self) -> HashMap<String, GeneResult> {
        let mut gene_results: HashMap<String, GeneResult> = HashMap::new();
        
        for (idx, region) in self.regions.iter().enumerate() {
            let rc = &self.region_counts[idx];
            
            gene_results.entry(region.gene.clone())
                .and_modify(|g| {
                    g.n_regions += 1;
                    g.total_bp += region.end - region.start;
                    g.merge(rc);
                })
                .or_insert_with(|| GeneResult {
                    gene: region.gene.clone(),
                    n_regions: 1,
                    total_bp: region.end - region.start,
                    ultra_short: rc.ultra_short,
                    core_short: rc.core_short,
                    mono_nucl: rc.mono_nucl,
                    di_nucl: rc.di_nucl,
                    long: rc.long,
                    ultra_long: rc.ultra_long,
                    total: rc.total,
                    unclassified: rc.unclassified,
                });
        }
        
        gene_results
    }
    
    /// Write gene-level output to TSV
    pub fn write_gene_output(&self, path: &Path) -> Result<()> {
        let gene_results = self.aggregate_to_genes();
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        
        // Header
        writeln!(writer, "gene\tn_regions\ttotal_bp\tultra_short\tcore_short\tmono_nucl\tdi_nucl\tlong\tultra_long\ttotal\tultra_short_ratio\tcore_short_ratio\tmono_nucl_ratio\tdi_nucl_ratio\tlong_ratio\tultra_long_ratio\tnormalized_depth")?;
        
        // Sort by gene name
        let mut genes: Vec<_> = gene_results.keys().collect();
        genes.sort();
        
        for gene in genes {
            let g = &gene_results[gene];
            let total = g.total.max(1e-9);
            
            // Ratios
            let us_r = g.ultra_short / total;
            let cs_r = g.core_short / total;
            let mn_r = g.mono_nucl / total;
            let dn_r = g.di_nucl / total;
            let lg_r = g.long / total;
            let ul_r = g.ultra_long / total;

            // Normalized depth (RPKM-like)
            let total_frags = self.total_fragments.max(1) as f64;
            let norm_depth = if g.total_bp > 0 && total_frags > 0.0 {
                (g.total * 1e9) / (g.total_bp as f64 * total_frags)
            } else {
                0.0
            };
            
            writeln!(writer, "{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                g.gene, g.n_regions, g.total_bp,
                g.ultra_short, g.core_short, g.mono_nucl, g.di_nucl, g.long, g.ultra_long, g.total,
                us_r, cs_r, mn_r, dn_r, lg_r, ul_r, norm_depth
            )?;
        }
        
        warn_if_unclassified(
            "GeneFSC (gene)",
            gene_results.values().map(|g| g.unclassified).sum(),
            gene_results.values().map(|g| g.total).sum(),
        );
        log::info!("GeneFSC: Wrote {} genes to {:?}", gene_results.len(), path);
        Ok(())
    }
    /// Write region-level (per-exon) output to TSV
    pub fn write_region_output(&self, path: &Path) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        
        // Header
        // strand/is_e1/is_alt_e1 are carried through so E1 selection happens
        // downstream on a fact rather than a coordinate guess.
        writeln!(writer, "chrom\tstart\tend\tgene\tregion_name\tregion_bp\tstrand\tis_e1\tis_alt_e1\tultra_short\tcore_short\tmono_nucl\tdi_nucl\tlong\tultra_long\ttotal\tultra_short_ratio\tcore_short_ratio\tmono_nucl_ratio\tdi_nucl_ratio\tlong_ratio\tultra_long_ratio\tnormalized_depth")?;
        
        let total_frags = self.total_fragments.max(1) as f64;
        
        for (idx, region) in self.regions.iter().enumerate() {
            let rc = &self.region_counts[idx];
            let region_bp = region.end - region.start;
            let total = rc.total.max(1e-9);
            
            // Ratios
            let us_r = rc.ultra_short / total;
            let cs_r = rc.core_short / total;
            let mn_r = rc.mono_nucl / total;
            let dn_r = rc.di_nucl / total;
            let lg_r = rc.long / total;
            let ul_r = rc.ultra_long / total;

            // Normalized depth (RPKM-like)
            let norm_depth = if region_bp > 0 && total_frags > 0.0 {
                (rc.total * 1e9) / (region_bp as f64 * total_frags)
            } else {
                0.0
            };
            
            // "." for an unknown flag, never "0": a legacy gene BED cannot
            // say whether a region is E1, and writing 0 would assert it is not.
            let fmt_flag = |v: Option<bool>| match v {
                Some(true) => "1",
                Some(false) => "0",
                None => ".",
            };
            writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                region.chrom, region.start, region.end, region.gene, region.name, region_bp,
                region.strand, fmt_flag(region.is_e1), fmt_flag(region.is_alt_e1),
                rc.ultra_short, rc.core_short, rc.mono_nucl, rc.di_nucl, rc.long, rc.ultra_long, rc.total,
                us_r, cs_r, mn_r, dn_r, lg_r, ul_r, norm_depth
            )?;
        }
        
        warn_if_unclassified(
            "GeneFSC (region)",
            self.region_counts.iter().map(|r| r.unclassified).sum(),
            self.region_counts.iter().map(|r| r.total).sum(),
        );
        log::info!("GeneFSC: Wrote {} regions to {:?}", self.regions.len(), path);
        Ok(())
    }
}

impl FragmentConsumer for GeneFscConsumer {
    fn name(&self) -> &str {
        "GeneFSC"
    }
    
    fn consume(&mut self, fragment: &Fragment) {
        // Filter by length; same bounds the channels are defined over, so the
        // six channels partition exactly what `total` counts.
        if fragment.length < MIN_FRAGMENT_LEN || fragment.length > MAX_FRAGMENT_LEN {
            return;
        }
        
        self.total_fragments += 1;
        
        if let Some(tree) = self.trees.get(&fragment.chrom_id) {
            // Use fragment midpoint for region assignment
            let mid = ((fragment.start + fragment.end) / 2) as i32;
            
            // Query overlapping regions
            tree.query(mid, mid, |node| {
                let region_idx = node.metadata.to_owned();
                
                // Calculate GC weight
                let gc_pct = (fragment.gc * 100.0).round() as u8;
                let weight = self.factors.as_ref()
                    .map(|f| f.get_factor(fragment.length, gc_pct))
                    .unwrap_or(1.0);
                
                self.region_counts[region_idx].add_fragment(fragment.length, weight);
            });
        }
    }
    
    fn merge(&mut self, other: Self) {
        self.total_fragments += other.total_fragments;
        for (i, rc) in other.region_counts.iter().enumerate() {
            self.region_counts[i].merge(rc);
        }
    }
}

/// Aggregate fragment counts by gene from BED.gz file.
/// 
/// Reads fragments from BED.gz, assigns to gene regions, counts by size channel,
/// applies GC correction if provided, and outputs gene-level or region-level FSC.
/// 
/// # Arguments
/// * `bed_path` - Path to fragment BED.gz file
/// * `gene_bed_path` - Path to gene regions BED (chrom, start, end, gene, [name])
/// * `output_path` - Path to output TSV
/// * `gc_factors_path` - Optional path to GC correction factors TSV
/// * `aggregate_by` - "gene" for gene-level aggregation, "region" for per-exon output
/// 
/// # Returns
/// Number of rows in output (genes or regions)
#[pyfunction]
#[pyo3(signature = (bed_path, gene_bed_path, output_path, gc_factors_path=None, aggregate_by="gene"))]
pub fn aggregate_by_gene(
    bed_path: PathBuf,
    gene_bed_path: PathBuf,
    output_path: PathBuf,
    gc_factors_path: Option<PathBuf>,
    aggregate_by: &str,
) -> PyResult<usize> {
    log::info!("GeneFSC: Starting {} aggregation - {:?}", aggregate_by, bed_path);
    
    // 1. Load gene regions
    let gene_regions = parse_gene_bed(&gene_bed_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(
            format!("Failed to parse gene BED: {}", e)
        ))?;
    
    if gene_regions.is_empty() {
        log::warn!("GeneFSC: No gene regions found");
        return Ok(0);
    }
    
    let n_regions = gene_regions.len();
    
    // 2. Load GC correction factors — Parquet-first auto-detection via load()
    let factors = if let Some(ref path) = gc_factors_path {
        // CorrectionFactors::load() tries .parquet first, then TSV/legacy CSV
        match CorrectionFactors::load(path) {
            Ok(f) => {
                log::info!("GeneFSC: Loaded {} GC correction factor bins from {:?}", f.data.len(), path);
                Some(Arc::new(f))
            }
            Err(e) => {
                log::warn!("GeneFSC: Failed to load GC factors from {:?}: {}", path, e);
                None
            }
        }
    } else {
        None
    };
    
    // 3. Create consumer and process
    let mut chrom_map = ChromosomeMap::new();
    let consumer = GeneFscConsumer::new(gene_regions, &mut chrom_map, factors);
    
    let engine = FragmentAnalyzer::new(consumer, 100_000);
    let final_consumer = engine.process_file(&bed_path, &mut chrom_map, false)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    
    // 4. Write output based on aggregate_by mode
    let output_count = if aggregate_by == "region" {
        final_consumer.write_region_output(&output_path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        n_regions
    } else {
        let gene_count = final_consumer.aggregate_to_genes().len();
        final_consumer.write_gene_output(&output_path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        gene_count
    };
    
    log::info!("GeneFSC: Complete - {} {}, {} fragments", 
        output_count, aggregate_by, final_consumer.total_fragments);
    
    Ok(output_count)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every length in the accepted range lands in exactly one channel.
    ///
    /// This is the property the ratio columns depend on: `channel / total` is
    /// only meaningful if the channels partition what `total` counts.
    #[test]
    fn every_accepted_length_has_exactly_one_channel() {
        for len in MIN_FRAGMENT_LEN..=MAX_FRAGMENT_LEN {
            assert!(
                SizeChannel::of(len).is_some(),
                "{len}bp is inside the accepted range but matches no channel"
            );
        }
    }

    /// Lengths the callers filter out must not be silently absorbed.
    #[test]
    fn lengths_outside_the_range_have_no_channel() {
        for len in [0, 1, 64, MAX_FRAGMENT_LEN + 1, 5_000] {
            assert!(
                SizeChannel::of(len).is_none(),
                "{len}bp is outside the accepted range but was assigned a channel"
            );
        }
    }

    /// Boundaries, spelled out.
    ///
    /// Both sides of each cut are pinned, because an off-by-one here silently
    /// moves fragments between channels rather than failing.
    #[test]
    fn channel_boundaries_are_exact() {
        let cases = [
            (65, SizeChannel::UltraShort),
            (100, SizeChannel::UltraShort),
            (101, SizeChannel::CoreShort),
            (149, SizeChannel::CoreShort),
            (150, SizeChannel::MonoNucl),
            (220, SizeChannel::MonoNucl),
            (221, SizeChannel::DiNucl),
            (260, SizeChannel::DiNucl),
            (261, SizeChannel::Long),
            (400, SizeChannel::Long),
            (401, SizeChannel::UltraLong),
            (1000, SizeChannel::UltraLong),
        ];
        for (len, expected) in cases {
            assert_eq!(
                SizeChannel::of(len),
                Some(expected),
                "{len}bp landed in the wrong channel"
            );
        }
    }

    /// The regression this module's shared classifier exists to prevent.
    ///
    /// The gene path used to carry its own boundaries, and they had drifted:
    /// 150-259 went to `mono_nucl` and 260-399 to `di_nucl`, so gene `di_nucl`
    /// held the genome path's `long` range. These four lengths are exactly
    /// where the two disagreed.
    #[test]
    fn gene_and_genome_bands_agree_where_they_used_to_diverge() {
        let mut gene = GeneResult::default();
        for len in [100, 230, 300, 500] {
            gene.add_fragment(len, 1.0);
        }

        assert_eq!(gene.ultra_short, 1.0, "100bp belongs to ultra_short");
        assert_eq!(gene.core_short, 0.0, "100bp used to fall into core_short");
        assert_eq!(gene.di_nucl, 1.0, "230bp belongs to di_nucl");
        assert_eq!(gene.mono_nucl, 0.0, "230bp used to fall into mono_nucl");
        assert_eq!(gene.long, 1.0, "300bp belongs to long");
        assert_eq!(gene.ultra_long, 1.0, "500bp belongs to ultra_long");
        assert_eq!(gene.total, 4.0);
    }

    /// The six gene channels partition `total`, so the six ratios sum to 1.
    #[test]
    fn gene_channels_partition_the_total() {
        let mut gene = GeneResult::default();
        for len in (MIN_FRAGMENT_LEN..=MAX_FRAGMENT_LEN).step_by(7) {
            gene.add_fragment(len, 1.0);
        }
        let summed = gene.ultra_short
            + gene.core_short
            + gene.mono_nucl
            + gene.di_nucl
            + gene.long
            + gene.ultra_long;
        assert_eq!(
            summed, gene.total,
            "the six channels must sum to total, or every ratio is wrong"
        );
    }

    /// A fragment matching no channel is counted, not dropped and not logged.
    ///
    /// The loop cannot log here -- pyo3-log routes into Python and deadlocks
    /// under the GIL when it fires per fragment -- so the signal has to survive
    /// as a number until `warn_if_unclassified` reports it at write time.
    #[test]
    fn unclassified_fragments_are_counted_rather_than_dropped() {
        let mut gene = GeneResult::default();
        gene.add_fragment(30, 1.0); // below MIN_FRAGMENT_LEN
        gene.add_fragment(2_000, 1.0); // above MAX_FRAGMENT_LEN
        gene.add_fragment(185, 1.0); // ordinary

        assert_eq!(gene.unclassified, 2.0, "out-of-range weight must be retained");
        assert_eq!(gene.mono_nucl, 1.0);
        assert_eq!(gene.total, 3.0, "total counts everything handed to it");
    }

    /// `merge` must carry the sixth channel; forgetting it would lose
    /// ultra-long counts only when a gene spans more than one region.
    #[test]
    fn merge_accumulates_every_channel() {
        let mut a = GeneResult::default();
        let mut b = GeneResult::default();
        for len in [80, 120, 200, 240, 300, 600] {
            a.add_fragment(len, 1.0);
            b.add_fragment(len, 2.0);
        }
        a.merge(&b);

        for (name, value) in [
            ("ultra_short", a.ultra_short),
            ("core_short", a.core_short),
            ("mono_nucl", a.mono_nucl),
            ("di_nucl", a.di_nucl),
            ("long", a.long),
            ("ultra_long", a.ultra_long),
        ] {
            assert_eq!(value, 3.0, "{name} did not accumulate across merge");
        }
        assert_eq!(a.total, 18.0);
        assert_eq!(a.unclassified, 0.0, "every length used here is in range");

        // And the counter itself survives a merge, or a divergence detected in
        // one region would vanish when its gene rolled up.
        let mut c = GeneResult::default();
        let mut d = GeneResult::default();
        c.add_fragment(10, 1.0);
        d.add_fragment(10, 1.0);
        c.merge(&d);
        assert_eq!(c.unclassified, 2.0, "unclassified did not survive merge");
    }
}
