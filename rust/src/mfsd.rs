//! Mutant Fragment Size Distribution (mFSD) calculation
//!
//! Analyzes fragment size distributions around variant positions to detect
//! tumor-derived cfDNA. Supports all small variant types with CIGAR-aware extraction.
//!
//! ## Variant Types Supported
//! - **SNV**: Single nucleotide variants
//! - **MNV**: Multi-nucleotide variants
//! - **Insertion**: Pure insertions (A→ATG)
//! - **Deletion**: Pure deletions (ATG→A)
//! - **Complex**: Combined substitution + indel
//!
//! ## Fragment Classification (4-way)
//! - **REF**: Supports reference allele (healthy baseline)
//! - **ALT**: Supports alternate allele (tumor signal)
//! - **NonREF**: Non-REF, non-ALT (sequencing errors, subclones)
//! - **N**: Contains N at variant position (low quality)
//!
//! ## Output
//! - {sample}.mFSD.tsv: Summary statistics per variant (39 columns)
//! - {sample}.mFSD.distributions.tsv: Per-size counts (optional)

use pyo3::prelude::*;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::Cigar;
use rust_htslib::faidx;
use rayon::prelude::*;
use indicatif::{ProgressBar, ProgressStyle};
use std::time::Duration;
use log::{info, warn, debug};

// ============================================================================
// PHASE 1: Data Structures
// ============================================================================

/// Classification of variant type based on REF/ALT alleles
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VariantType {
    Snv,       // Single nucleotide variant (A>T)
    Mnv,       // Multi-nucleotide variant (AT>GC)
    Insertion, // Insertion (A>ATG)
    Deletion,  // Deletion (ATG>A)
    Complex,   // Complex (ATG>CT)
}

impl VariantType {
    pub fn as_str(&self) -> &'static str {
        match self {
            VariantType::Snv => "SNV",
            VariantType::Mnv => "MNV",
            VariantType::Insertion => "INS",
            VariantType::Deletion => "DEL",
            VariantType::Complex => "COMPLEX",
        }
    }
}

/// Classification of a fragment's allele support
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FragmentClass {
    Ref,    // Supports reference allele
    Alt,    // Supports alternate allele
    NonRef, // Non-REF, non-ALT, non-N (errors, subclones)
    N,      // Contains N at variant position
}

/// Parsed variant with type classification
#[derive(Debug, Clone)]
struct Variant {
    chrom: String,
    pos: i64,           // 0-based
    ref_allele: String,
    alt_allele: String,
    var_type: VariantType,
}

/// Result accumulator for a single variant - 4-way classification
/// Supports both unweighted counts and GC-weighted counts
#[derive(Debug, Clone, Default)]
struct VariantResult {
    ref_lengths: Vec<f64>,
    alt_lengths: Vec<f64>,
    nonref_lengths: Vec<f64>,
    n_lengths: Vec<f64>,
    // Weighted counts (sum of weights, not raw counts)
    ref_weighted: f64,
    alt_weighted: f64,
    nonref_weighted: f64,
    n_weighted: f64,
    // Debug: filter counters
    skipped_proper_pair: usize,
    skipped_size: usize,
    skipped_extreme: usize,  // Discordant reads with TLEN > 10000
    // Duplex tag tracking (for --duplex warning)
    duplex_tags_found: usize,
}

impl VariantResult {
    fn new() -> Self {
        Self::default()
    }
    
    fn total_count(&self) -> usize {
        self.ref_lengths.len() + self.alt_lengths.len() + 
        self.nonref_lengths.len() + self.n_lengths.len()
    }
    
    /// Add a fragment to the result with optional GC correction weight
    fn add_fragment(&mut self, class: FragmentClass, frag_len: f64, weight: f64) {
        match class {
            FragmentClass::Ref => {
                self.ref_lengths.push(frag_len);
                self.ref_weighted += weight;
            },
            FragmentClass::Alt => {
                self.alt_lengths.push(frag_len);
                self.alt_weighted += weight;
            },
            FragmentClass::NonRef => {
                self.nonref_lengths.push(frag_len);
                self.nonref_weighted += weight;
            },
            FragmentClass::N => {
                self.n_lengths.push(frag_len);
                self.n_weighted += weight;
            },
        }
    }
}

/// Compute GC content from a sequence
fn compute_gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.5; // Default
    }
    let gc_count = seq.iter().filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c').count();
    gc_count as f64 / seq.len() as f64
}

// ============================================================================
// DUPLEX CONSENSUS WEIGHTING
// ============================================================================
//
// Duplex sequencing produces high-confidence consensus reads from multiple
// raw observations. The family size (number of raw reads) indicates confidence.
// We apply log-weighting to prioritize high-confidence duplex families.
//
// Supported formats:
// 1. fgbio/picard (XS2): Uses SAM tags aD/bD/cD for depths
//    - cD = duplex consensus depth (primary weight source)
//    - See: https://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html
//
// 2. Marianas (XS1): Encodes family size in read name
//    - Format: Marianas:UMI:chr:start:posCount:negCount:chr2:start2:pos2:neg2
//    - Family size = posCount + negCount
//    - See: https://cmo-ci.gitbook.io/marianas/read-name-information
// ============================================================================

/// Get duplex consensus weight from BAM record
/// 
/// Returns a log-scaled weight based on duplex family size.
/// Higher family size = more sequencing evidence = higher weight.
/// 
/// Weight formula: ln(family_size).max(1.0)
/// - Family size 1: weight 1.0
/// - Family size 3: weight 1.1
/// - Family size 10: weight 2.3
/// - Family size 50: weight 3.9
fn get_duplex_weight(record: &bam::Record) -> f64 {
    use rust_htslib::bam::record::Aux;
    
    // Method 1: fgbio/picard cD tag (duplex consensus depth)
    // cD = maximum depth of raw reads at any point in the duplex consensus
    if let Ok(aux) = record.aux(b"cD") {
        let depth = match aux {
            Aux::I8(v) => v as f64,
            Aux::U8(v) => v as f64,
            Aux::I16(v) => v as f64,
            Aux::U16(v) => v as f64,
            Aux::I32(v) => v as f64,
            Aux::U32(v) => v as f64,
            _ => 0.0,
        };
        if depth > 0.0 {
            return depth.ln().max(1.0);
        }
    }
    
    // Method 2: Marianas read name format
    // Format: Marianas:UMI:chr:start:posCount:negCount:chr2:start2:pos2:neg2
    if let Ok(qname) = std::str::from_utf8(record.qname()) {
        if qname.starts_with("Marianas:") {
            if let Some(family_size) = parse_marianas_family(qname) {
                if family_size > 0 {
                    return (family_size as f64).ln().max(1.0);
                }
            }
        }
    }
    
    1.0  // Default: no duplex weighting
}

/// Parse Marianas read name to extract family size
/// 
/// Read name format: Marianas:UMI:contig:start:posCount:negCount:contig2:start2:pos2:neg2
/// Family size = posCount + negCount (for read1, fields 4 and 5)
/// 
/// Example: "Marianas:ACT+TTA:2:48033828:4:3:2:48033899:4:3"
/// → posCount=4, negCount=3, family_size=7
fn parse_marianas_family(qname: &str) -> Option<u32> {
    let parts: Vec<&str> = qname.split(':').collect();
    if parts.len() >= 6 {
        let pos_count: u32 = parts[4].parse().ok()?;
        let neg_count: u32 = parts[5].parse().ok()?;
        Some(pos_count + neg_count)
    } else {
        None
    }
}


/// Fetch sequence from reference FASTA for a genomic region
fn fetch_reference_gc(
    fasta: &faidx::Reader,
    chrom: &str,
    start: u64,
    end: u64,
) -> Option<f64> {
    // Try with and without chr prefix
    let chroms_to_try = if chrom.starts_with("chr") {
        vec![chrom.to_string(), chrom.trim_start_matches("chr").to_string()]
    } else {
        vec![chrom.to_string(), format!("chr{}", chrom)]
    };
    
    for chr_name in chroms_to_try {
        if let Ok(seq) = fasta.fetch_seq(&chr_name, start as usize, end as usize) {
            return Some(compute_gc_content(&seq));
        }
    }
    None
}

// ============================================================================
// PHASE 1b: MAF Header Parsing & Input Validation
// ============================================================================

/// Column index map resolved from MAF header.
///
/// MAF files can have variable column ordering (e.g., cBioPortal adds a
/// `Consequence` column at index 8, shifting all subsequent indices).
/// This struct resolves required columns by name from the header line,
/// making the parser resilient to any column layout.
#[derive(Debug, Clone)]
struct MafColumnMap {
    chrom: usize,
    pos: usize,
    ref_allele: usize,
    alt_allele: usize,
    max_idx: usize, // Maximum required index (for bounds checking)
}

/// Parse MAF header line to resolve column indices by name.
///
/// Required columns: Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2.
/// Returns None with a warning if any required column is missing.
fn parse_maf_header(header_line: &str) -> Option<MafColumnMap> {
    let fields: Vec<&str> = header_line.split('\t').collect();
    let col_map: std::collections::HashMap<&str, usize> = fields.iter()
        .enumerate()
        .map(|(i, &name)| (name, i))
        .collect();

    let chrom = match col_map.get("Chromosome") {
        Some(&idx) => idx,
        None => {
            warn!("MAF header missing required column: Chromosome");
            return None;
        }
    };

    let pos = match col_map.get("Start_Position") {
        Some(&idx) => idx,
        None => {
            warn!("MAF header missing required column: Start_Position");
            return None;
        }
    };

    let ref_allele = match col_map.get("Reference_Allele") {
        Some(&idx) => idx,
        None => {
            warn!("MAF header missing required column: Reference_Allele");
            return None;
        }
    };

    let alt_allele = match col_map.get("Tumor_Seq_Allele2") {
        Some(&idx) => idx,
        None => {
            warn!("MAF header missing required column: Tumor_Seq_Allele2");
            return None;
        }
    };

    let max_idx = *[chrom, pos, ref_allele, alt_allele].iter().max().unwrap();

    debug!("MAF columns resolved: Chromosome={}, Start_Position={}, Reference_Allele={}, Tumor_Seq_Allele2={}",
        chrom, pos, ref_allele, alt_allele);

    Some(MafColumnMap { chrom, pos, ref_allele, alt_allele, max_idx })
}

/// Validate that an allele string contains only valid nucleotide characters.
///
/// Valid characters: A, C, G, T, N (case-insensitive) and `-` (MAF convention
/// for empty alleles in insertions/deletions).
///
/// Rejects non-nucleotide strings like "SNP", "DEL", "Missense_Mutation", etc.
/// which would indicate a MAF column mapping error.
fn is_valid_allele(allele: &str) -> bool {
    !allele.is_empty() && allele.chars().all(|c| "ACGTNacgtn-".contains(c))
}

// ============================================================================
// PHASE 2: Variant Type Classification
// ============================================================================

/// Classify variant type from REF and ALT alleles
fn classify_variant(ref_allele: &str, alt_allele: &str) -> VariantType {
    let ref_len = ref_allele.len();
    let alt_len = alt_allele.len();
    
    if ref_len == 1 && alt_len == 1 {
        VariantType::Snv
    } else if ref_len == alt_len {
        VariantType::Mnv
    } else if ref_len < alt_len && alt_allele.starts_with(ref_allele) {
        // Pure insertion: A -> ATG (ref is prefix of alt)
        VariantType::Insertion
    } else if ref_len > alt_len && ref_allele.starts_with(alt_allele) {
        // Pure deletion: ATG -> A (alt is prefix of ref)
        VariantType::Deletion
    } else {
        // Complex: substitution + indel
        VariantType::Complex
    }
}

// ============================================================================
// PHASE 3: CIGAR-aware Sequence Extraction
// ============================================================================

/// Extract sequence from a read at a variant position
/// Returns: (extracted_sequence, has_n_base, spans_variant)
fn extract_sequence_at_variant(
    record: &bam::Record,
    var: &Variant,
) -> Option<(String, bool)> {
    let read_start = record.pos() as i64;
    let seq = record.seq();
    
    // For SNV/MNV: extract `ref_len` bases starting at variant position
    // For insertions: check if read has insertion at position
    // For deletions: check if read has deletion at position
    
    match var.var_type {
        VariantType::Snv => {
            extract_snv_sequence(record, var.pos, read_start, &seq)
        },
        VariantType::Mnv => {
            extract_mnv_sequence(record, var.pos, var.ref_allele.len(), read_start, &seq)
        },
        VariantType::Insertion => {
            extract_insertion_sequence(record, var, read_start, &seq)
        },
        VariantType::Deletion => {
            extract_deletion_sequence(record, var, read_start, &seq)
        },
        VariantType::Complex => {
            // For complex variants, extract the full REF-length region
            extract_mnv_sequence(record, var.pos, var.ref_allele.len(), read_start, &seq)
        },
    }
}

/// Extract single base for SNV
fn extract_snv_sequence(
    record: &bam::Record,
    var_pos: i64,
    read_start: i64,
    seq: &rust_htslib::bam::record::Seq,
) -> Option<(String, bool)> {
    if var_pos < read_start { return None; }
    
    let target_offset = (var_pos - read_start) as usize;
    let mut ref_offset = 0usize;
    let mut query_offset = 0usize;
    
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len = *len as usize;
                if target_offset >= ref_offset && target_offset < ref_offset + len {
                    let dist = target_offset - ref_offset;
                    let q_idx = query_offset + dist;
                    if q_idx < seq.len() {
                        let base = seq[q_idx] as char;
                        let has_n = base == 'N' || base == 'n';
                        return Some((base.to_string().to_uppercase(), has_n));
                    }
                    return None;
                }
                ref_offset += len;
                query_offset += len;
            },
            Cigar::Ins(len) => { query_offset += *len as usize; },
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                let len = *len as usize;
                // If variant position falls within deletion, read doesn't span
                if target_offset >= ref_offset && target_offset < ref_offset + len {
                    return None;
                }
                ref_offset += len;
            },
            Cigar::SoftClip(len) => { query_offset += *len as usize; },
            Cigar::HardClip(_) | Cigar::Pad(_) => {},
        }
    }
    None
}

/// Extract multiple bases for MNV/Complex
fn extract_mnv_sequence(
    record: &bam::Record,
    var_pos: i64,
    ref_len: usize,
    read_start: i64,
    seq: &rust_htslib::bam::record::Seq,
) -> Option<(String, bool)> {
    if var_pos < read_start { return None; }
    
    let mut extracted = String::with_capacity(ref_len);
    let mut has_n = false;
    
    for offset in 0..ref_len as i64 {
        let pos = var_pos + offset;
        if let Some((base, is_n)) = extract_snv_sequence(record, pos, read_start, seq) {
            extracted.push_str(&base);
            if is_n { has_n = true; }
        } else {
            return None; // Read doesn't span this position
        }
    }
    
    Some((extracted, has_n))
}

/// Check for insertion at variant position
fn extract_insertion_sequence(
    record: &bam::Record,
    var: &Variant,
    read_start: i64,
    seq: &rust_htslib::bam::record::Seq,
) -> Option<(String, bool)> {
    if var.pos < read_start { return None; }
    
    let target_offset = (var.pos - read_start) as usize;
    let expected_ins_len = var.alt_allele.len() - var.ref_allele.len();
    
    let mut ref_offset = 0usize;
    let mut query_offset = 0usize;
    
    // First, find the anchor base
    let mut anchor_found = false;
    let mut anchor_query_pos = 0usize;
    
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len = *len as usize;
                if !anchor_found && target_offset >= ref_offset && target_offset < ref_offset + len {
                    let dist = target_offset - ref_offset;
                    anchor_query_pos = query_offset + dist;
                    anchor_found = true;
                }
                ref_offset += len;
                query_offset += len;
            },
            Cigar::Ins(len) => {
                let len = *len as usize;
                // Check if this insertion is right after our anchor position
                if anchor_found && len == expected_ins_len {
                    // Extract anchor + inserted bases
                    let mut result = String::new();
                    let mut has_n = false;
                    
                    // Anchor base
                    if anchor_query_pos < seq.len() {
                        let b = seq[anchor_query_pos] as char;
                        result.push(b.to_ascii_uppercase());
                        if b == 'N' || b == 'n' { has_n = true; }
                    }
                    
                    // Inserted bases
                    for i in 0..len {
                        let idx = query_offset + i;
                        if idx < seq.len() {
                            let b = seq[idx] as char;
                            result.push(b.to_ascii_uppercase());
                            if b == 'N' || b == 'n' { has_n = true; }
                        }
                    }
                    
                    return Some((result, has_n));
                }
                query_offset += len;
            },
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_offset += *len as usize;
            },
            Cigar::SoftClip(len) => { query_offset += *len as usize; },
            Cigar::HardClip(_) | Cigar::Pad(_) => {},
        }
    }
    
    // No insertion found - extract just the REF base to compare
    if anchor_found && anchor_query_pos < seq.len() {
        let b = seq[anchor_query_pos] as char;
        let has_n = b == 'N' || b == 'n';
        Some((b.to_ascii_uppercase().to_string(), has_n))
    } else {
        None
    }
}

/// Check for deletion at variant position
fn extract_deletion_sequence(
    record: &bam::Record,
    var: &Variant,
    read_start: i64,
    seq: &rust_htslib::bam::record::Seq,
) -> Option<(String, bool)> {
    if var.pos < read_start { return None; }
    
    let target_offset = (var.pos - read_start) as usize;
    let expected_del_len = var.ref_allele.len() - var.alt_allele.len();
    
    let mut ref_offset = 0usize;
    let mut query_offset = 0usize;
    
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len = *len as usize;
                ref_offset += len;
                query_offset += len;
            },
            Cigar::Ins(len) => { query_offset += *len as usize; },
            Cigar::Del(len) => {
                let len = *len as usize;
                // Check if this deletion starts at our variant position
                if ref_offset == target_offset + 1 && len == expected_del_len {
                    // This read has the deletion - it supports ALT
                    // Return the anchor base only (represents ALT)
                    let anchor_idx = query_offset.saturating_sub(1);
                    if anchor_idx < seq.len() {
                        let b = seq[anchor_idx] as char;
                        let has_n = b == 'N' || b == 'n';
                        return Some((format!("DEL{}", len), has_n));
                    }
                }
                ref_offset += len;
            },
            Cigar::RefSkip(len) => { ref_offset += *len as usize; },
            Cigar::SoftClip(len) => { query_offset += *len as usize; },
            Cigar::HardClip(_) | Cigar::Pad(_) => {},
        }
    }
    
    // No deletion found - extract REF-length sequence
    extract_mnv_sequence(record, var.pos, var.ref_allele.len(), read_start, seq)
}

// ============================================================================
// PHASE 4: Fragment Classification
// ============================================================================

/// Classify a fragment based on extracted sequence
fn classify_fragment(extracted: &str, var: &Variant) -> FragmentClass {
    let ref_upper = var.ref_allele.to_uppercase();
    let alt_upper = var.alt_allele.to_uppercase();
    
    match var.var_type {
        VariantType::Snv | VariantType::Mnv | VariantType::Complex => {
            if extracted == ref_upper {
                FragmentClass::Ref
            } else if extracted == alt_upper {
                FragmentClass::Alt
            } else {
                FragmentClass::NonRef
            }
        },
        VariantType::Insertion => {
            if extracted == alt_upper {
                FragmentClass::Alt
            } else if extracted == ref_upper || extracted.len() == 1 {
                // Single base = no insertion = REF
                FragmentClass::Ref
            } else {
                FragmentClass::NonRef
            }
        },
        VariantType::Deletion => {
            let expected_del_len = var.ref_allele.len() - var.alt_allele.len();
            if extracted == format!("DEL{}", expected_del_len) {
                FragmentClass::Alt
            } else if extracted == ref_upper {
                FragmentClass::Ref
            } else {
                FragmentClass::NonRef
            }
        },
    }
}

// ============================================================================
// PHASE 5: Statistics Calculation
// ============================================================================

const MIN_FOR_KS: usize = 2;

// ============================================================================
// LOG-LIKELIHOOD RATIO (LLR) SCORING
// ============================================================================
//
// For duplex/panel mode with very low fragment counts (N<5), traditional
// statistical tests (KS, t-test) are unreliable. We use a probabilistic
// approach that answers: "Are these fragments more likely tumor or healthy?"
//
// Model: Gaussian size distributions for tumor vs healthy cfDNA
// Default (Human): Healthy=167bp, Tumor=145bp
// Canine (Favaro et al.): Healthy=153bp, Tumor=130bp
// ssDNA library prep: Shift both by -10bp
//
// LLR > 0 → fragments are tumor-like
// LLR < 0 → fragments are healthy-like
// ============================================================================

/// LLR Model parameters for cross-species/assay support
/// 
/// Default: Human cfDNA (167bp healthy, 145bp tumor)
/// 
/// # Use Cases
/// - Human cfDNA: LLRModelParams::default() or ::human()
/// - Canine cfDNA: LLRModelParams::canine() - peaks at 153bp (Favaro et al.)
/// - ssDNA library: LLRModelParams::ssdna() - 10bp shorter due to ligation
/// - FFPE samples: Consider wider sigma values
#[derive(Clone, Debug)]
pub struct LLRModelParams {
    pub healthy_mu: f64,
    pub healthy_sigma: f64,
    pub tumor_mu: f64,
    pub tumor_sigma: f64,
}

impl Default for LLRModelParams {
    fn default() -> Self {
        Self::human()
    }
}

impl LLRModelParams {
    /// Human cfDNA defaults (most common use case)
    /// Peak at 167bp for healthy mono-nucleosome, 145bp for tumor-derived
    pub fn human() -> Self {
        Self {
            healthy_mu: 167.0,
            healthy_sigma: 30.0,
            tumor_mu: 145.0,
            tumor_sigma: 35.0,
        }
    }
    
    /// Canine cfDNA model (Favaro et al.)
    /// Dog cfDNA peaks at 153bp (shorter than human due to smaller nucleosomes)
    pub fn canine() -> Self {
        Self {
            healthy_mu: 153.0,
            healthy_sigma: 25.0,
            tumor_mu: 130.0,
            tumor_sigma: 30.0,
        }
    }
    
    /// ssDNA ligation library prep (peaks shifted -10bp)
    /// Single-strand library preps produce shorter fragments
    pub fn ssdna() -> Self {
        Self {
            healthy_mu: 157.0,
            healthy_sigma: 30.0,
            tumor_mu: 135.0,
            tumor_sigma: 35.0,
        }
    }
    
    /// Custom model from explicit parameters
    pub fn custom(healthy_mu: f64, healthy_sigma: f64, tumor_mu: f64, tumor_sigma: f64) -> Self {
        Self { healthy_mu, healthy_sigma, tumor_mu, tumor_sigma }
    }
}

/// Probability density function for healthy cfDNA fragment sizes (parameterized)
fn prob_healthy_param(length: f64, model: &LLRModelParams) -> f64 {
    let z = (length - model.healthy_mu) / model.healthy_sigma;
    let norm_const = model.healthy_sigma * (2.0 * std::f64::consts::PI).sqrt();
    (-0.5 * z * z).exp() / norm_const
}

/// Probability density function for tumor-derived cfDNA fragment sizes (parameterized)
fn prob_tumor_param(length: f64, model: &LLRModelParams) -> f64 {
    let z = (length - model.tumor_mu) / model.tumor_sigma;
    let norm_const = model.tumor_sigma * (2.0 * std::f64::consts::PI).sqrt();
    (-0.5 * z * z).exp() / norm_const
}

/// Calculate Log-Likelihood Ratio for a set of fragment lengths (parameterized)
/// 
/// Uses configurable model parameters for cross-species/assay support.
/// 
/// Returns: sum of log(P_tumor / P_healthy) for each fragment
fn calc_log_likelihood_ratio_param(lengths: &[f64], model: &LLRModelParams) -> f64 {
    if lengths.is_empty() {
        return 0.0;
    }
    
    lengths.iter().map(|&len| {
        let p_h = prob_healthy_param(len, model).max(1e-15);
        let p_t = prob_tumor_param(len, model).max(1e-15);
        p_t.ln() - p_h.ln()  // Positive = tumor-like
    }).sum()
}

/// Calculate Log-Likelihood Ratio for a set of fragment lengths
/// 
/// Returns: sum of log(P_tumor / P_healthy) for each fragment
/// 
/// Interpretation:
/// - LLR > +1.0: Strong tumor signature
/// - LLR > +0.5: Weak tumor signature  
/// - LLR between -0.5 and +0.5: Inconclusive
/// - LLR < -0.5: Healthy-like
/// 
/// Works well even with N=1-5 fragments (unlike KS test)
/// Uses default human model (167bp healthy, 145bp tumor)
fn calc_log_likelihood_ratio(lengths: &[f64]) -> f64 {
    calc_log_likelihood_ratio_param(lengths, &LLRModelParams::human())
}

/// Calculate mean of a vector, returns 0.0 if empty
fn calc_mean(v: &[f64]) -> f64 {
    if v.is_empty() { 0.0 } else { v.iter().sum::<f64>() / v.len() as f64 }
}

/// Two-sample Kolmogorov-Smirnov test
/// Returns (D statistic, p-value)
fn ks_test(a: &[f64], b: &[f64]) -> (f64, f64) {
    if a.len() < MIN_FOR_KS || b.len() < MIN_FOR_KS {
        return (f64::NAN, 1.0);
    }
    
    let mut a_sorted = a.to_vec();
    let mut b_sorted = b.to_vec();
    a_sorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
    b_sorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
    
    let n_a = a_sorted.len();
    let n_b = b_sorted.len();
    let mut i = 0;
    let mut j = 0;
    let mut d_max: f64 = 0.0;
    let mut cdf_a: f64 = 0.0;
    let mut cdf_b: f64 = 0.0;
    
    while i < n_a && j < n_b {
        let v_a = a_sorted[i];
        let v_b = b_sorted[j];
        if v_a < v_b {
            cdf_a += 1.0 / n_a as f64;
            i += 1;
        } else if v_b < v_a {
            cdf_b += 1.0 / n_b as f64;
            j += 1;
        } else {
            cdf_a += 1.0 / n_a as f64;
            cdf_b += 1.0 / n_b as f64;
            i += 1;
            j += 1;
        }
        d_max = d_max.max((cdf_a - cdf_b).abs());
    }
    
    // P-value approximation
    let m = n_a as f64;
    let n = n_b as f64;
    let en = (m * n) / (m + n);
    let lambda = (en.sqrt() + 0.12 + 0.11 / en.sqrt()) * d_max;
    
    let mut p_val = 0.0;
    for k in 1..=100 {
        let term = (-1.0_f64).powi(k - 1) * (-2.0 * (k as f64).powi(2) * lambda.powi(2)).exp();
        p_val += term;
    }
    let p_val = (2.0 * p_val).clamp(0.0, 1.0);
    
    (d_max, p_val)
}

/// Confidence level based on count
fn confidence_level(n: usize) -> &'static str {
    if n >= 5 { "HIGH" }
    else if n >= 1 { "LOW" }
    else { "NONE" }
}

// ============================================================================
// PHASE 6: Main Entry Point
// ============================================================================

/// Calculate Mutant Fragment Size Distribution (mFSD) - Enhanced Version
/// 
/// Supports: SNV, MNV, Insertion, Deletion, Complex variants
/// Classification: REF, ALT, NonREF, N (4-way)
/// 
/// # Arguments
/// * `bam_path` - Path to the input BAM file
/// * `input_file` - Path to VCF/MAF file containing variants
/// * `output_file` - Path to output TSV
/// * `input_format` - "vcf" or "maf" (or "auto")
/// * `map_quality` - Minimum mapping quality
/// * `output_distributions` - If true, write per-variant size distributions
/// * `reference_path` - Optional path to reference FASTA (for GC computation)
/// * `correction_factors_path` - Optional path to pre-computed correction_factors.csv
#[pyfunction]
#[pyo3(signature = (bam_path, input_file, output_file, input_format, map_quality, min_frag_len=65, max_frag_len=400, output_distributions=false, reference_path=None, correction_factors_path=None, require_proper_pair=false, duplex_mode=false, silent=false))]
pub fn calculate_mfsd(
    bam_path: PathBuf,
    input_file: PathBuf,
    output_file: PathBuf,
    input_format: String,
    map_quality: u8,
    min_frag_len: i64,
    max_frag_len: i64,
    output_distributions: bool,
    reference_path: Option<PathBuf>,
    correction_factors_path: Option<PathBuf>,
    require_proper_pair: bool,
    duplex_mode: bool,
    silent: bool,
) -> PyResult<()> {
    use crate::gc_correction::CorrectionFactors;
    use std::sync::Arc;
    
    // Load correction factors if provided
    let factors: Option<Arc<CorrectionFactors>> = if let Some(ref factors_path) = correction_factors_path {
        match CorrectionFactors::load_csv(factors_path) {
            Ok(f) => {
                info!("Loaded GC correction factors from {:?}", factors_path);
                Some(Arc::new(f))
            },
            Err(e) => {
                warn!("Failed to load correction factors: {}. Proceeding without GC correction.", e);
                None
            }
        }
    } else {
        None
    };
    
    // Load reference FASTA if provided (for GC content computation)
    let _fasta: Option<faidx::Reader> = if let Some(ref ref_path) = reference_path {
        match faidx::Reader::from_path(ref_path) {
            Ok(f) => {
                info!("Loaded reference FASTA: {:?}", ref_path);
                Some(f)
            },
            Err(e) => {
                warn!("Failed to load reference FASTA: {}. GC correction disabled.", e);
                None
            }
        }
    } else {
        None
    };
    // 1. Parse Variants
    let mut variants = Vec::new();
    let file = File::open(&input_file)?;
    let reader = BufReader::new(file);
    
    let is_vcf = input_format == "vcf" || 
        (input_format == "auto" && input_file.extension().map_or(false, |e| e == "vcf" || e == "gz"));
    
    info!("Parsing variants from {:?}...", input_file);
    
    // MAF header-based column resolution (parsed from first non-comment data line)
    let mut maf_cols: Option<MafColumnMap> = None;
    let mut skipped_invalid_alleles: usize = 0;
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') { continue; }
        
        let fields: Vec<&str> = line.split('\t').collect();
        
        let (chrom, pos, ref_allele, alt_allele) = if is_vcf {
            if fields.len() < 5 { continue; }
            let pos: i64 = fields[1].parse().unwrap_or(0) - 1;
            (fields[0].to_string(), pos, fields[3].to_string(), fields[4].to_string())
        } else {
            // MAF Parsing: header-based column lookup
            
            // Detect and parse MAF header line
            if maf_cols.is_none() {
                if fields[0] == "Hugo_Symbol" || fields[0].starts_with("Hugo_Symbol") {
                    maf_cols = parse_maf_header(&line);
                    if maf_cols.is_none() {
                        return Err(pyo3::exceptions::PyValueError::new_err(
                            "Failed to parse MAF header: missing required columns (Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2)"
                        ));
                    }
                    info!("MAF header parsed: {} columns detected", fields.len());
                    continue; // Header line, not a data row
                } else {
                    // No header found — fall back to error
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        "MAF file missing header line (expected first non-comment line to start with Hugo_Symbol)"
                    ));
                }
            }
            
            let cols = maf_cols.as_ref().unwrap();
            
            // Bounds check using resolved max index
            if fields.len() <= cols.max_idx { continue; }
            
            let pos: i64 = fields[cols.pos].parse().unwrap_or(0) - 1;
            let ref_allele = fields[cols.ref_allele].to_string();
            let alt_allele = fields[cols.alt_allele].to_string();
            
            // Allele validation — catch column mapping errors early
            if !is_valid_allele(&ref_allele) || !is_valid_allele(&alt_allele) {
                if skipped_invalid_alleles == 0 {
                    // Log first occurrence with details for debugging
                    warn!("Invalid allele at {}:{}: REF='{}' ALT='{}' — skipping (possible MAF column mismatch)",
                        fields[cols.chrom], fields[cols.pos], ref_allele, alt_allele);
                }
                skipped_invalid_alleles += 1;
                continue;
            }
            
            // Handle MAF dash convention: '-' represents empty allele for indels
            // Convert to VCF-style representation for downstream processing
            let ref_allele = if ref_allele == "-" { String::new() } else { ref_allele };
            let alt_allele = if alt_allele == "-" { String::new() } else { alt_allele };
            
            (fields[cols.chrom].to_string(), pos, ref_allele, alt_allele)
        };
        
        let var_type = classify_variant(&ref_allele, &alt_allele);
        
        variants.push(Variant {
            chrom,
            pos,
            ref_allele,
            alt_allele,
            var_type,
        });
    }
    
    // Log invalid allele summary
    if skipped_invalid_alleles > 0 {
        warn!("Skipped {} variants with invalid alleles (non-nucleotide characters in REF/ALT)", 
            skipped_invalid_alleles);
    }
    
    let total_vars = variants.len();
    info!("Found {} variants. Processing in parallel...", total_vars);
    
    // Debug: variant type breakdown
    let snv_count = variants.iter().filter(|v| v.var_type == VariantType::Snv).count();
    let mnv_count = variants.iter().filter(|v| v.var_type == VariantType::Mnv).count();
    let ins_count = variants.iter().filter(|v| v.var_type == VariantType::Insertion).count();
    let del_count = variants.iter().filter(|v| v.var_type == VariantType::Deletion).count();
    let complex_count = variants.iter().filter(|v| v.var_type == VariantType::Complex).count();
    debug!("Variant types: {} SNV, {} MNV, {} INS, {} DEL, {} Complex", 
        snv_count, mnv_count, ins_count, del_count, complex_count);
    
    // Log filter settings
    info!("Fragment filters: length {}-{}bp, proper_pair={}, duplex_mode={}", min_frag_len, max_frag_len, require_proper_pair, duplex_mode);

    // Progress Bar (hidden when called from wrapper)
    let pb = if silent {
        ProgressBar::hidden()
    } else {
        let pb = ProgressBar::new(total_vars as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"));
        pb.enable_steady_tick(Duration::from_millis(100));
        pb
    };

    // 2. Process Variants (Parallel)
    let results: Vec<(Variant, VariantResult)> = variants.par_iter()
        .map(|var| {
            let mut result = VariantResult::new();
            
            // Thread-local BAM reader
            let mut bam = match bam::IndexedReader::from_path(&bam_path) {
                Ok(b) => b,
                Err(_) => return (var.clone(), result),
            };
            
            // Thread-local FASTA reader (for GC correction)
            let thread_fasta: Option<faidx::Reader> = if factors.is_some() {
                reference_path.as_ref().and_then(|p| faidx::Reader::from_path(p).ok())
            } else {
                None
            };

            // Get chromosome tid
            let tid = match bam.header().tid(var.chrom.as_bytes()) {
                Some(t) => t,
                None => {
                    let alt_name = if var.chrom.starts_with("chr") {
                        var.chrom.trim_start_matches("chr").to_string()
                    } else {
                        format!("chr{}", var.chrom)
                    };
                    match bam.header().tid(alt_name.as_bytes()) {
                        Some(t) => t,
                        None => {
                            warn!("Chromosome {} not found in BAM for variant at pos {}", var.chrom, var.pos + 1);
                            return (var.clone(), result);
                        }
                    }
                }
            };
            
            // Fetch region - extend to cover variant span
            let var_end = var.pos + var.ref_allele.len().max(var.alt_allele.len()) as i64;
            if bam.fetch((tid, var.pos as u64, var_end as u64 + 1)).is_err() {
                return (var.clone(), result);
            }
            
            // Process reads - ONLY R1 to count fragments, not reads
            for record_res in bam.records() {
                let record = match record_res {
                    Ok(r) => r,
                    Err(_) => continue,
                };
                
                // Filters
                if record.mapq() < map_quality { continue; }
                if record.is_duplicate() { continue; }
                if record.is_unmapped() { continue; }
                if record.is_secondary() { continue; }
                if record.is_supplementary() { continue; }
                
                // Only process R1 for fragment counting (avoid double-counting)
                if record.is_paired() && !record.is_first_in_template() { continue; }
                
                // Proper pair filter (optional - disabled by default for duplex BAMs)
                if require_proper_pair && record.is_paired() && !record.is_proper_pair() {
                    result.skipped_proper_pair += 1;
                    continue;
                }
                
                // Get fragment length first to apply size filter early
                let tlen = record.insert_size().abs();
                let frag_len = if tlen == 0 { record.seq().len() as i64 } else { tlen };
                
                // Fragment size filter - skip discordant/chimeric reads with extreme TLEN
                if frag_len < min_frag_len || frag_len > max_frag_len {
                    result.skipped_size += 1;
                    // Track extreme outliers (discordant reads) for summary logging
                    if frag_len > 10000 {
                        result.skipped_extreme += 1;
                    }
                    continue;
                }
                
                let frag_len_f64 = frag_len as f64;
                
                // Extract sequence at variant
                let (extracted, has_n) = match extract_sequence_at_variant(&record, var) {
                    Some((s, n)) => (s, n),
                    None => continue, // Read doesn't span variant
                };
                
                // Compute GC correction weight
                let gc_weight = if let Some(ref factors_arc) = factors {
                    // Get fragment coordinates for GC lookup
                    let frag_start = record.pos() as u64;
                    let frag_end = if record.insert_size() > 0 {
                        frag_start + record.insert_size() as u64
                    } else {
                        frag_start + record.seq().len() as u64
                    };
                    
                    // Get GC from thread-local FASTA reader
                    let gc_pct = if let Some(ref fasta) = thread_fasta {
                        fetch_reference_gc(fasta, &var.chrom, frag_start, frag_end)
                            .map(|gc| (gc * 100.0).round() as u8)
                            .unwrap_or(50) // Default 50% GC
                    } else {
                        50 // Default if no reference
                    };
                    
                    factors_arc.get_factor(frag_len as u64, gc_pct)
                } else {
                    1.0 // No GC correction
                };
                
                // Compute duplex consensus weight (for fgbio/Marianas duplex BAMs)
                // Only applied when --duplex flag is set
                let duplex_weight = if duplex_mode {
                    get_duplex_weight(&record)
                } else {
                    1.0  // No duplex weighting for non-duplex data
                };
                
                // Track duplex tags found (weight != 1.0 means a tag was found)
                if duplex_mode && duplex_weight != 1.0 {
                    result.duplex_tags_found += 1;
                }
                
                // Combined weight: GC correction * duplex confidence
                let weight = gc_weight * duplex_weight;
                
                // Classify and add with combined weight
                let class = if has_n {
                    FragmentClass::N
                } else {
                    classify_fragment(&extracted, var)
                };
                
                result.add_fragment(class, frag_len_f64, weight);
            }
            
            pb.inc(1);
            (var.clone(), result)
        })
        .collect();
        
    pb.finish_with_message("Done!");

    // Calculate summary statistics
    let mut total_ref = 0usize;
    let mut total_alt = 0usize;
    let mut total_nonref = 0usize;
    let mut total_n = 0usize;
    let mut variants_with_alt = 0usize;
    let mut variants_no_coverage = 0usize;
    let mut variants_suspicious = 0usize; // 0 REF + 0 ALT + high NonREF = likely parsing error
    let mut total_skipped_proper_pair = 0usize;
    let mut total_skipped_size = 0usize;
    let mut total_skipped_extreme = 0usize;
    
    for (var, res) in &results {
        total_ref += res.ref_lengths.len();
        total_alt += res.alt_lengths.len();
        total_nonref += res.nonref_lengths.len();
        total_n += res.n_lengths.len();
        total_skipped_proper_pair += res.skipped_proper_pair;
        total_skipped_size += res.skipped_size;
        total_skipped_extreme += res.skipped_extreme;
        if !res.alt_lengths.is_empty() {
            variants_with_alt += 1;
        }
        if res.total_count() == 0 {
            variants_no_coverage += 1;
        }
        // Sanity check: 0 REF + 0 ALT but many NonREF strongly suggests allele parsing error
        if res.ref_lengths.is_empty() && res.alt_lengths.is_empty() && res.nonref_lengths.len() > 10 {
            if variants_suspicious == 0 {
                warn!("Suspicious: variant {}:{} has {} NonREF but 0 REF/0 ALT — possible allele parsing error (REF='{}', ALT='{}')",
                    var.chrom, var.pos + 1, res.nonref_lengths.len(), var.ref_allele, var.alt_allele);
            }
            variants_suspicious += 1;
        }
    }
    
    info!("Summary: {} REF, {} ALT, {} NonREF, {} N fragments", total_ref, total_alt, total_nonref, total_n);
    info!("Variants with ALT support: {}/{} ({:.1}%)", 
        variants_with_alt, results.len(), 
        if results.len() > 0 { variants_with_alt as f64 / results.len() as f64 * 100.0 } else { 0.0 });
    
    // Log filter statistics
    if total_skipped_proper_pair > 0 || total_skipped_size > 0 {
        info!("Filtered: {} improper-pair, {} out-of-size-range ({}-{}bp)", 
            total_skipped_proper_pair, total_skipped_size, min_frag_len, max_frag_len);
    }
    
    // Log extreme outliers if any (discordant reads with TLEN > 10000bp)
    if total_skipped_extreme > 0 {
        debug!("Skipped {} discordant reads with TLEN > 10000bp", total_skipped_extreme);
    }
    
    if variants_no_coverage > 0 {
        warn!("{} variants had no fragment coverage", variants_no_coverage);
    }
    
    // Sanity: warn if many variants have 0 REF + 0 ALT (suggests systemic allele parsing issue)
    if variants_suspicious > 0 {
        warn!("{}/{} variants have 0 REF + 0 ALT but high NonREF — check MAF allele columns",
            variants_suspicious, results.len());
    }
    
    // Duplex tag warning: alert user if --duplex was set but no tags were found
    if duplex_mode {
        let total_fragments = total_ref + total_alt + total_nonref + total_n;
        let total_duplex_tags: usize = results.iter().map(|(_, res)| res.duplex_tags_found).sum();
        if total_duplex_tags == 0 && total_fragments > 0 {
            warn!("--duplex mode enabled but NO cD/Marianas tags found in {} fragments. \
                   Ensure BAM was processed by fgbio/Marianas consensus caller.", total_fragments);
        } else {
            info!("Duplex tags found in {}/{} ({:.1}%) fragments", 
                total_duplex_tags, total_fragments,
                if total_fragments > 0 { total_duplex_tags as f64 / total_fragments as f64 * 100.0 } else { 0.0 });
        }
    }

    // 3. Write Main Output
    info!("Writing output to {:?}...", output_file);
    let mut out_file = File::create(&output_file)?;
    
    // Header (42 columns - added 5 GC-weighted columns)
    writeln!(out_file, "{}", [
        // Variant info (5)
        "Chrom", "Pos", "Ref", "Alt", "VarType",
        // Counts (5)
        "REF_Count", "ALT_Count", "NonREF_Count", "N_Count", "Total_Count",
        // GC-Weighted Counts (5)
        "REF_Weighted", "ALT_Weighted", "NonREF_Weighted", "N_Weighted", "VAF_GC_Corrected",
        // Log-Likelihood Ratio (2) - for duplex/panel with low N
        "ALT_LLR", "REF_LLR",
        // Mean sizes (4)
        "REF_MeanSize", "ALT_MeanSize", "NonREF_MeanSize", "N_MeanSize",
        // Primary: ALT vs REF (3)
        "Delta_ALT_REF", "KS_ALT_REF", "KS_Pval_ALT_REF",
        // Secondary: ALT vs NonREF (3)
        "Delta_ALT_NonREF", "KS_ALT_NonREF", "KS_Pval_ALT_NonREF",
        // REF vs NonREF (3)
        "Delta_REF_NonREF", "KS_REF_NonREF", "KS_Pval_REF_NonREF",
        // ALT vs N (3)
        "Delta_ALT_N", "KS_ALT_N", "KS_Pval_ALT_N",
        // Tertiary: REF vs N (3)
        "Delta_REF_N", "KS_REF_N", "KS_Pval_REF_N",
        // NonREF vs N (3)
        "Delta_NonREF_N", "KS_NonREF_N", "KS_Pval_NonREF_N",
        // Derived (5)
        "VAF_Proxy", "Error_Rate", "N_Rate", "Size_Ratio", "Quality_Score",
        // Quality flags (2)
        "ALT_Confidence", "KS_Valid",
    ].join("\t"))?;
    
    // Optional: distributions file
    let mut dist_file = if output_distributions {
        let dist_path = output_file.with_extension("distributions.tsv");
        info!("Writing distributions to {:?}...", dist_path);
        let f = File::create(&dist_path)?;
        let mut f = std::io::BufWriter::new(f);
        writeln!(f, "Chrom\tPos\tRef\tAlt\tCategory\tSize\tCount")?;
        Some(f)
    } else {
        None
    };
    
    for (var, res) in &results {
        // Counts
        let n_ref = res.ref_lengths.len();
        let n_alt = res.alt_lengths.len();
        let n_nonref = res.nonref_lengths.len();
        let n_n = res.n_lengths.len();
        let n_total = res.total_count();
        
        // Means
        let mean_ref = calc_mean(&res.ref_lengths);
        let mean_alt = calc_mean(&res.alt_lengths);
        let mean_nonref = calc_mean(&res.nonref_lengths);
        let mean_n = calc_mean(&res.n_lengths);
        
        // Pairwise comparisons
        let (ks_alt_ref, pval_alt_ref) = ks_test(&res.alt_lengths, &res.ref_lengths);
        let delta_alt_ref = if n_alt > 0 && n_ref > 0 { mean_alt - mean_ref } else { f64::NAN };
        
        let (ks_alt_nonref, pval_alt_nonref) = ks_test(&res.alt_lengths, &res.nonref_lengths);
        let delta_alt_nonref = if n_alt > 0 && n_nonref > 0 { mean_alt - mean_nonref } else { f64::NAN };
        
        let (ks_ref_nonref, pval_ref_nonref) = ks_test(&res.ref_lengths, &res.nonref_lengths);
        let delta_ref_nonref = if n_ref > 0 && n_nonref > 0 { mean_ref - mean_nonref } else { f64::NAN };
        
        let (ks_alt_n, pval_alt_n) = ks_test(&res.alt_lengths, &res.n_lengths);
        let delta_alt_n = if n_alt > 0 && n_n > 0 { mean_alt - mean_n } else { f64::NAN };
        
        let (ks_ref_n, pval_ref_n) = ks_test(&res.ref_lengths, &res.n_lengths);
        let delta_ref_n = if n_ref > 0 && n_n > 0 { mean_ref - mean_n } else { f64::NAN };
        
        let (ks_nonref_n, pval_nonref_n) = ks_test(&res.nonref_lengths, &res.n_lengths);
        let delta_nonref_n = if n_nonref > 0 && n_n > 0 { mean_nonref - mean_n } else { f64::NAN };
        
        // Derived metrics
        let vaf_proxy = if n_alt + n_ref > 0 { n_alt as f64 / (n_alt + n_ref) as f64 } else { 0.0 };
        let error_rate = if n_total > 0 { n_nonref as f64 / n_total as f64 } else { 0.0 };
        let n_rate = if n_total > 0 { n_n as f64 / n_total as f64 } else { 0.0 };
        let size_ratio = if mean_ref > 0.0 { mean_alt / mean_ref } else { f64::NAN };
        let quality_score = 1.0 - n_rate - error_rate;
        
        // Quality flags
        let alt_confidence = confidence_level(n_alt);
        let ks_valid = n_alt >= MIN_FOR_KS && n_ref >= MIN_FOR_KS;
        
        let w_ref = res.ref_weighted;
        let w_alt = res.alt_weighted;
        let w_nonref = res.nonref_weighted;
        let w_n = res.n_weighted;
        let vaf_gc_corrected = if w_alt + w_ref > 0.0 { w_alt / (w_alt + w_ref) } else { 0.0 };
        
        // LLR scores - for low-N duplex/panel mode
        // Positive = tumor-like, Negative = healthy-like
        let alt_llr = calc_log_likelihood_ratio(&res.alt_lengths);
        let ref_llr = calc_log_likelihood_ratio(&res.ref_lengths);
        
        // Format NaN as "NA"
        let fmt = |v: f64| if v.is_nan() { "NA".to_string() } else { format!("{:.4}", v) };
        
        writeln!(out_file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            // Variant info (5)
            var.chrom, var.pos + 1, var.ref_allele, var.alt_allele, var.var_type.as_str(),
            // Raw Counts (5)
            n_ref, n_alt, n_nonref, n_n, n_total,
            // GC-Weighted Counts (5)
            fmt(w_ref), fmt(w_alt), fmt(w_nonref), fmt(w_n), fmt(vaf_gc_corrected),
            // LLR scores (2)
            fmt(alt_llr), fmt(ref_llr),
            // Mean sizes (4)
            fmt(mean_ref), fmt(mean_alt), fmt(mean_nonref), fmt(mean_n),
            // Primary: ALT vs REF (3)
            fmt(delta_alt_ref), fmt(ks_alt_ref), fmt(pval_alt_ref),
            // ALT vs NonREF (3)
            fmt(delta_alt_nonref), fmt(ks_alt_nonref), fmt(pval_alt_nonref),
            // REF vs NonREF (3)
            fmt(delta_ref_nonref), fmt(ks_ref_nonref), fmt(pval_ref_nonref),
            // ALT vs N (3)
            fmt(delta_alt_n), fmt(ks_alt_n), fmt(pval_alt_n),
            // REF vs N (3)
            fmt(delta_ref_n), fmt(ks_ref_n), fmt(pval_ref_n),
            // NonREF vs N (3)
            fmt(delta_nonref_n), fmt(ks_nonref_n), fmt(pval_nonref_n),
            // Derived (5)
            fmt(vaf_proxy), fmt(error_rate), fmt(n_rate), fmt(size_ratio), fmt(quality_score),
            // Quality flags (2)
            alt_confidence, ks_valid,
        )?;
        
        // Write distributions if requested
        if let Some(ref mut df) = dist_file {
            // Group sizes by count for compact output
            let write_dist = |df: &mut std::io::BufWriter<File>, category: &str, lengths: &[f64]| -> std::io::Result<()> {
                use std::collections::HashMap;
                let mut counts: HashMap<i64, u64> = HashMap::new();
                for &len in lengths {
                    *counts.entry(len as i64).or_default() += 1;
                }
                for (size, count) in counts {
                    writeln!(df, "{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                        var.chrom, var.pos + 1, var.ref_allele, var.alt_allele, category, size, count)?;
                }
                Ok(())
            };
            
            write_dist(df, "REF", &res.ref_lengths)?;
            write_dist(df, "ALT", &res.alt_lengths)?;
            write_dist(df, "NonREF", &res.nonref_lengths)?;
            write_dist(df, "N", &res.n_lengths)?;
        }
    }

    info!("Done! Processed {} variants.", results.len());
    Ok(())
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ======================================================================
    // MAF Header Parsing Tests
    // ======================================================================

    #[test]
    fn test_parse_maf_header_standard() {
        // Standard GDC MAF without Consequence column
        let header = "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2";
        let cols = parse_maf_header(header).expect("Should parse standard MAF header");
        assert_eq!(cols.chrom, 4, "Chromosome should be at index 4");
        assert_eq!(cols.pos, 5, "Start_Position should be at index 5");
        assert_eq!(cols.ref_allele, 10, "Reference_Allele should be at index 10");
        assert_eq!(cols.alt_allele, 12, "Tumor_Seq_Allele2 should be at index 12");
    }

    #[test]
    fn test_parse_maf_header_cbio() {
        // cBioPortal MAF with extra Consequence column at index 8
        // This is the exact format that caused the original bug
        let header = "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tConsequence\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2";
        let cols = parse_maf_header(header).expect("Should parse cBioPortal MAF header");
        assert_eq!(cols.chrom, 4, "Chromosome should be at index 4");
        assert_eq!(cols.pos, 5, "Start_Position should be at index 5");
        assert_eq!(cols.ref_allele, 11, "Reference_Allele should be at index 11 (shifted by Consequence)");
        assert_eq!(cols.alt_allele, 13, "Tumor_Seq_Allele2 should be at index 13 (shifted by Consequence)");
    }

    #[test]
    fn test_parse_maf_header_minimal() {
        // Minimal MAF with only required columns (in non-standard order)
        let header = "Reference_Allele\tChromosome\tStart_Position\tTumor_Seq_Allele2";
        let cols = parse_maf_header(header).expect("Should parse minimal MAF header");
        assert_eq!(cols.chrom, 1);
        assert_eq!(cols.pos, 2);
        assert_eq!(cols.ref_allele, 0);
        assert_eq!(cols.alt_allele, 3);
    }

    #[test]
    fn test_parse_maf_header_missing_col() {
        // Missing Tumor_Seq_Allele2
        let header = "Hugo_Symbol\tChromosome\tStart_Position\tReference_Allele";
        assert!(parse_maf_header(header).is_none(), "Should return None for missing required column");
    }

    // ======================================================================
    // Allele Validation Tests
    // ======================================================================

    #[test]
    fn test_validate_allele_valid() {
        assert!(is_valid_allele("A"), "Single base A");
        assert!(is_valid_allele("ATG"), "Multi-base ATG");
        assert!(is_valid_allele("a"), "Lowercase base");
        assert!(is_valid_allele("N"), "Ambiguous base N");
        assert!(is_valid_allele("ACGTN"), "All valid bases");
    }

    #[test]
    fn test_validate_allele_invalid() {
        assert!(!is_valid_allele("SNP"), "Variant_Type value SNP");
        assert!(!is_valid_allele("DEL"), "Variant_Type value DEL");
        assert!(!is_valid_allele("123"), "Numeric string");
        assert!(!is_valid_allele(""), "Empty string");
        assert!(!is_valid_allele("Missense_Mutation"), "Classification value");
        assert!(!is_valid_allele("A T"), "Contains space");
    }

    #[test]
    fn test_validate_allele_maf_dash() {
        // MAF uses '-' for empty alleles in indels
        assert!(is_valid_allele("-"), "MAF dash for insertion REF or deletion ALT");
        assert!(is_valid_allele("A-T"), "Should not appear but is technically valid chars");
    }

    // ======================================================================
    // Variant Classification Tests
    // ======================================================================

    #[test]
    fn test_classify_variant_snv() {
        assert_eq!(classify_variant("A", "T"), VariantType::Snv);
        assert_eq!(classify_variant("C", "G"), VariantType::Snv);
    }

    #[test]
    fn test_classify_variant_mnv() {
        assert_eq!(classify_variant("AT", "GC"), VariantType::Mnv);
        assert_eq!(classify_variant("ATG", "CCT"), VariantType::Mnv);
    }

    #[test]
    fn test_classify_variant_insertion() {
        assert_eq!(classify_variant("A", "ATG"), VariantType::Insertion);
        assert_eq!(classify_variant("T", "TCCG"), VariantType::Insertion);
    }

    #[test]
    fn test_classify_variant_deletion() {
        assert_eq!(classify_variant("ATG", "A"), VariantType::Deletion);
        assert_eq!(classify_variant("TCCG", "T"), VariantType::Deletion);
    }

    #[test]
    fn test_classify_variant_complex() {
        assert_eq!(classify_variant("ATG", "CT"), VariantType::Complex);
        assert_eq!(classify_variant("AT", "GCT"), VariantType::Complex);
    }

    // ======================================================================
    // Fragment Classification Tests
    // ======================================================================

    #[test]
    fn test_classify_fragment_snv_ref() {
        let var = Variant {
            chrom: "chr1".to_string(), pos: 1000,
            ref_allele: "A".to_string(), alt_allele: "T".to_string(),
            var_type: VariantType::Snv,
        };
        assert_eq!(classify_fragment("A", &var), FragmentClass::Ref);
    }

    #[test]
    fn test_classify_fragment_snv_alt() {
        let var = Variant {
            chrom: "chr1".to_string(), pos: 1000,
            ref_allele: "A".to_string(), alt_allele: "T".to_string(),
            var_type: VariantType::Snv,
        };
        assert_eq!(classify_fragment("T", &var), FragmentClass::Alt);
    }

    #[test]
    fn test_classify_fragment_snv_nonref() {
        let var = Variant {
            chrom: "chr1".to_string(), pos: 1000,
            ref_allele: "A".to_string(), alt_allele: "T".to_string(),
            var_type: VariantType::Snv,
        };
        assert_eq!(classify_fragment("G", &var), FragmentClass::NonRef);
    }

    #[test]
    fn test_classify_fragment_deletion_alt() {
        let var = Variant {
            chrom: "chr1".to_string(), pos: 1000,
            ref_allele: "ATG".to_string(), alt_allele: "A".to_string(),
            var_type: VariantType::Deletion,
        };
        assert_eq!(classify_fragment("DEL2", &var), FragmentClass::Alt);
    }
}
