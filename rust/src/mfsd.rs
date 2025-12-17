use pyo3::prelude::*;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::Cigar;
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
#[derive(Debug, Clone, Default)]
struct VariantResult {
    ref_lengths: Vec<f64>,
    alt_lengths: Vec<f64>,
    nonref_lengths: Vec<f64>,
    n_lengths: Vec<f64>,
}

impl VariantResult {
    fn new() -> Self {
        Self::default()
    }
    
    fn total_count(&self) -> usize {
        self.ref_lengths.len() + self.alt_lengths.len() + 
        self.nonref_lengths.len() + self.n_lengths.len()
    }
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
#[pyfunction]
#[pyo3(signature = (bam_path, input_file, output_file, input_format, map_quality, output_distributions=false))]
pub fn calculate_mfsd(
    bam_path: PathBuf,
    input_file: PathBuf,
    output_file: PathBuf,
    input_format: String,
    map_quality: u8,
    output_distributions: bool,
) -> PyResult<()> {
    // 1. Parse Variants
    let mut variants = Vec::new();
    let file = File::open(&input_file)?;
    let reader = BufReader::new(file);
    
    let is_vcf = input_format == "vcf" || 
        (input_format == "auto" && input_file.extension().map_or(false, |e| e == "vcf" || e == "gz"));
    
    info!("Parsing variants from {:?}...", input_file);
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') { continue; }
        
        let fields: Vec<&str> = line.split('\t').collect();
        
        let (chrom, pos, ref_allele, alt_allele) = if is_vcf {
            if fields.len() < 5 { continue; }
            let pos: i64 = fields[1].parse().unwrap_or(0) - 1;
            (fields[0].to_string(), pos, fields[3].to_string(), fields[4].to_string())
        } else {
            // MAF Parsing
            if fields.len() < 13 { continue; }
            if fields[0] == "Hugo_Symbol" || fields[0].starts_with("Hugo_Symbol") { continue; }
            let pos: i64 = fields[5].parse().unwrap_or(0) - 1;
            (fields[4].to_string(), pos, fields[10].to_string(), fields[12].to_string())
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

    // Progress Bar
    let pb = ProgressBar::new(total_vars as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"));
    pb.enable_steady_tick(Duration::from_millis(100));

    // 2. Process Variants (Parallel)
    let results: Vec<(Variant, VariantResult)> = variants.par_iter()
        .map(|var| {
            let mut result = VariantResult::new();
            
            // Thread-local BAM reader
            let mut bam = match bam::IndexedReader::from_path(&bam_path) {
                Ok(b) => b,
                Err(_) => return (var.clone(), result),
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
                
                // Extract sequence at variant
                let (extracted, has_n) = match extract_sequence_at_variant(&record, var) {
                    Some((s, n)) => (s, n),
                    None => continue, // Read doesn't span variant
                };
                
                // Get fragment length
                let tlen = record.insert_size().abs();
                let frag_len = if tlen == 0 { record.seq().len() as i64 } else { tlen };
                let frag_len = frag_len as f64;
                
                // Classify
                if has_n {
                    result.n_lengths.push(frag_len);
                } else {
                    match classify_fragment(&extracted, var) {
                        FragmentClass::Ref => result.ref_lengths.push(frag_len),
                        FragmentClass::Alt => result.alt_lengths.push(frag_len),
                        FragmentClass::NonRef => result.nonref_lengths.push(frag_len),
                        FragmentClass::N => result.n_lengths.push(frag_len),
                    }
                }
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
    
    for (_, res) in &results {
        total_ref += res.ref_lengths.len();
        total_alt += res.alt_lengths.len();
        total_nonref += res.nonref_lengths.len();
        total_n += res.n_lengths.len();
        if !res.alt_lengths.is_empty() {
            variants_with_alt += 1;
        }
        if res.total_count() == 0 {
            variants_no_coverage += 1;
        }
    }
    
    info!("Summary: {} REF, {} ALT, {} NonREF, {} N fragments", total_ref, total_alt, total_nonref, total_n);
    info!("Variants with ALT support: {}/{} ({:.1}%)", 
        variants_with_alt, results.len(), 
        if results.len() > 0 { variants_with_alt as f64 / results.len() as f64 * 100.0 } else { 0.0 });
    if variants_no_coverage > 0 {
        warn!("{} variants had no fragment coverage", variants_no_coverage);
    }

    // 3. Write Main Output
    info!("Writing output to {:?}...", output_file);
    let mut out_file = File::create(&output_file)?;
    
    // Header (37 columns)
    writeln!(out_file, "{}", [
        // Variant info (5)
        "Chrom", "Pos", "Ref", "Alt", "VarType",
        // Counts (5)
        "REF_Count", "ALT_Count", "NonREF_Count", "N_Count", "Total_Count",
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
        
        // Format NaN as "NA"
        let fmt = |v: f64| if v.is_nan() { "NA".to_string() } else { format!("{:.4}", v) };
        
        writeln!(out_file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            // Variant info (5)
            var.chrom, var.pos + 1, var.ref_allele, var.alt_allele, var.var_type.as_str(),
            // Counts (5)
            n_ref, n_alt, n_nonref, n_n, n_total,
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
