use pyo3::prelude::*;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::Cigar;
use rayon::prelude::*;
use indicatif::{ProgressBar, ProgressStyle};
use std::time::Duration;
use log::info;

/// Calculate Mutant Fragment Size Distribution (mFSD)
/// 
/// # Arguments
/// * `bam_path` - Path to the input BAM file
/// * `input_file` - Path to VCF/MAF file containing variants
/// * `output_file` - Path to output TSV
/// * `input_format` - "vcf" or "maf" (or "auto")
/// * `map_quality` - Minimum mapping quality
#[pyfunction]
pub fn calculate_mfsd(
    bam_path: PathBuf,
    input_file: PathBuf,
    output_file: PathBuf,
    input_format: String,
    map_quality: u8,
) -> PyResult<()> {
    // 1. Parse Variants
    struct Variant {
        chrom: String,
        pos: i64, // 0-based
        ref_base: String,
        alt_base: String,
    }
    
    // SAFETY: We need to share BAM path across threads. PathBuf is Clone/Sync.
    // Variants need to be Sync (they are).

    let mut variants = Vec::new();
    let file = File::open(&input_file)?;
    let reader = BufReader::new(file);
    
    let is_vcf = input_format == "vcf" || (input_format == "auto" && input_file.extension().map_or(false, |e| e == "vcf" || e == "gz"));
    
    info!("Parsing variants from {:?}...", input_file);
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') { continue; }
        
        let fields: Vec<&str> = line.split('\t').collect();
        
        if is_vcf {
             if fields.len() < 5 { continue; }
             let pos: i64 = fields[1].parse().unwrap_or(0) - 1; // VCF is 1-based, convert to 0-based
             variants.push(Variant {
                 chrom: fields[0].to_string(),
                 pos,
                 ref_base: fields[3].to_string(),
                 alt_base: fields[4].to_string(),
             });
        } else {
            // MAF Parsing
            if fields.len() < 13 { continue; }
            if fields[0] == "Hugo_Symbol" || fields[0].starts_with("Hugo_Symbol") { continue; }
            
            let pos: i64 = fields[5].parse().unwrap_or(0) - 1; // MAF is 1-based, convert to 0-based
            variants.push(Variant {
                chrom: fields[4].to_string(),
                pos,
                ref_base: fields[10].to_string(),
                alt_base: fields[12].to_string(),
            });
        }
    }
    
    let total_vars = variants.len();
    info!("Found {} variants. Processing in parallel...", total_vars);
    
    // Progress Bar
    let pb = ProgressBar::new(total_vars as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"));
    pb.enable_steady_tick(Duration::from_millis(100));

    // 2. Process Variants (Parallel)
    // We collect results into a Vec<String> to write sequentially at the end
    // Or we can use a parallel iterator and write to a Mutex<File>, but collecting strings is safer/faster for small outputs.
    // If output is massive, we should stream, but mFSD output is one line per variant (small).
    
    let results: Vec<String> = variants.par_iter()
        .map(|var| {
            // Thread-local BAM reader
            // In efficient implementation we might use a pool, but opening a fresh reader per variant 
            // is acceptable if variant count isn't massive (thousands is fine).
            // Opening BAM is relatively cheap (reading header), seeking is the cost.
            // If optimization is needed, use `thread_local!` or a specialized pool.
            // For now, simple approach: new reader per thread (or per task if many tasks).
            // Actually, best practice with Rayon `map_init` or `fold` isn't easy here.
            
            // Let's rely on OS file handle caching.
            let mut bam = match bam::IndexedReader::from_path(&bam_path) {
                Ok(b) => b,
                Err(e) => return format!("{}\t{}\t{}\t{}\t0\t0\t0.00\t0.00\t0.00\t0.0000\t1.0000\t# Error: {}", 
                    var.chrom, var.pos+1, var.ref_base, var.alt_base, e),
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
                        None => return format!("{}\t{}\t{}\t{}\t0\t0\t0.00\t0.00\t0.00\t0.0000\t1.0000",
                            var.chrom, var.pos+1, var.ref_base, var.alt_base),
                    }
                }
            };
            
            // Fetch region
            let start = var.pos as u64;
            let end = var.pos as u64 + 1;
            if let Err(_) = bam.fetch((tid, start, end)) {
                return format!("{}\t{}\t{}\t{}\t0\t0\t0.00\t0.00\t0.00\t0.0000\t1.0000",
                            var.chrom, var.pos+1, var.ref_base, var.alt_base);
            }
                 
            let mut mut_lengths = Vec::new();
            let mut wt_lengths = Vec::new();
            
            // Process reads
            for result in bam.records() {
                let record = match result {
                    Ok(r) => r,
                    Err(_) => continue,
                };
                
                // Filters
                if record.mapq() < map_quality { continue; }
                if record.is_duplicate() { continue; }
                if record.is_unmapped() { continue; }
                if record.is_secondary() { continue; }
                if record.is_supplementary() { continue; }

                // Get base
                let read_start = record.pos() as i64; 
                if var.pos < read_start { continue; }
                let target_offset = (var.pos - read_start) as usize;
                
                // Walk CIGAR
                let mut ref_offset = 0usize;
                let mut query_offset = 0usize;
                let mut found = false;
                let mut base_char = 'N';
                
                for op in record.cigar().iter() {
                    match op {
                        Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                            let len = *len as usize;
                            if target_offset >= ref_offset && target_offset < ref_offset + len {
                                let dist = target_offset - ref_offset;
                                let q_idx = query_offset + dist;
                                let seq = record.seq();
                                if q_idx < seq.len() {
                                    base_char = seq[q_idx] as char;
                                }
                                found = true;
                                break;
                            }
                            ref_offset += len;
                            query_offset += len;
                        },
                        Cigar::Ins(len) => { query_offset += *len as usize; },
                        Cigar::Del(len) | Cigar::RefSkip(len) => {
                            let len = *len as usize;
                            if target_offset >= ref_offset && target_offset < ref_offset + len { break; }
                            ref_offset += len;
                        },
                        Cigar::SoftClip(len) => { query_offset += *len as usize; },
                        Cigar::HardClip(_) | Cigar::Pad(_) => {},
                    }
                }
                
                if found {
                    let tlen = record.insert_size().abs();
                    let len = if tlen == 0 { record.seq().len() as i64 } else { tlen };
                    
                    let b_str = base_char.to_string().to_uppercase();
                    let ref_upper = var.ref_base.to_uppercase();
                    let alt_upper = var.alt_base.to_uppercase();
                    
                    if b_str == alt_upper {
                        mut_lengths.push(len as f64);
                    } else if b_str == ref_upper {
                        wt_lengths.push(len as f64);
                    }
                }
            }
            
            // Statistics
            let n_mut = mut_lengths.len();
            let n_wt = wt_lengths.len();
            let mut mut_mean = 0.0;
            let mut wt_mean = 0.0;
            let mut delta = 0.0;
            let mut ks_stat = 0.0;
            let mut ks_pval = 1.0;
            
            if n_mut > 0 { mut_mean = mut_lengths.iter().sum::<f64>() / n_mut as f64; }
            if n_wt > 0 { wt_mean = wt_lengths.iter().sum::<f64>() / n_wt as f64; }
            if n_mut > 0 && n_wt > 0 { delta = wt_mean - mut_mean; }
            
            // KS Test
            if n_mut > 0 && n_wt > 0 {
                mut_lengths.sort_by(|a, b| a.partial_cmp(b).unwrap());
                wt_lengths.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let mut i = 0;
                let mut j = 0;
                let mut d_max: f64 = 0.0;
                let mut cdf1: f64 = 0.0;
                let mut cdf2: f64 = 0.0;
                
                while i < n_mut && j < n_wt {
                    let v1 = mut_lengths[i];
                    let v2 = wt_lengths[j];
                    if v1 < v2 { cdf1 += 1.0 / n_mut as f64; i += 1; }
                    else if v2 < v1 { cdf2 += 1.0 / n_wt as f64; j += 1; }
                    else { cdf1 += 1.0 / n_mut as f64; cdf2 += 1.0 / n_wt as f64; i += 1; j += 1; }
                    d_max = d_max.max((cdf1 - cdf2).abs());
                }
                ks_stat = d_max;
                let m = n_mut as f64;
                let n = n_wt as f64;
                let en = (m * n) / (m + n);
                let lambda = (en.sqrt() + 0.12 + 0.11 / en.sqrt()) * ks_stat;
                let mut p_val = 0.0;
                for k in 1..101 {
                    let term = (-1.0_f64).powi(k - 1) * (-2.0 * k as f64 * k as f64 * lambda * lambda).exp();
                    p_val += term;
                }
                ks_pval = (2.0 * p_val).clamp(0.0, 1.0);
            }
            
            pb.inc(1);
            
            format!("{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.4}\t{:.4}",
                var.chrom, var.pos+1, var.ref_base, var.alt_base,
                n_mut, n_wt, mut_mean, wt_mean, delta, ks_stat, ks_pval)
        })
        .collect();
        
    pb.finish_with_message("Done!");

    // 3. Write Output
    info!("Writing output to {:?}...", output_file);
    let mut out_file = File::create(&output_file)?;
    writeln!(out_file, "Chrom\tPos\tRef\tAlt\tMut_Count\tWT_Count\tMut_MeanSize\tWT_MeanSize\tDelta_Size\tKS_Stat\tKS_Pval")?;
    for line in results {
        writeln!(out_file, "{}", line)?;
    }

    Ok(())
}
