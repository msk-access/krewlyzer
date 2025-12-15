use pyo3::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader}; // Needed for blacklist reading
use rust_htslib::bam::{self, Read};
use rust_htslib::faidx;
use rayon::prelude::*;
use indicatif::{ProgressBar, ProgressStyle};
use std::time::Duration;
use log::{info};

// Config struct moved to avoid duplication

struct Chunk {
    tid: u32,
    chrom: String,
    start: u64,
    end: u64,
}

// Result from a single chunk
struct ChunkResult {
    fragments: Vec<String>,
    end_motifs: HashMap<String, u64>,
    bp_motifs: HashMap<String, u64>,
    count: u64,
}

impl ChunkResult {
    fn new() -> Self {
        Self {
            fragments: Vec::new(),
            end_motifs: HashMap::new(),
            bp_motifs: HashMap::new(),
            count: 0,
        }
    }
}

// Blacklist helpers
fn load_exclude_regions(path: &str) -> HashSet<(String, u64, u64)> {
    let mut regions = HashSet::new();
    if let Ok(file) = File::open(path) {
        let reader = BufReader::new(file);
        for line in reader.lines().flatten() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 {
                if let (Ok(start), Ok(end)) = (fields[1].parse(), fields[2].parse()) {
                    regions.insert((fields[0].to_string(), start, end));
                }
            }
        }
    }
    regions
}

fn overlaps_exclude(chrom: &str, start: u64, end: u64, exclude: &HashSet<(String, u64, u64)>) -> bool {
    for (ex_chrom, ex_start, ex_end) in exclude {
        if chrom == ex_chrom && start < *ex_end && end > *ex_start {
            return true;
        }
    }
    false
}

/// Calculate GC content
fn calculate_gc(seq: &[u8]) -> f64 {
    let mut gc = 0;
    let mut valid = 0;
    for &b in seq {
        match b {
            b'G' | b'g' | b'C' | b'c' => { gc += 1; valid += 1; },
            b'A' | b'a' | b'T' | b't' => { valid += 1; },
            _ => {}
        }
    }
    if valid == 0 { 0.0 } else { gc as f64 / valid as f64 }
}

/// Reverse complement
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|b| match b {
        b'A' => b'T', b'a' => b'T',
        b'T' => b'A', b't' => b'A',
        b'G' => b'C', b'g' => b'C',
        b'C' => b'G', b'c' => b'G',
        x => *x,
    }).collect()
}

#[derive(Clone)]
pub struct UnifiedConfig {
    pub mapq: u8,
    pub min_len: u32,
    pub max_len: u32,
    pub kmer: usize,
    pub count_motifs: bool,
    pub output_bed: bool,
    pub skip_duplicates: bool,
    pub require_proper_pair: bool,
}

/// Unified Parallel Engine
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (bam_path, fasta_path, mapq=20, min_len=65, max_len=400, kmer=4, threads=0, output_bed_path=None, output_motif_prefix=None, exclude_path=None, skip_duplicates=true, require_proper_pair=true))]
pub fn process_bam_parallel(
    bam_path: String,
    fasta_path: String,
    mapq: u8,
    min_len: u32,
    max_len: u32,
    kmer: usize,
    threads: usize,
    output_bed_path: Option<String>,
    output_motif_prefix: Option<String>,
    exclude_path: Option<String>,
    skip_duplicates: bool,
    require_proper_pair: bool,
) -> PyResult<(u64, HashMap<String, u64>, HashMap<String, u64>)> {
    
    // 1. Configure Global Thread Pool if needed
    if threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new().num_threads(threads).build_global();
    }
    
    // Load blacklist (Main thread)
    let exclude_regions = match exclude_path {
        Some(p) => load_exclude_regions(&p),
        None => HashSet::new(),
    };
    // To share with Rayon, we wrap in Arc
    let exclude_arc = std::sync::Arc::new(exclude_regions);

    // Load available chromosomes from FASTA to avoid fetch crashes
    let valid_chroms = if let Ok(fa) = faidx::Reader::from_path(&fasta_path) {
        let n = fa.n_seqs();
        let mut s = HashSet::new();
        for i in 0..n {
             if let Ok(name) = fa.seq_name(i as i32) {
                 s.insert(name);
             }
        }
        Some(s)
    } else {
        None
    };
    let valid_chroms_arc = std::sync::Arc::new(valid_chroms);

    let config = UnifiedConfig {
        mapq,
        min_len,
        max_len,
        kmer,
        count_motifs: output_motif_prefix.is_some(),
        output_bed: output_bed_path.is_some(),
        skip_duplicates,
        require_proper_pair,
    };
    
    // 2. Scan Header for Chunks
    // We need to read the BAM header once to define chunks
    let bam_reader = bam::IndexedReader::from_path(&bam_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open BAM: {}", e)))?;
    
    let header = bam_reader.header();
    let mut chunks = Vec::new();
    let chunk_size = 10_000_000; // 10Mb chunks
    
    for (tid, name_bytes) in header.target_names().iter().enumerate() {
        let name = String::from_utf8_lossy(name_bytes).to_string();
        let len = header.target_len(tid as u32).unwrap_or(0);
        
        let mut start = 0;
        while start < len {
            let end = (start + chunk_size).min(len);
            chunks.push(Chunk {
                tid: tid as u32,
                chrom: name.clone(),
                start,
                end,
            });
            start = end;
        }
    }
    
    info!("Split genome into {} chunks. Processing...", chunks.len());
    
    // Progress Bar
    let pb = ProgressBar::new(chunks.len() as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"));
    pb.enable_steady_tick(Duration::from_millis(100));

    // 3. Parallel Processing
    let results: Vec<ChunkResult> = chunks.par_iter().map(|chunk| {
        let mut result = ChunkResult::new();
        
        // Thread-local Readers
        let mut local_bam = match bam::IndexedReader::from_path(&bam_path) {
            Ok(b) => b,
            Err(_) => return result, // Skip chunk on error
        };
        
        // Only open FASTA if needed (for GC or Breakpoint Motif)
        let local_fasta = if config.output_bed || (config.count_motifs && kmer > 0) {
            faidx::Reader::from_path(&fasta_path).ok()
        } else { None };

        // Fetch Region
        if local_bam.fetch((chunk.tid, chunk.start, chunk.end)).is_err() {
            return result;
        }

        // Check if chromosome is in reference
        let chrom_valid = if let Some(ref set) = *valid_chroms_arc {
             set.contains(&chunk.chrom)
        } else { false };
        
        // Only fetch if chrom is valid AND we have a reader
        let can_fetch_seq = chrom_valid && local_fasta.is_some();
        
        let context_len = (config.kmer as f64 / 2.0).ceil() as usize;

        for record_res in local_bam.records() {
            let record = match record_res { Ok(r) => r, Err(_) => continue };
            
            // Standard Filters
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() 
               || record.is_quality_check_failed() { continue; }
            
            if config.skip_duplicates && record.is_duplicate() { continue; }
            if record.mapq() < config.mapq { continue; }
            
            if config.require_proper_pair {
                 if !record.is_paired() || !record.is_proper_pair() { continue; }
            }
            
            // Note: We iterate ALL reads in the chunk.
            // But chunks overlap (boundary issues).
            // Logic: process read only if its start position is within [chunk.start, chunk.end)
            if (record.pos() as u64) < chunk.start || (record.pos() as u64) >= chunk.end {
                continue;
            }

            // --- Extract Logic (BED Output) ---
            // Only process R1 for Fragment definition to avoid double counting
            if config.output_bed && record.is_first_in_template() {
                 let tlen = record.insert_size().abs() as i64;
                 if tlen >= config.min_len as i64 && tlen <= config.max_len as i64 {
                     let start = record.pos() as u64; // 0-based
                     let end = start + tlen as u64;
                     
                     // Blacklist Check
                     if !exclude_arc.is_empty() && overlaps_exclude(&chunk.chrom, start, end, &exclude_arc) {
                         continue;
                     }
                     
                     // Calculate GC
                     let gc = if can_fetch_seq {
                         if let Some(ref fa) = local_fasta {
                             match fa.fetch_seq(&chunk.chrom, start as usize, (end - 1) as usize) {
                                 Ok(seq) => calculate_gc(&seq),
                                 Err(_) => 0.0,
                             }
                         } else { 0.0 }
                     } else { 0.0 };
                     
                     result.fragments.push(format!("{}\t{}\t{}\t{:.4}", chunk.chrom, start, end, gc));
                     result.count += 1;
                 }
            }
            
            // --- Motif Logic ---
            if config.count_motifs {
                // ... (Motif logic needs independent blacklist check? Or relies on R1 check?)
                // Motif tool counts *individual fragment ends*.
                // If the fragment is blacklisted, we should skip counting its motifs.
                // But since we process R1 and R2 *independently*, we need to know the fragment coordinates involved.
                // TLEN gives us length. pos gives us one end.
                // R1: [pos, pos+tlen]
                // R2: [mate_pos, mate_pos+tlen] -> [pos, pos+tlen] (if proper pair)
                // We can reconstruct fragment coordinates from EITHER read using TLEN.
                // But we must be careful about `tlen` sign.
                // R1 (fwd): tlen > 0. Frag: [pos, pos+tlen]
                // R2 (rev): tlen < 0. Frag: [mate_pos, mate_pos+abs(tlen)] => [end-abs(tlen), end]. 
                // Where `end` of R2 ALIGNMENT is the end of the fragment?
                // R2 alignment: [pos, end]. The fragment end is `end`. The fragment start is `end - abs(tlen)`.
                
                let tlen = record.insert_size();
                let abs_len = tlen.abs() as i64;
                if abs_len < config.min_len as i64 || abs_len > config.max_len as i64 { continue; }
                
                // Reconstruct fragment coords for blacklist check
                let (frag_start, frag_end) = if tlen > 0 {
                    (record.pos() as u64, (record.pos() + tlen) as u64)
                } else {
                    // tlen < 0. record is R2 (usually).
                    // The mate is at `record.mate_pos()`.
                    // But checking blacklist requires `start, end`.
                    // Frag End is `record.cigar().end_pos()`.
                    // Frag Start is `Frag End - abs_len`.
                    let e = record.cigar().end_pos() as u64;
                    (e.saturating_sub(abs_len as u64), e)
                };
                
                if !exclude_arc.is_empty() && overlaps_exclude(&chunk.chrom, frag_start, frag_end, &exclude_arc) {
                     continue;
                }
                
                // Get sequence
                let seq = record.seq();
                if seq.len() < config.kmer { continue; }
                
                let mut motif_seq: Vec<u8> = Vec::new();
                let genomic_coords; 
                
                if !record.is_reverse() {
                     for i in 0..config.kmer { motif_seq.push(seq[i]); }
                     genomic_coords = (record.pos(), record.pos());
                } else {
                     let len = seq.len();
                     let raw_end = &seq.as_bytes()[len-config.kmer..];
                     motif_seq = reverse_complement(raw_end); 
                     genomic_coords = (record.cigar().end_pos(), record.cigar().end_pos());
                }
                
                // Store End Motif
                let motif_str = String::from_utf8_lossy(&motif_seq).to_string();
                *result.end_motifs.entry(motif_str.clone()).or_default() += 1;
                
                // Breakpoint Motif ... (retained)
                // Needs FASTA fetch
                if can_fetch_seq { // Safe guard
                    if let Some(ref fa) = local_fasta {
                    // Logic for Breakpoint context needs careful coordinate handling matching `motif.rs`
                    // Simplified for now: strictly context + motif
                    // If FWD: Ref[pos-p..pos] + Frag[0..p]
                    // If REV: Frag[end-p..end] (RC) + Ref[end..end+p]
                    
                    // Implementation skipped for brevity in this step, focusing on framework.
                    // But we MUST implement it if we claim parity.
                    // For now, let's placeholder "BP logic here".
                    // Actually, "krewlyzer motif" relies on BP motif. I should include it.
                    
                    let p = context_len;
                    let bp_str = if !record.is_reverse() {
                        let pos = genomic_coords.0 as usize;
                        let ctx_start = pos.saturating_sub(p);
                        let ref_seq = fa.fetch_seq(&chunk.chrom, ctx_start, pos.saturating_sub(1)).unwrap_or_default();
                        let frag_part = &motif_seq[0..p.min(motif_seq.len())];
                        let full = [ref_seq, frag_part.to_vec()].concat();
                        String::from_utf8_lossy(&full).to_string()
                    } else {
                        // Reverse
                        let pos = genomic_coords.0 as usize; // end of mapping
                        let ref_seq = fa.fetch_seq(&chunk.chrom, pos, pos + p - 1).unwrap_or_default();
                        let frag_part = &motif_seq[0..p.min(motif_seq.len())]; // acts as RC of end
                        let full = [frag_part.to_vec(), ref_seq].concat();
                        String::from_utf8_lossy(&full).to_string()
                    };
                    
                     *result.bp_motifs.entry(bp_str).or_default() += 1;
                }
                } // End safe guard
            }
        }
        
        pb.inc(1);
        result
    }).collect();
    
    pb.finish_with_message("Done processing chunks.");
    
    // 4. Merge Results
    info!("Merging results...");
    let mut total_count = 0;
    let mut final_end_motifs = HashMap::new();
    let mut final_bp_motifs = HashMap::new();
    
    // Write BED if requested
    if let Some(path) = output_bed_path {
        info!("Writing BED to {}...", path);
        let mut f = File::create(path).map(BufWriter::new)?;
        for res in &results {
            for line in &res.fragments {
                writeln!(f, "{}", line)?;
            }
            total_count += res.count;
        }
    } else {
        // Just sum counts if not writing output
         for res in &results {
            total_count += res.count;
        }
    }
    
    // Merge Motifs
    if config.count_motifs {
        for res in results {
            for (k, v) in res.end_motifs {
                *final_end_motifs.entry(k).or_default() += v;
            }
            for (k, v) in res.bp_motifs {
                *final_bp_motifs.entry(k).or_default() += v;
            }
        }
    }
    
    Ok((total_count, final_end_motifs, final_bp_motifs))
}
