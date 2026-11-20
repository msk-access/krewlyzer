//! BAM/CRAM Fragment Extraction with Motif Counting
//!
//! Extracts cfDNA fragments from BAM/CRAM files in a single pass:
//! - Fragment coordinates → BED.gz output
//! - End motif (EM) k-mer counts
//! - Breakpoint motif (BPM) k-mer counts
//! - GC content observations for bias correction
//!
//! For panel data, separates on-target and off-target counts.
//! Off-target is used for unbiased GC model training.

use pyo3::prelude::*;
use std::fs::File;
use std::io::Write;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use rust_htslib::bam::{self, Read};
use rust_htslib::faidx;
use rayon::prelude::*;
use indicatif::{ProgressBar, ProgressStyle};
use std::time::Duration;
use log::{info, debug};
use noodles::bgzf;

// Use shared utilities to avoid code duplication
use crate::motif_utils::{reverse_complement, calculate_gc};
use crate::contigs::{norm_chrom, ContigNames};

struct Chunk {
    tid: u32,
    chrom: String,
    /// The reference FASTA's name for `chrom`, used by every sequence fetch; `None` when the
    /// reference lacks the contig (a decoy or alt), so its fragments get no GC.
    fasta_name: Option<String>,
    start: u64,
    end: u64,
}

// Result from a single chunk
struct ChunkResult {
    /// (fragment start, BED line). The start is kept alongside so the
    /// writer can restore coordinate order: a reverse-oriented fragment
    /// begins before the read that produced it, so read order is no longer
    /// fragment order and tabix would reject the output.
    fragments: Vec<(u64, String)>,
    // Off-target motifs (primary - unbiased)
    end_motifs: HashMap<String, u64>,
    bp_motifs: HashMap<String, u64>,
    // On-target motifs (for comparison - PCR biased)
    end_motifs_on: HashMap<String, u64>,
    bp_motifs_on: HashMap<String, u64>,
    count: u64,
    // GC observations for correction factor computation
    // Key: (length_bin, gc_percent), Value: count
    gc_observations: HashMap<(u8, u8), u64>,         // Off-target (unbiased)
    gc_observations_ontarget: HashMap<(u8, u8), u64>, // On-target (for panel mode)
    // Fragments whose contig the reference lacks (their GC cannot be computed), and fragments
    // whose GC fetch failed on a contig it has. Neither enters the GC observations.
    n_missing_contig: u64,
    n_gc_fetch_failed: u64,
}

impl ChunkResult {
    fn new() -> Self {
        Self {
            fragments: Vec::new(),
            end_motifs: HashMap::new(),
            bp_motifs: HashMap::new(),
            end_motifs_on: HashMap::new(),
            bp_motifs_on: HashMap::new(),
            count: 0,
            gc_observations: HashMap::new(),
            gc_observations_ontarget: HashMap::new(),
            n_missing_contig: 0,
            n_gc_fetch_failed: 0,
        }
    }
}

/// Regions (chrom, start, end) from a BED, contig names without a leading `chr`, so `chr1` and `1`
/// match. Used for targets and for the blacklist; queries go through [`overlaps_region`], which
/// normalises the same way. The bundled hg19 blacklist names `chr1` while MSK-ACCESS BAMs name `1`:
/// with exact matching it matched nothing and was never applied.
fn load_regions(path: &str) -> HashSet<(String, u64, u64)> {
    let mut regions = HashSet::new();
    let path_buf = std::path::PathBuf::from(path);
    if let Ok(reader) = crate::bed::get_reader(&path_buf) {
        for line in reader.lines().map_while(Result::ok) {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 {
                if let (Ok(start), Ok(end)) = (fields[1].parse(), fields[2].parse()) {
                    regions.insert((norm_chrom(fields[0]).to_string(), start, end));
                }
            }
        }
    }
    regions
}

fn overlaps_region(chrom: &str, start: u64, end: u64, exclude: &HashSet<(String, u64, u64)>) -> bool {
    let chrom = norm_chrom(chrom);
    for (ex_chrom, ex_start, ex_end) in exclude {
        if chrom == ex_chrom && start < *ex_end && end > *ex_start {
            return true;
        }
    }
    false
}

// NOTE: reverse_complement and calculate_gc moved to motif_utils.rs

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
#[pyo3(signature = (bam_path, fasta_path, mapq=20, min_len=65, max_len=400, kmer=4, threads=0, output_bed_path=None, output_motif_prefix=None, exclude_path=None, target_regions_path=None, skip_duplicates=true, require_proper_pair=true, silent=false))]
pub fn process_bam_parallel(
    py: Python<'_>,
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
    target_regions_path: Option<String>,
    skip_duplicates: bool,
    require_proper_pair: bool,
    silent: bool,
) -> PyResult<(
    u64,                      // fragment count
    HashMap<String, u64>,     // end_motifs (off-target)
    HashMap<String, u64>,     // bp_motifs (off-target)
    HashMap<(u8, u8), u64>,   // gc_observations (off-target)
    HashMap<String, u64>,     // end_motifs_on (on-target)
    HashMap<String, u64>,     // bp_motifs_on (on-target)
    HashMap<(u8, u8), u64>,   // gc_observations_ontarget (on-target)
    HashMap<String, u64>,     // contig counts: fragments_on_missing_contig, gc_fetch_failed
)> {
    
    // 1. Configure Global Thread Pool if needed
    if threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new().num_threads(threads).build_global();
    }
    
    // Load blacklist (Main thread)
    let exclude_regions = match exclude_path {
        Some(p) => load_regions(&p),
        None => HashSet::new(),
    };
    // To share with Rayon, we wrap in Arc
    let exclude_arc = std::sync::Arc::new(exclude_regions);

    // Load target regions for off-target GC model (panel data like MSK-ACCESS)
    let target_regions = match target_regions_path {
        Some(ref p) => {
            let regions = load_regions(p);
            info!("Panel mode: GC model will use off-target reads only ({} target regions)", regions.len());
            regions
        },
        None => HashSet::new(),
    };
    let target_arc = std::sync::Arc::new(target_regions);
    let is_panel_mode = !target_arc.is_empty();

    // The reference's contigs (names and lengths), for the one BAM-to-FASTA rule in contigs.rs.
    // A reference that cannot be opened is an error: it used to leave every fragment at GC 0.
    let contig_names = faidx::Reader::from_path(&fasta_path)
        .map(|fa| ContigNames::from_faidx(&fa))
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
            format!("Failed to open reference FASTA {}: {}", fasta_path, e)))?;

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
    let bam_contigs: Vec<(String, u64)> = header.target_names().iter().enumerate()
        .map(|(tid, name)| (String::from_utf8_lossy(name).to_string(),
                            header.target_len(tid as u32).unwrap_or(0)))
        .collect();

    // Resolve every BAM contig before any read is processed. A primary contig the reference lacks,
    // or any contig whose length differs (another assembly), is an error; other absent contigs
    // (decoys, alts) are counted per fragment and reported by the caller.
    let fasta_names = contig_names.resolve_header(&bam_contigs)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("{} against reference {}: {}", bam_path, fasta_path, e)))?;
    let absent: Vec<&str> = bam_contigs.iter().zip(&fasta_names)
        .filter(|(_, f)| f.is_none()).map(|((n, _), _)| n.as_str()).collect();
    if !absent.is_empty() {
        info!("{} BAM contigs are not in the reference (decoys, alts, unplaced); their fragments get no GC: {}{}",
              absent.len(), absent[..absent.len().min(5)].join(", "), if absent.len() > 5 { ", ..." } else { "" });
    }

    let mut chunks = Vec::new();
    let chunk_size = 10_000_000; // 10Mb chunks

    for (tid, ((name, len), fasta_name)) in bam_contigs.iter().zip(fasta_names).enumerate() {
        let mut start = 0;
        while start < *len {
            let end = (start + chunk_size).min(*len);
            chunks.push(Chunk {
                tid: tid as u32,
                chrom: name.clone(),
                fasta_name: fasta_name.clone(),
                start,
                end,
            });
            start = end;
        }
    }
    
    info!("Split genome into {} chunks. Processing...", chunks.len());
    
    // Progress Bar (hidden when called from wrapper)
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
    let results: Vec<ChunkResult> = py.allow_threads(|| {
        chunks.par_iter().map(|chunk| {
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

        // The reference's name for this contig, for every fetch below; `None` when it lacks it.
        let fasta_name = chunk.fasta_name.as_deref();

        let context_len = (config.kmer as f64 / 2.0).ceil() as usize;

        for record_res in local_bam.records() {
            let record = match record_res { Ok(r) => r, Err(_) => continue };
            
            // Standard Filters
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() 
               || record.is_quality_check_failed() { continue; }
            
            if config.skip_duplicates && record.is_duplicate() { continue; }
            if record.mapq() < config.mapq { continue; }
            
            if config.require_proper_pair
                 && (!record.is_paired() || !record.is_proper_pair()) { continue; }
            
            // Note: We iterate ALL reads in the chunk.
            // But chunks overlap (boundary issues).
            // Logic: process read only if its start position is within [chunk.start, chunk.end)
            if (record.pos() as u64) < chunk.start || (record.pos() as u64) >= chunk.end {
                continue;
            }

            // --- The fragment, from either mate ---
            let tlen_signed = record.insert_size();
            let tlen = tlen_signed.unsigned_abs();
            if tlen < config.min_len as u64 || tlen > config.max_len as u64 { continue; }

            // The fragment spans the outermost bases of the pair, which is not `pos() + |tlen|`
            // when this mate is the rightmost one: there `pos()` is near the fragment's 3' end and
            // the interval lands ~(tlen - read_length) too far right. Consensus BAMs normalise R1
            // to forward so this is rare there (<1%), but uncollapsed and WGS input run ~50%.
            //
            // Branch on the TLEN sign, not on is_reverse(). htslib gives the rightmost segment of
            // a pair a negative TLEN, and "am I the rightmost mate" is exactly the question here.
            // The two agree for well-formed FR pairs and diverge for malformed ones -- real data
            // contains records flagged forward whose mate lies to their left.
            let (start, end) = if tlen_signed >= 0 {
                let s = record.pos() as u64;
                (s, s + tlen)
            } else {
                let e = record.cigar().end_pos() as u64;
                (e.saturating_sub(tlen), e)
            };

            if !exclude_arc.is_empty() && overlaps_region(&chunk.chrom, start, end, &exclude_arc) {
                continue;
            }
            let is_on_target = is_panel_mode && overlaps_region(&chunk.chrom, start, end, &target_arc);

            // --- Fragment count, contig counts, BED and GC: R1 only, so each pair counts once ---
            if record.is_first_in_template() {
                result.count += 1;
                if fasta_name.is_none() { result.n_missing_contig += 1; }

                if config.output_bed {
                    // None when the reference lacks the contig or the fetch fails. Such a fragment
                    // is still written (GC 0.0 in the BED, as before) but enters no GC observation:
                    // counted at 0% it would bias the model's lowest bin.
                    let gc = match (fasta_name, &local_fasta) {
                        (Some(name), Some(fa)) => match fa.fetch_seq(name, start as usize, (end - 1) as usize) {
                            Ok(seq) => Some(calculate_gc(&seq)),
                            Err(_) => { result.n_gc_fetch_failed += 1; None }
                        },
                        // The reference opened on the main thread but not in this one.
                        (Some(_), None) => { result.n_gc_fetch_failed += 1; None }
                        (None, _) => None,
                    };

                    result
                        .fragments
                        .push((start, format!("{}\t{}\t{}\t{:.4}", chunk.chrom, start, end, gc.unwrap_or(0.0))));

                    // GC observation for correction factor computation. Length bins: 5bp bins from
                    // 60-400bp = 68 bins (matches gc_correction.rs). Panel mode: separate on-target
                    // and off-target for dual GC models.
                    if let Some(gc) = gc {
                        // 5bp bins: (length - 60) / 5 -> 0..67
                        let length_bin = ((tlen.saturating_sub(60)) / 5).min(67) as u8;
                        let gc_pct = (gc * 100.0).round() as u8;
                        let observations = if is_on_target {
                            &mut result.gc_observations_ontarget
                        } else {
                            &mut result.gc_observations
                        };
                        *observations.entry((length_bin, gc_pct)).or_insert(0) += 1;
                    }
                }
            }
            
            // --- Motif Logic: every fragment end, so both mates ---
            if config.count_motifs {
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
                
                // Store End Motif (route to on-target or off-target)
                let motif_str = String::from_utf8_lossy(&motif_seq).to_string();
                if is_on_target {
                    *result.end_motifs_on.entry(motif_str.clone()).or_default() += 1;
                } else {
                    *result.end_motifs.entry(motif_str.clone()).or_default() += 1;
                }
                
                // Breakpoint motif: reference context + fragment end, so it needs the FASTA.
                // If FWD: Ref[pos-p..pos] + Frag[0..p]
                // If REV: Frag[end-p..end] (RC) + Ref[end..end+p]
                if let (Some(name), Some(fa)) = (fasta_name, &local_fasta) {
                    let p = context_len;
                    let bp_str = if !record.is_reverse() {
                        let pos = genomic_coords.0 as usize;
                        let ctx_start = pos.saturating_sub(p);
                        let ref_seq = fa.fetch_seq(name, ctx_start, pos.saturating_sub(1)).unwrap_or_default();
                        let frag_part = &motif_seq[0..p.min(motif_seq.len())];
                        let full = [ref_seq, frag_part.to_vec()].concat();
                        String::from_utf8_lossy(&full).to_string()
                    } else {
                        // Reverse
                        let pos = genomic_coords.0 as usize; // end of mapping
                        let ref_seq = fa.fetch_seq(name, pos, pos + p - 1).unwrap_or_default();
                        let frag_part = &motif_seq[0..p.min(motif_seq.len())]; // acts as RC of end
                        let full = [frag_part.to_vec(), ref_seq].concat();
                        String::from_utf8_lossy(&full).to_string()
                    };
                    
                    // Route BP motif to on-target or off-target
                    if is_on_target {
                        *result.bp_motifs_on.entry(bp_str).or_default() += 1;
                    } else {
                        *result.bp_motifs.entry(bp_str).or_default() += 1;
                    }
                }
            }
        }
        
        pb.inc(1);
        result
    }).collect()
    }); // end py.allow_threads
    
    pb.finish_with_message("Done processing chunks.");
    
    // 4. Merge Results
    info!("Merging results...");
    let mut total_count = 0;
    let (mut n_missing_contig, mut n_gc_fetch_failed) = (0u64, 0u64);
    let mut final_end_motifs = HashMap::new();
    let mut final_bp_motifs = HashMap::new();
    let mut final_gc_observations: HashMap<(u8, u8), u64> = HashMap::new();
    let mut final_gc_observations_ontarget: HashMap<(u8, u8), u64> = HashMap::new();

    for res in &results {
        total_count += res.count;
        n_missing_contig += res.n_missing_contig;
        n_gc_fetch_failed += res.n_gc_fetch_failed;
        for ((len_bin, gc_pct), count) in &res.gc_observations {
            *final_gc_observations.entry((*len_bin, *gc_pct)).or_insert(0) += count;
        }
        for ((len_bin, gc_pct), count) in &res.gc_observations_ontarget {
            *final_gc_observations_ontarget.entry((*len_bin, *gc_pct)).or_insert(0) += count;
        }
    }

    // Write BED if requested
    if let Some(path) = output_bed_path {
        info!("Writing BED to {}...", path);
        let file = File::create(path)?;
        let mut writer = bgzf::io::Writer::new(file);

        // Emit each chromosome's fragments in coordinate order.
        //
        // Chunks arrive in header order, and every chunk of a chromosome is
        // contiguous, so buffering one chromosome at a time is enough -- no
        // need to hold the genome. Sorting is required because a fragment whose
        // R1 is the rightmost mate starts before the read that produced it, by
        // up to (max_len - read_length). That is small next to a 10Mb chunk but
        // it does cross chunk boundaries, so a per-chunk sort would not be
        // enough. Unsorted output makes tabix indexing fail outright, which
        // takes every downstream feature with it.
        let mut pending: Vec<(u64, &String)> = Vec::new();
        let mut pending_chrom: Option<&str> = None;

        let flush = |buf: &mut Vec<(u64, &String)>,
                         w: &mut bgzf::io::Writer<File>|
         -> std::io::Result<()> {
            buf.sort_by_key(|(start, _)| *start);
            for (_, line) in buf.iter() {
                writeln!(w, "{}", line)?;
            }
            buf.clear();
            Ok(())
        };

        for (chunk, res) in chunks.iter().zip(results.iter()) {
            if pending_chrom != Some(chunk.chrom.as_str()) {
                flush(&mut pending, &mut writer)?;
                pending_chrom = Some(chunk.chrom.as_str());
            }
            pending.extend(res.fragments.iter().map(|(s, line)| (*s, line)));
        }
        flush(&mut pending, &mut writer)?;
        writer.finish()?;  // Ensure proper EOF marker
    }
    
    debug!("Total GC observation bins: {} off-target, {} on-target", 
           final_gc_observations.len(), final_gc_observations_ontarget.len());
    info!("Extracted {} fragments", total_count);
    
    // Merge Motifs (off-target and on-target separately)
    let mut final_end_motifs_on = HashMap::new();
    let mut final_bp_motifs_on = HashMap::new();
    
    if config.count_motifs {
        for res in results {
            // Off-target motifs
            for (k, v) in res.end_motifs {
                *final_end_motifs.entry(k).or_default() += v;
            }
            for (k, v) in res.bp_motifs {
                *final_bp_motifs.entry(k).or_default() += v;
            }
            // On-target motifs
            for (k, v) in res.end_motifs_on {
                *final_end_motifs_on.entry(k).or_default() += v;
            }
            for (k, v) in res.bp_motifs_on {
                *final_bp_motifs_on.entry(k).or_default() += v;
            }
        }
        
        // Log on-target counts if present
        let on_count: u64 = final_end_motifs_on.values().sum();
        let off_count: u64 = final_end_motifs.values().sum();
        if on_count > 0 {
            info!("Motif: {} off-target, {} on-target fragment ends", off_count, on_count);
        }
    }
    
    // Reported here; the policy (when to warn) is check_reference_contigs in core/gc_assets.py.
    info!("Contigs: {} fragments on contigs the reference lacks, {} failed GC fetches; neither enters a GC model",
          n_missing_contig, n_gc_fetch_failed);
    let guard: HashMap<String, u64> = [
        ("fragments_on_missing_contig".to_string(), n_missing_contig),
        ("gc_fetch_failed".to_string(), n_gc_fetch_failed),
    ].into_iter().collect();
    Ok((total_count, final_end_motifs, final_bp_motifs, final_gc_observations, final_end_motifs_on, final_bp_motifs_on, final_gc_observations_ontarget, guard))
}
