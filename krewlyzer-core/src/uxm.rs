use pyo3::prelude::*;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::Aux;
use rayon::prelude::*;
use indicatif::{ProgressBar, ProgressStyle};
use std::time::Duration;
use log::info;

/// Calculate UXM (Fragment-level Methylation)
/// 
/// # Arguments
/// * `bam_path` - Path to the input BAM file
/// * `marker_path` - Path to marker BED file
/// * `output_file` - Path to output TSV
/// * `map_quality` - Min mapping quality
/// * `min_cpg` - Min CpG count per fragment
/// * `methy_threshold` - Threshold for 'M'
/// * `unmethy_threshold` - Threshold for 'U'
/// * `pe_type` - "SE" or "PE"
#[pyfunction]
pub fn calculate_uxm(
    bam_path: PathBuf,
    marker_path: PathBuf,
    output_file: PathBuf,
    map_quality: u8,
    min_cpg: u32,
    methy_threshold: f64,
    unmethy_threshold: f64,
    _pe_type: String, // Prefixed with _ as currently unused
) -> PyResult<()> {
    // 1. Load Markers
    struct Marker {
        chrom: String,
        start: u64,
        end: u64,
    }
    // SAFETY: Share BAM path strings.
    
    let mut markers = Vec::new();
    let file = File::open(&marker_path)?;
    let reader = BufReader::new(file);
    
    info!("Loading markers from {:?}...", marker_path);

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 { continue; }
        
        markers.push(Marker {
            chrom: fields[0].to_string(),
            start: fields[1].parse().unwrap_or(0),
            end: fields[2].parse().unwrap_or(0),
        });
    }
    
    let total_markers = markers.len();
    info!("Found {} markers. Processing in parallel...", total_markers);
    
    // Progress Bar
    let pb = ProgressBar::new(total_markers as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"));
    pb.enable_steady_tick(Duration::from_millis(100));

    // 2. Process Markers (Parallel)
    let results: Vec<String> = markers.par_iter()
        .map(|marker| {
            let region_str = format!("{}:{}-{}", marker.chrom, marker.start, marker.end);
            
             // Thread-local BAM reader
            let mut bam = match bam::IndexedReader::from_path(&bam_path) {
                Ok(b) => b,
                Err(e) => return format!("{}\t0\t0\t0\t# Error: {}", region_str, e),
            };

            // Get chromosome tid
            let tid = match bam.header().tid(marker.chrom.as_bytes()) {
                Some(t) => t,
                None => return format!("{}\t0\t0\t0", region_str),
            };
            
            // Fetch region
            if bam.fetch((tid, marker.start, marker.end)).is_err() {
                return format!("{}\t0\t0\t0", region_str);
            }
            
            let mut u_frag = 0;
            let mut x_frag = 0;
            let mut m_frag = 0;
            
            for result in bam.records() {
                let record = match result {
                    Ok(r) => r,
                    Err(_) => continue,
                };

                // Apply filters
                if record.mapq() < map_quality { continue; }
                if record.is_duplicate() { continue; }
                if record.is_unmapped() { continue; }
                if record.is_secondary() { continue; }
                if record.is_supplementary() { continue; }
                
                // Get XM tag
                let mut num_meth = 0u32;
                let mut num_unmeth = 0u32;
                
                if let Ok(aux) = record.aux(b"XM") {
                    match aux {
                        Aux::String(s) => {
                            for &b in s.as_bytes() {
                                if b == b'Z' { num_meth += 1; }
                                else if b == b'z' { num_unmeth += 1; }
                            }
                        },
                        _ => {}
                    }
                }
                
                let total_cpg = num_meth + num_unmeth;
                if total_cpg < min_cpg { continue; }
                
                let ratio = num_meth as f64 / total_cpg as f64;
                if ratio >= methy_threshold { m_frag += 1; }
                else if ratio <= unmethy_threshold { u_frag += 1; }
                else { x_frag += 1; }
            }
            
            pb.inc(1);
            
            let total = u_frag + x_frag + m_frag;
            if total == 0 {
                format!("{}\t0\t0\t0", region_str)
            } else {
                let u_rat = u_frag as f64 / total as f64;
                let x_rat = x_frag as f64 / total as f64;
                let m_rat = m_frag as f64 / total as f64;
                format!("{}\t{:.6}\t{:.6}\t{:.6}", region_str, u_rat, x_rat, m_rat)
            }
        })
        .collect();
        
    pb.finish_with_message("Done!");

    // 3. Write Output
    info!("Writing output to {:?}...", output_file);
    let mut out_file = File::create(&output_file)?;
    writeln!(out_file, "region\tU\tX\tM")?;
    for line in results {
        writeln!(out_file, "{}", line)?;
    }

    Ok(())
}
