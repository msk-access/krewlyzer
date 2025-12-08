//! Windowed Protection Score (WPS) calculation
//!
//! Calculates WPS, Coverage, and Starts for transcript regions.
//! Writes output directly to gzipped TSV files per region.

use std::path::Path;
use std::io::{BufRead, BufReader, Write};
use std::fs::File;
use anyhow::{Result, Context};
use rayon::prelude::*;
use pyo3::prelude::*;
use flate2::write::GzEncoder;
use flate2::Compression;

/// Transcript region from TSV
#[derive(Debug, Clone)]
struct Region {
    id: String,
    chrom: String,
    start: u64,
    end: u64,
    strand: String,
}

/// Fragment range (0-based start, 0-based exclusive end)
#[derive(Debug, Clone, Copy)]
struct Fragment {
    // We use u32 to save memory, assuming chromosomes < 4GB
    // We store chrom as index or hash? For simplicity string for now, or match on string in loop
    // To save memory for 100M fragments, we should map chrom string to u8/u16 ID.
    chrom_id: u8, 
    start: u32,
    end: u32,
}

/// Parse a BGZF BED file into a vector of Fragments (lighter memory usage)
/// Returns (fragments, chrom_map)
fn parse_fragments_wps(bedgz_path: &Path) -> Result<(Vec<Fragment>, Vec<String>)> {
    use noodles_bgzf as bgzf;
    
    let file = File::open(bedgz_path)
        .with_context(|| format!("Failed to open BED.gz file: {:?}", bedgz_path))?;
    let bgzf_reader = bgzf::Reader::new(file);
    let reader = BufReader::new(bgzf_reader);
    
    let mut fragments = Vec::with_capacity(1_000_000);
    let mut chrom_map: Vec<String> = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 { continue; }
        
        let chrom = fields[0]; // normalized later?
        let start: u32 = fields[1].parse().unwrap_or(0);
        let end: u32 = fields[2].parse().unwrap_or(0);
        
        // Find or add chrom ID
        let chrom_id = if let Some(pos) = chrom_map.iter().position(|c| c == chrom) {
            pos as u8
        } else {
            if chrom_map.len() >= 255 {
                // Warning or error? access-solid usually has standard chroms
                // Just fallback to 255? or error
                eprintln!("Warning: Too many unique chromosomes (>255), skipping {}", chrom);
                continue; 
            }
            chrom_map.push(chrom.to_string());
            (chrom_map.len() - 1) as u8
        };
        
        fragments.push(Fragment { chrom_id, start, end });
    }
    
    Ok((fragments, chrom_map))
}

/// Parse TSV transcript file
fn parse_regions(tsv_path: &Path) -> Result<Vec<Region>> {
    let file = File::open(tsv_path).with_context(|| "Failed to open regions file")?;
    let reader = BufReader::new(file);
    let mut regions = Vec::new();
    
    let valid_chroms: Vec<String> = (1..=22).map(|i| i.to_string()).chain(vec!["X".to_string(), "Y".to_string()]).collect();
    
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() { continue; }
        
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 5 { continue; }
        
        let id = fields[0].to_string();
        let chrom_raw = fields[1];
        let start_str = fields[2];
        let end_str = fields[3];
        let strand = fields[4].to_string();
        
        // Normalize chrom: remove chr prefix
        let chrom_norm = chrom_raw.trim_start_matches("chr").to_string();
        
        // Filter valid chroms (1-22, X, Y) if desired (wps.py does this)
        if !valid_chroms.iter().any(|c| c == &chrom_norm) {
            continue;
        }
        
        let start: u64 = start_str.parse::<f64>().unwrap_or(0.0) as u64; // wps.py does int(float(str))
        let end: u64 = end_str.parse::<f64>().unwrap_or(0.0) as u64;
        
        if start < 1 { continue; }
        
        regions.push(Region { id, chrom: chrom_norm, start, end, strand });
    }
    
    Ok(regions)
}

#[pyfunction]
#[pyo3(signature = (bedgz_path, tsv_path, output_dir, file_stem, empty=false, protect=120, min_size=120, max_size=180, total_fragments=None))]
pub fn calculate_wps(
    _py: Python<'_>,
    bedgz_path: &str,
    tsv_path: &str,
    output_dir: &str,
    file_stem: &str,
    empty: bool,
    protect: i32,
    min_size: u32,
    max_size: u32,
    total_fragments: Option<u64>,
) -> PyResult<usize> {
    let bed_path = Path::new(bedgz_path);
    let tsv = Path::new(tsv_path);
    let out = Path::new(output_dir);
    
    // Parse inputs
    let (fragments, chrom_map) = parse_fragments_wps(bed_path)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to parse fragments: {}", e)))?;
        
    let regions = parse_regions(tsv)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to parse regions: {}", e)))?;

    let protection = (protect as u32) / 2;

    // Process parallel
    let processed_count = regions.par_iter().map(|region| {
        // Resolve region chrom to ID
        let region_chrom_norm = region.chrom.trim_start_matches("chr");
        let chrom_id_opt = chrom_map.iter().position(|c| c.trim_start_matches("chr") == region_chrom_norm);
        
        if chrom_id_opt.is_none() {
            return 0; // Chromosome not in BED file
        }
        let chrom_id = chrom_id_opt.unwrap() as u8;
        
        // Region coordinates (1-based inclusive)
        let r_start = region.start;
        let r_end = region.end;
        let length = (r_end - r_start + 1) as usize;
        
        let mut cov_arr = vec![0u32; length];
        let mut start_arr = vec![0u32; length];
        let mut gcount_arr = vec![0u32; length];
        let mut total_arr = vec![0u32; length];
        
        // Fetch candidates (optimization: fragments sorted ideally, but we iterate all matching chrom)
        // If slow, optimize by sorting fragments and binary searching.
        // For standard targeted panels, fragments vector is ~1M, filtering by chrom reduces to ~50k. Linear scan is fast enough.
        
        // Define fetch bounds: region expanded by protection
        let fetch_start = r_start.saturating_sub(protection as u64 + 1);
        let fetch_end = r_end + protection as u64;
        
        for frag in &fragments {
            if frag.chrom_id != chrom_id { continue; }
            
            // Overlap check (fragment 0-based [start, end))
            // Convert to 1-based [start+1, end] for logic comparison with fetch window
            let f_start_1 = frag.start as u64 + 1;
            let f_end_1 = frag.end as u64;
            
            // Check if fragment overlaps the extended fetch window
            if f_end_1 <= fetch_start || f_start_1 > fetch_end {
                continue;
            }
            
            let len = f_end_1 - f_start_1 + 1;
            if len < min_size as u64 || len > max_size as u64 {
                continue;
            }
            
            // 1. Coverage
            // Overlap of fragment [f_start_1, f_end_1] with region [r_start, r_end]
            let ov_start = r_start.max(f_start_1);
            let ov_end = r_end.min(f_end_1);
            
            if ov_start <= ov_end {
                let s = (ov_start - r_start) as usize;
                let e = (ov_end - r_start + 1) as usize;
                for i in s..e { cov_arr[i] += 1; }
            }
            
            // 2. Starts (ends)
            if f_start_1 >= r_start && f_start_1 <= r_end {
                start_arr[(f_start_1 - r_start) as usize] += 1;
            }
            if f_end_1 >= r_start && f_end_1 <= r_end {
                start_arr[(f_end_1 - r_start) as usize] += 1;
            }
            
            // 3. WPS
            // gcount (spanning): window [k-P, k+P] inside read
            // Valid range for k (relative to read): [f_start_1 + P, f_end_1 - P]
            if (f_end_1 - f_start_1) >= (2 * protection as u64) {
                 let g_start = f_start_1 + protection as u64;
                 let g_end = f_end_1 - protection as u64;
                 
                 let g_ov_start = r_start.max(g_start);
                 let g_ov_end = r_end.min(g_end);
                 
                 if g_ov_start <= g_ov_end {
                     let s = (g_ov_start - r_start) as usize;
                     let e = (g_ov_end - r_start + 1) as usize;
                     for i in s..e { gcount_arr[i] += 1; }
                 }
            }
            
            // total (overlapping): window [k-P, k+P] overlaps read
            // Valid range for k: [f_start_1 - P, f_end_1 + P]
            // Note: f_start_1 - P might be < 1.
            let t_start = f_start_1.saturating_sub(protection as u64);
            let t_end = f_end_1 + protection as u64;
            
            let t_ov_start = r_start.max(t_start);
            let t_ov_end = r_end.min(t_end);
            
            if t_ov_start <= t_ov_end {
                let s = (t_ov_start - r_start) as usize;
                let e = (t_ov_end - r_start + 1) as usize;
                for i in s..e { total_arr[i] += 1; }
            }
        }
        
        // Compute WPS
        let mut wps_arr = Vec::with_capacity(length);
        for i in 0..length {
            let val = 2 * (gcount_arr[i] as i32) - (total_arr[i] as i32);
            wps_arr.push(val);
        }
        
        // Write logic
        // Skip empty?
        let total_cov: u32 = cov_arr.iter().sum();
        if total_cov == 0 && !empty {
            return 0;
        }
        
        // Prepare output data
        struct Row {
            pos: u64,
            cov: u32,
            starts: u32,
            wps: i32,
            wps_norm: f64,
        }
        
        let mut rows = Vec::with_capacity(length);
        for i in 0..length {
            let wps_norm = if let Some(total) = total_fragments {
                wps_arr[i] as f64 / (total as f64 / 1_000_000.0)
            } else {
                0.0
            };
            
            rows.push(Row {
                pos: r_start + i as u64,
                cov: cov_arr[i],
                starts: start_arr[i],
                wps: wps_arr[i],
                wps_norm,
            });
        }
        
        // Handle strand: if "-", reverse rows
        if region.strand == "-" {
            rows.reverse();
        }
        
        // Write file
        // Pattern: file_stem.id.WPS.tsv.gz
        let filename = format!("{}.{}.WPS.tsv.gz", file_stem, region.id);
        let path = out.join(filename);
        
        if let Ok(file_out) = File::create(&path) {
            let mut gz = GzEncoder::new(file_out, Compression::default());
            
            for row in rows {
                if let Some(_) = total_fragments {
                    writeln!(gz, "{}\t{}\t{}\t{}\t{}\t{:.6}", region.chrom, row.pos, row.cov, row.starts, row.wps, row.wps_norm).ok();
                } else {
                    writeln!(gz, "{}\t{}\t{}\t{}\t{}", region.chrom, row.pos, row.cov, row.starts, row.wps).ok();
                }
            }
            gz.finish().ok();
        } else {
             return 0; // Error creating file
        }
        
        1 // Processed
    }).sum();
    
    Ok(processed_count)
}
