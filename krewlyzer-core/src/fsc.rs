//! Fragment Size Coverage (FSC) calculation
//!
//! Counts fragments in genomic bins by size category:
//! - Short: 65-150bp
//! - Intermediate: 151-260bp  
//! - Long: 261-400bp

use std::path::Path;
use std::io::{BufRead, BufReader};
use std::fs::File;
use anyhow::{Result, Context};
use rayon::prelude::*;
use pyo3::prelude::*;
use numpy::{PyArray1, IntoPyArray};

use crate::bed::Region;

/// Result of FSC calculation for a single bin
#[derive(Debug, Clone, Default)]
pub struct BinResult {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub short_count: u32,
    pub intermediate_count: u32,
    pub long_count: u32,
    pub total_count: u32,
    pub mean_gc: f64,
}

/// Parse a BED file to get regions (bins)
pub fn parse_bin_file(bin_path: &Path) -> Result<Vec<Region>> {
    let file = File::open(bin_path)
        .with_context(|| format!("Failed to open bin file: {:?}", bin_path))?;
    let reader = BufReader::new(file);
    
    let mut regions = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        
        let chrom = fields[0].to_string();
        let start: u64 = fields[1].parse().with_context(|| "Invalid start position")?;
        let end: u64 = fields[2].parse().with_context(|| "Invalid end position")?;
        
        regions.push(Region::new(chrom, start, end));
    }
    
    Ok(regions)
}

/// Parse a BGZF BED file using noodles-bgzf with libdeflate
fn parse_bedgz_fragments(bedgz_path: &Path) -> Result<Vec<(String, u64, u64, f64)>> {
    use noodles_bgzf as bgzf;
    
    let file = File::open(bedgz_path)
        .with_context(|| format!("Failed to open BED.gz file: {:?}", bedgz_path))?;
    let bgzf_reader = bgzf::Reader::new(file);
    let reader = BufReader::new(bgzf_reader);
    
    let mut fragments = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }
        
        let chrom = fields[0].to_string();
        let start: u64 = fields[1].parse().unwrap_or(0);
        let end: u64 = fields[2].parse().unwrap_or(0);
        let gc: f64 = fields[3].parse().unwrap_or(0.0);
        
        fragments.push((chrom, start, end, gc));
    }
    
    Ok(fragments)
}

/// Count fragments in regions (sequential but memory-efficient for targeted panels)
pub fn count_fragments_sequential(
    bedgz_path: &Path,
    regions: &[Region],
) -> Result<Vec<BinResult>> {
    // Load all fragments (fine for targeted panels, ~few million fragments)
    let fragments = parse_bedgz_fragments(bedgz_path)?;
    
    // Create results for each region
    let results: Vec<BinResult> = regions
        .par_iter()
        .map(|region| {
            let mut result = BinResult {
                chrom: region.chrom.clone(),
                start: region.start,
                end: region.end,
                ..Default::default()
            };
            
            let mut gc_sum = 0.0;
            let mut gc_count = 0u32;
            
            // Normalize chromosome names (handle chr prefix)
            let region_chrom_normalized = region.chrom.trim_start_matches("chr").to_string();
            
            for (frag_chrom, frag_start, frag_end, gc) in &fragments {
                let frag_chrom_normalized = frag_chrom.trim_start_matches("chr");
                
                // Check chromosome match
                if frag_chrom_normalized != region_chrom_normalized {
                    continue;
                }
                
                // Check overlap: fragment overlaps with region
                if *frag_end <= region.start || *frag_start >= region.end {
                    continue;
                }
                
                let length = frag_end.saturating_sub(*frag_start) as u32;
                
                // Count by size category
                if length >= 65 && length <= 400 {
                    result.total_count += 1;
                    gc_sum += gc;
                    gc_count += 1;
                    
                    if length <= 150 {
                        result.short_count += 1;
                    } else if length <= 260 {
                        result.intermediate_count += 1;
                    } else {
                        result.long_count += 1;
                    }
                }
            }
            
            result.mean_gc = if gc_count > 0 { gc_sum / gc_count as f64 } else { f64::NAN };
            
            result
        })
        .collect();
    
    Ok(results)
}

/// Python-exposed function to calculate FSC
#[pyfunction]
#[pyo3(signature = (bedgz_path, bin_path))]
pub fn count_fragments_by_bins(
    py: Python<'_>,
    bedgz_path: &str,
    bin_path: &str,
) -> PyResult<(
    Py<PyArray1<u32>>,  // short counts
    Py<PyArray1<u32>>,  // intermediate counts
    Py<PyArray1<u32>>,  // long counts
    Py<PyArray1<u32>>,  // total counts
    Py<PyArray1<f64>>,  // mean GC
)> {
    let bedgz = Path::new(bedgz_path);
    let bins = Path::new(bin_path);
    
    let regions = parse_bin_file(bins)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    
    let results = count_fragments_sequential(bedgz, &regions)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    
    let shorts: Vec<u32> = results.iter().map(|r| r.short_count).collect();
    let intermediates: Vec<u32> = results.iter().map(|r| r.intermediate_count).collect();
    let longs: Vec<u32> = results.iter().map(|r| r.long_count).collect();
    let totals: Vec<u32> = results.iter().map(|r| r.total_count).collect();
    let gcs: Vec<f64> = results.iter().map(|r| r.mean_gc).collect();
    
    Ok((
        shorts.into_pyarray(py).into(),
        intermediates.into_pyarray(py).into(),
        longs.into_pyarray(py).into(),
        totals.into_pyarray(py).into(),
        gcs.into_pyarray(py).into(),
    ))
}
