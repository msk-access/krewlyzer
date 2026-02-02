//! BED file parsing and fragment counting
//!
//! Handles both standard gzip and BGZF-compressed BED files for fragment analysis.
//! Uses noodles::bgzf for BGZF files (optimal for tabix-indexed files) with
//! flate2::MultiGzDecoder as fallback for standard gzip.

use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use anyhow::{Context, Result};

/// Check if a file appears to be BGZF format by examining the header.
/// 
/// BGZF files have the gzip magic bytes (0x1f 0x8b) followed by specific
/// flags indicating extra fields, and a "BC" subfield identifier.
fn is_bgzf_file(path: &Path) -> bool {
    if let Ok(mut file) = File::open(path) {
        let mut header = [0u8; 18];
        if file.read_exact(&mut header).is_ok() {
            // Check gzip magic (0x1f 0x8b)
            if header[0] != 0x1f || header[1] != 0x8b {
                return false;
            }
            // Check FEXTRA flag (bit 2 of FLG byte at position 3)
            if header[3] & 0x04 == 0 {
                return false;
            }
            // Check for BC subfield at position 12-13
            if header[12] == b'B' && header[13] == b'C' {
                return true;
            }
        }
    }
    false
}

/// Open a file and return a buffered reader, handling compression transparently.
/// 
/// For .gz files:
/// - First checks if the file is BGZF format (via magic bytes)
/// - Uses noodles::bgzf::io::Reader for BGZF files (optimal for tabix-indexed)
/// - Falls back to flate2::MultiGzDecoder for standard gzip files
/// 
/// This is critical for tabix-indexed BED.gz files which use BGZF compression.
pub fn get_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let is_gz = path.extension()
        .map(|ext| ext == "gz")
        .unwrap_or(false);
        
    if is_gz {
        if is_bgzf_file(path) {
            // Use noodles::bgzf for proper BGZF reading (supports seeking, multi-block)
            let file = File::open(path)
                .with_context(|| format!("Failed to open BGZF file: {:?}", path))?;
            let bgzf_reader = noodles::bgzf::io::Reader::new(file);
            Ok(Box::new(BufReader::new(bgzf_reader)))
        } else {
            // Fallback to MultiGzDecoder for standard gzip files
            let file = File::open(path)
                .with_context(|| format!("Failed to open gzip file: {:?}", path))?;
            use flate2::read::MultiGzDecoder;
            Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
        }
    } else {
        let file = File::open(path)
            .with_context(|| format!("Failed to open file: {:?}", path))?;
        Ok(Box::new(BufReader::new(file)))
    }
}
/// Count fragments in a region by size category
#[derive(Debug, Default, Clone)]
pub struct FragmentCounts {
    pub short: u32,        // 65-150bp
    pub intermediate: u32, // 151-260bp
    pub long: u32,         // 261-400bp
    pub total: u32,        // 65-400bp
}

impl FragmentCounts {
    pub fn new() -> Self {
        Self::default()
    }
    
    pub fn add_fragment(&mut self, length: u32) {
        if length >= 65 && length <= 400 {
            self.total += 1;
            
            if length <= 150 {
                self.short += 1;
            } else if length <= 260 {
                self.intermediate += 1;
            } else {
                self.long += 1;
            }
        }
    }
}

/// Region definition (from BIN file)
#[derive(Debug, Clone)]
pub struct Region {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl Region {
    pub fn new(chrom: String, start: u64, end: u64) -> Self {
        Self { chrom, start, end }
    }
}

/// Unified Fragment definition for the engine
#[derive(Debug, Clone, Copy)]
pub struct Fragment {
    pub chrom_id: u32,
    pub start: u64,
    pub end: u64,
    pub length: u64,
    pub gc: f64,
    /// GC correction weight (default 1.0)
    pub weight: f64,
}

/// Helper to map string chromosomes to integer IDs
#[derive(Debug, Default, Clone)]
pub struct ChromosomeMap {
    pub name_to_id: std::collections::HashMap<String, u32>,
    pub id_to_name: Vec<String>,
}

impl ChromosomeMap {
    pub fn new() -> Self {
        Self::default()
    }
    
    pub fn get_id(&mut self, chrom: &str) -> u32 {
        if let Some(&id) = self.name_to_id.get(chrom) {
            id
        } else {
            let id = self.id_to_name.len() as u32;
            self.id_to_name.push(chrom.to_string());
            self.name_to_id.insert(chrom.to_string(), id);
            id
        }
    }
    
    pub fn get_name(&self, id: u32) -> Option<&str> {
        if (id as usize) < self.id_to_name.len() {
            Some(&self.id_to_name[id as usize])
        } else {
            None
        }
    }
}

// Note: Fragment counting functions are implemented in:
//   - fsc.rs: calc_fsc() - fragment size coverage with GC correction
//   - fsr.rs: calc_fsr() - fragment size ratios
//   - fsd.rs: calculate_fsd() - fragment size distributions
