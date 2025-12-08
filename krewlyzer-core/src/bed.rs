//! BED file parsing and fragment counting
//!
//! Handles tabix-indexed BED.gz files for fragment analysis.

use std::path::Path;
use rayon::prelude::*;

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

// TODO: Phase 1 implementation
// - parse_bed_regions() - load regions from BIN file
// - count_fragments_in_regions() - parallel counting using rayon
// - calc_fsc() - fragment size coverage with GC correction
// - calc_fsr() - fragment size ratios
