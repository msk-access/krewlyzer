//! Read filtering logic
//!
//! Handles BAM/SAM read filtering based on flags and quality thresholds.

use crate::ReadFilters;

/// Check if a read passes all filters
pub fn passes_filters(
    flags: u16,
    mapq: u8,
    fragment_length: u32,
    filters: &ReadFilters,
) -> bool {
    // SAM flags (bit positions)
    const FLAG_UNMAPPED: u16 = 0x4;
    const FLAG_SECONDARY: u16 = 0x100;
    const FLAG_FAILED_QC: u16 = 0x200;
    const FLAG_DUPLICATE: u16 = 0x400;
    const FLAG_SUPPLEMENTARY: u16 = 0x800;
    const FLAG_PROPER_PAIR: u16 = 0x2;
    
    // Skip unmapped reads
    if flags & FLAG_UNMAPPED != 0 {
        return false;
    }
    
    // Skip secondary alignments
    if flags & FLAG_SECONDARY != 0 {
        return false;
    }
    
    // Skip supplementary alignments
    if flags & FLAG_SUPPLEMENTARY != 0 {
        return false;
    }
    
    // Skip failed QC
    if flags & FLAG_FAILED_QC != 0 {
        return false;
    }
    
    // Skip duplicates if configured
    if filters.skip_duplicates && flags & FLAG_DUPLICATE != 0 {
        return false;
    }
    
    // Require proper pair if configured
    if filters.require_proper_pair && flags & FLAG_PROPER_PAIR == 0 {
        return false;
    }
    
    // Check mapping quality
    if mapq < filters.mapq {
        return false;
    }
    
    // Check fragment length
    if fragment_length < filters.min_length || fragment_length > filters.max_length {
        return false;
    }
    
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_passes_filters_basic() {
        let filters = ReadFilters {
            mapq: 20,
            min_length: 65,
            max_length: 400,
            skip_duplicates: true,
            require_proper_pair: true,
        };
        
        // Good read: proper pair, mapq=30, length=150
        assert!(passes_filters(0x2, 30, 150, &filters));
        
        // Bad: unmapped
        assert!(!passes_filters(0x4, 30, 150, &filters));
        
        // Bad: duplicate
        assert!(!passes_filters(0x402, 30, 150, &filters));
        
        // Bad: low mapq
        assert!(!passes_filters(0x2, 10, 150, &filters));
        
        // Bad: too short
        assert!(!passes_filters(0x2, 30, 50, &filters));
    }
}
