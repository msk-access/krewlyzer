//! Shared Motif Utilities for krewlyzer
//!
//! Contains common functions for DNA sequence manipulation and motif analysis.
//! Used by both `extract_motif.rs` (global MDS) and `region_mds.rs` (per-region MDS).

/// Reverse complement a DNA sequence.
///
/// Converts A<->T and G<->C, reversing the sequence.
/// Non-ACGT characters are passed through unchanged.
///
/// # Example
/// ```ignore
/// let seq = b"ACGT";
/// assert_eq!(reverse_complement(seq), vec![b'A', b'C', b'G', b'T']);
/// ```
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'a' => b'T',
            b'T' => b'A',
            b't' => b'A',
            b'G' => b'C',
            b'g' => b'C',
            b'C' => b'G',
            b'c' => b'G',
            x => *x,
        })
        .collect()
}

/// Convert a nucleotide to its 2-bit encoding index.
///
/// A/a -> 0, C/c -> 1, G/g -> 2, T/t -> 3.
/// Returns None for non-ACGT characters.
#[inline]
fn base_to_index(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Convert a 4-mer sequence to an index (0-255).
///
/// Uses 2-bit encoding: A=0, C=1, G=2, T=3.
/// Index = base0 * 64 + base1 * 16 + base2 * 4 + base3.
///
/// # Returns
/// Some(index) if all bases are valid ACGT, None otherwise.
pub fn kmer4_to_index(kmer: &[u8]) -> Option<usize> {
    if kmer.len() != 4 {
        return None;
    }
    let b0 = base_to_index(kmer[0])?;
    let b1 = base_to_index(kmer[1])?;
    let b2 = base_to_index(kmer[2])?;
    let b3 = base_to_index(kmer[3])?;
    Some(b0 * 64 + b1 * 16 + b2 * 4 + b3)
}

/// Calculate Motif Diversity Score (MDS) from a 4-mer histogram.
///
/// MDS is the Shannon entropy of the 4-mer distribution, normalized
/// by the maximum possible entropy (log2(256) = 8 bits).
///
/// Range: 0.0 (all fragments have same motif) to 1.0 (uniform distribution).
///
/// Formula: MDS = H / log2(256) where H = -Σ p(x) * log2(p(x))
///
/// # Arguments
/// * `counts` - Array of 256 counts, one per possible 4-mer
///
/// # Returns
/// MDS value in range [0.0, 1.0]
pub fn calculate_mds(counts: &[u64; 256]) -> f64 {
    let total: u64 = counts.iter().sum();
    if total == 0 {
        return 0.0;
    }

    let total_f = total as f64;
    let mut entropy = 0.0;

    for &count in counts {
        if count > 0 {
            let p = count as f64 / total_f;
            entropy -= p * p.log2();
        }
    }

    // Normalize by maximum entropy (log2(256) = 8)
    entropy / 8.0
}

/// Calculate GC content of a DNA sequence.
///
/// Returns the fraction of G+C bases among valid ACGT bases.
/// Non-ACGT characters are ignored.
pub fn calculate_gc(seq: &[u8]) -> f64 {
    let mut gc = 0;
    let mut valid = 0;
    for &b in seq {
        match b {
            b'G' | b'g' | b'C' | b'c' => {
                gc += 1;
                valid += 1;
            }
            b'A' | b'a' | b'T' | b't' => {
                valid += 1;
            }
            _ => {}
        }
    }
    if valid == 0 {
        0.0
    } else {
        gc as f64 / valid as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), vec![b'A', b'C', b'G', b'T']);
        assert_eq!(reverse_complement(b"AAAA"), vec![b'T', b'T', b'T', b'T']);
        assert_eq!(reverse_complement(b""), Vec::<u8>::new());
    }

    #[test]
    fn test_kmer4_to_index() {
        // AAAA = 0*64 + 0*16 + 0*4 + 0 = 0
        assert_eq!(kmer4_to_index(b"AAAA"), Some(0));
        // TTTT = 3*64 + 3*16 + 3*4 + 3 = 255
        assert_eq!(kmer4_to_index(b"TTTT"), Some(255));
        // ACGT = 0*64 + 1*16 + 2*4 + 3 = 27
        assert_eq!(kmer4_to_index(b"ACGT"), Some(27));
        // Invalid: contains N
        assert_eq!(kmer4_to_index(b"ACGN"), None);
        // Invalid: wrong length
        assert_eq!(kmer4_to_index(b"ACG"), None);
    }

    #[test]
    fn test_calculate_mds() {
        // All zeros -> 0.0
        let zeros = [0u64; 256];
        assert_eq!(calculate_mds(&zeros), 0.0);

        // Single motif -> 0.0 entropy
        let mut single = [0u64; 256];
        single[0] = 100;
        assert_eq!(calculate_mds(&single), 0.0);

        // Uniform distribution -> 1.0 (max entropy)
        let mut uniform = [1u64; 256];
        let mds = calculate_mds(&uniform);
        assert!((mds - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_calculate_gc() {
        assert_eq!(calculate_gc(b"GGCC"), 1.0);
        assert_eq!(calculate_gc(b"AATT"), 0.0);
        assert_eq!(calculate_gc(b"ACGT"), 0.5);
        assert_eq!(calculate_gc(b""), 0.0);
    }
}
