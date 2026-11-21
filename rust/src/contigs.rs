//! One rule for matching a BAM contig to the reference FASTA.
//!
//! Every place a BAM contig meets the reference goes through [`ContigNames`]: extract's GC and
//! breakpoint-context fetches, mFSD (variant to BAM, and BAM to FASTA) and the expected-GC
//! builder. Before it, each site had its own rule and none checked lengths. Extract matched
//! exactly, so a `chr1`-named BAM read with a `1`-named FASTA silently gave every fragment GC 0.0;
//! the others tried the second spelling after a miss, which crashed instead, because rust-htslib's
//! `fetch_seq` answers a name the FASTA lacks with a `Vec` over a null pointer, not `Err`. Fetching
//! only under a name taken from the FASTA's own index (here) is what keeps that from happening.
//!
//! Resolution, per BAM contig (name and `@SQ` length):
//! 1. the exact name, if the FASTA has it;
//! 2. otherwise the name without a leading `chr`, matched against FASTA names without one;
//! 3. the FASTA length must equal the BAM length -- a match by name with another length is a
//!    different sequence, which no name rule can detect. On a primary contig (1-22, X, Y) that
//!    means another assembly (hg38 against hg19) and is refused; on any other contig (a chrM
//!    of 16,571 bp against rCRS at 16,569) the contig is left unresolved and its fragments get no GC.
//!
//! Exact-first means a FASTA holding both `chr1` and `1` never makes a BAM contig ambiguous: either
//! spelling matches one of them exactly. `chrM` and `MT` are not special-cased: UCSC hg19 `chrM`
//! (16,571 bp) and b37 `MT` (16,569 bp) are different sequences, and the length check keeps them
//! apart.

use std::collections::HashMap;

/// A contig name without a leading `chr`, so `chr1` and `1` compare equal.
pub fn norm_chrom(c: &str) -> &str {
    c.strip_prefix("chr").unwrap_or(c)
}

/// True for 1-22, X and Y in either spelling: the contigs a reference must have for the sample.
pub fn is_primary(name: &str) -> bool {
    let k = norm_chrom(name);
    k == "X" || k == "Y" || k.parse::<u8>().is_ok_and(|n| (1..=22).contains(&n))
}

/// How one BAM contig resolves against the FASTA.
#[derive(Debug, PartialEq, Eq)]
pub enum Resolution {
    /// The FASTA's name for it (the same spelling or the other `chr` one), with the same length.
    Resolved(String),
    /// No FASTA contig under either spelling.
    Absent,
    /// A FASTA contig matches by name but has another length: a different assembly.
    LengthMismatch { fasta: String, bam_len: u64, fasta_len: u64 },
}

impl Resolution {
    /// The matching contig's own name (the FASTA's spelling, to fetch with), if it resolved.
    pub fn name(&self) -> Option<&str> {
        match self {
            Resolution::Resolved(n) => Some(n),
            _ => None,
        }
    }
}

/// The reference FASTA's contigs (name -> length), indexed for [`Resolution`].
#[derive(Debug, Default)]
pub struct ContigNames {
    lengths: HashMap<String, u64>,
    by_key: HashMap<String, String>,
}

impl ContigNames {
    pub fn new<I: IntoIterator<Item = (String, u64)>>(fasta_contigs: I) -> Self {
        let mut names = ContigNames::default();
        for (name, len) in fasta_contigs {
            // With both spellings present, the fallback key keeps the first; exact-first
            // resolution means it is only consulted for a name the FASTA does not hold exactly.
            names.by_key.entry(norm_chrom(&name).to_string()).or_insert_with(|| name.clone());
            names.lengths.insert(name, len);
        }
        names
    }

    /// From an open FASTA index: every sequence's name and length.
    pub fn from_faidx(fa: &rust_htslib::faidx::Reader) -> Self {
        let n = fa.n_seqs();
        let mut contigs = Vec::with_capacity(n as usize);
        for i in 0..n {
            if let Ok(name) = fa.seq_name(i as i32) {
                let len = fa.fetch_seq_len(&name);
                contigs.push((name, len));
            }
        }
        Self::new(contigs)
    }

    /// The length of a contig, by the FASTA's own name for it (as [`Resolution::name`] gives it).
    pub fn length(&self, fasta_name: &str) -> Option<u64> {
        self.lengths.get(fasta_name).copied()
    }

    /// Resolve one BAM contig. `bam_len` of 0 means unknown and skips the length check.
    pub fn resolve(&self, bam_name: &str, bam_len: u64) -> Resolution {
        let fasta = if self.lengths.contains_key(bam_name) {
            bam_name
        } else if let Some(other) = self.by_key.get(norm_chrom(bam_name)) {
            other.as_str()
        } else {
            return Resolution::Absent;
        };
        let fasta_len = self.lengths[fasta];
        if bam_len != 0 && fasta_len != bam_len {
            return Resolution::LengthMismatch { fasta: fasta.to_string(), bam_len, fasta_len };
        }
        Resolution::Resolved(fasta.to_string())
    }

    /// Resolve every BAM contig and refuse a reference that cannot serve the sample: a primary
    /// contig (1-22, X, Y) absent from the FASTA or present with another length. Other contigs
    /// that do not resolve (decoys, alts, unplaced; a non-primary length mismatch, warned about
    /// here) are `None`, and their fragments are counted by the caller. Returns the FASTA name per
    /// BAM contig, in BAM order.
    pub fn resolve_header(&self, bam: &[(String, u64)]) -> Result<Vec<Option<String>>, String> {
        let mut out = Vec::with_capacity(bam.len());
        for (name, len) in bam {
            match self.resolve(name, *len) {
                Resolution::Resolved(n) => out.push(Some(n)),
                Resolution::Absent if is_primary(name) => {
                    return Err(format!(
                        "contig {name} is in the BAM but not in the reference FASTA under either \
                         spelling; this reference cannot serve the sample"
                    ))
                }
                Resolution::Absent => out.push(None),
                Resolution::LengthMismatch { fasta, bam_len, fasta_len } if !is_primary(name) => {
                    log::warn!("contig {name} is {bam_len} bp in the BAM but {fasta} is {fasta_len} bp \
                                in the reference FASTA; its fragments get no GC");
                    out.push(None)
                }
                Resolution::LengthMismatch { fasta, bam_len, fasta_len } => {
                    return Err(format!(
                        "contig {name} is {bam_len} bp in the BAM but {fasta} is {fasta_len} bp in \
                         the reference FASTA: a different assembly"
                    ))
                }
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b37() -> ContigNames {
        ContigNames::new([
            ("1".to_string(), 249_250_621),
            ("X".to_string(), 155_270_560),
            ("MT".to_string(), 16_569),
            ("hs37d5".to_string(), 35_477_943),
        ])
    }

    #[test]
    fn exact_names_resolve_as_themselves() {
        assert_eq!(b37().resolve("1", 249_250_621), Resolution::Resolved("1".into()));
    }

    #[test]
    fn the_other_chr_spelling_resolves_in_both_directions() {
        assert_eq!(b37().resolve("chr1", 249_250_621), Resolution::Resolved("1".into()));
        let ucsc = ContigNames::new([("chr1".to_string(), 249_250_621)]);
        assert_eq!(ucsc.resolve("1", 249_250_621), Resolution::Resolved("chr1".into()));
    }

    #[test]
    fn a_name_match_with_another_length_is_a_different_assembly() {
        // hg38 chr1 is 248,956,422 bp; b37 1 is 249,250,621 bp.
        let r = b37().resolve("chr1", 248_956_422);
        assert!(matches!(r, Resolution::LengthMismatch { .. }), "{r:?}");
        assert!(b37().resolve_header(&[("chr1".into(), 248_956_422)]).is_err());
        // The same spelling with another length too: hg38 chr1 against hg19 chr1.
        let hg19 = ContigNames::new([("chr1".to_string(), 249_250_621)]);
        assert!(hg19.resolve_header(&[("chr1".into(), 248_956_422)]).is_err());
    }

    #[test]
    fn chrm_and_mt_stay_apart_by_length() {
        // UCSC hg19 chrM is 16,571 bp; b37 MT is 16,569 bp: "MT" never matches "chrM" by name, and
        // a b37 BAM's MT against a hg19 FASTA's chrM is Absent, not a silent rename.
        assert_eq!(b37().resolve("chrM", 16_571), Resolution::Absent);
        assert_eq!(b37().resolve_header(&[("chrM".into(), 16_571)]).unwrap(), vec![None]);
    }

    #[test]
    fn another_length_on_a_non_primary_contig_is_unresolved_not_refused() {
        // UCSC hg19 chrM (16,571 bp) against an hg19 build carrying rCRS chrM (16,569 bp) under
        // the same name: its fragments get no GC, the primary contigs still serve the sample.
        let rcrs = ContigNames::new([("chr1".to_string(), 249_250_621), ("chrM".to_string(), 16_569)]);
        let bam = [("chr1".to_string(), 249_250_621), ("chrM".to_string(), 16_571)];
        assert_eq!(rcrs.resolve_header(&bam).unwrap(), vec![Some("chr1".to_string()), None]);
    }

    #[test]
    fn exact_wins_when_the_fasta_holds_both_spellings() {
        let both = ContigNames::new([("1".to_string(), 10), ("chr1".to_string(), 20)]);
        assert_eq!(both.resolve("1", 10), Resolution::Resolved("1".into()));
        assert_eq!(both.resolve("chr1", 20), Resolution::Resolved("chr1".into()));
    }

    #[test]
    fn a_missing_primary_contig_is_refused_and_a_missing_decoy_is_not() {
        let bam = [("1".to_string(), 249_250_621), ("GL000192.1".to_string(), 547_496)];
        assert_eq!(b37().resolve_header(&bam).unwrap(), vec![Some("1".to_string()), None]);
        assert!(b37().resolve_header(&[("22".into(), 51_304_566)]).is_err());
    }

    #[test]
    fn primary_contigs_are_1_to_22_x_and_y_in_either_spelling() {
        for n in ["1", "chr22", "X", "chrY"] {
            assert!(is_primary(n), "{n}");
        }
        for n in ["23", "MT", "chrM", "hs37d5", "chrUn_gl000220", "1_random"] {
            assert!(!is_primary(n), "{n}");
        }
    }
}
