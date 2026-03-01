"""
Motif processor for End Motif (EDM), Breakpoint Motif (BPM), and Motif Diversity Score (MDS).

Shared processing logic for writing motif output files.
Used by both standalone motif.py and run-all wrapper.py.
"""

from pathlib import Path
from typing import Dict
import numpy as np
import itertools
import logging

logger = logging.getLogger("core.motif_processor")


def _write_motif_helper(
    counts: Dict[str, int],
    output_path: Path,
    motif_type: str,
    kmer: int = 4,
    include_header: bool = True,
) -> int:
    """
    Write motif frequencies to TSV file (shared helper).

    Only valid ACGT k-mers are included (256 for k=4).

    Args:
        counts: Dictionary of k-mer -> count
        output_path: Path to write output TSV
        motif_type: Type name for logging (e.g., "End Motif", "Breakpoint Motif")
        kmer: K-mer size for initializing all possible k-mers
        include_header: Whether to include header line

    Returns:
        Total fragment count
    """
    bases = ["A", "C", "T", "G"]
    valid_kmers = set("".join(i) for i in itertools.product(bases, repeat=kmer))

    # Initialize all valid k-mers to 0, then update with counts
    all_kmers = {k: 0 for k in sorted(valid_kmers)}
    for k, v in counts.items():
        if k in valid_kmers:
            all_kmers[k] = v

    total = sum(all_kmers.values())
    logger.info(
        f"Writing {motif_type}: {output_path} ({len(all_kmers)} k-mers, {total:,} total)"
    )

    with open(output_path, "w") as f:
        if include_header:
            f.write("Motif\tFrequency\n")
        for k, v in all_kmers.items():
            freq = v / total if total else 0
            f.write(f"{k}\t{freq:.6f}\n")

    return total


def write_end_motif(
    em_counts: Dict[str, int],
    output_path: Path,
    kmer: int = 4,
    include_header: bool = True,
) -> int:
    """Write End Motif frequencies to TSV file."""
    return _write_motif_helper(
        em_counts, output_path, "End Motif", kmer, include_header
    )


def write_breakpoint_motif(
    bpm_counts: Dict[str, int],
    output_path: Path,
    kmer: int = 4,
    include_header: bool = True,
) -> int:
    """Write Breakpoint Motif frequencies to TSV file."""
    return _write_motif_helper(
        bpm_counts, output_path, "Breakpoint Motif", kmer, include_header
    )


def compute_mds(em_counts: Dict[str, int], kmer: int = 4) -> float:
    """
    Compute Motif Diversity Score from k-mer counts.

    MDS is the normalized Shannon entropy of the k-mer distribution,
    measuring diversity of fragment end sequences. Higher values indicate
    more diverse (healthy) samples; lower values may indicate tumor-derived
    cell-free DNA with altered fragmentation patterns.

    Only valid ACGT k-mers are used (256 for k=4).

    Args:
        em_counts: Dictionary of k-mer -> count
        kmer: K-mer size (default 4)

    Returns:
        MDS value normalized to [0.0, 1.0]
    """
    bases = ["A", "C", "T", "G"]
    valid_kmers = set("".join(i) for i in itertools.product(bases, repeat=kmer))

    # Filter to valid ACGT k-mers only
    all_kmers = {k: 0 for k in valid_kmers}
    for k, v in em_counts.items():
        if k in valid_kmers:
            all_kmers[k] = v

    total = sum(all_kmers.values())
    n_kmers = len(all_kmers)  # 256 for k=4

    if total == 0:
        return 0.0

    freq = np.array(list(all_kmers.values())) / total

    # MDS = -sum(p * log2(p)) / log2(N) where N = 4^k
    entropy = -np.sum(freq * np.log2(freq + 1e-12))
    mds = entropy / np.log2(n_kmers)

    return mds


def write_mds(
    em_counts: Dict[str, int],
    output_path: Path,
    sample_name: str = None,
    kmer: int = 4,
    include_header: bool = True,
) -> float:
    """
    Calculate and write Motif Diversity Score to TSV file.

    MDS is calculated as normalized Shannon entropy of end motif frequencies.
    Only valid ACGT k-mers are used (256 for k=4).

    Args:
        em_counts: Dictionary of k-mer -> count (End Motif)
        output_path: Path to write output TSV
        sample_name: Optional sample name for output
        kmer: K-mer size
        include_header: Whether to include header line

    Returns:
        MDS value (0.0 to 1.0, normalized by log2(4^k))
    """
    # Use the canonical compute_mds function
    mds = compute_mds(em_counts, kmer)

    logger.info(f"Writing MDS: {output_path} (MDS={mds:.6f})")

    with open(output_path, "w") as f:
        if include_header:
            f.write("Sample\tMDS\n")
        if sample_name:
            f.write(f"{sample_name}\t{mds:.6f}\n")
        else:
            f.write(f"{mds:.6f}\n")

    return mds


def process_motif_outputs(
    em_counts: Dict[str, int],
    bpm_counts: Dict[str, int],
    edm_output: Path,
    bpm_output: Path,
    mds_output: Path,
    sample_name: str = None,
    kmer: int = 4,
    include_headers: bool = True,
) -> tuple:
    """
    Write all motif outputs (EDM, BPM, MDS).

    This is the main entry point for both standalone and run-all.

    Args:
        em_counts: End motif k-mer counts
        bpm_counts: Breakpoint motif k-mer counts
        edm_output: Path for End Motif TSV
        bpm_output: Path for Breakpoint Motif TSV
        mds_output: Path for MDS TSV
        sample_name: Optional sample name for MDS output
        kmer: K-mer size (default 4)
        include_headers: Include headers in all output files (default True)

    Returns:
        Tuple of (total_em, total_bpm, mds_value)
    """
    total_em = write_end_motif(em_counts, edm_output, kmer, include_headers)
    total_bpm = write_breakpoint_motif(bpm_counts, bpm_output, kmer, include_headers)
    mds = write_mds(em_counts, mds_output, sample_name, kmer, include_headers)

    return total_em, total_bpm, mds


# =============================================================================
# 1-MER MOTIF ANALYSIS (Jagged Index / C-end Fraction)
# =============================================================================


def compute_1mer_from_kmer(em_counts: Dict[str, int]) -> Dict[str, int]:
    """
    Compute 1-mer (single base) counts from k-mer counts.

    Aggregates k-mer counts by their first base, which represents
    the fragment end base. This avoids needing Rust changes.

    Args:
        em_counts: Dictionary of k-mer -> count (e.g., 4-mers)

    Returns:
        Dictionary of 1-mer (A/C/G/T) -> count

    Example:
        >>> em_counts = {'ACGT': 100, 'ATTT': 50, 'CGAT': 75, 'GGGG': 25}
        >>> compute_1mer_from_kmer(em_counts)
        {'A': 150, 'C': 75, 'G': 25, 'T': 0}
    """
    counts_1mer: Dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}

    for kmer, count in em_counts.items():
        if len(kmer) >= 1:
            first_base = kmer[0].upper()
            if first_base in counts_1mer:
                counts_1mer[first_base] += count

    return counts_1mer


def compute_c_end_fraction(em_counts: Dict[str, int]) -> Dict[str, float]:
    """
    Compute C-end fraction and related metrics for jagged end analysis.

    The C-end fraction is elevated in healthy cfDNA due to apoptotic
    nuclease preferences. Lower C-fraction may indicate tumor-derived DNA
    with altered fragmentation patterns.

    Reference: Cristiano et al., Nature 2019 - DELFI approach

    Args:
        em_counts: Dictionary of k-mer -> count (4-mers from End Motif)

    Returns:
        Dictionary with:
        - 'a_fraction': A-end fraction
        - 'c_fraction': C-end fraction (key signal for jagged ends)
        - 'g_fraction': G-end fraction
        - 't_fraction': T-end fraction
        - 'entropy': Shannon entropy of 1-mer distribution
        - 'c_bias': C-fraction relative to expected 0.25

    Note:
        Healthy cfDNA typically has C-fraction of 0.28-0.32.
        Tumor-derived DNA may show reduced C-fraction.
    """
    counts_1mer = compute_1mer_from_kmer(em_counts)

    total = sum(counts_1mer.values())
    if total == 0:
        return {
            "a_fraction": 0.0,
            "c_fraction": 0.0,
            "g_fraction": 0.0,
            "t_fraction": 0.0,
            "entropy": 0.0,
            "c_bias": 0.0,
        }

    fractions = {
        "a_fraction": counts_1mer["A"] / total,
        "c_fraction": counts_1mer["C"] / total,
        "g_fraction": counts_1mer["G"] / total,
        "t_fraction": counts_1mer["T"] / total,
    }

    # Shannon entropy of 1-mer distribution
    freqs = np.array([counts_1mer[b] for b in "ACGT"]) / total
    entropy = -np.sum(freqs * np.log2(freqs + 1e-12))
    fractions["entropy"] = entropy

    # C-bias relative to expected 0.25
    fractions["c_bias"] = fractions["c_fraction"] - 0.25

    return fractions


def write_end_motif_1mer(
    em_counts: Dict[str, int],
    output_path: Path,
    sample_name: str = None,
    include_header: bool = True,
) -> Dict[str, float]:
    """
    Write 1-mer end motif frequencies and C-end metrics to TSV file.

    Creates {sample}.EndMotif1mer.tsv with single base frequencies
    and derived jagged end metrics.

    Args:
        em_counts: Dictionary of k-mer -> count (4-mers)
        output_path: Path to write output TSV
        sample_name: Optional sample name
        include_header: Whether to include header line

    Returns:
        Dictionary of C-end fraction and related metrics

    Output columns:
        base, count, fraction
        + summary row with c_fraction, entropy, c_bias
    """
    counts_1mer = compute_1mer_from_kmer(em_counts)
    metrics = compute_c_end_fraction(em_counts)

    total = sum(counts_1mer.values())
    logger.info(f"Writing 1-mer End Motif: {output_path} ({total:,} fragments)")

    with open(output_path, "w") as f:
        if include_header:
            f.write("base\tcount\tfraction\n")

        for base in "ACGT":
            count = counts_1mer[base]
            frac = count / total if total > 0 else 0.0
            f.write(f"{base}\t{count}\t{frac:.6f}\n")

        # Summary metrics
        f.write(f"# c_fraction\t{metrics['c_fraction']:.6f}\n")
        f.write(f"# entropy\t{metrics['entropy']:.6f}\n")
        f.write(f"# c_bias\t{metrics['c_bias']:.6f}\n")
        if sample_name:
            f.write(f"# sample\t{sample_name}\n")

    logger.debug(
        f"  C-end fraction: {metrics['c_fraction']:.4f} (bias: {metrics['c_bias']:+.4f})"
    )

    return metrics
