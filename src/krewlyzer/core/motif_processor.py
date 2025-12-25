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


def write_end_motif(
    em_counts: Dict[str, int],
    output_path: Path,
    kmer: int = 4,
    include_header: bool = True
) -> int:
    """
    Write End Motif frequencies to TSV file.
    
    Args:
        em_counts: Dictionary of k-mer -> count
        output_path: Path to write output TSV
        kmer: K-mer size for initializing all possible k-mers
        include_header: Whether to include header line (default True for consistency)
        
    Returns:
        Total fragment count
    """
    bases = ['A', 'C', 'T', 'G']
    all_kmers = {''.join(i): 0 for i in itertools.product(bases, repeat=kmer)}
    all_kmers.update(em_counts)
    
    total = sum(all_kmers.values())
    logger.info(f"Writing End Motif: {output_path} ({len(all_kmers)} k-mers, {total:,} total)")
    
    with open(output_path, 'w') as f:
        if include_header:
            f.write("Motif\tFrequency\n")
        for k, v in all_kmers.items():
            freq = v / total if total else 0
            f.write(f"{k}\t{freq:.6f}\n")
    
    return total


def write_breakpoint_motif(
    bpm_counts: Dict[str, int],
    output_path: Path,
    kmer: int = 4,
    include_header: bool = True
) -> int:
    """
    Write Breakpoint Motif frequencies to TSV file.
    
    Args:
        bpm_counts: Dictionary of k-mer -> count
        output_path: Path to write output TSV
        kmer: K-mer size for initializing all possible k-mers
        include_header: Whether to include header line (default True for consistency)
        
    Returns:
        Total fragment count
    """
    bases = ['A', 'C', 'T', 'G']
    all_kmers = {''.join(i): 0 for i in itertools.product(bases, repeat=kmer)}
    all_kmers.update(bpm_counts)
    
    total = sum(all_kmers.values())
    logger.info(f"Writing Breakpoint Motif: {output_path} ({len(all_kmers)} k-mers, {total:,} total)")
    
    with open(output_path, 'w') as f:
        if include_header:
            f.write("Motif\tFrequency\n")
        for k, v in all_kmers.items():
            freq = v / total if total else 0
            f.write(f"{k}\t{freq:.6f}\n")
    
    return total


def write_mds(
    em_counts: Dict[str, int],
    output_path: Path,
    sample_name: str = None,
    kmer: int = 4,
    include_header: bool = True
) -> float:
    """
    Calculate and write Motif Diversity Score to TSV file.
    
    MDS is calculated as normalized Shannon entropy of end motif frequencies.
    
    Args:
        em_counts: Dictionary of k-mer -> count (End Motif)
        output_path: Path to write output TSV
        sample_name: Optional sample name for output
        kmer: K-mer size
        include_header: Whether to include header line
        
    Returns:
        MDS value
    """
    bases = ['A', 'C', 'T', 'G']
    all_kmers = {''.join(i): 0 for i in itertools.product(bases, repeat=kmer)}
    all_kmers.update(em_counts)
    
    total = sum(all_kmers.values())
    
    # Calculate normalized Shannon entropy
    if total > 0:
        freq = np.array(list(all_kmers.values())) / total
    else:
        freq = np.zeros(len(all_kmers))
    
    # MDS = -sum(p * log2(p)) / log2(N)
    # Add small epsilon to avoid log(0)
    mds = -np.sum(freq * np.log2(freq + 1e-12)) / np.log2(len(freq))
    
    logger.info(f"Writing MDS: {output_path} (MDS={mds:.6f})")
    
    with open(output_path, 'w') as f:
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
    include_headers: bool = True
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
