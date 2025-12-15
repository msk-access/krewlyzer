"""
Motif feature extraction from BAM files.

Extracts End Motif (EDM), Breakpoint Motif (BPM), and Motif Diversity Score (MDS).
Uses Rust backend for accelerated extraction.

Note: For fragment extraction (BED.gz), use `krewlyzer extract` instead.
"""

import typer
from pathlib import Path
from typing import Optional
import logging
import numpy as np
import itertools

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("motif")

# Rust backend is required
from krewlyzer import _core


def motif(
    bam_input: Path = typer.Argument(..., help="Input BAM file"),
    genome_reference: Path = typer.Option(..., '-g', '--genome', help="Reference genome FASTA (indexed)"),
    output: Path = typer.Option(..., '-o', '--output', help="Output directory"),
    kmer: int = typer.Option(4, '-k', '--kmer', help="K-mer size for motif extraction"),
    chromosomes: Optional[str] = typer.Option(None, '--chromosomes', help="Comma-separated chromosomes to process"),
    sample_name: Optional[str] = typer.Option(None, '--sample-name', '-s', help="Sample name for output files (default: derived from BAM filename)"),
    threads: int = typer.Option(0, '--threads', '-t', help="Number of threads (0=all cores)")
):
    """
    Extract k-mer motif features from a BAM file.
    
    Output:
    - {sample}.EndMotif: End motif k-mer frequencies
    - {sample}.BreakPointMotif: Breakpoint motif k-mer frequencies
    - {sample}.MDS: Motif Diversity Score
    
    Note: For fragment extraction (BED.gz), use `krewlyzer extract` instead.
    """
    # Configure Rust thread pool
    if threads > 0:
        try:
            _core.configure_threads(threads)
            logger.info(f"Configured {threads} threads")
        except Exception as e:
            logger.warning(f"Could not configure threads: {e}")
    
    # Input validation
    if not bam_input.exists():
        logger.error(f"BAM file not found: {bam_input}")
        raise typer.Exit(1)
    
    if not str(bam_input).endswith('.bam'):
        logger.error(f"Input must be a .bam file: {bam_input}")
        raise typer.Exit(1)
    
    if not genome_reference.exists():
        logger.error(f"Reference genome not found: {genome_reference}")
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from BAM filename)
    if sample_name is None:
        sample_name = bam_input.stem.replace('.bam', '')
    
    # Output file paths
    edm_output = output / f"{sample_name}.EndMotif.tsv"
    bpm_output = output / f"{sample_name}.BreakPointMotif.tsv"
    mds_output = output / f"{sample_name}.MDS.tsv"
    
    try:
        logger.info(f"Extracting motif features from {bam_input.name}")
        
        # Initialize motif dictionaries
        bases = ['A', 'C', 'T', 'G']
        End_motif = {''.join(i): 0 for i in itertools.product(bases, repeat=kmer)}
        Breakpoint_motif = {''.join(i): 0 for i in itertools.product(bases, repeat=kmer)}
        
        # Parse chromosomes
        chroms = chromosomes.split(',') if chromosomes else None
        
        # Call Unified Rust Engine (Extract + Motif)
        # We pass output_motif_prefix="enable" to trigger motif counting in Rust
        # We pass output_bed_path=None to skip writing BED (unless debugging)
        
        fragment_count, em_counts, bpm_counts = _core.extract_motif.process_bam_parallel(
            str(bam_input),
            str(genome_reference),
            20,    # Default mapQ
            65,    # Default min length
            400,   # Default max length
            kmer,
            threads,
            None,  # output_bed_path
            "enable" # output_motif_prefix (triggers counting)
        )
        
        End_motif.update(em_counts)
        Breakpoint_motif.update(bpm_counts)
        logger.info(f"Processed {fragment_count:,} fragments")
        
        # Write End Motif
        logger.info(f"Writing End Motif: {edm_output}")
        total_em = sum(End_motif.values())
        with open(edm_output, 'w') as f:
            for k, v in End_motif.items():
                f.write(f"{k}\t{v/total_em if total_em else 0}\n")
        
        # Write Breakpoint Motif
        logger.info(f"Writing Breakpoint Motif: {bpm_output}")
        total_bpm = sum(Breakpoint_motif.values())
        with open(bpm_output, 'w') as f:
            for k, v in Breakpoint_motif.items():
                f.write(f"{k}\t{v/total_bpm if total_bpm else 0}\n")
        
        # Write MDS
        logger.info(f"Writing MDS: {mds_output}")
        freq = np.array(list(End_motif.values())) / total_em if total_em else np.zeros(len(End_motif))
        mds = -np.sum(freq * np.log2(freq + 1e-12)) / np.log2(len(freq))
        with open(mds_output, 'w') as f:
            f.write(f"{mds}\n")
        
        logger.info(f"Motif extraction complete")

    except Exception as e:
        logger.error(f"Motif extraction failed: {e}")
        raise typer.Exit(1)
