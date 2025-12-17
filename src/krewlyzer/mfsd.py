"""
Mutant Fragment Size Distribution (mFSD) calculation.

Calculates mFSD features for a single sample.
Uses Rust backend for accelerated computation.

Supports all small variant types: SNV, MNV, Insertion, Deletion, Complex.
4-way fragment classification: REF, ALT, NonREF, N.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("mfsd")

# Rust backend is required
from krewlyzer import _core


def mfsd(
    bam_input: Path = typer.Argument(..., help="Input BAM file"),
    input_file: Path = typer.Option(..., "--input-file", "-i", help="VCF or MAF file with variants"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file (default: derived from input filename)"),
    mapq: int = typer.Option(20, "--mapq", "-q", help="Minimum mapping quality"),
    output_distributions: bool = typer.Option(False, "--output-distributions", "-d", help="Output per-variant size distributions"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose/debug logging"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)")
):
    """
    Calculate Mutant Fragment Size Distribution (mFSD) features for a single sample.
    
    Supports all small variant types: SNV, MNV, Insertion, Deletion, Complex.
    
    Fragment classification:
    - REF: Supports reference allele (healthy cfDNA baseline)
    - ALT: Supports alternate allele (tumor signal)
    - NonREF: Non-REF, non-ALT (sequencing errors, subclones)
    - N: Contains N at variant position (low quality)
    
    Input: BAM file and VCF/MAF file with variants
    Output: 
    - {sample}.mFSD.tsv: Summary statistics (39 columns)
    - {sample}.mFSD.distributions.tsv: Per-size fragment counts (optional, with -d)
    """
    # Set log level based on verbose flag
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    # Configure Rust thread pool
    if threads > 0:
        try:
            _core.configure_threads(threads)
            logger.info(f"Configured {threads} threads for parallel processing")
        except Exception as e:
            logger.warning(f"Could not configure threads: {e}")
    
    # Input validation
    if not bam_input.exists():
        logger.error(f"BAM file not found: {bam_input}")
        raise typer.Exit(1)
    
    if not str(bam_input).endswith('.bam'):
        logger.error(f"Input must be a .bam file: {bam_input}")
        raise typer.Exit(1)
    
    if not input_file.exists():
        logger.error(f"Variant file not found: {input_file}")
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bam_input.stem.replace('.bam', '')
    
    output_file = output / f"{sample_name}.mFSD.tsv"
    
    try:
        logger.info(f"Processing {bam_input.name}")
        
        # Detect input type
        input_type = "vcf" if str(input_file).endswith(".vcf") or str(input_file).endswith(".vcf.gz") else "maf"
        logger.info(f"Detected input type: {input_type}")
        
        # Call Rust backend (new signature with output_distributions)
        _core.mfsd.calculate_mfsd(
            str(bam_input),
            str(input_file),
            str(output_file),
            input_type,
            mapq,
            output_distributions
        )
        
        logger.info(f"mFSD complete: {output_file}")
        if output_distributions:
            dist_file = output_file.with_suffix('.distributions.tsv')
            logger.info(f"Distributions: {dist_file}")

    except Exception as e:
        logger.error(f"mFSD calculation failed: {e}")
        raise typer.Exit(1)
