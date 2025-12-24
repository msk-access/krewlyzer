"""
Fragment-level Methylation analysis (UXM) calculation.

Calculates UXM features for a single sample.
Uses Rust backend for accelerated computation.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("uxm")

# Rust backend is required
from krewlyzer import _core


def uxm(
    bam_input: Path = typer.Argument(..., help="Input bisulfite BAM file"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file (default: derived from input filename)"),
    mark_input: Optional[Path] = typer.Option(None, "--mark-input", "-m", help="Path to genomic marker file"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging")
):
    """
    Calculate Fragment-level Methylation (UXM) features for a single sample.
    
    Input: Bisulfite BAM file
    Output: {sample}.UXM.tsv file with fragment-level methylation scores
    """
    from .assets import AssetManager
    
    # Configure Rust thread pool
    if threads > 0:
        try:
            _core.configure_threads(threads)
            logger.info(f"Configured {threads} threads for parallel processing")
        except Exception as e:
            logger.warning(f"Could not configure threads: {e}")
    
    # Input validation
    if not bam_input.exists():
        logger.error(f"Input BAM not found: {bam_input}")
        raise typer.Exit(1)
    
    if not str(bam_input).endswith('.bam'):
        logger.error(f"Input must be a .bam file: {bam_input}")
        raise typer.Exit(1)
    
    # Initialize Asset Manager
    try:
        assets = AssetManager(genome)
        logger.info(f"Genome: {assets.raw_genome} → {assets.genome_dir}")
    except ValueError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    
    # Default marker file via AssetManager
    if mark_input is None:
        try:
            mark_input = assets.resolve("methylation_markers")
            logger.info(f"Using default marker file: {mark_input}")
        except FileNotFoundError as e:
            logger.error(f"Methylation markers not found: {e}")
            logger.error(f"Please provide --mark-input or ensure markers exist for {genome}")
            raise typer.Exit(1)
    
    if not mark_input.exists():
        logger.error(f"Marker file not found: {mark_input}")
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bam_input.stem.replace('.bam', '')
    
    output_file = output / f"{sample_name}.UXM.tsv"
    
    try:
        logger.info(f"Processing {bam_input.name}")
        
        # Call Rust backend with all parameters
        _core.uxm.calculate_uxm(
            str(bam_input),
            str(mark_input),
            str(output_file),
            20,    # map_quality
            1,     # min_cpg
            0.5,   # methy_threshold
            0.5,   # unmethy_threshold
            "SE"   # pe_type (single-end default)
        )
        
        logger.info(f"UXM complete: {output_file}")

    except Exception as e:
        logger.error(f"UXM calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
