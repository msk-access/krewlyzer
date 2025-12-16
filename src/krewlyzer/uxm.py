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
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)")
):
    """
    Calculate Fragment-level Methylation (UXM) features for a single sample.
    
    Input: Bisulfite BAM file
    Output: {sample}.UXM.tsv file with fragment-level methylation scores
    """
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
    
    # Default marker file
    if mark_input is None:
        pkg_dir = Path(__file__).parent
        mark_input = pkg_dir / "data" / "methylation-markers" / "uxm_markers_hg19.bed"
        if mark_input.exists():
            logger.info(f"Using default marker file: {mark_input}")
        else:
            logger.error("No default marker file found. Please provide --mark-input")
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
        raise typer.Exit(1)
