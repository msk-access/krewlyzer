"""
Orientation-aware cfDNA Fragmentation (OCF) calculation.

Calculates OCF features for a single sample.
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
logger = logging.getLogger("ocf")

# Rust backend is required
from krewlyzer import _core


def ocf(
    bedgz_input: Path = typer.Argument(..., help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output files (default: derived from input filename)"),
    ocr_input: Optional[Path] = typer.Option(None, "--ocr-input", "-r", help="Path to open chromatin regions file"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)")
):
    """
    Calculate Orientation-aware cfDNA Fragmentation (OCF) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.OCF.tsv (summary) and {sample}.OCF.sync.tsv (detailed sync data)
    """
    # Configure Rust thread pool
    if threads > 0:
        try:
            _core.configure_threads(threads)
            logger.info(f"Configured {threads} threads for parallel processing")
        except Exception as e:
            logger.warning(f"Could not configure threads: {e}")
    
    # Input validation
    if not bedgz_input.exists():
        logger.error(f"Input file not found: {bedgz_input}")
        raise typer.Exit(1)
    
    if not str(bedgz_input).endswith('.bed.gz'):
        logger.error(f"Input must be a .bed.gz file: {bedgz_input}")
        raise typer.Exit(1)
    
    # Default OCR file
    if ocr_input is None:
        pkg_dir = Path(__file__).parent
        ocr_input = pkg_dir / "data" / "OpenChromatinRegion" / "7specificTissue.all.OC.bed"
        logger.info(f"Using default OCR file: {ocr_input}")
    
    if not ocr_input.exists():
        logger.error(f"OCR file not found: {ocr_input}")
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    # Use a subdirectory for Rust output to avoid collisions/hardcoded names
    # Rust writes 'all.ocf.csv' and 'all.sync.tsv'
    sample_dir = output / sample_name
    sample_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        logger.info(f"Processing {bedgz_input.name}")
        
        # Call Rust backend
        _core.ocf.calculate_ocf(
            str(bedgz_input),
            str(ocr_input),
            str(sample_dir)
        )
        
        # Rename/Move files to krewlyzer standard: {output}/{sample}.{EXT}
        import shutil
        
        # Rust hardcoded outputs
        rust_ocf = sample_dir / "all.ocf.csv"
        rust_sync = sample_dir / "all.sync.tsv"
        
        # Standardized outputs
        final_ocf = output / f"{sample_name}.OCF.csv"
        final_sync = output / f"{sample_name}.OCF.sync.tsv"
        
        if rust_ocf.exists():
            shutil.move(str(rust_ocf), str(final_ocf))
        if rust_sync.exists():
            shutil.move(str(rust_sync), str(final_sync))
            
        # Clean up temporary sample dir if empty
        try:
            sample_dir.rmdir()
        except OSError:
            pass # Directory not empty or other error, leave it
        
        logger.info(f"OCF complete: {final_ocf}, {final_sync}")

    except Exception as e:
        logger.error(f"OCF calculation failed: {e}")
        raise typer.Exit(1)
