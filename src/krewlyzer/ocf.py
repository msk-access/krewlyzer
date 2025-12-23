"""
Orientation-aware cfDNA Fragmentation (OCF) calculation.

Calculates OCF features for a single sample.
Uses Rust backend via unified pipeline for accelerated computation with GC correction.
"""

import typer
from pathlib import Path
from typing import Optional
import logging
import shutil

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
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)")
):
    """
    Calculate Orientation-aware cfDNA Fragmentation (OCF) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.OCF.tsv (summary) and {sample}.OCF.sync.tsv (detailed sync data)
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
    if not bedgz_input.exists():
        logger.error(f"Input file not found: {bedgz_input}")
        raise typer.Exit(1)
    
    if not str(bedgz_input).endswith('.bed.gz'):
        logger.error(f"Input must be a .bed.gz file: {bedgz_input}")
        raise typer.Exit(1)
    
    # Initialize Asset Manager
    try:
        assets = AssetManager(genome)
        logger.info(f"Genome: {assets.raw_genome} -> {assets.genome_dir}")
    except ValueError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    
    # Default OCR file from assets
    if ocr_input is None:
        try:
            ocr_input = assets.resolve("ocf_regions")
            logger.info(f"Using default OCR file: {ocr_input}")
        except FileNotFoundError as e:
            logger.error(str(e))
            raise typer.Exit(1)
    
    if not ocr_input.exists():
        logger.error(f"OCR file not found: {ocr_input}")
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    # Use a subdirectory for Rust output to avoid collisions/hardcoded names
    # Rust writes 'all.ocf.tsv' and 'all.sync.tsv'
    sample_dir = output / f"{sample_name}_ocf_tmp"
    sample_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        logger.info(f"Processing {bedgz_input.name}")
        
        # Resolve GC correction assets
        gc_ref = None
        valid_regions = None
        factors_out = None
        
        if gc_correct:
            try:
                gc_ref = assets.resolve("gc_reference")
                valid_regions = assets.resolve("valid_regions")
                factors_out = output / f"{sample_name}.correction_factors.csv"
                logger.info(f"GC correction enabled using bundled assets for {genome}")
            except FileNotFoundError as e:
                logger.warning(f"GC correction assets not found: {e}")
                logger.warning("Proceeding without GC correction. Use --no-gc-correct to suppress this warning.")
                gc_correct = False
        
        # Check for pre-computed correction factors (from extract step)
        factors_input = None
        if gc_correct:
            potential_factors = bedgz_input.parent / f"{bedgz_input.stem.replace('.bed', '')}.correction_factors.csv"
            if potential_factors.exists():
                factors_input = potential_factors
                logger.info(f"Using pre-computed correction factors: {factors_input}")
                gc_ref = None
                valid_regions = None
                factors_out = None
        
        # Call Unified Pipeline (OCF only)
        logger.info("Running unified pipeline for OCF...")
        _core.run_unified_pipeline(
            str(bedgz_input),
            # GC Correction (compute)
            str(gc_ref) if gc_ref else None,
            str(valid_regions) if valid_regions else None,
            str(factors_out) if factors_out else None,
            # GC Correction (load pre-computed)
            str(factors_input) if factors_input else None,
            # FSC - disabled
            None, None,
            # WPS - disabled
            None, None, False,
            # FSD - disabled
            None, None,
            # OCF - enabled
            str(ocr_input), str(sample_dir)
        )
        
        # Rename/Move files to krewlyzer standard: {output}/{sample}.{EXT}
        # Rust hardcoded outputs
        rust_ocf = sample_dir / "all.ocf.tsv"
        rust_sync = sample_dir / "all.sync.tsv"
        
        # Standardized outputs
        final_ocf = output / f"{sample_name}.OCF.tsv"
        final_sync = output / f"{sample_name}.OCF.sync.tsv"
        
        if rust_ocf.exists():
            shutil.move(str(rust_ocf), str(final_ocf))
        if rust_sync.exists():
            shutil.move(str(rust_sync), str(final_sync))
            
        # Clean up temporary sample dir if empty
        try:
            sample_dir.rmdir()
        except OSError:
            pass  # Directory not empty or other error, leave it
        
        logger.info(f"OCF complete: {final_ocf}, {final_sync}")

    except Exception as e:
        logger.error(f"OCF calculation failed: {e}")
        raise typer.Exit(1)
