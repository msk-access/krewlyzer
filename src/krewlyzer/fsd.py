"""
Fragment Size Distribution (FSD) calculation.

Calculates FSD features showing fragment size histograms per chromosome arm.
Key biomarker for cancer detection based on cfDNA fragmentation patterns.

Uses Rust backend via unified pipeline for accelerated computation with GC correction.

Output Structure:
    - Rows: Chromosome arms (e.g., chr1p, chr1q, chr2p, ...)
    - Columns: Fragment counts per size (65bp to 400bp)
    - With --pon-model: Z-scores vs healthy baseline per arm

Panel Mode: Generates both off-target (.FSD.tsv) and on-target (.FSD.ontarget.tsv) outputs.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

import pandas as pd
from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("fsd")

# Rust backend is required
from krewlyzer import _core


def fsd(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    arms_file: Optional[Path] = typer.Option(None, "--arms-file", "-a", help="Path to chromosome arms BED file"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED (for panel data: generates on/off-target FSD)"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for z-score computation"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    format: Optional[str] = typer.Option(None, "--format", "-f", help="Output format override: tsv, parquet, json (default: tsv)")
):
    """
    Calculate fragment size distribution (FSD) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSD.tsv file with fragment size histogram per chromosome arm
    
    With --pon-model: Additional z-score columns comparing to PON baseline
    """
    from .assets import AssetManager
    from .core.pon_integration import load_pon_model
    from .core.fsd_processor import apply_fsd_pon
    
    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("core.fsd_processor").setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")
    
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
    
    # Default arms file from assets
    if arms_file is None:
        try:
            arms_file = assets.resolve("arms")
            logger.info(f"Using default arms file: {arms_file}")
        except FileNotFoundError as e:
            logger.error(str(e))
            raise typer.Exit(1)
    
    if not arms_file.exists():
        logger.error(f"Arms file not found: {arms_file}")
        raise typer.Exit(1)
    
    logger.debug(f"Arms file: {arms_file}")
    
    # Load PON model if provided
    pon = None
    if pon_model:
        pon = load_pon_model(pon_model)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    output_file = output / f"{sample_name}.FSD.tsv"
    logger.debug(f"Sample: {sample_name}, Output: {output_file}")
    
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
                factors_out = output / f"{sample_name}.correction_factors.tsv"
                logger.info(f"GC correction enabled using bundled assets for {genome}")
            except FileNotFoundError as e:
                logger.warning(f"GC correction assets not found: {e}")
                logger.warning("Proceeding without GC correction. Use --no-gc-correct to suppress this warning.")
                gc_correct = False
        
        # Check for pre-computed correction factors (from extract step)
        factors_input = None
        if gc_correct:
            potential_factors = bedgz_input.parent / f"{bedgz_input.stem.replace('.bed', '')}.correction_factors.tsv"
            if potential_factors.exists():
                factors_input = potential_factors
                logger.info(f"Using pre-computed correction factors: {factors_input}")
                gc_ref = None
                valid_regions = None
                factors_out = None
        
        # Call Unified Pipeline (FSD only)
        logger.info("Running unified pipeline for FSD...")
        
        is_panel_mode = target_regions and target_regions.exists()
        if is_panel_mode:
            logger.info(f"Panel mode: on/off-target split enabled (targets: {target_regions.name})")
        
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
            None, None,
            # WPS Background - disabled
            None, None, False,
            # FSD - enabled
            str(arms_file), str(output_file),
            # OCF - disabled
            None, None,
            # Target regions for on/off-target split
            str(target_regions) if is_panel_mode else None,
            50,  # bait_padding
            False  # silent
        )
        
        # If PON provided, apply z-scores using shared processor
        if pon:
            apply_fsd_pon(output_file, pon)
        
        logger.info(f"✅ FSD complete: {output_file}")

    except Exception as e:
        logger.error(f"FSD calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
