"""
Windowed Protection Score (WPS) calculation.

Calculates unified WPS features (Long, Short, Ratio) for a single sample.
Uses Rust backend via unified pipeline for accelerated computation with GC correction.

Dual-Stream Processing:
- WPS-Nuc (Nucleosome): Weighted fragments (1.0 for [160,175], 0.5 for [120,159]∪[176,180])
- WPS-TF (Transcription Factor): Binary [35,80] fragments
"""

import typer
from pathlib import Path
from typing import Optional
import logging
import json

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("wps")

# Rust backend is required
from krewlyzer import _core


def wps(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    tsv_input: Optional[Path] = typer.Option(None, "--tsv-input", "-W", help="Path to transcript/region TSV file (legacy)"),
    wps_anchors: Optional[Path] = typer.Option(None, "--wps-anchors", help="WPS anchors BED (merged TSS+CTCF) for dual-stream profiling"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Panel capture BED (enables bait edge masking)"),
    bait_padding: int = typer.Option(50, "--bait-padding", help="Bait edge padding in bp (default 50, use 15-20 for small exon panels)"),
    background: Optional[Path] = typer.Option(None, "--background", "-B", help="Background Alu BED for hierarchical stacking (auto-loaded if not specified)"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for z-score computation"),
    empty: bool = typer.Option(False, "--empty/--no-empty", help="Include regions with no coverage"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging")
):
    """
    Calculate unified Windowed Protection Score (WPS) features for a single sample.
    
    Dual-stream weighted fragment processing:
    - WPS-Nuc (Nucleosome): 120bp window, weighted fragments
      * Primary [160,175bp]: weight 1.0
      * Secondary [120,159]∪[176,180bp]: weight 0.5
    - WPS-TF (Transcription Factor): 16bp window, fragments [35,80bp]
    
    Input: .bed.gz file from extract step
    Output: {sample}.WPS.tsv.gz file with columns:
        - gene_id, chrom, pos
        - cov_long, cov_short (coverage)
        - wps_long, wps_short, wps_ratio (raw WPS)
        - wps_long_norm, wps_short_norm, wps_ratio_norm (normalized)
    
    With --pon-model: Adds wps_long_z and wps_short_z columns (z-scores vs PON)
    """
    from .assets import AssetManager
    from .core.pon_integration import load_pon_model
    from .core.wps_processor import apply_wps_pon
    
    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("core.wps_processor").setLevel(logging.DEBUG)
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
    
    # Resolve WPS anchor file (prefer wps_anchors, fallback to tsv_input for legacy)
    regions_file = None
    if wps_anchors is not None:
        if wps_anchors.exists():
            regions_file = wps_anchors
            logger.info(f"WPS anchors: {wps_anchors}")
        else:
            logger.error(f"WPS anchors file not found: {wps_anchors}")
            raise typer.Exit(1)
    elif tsv_input is not None:
        if tsv_input.exists():
            regions_file = tsv_input
            logger.info(f"Using legacy transcript file: {tsv_input}")
        else:
            logger.error(f"Transcript file not found: {tsv_input}")
            raise typer.Exit(1)
    else:
        # Try new wps_anchors first, fallback to legacy transcript_anno
        try:
            regions_file = assets.resolve("wps_anchors")
            logger.info(f"Using default WPS anchors: {regions_file}")
        except FileNotFoundError:
            try:
                regions_file = assets.resolve("transcript_anno")
                logger.info(f"Using legacy transcript file: {regions_file}")
            except FileNotFoundError as e:
                logger.error(f"No WPS anchor or transcript file found: {e}")
                raise typer.Exit(1)
    
    # Log panel mode if target_regions provided
    if target_regions is not None:
        if target_regions.exists():
            logger.info(f"Panel mode: bait edge masking enabled ({target_regions.name})")
        else:
            logger.warning(f"Target regions file not found: {target_regions}")
            target_regions = None
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    # Output file path (Parquet only)
    output_file = output / f"{sample_name}.WPS.parquet"
    
    try:
        logger.info(f"Processing {bedgz_input.name}")
        logger.info("WPS dual-stream: Nuc (weighted 120-180bp), TF (35-80bp)")
        
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
            # Look for existing correction_factors.csv next to input BED.gz
            potential_factors = bedgz_input.parent / f"{bedgz_input.stem.replace('.bed', '')}.correction_factors.csv"
            if potential_factors.exists():
                factors_input = potential_factors
                logger.info(f"Using pre-computed correction factors: {factors_input}")
                # Don't recompute, just load
                gc_ref = None
                valid_regions = None
                factors_out = None
        
        # Output paths (Parquet only)
        output_file = output / f"{file_stem}.WPS.parquet"
        output_bg_file = output / f"{file_stem}.WPS_background.parquet"
        
        # Load background Alu regions for hierarchical stacking
        # Explicit --background takes precedence over auto-load
        if background and background.exists():
            bg_regions = background
            logger.info(f"Using explicit background: {bg_regions}")
        else:
            try:
                bg_regions = assets.resolve("wps_background")
                logger.info(f"Auto-loaded background ({genome}): {bg_regions}")
            except FileNotFoundError:
                bg_regions = None
                logger.warning("Background Alu regions not found - skipping background WPS")
        
        # Call Unified Pipeline (WPS with foreground + background)
        logger.info("Running unified pipeline for WPS...")
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
            # WPS Foreground
            str(regions_file), str(output_file),
            # WPS Background (Alu stacking)
            str(bg_regions) if bg_regions else None, 
            str(output_bg_file) if bg_regions else None,
            empty,
            # FSD - disabled
            None, None,
            # OCF - disabled
            None, None
        )
        
        logger.info(f"WPS complete: {output_file}")
        if bg_regions and output_bg_file.exists():
            logger.info(f"WPS background: {output_bg_file}")
        
        # Post-processing: smoothing, FFT periodicity
        from .core.wps_processor import post_process_wps
        wps_result = post_process_wps(
            wps_parquet=output_file,
            wps_background_parquet=output_bg_file if output_bg_file.exists() else None,
            pon_baseline_parquet=None,  # TODO: extract from pon_model if provided
            smooth=True,
            extract_periodicity=True
        )
        if wps_result.get("smoothed"):
            logger.info("Applied Savitzky-Golay smoothing to WPS profiles")
        if wps_result.get("periodicity_extracted"):
            logger.info(f"Extracted FFT periodicity (NRL score: {wps_result.get('periodicity_score', 0):.3f})")

    except Exception as e:
        logger.error(f"WPS calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)

