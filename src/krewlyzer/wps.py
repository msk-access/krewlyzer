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
    assay: Optional[str] = typer.Option(None, "--assay", "-A", help="Assay code (xs1, xs2) for dual WPS output (genome-wide + panel)"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Panel capture BED (enables bait edge masking)"),
    bait_padding: int = typer.Option(50, "--bait-padding", help="Bait edge padding in bp (default 50, use 15-20 for small exon panels)"),
    background: Optional[Path] = typer.Option(None, "--background", "-B", help="Background Alu BED for hierarchical stacking (auto-loaded if not specified)"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for z-score computation"),
    empty: bool = typer.Option(False, "--empty/--no-empty", help="Include regions with no coverage"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
    format: Optional[str] = typer.Option(None, "--format", "-f", help="Output format override: tsv, parquet, json (default: parquet)")
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
    from .core.wps_processor import subtract_pon_baseline
    
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
    
    # Resolve assay for dual WPS output
    resolved_assay = assay if isinstance(assay, str) else None
    panel_anchors = None
    if resolved_assay:
        try:
            panel_anchors = assets.get_wps_anchors(resolved_assay)
            logger.info(f"Dual WPS enabled: panel anchors ({resolved_assay}) will generate WPS.panel.parquet")
        except Exception:
            logger.warning(f"Panel anchors not found for {resolved_assay}, dual WPS disabled")
    
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
        
        # Resolve GC correction assets (centralized helper)
        from .core.gc_assets import resolve_gc_assets
        gc = resolve_gc_assets(assets, output, sample_name, bedgz_input, gc_correct, genome)
        gc_ref = gc.gc_ref
        valid_regions = gc.valid_regions
        factors_out = gc.factors_out
        factors_input = gc.factors_input
        gc_correct = gc.gc_correct_enabled
        
        # Output paths (Parquet only)
        output_file = output / f"{sample_name}.WPS.parquet"
        output_bg_file = output / f"{sample_name}.WPS_background.parquet"
        
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
        
        # === DUAL WPS: Run panel-specific WPS if --assay was provided ===
        output_panel_file = None
        if panel_anchors and panel_anchors.exists():
            output_panel_file = output / f"{sample_name}.WPS.panel.parquet"
            logger.info(f"Running panel-specific WPS ({panel_anchors.name})...")
            _core.run_unified_pipeline(
                str(bedgz_input),
                None, None, None,  # GC already computed
                str(factors_out) if factors_out and factors_out.exists() else (str(factors_input) if factors_input else None),
                None, None,  # No FSC
                str(panel_anchors), str(output_panel_file),  # Panel WPS
                None, None,  # No background for panel
                empty,
                None, None,  # No FSD
                None, None,  # No OCF
                str(target_regions) if target_regions else None,
                bait_padding
            )
            logger.info(f"Panel WPS complete: {output_panel_file}")
        
        # Load PON model if provided and extract WPS baseline
        pon_wps_baseline_path = None
        if pon_model and pon_model.exists():
            try:
                pon = load_pon_model(pon_model)
                if pon and pon.wps_baseline and pon.wps_baseline.regions is not None:
                    pon_wps_baseline_path = output / f".{sample_name}.pon_wps_baseline.parquet"
                    pon.wps_baseline.regions.to_parquet(pon_wps_baseline_path, index=False)
                    logger.info(f"Loaded PON WPS baseline: {len(pon.wps_baseline.regions)} regions")
            except Exception as e:
                logger.warning(f"Failed to load PON model: {e}")
        
        # Post-processing: smoothing, PON subtraction, FFT periodicity
        from .core.wps_processor import post_process_wps
        wps_result = post_process_wps(
            wps_parquet=output_file,
            wps_background_parquet=output_bg_file if output_bg_file.exists() else None,
            pon_baseline_parquet=pon_wps_baseline_path,
            smooth=True,
            extract_periodicity=True
        )
        if wps_result.get("smoothed"):
            logger.info("Applied Savitzky-Golay smoothing to WPS profiles")
        if wps_result.get("pon_subtracted"):
            logger.info("Subtracted PON baseline (added *_delta and *_z columns)")
        if wps_result.get("periodicity_extracted"):
            score = wps_result.get('periodicity_score')
            if score is not None:
                logger.info(f"Extracted FFT periodicity (NRL score: {score:.3f})")
            else:
                logger.info("Extracted FFT periodicity")

    except Exception as e:
        logger.error(f"WPS calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)

