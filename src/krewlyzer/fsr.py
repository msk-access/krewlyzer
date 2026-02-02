"""
Fragment Size Ratio (FSR) calculation.

Calculates FSR features for a single sample with comprehensive counts and ratios.
This CLI is a thin wrapper around the unified processor.

Fragment Size Bins (as defined by Rust backend):
    - ultra_short: 65-99bp (TF footprints, highly tumor-specific)
    - core_short: 100-149bp (tumor-enriched, primary biomarker)
    - mono_nucl: 150-259bp (mono-nucleosomal fragments)
    - di_nucl: 260-399bp (di-nucleosomal fragments)
    - long: 400+bp (di-nucleosomal and larger)
    
Output Columns:
    Counts: short_count, long_count, total_count
    Ratios: short_long_ratio, short_long_log2
    Fractions: short_frac, long_frac
"""

import typer
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("fsr")

# Import asset resolution and startup banner
from .core.asset_resolution import resolve_target_regions, resolve_pon_model
from .core.logging import log_startup_banner, ResolvedAsset
from . import __version__


def fsr(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for hybrid GC correction"),
    pon_variant: str = typer.Option("all_unique", "--pon-variant", help="PON variant: 'all_unique' (default, max coverage) or 'duplex' (highest accuracy)"),
    skip_pon: bool = typer.Option(False, "--skip-pon", help="Skip PON z-score normalization"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED (for panel data: generates on/off-target FSR)"),
    skip_target_regions: bool = typer.Option(False, "--skip-target-regions", help="Disable panel mode even when --assay has bundled targets"),
    assay: Optional[str] = typer.Option(None, "--assay", "-A", help="Assay code (xs1/xs2) for auto-loading bundled assets"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
):
    """
    Calculate Fragment Size Ratio (FSR) features for a single sample.
    
    Computes fragment size distributions and ratios for cancer biomarker analysis.
    Uses Rust backend for performance with GC-corrected counts.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSR.tsv with short/long ratios per window
    """
    from .core.unified_processor import run_features
    from .assets import AssetManager
    
    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("krewlyzer.core.unified_processor").setLevel(logging.DEBUG)
    
    # Input validation
    if not bedgz_input.exists():
        logger.error(f"Input file not found: {bedgz_input}")
        raise typer.Exit(1)
    
    # Derive sample name
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    # Initialize AssetManager
    assets = AssetManager(genome)
    
    # ═══════════════════════════════════════════════════════════════════
    # ASSET RESOLUTION
    # ═══════════════════════════════════════════════════════════════════
    try:
        resolved_pon_path, pon_source = resolve_pon_model(
            explicit_path=pon_model, assay=assay, skip_pon=skip_pon,
            assets=assets, variant=pon_variant, log=logger
        )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)
    
    try:
        resolved_target_path, target_source = resolve_target_regions(
            explicit_path=target_regions, assay=assay, skip_target_regions=skip_target_regions,
            assets=assets, log=logger
        )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)
    
    is_panel_mode = resolved_target_path is not None and resolved_target_path.exists()
    
    # ═══════════════════════════════════════════════════════════════════
    # STARTUP BANNER
    # ═══════════════════════════════════════════════════════════════════
    log_startup_banner(
        tool_name="fsr", version=__version__,
        inputs={"Fragments": str(bedgz_input.name), "Output": str(output)},
        config={"Genome": f"{assets.raw_genome} → {assets.genome_dir}", "Assay": assay or "None", "Mode": "Panel" if is_panel_mode else "WGS"},
        assets=[ResolvedAsset("PON", resolved_pon_path, pon_source), ResolvedAsset("Targets", resolved_target_path, target_source)],
        logger=logger
    )
    
    try:
        outputs = run_features(
            bed_path=bedgz_input,
            output_dir=output,
            sample_name=sample_name,
            genome=genome,
            enable_fsr=True,
            target_regions=resolved_target_path,
            assay=assay,
            pon_model=resolved_pon_path,
            skip_pon_zscore=skip_pon,
            fsc_bins=bin_input,
            fsc_windows=windows,
            fsc_continue_n=continue_n,
            gc_correct=gc_correct,
            threads=threads,
            verbose=verbose,
        )
        
        if outputs.fsr and outputs.fsr.exists():
            logger.info(f"✅ FSR complete: {outputs.fsr}")
        if outputs.fsr_ontarget and outputs.fsr_ontarget.exists():
            logger.info(f"✅ FSR on-target: {outputs.fsr_ontarget}")
            
    except FileNotFoundError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    except RuntimeError as e:
        logger.error(f"FSR calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
