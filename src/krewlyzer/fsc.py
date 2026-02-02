"""
Fragment Size Coverage (FSC) calculation.

Calculates FSC features for a single sample showing fragment coverage depth
across genomic windows, stratified by fragment size bin.

This CLI is a thin wrapper around the unified processor.
See krewlyzer.core.unified_processor for the core implementation.

Fragment Size Bins (from Rust backend):
    - ultra_short: 65-99bp (TF footprints)
    - core_short: 100-149bp (tumor-enriched)
    - mono_nucl: 150-259bp (mono-nucleosomal)
    - di_nucl: 260-399bp (di-nucleosomal)
    - long: 400+bp

Output: {sample}.FSC.tsv with per-window coverage counts and z-scores.
Panel mode: Additional {sample}.FSC.ontarget.tsv for on-target regions.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("fsc")

# Import asset resolution functions
from .core.asset_resolution import resolve_target_regions, resolve_pon_model

# Import startup banner logging
from .core.logging import log_startup_banner, ResolvedAsset
from . import __version__


def fsc(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for hybrid GC correction"),
    pon_variant: str = typer.Option("all_unique", "--pon-variant", help="PON variant: 'all_unique' (default, max coverage) or 'duplex' (highest accuracy)"),
    skip_pon: bool = typer.Option(False, "--skip-pon", help="Skip PON z-score normalization"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED (for panel data: generates on/off-target FSC)"),
    skip_target_regions: bool = typer.Option(False, "--skip-target-regions", help="Disable panel mode even when --assay has bundled targets"),
    assay: Optional[str] = typer.Option(None, "--assay", "-A", help="Assay type (xs1/xs2) for gene-centric FSC aggregation"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size (default: 100000)"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
):
    """
    Calculate fragment size coverage (FSC) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSC.tsv file with z-scored fragment size coverage per window
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
    
    # Validate user-provided override files
    from .core.asset_validation import validate_file, FileSchema
    if bin_input and bin_input.exists():
        logger.debug(f"Validating user-provided bin file: {bin_input}")
        validate_file(bin_input, FileSchema.BED3)
    if target_regions and target_regions.exists():
        logger.debug(f"Validating user-provided target regions: {target_regions}")
        validate_file(target_regions, FileSchema.BED3)
    
    # Derive sample name
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    # Initialize AssetManager
    assets = AssetManager(genome)
    
    # ═══════════════════════════════════════════════════════════════════
    # ASSET RESOLUTION
    # ═══════════════════════════════════════════════════════════════════
    # Resolve PON model
    try:
        resolved_pon_path, pon_source = resolve_pon_model(
            explicit_path=pon_model,
            assay=assay,
            skip_pon=skip_pon,
            assets=assets,
            variant=pon_variant,
            log=logger
        )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)
    
    # Resolve target regions
    try:
        resolved_target_path, target_source = resolve_target_regions(
            explicit_path=target_regions,
            assay=assay,
            skip_target_regions=skip_target_regions,
            assets=assets,
            log=logger
        )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)
    
    is_panel_mode = resolved_target_path is not None and resolved_target_path.exists()
    
    # ═══════════════════════════════════════════════════════════════════
    # STARTUP BANNER
    # ═══════════════════════════════════════════════════════════════════
    log_startup_banner(
        tool_name="fsc",
        version=__version__,
        inputs={
            "Fragments": str(bedgz_input.name),
            "Output": str(output),
        },
        config={
            "Genome": f"{assets.raw_genome} → {assets.genome_dir}",
            "Assay": assay or "None",
            "Mode": "Panel" if is_panel_mode else "WGS",
        },
        assets=[
            ResolvedAsset("PON", resolved_pon_path, pon_source),
            ResolvedAsset("Targets", resolved_target_path, target_source),
        ],
        logger=logger
    )
    
    try:
        # Call unified processor with FSC enabled
        outputs = run_features(
            bed_path=bedgz_input,
            output_dir=output,
            sample_name=sample_name,
            genome=genome,
            enable_fsc=True,
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
        
        # Report results
        if outputs.fsc and outputs.fsc.exists():
            logger.info(f"✅ FSC complete: {outputs.fsc}")
        if outputs.fsc_ontarget and outputs.fsc_ontarget.exists():
            logger.info(f"✅ FSC on-target: {outputs.fsc_ontarget}")
        if outputs.fsc_gene and outputs.fsc_gene.exists():
            logger.info(f"✅ FSC gene: {outputs.fsc_gene}")
            
    except FileNotFoundError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    except RuntimeError as e:
        logger.error(f"FSC calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
