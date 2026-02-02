"""
Fragment Size Distribution (FSD) calculation.

Calculates FSD features showing fragment size histograms per chromosome arm.
This CLI is a thin wrapper around the unified processor.

Key biomarker for cancer detection based on cfDNA fragmentation patterns.

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

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("fsd")

# Import asset resolution and startup banner
from .core.asset_resolution import resolve_target_regions, resolve_pon_model
from .core.logging import log_startup_banner, ResolvedAsset
from . import __version__


def fsd(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    arms_file: Optional[Path] = typer.Option(None, "--arms-file", "-a", help="Path to chromosome arms BED file"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED (for panel data: generates on/off-target FSD)"),
    skip_target_regions: bool = typer.Option(False, "--skip-target-regions", help="Disable panel mode even when --assay has bundled targets"),
    assay: Optional[str] = typer.Option(None, "--assay", "-A", help="Assay code (xs1/xs2) for auto-loading bundled assets"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for z-score computation"),
    pon_variant: str = typer.Option("all_unique", "--pon-variant", help="PON variant: 'all_unique' (default, max coverage) or 'duplex' (highest accuracy)"),
    skip_pon: bool = typer.Option(False, "--skip-pon", help="Skip PON z-score normalization"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
):
    """
    Calculate fragment size distribution (FSD) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSD.tsv file with fragment size histogram per chromosome arm
    
    With --pon-model: Additional z-score columns comparing to PON baseline
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
    if arms_file and arms_file.exists():
        logger.debug(f"Validating user-provided arms file: {arms_file}")
        validate_file(arms_file, FileSchema.ARMS_BED)
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
        tool_name="fsd", version=__version__,
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
            enable_fsd=True,
            target_regions=resolved_target_path,
            assay=assay,
            pon_model=resolved_pon_path,
            skip_pon_zscore=skip_pon,
            fsd_arms=arms_file,
            gc_correct=gc_correct,
            threads=threads,
            verbose=verbose,
        )
        
        if outputs.fsd and outputs.fsd.exists():
            logger.info(f"✅ FSD complete: {outputs.fsd}")
        if outputs.fsd_ontarget and outputs.fsd_ontarget.exists():
            logger.info(f"✅ FSD on-target: {outputs.fsd_ontarget}")
            
    except FileNotFoundError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    except RuntimeError as e:
        logger.error(f"FSD calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
