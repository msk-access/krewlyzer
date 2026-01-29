"""
Orientation-aware cfDNA Fragmentation (OCF) calculation.

Calculates OCF features showing fragment orientation patterns around open chromatin regions.
This CLI is a thin wrapper around the unified processor.

OCF measures tissue-of-origin based on end-motif orientation at regulatory elements.
OCF Score = (Upstream - Downstream) / (Upstream + Downstream)

Output Files (WGS mode: 2 files, Panel mode: 6 files):
    - {sample}.OCF.tsv: Summary scores (all fragments)
    - {sample}.OCF.sync.tsv: Detailed sync data (all fragments)
    - {sample}.OCF.ontarget.tsv: Summary scores (on-target only)
    - {sample}.OCF.ontarget.sync.tsv: Detailed sync (on-target)
    - {sample}.OCF.offtarget.tsv: Summary scores (off-target only)
    - {sample}.OCF.offtarget.sync.tsv: Detailed sync (off-target)
"""

import typer
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("ocf")

# Import asset resolution and startup banner
from .core.asset_resolution import resolve_target_regions, resolve_pon_model
from .core.logging import log_startup_banner, ResolvedAsset
from . import __version__


def ocf(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output files"),
    ocr_input: Optional[Path] = typer.Option(None, "--ocr-input", "-r", help="Path to open chromatin regions file"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED (for panel data: generates on/off-target OCF)"),
    skip_target_regions: bool = typer.Option(False, "--skip-target-regions", help="Disable panel mode even when --assay has bundled targets"),
    assay: Optional[str] = typer.Option(None, "--assay", "-A", help="Assay code (xs1/xs2) for auto-loading bundled assets"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for z-score computation"),
    skip_pon: bool = typer.Option(False, "--skip-pon", help="Skip PON z-score normalization"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
):
    """
    Calculate Orientation-aware cfDNA Fragmentation (OCF) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.OCF.tsv (summary) and {sample}.OCF.sync.tsv (detailed sync data)
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
    if ocr_input and ocr_input.exists():
        logger.debug(f"Validating user-provided OCR file: {ocr_input}")
        validate_file(ocr_input, FileSchema.REGION_BED)
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
            assets=assets, log=logger
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
        tool_name="ocf", version=__version__,
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
            enable_ocf=True,
            target_regions=resolved_target_path,
            assay=assay,
            pon_model=resolved_pon_path,
            skip_pon_zscore=skip_pon,
            ocf_regions=ocr_input,
            gc_correct=gc_correct,
            threads=threads,
            verbose=verbose,
        )
        
        if outputs.ocf and outputs.ocf.exists():
            logger.info(f"✅ OCF complete: {outputs.ocf}")
        if outputs.ocf_sync and outputs.ocf_sync.exists():
            logger.info(f"  OCF sync: {outputs.ocf_sync}")
        if outputs.ocf_ontarget and outputs.ocf_ontarget.exists():
            logger.info(f"  OCF on-target: {outputs.ocf_ontarget}")
        if outputs.ocf_offtarget and outputs.ocf_offtarget.exists():
            logger.info(f"  OCF off-target: {outputs.ocf_offtarget}")
            
    except FileNotFoundError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    except RuntimeError as e:
        logger.error(f"OCF calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
