"""
Windowed Protection Score (WPS) calculation.

Calculates WPS features for a single sample measuring nucleosome positioning.
This CLI is a thin wrapper around the unified processor.

Dual-stream weighted fragment processing:
- WPS-Nuc (Nucleosome): 120bp window, weighted fragments
  * Primary [160,175bp]: weight 1.0
  * Secondary [120,159]∪[176,180bp]: weight 0.5
- WPS-TF (Transcription Factor): 16bp window, fragments [35,80bp]

Output: {sample}.WPS.parquet with columns:
    - gene_id, chrom, pos
    - cov_long, cov_short (coverage)
    - wps_long, wps_short, wps_ratio (raw WPS)
    
With --pon-model: Adds wps_long_z and wps_short_z columns (z-scores vs PON)
With --assay: Adds {sample}.WPS.panel.parquet for panel-specific anchors
"""

import typer
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("wps")


def wps(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    wps_anchors: Optional[Path] = typer.Option(None, "--wps-anchors", help="WPS anchors BED (merged TSS+CTCF) for dual-stream profiling"),
    assay: Optional[str] = typer.Option(None, "--assay", "-A", help="Assay code (xs1, xs2) for dual WPS output (genome-wide + panel)"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Panel capture BED (enables bait edge masking)"),
    bait_padding: int = typer.Option(50, "--bait-padding", help="Bait edge padding in bp (default 50, use 15-20 for small exon panels)"),
    background: Optional[Path] = typer.Option(None, "--background", "-B", help="Background Alu BED for hierarchical stacking (auto-loaded if not specified)"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for z-score computation"),
    skip_pon: bool = typer.Option(False, "--skip-pon", help="Skip PON z-score normalization"),
    empty: bool = typer.Option(False, "--empty/--no-empty", help="Include regions with no coverage"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
    format: Optional[str] = typer.Option(None, "--format", "-f", help="Output format override: tsv, parquet, json (default: parquet)")
):
    """
    Calculate unified Windowed Protection Score (WPS) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.WPS.parquet with nucleosome positioning scores per region
    """
    from .core.unified_processor import run_features
    
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
    if wps_anchors and wps_anchors.exists():
        logger.debug(f"Validating user-provided WPS anchors: {wps_anchors}")
        validate_file(wps_anchors, FileSchema.WPS_ANCHORS)
    if background and background.exists():
        logger.debug(f"Validating user-provided WPS background: {background}")
        validate_file(background, FileSchema.WPS_ANCHORS)
    if target_regions and target_regions.exists():
        logger.debug(f"Validating user-provided target regions: {target_regions}")
        validate_file(target_regions, FileSchema.BED3)
    
    # Derive sample name
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    try:
        # Call unified processor with WPS enabled
        outputs = run_features(
            bed_path=bedgz_input,
            output_dir=output,
            sample_name=sample_name,
            genome=genome,
            enable_wps=True,
            target_regions=target_regions,
            assay=assay,
            pon_model=pon_model,
            skip_pon_zscore=skip_pon,
            wps_anchors=wps_anchors,
            wps_background=background,
            wps_bait_padding=bait_padding,
            wps_empty=empty,
            gc_correct=gc_correct,
            threads=threads,
            verbose=verbose,
        )
        
        # Report results
        if outputs.wps and outputs.wps.exists():
            logger.info(f"✅ WPS complete: {outputs.wps}")
        if outputs.wps_background and outputs.wps_background.exists():
            logger.info(f"✅ WPS background: {outputs.wps_background}")
        if outputs.wps_panel and outputs.wps_panel.exists():
            logger.info(f"✅ WPS panel: {outputs.wps_panel}")
            
    except FileNotFoundError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    except RuntimeError as e:
        logger.error(f"WPS calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
