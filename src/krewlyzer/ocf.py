"""
Orientation-aware cfDNA Fragmentation (OCF) calculation.

Calculates OCF features showing fragment orientation patterns around open chromatin regions.
This CLI is a thin wrapper around the unified processor.

OCF measures tissue-of-origin based on end-motif orientation at regulatory elements.

Output Files:
    - {sample}.OCF.tsv: Per-region OCF scores and fragment counts
    - {sample}.OCF.sync.tsv: Synchronized regions for panel comparison
    - Panel mode: {sample}.OCF.ontarget.tsv for on-target regions

OCF Score: Measures strand asymmetry of fragment ends = (U-D)/(U+D) where U=upstream, D=downstream.
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


def ocf(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output files"),
    ocr_input: Optional[Path] = typer.Option(None, "--ocr-input", "-r", help="Path to open chromatin regions file"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED (for panel data: generates on/off-target OCF)"),
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
    
    try:
        # Call unified processor with OCF enabled
        outputs = run_features(
            bed_path=bedgz_input,
            output_dir=output,
            sample_name=sample_name,
            genome=genome,
            enable_ocf=True,
            target_regions=target_regions,
            pon_model=pon_model,
            skip_pon_zscore=skip_pon,
            ocf_regions=ocr_input,
            gc_correct=gc_correct,
            threads=threads,
            verbose=verbose,
        )
        
        # Report results
        if outputs.ocf and outputs.ocf.exists():
            logger.info(f"✅ OCF complete: {outputs.ocf}")
        if outputs.ocf_sync and outputs.ocf_sync.exists():
            logger.info(f"✅ OCF sync: {outputs.ocf_sync}")
        if outputs.ocf_ontarget and outputs.ocf_ontarget.exists():
            logger.info(f"✅ OCF on-target: {outputs.ocf_ontarget}")
            
    except FileNotFoundError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    except RuntimeError as e:
        logger.error(f"OCF calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
