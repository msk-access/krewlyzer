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


def fsc(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for hybrid GC correction"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED (for panel data: generates on/off-target FSC)"),
    assay: Optional[str] = typer.Option(None, "--assay", "-A", help="Assay type (xs1/xs2) for gene-centric FSC aggregation"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size (default: 100000)"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    format: Optional[str] = typer.Option(None, "--format", "-f", help="Output format override: tsv, parquet, json (default: tsv)")
):
    """
    Calculate fragment size coverage (FSC) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSC.tsv file with z-scored fragment size coverage per window
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
    
    # Derive sample name
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    try:
        # Call unified processor with FSC enabled
        outputs = run_features(
            bed_path=bedgz_input,
            output_dir=output,
            sample_name=sample_name,
            genome=genome,
            enable_fsc=True,
            target_regions=target_regions,
            assay=assay,
            pon_model=pon_model,
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
