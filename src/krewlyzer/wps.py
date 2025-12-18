"""
Windowed Protection Score (WPS) calculation.

Calculates unified WPS features (Long, Short, Ratio) for a single sample.
Uses Rust backend for accelerated computation.
"""

import typer
from pathlib import Path
from typing import Optional
import logging
import json

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("wps")

# Rust backend is required
from krewlyzer import _core


def wps(
    bedgz_input: Path = typer.Argument(..., help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file (default: derived from input filename)"),
    tsv_input: Optional[Path] = typer.Option(None, "--tsv-input", "-t", help="Path to transcript/region TSV file"),
    reference: Optional[Path] = typer.Option(None, "--reference", "-r", help="Reference FASTA for GC computation (required for GC correction)"),
    empty: bool = typer.Option(False, "--empty/--no-empty", help="Include regions with no coverage"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction using LOESS (default: True, requires --reference)"),
    threads: int = typer.Option(0, "--threads", "-p", help="Number of threads (0=all cores)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging")
):
    """
    Calculate unified Windowed Protection Score (WPS) features for a single sample.
    
    Calculates both Long WPS (nucleosome, 120-180bp) and Short WPS (TF, 35-80bp)
    in a single pass, plus their ratio and normalized versions.
    
    Input: .bed.gz file from extract step
    Output: {sample}.WPS.tsv.gz file with columns:
        - gene_id, chrom, pos
        - cov_long, cov_short (coverage)
        - wps_long, wps_short, wps_ratio (raw WPS)
        - wps_long_norm, wps_short_norm, wps_ratio_norm (normalized per million)
    
    GC Correction:
        By default, computes region GC content from the reference FASTA.
        Use --no-gc-correct to disable, or provide --reference to enable.
    """
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
    
    # GC correction requires reference
    if gc_correct and reference is None:
        logger.warning("GC correction enabled but --reference not provided. Disabling GC correction.")
        logger.warning("Use --reference (-r) to enable GC correction, or --no-gc-correct to suppress this warning.")
        gc_correct = False
    
    if reference and not reference.exists():
        logger.error(f"Reference FASTA not found: {reference}")
        raise typer.Exit(1)
    
    # Default transcript file
    if tsv_input is None:
        pkg_dir = Path(__file__).parent
        tsv_input = pkg_dir / "data" / "TranscriptAnno" / "transcriptAnno-hg19-1kb.tsv"
        logger.info(f"Using default transcript file: {tsv_input}")
    
    if not tsv_input.exists():
        logger.error(f"Transcript file not found: {tsv_input}")
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    try:
        logger.info(f"Processing {bedgz_input.name}")
        logger.info("Calculating unified WPS (Long: 120-180bp, Short: 35-80bp)")
        
        if gc_correct:
            logger.info(f"GC correction enabled using reference: {reference.name}")
        
        # Get metadata if available
        total_fragments = None
        metadata_file = str(bedgz_input).replace('.bed.gz', '.metadata.json')
        if Path(metadata_file).exists():
            try:
                with open(metadata_file, 'r') as f:
                    meta = json.load(f)
                    total_fragments = meta.get('total_fragments')
                    if total_fragments:
                        logger.info(f"Loaded metadata: {total_fragments:,} fragments")
            except Exception as e:
                logger.warning(f"Could not load metadata: {e}")
        else:
            logger.warning(f"Metadata file not found: {metadata_file}")
            logger.warning("WPS normalized columns (wps_*_norm) will NOT be comparable across samples!")
            logger.warning("Run 'krewlyzer extract' first to generate metadata.json")

        
        # Call Rust backend (unified calculation with GC correction support)
        reference_path = str(reference) if reference else None
        count = _core.calculate_wps(
            str(bedgz_input),
            str(tsv_input),
            str(output),
            sample_name,
            empty,
            total_fragments,
            reference_path,
            gc_correct,
            verbose
        )
        
        logger.info(f"WPS complete: processed {count} regions → {output}/{sample_name}.WPS.tsv.gz")
        logger.info("Output columns: gene_id, chrom, pos, cov_long, cov_short, wps_long, wps_short, wps_ratio, wps_long_norm, wps_short_norm, wps_ratio_norm")

    except Exception as e:
        logger.error(f"WPS calculation failed: {e}")
        raise typer.Exit(1)

