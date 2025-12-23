"""
Windowed Protection Score (WPS) calculation.

Calculates unified WPS features (Long, Short, Ratio) for a single sample.
Uses Rust backend via unified pipeline for accelerated computation with GC correction.
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
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    tsv_input: Optional[Path] = typer.Option(None, "--tsv-input", "-T", help="Path to transcript/region TSV file"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for z-score computation"),
    empty: bool = typer.Option(False, "--empty/--no-empty", help="Include regions with no coverage"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
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
        - wps_long_norm, wps_short_norm, wps_ratio_norm (normalized)
    
    With --pon-model: Adds wps_long_z and wps_short_z columns (z-scores vs PON)
    """
    from .assets import AssetManager
    
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
    
    # Default transcript file from assets
    if tsv_input is None:
        try:
            tsv_input = assets.resolve("transcript_anno")
            logger.info(f"Using default transcript file: {tsv_input}")
        except FileNotFoundError as e:
            logger.error(str(e))
            raise typer.Exit(1)
    
    if not tsv_input.exists():
        logger.error(f"Transcript file not found: {tsv_input}")
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    # Output file path
    output_file = output / f"{sample_name}.WPS.tsv.gz"
    
    try:
        logger.info(f"Processing {bedgz_input.name}")
        logger.info("Calculating unified WPS (Long: 120-180bp, Short: 35-80bp)")
        
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
        
        # Call Unified Pipeline (WPS only)
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
            # WPS - enabled
            str(tsv_input), str(output_file), empty,
            # FSD - disabled
            None, None,
            # OCF - disabled
            None, None
        )
        
        logger.info(f"WPS complete: {output_file}")
        
        # If PON provided, compute z-scores for each region
        if pon_model:
            from krewlyzer.pon.model import PonModel
            import pandas as pd
            import gzip
            
            try:
                pon = PonModel.load(pon_model)
                logger.info(f"Loaded PON model: {pon.assay} (n={pon.n_samples})")
                
                if pon.wps_baseline:
                    logger.info("Computing z-scores against PON baseline...")
                    
                    # Read WPS output
                    df = pd.read_csv(output_file, sep="\t", compression="gzip")
                    
                    # Aggregate per-region means
                    region_stats = df.groupby("gene_id").agg({
                        "wps_long": "mean",
                        "wps_short": "mean"
                    }).reset_index()
                    
                    # Merge with PON baseline
                    pon_df = pon.wps_baseline.regions
                    merged = region_stats.merge(pon_df, left_on="gene_id", right_on="region_id", how="left")
                    
                    # Compute z-scores
                    merged["wps_long_z"] = (merged["wps_long"] - merged["wps_long_mean"]) / merged["wps_long_std"].replace(0, 1)
                    merged["wps_short_z"] = (merged["wps_short"] - merged["wps_short_mean"]) / merged["wps_short_std"].replace(0, 1)
                    
                    # Save summary file
                    summary_file = output / f"{sample_name}.WPS_summary.tsv"
                    summary_cols = ["gene_id", "wps_long", "wps_short", "wps_long_z", "wps_short_z"]
                    merged[summary_cols].to_csv(summary_file, sep="\t", index=False)
                    logger.info(f"Saved z-score summary: {summary_file}")
                    
            except Exception as e:
                logger.warning(f"Could not compute PON z-scores: {e}")

    except Exception as e:
        logger.error(f"WPS calculation failed: {e}")
        raise typer.Exit(1)
