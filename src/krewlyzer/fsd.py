"""
Fragment Size Distribution (FSD) calculation.

Calculates FSD features for a single sample.
Uses Rust backend via unified pipeline for accelerated computation with GC correction.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

import pandas as pd
from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("fsd")

# Rust backend is required
from krewlyzer import _core


def fsd(
    bedgz_input: Path = typer.Argument(..., help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    arms_file: Optional[Path] = typer.Option(None, "--arms-file", "-a", help="Path to chromosome arms BED file"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for z-score computation"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)")
):
    """
    Calculate fragment size distribution (FSD) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSD.tsv file with fragment size histogram per chromosome arm
    
    With --pon-model: Additional z-score columns comparing to PON baseline
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
    
    # Default arms file from assets
    if arms_file is None:
        try:
            arms_file = assets.resolve("arms")
            logger.info(f"Using default arms file: {arms_file}")
        except FileNotFoundError as e:
            logger.error(str(e))
            raise typer.Exit(1)
    
    if not arms_file.exists():
        logger.error(f"Arms file not found: {arms_file}")
        raise typer.Exit(1)
    
    # Load PON model if provided
    pon = None
    if pon_model:
        from krewlyzer.pon.model import PonModel
        try:
            pon = PonModel.load(pon_model)
            logger.info(f"Loaded PON model: {pon.assay} (n={pon.n_samples})")
        except Exception as e:
            logger.warning(f"Could not load PON model: {e}")
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    output_file = output / f"{sample_name}.FSD.tsv"
    
    try:
        logger.info(f"Processing {bedgz_input.name}")
        
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
            potential_factors = bedgz_input.parent / f"{bedgz_input.stem.replace('.bed', '')}.correction_factors.csv"
            if potential_factors.exists():
                factors_input = potential_factors
                logger.info(f"Using pre-computed correction factors: {factors_input}")
                gc_ref = None
                valid_regions = None
                factors_out = None
        
        # Call Unified Pipeline (FSD only)
        logger.info("Running unified pipeline for FSD...")
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
            # WPS - disabled
            None, None, False,
            # FSD - enabled
            str(arms_file), str(output_file),
            # OCF - disabled
            None, None
        )
        
        # If PON provided, add z-score columns
        if pon and pon.fsd_baseline:
            logger.info("Computing z-scores against PON baseline...")
            df = pd.read_csv(output_file, sep="\t")
            
            # Compute z-scores for each arm-size combination
            z_scores = []
            for _, row in df.iterrows():
                arm = row.iloc[0]  # First column is arm name
                size_bin = int(row.name) if isinstance(row.name, int) else None
                
                # For now, compute overall z-score per arm
                # Future: Per-size-bin z-scores
                if arm in pon.fsd_baseline.arms and size_bin:
                    expected = pon.fsd_baseline.get_expected(arm, size_bin)
                    std = pon.fsd_baseline.get_std(arm, size_bin)
                    if std > 0:
                        z = (row.iloc[1] - expected) / std  # Assume 2nd col is count
                    else:
                        z = 0.0
                else:
                    z = None
                z_scores.append(z)
            
            # Add z-score column and save
            df["pon_zscore"] = z_scores
            df.to_csv(output_file, sep="\t", index=False)
            logger.info(f"Added z-score column from PON baseline")
        
        logger.info(f"FSD complete: {output_file}")

    except Exception as e:
        logger.error(f"FSD calculation failed: {e}")
        raise typer.Exit(1)
