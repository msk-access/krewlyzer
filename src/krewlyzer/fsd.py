"""
Fragment Size Distribution (FSD) calculation.

Calculates FSD features for a single sample.
Uses Rust backend for accelerated computation.
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
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for z-score computation"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)")
):
    """
    Calculate fragment size distribution (FSD) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSD.tsv file with fragment size histogram per chromosome arm
    
    With --pon-model: Additional z-score columns comparing to PON baseline
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
    
    # Default arms file
    if arms_file is None:
        pkg_dir = Path(__file__).parent
        arms_file = pkg_dir / "data" / "ChormosomeArms" / "hg19.arms.bed.gz"
        logger.info(f"Using default arms file: {arms_file}")
    
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
        
        # Call Rust backend
        _core.fsd.calculate_fsd(
            str(bedgz_input),
            str(arms_file),
            str(output_file)
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
