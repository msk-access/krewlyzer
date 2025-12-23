"""
Fragment Size Coverage (FSC) calculation.

Calculates FSC features for a single sample.
Uses Rust backend for accelerated computation with GC correction.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

import numpy as np
import pandas as pd
from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("fsc")

# Rust backend is required
from krewlyzer import _core


def fsc(
    bedgz_input: Path = typer.Argument(..., help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for hybrid GC correction"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size (default: 100000)"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)")
):
    """
    Calculate fragment size coverage (FSC) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSC.tsv file with z-scored fragment size coverage per window
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
    
    # Initialize Asset Manager
    try:
        assets = AssetManager(genome)
    except ValueError as e:
        logger.error(str(e))
        raise typer.Exit(1)
        
    # Resolve bin file
    if bin_input is None:
        try:
            bin_input = assets.resolve("bins_100kb")
            logger.info(f"Using default bin file: {bin_input}")
        except FileNotFoundError as e:
            logger.error(str(e))
            raise typer.Exit(1)
            
    if not bin_input.exists():
        logger.error(f"Bin file not found: {bin_input}")
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    output_file = output / f"{sample_name}.FSC.tsv"
    fsc_counts_file = output / f"{sample_name}.fsc_counts.tsv"
    
    # Load PON model if provided
    pon = None
    if pon_model:
        from krewlyzer.pon.model import PonModel
        try:
            pon = PonModel.load(pon_model)
            logger.info(f"Loaded PON model: {pon.assay} (n={pon.n_samples})")
        except Exception as e:
            logger.warning(f"Could not load PON model: {e}")
            logger.warning("Falling back to within-sample correction only")

    try:
        logger.info(f"Processing {bedgz_input.name}")
        
        # Resolve GC correction assets
        gc_ref = None
        valid_regions = None
        factors_out = None
        
        if gc_correct:
            try:
                gc_ref = assets.resolve("gc_correction")
                valid_regions = assets.resolve("valid_regions")
                factors_out = output / f"{sample_name}.correction_factors.csv"
            except FileNotFoundError as e:
                logger.warning(f"GC correction assets not found: {e}")
                logger.warning("Proceeding without GC correction")
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
        
        # Call Unified Pipeline (FSC only)
        logger.info("Running fragment counting...")
        _core.run_unified_pipeline(
            str(bedgz_input),
            # GC Correction (compute)
            str(gc_ref) if gc_ref else None,
            str(valid_regions) if valid_regions else None,
            str(factors_out) if factors_out else None,
            # GC Correction (load pre-computed)
            str(factors_input) if factors_input else None,
            # FSC
            str(bin_input), str(fsc_counts_file),
            # Others
            None, None, False, # WPS
            None, None,        # FSD
            None, None         # OCF
        )
        
        logger.info(f"Reading counts from {fsc_counts_file}")
        
        # Load the FSC counts TSV
        # Format: chrom, start, end, ultra_short, short, intermediate, long, total, mean_gc
        # Using pandas
        df_counts = pd.read_csv(fsc_counts_file, sep='\t')
        
        # Extract arrays for downstream logic
        short = df_counts['short'].values
        intermediate = df_counts['intermediate'].values
        long = df_counts['long'].values
        total = df_counts['total'].values
        gc_arr = df_counts['mean_gc'].values
        
        # Step 2: Apply PON correction if available (overlay on corrected counts)
        if pon and pon.gc_bias:
            logger.info("Applying PON-based normalization...")
            for i, gc in enumerate(gc_arr):
                short_exp = pon.gc_bias.get_expected(gc, "short")
                inter_exp = pon.gc_bias.get_expected(gc, "intermediate")
                long_exp = pon.gc_bias.get_expected(gc, "long")
                
                # Correct: observed / expected
                # Note: 'observed' here is already GC-corrected by factors if gc_correct=True
                if short_exp > 0:
                    short[i] /= short_exp
                if inter_exp > 0:
                    intermediate[i] /= inter_exp
                if long_exp > 0:
                    long[i] /= long_exp
            
            # Recalculate total after PON
            total = short + intermediate + long
            
        logger.info(f"Processed {len(total)} bins")
        
        # Create DataFrame for aggregation (reusing loaded columns)
        df = pd.DataFrame({
            'chrom': df_counts['chrom'],
            'start': df_counts['start'],
            'end': df_counts['end'],
            'shorts': short,
            'intermediates': intermediate,
            'longs': long,
            'totals': total
        })
        
        # Aggregation into windows
        results = []
        for chrom, group in df.groupby('chrom', sort=False):
            n_bins = len(group)
            n_windows = n_bins // continue_n
            
            if n_windows == 0:
                continue
            
            trunc_len = n_windows * continue_n
            
            shorts_mat = group['shorts'].values[:trunc_len].reshape(n_windows, continue_n)
            inter_mat = group['intermediates'].values[:trunc_len].reshape(n_windows, continue_n)
            longs_mat = group['longs'].values[:trunc_len].reshape(n_windows, continue_n)
            totals_mat = group['totals'].values[:trunc_len].reshape(n_windows, continue_n)
            
            sum_shorts = shorts_mat.sum(axis=1)
            sum_inter = inter_mat.sum(axis=1)
            sum_longs = longs_mat.sum(axis=1)
            sum_totals = totals_mat.sum(axis=1)
            
            window_starts = np.arange(n_windows) * continue_n * windows
            window_ends = (np.arange(n_windows) + 1) * continue_n * windows - 1
            
            results.append(pd.DataFrame({
                'chrom': chrom,
                'start': window_starts,
                'end': window_ends,
                'short_sum': sum_shorts,
                'inter_sum': sum_inter,
                'long_sum': sum_longs,
                'total_sum': sum_totals
            }))
        
        if not results:
            logger.warning("No valid windows found.")
            return
        
        final_df = pd.concat(results, ignore_index=True)
        
        # Calculate Z-scores globally
        final_df['short_z'] = (final_df['short_sum'] - final_df['short_sum'].mean()) / final_df['short_sum'].std()
        final_df['inter_z'] = (final_df['inter_sum'] - final_df['inter_sum'].mean()) / final_df['inter_sum'].std()
        final_df['long_z'] = (final_df['long_sum'] - final_df['long_sum'].mean()) / final_df['long_sum'].std()
        final_df['total_z'] = (final_df['total_sum'] - final_df['total_sum'].mean()) / final_df['total_sum'].std()
        
        # Write output
        with open(output_file, 'w') as f:
            f.write("region\tshort-fragment-zscore\titermediate-fragment-zscore\tlong-fragment-zscore\ttotal-fragment-zscore\n")
            for _, row in final_df.iterrows():
                region = f"{row['chrom']}:{int(row['start'])}-{int(row['end'])}"
                f.write(f"{region}\t{row['short_z']:.4f}\t{row['inter_z']:.4f}\t{row['long_z']:.4f}\t{row['total_z']:.4f}\n")
        
        logger.info(f"FSC complete: {output_file}")

    except Exception as e:
        logger.error(f"FSC calculation failed: {e}")
        raise typer.Exit(1)
