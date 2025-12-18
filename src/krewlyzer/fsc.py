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
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file (default: derived from input filename)"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file (default: hg19_window_100kb.bed)"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size (default: 100000)"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number (default: 50)"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction using LOESS (default: True)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)")
):
    """
    Calculate fragment size coverage (FSC) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSC.txt file with z-scored fragment size coverage per window
    
    GC Correction:
        By default, per-fragment-type LOESS correction is applied to remove GC bias.
        Use --no-gc-correct to disable.
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
    
    # Default bin input
    if bin_input is None:
        pkg_dir = Path(__file__).parent
        bin_input = pkg_dir / "data" / "ChormosomeBins" / "hg19_window_100kb.bed"
        logger.info(f"Using default bin file: {bin_input}")
    
    if not bin_input.exists():
        logger.error(f"Bin file not found: {bin_input}")
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    output_file = output / f"{sample_name}.FSC.tsv"
    
    try:
        logger.info(f"Processing {bedgz_input.name}")
        
        # Count fragments using Rust backend with optional GC correction
        if gc_correct:
            logger.info("Counting fragments with GC correction (LOESS per fragment type)...")
            short, intermediate, long, gc_values = _core.count_fragments_gc_corrected(
                str(bedgz_input),
                str(bin_input),
                verbose
            )
            # Compute total from corrected values
            total = np.array(short) + np.array(intermediate) + np.array(long)
            short = np.array(short)
            intermediate = np.array(intermediate)
            long = np.array(long)
            logger.info(f"GC correction applied to short/intermediate/long counts")
        else:
            logger.info("Counting fragments (no GC correction)...")
            _, short, intermediate, long, total, gc_values = _core.count_fragments_by_bins(
                str(bedgz_input),
                str(bin_input)
            )
            short = np.array(short)
            intermediate = np.array(intermediate)
            long = np.array(long)
            total = np.array(total)
        
        logger.info(f"Processed {len(total)} bins")
        
        # Load bin file for coordinates
        bins_df = pd.read_csv(bin_input, sep='\t', header=None, usecols=[0, 1, 2], 
                            names=['chrom', 'start', 'end'], dtype={'chrom': str, 'start': int, 'end': int})
        
        # Create DataFrame with results (already GC-corrected if enabled)
        df = pd.DataFrame({
            'chrom': bins_df['chrom'],
            'start': bins_df['start'],
            'end': bins_df['end'],
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
