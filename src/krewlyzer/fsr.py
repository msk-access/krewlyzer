"""
Fragment Size Ratio (FSR) calculation.

Calculates FSR features for a single sample with comprehensive counts and ratios.
Uses Rust backend for accelerated computation.
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
logger = logging.getLogger("fsr")

# Rust backend is required
from krewlyzer import _core


def fsr(
    bedgz_input: Path = typer.Argument(..., help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file (default: derived from input filename)"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file (default: hg19_window_100kb.bed)"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size (default: 100000)"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number (default: 50)"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction using LOESS (default: True)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging")
):
    """
    Calculate Fragment Size Ratio (FSR) features for a single sample.
    
    Outputs comprehensive fragment counts and ratios for cancer biomarker analysis.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSR.txt with columns:
        Counts: ultra_short_count, short_count, inter_count, long_count, total_count
        Ratios: short_ratio, inter_ratio, long_ratio, short_long_ratio, ultra_short_ratio
    
    Size categories (matching cfDNAFE):
        - Ultra-short: 65-100bp (TF footprints)
        - Short: 65-149bp (tumor-enriched)
        - Intermediate: 151-259bp (nucleosome dynamics)
        - Long: 261-399bp (healthy cell contribution)
        - Total: 65-399bp
    
    GC Correction:
        By default, per-fragment-type LOESS correction is applied to remove GC bias
        before calculating ratios. Use --no-gc-correct to disable.
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
    
    output_file = output / f"{sample_name}.FSR.tsv"
    
    try:
        logger.info(f"Processing {bedgz_input.name}")
        
        if gc_correct:
            # Use GC-corrected counts from Rust backend
            logger.info("Counting fragments with GC correction (LOESS per fragment type)...")
            short, intermediate, long, gc_values = _core.count_fragments_gc_corrected(
                str(bedgz_input),
                str(bin_input),
                verbose
            )
            # Note: GC-corrected version doesn't return ultra_short or total (derived from corrected)
            # For ultra_short ratio, we need raw counts - use non-corrected for ultra_short only
            ultra_short_raw, _, _, _, total_raw, _ = _core.count_fragments_by_bins(
                str(bedgz_input),
                str(bin_input)
            )
            ultra_short = np.array(ultra_short_raw)
            # Total from corrected values
            total = np.array(short) + np.array(intermediate) + np.array(long)
            short = np.array(short)
            intermediate = np.array(intermediate)
            long = np.array(long)
            logger.info(f"GC correction applied to short/intermediate/long counts")
        else:
            # Use raw counts (no GC correction)
            logger.info("Counting fragments (no GC correction)...")
            ultra_short, short, intermediate, long, total, gc_values = _core.count_fragments_by_bins(
                str(bedgz_input),
                str(bin_input)
            )
            ultra_short = np.array(ultra_short)
            short = np.array(short)
            intermediate = np.array(intermediate)
            long = np.array(long)
            total = np.array(total)
        
        logger.info(f"Processed {len(total)} bins")
        
        # Load bin file for coordinates
        bins_df = pd.read_csv(bin_input, sep='\t', header=None, usecols=[0, 1, 2], 
                            names=['chrom', 'start', 'end'], dtype={'chrom': str, 'start': int, 'end': int})
        
        # Create DataFrame with counts (GC-corrected if gc_correct=True)
        df = pd.DataFrame({
            'chrom': bins_df['chrom'],
            'start': bins_df['start'],
            'end': bins_df['end'],
            'ultra_short': ultra_short,
            'short': short,
            'intermediate': intermediate,
            'long': long,
            'total': total
        })
        
        # Aggregation into windows
        results = []
        for chrom, group in df.groupby('chrom', sort=False):
            n_bins = len(group)
            n_windows = n_bins // continue_n
            
            if n_windows == 0:
                continue
            
            trunc_len = n_windows * continue_n
            
            # Sum counts within aggregated windows
            ultra_short_sums = group['ultra_short'].values[:trunc_len].reshape(n_windows, continue_n).sum(axis=1)
            short_sums = group['short'].values[:trunc_len].reshape(n_windows, continue_n).sum(axis=1)
            inter_sums = group['intermediate'].values[:trunc_len].reshape(n_windows, continue_n).sum(axis=1)
            long_sums = group['long'].values[:trunc_len].reshape(n_windows, continue_n).sum(axis=1)
            total_sums = group['total'].values[:trunc_len].reshape(n_windows, continue_n).sum(axis=1)
            
            window_starts = np.arange(n_windows) * continue_n * windows
            window_ends = (np.arange(n_windows) + 1) * continue_n * windows - 1
            
            for i in range(n_windows):
                region = f"{chrom}:{window_starts[i]}-{window_ends[i]}"
                
                ultra_count = int(ultra_short_sums[i])
                short_count = int(short_sums[i])
                inter_count = int(inter_sums[i])
                long_count = int(long_sums[i])
                total_count = int(total_sums[i])
                
                # Calculate ratios (avoid division by zero)
                if total_count > 0:
                    short_ratio = short_count / total_count
                    inter_ratio = inter_count / total_count
                    long_ratio = long_count / total_count
                    ultra_ratio = ultra_count / total_count
                else:
                    short_ratio = inter_ratio = long_ratio = ultra_ratio = 0.0
                
                # Short/Long ratio (primary cancer biomarker)
                if long_count > 0:
                    short_long_ratio = short_count / long_count
                else:
                    short_long_ratio = short_count if short_count > 0 else 0.0
                
                results.append({
                    'region': region,
                    'ultra_short_count': ultra_count,
                    'short_count': short_count,
                    'inter_count': inter_count,
                    'long_count': long_count,
                    'total_count': total_count,
                    'short_ratio': short_ratio,
                    'inter_ratio': inter_ratio,
                    'long_ratio': long_ratio,
                    'short_long_ratio': short_long_ratio,
                    'ultra_short_ratio': ultra_ratio
                })
        
        # Write output
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
        
        logger.info(f"FSR complete: {len(results_df)} windows → {output_file}")
        logger.info("Output columns: counts (ultra_short, short, inter, long, total) + ratios (short, inter, long, short_long, ultra_short)")

    except Exception as e:
        logger.error(f"FSR calculation failed: {e}")
        raise typer.Exit(1)
