import typer
from pathlib import Path
from typing import Optional
import logging
import sys

import pysam

import numpy as np
from rich.console import Console
from rich.logging import RichHandler

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("fsr")

from .helpers import gc_correct

import pandas as pd

def _calc_fsr(
    bedgz_input: str | Path, 
    bin_input: str | Path, 
    windows: int, 
    continue_n: int, 
    output_file: str | Path
):
    """
    Internal: Calculate fragment size ratio (FSR) for a single .bed.gz file.
    Optimized with vectorized operations.
    """
    try:
        logger.info(f"Processing {bedgz_input} with bins from {bin_input}")
        
        # Load bins
        try:
            bins_df = pd.read_csv(bin_input, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'], dtype={'chrom': str, 'start': int, 'end': int})
        except Exception as e:
            logger.error(f"Could not load bins from {bin_input}: {e}")
            raise typer.Exit(1)
            
        try:
            tbx = pysam.TabixFile(filename=bedgz_input, mode="r")
        except Exception as e:
            logger.error(f"Could not open {bedgz_input} as Tabix file: {e}")
            raise typer.Exit(1)

        shorts_ratios = []
        ultra_shorts_ratios = []
        inter_ratios = []
        longs_ratios = []
        
        # Iterate over bins
        for _, bin_row in bins_df.iterrows():
            chrom = bin_row['chrom']
            start = bin_row['start']
            end = bin_row['end']
            
            try:
                rows = list(tbx.fetch(chrom, start, end, parser=pysam.asTuple()))
            except ValueError:
                rows = []
            except Exception as e:
                logger.error(f"Error fetching {chrom}:{start}-{end}: {e}")
                raise typer.Exit(1)
            
            if not rows:
                shorts_ratios.append(0)
                ultra_shorts_ratios.append(0)
                inter_ratios.append(0)
                longs_ratios.append(0)
                continue
                
            try:
                # Vectorized parsing
                _, starts, ends, _ = zip(*rows)
                starts = np.array(starts, dtype=int)
                ends = np.array(ends, dtype=int)
                lengths = ends - starts
                
                # Filter 65-400
                mask = (lengths >= 65) & (lengths <= 400)
                valid_lengths = lengths[mask]
                
                total = len(valid_lengths)
                
                if total == 0:
                    shorts_ratios.append(0)
                    ultra_shorts_ratios.append(0)
                    inter_ratios.append(0)
                    longs_ratios.append(0)
                else:
                    shorts = np.sum((valid_lengths >= 65) & (valid_lengths <= 150))
                    ultra_shorts = np.sum((valid_lengths >= 65) & (valid_lengths <= 100))
                    intermediates = np.sum((valid_lengths >= 151) & (valid_lengths <= 260))
                    longs = np.sum((valid_lengths >= 261) & (valid_lengths <= 400))
                    
                    shorts_ratios.append(shorts / total)
                    ultra_shorts_ratios.append(ultra_shorts / total)
                    inter_ratios.append(intermediates / total)
                    longs_ratios.append(longs / total)
                    
            except Exception as e:
                logger.error(f"Error processing data in bin {chrom}:{start}-{end}: {e}")
                raise typer.Exit(1)

        # Aggregation into windows
        df = pd.DataFrame({
            'chrom': bins_df['chrom'],
            'start': bins_df['start'],
            'end': bins_df['end'],
            'short_r': shorts_ratios,
            'ultra_short_r': ultra_shorts_ratios,
            'inter_r': inter_ratios,
            'long_r': longs_ratios
        })
        
        results = []
        
        for chrom, group in df.groupby('chrom', sort=False):
            n_bins = len(group)
            n_windows = n_bins // continue_n
            
            if n_windows == 0:
                continue
                
            trunc_len = n_windows * continue_n
            
            short_mat = group['short_r'].values[:trunc_len].reshape(n_windows, continue_n)
            ultra_short_mat = group['ultra_short_r'].values[:trunc_len].reshape(n_windows, continue_n)
            inter_mat = group['inter_r'].values[:trunc_len].reshape(n_windows, continue_n)
            long_mat = group['long_r'].values[:trunc_len].reshape(n_windows, continue_n)
            
            # Mean of ratios
            mean_short = short_mat.mean(axis=1)
            mean_ultra_short = ultra_short_mat.mean(axis=1)
            mean_inter = inter_mat.mean(axis=1)
            mean_long = long_mat.mean(axis=1)
            
            window_starts = np.arange(n_windows) * continue_n * windows
            window_ends = (np.arange(n_windows) + 1) * continue_n * windows - 1
            
            results.append(pd.DataFrame({
                'chrom': chrom,
                'start': window_starts,
                'end': window_ends,
                'short_mean': mean_short,
                'ultra_short_mean': mean_ultra_short,
                'inter_mean': mean_inter,
                'long_mean': mean_long
            }))
            
        if not results:
            logger.warning("No valid windows found.")
            return

        final_df = pd.concat(results, ignore_index=True)
        
        # Write output
        with open(output_file, 'w') as f:
            f.write("region\tshort-ratio\tultra-short-ratio\titermediate-ratio\tlong-ratio\n")
            for _, row in final_df.iterrows():
                region = f"{row['chrom']}:{int(row['start'])}-{int(row['end'])}"
                f.write(f"{region}\t{row['short_mean']:.4f}\t{row['ultra_short_mean']:.4f}\t{row['inter_mean']:.4f}\t{row['long_mean']:.4f}\n")
                
        logger.info(f"FSR calculation complete. Results written to {output_file}")

    except Exception as e:
        logger.error(f"Fatal error in _calc_fsr: {e}")
        raise typer.Exit(1)

def _run_fsr_file(bedgz_file, output_dir, bin_input, windows, continue_n):
    """Module-level worker function for parallel FSR calculation."""
    output_file = Path(output_dir) / (Path(bedgz_file).stem.replace('.bed', '') + '.FSR.txt')
    _calc_fsr(str(bedgz_file), str(bin_input), windows, continue_n, str(output_file))
    return str(output_file)


def fsr(
    bedgz_path: Path = typer.Argument(..., help="Folder containing .bed.gz files (should be the output directory from motif.py)"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file (default: data/ChormosomeBins/hg19_window_100kb.bed)"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size (default: 100000)"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number (default: 50)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output folder for results"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of parallel processes (default: 1)")
):
    """
    Calculate fragment size ratio (FSR) features for all .bed.gz files in a folder.
    The input folder should be the output directory produced by motif.py, containing the .bed.gz files.
    Output files are written to the output directory, one per .bed.gz file.
    """
    # Input checks
    if not bedgz_path.exists():
        logger.error(f"Input directory not found: {bedgz_path}")
        raise typer.Exit(1)
    if bin_input and not bin_input.exists():
        logger.error(f"Bin input file not found: {bin_input}")
        raise typer.Exit(1)
    try:
        output.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Could not create output directory {output}: {e}")
        raise typer.Exit(1)
    if not output.exists():
        output.mkdir(parents=True, exist_ok=True)
    bedgz_files = [f for f in bedgz_path.iterdir() if f.suffixes == ['.bed', '.gz']]
    if not bedgz_files:
        logger.error("No .bed.gz files found in the specified folder.")
        raise typer.Exit(1)
    if bin_input is None:
        # Use package-level data (inside krewlyzer/data/)
        pkg_dir = Path(__file__).parent
        bin_input = pkg_dir / "data" / "ChormosomeBins" / "hg19_window_100kb.bed"
        logger.info(f"No bin_input specified. Using default: {bin_input}")
    if not bin_input.exists():
        logger.error(f"Bin input file does not exist: {bin_input}")
        raise typer.Exit(1)
    logger.info(f"Calculating FSR for {len(bedgz_files)} files...")
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from functools import partial
    logger.info(f"Starting parallel FSR calculation using {threads} processes...")
    
    worker = partial(_run_fsr_file, output_dir=str(output), bin_input=str(bin_input),
                     windows=windows, continue_n=continue_n)
    
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(worker, str(bedgz_file)): bedgz_file for bedgz_file in bedgz_files}
        for future in as_completed(futures):
            bedgz_file = futures[future]
            try:
                result = future.result()
                logger.info(f"FSR calculated: {result}")
            except Exception as exc:
                logger.error(f"FSR calculation failed for {bedgz_file}: {exc}")
    logger.info(f"FSR features calculated for {len(bedgz_files)} files.")

