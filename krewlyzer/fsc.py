import typer
from pathlib import Path
from typing import Optional
import logging
import sys

import pysam

import numpy as np
import pandas as pd
from skmisc.loess import loess
from rich.console import Console
from rich.logging import RichHandler

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("fsc")


from .helpers import gc_correct


def _calc_fsc(
    bedgz_input: str | Path, 
    bin_input: str | Path, 
    windows: int, 
    continue_n: int, 
    output_file: str | Path
):
    """
    Internal: Calculate fragment size coverage (FSC) for a single .bed.gz file.
    Optimized with vectorized operations.
    """
    try:
        logger.info(f"Processing {bedgz_input} with bins from {bin_input}")
        
        # Load bins
        try:
            # Use pandas for faster loading if possible, but is fine for iteration
            # bins = pybedtools.BedTool(bin_input)
            # Actually, reading bins into a dataframe might be easier for grouping
            bins_df = pd.read_csv(bin_input, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'], dtype={'chrom': str, 'start': int, 'end': int})
        except Exception as e:
            logger.error(f"Could not load bins from {bin_input}: {e}")
            raise typer.Exit(1)
            
        try:
            tbx = pysam.TabixFile(filename=bedgz_input, mode="r")
        except Exception as e:
            logger.error(f"Could not open {bedgz_input} as Tabix file: {e}")
            raise typer.Exit(1)

        shorts_data = []
        intermediates_data = []
        longs_data = []
        totals_data = []
        bingc = []
        
        # Iterate over bins
        # To optimize, we can process by chromosome to avoid random seeking if bins are sorted?
        # But Tabix is good at random access.
        
        for _, bin_row in bins_df.iterrows():
            chrom = bin_row['chrom']
            start = bin_row['start']
            end = bin_row['end']
            
            try:
                # Fetch reads in bin
                # parser=pysam.asTuple() is slightly faster than splitting string manually
                rows = list(tbx.fetch(chrom, start, end, parser=pysam.asTuple()))
            except ValueError:
                # Region not in file (e.g. chromosome not present)
                rows = []
            except Exception as e:
                logger.error(f"Error fetching {chrom}:{start}-{end}: {e}")
                raise typer.Exit(1)
            
            if not rows:
                bingc.append(np.nan)
                shorts_data.append(0)
                intermediates_data.append(0)
                longs_data.append(0)
                totals_data.append(0)
                continue
                
            # Vectorize parsing
            # rows is a list of tuples: (chrom, start, end, gc)
            # We need start (idx 1), end (idx 2), gc (idx 3)
            # Note: BED is 0-based start, 1-based end? 
            # motif.py writes: rstart, rend, gc.
            # pysam.TabixFile returns what's in the file.
            # motif.py writes: f"{read1.reference_name}\t{rstart}\t{rend}\t{gc}\n"
            # So col 1 is start, col 2 is end.
            
            try:
                # Extract columns
                # This is still a Python loop but faster than append in loop
                # We can use zip to transpose
                _, starts, ends, gcs = zip(*rows)
                starts = np.array(starts, dtype=int)
                ends = np.array(ends, dtype=int)
                gcs = np.array(gcs, dtype=float)
                
                lengths = ends - starts
                
                # Filter by length (65-400)
                # We only care about reads with length in range [65, 400] for GC calculation?
                # Original code: if 65 <= len <= 400: gc.append(...)
                
                mask = (lengths >= 65) & (lengths <= 400)
                valid_gcs = gcs[mask]
                
                if len(valid_gcs) == 0:
                    bingc.append(np.nan)
                else:
                    bingc.append(np.mean(valid_gcs))
                
                # Counts
                # np.bincount requires non-negative integers
                # We can filter lengths to be within [0, 400] for bincount safety, though mask handles 65-400
                
                # We need counts for specific ranges
                # shorts: 65-150
                # intermediates: 151-260
                # longs: 261-400
                # totals: 65-400
                
                # Using histogram might be cleaner or just boolean indexing
                shorts = np.sum((lengths >= 65) & (lengths <= 150))
                intermediates = np.sum((lengths >= 151) & (lengths <= 260))
                longs = np.sum((lengths >= 261) & (lengths <= 400))
                totals = np.sum(mask) # 65-400
                
                shorts_data.append(shorts)
                intermediates_data.append(intermediates)
                longs_data.append(longs)
                totals_data.append(totals)
                
            except Exception as e:
                logger.error(f"Error processing data in bin {chrom}:{start}-{end}: {e}")
                raise typer.Exit(1)

        # GC Correction
        try:
            correct_shorts = gc_correct(shorts_data, bingc)
            correct_intermediates = gc_correct(intermediates_data, bingc)
            correct_longs = gc_correct(longs_data, bingc)
            correct_totals = gc_correct(totals_data, bingc)
        except Exception as e:
            logger.error(f"GC correction failed: {e}")
            raise typer.Exit(1)
            
        # Aggregation into windows
        # The original logic aggregates 'continue_n' bins into one window.
        # It assumes bins are contiguous and ordered by chromosome.
        # It resets when chromosome changes.
        
        # We can do this more pandas-style
        df = pd.DataFrame({
            'chrom': bins_df['chrom'],
            'start': bins_df['start'], # bin start
            'end': bins_df['end'], # bin end
            'shorts': correct_shorts,
            'intermediates': correct_intermediates,
            'longs': correct_longs,
            'totals': correct_totals
        })
        
        results = []
        
        # Group by chromosome
        for chrom, group in df.groupby('chrom', sort=False):
            # Rolling window or block aggregation?
            # Original code:
            # num = chrom.count(chrom[step]) -> count of bins in this chrom
            # continues_bin = num // continue_n
            # for i in range(continues_bin):
            #    combine...
            # This is block aggregation (non-overlapping blocks of size continue_n)
            
            n_bins = len(group)
            n_windows = n_bins // continue_n
            
            # We can reshape the array to (n_windows, continue_n) and sum along axis 1
            # But we need to handle the remainder (last_bin) which is ignored in original code?
            # Original: "if last_bin != 0: step += last_bin; start = 0" -> It skips the remainder?
            # "for i in range(continues_bin): ... step += continue_n"
            # It seems it drops the last partial window.
            
            if n_windows == 0:
                continue
                
            # Truncate to multiple of continue_n
            trunc_len = n_windows * continue_n
            
            shorts_mat = group['shorts'].values[:trunc_len].reshape(n_windows, continue_n)
            inter_mat = group['intermediates'].values[:trunc_len].reshape(n_windows, continue_n)
            longs_mat = group['longs'].values[:trunc_len].reshape(n_windows, continue_n)
            totals_mat = group['totals'].values[:trunc_len].reshape(n_windows, continue_n)
            
            sum_shorts = shorts_mat.sum(axis=1)
            sum_inter = inter_mat.sum(axis=1)
            sum_longs = longs_mat.sum(axis=1)
            sum_totals = totals_mat.sum(axis=1)
            
            # Calculate window coordinates
            # bin_start = start * windows
            # bin_end = (start + continue_n) * windows - 1
            # 'windows' param is actually the bin size? Or the window size?
            # CLI says: --windows "-w" help="Window size (default: 100000)"
            # CLI says: --bin-input "-b" help="Path to bin file (default: .../hg19_window_100kb.bed)"
            # If bin file has 100kb bins, and 'windows' is 100000 (100kb).
            # And continue_n is 50.
            # Then the output window is 50 * 100kb = 5Mb?
            # Original code: bin_start = start * windows
            # start increments by continue_n.
            # So yes, it seems to be creating larger windows from smaller bins.
            
            # Let's trust the logic:
            # start index 0, 1, 2... corresponding to blocks of continue_n
            
            window_starts = np.arange(n_windows) * continue_n * windows
            window_ends = (np.arange(n_windows) + 1) * continue_n * windows - 1
            
            # Z-scores
            # Calculated per chromosome or globally?
            # Original code:
            # short_z = (short_s - np.mean(short_s)) / np.std(short_s)
            # It accumulates `short_s` lists across ALL chromosomes in the loop `while step < length`.
            # Then calculates Z-score on the full list.
            # So Z-score is global.
            
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
                
        logger.info(f"FSC calculation complete. Results written to {output_file}")

    except Exception as e:
        logger.error(f"Fatal error in _calc_fsc: {e}")
        raise typer.Exit(1)


def fsc(
    bedgz_path: Path = typer.Argument(..., help="Folder containing .bed.gz files (should be the output directory from motif.py)"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file (default: data/ChormosomeBins/hg19_window_100kb.bed)"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size (default: 100000)"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number (default: 50)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output folder for results"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of parallel processes (default: 1)")
):
    """
    Calculate fragment size coverage (FSC) features for all .bed.gz files in a folder.
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
        logger.error(f"Output directory not found: {output}")
        raise typer.Exit(1)
    if not output.is_dir():
        logger.error(f"Output path is not a directory: {output}")
        raise typer.Exit(1)

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
    logger.info(f"Calculating FSC for {len(bedgz_files)} files...")
    from concurrent.futures import ProcessPoolExecutor, as_completed
    logger.info(f"Starting parallel FSC calculation using {threads} processes...")
    def run_fsc_file(bedgz_file):
        output_file = output / (bedgz_file.stem.replace('.bed', '') + '.FSC.txt')
        _calc_fsc(str(bedgz_file), str(bin_input), windows, continue_n, str(output_file))
        return str(output_file)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(run_fsc_file, bedgz_file): bedgz_file for bedgz_file in bedgz_files}
        for future in as_completed(futures):
            bedgz_file = futures[future]
            try:
                result = future.result()
                logger.info(f"FSC calculated: {result}")
            except Exception as exc:
                logger.error(f"FSC calculation failed for {bedgz_file}: {exc}")
    logger.info(f"FSC features calculated for {len(bedgz_files)} files.")
