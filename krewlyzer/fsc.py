import typer
from pathlib import Path
from typing import Optional
import logging
import sys

import pysam
import pybedtools
import numpy as np
import pandas as pd
from skmisc.loess import loess
from rich.console import Console
from rich.logging import RichHandler

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("fsc")


from .helpers import gc_correct


def _calc_fsc(bedgz_input, bin_input, windows, continue_n, output_file):
    """
    Internal: Calculate fragment size coverage (FSC) for a single .bed.gz file.
    Handles errors and logs all steps. Raises typer.Exit(1) on fatal errors.
    """
    try:
        logger.info(f"input file: {bedgz_input}, {bin_input}")
        try:
            inputbed = pysam.Tabixfile(filename=bedgz_input, mode="r")
        except Exception as e:
            logger.error(f"Could not open {bedgz_input} as Tabix file: {e}")
            raise typer.Exit(1)
        try:
            bins = pybedtools.BedTool(bin_input)
        except Exception as e:
            logger.error(f"Could not load bins from {bin_input}: {e}")
            raise typer.Exit(1)
        length = len(bins)
        shorts_data, intermediates_data, longs_data, totals_data, bingc = [], [], [], [], []
        chrom = []
        logger.info(f"output file: {output_file}")
        for idx in range(length):
            bin = bins[idx]
            try:
                chrom.append(bin.chrom)
                inputbed.fetch(bin.chrom, bin.start, bin.end)
            except ValueError:
                bingc.append(np.nan)
                shorts_data.append(0)
                intermediates_data.append(0)
                longs_data.append(0)
                totals_data.append(0)
            except Exception as e:
                logger.error(f"Error fetching bin {bin}: {e}")
                raise typer.Exit(1)
            else:
                bin_data = []
                gc = []
                try:
                    for read in inputbed.fetch(bin.chrom, bin.start, bin.end):
                        bin_data.append(int(read.split("\t")[2]) - int(read.split("\t")[1]))
                        if 65 <= int(read.split("\t")[2]) - int(read.split("\t")[1]) <= 400:
                            gc.append(float(read.split("\t")[3]))
                    count = np.bincount(bin_data, minlength=401)
                except Exception as e:
                    logger.error(f"Error processing reads in bin {bin}: {e}")
                    raise typer.Exit(1)
                if len(gc) == 0:
                    bingc.append(np.nan)
                else:
                    bingc.append(np.mean(gc))
                shorts = sum(count[65:150])
                intermediates = sum(count[151:260])
                longs = sum(count[261:400])
                totals = sum(count[65:400])
                shorts_data.append(shorts)
                intermediates_data.append(intermediates)
                longs_data.append(longs)
                totals_data.append(totals)
        try:
            correct_shorts = gc_correct(shorts_data, bingc)
            correct_intermediates = gc_correct(intermediates_data, bingc)
            correct_longs = gc_correct(longs_data, bingc)
            correct_totals = gc_correct(totals_data, bingc)
        except Exception as e:
            logger.error(f"GC correction failed: {e}")
            raise typer.Exit(1)
        start = 0
        step = 0
        short_s, intermediate_s, long_s, total_s = [], [], [], []
        region = []
        try:
            with open(output_file, 'w') as fscfile:
                fscfile.write(
                    "region\tshort-fragment-zscore\titermediate-fragment-zscore\tlong-fragment-zscore\ttotal-fragment-zscore\n"
                )
                while step < length:
                    num = chrom.count(chrom[step])
                    continues_bin = num // continue_n
                    last_bin = num % continue_n
                    for _ in range(continues_bin):
                        bin_start = start * windows
                        bin_end = (start + continue_n) * windows - 1
                        combine_shorts = correct_shorts[step: step + continue_n]
                        combine_intermediates = correct_intermediates[step: step + continue_n]
                        combine_longs = correct_longs[step: step + continue_n]
                        combine_totals = correct_totals[step: step + continue_n]
                        short_s.append(np.sum(combine_shorts))
                        intermediate_s.append(np.sum(combine_intermediates))
                        long_s.append(np.sum(combine_longs))
                        total_s.append(np.sum(combine_totals))
                        region.append(f"{chrom[step]}:{bin_start}-{bin_end}")
                        step += continue_n
                        start += continue_n
                    if last_bin != 0:
                        step += last_bin
                        start = 0
                try:
                    short_z = (np.array(short_s) - np.mean(short_s)) / np.std(short_s)
                    intermediate_z = (np.array(intermediate_s) - np.mean(intermediate_s)) / np.std(intermediate_s)
                    long_z = (np.array(long_s) - np.mean(long_s)) / np.std(long_s)
                    total_z = (np.array(total_s) - np.mean(total_s)) / np.std(total_s)
                except Exception as e:
                    logger.error(f"Error calculating z-scores: {e}")
                    raise typer.Exit(1)
                for j in range(len(region)):
                    temp_str = f"{region[j]}\t{short_z[j]}\t{intermediate_z[j]}\t{long_z[j]}\t{total_z[j]}\n"
                    fscfile.write(temp_str)
        except Exception as e:
            logger.error(f"Error writing FSC output file: {e}")
            raise typer.Exit(1)
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
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads (default: 1)")
):
    """
    Calculate fragment size coverage (FSC) features for all .bed.gz files in a folder.
    The input folder should be the output directory produced by motif.py, containing the .bed.gz files.
    Output files are written to the output directory, one per .bed.gz file.
    """
    if not output.exists():
        output.mkdir(parents=True, exist_ok=True)
    bedgz_files = [f for f in bedgz_path.iterdir() if f.suffixes == ['.bed', '.gz']]
    if not bedgz_files:
        logger.error("No .bed.gz files found in the specified folder.")
        raise typer.Exit(1)
    if bin_input is None:
        # Use package-relative default
        bin_input = Path(__file__).parent / "data" / "ChormosomeBins" / "hg19_window_100kb.bed"
        logger.info(f"No bin_input specified. Using default: {bin_input}")
    if not bin_input.exists():
        logger.error(f"Bin input file does not exist: {bin_input}")
        raise typer.Exit(1)
    logger.info(f"Calculating FSC for {len(bedgz_files)} files...")
    for bedgz_file in bedgz_files:
        output_file = output / (bedgz_file.stem.replace('.bed', '') + '.FSC.txt')
        _calc_fsc(str(bedgz_file), str(bin_input), windows, continue_n, str(output_file))
    logger.info(f"FSC features calculated for {len(bedgz_files)} files.")
