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


def _calc_fsr(bedgz_input, bin_input, windows, continue_n, output_file):
    """
    Internal: Calculate fragment size ratio (FSR) for a single .bed.gz file.
    Writes region-based ratios for short, intermediate, and long fragments.
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
                if totals == 0:
                    shorts_data.append(0)
                    intermediates_data.append(0)
                    longs_data.append(0)
                else:
                    shorts_data.append(shorts / totals)
                    intermediates_data.append(intermediates / totals)
                    longs_data.append(longs / totals)
        start = 0
        step = 0
        try:
            with open(output_file, 'w') as fsrfile:
                fsrfile.write("region\tshort-ratio\titermediate-ratio\tlong-ratio\n")
                while step < length:
                    num = chrom.count(chrom[step])
                    continues_bin = num // continue_n
                    last_bin = num % continue_n
                    for _ in range(continues_bin):
                        bin_start = start * windows
                        bin_end = (start + continue_n) * windows - 1
                        combine_shorts = shorts_data[step: step + continue_n]
                        combine_intermediates = intermediates_data[step: step + continue_n]
                        combine_longs = longs_data[step: step + continue_n]
                        tmp_array = np.zeros(3)
                        tmp_array[0] = np.mean(combine_shorts)
                        tmp_array[1] = np.mean(combine_intermediates)
                        tmp_array[2] = np.mean(combine_longs)
                        region = f"{chrom[step]}:{bin_start}-{bin_end}"
                        temp_str = f"{region}\t" + "\t".join(map(str, tmp_array)) + "\n"
                        fsrfile.write(temp_str)
                        step += continue_n
                        start += continue_n
                    if last_bin != 0:
                        step += last_bin
                        start = 0
        except Exception as e:
            logger.error(f"Error writing FSR output file: {e}")
            raise typer.Exit(1)
        logger.info(f"FSR calculation complete. Results written to {output_file}")
    except Exception as e:
        logger.error(f"Fatal error in _calc_fsr: {e}")
        raise typer.Exit(1)

def _calc_fsd(bedgz_input, arms_file, output_file):
    """
    Internal: Calculate fragment size distribution (FSD) for a single .bed.gz file.
    Writes region-based fragment size distributions in 5bp bins from 65-399bp.
    """
    try:
        logger.info(f"input file: {bedgz_input}, {arms_file}")
        try:
            inputbed = pysam.Tabixfile(filename=bedgz_input, mode="r")
        except Exception as e:
            logger.error(f"Could not open {bedgz_input} as Tabix file: {e}")
            raise typer.Exit(1)
        try:
            bins = pybedtools.BedTool(arms_file)
        except Exception as e:
            logger.error(f"Could not load bins from {arms_file}: {e}")
            raise typer.Exit(1)
        length = len(bins)
        interval_data = []
        region = []
        logger.info(f"output file: {output_file}")
        for idx in range(length):
            bin = bins[idx]
            region.append(f"{bin.chrom}:{bin.start}-{bin.end}")
            try:
                inputbed.fetch(bin.chrom, bin.start, bin.end)
            except ValueError:
                interval_data.append([0] * 67)
                continue
            except Exception as e:
                logger.error(f"Error fetching bin {bin}: {e}")
                raise typer.Exit(1)
            else:
                bin_data = []
                try:
                    for read in inputbed.fetch(bin.chrom, bin.start, bin.end):
                        bin_data.append(int(read.split("\t")[2]) - int(read.split("\t")[1]))
                    count = np.bincount(bin_data, minlength=401)
                    step_size = 5
                    start_bin = 65
                    end_bin = 400
                    bin_len = int((end_bin - start_bin) / step_size)
                    temp_bin = []
                    for bin_id in range(bin_len):
                        temp_bin.append(np.sum(count[(start_bin + step_size * bin_id):(start_bin + step_size * (bin_id + 1))]))
                    interval_data.append(temp_bin)
                except Exception as e:
                    logger.error(f"Error processing reads in bin {bin}: {e}")
                    interval_data.append([0] * 67)
                    continue
        try:
            with open(output_file, 'w') as fsdfile:
                sbin = np.arange(65, 400, 5)
                head_str = 'region' + '\t' + '\t'.join([f"{s}-{s+4}" for s in sbin]) + '\n'
                fsdfile.write(head_str)
                for i in range(length):
                    arms = interval_data[i]
                    score = np.zeros(67)
                    if np.sum(arms) != 0:
                        score = np.array(arms) / np.sum(arms)
                    temp_str = region[i] + '\t' + '\t'.join(map(str, score)) + '\n'
                    fsdfile.write(temp_str)
        except Exception as e:
            logger.error(f"Error writing FSD output file: {e}")
            raise typer.Exit(1)
        logger.info(f"FSD calculation complete. Results written to {output_file}")
    except Exception as e:
        logger.error(f"Fatal error in _calc_fsd: {e}")
        raise typer.Exit(1)


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
    if not output.is_writable():
        logger.error(f"Output directory is not writable: {output}")
        raise typer.Exit(1)

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


def fsr(
    bedgz_path: Path = typer.Argument(..., help="Folder containing .bed.gz files (should be the output directory from motif.py)"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file (default: data/ChormosomeBins/hg19_window_100kb.bed)"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size (default: 100000)"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number (default: 50)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output folder for results"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads (default: 1)")
):
    """
    Calculate fragment size ratio (FSR) features for all .bed.gz files in a folder.
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
        bin_input = Path(__file__).parent / "data" / "ChormosomeBins" / "hg19_window_100kb.bed"
        logger.info(f"No bin_input specified. Using default: {bin_input}")
    if not bin_input.exists():
        logger.error(f"Bin input file does not exist: {bin_input}")
        raise typer.Exit(1)
    logger.info(f"Calculating FSR for {len(bedgz_files)} files...")
    from concurrent.futures import ProcessPoolExecutor, as_completed
    logger.info(f"Starting parallel FSR calculation using {threads} processes...")
    def run_fsr_file(bedgz_file):
        output_file = output / (bedgz_file.stem.replace('.bed', '') + '.FSR.txt')
        _calc_fsr(str(bedgz_file), str(bin_input), windows, continue_n, str(output_file))
        return str(output_file)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(run_fsr_file, bedgz_file): bedgz_file for bedgz_file in bedgz_files}
        for future in as_completed(futures):
            bedgz_file = futures[future]
            try:
                result = future.result()
                logger.info(f"FSR calculated: {result}")
            except Exception as exc:
                logger.error(f"FSR calculation failed for {bedgz_file}: {exc}")
    logger.info(f"FSR features calculated for {len(bedgz_files)} files.")


def fsd(
    bedgz_path: Path = typer.Argument(..., help="Folder containing .bed.gz files (should be the output directory from motif.py)"),
    arms_file: Path = typer.Option(..., "--arms-file", "-a", help="Path to arms/region file (BED format)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output folder for results"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads (default: 1)")
):
    """
    Calculate fragment size distribution (FSD) features for all .bed.gz files in a folder.
    The input folder should be the output directory produced by motif.py, containing the .bed.gz files.
    Output files are written to the output directory, one per .bed.gz file.
    """
    if not output.exists():
        output.mkdir(parents=True, exist_ok=True)
    bedgz_files = [f for f in bedgz_path.iterdir() if f.suffixes == ['.bed', '.gz']]
    if not bedgz_files:
        logger.error("No .bed.gz files found in the specified folder.")
        raise typer.Exit(1)
    if not arms_file.exists():
        logger.error(f"Arms/region file does not exist: {arms_file}")
        raise typer.Exit(1)
    logger.info(f"Calculating FSD for {len(bedgz_files)} files...")
    from concurrent.futures import ProcessPoolExecutor, as_completed
    logger.info(f"Starting parallel FSD calculation using {threads} processes...")
    def run_fsd_file(bedgz_file):
        output_file = output / (bedgz_file.stem.replace('.bed', '') + '.FSD.txt')
        _calc_fsd(str(bedgz_file), str(arms_file), str(output_file))
        return str(output_file)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(run_fsd_file, bedgz_file): bedgz_file for bedgz_file in bedgz_files}
        for future in as_completed(futures):
            bedgz_file = futures[future]
            try:
                result = future.result()
                logger.info(f"FSD calculated: {result}")
            except Exception as exc:
                logger.error(f"FSD calculation failed for {bedgz_file}: {exc}")
    logger.info(f"FSD features calculated for {len(bedgz_files)} files.")
