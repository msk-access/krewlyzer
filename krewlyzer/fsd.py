import typer
from pathlib import Path
import logging
import pysam
import pybedtools
import numpy as np
from rich.console import Console
from rich.logging import RichHandler

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("fsd")

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
    # Input checks
    if not bedgz_path.exists():
        logger.error(f"Input directory not found: {bedgz_path}")
        raise typer.Exit(1)
    if not arms_file.exists():
        logger.error(f"Arms/region file not found: {arms_file}")
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
