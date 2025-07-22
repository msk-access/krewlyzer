import typer
from pathlib import Path
import logging
import pysam
import pybedtools
import numpy as np
from collections import defaultdict
import gzip
from rich.console import Console
from rich.logging import RichHandler

from .helpers import max_core, commonError

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("wps")


def _calc_wps(bedgz_input, tsv_input, output_file_pattern, empty=False, protect_input=120, min_size=120, max_size=180):
    """
    Calculate Windowed Protection Score (WPS) for a single .bed.gz file and transcript region file.
    Output is gzipped TSV per region.
    """
    try:
        bedgzfile = str(bedgz_input)
        tbx = pysam.TabixFile(bedgzfile)
        protection = protect_input // 2
        with open(tsv_input, 'r') as infile:
            prefix = "chr"
            valid_chroms = set(map(str, list(range(1, 23)) + ["X"]))
            logger.info(f"input file: {bedgz_input}, {tsv_input}")
            for line in infile:
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 5:
                    continue
                cid, chrom, start, end, strand = parts[:5]
                chrom = chrom.replace("chr", "")
                if chrom not in valid_chroms:
                    continue
                region_start, region_end = int(float(start)), int(float(end))
                if region_start < 1:
                    continue
                pos_range = defaultdict(lambda: [0, 0])
                try:
                    from bx.intervals.intersection import Intersecter, Interval
                    filtered_reads = Intersecter()
                    for row in tbx.fetch(prefix + chrom, region_start - protection, region_end + protection):
                        tmp_row = row.split()
                        rstart = int(tmp_row[1])
                        rend = int(tmp_row[2])
                        lseq = rend - rstart
                        if lseq < min_size or lseq > max_size:
                            continue
                        filtered_reads.add_interval(Interval(rstart, rend))
                        for i in range(rstart, rend):
                            if region_start <= i <= region_end:
                                pos_range[i][0] += 1
                        if region_start <= rstart <= region_end:
                            pos_range[rstart][1] += 1
                        if region_start <= rend <= region_end:
                            pos_range[rend][1] += 1
                except Exception as e:
                    logger.error(f"Error fetching region {chrom}:{region_start}-{region_end}: {e}")
                    continue
                filename = output_file_pattern % cid
                with gzip.open(filename, 'wt') as outfile:
                    cov_sites = 0
                    out_lines = []
                    for pos in range(region_start, region_end + 1):
                        rstart, rend = pos - protection, pos + protection
                        gcount, bcount = 0, 0
                        for read in filtered_reads.find(rstart, rend):
                            if (read.start > rstart) or (read.end < rend):
                                bcount += 1
                            else:
                                gcount += 1
                        cov_count, start_count = pos_range[pos]
                        cov_sites += cov_count
                        out_lines.append(f"{chrom}\t{pos}\t{cov_count}\t{start_count}\t{gcount - bcount}\n")
                    if strand == "-":
                        out_lines = out_lines[::-1]
                    for line in out_lines:
                        outfile.write(line)
                if cov_sites == 0 and not empty:
                    import os
                    os.remove(filename)
        logger.info(f"WPS calculation complete. Results written to pattern: {output_file_pattern}")
    except Exception as e:
        logger.error(f"Fatal error in _calc_wps: {e}")
        raise typer.Exit(1)


def wps(
    bedgz_path: Path = typer.Argument(..., help="Folder containing .bed.gz files (should be the output directory from motif.py)"),
    tsv_input: Path = typer.Option(None, "--tsv-input", "-t", help="Path to transcript/region file (TSV format)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output folder for results"),
    wpstype: str = typer.Option('L', "--wpstype", "-w", help="WPS type: 'L' for long (default), 'S' for short"),
    empty: bool = typer.Option(False, "--empty", help="Keep files of empty blocks (default: False)"),
    threads: int = typer.Option(1, "--threads", "-p", help="Number of threads (default: 1)")
):
    """
    Calculate Windowed Protection Score (WPS) features for all .bed.gz files in a folder.
    """
    # Input checks
    if not bedgz_path.exists():
        logger.error(f"Input directory not found: {bedgz_path}")
        raise typer.Exit(1)
    if tsv_input and not tsv_input.exists():
        logger.error(f"Transcript region file not found: {tsv_input}")
        raise typer.Exit(1)
    try:
        output.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Could not create output directory {output}: {e}")
        raise typer.Exit(1)
    try:
        output.touch()
    except Exception as e:
        logger.error(f"Output directory {output} is not writable: {e}")
        raise typer.Exit(1)
    try:
        bedgz_files = list(Path(bedgz_path).glob("*.bed.gz"))
        if not bedgz_files:
            logger.error("No .bed.gz files found in the specified folder.")
            raise typer.Exit(1)
        if tsv_input is None:
            # Default to package data transcriptAnno-hg19-1kb.tsv
            tsv_input = Path(__file__).parent.parent / "data" / "TranscriptAnno" / "transcriptAnno-hg19-1kb.tsv"
            logger.info(f"No tsv_input specified. Using default: {tsv_input}")
        if not tsv_input.exists():
            logger.error(f"Transcript/region file does not exist: {tsv_input}")
            raise typer.Exit(1)
        if wpstype == 'L':
            protect_input = 120
            min_size = 120
            max_size = 180
        else:
            protect_input = 16
            min_size = 35
            max_size = 80
        output.mkdir(parents=True, exist_ok=True)
        logger.info(f"Calculating WPS for {len(bedgz_files)} files...")
        from concurrent.futures import ProcessPoolExecutor, as_completed
        import traceback
        def wps_task(bedgz_file):
            try:
                output_file_pattern = str(output / (bedgz_file.stem.replace('.bed', '') + ".%s.WPS.tsv.gz"))
                _calc_wps(
                    bedgz_input=str(bedgz_file),
                    tsv_input=str(tsv_input),
                    output_file_pattern=output_file_pattern,
                    empty=empty,
                    protect_input=protect_input,
                    min_size=min_size,
                    max_size=max_size
                )
                return None
            except Exception as exc:
                return traceback.format_exc()
        n_procs = max_core(threads) if threads else 1
        logger.info(f"Calculating WPS for {len(bedgz_files)} files using {n_procs} processes...")
        with ProcessPoolExecutor(max_workers=n_procs) as executor:
            futures = {executor.submit(wps_task, bedgz_file): bedgz_file for bedgz_file in bedgz_files}
            for future in as_completed(futures):
                exc = future.result()
                if exc:
                    logger.error(f"WPS calculation failed for {futures[future]}:\n{exc}")
        logger.info(f"WPS features calculated for {len(bedgz_files)} files.")
    except Exception as e:
        logger.error(f"Fatal error in wps CLI: {e}")
        raise typer.Exit(1)
