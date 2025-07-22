import typer
from pathlib import Path
from typing import Optional
import logging
import pysam
import pandas as pd
from collections import defaultdict
from functools import partial
from rich.console import Console
from rich.logging import RichHandler
from concurrent.futures import ProcessPoolExecutor, as_completed
import os

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("ocf")


def calc_ocf(bedgz_file: Path, ocr_file: Path, output_dir: Path):
    """
    Calculate OCF for a single .bed.gz file and OCR region file.
    Output is per-region .sync.end files and a summary all.ocf.csv.
    """
    try:
        tbx = pysam.TabixFile(str(bedgz_file))
        regions = pd.read_csv(ocr_file, sep="\t", header=None, names=["chr", "start", "end", "description"])
        leftPOS = defaultdict(partial(defaultdict, int))
        rightPOS = defaultdict(partial(defaultdict, int))
        total = defaultdict(lambda: [0, 0])
        for _, region in regions.iterrows():
            region_Chr, region_Start, region_End, region_Label = (
                region["chr"], region["start"], region["end"], region["description"])
            try:
                fetched_reads = tbx.fetch(region_Chr, region_Start, region_End)
            except ValueError:
                continue
            for row in fetched_reads:
                tmp_row = row.split()
                rstart = int(tmp_row[1])
                rend = int(tmp_row[2])
                if rstart >= region_Start:
                    s = rstart - region_Start
                    leftPOS[region_Label][s] += 1
                    total[region_Label][0] += 1
                if rend <= region_End:
                    e = rend - region_Start + 1
                    rightPOS[region_Label][e] += 1
                    total[region_Label][1] += 1
        Labels = []
        ocf = []
        outputfile = output_dir / 'all.ocf.csv'
        for label in total.keys():
            output = output_dir / f'{label}.sync.end'
            Labels.append(label)
            le = leftPOS[label]
            re = rightPOS[label]
            ts = total[label][0] / 10000 if total[label][0] else 1
            te = total[label][1] / 10000 if total[label][1] else 1
            num = 2000
            with open(output, 'w') as output_write:
                for k in range(num):
                    l = le[k]
                    r = re[k]
                    output_write.write(
                        f"{k - 1000}\t{l}\t{l / ts}\t{r}\t{r / te}\n")
            # OCF calculation
            with open(output, 'r') as o:
                peak = 60
                bin = 10
                trueends = 0
                background = 0
                for line in o.readlines():
                    loc, left, Left, right, Right = line.split()
                    loc = int(loc)
                    if -peak - bin <= loc <= -peak + bin:
                        trueends += float(Right)
                        background += float(Left)
                    elif peak - bin <= loc <= peak + bin:
                        trueends += float(Left)
                        background += float(Right)
                ocf.append(trueends - background)
        import pandas as pd
        ocf_df = pd.DataFrame({"tissue": Labels, "OCF": ocf})
        ocf_df.to_csv(outputfile, sep="\t", index=None)
        logger.info(f"OCF calculation complete for {bedgz_file}. Results in {output_dir}.")
    except Exception as e:
        logger.error(f"Fatal error in calc_ocf: {e}")
        raise typer.Exit(1)


def ocf(
    bedgz_path: Path = typer.Argument(..., help="Folder containing .bed.gz files (should be the output directory from motif.py)"),
    ocr_input: Optional[Path] = typer.Option(None, "--ocr-input", "-r", help="Path to open chromatin region BED file (default: packaged tissue file)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output folder for results"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of parallel processes (default: 1)")
):
    """
    Calculate orientation-aware cfDNA fragmentation (OCF) features for all .bed.gz files in a folder.
    """
    # Input checks
    if not bedgz_path.exists():
        logger.error(f"Input directory not found: {bedgz_path}")
        raise typer.Exit(1)
    if ocr_input and not ocr_input.exists():
        logger.error(f"OCR region BED file not found: {ocr_input}")
        raise typer.Exit(1)
    try:
        output.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Could not create output directory {output}: {e}")
        raise typer.Exit(1)
    # Set default OCR file if not provided
    if ocr_input is None:
        pkg_dir = Path(__file__).parent
        ocr_input = pkg_dir / "data/OpenChromatinRegion/7specificTissue.all.OC.bed"
    bedgz_files = [f for f in Path(bedgz_path).glob("*.bed.gz")]
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    def run_ocf_file(bedgz_file):
        sample_dir = output / bedgz_file.stem.replace('.bed', '')
        sample_dir.mkdir(exist_ok=True)
        calc_ocf(bedgz_file, ocr_input, sample_dir)
        return str(sample_dir)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(run_ocf_file, bedgz_file): bedgz_file for bedgz_file in bedgz_files}
        for future in as_completed(futures):
            bedgz_file = futures[future]
            try:
                result = future.result()
                logger.info(f"OCF calculated: {result}")
            except Exception as exc:
                logger.error(f"OCF calculation failed for {bedgz_file}: {exc}")
    logger.info(f"OCF features calculated for {len(bedgz_files)} files.")
