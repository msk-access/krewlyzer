import typer
from pathlib import Path
from typing import Optional
import logging
import pysam
import pybedtools
import numpy as np
from rich.console import Console
from rich.logging import RichHandler
from concurrent.futures import ProcessPoolExecutor, as_completed
import os

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("uxm")

def calc_uxm(
    bam_file: Path,
    mark_file: Path,
    output_file: Path,
    map_quality: int,
    min_cpg: int,
    methy_threshold: float,
    unmethy_threshold: float,
    pe_type: str = "PE"
):
    """
    Calculate UXM fragment-level methylation for a single BAM and marker file.
    Output is a .UXM.tsv file with region, U, X, M proportions.
    """
    try:
        bai = str(bam_file) + ".bai"
        if not os.path.exists(bai):
            pysam.sort("-o", str(bam_file), str(bam_file))
            pysam.index(str(bam_file))
            logger.warning(f"Index file {bai} did not exist. Sorted and indexed BAM.")
        input_file = pysam.AlignmentFile(str(bam_file))
        marks = pybedtools.BedTool(str(mark_file))
        res = []
        for mark in marks:
            region = f"{mark.chrom}:{mark.start}-{mark.end}"
            try:
                input_file.fetch(mark.chrom, mark.start, mark.end)
            except ValueError:
                res.append(f"{region}\t0\t0\t0")
                continue
            Ufragment = 0
            Xfragment = 0
            Mfragment = 0
            if pe_type == "PE":
                from krewlyzer.helpers import read_pair_generator
                region_string = f"{mark.chrom}:{mark.start}-{mark.end}"
                for read1, read2 in read_pair_generator(input_file, region_string):
                    if read1 is None or read2 is None:
                        continue
                    if read1.mapping_quality < map_quality or read2.mapping_quality < map_quality:
                        continue
                    try:
                        m1 = read1.get_tag("XM")
                        m2 = read2.get_tag("XM")
                    except KeyError:
                        continue
                    read1Start = read1.reference_start
                    read1End = read1.reference_end
                    read2Start = read2.reference_start
                    read2End = read2.reference_end
                    # cfDNAFE logic for overlap
                    if not read1.is_reverse:  # read1 is forward, read2 is reverse
                        if read2Start < read1End:
                            overlap = read1End - read2Start
                            num_methylated = m1.count("Z") + m2[overlap:].count("Z")
                            num_unmethylated = m1.count("z") + m2[overlap:].count("z")
                        else:
                            num_methylated = m1.count("Z") + m2.count("Z")
                            num_unmethylated = m1.count("z") + m2.count("z")
                    else:  # read1 is reverse, read2 is forward
                        if read1Start < read2End:
                            overlap = read2End - read1Start
                            num_methylated = m2.count("Z") + m1[overlap:].count("Z")
                            num_unmethylated = m2.count("z") + m1[overlap:].count("z")
                        else:
                            num_methylated = m1.count("Z") + m2.count("Z")
                            num_unmethylated = m1.count("z") + m2.count("z")
                    if num_methylated + num_unmethylated < min_cpg:
                        continue
                    ratio = num_methylated / (num_methylated + num_unmethylated)
                    if ratio >= methy_threshold:
                        Mfragment += 1
                    elif ratio <= unmethy_threshold:
                        Ufragment += 1
                    else:
                        Xfragment += 1
            elif pe_type == "SE":
                for read in input_file.fetch(mark.chrom, mark.start, mark.end):
                    if read.mapping_quality < map_quality:
                        continue
                    try:
                        m = read.get_tag("XM")
                    except KeyError:
                        continue
                    num_methylated = m.count("Z")
                    num_unmethylated = m.count("z")
                    if num_methylated + num_unmethylated < min_cpg:
                        continue
                    ratio = num_methylated / (num_methylated + num_unmethylated)
                    if ratio >= methy_threshold:
                        Mfragment += 1
                    elif ratio <= unmethy_threshold:
                        Ufragment += 1
                    else:
                        Xfragment += 1
            else:
                logger.error("type must be SE or PE")
                raise typer.Exit(1)
            total = Mfragment + Ufragment + Xfragment
            if total == 0:
                res.append(f"{region}\t0\t0\t0")
            else:
                tmp_array = np.zeros(3)
                tmp_array[0] = Ufragment / total
                tmp_array[1] = Xfragment / total
                tmp_array[2] = Mfragment / total
                res.append(f"{region}\t" + "\t".join(map(str, tmp_array)))
        with open(output_file, 'w') as f:
            f.write('region\tU\tX\tM\n')
            for i in res:
                f.write(i + '\n')
        logger.info(f"UXM calculation complete for {bam_file}. Results in {output_file}.")
    except Exception as e:
        logger.error(f"Fatal error in calc_uxm: {e}")
        raise typer.Exit(1)

def uxm(
    bam_path: Path = typer.Argument(..., help="Folder containing .bam files for UXM calculation."),
    mark_input: Optional[Path] = typer.Option(None, "--mark-input", "-m", help="Marker BED file (default: packaged atlas)", show_default=False),
    output: Path = typer.Option(..., "--output", "-o", help="Output folder for results"),
    map_quality: int = typer.Option(30, "--map-quality", "-q", help="Minimum mapping quality"),
    min_cpg: int = typer.Option(4, "--min-cpg", "-c", help="Minimum CpG count per fragment"),
    methy_threshold: float = typer.Option(0.75, "--methy-threshold", "-tM", help="Methylation threshold for M fragments"),
    unmethy_threshold: float = typer.Option(0.25, "--unmethy-threshold", "-tU", help="Unmethylation threshold for U fragments"),
    pe_type: str = typer.Option("SE", "--type", help="Fragment type: SE or PE (default: SE)"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of parallel processes (default: 1)")
):
    """
    Calculate fragment-level methylation (UXM) features for all BAM files in a folder.
    """
    # Input checks
    if not bam_path.exists():
        logger.error(f"Input BAM directory not found: {bam_path}")
        raise typer.Exit(1)
    if mark_input and not mark_input.exists():
        logger.error(f"Marker BED file not found: {mark_input}")
        raise typer.Exit(1)
    try:
        output.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Could not create output directory {output}: {e}")
        raise typer.Exit(1)
    if mark_input is None:
        pkg_dir = Path(__file__).parent
        mark_input = pkg_dir / "data/MethMark/Atlas.U25.l4.hg19.bed"
    bam_files = [f for f in Path(bam_path).glob("*.bam")]
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    def run_uxm_file(bam_file):
        sample_prefix = bam_file.stem.replace('.bam', '')
        output_file = output / f"{sample_prefix}.UXM.tsv"
        calc_uxm(
            bam_file,
            mark_input,
            output_file,
            map_quality,
            min_cpg,
            methy_threshold,
            unmethy_threshold,
            pe_type
        )
        return str(output_file)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(run_uxm_file, bam_file): bam_file for bam_file in bam_files}
        for future in as_completed(futures):
            bam_file = futures[future]
            try:
                result = future.result()
                logger.info(f"UXM calculated: {result}")
            except Exception as exc:
                logger.error(f"UXM calculation failed for {bam_file}: {exc}")
    logger.info(f"UXM features calculated for {len(bam_files)} files.")
