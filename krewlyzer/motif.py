# motif.py: Extracts motif-based features from BAM files
import typer
from pathlib import Path
from typing import Optional
import os
import pysam
import numpy as np
import pandas as pd
import math
import pybedtools
from collections import defaultdict
from .helpers import (
    reverse_seq,
    get_End_motif,
    get_Breakpoint_motif,
    GCcontent,
    read_pair_generator,
    maxCore,
    rmEndString,
    calc_MDS
)
from rich.progress import Progress
from rich.console import Console
from rich.logging import RichHandler
import logging

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("motif")


def motif(
    bam_path: Path = typer.Argument(..., help="Path to input BAM file or directory of BAM files (GRCh37 aligned)"),
    genome_reference: Path = typer.Option(..., '-g', help="Path to genome reference file (GRCh37/hg19)"),
    output: Path = typer.Option(..., '-o', help="Output directory"),
    blacklist: Optional[Path] = typer.Option(None, '-b', help="Path to blacklist regions file"),
    map_quality: int = typer.Option(20, '-m', help="Minimum mapping quality"),
    min_length: int = typer.Option(65, '--minlen', help="Minimum fragment length"),
    max_length: int = typer.Option(400, '--maxlen', help="Maximum fragment length"),
    kmer: int = typer.Option(3, '-k', help="K-mer size for motif extraction"),
    chromosomes: Optional[str] = typer.Option(None, '--chromosomes', help="Comma-separated list of chromosomes to process"),
    verbose: bool = typer.Option(False, '--verbose', help="Enable verbose logging"),
    threads: int = typer.Option(1, '--threads', help="Number of parallel processes (default: 1)")
):
    """
    Extract motif-based features from BAM files.
    """
    # Input checks
    if not bam_path.exists():
        logger.error(f"Input BAM file or directory not found: {bam_path}")
        raise typer.Exit(1)
    if not genome_reference.exists() or not genome_reference.is_file():
        logger.error(f"Reference genome file not found: {genome_reference}")
        raise typer.Exit(1)
    try:
        output.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Could not create output directory {output}: {e}")
        raise typer.Exit(1)
    """
    Extracts end motif, breakpoint motif, and Motif-Diversity Score (MDS) from one or more BAM files.
    If a directory is provided, all BAM files in the directory will be processed in parallel using multiple processes.
    Output files are written to the output directory, with EDM, BPM, and MDS subfolders.
    """
    import concurrent.futures
    if verbose:
        logger.setLevel(logging.DEBUG)
    logger.info(f"Reference genome: {genome_reference}")
    logger.info(f"Output directory: {output}")
    if bam_path.is_dir():
        bam_files = sorted([f for f in bam_path.iterdir() if f.suffix == '.bam'])
        if not bam_files:
            logger.error(f"No BAM files found in directory: {bam_path}")
            raise typer.Exit(1)
        logger.info(f"Processing {len(bam_files)} BAM files in parallel using {threads} processes...")
        def run_motif_for_bam(bam_file):
            logger.info(f"Processing BAM: {bam_file}")
            motif_process(
                str(bam_file),
                str(blacklist) if blacklist else None,
                str(output / (bam_file.stem + '.bed')),
                str(genome_reference),
                chromosomes.split(',') if chromosomes else None,
                map_quality,
                kmer,
                fragFilter=True,
                minLen=min_length,
                maxLen=max_length
            )
            return str(bam_file)
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(run_motif_for_bam, bam_file): bam_file for bam_file in bam_files}
            for future in concurrent.futures.as_completed(futures):
                bam_file = futures[future]
                try:
                    result = future.result()
                    logger.info(f"Motif extraction complete for: {result}")
                except Exception as exc:
                    logger.error(f"Motif extraction failed for {bam_file}: {exc}")
        logger.info(f"All BAM files processed.")
    else:
        logger.info(f"Processing BAM: {bam_path}")
        motif_process(
            str(bam_path),
            str(blacklist) if blacklist else None,
            str(output / (bam_path.stem + '.bed')),
            str(genome_reference),
            chromosomes.split(',') if chromosomes else None,
            map_quality,
            kmer,
            fragFilter=True,
            minLen=min_length,
            maxLen=max_length
        )
        logger.info("End motif, Breakpoint motif, and MDS extraction complete.")


def motif_process(
    bamInput,
    blacklistInput,
    bedOutput,
    genome_reference,
    CHR,
    mapQuality,
    k_mer,
    fragFilter=False,
    minLen=None,
    maxLen=None
):
    """
    Main motif feature extraction process with rich logging and consistent CLI output.
    """
    from rich.table import Table
    from rich.panel import Panel
    from rich.text import Text
    bedOutput_path = os.path.abspath(bedOutput)
    EDM_output_path = os.path.join(os.path.dirname(bedOutput_path), 'EDM')
    BPM_output_path = os.path.join(os.path.dirname(bedOutput_path), 'BPM')
    MDS_output_path = os.path.join(os.path.dirname(bedOutput_path), 'MDS')
    try:
        os.makedirs(EDM_output_path, exist_ok=True)
        os.makedirs(BPM_output_path, exist_ok=True)
        os.makedirs(MDS_output_path, exist_ok=True)
    except Exception as e:
        logger.error(f"Failed to create output directories: {e}")
        raise typer.Exit(1)
    bases = ['A', 'C', 'T', 'G']
    End_motif = {''.join(i): 0 for i in itertools.product(bases, repeat=k_mer)}
    Breakpoint_motif = {''.join(i): 0 for i in itertools.product(bases, repeat=k_mer)}
    try:
        bamfile = pysam.AlignmentFile(bamInput, 'rb')
    except Exception as e:
        logger.error(f"Failed to open BAM file: {e}")
        raise typer.Exit(1)
    try:
        genome = pysam.FastaFile(genome_reference)
    except Exception as e:
        logger.error(f"Failed to open genome FASTA: {e}")
        raise typer.Exit(1)
    temp_bed = bedOutput + '.tmp'
    try:
        bedWrite = open(temp_bed, 'w')
    except Exception as e:
        logger.error(f"Failed to open temp BED for writing: {e}")
        raise typer.Exit(1)
    chroms = CHR if CHR else list(bamfile.references)
    logger.info("Extracting motif features from BAM file...")
    total_pairs = bamfile.mapped // 2 if bamfile.mapped else 1000000
    motif_errors = 0
    with Progress(console=console, transient=True) as progress:
        task = progress.add_task("Processing fragments", total=total_pairs)
        for idx, pair in enumerate(read_pair_generator(bamfile)):
            try:
                read1, read2 = pair
                if read1.mapping_quality < mapQuality or read2.mapping_quality < mapQuality or read1.reference_name not in chroms:
                    continue
                read1Start = read1.reference_start
                read1End = read1.reference_end
                read2Start = read2.reference_start
                read2End = read2.reference_end
                if not read1.is_reverse:
                    rstart = read1Start
                    rend = read2End
                    forward_end5 = read1.query_sequence[:k_mer].upper()
                    forward_end3 = read2.query_sequence[-k_mer:].upper()
                else:
                    rstart = read2Start
                    rend = read1End
                    forward_end5 = read2.query_sequence[:k_mer].upper()
                    forward_end3 = read1.query_sequence[-k_mer:].upper()
                if (rstart < 0) or (rend < 0) or (rstart >= rend):
                    continue
                if fragFilter:
                    readLen = rend - rstart
                    if (minLen and readLen < minLen) or (maxLen and readLen > maxLen):
                        continue
                gc = GCcontent(genome.fetch(read1.reference_name, rstart, rend))
                bedWrite.write(f"{read1.reference_name}\t{rstart+1}\t{rend+1}\t{gc}\n")
                End_motif = get_End_motif(End_motif, forward_end3, forward_end3)
                pos = math.ceil(k_mer / 2)
                try:
                    if k_mer % 2 == 0:
                        ref_seq1 = genome.fetch(read1.reference_name, rstart - pos, rstart).upper()
                        ref_seq2 = genome.fetch(read2.reference_name, rend, rend + pos).upper()
                        Breakpoint_motif = get_Breakpoint_motif(Breakpoint_motif, ref_seq1 + forward_end5[:pos], forward_end3[-pos:] + ref_seq2)
                    else:
                        ref_seq1 = genome.fetch(read1.reference_name, rstart - pos + 1, rstart).upper()
                        ref_seq2 = genome.fetch(read2.reference_name, rend, rend + pos - 1).upper()
                        Breakpoint_motif = get_Breakpoint_motif(Breakpoint_motif, ref_seq1 + forward_end5[:pos], forward_end3[-pos:] + ref_seq2)
                except Exception as e:
                    motif_errors += 1
                    logger.warning(f"Motif extraction failed for fragment at {read1.reference_name}:{rstart}-{rend}: {e}")
                    continue
                if idx % 10000 == 0:
                    progress.update(task, advance=10000)
            except Exception as e:
                motif_errors += 1
                logger.error(f"Unexpected error during fragment processing: {e}")
                continue
        progress.update(task, completed=total_pairs)
    bedWrite.close()
    logger.info("Filtering and sorting fragments with blacklist (if provided)...")
    try:
        bedData = pybedtools.BedTool(temp_bed)
        if blacklistInput:
            black_reigon = pybedtools.BedTool(blacklistInput)
            bedData = bedData.subtract(black_reigon, A=True)
        bedData.sort(output=bedOutput)
        os.remove(temp_bed)
    except Exception as e:
        logger.error(f"Error during BED filtering/sorting: {e}")
        raise typer.Exit(1)
    # Write EndMotif
    edm_file = os.path.join(EDM_output_path, Path(bedOutput).stem + '.EndMotif')
    logger.info(f"Writing End Motif frequencies to {edm_file}")
    try:
        with open(edm_file, 'w') as f:
            total = sum(End_motif.values())
            for k, v in End_motif.items():
                f.write(f"{k}\t{v/total if total else 0}\n")
    except Exception as e:
        logger.error(f"Failed to write End Motif output: {e}")
        raise typer.Exit(1)
    # Write BreakPointMotif
    bpm_file = os.path.join(BPM_output_path, Path(bedOutput).stem + '.BreakPointMotif')
    logger.info(f"Writing Breakpoint Motif frequencies to {bpm_file}")
    try:
        with open(bpm_file, 'w') as f:
            total = sum(Breakpoint_motif.values())
            for k, v in Breakpoint_motif.items():
                f.write(f"{k}\t{v/total if total else 0}\n")
    except Exception as e:
        logger.error(f"Failed to write Breakpoint Motif output: {e}")
        raise typer.Exit(1)
    # Write MDS (Motif Diversity Score)
    mds_file = os.path.join(MDS_output_path, Path(bedOutput).stem + '.MDS')
    logger.info(f"Writing Motif Diversity Score to {mds_file}")
    try:
        df = pd.read_csv(edm_file, sep='\t', header=None, names=['motif', 'frequency'])
        freq = df['frequency'].values
        mds = -np.sum(freq * np.log2(freq + 1e-12)) / np.log2(len(freq))
        with open(mds_file, 'w') as f:
            f.write(f"{mds}\n")
    except Exception as e:
        logger.error(f"Failed to write MDS output: {e}")
        raise typer.Exit(1)
    # Print summary
    summary_table = Table(title="Motif Extraction Summary", show_header=True, header_style="bold magenta")
    summary_table.add_column("Output Type", style="bold")
    summary_table.add_column("File Path")
    summary_table.add_row("End Motif (EDM)", edm_file)
    summary_table.add_row("Breakpoint Motif (BPM)", bpm_file)
    summary_table.add_row("Motif Diversity Score (MDS)", mds_file)
    console.print(Panel(summary_table, title="[green]Extraction Complete", subtitle=f"Motif errors: {motif_errors}", expand=False))
    logger.info("Motif feature extraction complete.")
