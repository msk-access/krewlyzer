"""
Fragment-level Methylation analysis (UXM) calculation.

Analyzes CpG methylation patterns in cell-free DNA fragments using bisulfite BAM.
UXM (Unmethylated-to-methylated Xition) scores measure epigenetic tissue signatures.

Uses Rust backend for accelerated computation.

Input: Bisulfite-converted BAM file (from WGBS or capture panel)
Output: {sample}.UXM.tsv - Fragment-level methylation scores at marker regions

Markers: Genomic regions where methylation status differs between tissues.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

from .core.output_utils import (
    read_exact_table,
    write_table,
)  # Unified I/O for Parquet support

console = Console(stderr=True)
# Guard: only configure root logging if not already configured (e.g. by pytest or parent app).
# This prevents double-logging when uxm is imported as a library rather than called as a CLI.
if not logging.root.handlers:
    logging.basicConfig(
        level="INFO",
        handlers=[RichHandler(console=console, show_time=True, show_path=False)],
        format="%(message)s",
    )
logger = logging.getLogger("uxm")

# Rust backend is required
from krewlyzer import _core

# ---------------------------------------------------------------------------
# Fragment methylation classification thresholds (Loyfer et al., Nature 2022)
#
# A fragment is classified by the fraction of its CpGs that are methylated:
#     ratio >= METHY_THRESHOLD    -> M (methylated)
#     ratio <= UNMETHY_THRESHOLD  -> U (unmethylated)
#     otherwise                   -> X (mixed / intermediate)
#
# These MUST straddle a gap. Passing an equal value for both (the historical
# 0.5 / 0.5) makes the X branch unreachable because the Rust backend tests
# `>= methy` first, so every fragment collapses into M or U and the published
# X column is identically zero.
# ---------------------------------------------------------------------------
METHY_THRESHOLD = 0.75
UNMETHY_THRESHOLD = 0.25


def uxm_stage(
    bam_input: Path,
    output: Path,
    sample_name: Optional[str] = None,
    mark_input: Optional[Path] = None,
    genome: str = "hg19",
    threads: int = 0,
    verbose: bool = False,
    *,
    output_format: str,
    compress: bool,
) -> Path:
    """Compute UXM for one sample and return the table written. Raises on failure.

    The body of the `uxm` command, with plain arguments. run-all used to call the typer command
    as a function and pass no `threads`, so `OptionInfo > 0` raised `TypeError` on every run and
    run-all's catch-all hid it. The command converts a failure to `typer.Exit`; this raises the
    real exception, so a caller can record what went wrong.
    """
    from .assets import AssetManager

    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    # Configure Rust thread pool
    if threads > 0:
        try:
            _core.configure_threads(threads)
            logger.info(f"Configured {threads} threads for parallel processing")
        except Exception as e:
            logger.warning(f"Could not configure threads: {e}")

    # Input validation
    if not bam_input.exists():
        raise FileNotFoundError(f"Input BAM not found: {bam_input}")

    if not str(bam_input).endswith(".bam"):
        raise ValueError(f"Input must be a .bam file: {bam_input}")

    # Validate user-provided override files
    from .core.asset_validation import validate_file, FileSchema

    if mark_input and mark_input.exists():
        logger.debug(f"Validating user-provided marker file: {mark_input}")
        validate_file(mark_input, FileSchema.BED3)

    assets = AssetManager(genome)  # ValueError for an unknown genome
    logger.info(f"Genome: {assets.raw_genome} → {assets.genome_dir}")

    # Default marker file via AssetManager
    if mark_input is None:
        try:
            mark_input = assets.resolve("methylation_markers")
        except FileNotFoundError as e:
            raise FileNotFoundError(
                f"Methylation markers not found for {genome}; provide --mark-input ({e})"
            ) from e
        logger.info(f"Using default marker file: {mark_input}")

    if not mark_input.exists():
        raise FileNotFoundError(f"Marker file not found: {mark_input}")

    # Create output directory
    output.mkdir(parents=True, exist_ok=True)

    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bam_input.stem.replace(".bam", "")

    output_file = output / f"{sample_name}.UXM.tsv"

    logger.info(f"Processing {bam_input.name}")

    # Call Rust backend with all parameters
    _core.uxm.calculate_uxm(
        str(bam_input),
        str(mark_input),
        str(output_file),
        20,  # map_quality
        1,  # min_cpg
        METHY_THRESHOLD,  # methy_threshold (ratio >= this -> M)
        UNMETHY_THRESHOLD,  # unmethy_threshold (ratio <= this -> U)
        "SE",  # pe_type (single-end default)
    )

    # Rust always writes TSV; convert to Parquet/both afterward if requested.
    # This two-step approach avoids re-implementing Arrow schema in Rust for UXM.
    logger.info(f"UXM complete: {output_file} (raw TSV from Rust)")
    logger.debug(
        f"UXM post-process: output_format={output_format!r}, compress={compress}"
    )
    if output_format != "tsv":
        logger.debug(f"UXM: converting {output_file.name} → format={output_format!r}")
        # read_exact_table: Rust wrote this exact TSV. read_table is
        # parquet-first and would prefer a stale `.parquet` sibling from an
        # earlier run into the same directory (c92ed86).
        df = read_exact_table(output_file)
        if df is not None:
            write_table(
                df,
                output_file.with_suffix(""),
                output_format=output_format,
                compress=compress,
            )
        else:
            logger.warning(
                f"UXM: could not read {output_file.name} for format conversion"
            )
    elif compress:
        # Re-write as gzip TSV (Rust writes uncompressed TSV, Python gzips it)
        logger.debug(f"UXM: compressing {output_file.name} → .tsv.gz")
        # read_exact_table: Rust wrote this exact TSV. read_table is
        # parquet-first and would prefer a stale `.parquet` sibling from an
        # earlier run into the same directory (c92ed86).
        df = read_exact_table(output_file)
        if df is not None:
            write_table(
                df, output_file.with_suffix(""), output_format="tsv", compress=True
            )
    return output_file


def uxm(
    bam_input: Path = typer.Option(
        ..., "--input", "-i", help="Input bisulfite BAM file"
    ),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(
        None,
        "--sample-name",
        "-s",
        help="Sample name for output file (default: derived from input filename)",
    ),
    mark_input: Optional[Path] = typer.Option(
        None, "--mark-input", "-m", help="Path to genomic marker file"
    ),
    genome: str = typer.Option(
        "hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"
    ),
    threads: int = typer.Option(
        0, "--threads", "-t", help="Number of threads (0=all cores)"
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable verbose logging"
    ),
    output_format: str = typer.Option(
        "tsv",
        "--output-format",
        "-F",
        help="Output format: tsv | parquet | both (default: tsv)",
    ),
    compress: bool = typer.Option(
        False,
        "--compress",
        help="Gzip-compress TSV output (only when format is tsv or both)",
    ),
):
    """
    Calculate Fragment-level Methylation (UXM) features for a single sample.

    Input: Bisulfite BAM file
    Output: {sample}.UXM.tsv file with fragment-level methylation scores
    """
    try:
        written = uxm_stage(
            bam_input,
            output,
            sample_name,
            mark_input,
            genome,
            threads=threads,
            verbose=verbose,
            output_format=output_format,
            compress=compress,
        )
    except Exception as e:
        logger.error(f"UXM calculation failed: {e}")
        if verbose:
            import traceback

            traceback.print_exc()
        raise typer.Exit(1)
    # Into the sample's record, if it has one: rerunning this command (at run-all's --output-format) is how a stage that failed
    # under run-all gets fixed. Success only -- a failure here exits non-zero, and may not have
    # touched the table a good earlier run left.
    from .core.output_provenance import Outcome, record_stage_outcomes

    record_stage_outcomes(
        output,
        written.name[: -len(".UXM.tsv")],  # the sample name the stage used
        {"uxm": Outcome("ran")},
        output_format=output_format,
        compress=compress,
        owns_record=False,
    )
