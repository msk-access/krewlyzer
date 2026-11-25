"""
Mutant Fragment Size Distribution (mFSD) calculation.

Calculates mFSD features for a single sample.
Uses Rust backend for accelerated computation with optional GC correction.

Supports all small variant types: SNV, MNV, Insertion, Deletion, Complex.
4-way fragment classification: REF, ALT, NonREF, N.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

from .core.output_utils import (
    convert_rust_tsv,
    table_paths,
)  # Unified I/O for Parquet support

console = Console(stderr=True)
# Guard: only configure root logging if not already configured (e.g. by pytest or parent app).
# This prevents double-logging when mfsd is imported as a library rather than called as a CLI.
if not logging.root.handlers:
    logging.basicConfig(
        level="INFO",
        handlers=[RichHandler(console=console, show_time=True, show_path=False)],
        format="%(message)s",
    )
logger = logging.getLogger("mfsd")

# Rust backend is required
from krewlyzer import _core


def _sample_name(bam_input: Path, sample_name: Optional[str]) -> str:
    """The name given, else the BAM's: one rule for the stage and the folder check."""
    return sample_name or bam_input.stem.replace(".bam", "")


def mfsd_stage(
    bam_input: Path,
    input_file: Path,
    output: Path,
    sample_name: Optional[str] = None,
    reference: Optional[Path] = None,
    correction_factors: Optional[Path] = None,
    mapq: int = 20,
    minlen: int = 65,
    maxlen: int = 1000,
    skip_duplicates: bool = True,
    require_proper_pair: bool = False,
    duplex: bool = False,
    output_distributions: bool = False,
    min_baseq: int = 20,
    verbose: bool = False,
    threads: int = 0,
    *,
    output_format: str,
    compress: bool,
) -> Path:
    """Compute mFSD for one sample and return the table written. Raises on failure.

    The body of the `mfsd` command, with plain arguments. The command converts every failure to
    `typer.Exit`, so a caller recording what went wrong would only ever see `Exit`; run-all calls
    this instead, and so can record the real class.
    """
    # Set log level based on verbose flag
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
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
        raise FileNotFoundError(f"BAM file not found: {bam_input}")

    if not str(bam_input).endswith(".bam"):
        raise ValueError(f"Input must be a .bam file: {bam_input}")

    if not input_file.exists():
        raise FileNotFoundError(f"Variant file not found: {input_file}")

    # Validate optional paths
    if reference and not reference.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference}")

    if correction_factors and not correction_factors.exists():
        raise FileNotFoundError(
            f"Correction factors file not found: {correction_factors}"
        )

    # Create output directory
    output.mkdir(parents=True, exist_ok=True)

    sample_name = _sample_name(bam_input, sample_name)

    output_file = output / f"{sample_name}.mFSD.tsv"

    logger.info(f"Processing {bam_input.name}")

    # Detect input type
    input_type = (
        "vcf"
        if str(input_file).endswith(".vcf") or str(input_file).endswith(".vcf.gz")
        else "maf"
    )
    logger.info(f"Detected input type: {input_type}")

    # Log GC correction status
    if correction_factors:
        logger.info(f"GC correction enabled using: {correction_factors}")
    else:
        logger.info("GC correction disabled (no --correction-factors provided)")

    # Log duplex mode
    if duplex:
        logger.info(
            "Duplex mode enabled: using fgbio cD tag or Marianas read names for family size weighting"
        )

    # Call Rust backend with optional reference and correction factors
    _core.mfsd.calculate_mfsd(
        str(bam_input),
        str(input_file),
        str(output_file),
        input_type,
        mapq,
        minlen,
        maxlen,
        output_distributions,
        str(reference) if reference else None,
        str(correction_factors) if correction_factors else None,
        require_proper_pair,
        duplex,
        min_baseq=min_baseq,
    )

    # Rust always writes TSV; it becomes exactly the run's formats here (Parquet, gzip), and the
    # raw TSV goes unless it is one of them. This two-step approach avoids re-implementing the
    # Arrow schema in Rust for mFSD.
    logger.info(f"mFSD complete: {output_file} (raw TSV from Rust)")
    logger.debug(
        f"mFSD post-process: output_format={output_format!r}, compress={compress}"
    )
    convert_rust_tsv(output_file, output_format=output_format, compress=compress)

    if output_distributions:
        dist_file = output_file.with_suffix(".distributions.tsv")
        logger.info(f"Distributions: {dist_file}")
        convert_rust_tsv(dist_file, output_format=output_format, compress=compress)
    return table_paths(
        output_file.with_suffix(""), output_format=output_format, compress=compress
    )[0]


def mfsd(
    bam_input: Path = typer.Option(..., "--input", "-i", help="Input BAM file"),
    input_file: Path = typer.Option(
        ..., "--input-file", "-V", help="VCF or MAF file with variants"
    ),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(
        None,
        "--sample-name",
        "-s",
        help="Sample name for output file (default: derived from input filename)",
    ),
    reference: Optional[Path] = typer.Option(
        None, "--reference", "-g", help="Reference genome FASTA (for GC correction)"
    ),
    correction_factors: Optional[Path] = typer.Option(
        None,
        "--correction-factors",
        "-F",
        help="Pre-computed correction_factors.tsv (from extract/run-all)",
    ),
    mapq: int = typer.Option(20, "--mapq", "-q", help="Minimum mapping quality"),
    minlen: int = typer.Option(
        65, "--minlen", help="Minimum fragment length (filters discordant reads)"
    ),
    maxlen: int = typer.Option(
        1000,
        "--maxlen",
        help="Maximum fragment length (default: 1000 for extended FSD range)",
    ),
    skip_duplicates: bool = typer.Option(
        True,
        "--skip-duplicates/--no-skip-duplicates",
        help="Skip duplicate reads (always enabled in Rust backend)",
    ),
    require_proper_pair: bool = typer.Option(
        False,
        "--require-proper-pair/--no-require-proper-pair",
        help="Require proper pairs (disable for duplex BAMs)",
    ),
    duplex: bool = typer.Option(
        False,
        "--duplex",
        "-D",
        help="Enable duplex weighting (fgbio cD tag or Marianas read names). Only use for duplex consensus BAMs.",
    ),
    output_distributions: bool = typer.Option(
        False,
        "--output-distributions",
        "-d",
        help="Output per-variant size distributions",
    ),
    min_baseq: int = typer.Option(
        20,
        "--min-baseq",
        "-Q",
        help="Minimum base quality at variant position (filters low-quality evidence)",
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable verbose/debug logging"
    ),
    threads: int = typer.Option(
        0, "--threads", "-t", help="Number of threads (0=all cores)"
    ),
    output_format: str = typer.Option(
        "tsv",
        "--output-format",
        help="Output format: tsv | parquet | both (default: tsv)",
    ),
    compress: bool = typer.Option(
        False,
        "--compress",
        help="Gzip-compress TSV output (only when format is tsv or both)",
    ),
):
    """
    Calculate Mutant Fragment Size Distribution (mFSD) features for a single sample.

    Supports all small variant types: SNV, MNV, Insertion, Deletion, Complex.

    Fragment classification:
    - REF: Supports reference allele (healthy cfDNA baseline)
    - ALT: Supports alternate allele (tumor signal)
    - NonREF: Non-REF, non-ALT (sequencing errors, subclones)
    - N: Contains N at variant position (low quality)

    GC Correction:
    - If --correction-factors is provided, fragment counts are weighted by GC correction factors
    - These factors can be computed automatically by run-all or standalone tools

    Input: BAM file and VCF/MAF file with variants
    Output:
    - {sample}.mFSD.tsv: Summary statistics (46 columns)
    - {sample}.mFSD.distributions.tsv: Per-size fragment counts (optional, with -d)
    """
    # Before any work: mFSD's tables depend on no assay or targets, but a run in another format
    # than the folder's leaves its earlier tables beside the new ones, and downstream reads the
    # Parquet. The filters too: mFSD reads the BAM, and run-all gives it the sample's mapq and
    # length range (`--mfsd-bam` included), so other ones would put a second fragment set under
    # the record. Not the other two: run-all always runs mFSD without the proper-pair filter, and
    # the Rust reader drops duplicates whatever `--skip-duplicates` says. run-all checks its own
    # folder; it calls `mfsd_stage`, not this command.
    from .core.output_provenance import (
        ConfigurationMismatch,
        filter_settings,
        refuse_a_different_configuration,
        refuse_tables_of_unrecorded_filters,
    )
    from .validate.contract import STAGE_TABLES

    try:
        refuse_a_different_configuration(
            output,
            _sample_name(bam_input, sample_name),
            filter_settings(mapq=mapq, minlen=minlen, maxlen=maxlen),
            output_format=output_format,
            compress=compress,
        )
        # Nor beside tables another step made with filters the record does not name.
        refuse_tables_of_unrecorded_filters(
            output,
            _sample_name(bam_input, sample_name),
            own=STAGE_TABLES["mfsd"],
            creates_the_record=False,
            output_format=output_format,
            compress=compress,
        )
    except ConfigurationMismatch as e:  # `UnrecordedFilters` is one
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)

    try:
        mfsd_stage(
            bam_input,
            input_file,
            output,
            sample_name,
            reference=reference,
            correction_factors=correction_factors,
            mapq=mapq,
            minlen=minlen,
            maxlen=maxlen,
            skip_duplicates=skip_duplicates,
            require_proper_pair=require_proper_pair,
            duplex=duplex,
            output_distributions=output_distributions,
            min_baseq=min_baseq,
            verbose=verbose,
            threads=threads,
            output_format=output_format,
            compress=compress,
        )
    except Exception as e:
        logger.error(f"mFSD calculation failed: {e}")
        raise typer.Exit(1)
    # Into the sample's record, if it has one: rerunning this command (at run-all's --output-format) is how a stage that failed
    # under run-all gets fixed. Success only -- a failure here exits non-zero, and may not have
    # touched the table a good earlier run left.
    from .core.output_provenance import Outcome, record_stage_outcomes

    record_stage_outcomes(
        output,
        _sample_name(bam_input, sample_name),
        {"mfsd": Outcome("ran")},
        output_format=output_format,
        compress=compress,
        owns_record=False,
    )
