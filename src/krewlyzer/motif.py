"""
Motif feature extraction from BAM files.

Extracts End Motif (EDM), Breakpoint Motif (BPM), and Motif Diversity Score (MDS)
from cfDNA fragment endpoints. These features capture tissue-of-origin signals.

Uses Rust backend for accelerated extraction.

Output Files:
    - {sample}.EndMotif.tsv: End motif 4-mer frequencies (256 columns)
    - {sample}.BreakPointMotif.tsv: Breakpoint motif frequencies
    - {sample}.MDS.tsv: Motif Diversity Score summary

Panel Mode: Generates separate on-target (.ontarget.tsv) files.

Note: For fragment extraction (BED.gz), use `krewlyzer extract` instead.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(
    level="INFO",
    handlers=[RichHandler(console=console, show_time=True, show_path=False)],
    format="%(message)s",
)
logger = logging.getLogger("motif")

# Rust backend is required
from krewlyzer import _core

# Import asset resolution and startup banner
from .core.asset_resolution import (
    resolve_exclude_regions,
    resolve_pon_model,
    resolve_target_regions,
)
from .core.logging import log_startup_banner, ResolvedAsset
from . import __version__


def motif(
    bam_input: Path = typer.Option(..., "--input", "-i", help="Input BAM file"),
    genome_reference: Path = typer.Option(
        ..., "-r", "--reference", help="Reference genome FASTA (indexed)"
    ),
    output: Path = typer.Option(..., "-o", "--output", help="Output directory"),
    genome: str = typer.Option(
        "hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"
    ),
    target_regions: Optional[Path] = typer.Option(
        None,
        "-T",
        "--target-regions",
        help="Target regions BED (for panel data: generates on/off-target motifs)",
    ),
    skip_target_regions: bool = typer.Option(
        False,
        "--skip-target-regions",
        help="Disable panel mode even when --assay has bundled targets",
    ),
    assay: Optional[str] = typer.Option(
        None,
        "--assay",
        "-A",
        help="Assay code (xs1/xs2) for auto-loading bundled assets",
    ),
    pon_model: Optional[Path] = typer.Option(
        None, "-P", "--pon-model", help="PON model for MDS z-score computation"
    ),
    pon_variant: str = typer.Option(
        "all_unique",
        "--pon-variant",
        help="PON variant: 'all_unique' (default, max coverage) or 'duplex' (highest accuracy)",
    ),
    skip_pon: bool = typer.Option(
        False, "--skip-pon", help="Skip PON z-score normalization"
    ),
    kmer: int = typer.Option(4, "-k", "--kmer", help="K-mer size for motif extraction"),
    # Read filters and blacklist: run-all's, so motif counts the fragments run-all's tables do
    exclude_regions: Optional[Path] = typer.Option(
        None,
        "-x",
        "--exclude-regions",
        help="Exclude regions BED file (default: the genome's bundled blacklist)",
    ),
    mapq: int = typer.Option(20, "--mapq", "-q", help="Minimum mapping quality"),
    minlen: int = typer.Option(65, "--minlen", help="Minimum fragment length"),
    maxlen: int = typer.Option(1000, "--maxlen", help="Maximum fragment length"),
    skip_duplicates: bool = typer.Option(
        True, "--skip-duplicates/--no-skip-duplicates", help="Skip duplicate reads"
    ),
    # (`--chromosomes` was removed: it was parsed and ignored, so a "chr1 only" run counted
    # motifs over the whole genome.)
    sample_name: Optional[str] = typer.Option(
        None,
        "--sample-name",
        "-s",
        help="Sample name for output files (default: derived from BAM filename)",
    ),
    require_proper_pair: bool = typer.Option(
        True,
        "--require-proper-pair/--no-require-proper-pair",
        help="Require proper pairs (disable for duplex/consensus BAMs)",
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
    threads: int = typer.Option(
        0, "--threads", "-t", help="Number of threads (0=all cores)"
    ),
):
    """
    Extract k-mer motif features from a BAM file.

    Output:
    - {sample}.EndMotif: End motif k-mer frequencies
    - {sample}.BreakPointMotif: Breakpoint motif k-mer frequencies
    - {sample}.MDS: Motif Diversity Score

    Note: For fragment extraction (BED.gz), use `krewlyzer extract` instead.
    """
    from .assets import AssetManager

    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("core.motif_processor").setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    # Configure Rust thread pool
    if threads > 0:
        try:
            _core.configure_threads(threads)
            logger.info(f"Configured {threads} threads")
        except Exception as e:
            logger.warning(f"Could not configure threads: {e}")

    # Input validation
    if not bam_input.exists():
        logger.error(f"BAM file not found: {bam_input}")
        raise typer.Exit(1)

    if not str(bam_input).endswith(".bam"):
        logger.error(f"Input must be a .bam file: {bam_input}")
        raise typer.Exit(1)

    if not genome_reference.exists():
        logger.error(f"Reference genome not found: {genome_reference}")
        raise typer.Exit(1)

    # Create output directory
    output.mkdir(parents=True, exist_ok=True)

    # Derive sample name (use provided or derive from BAM filename)
    if sample_name is None:
        sample_name = bam_input.stem.replace(".bam", "")

    # Initialize AssetManager
    assets = AssetManager(genome)

    # ═══════════════════════════════════════════════════════════════════
    # ASSET RESOLUTION
    # ═══════════════════════════════════════════════════════════════════
    try:
        resolved_pon_path, pon_source = resolve_pon_model(
            explicit_path=pon_model,
            assay=assay,
            skip_pon=skip_pon,
            assets=assets,
            variant=pon_variant,
            log=logger,
        )
        # Before any work: an unusable PON, or one other than the PON this sample's record
        # names, stops the run before a table is written.
        from .core.pon_integration import load_requested_pon

        requested_pon = load_requested_pon(
            resolved_pon_path,
            record_dir=output,
            sample=sample_name,
            output_format=output_format,
            compress=compress,
        )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)

    try:
        resolved_target_path, target_source = resolve_target_regions(
            explicit_path=target_regions,
            assay=assay,
            skip_target_regions=skip_target_regions,
            assets=assets,
            log=logger,
        )
        # The blacklist before the check, which compares it with the one the BED was made with.
        resolved_exclude_path, exclude_source = resolve_exclude_regions(
            exclude_regions, assets, logger
        )
        # Before any work: other settings than the folder's record would leave the first run's
        # extra tables looking current. A `ConfigurationMismatch` is a ValueError.
        from .core.output_provenance import (
            filter_settings,
            refuse_a_different_configuration,
        )

        # The filters and the blacklist too: motif reads the BAM, so with other ones its tables
        # would come from fragments the sample's BED does not hold.
        refuse_a_different_configuration(
            output,
            sample_name,
            {
                "genome": genome,
                "assay": assay,
                "target_regions": resolved_target_path,
                **filter_settings(
                    mapq=mapq,
                    minlen=minlen,
                    maxlen=maxlen,
                    skip_duplicates=skip_duplicates,
                    require_proper_pair=require_proper_pair,
                    blacklist=resolved_exclude_path,
                ),
            },
            output_format=output_format,
            compress=compress,
        )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)

    is_panel_mode = resolved_target_path is not None and resolved_target_path.exists()

    # ═══════════════════════════════════════════════════════════════════
    # STARTUP BANNER
    # ═══════════════════════════════════════════════════════════════════
    log_startup_banner(
        tool_name="motif",
        version=__version__,
        inputs={
            "BAM": str(bam_input.name),
            "Reference": str(genome_reference.name),
            "Output": str(output),
        },
        config={
            "Genome": f"{assets.raw_genome} → {assets.genome_dir}",
            "Assay": assay or "None",
            "Mode": "Panel" if is_panel_mode else "WGS",
            "K-mer": str(kmer),
        },
        assets=[
            ResolvedAsset("PON", resolved_pon_path, pon_source),
            ResolvedAsset("Targets", resolved_target_path, target_source),
            ResolvedAsset("Exclude regions", resolved_exclude_path, exclude_source),
        ],
        logger=logger,
    )

    # Pre-check BAM compatibility with current filters
    if require_proper_pair:
        from .core.bam_utils import check_bam_compatibility

        logger.info("Checking BAM read compatibility with filters...")
        compat = check_bam_compatibility(
            bam_input, require_proper_pair, skip_duplicates, mapq
        )

        if compat["pass_rate"] < 0.01 and compat["total_sampled"] > 100:
            # Critical: Almost no reads would pass!
            console.print("\n[bold red]⚠️  FILTER COMPATIBILITY WARNING[/bold red]\n")
            console.print(
                f"Only [bold]{compat['pass_rate']:.2%}[/bold] of sampled reads would pass current filters.\n"
            )

            for issue in compat["issues"]:
                console.print(f"  • {issue}")

            if compat["suggested_flags"]:
                suggested = " ".join(compat["suggested_flags"])
                console.print("\n[bold yellow]Suggested command:[/bold yellow]")
                console.print(
                    f"  krewlyzer motif {bam_input.name} -r {genome_reference.name} -o {output} [bold]{suggested}[/bold]\n"
                )

                logger.error(f"Filter mismatch detected. Re-run with: {suggested}")
                raise typer.Exit(1)

    try:
        logger.info(f"Extracting motif features from {bam_input.name}")

        # Use unified extract_sample() for extraction
        from .core.sample_processor import extract_sample, write_motif_tables

        result = extract_sample(
            input_path=bam_input,
            reference=genome_reference,
            mapq=mapq,
            minlen=minlen,
            maxlen=maxlen,
            kmer=kmer,
            threads=threads,
            skip_duplicates=skip_duplicates,
            require_proper_pair=require_proper_pair,
            target_regions=resolved_target_path,
            exclude_regions=resolved_exclude_path,
            write_bed=False,  # Motif command doesn't write BED
            output_dir=None,
            sample_name=sample_name,
            genome=genome,
        )

        logger.info(f"Processed {result.fragment_count:,} fragments")

        # Load PON model if resolved and not skipped
        pon = None
        if resolved_pon_path and resolved_pon_path.exists():
            # The guarded loader, not `PonModel.load`.
            #
            # `load_pon_model` is where the version floor lives. Loading
            # directly here meant `krewlyzer motif` would happily score
            # against a pre-0.9.0 model that run-all refuses -- invariant #6,
            # single-tool and run-all behaving differently.
            pon = requested_pon  # loaded, and checked, before the extraction

        # Write motif outputs using shared function
        _, motif_outcomes = write_motif_tables(
            result=result,
            output_dir=output,
            pon=pon,
            include_ontarget=is_panel_mode,
            output_format=output_format,
            compress=compress,
        )
        from .core.output_provenance import record_pon_used, record_stage_outcomes

        # What the MDS z-scores did, as run-all records it; a single tool, so only what it tried.
        record_stage_outcomes(
            output,
            sample_name,
            motif_outcomes,
            output_format=output_format,
            compress=compress,
            owns_record=False,
        )

        # Recorded as scored only if the MDS z-score reached the table: from the outcome, not
        # from the PON having a baseline, which also held when the z-score then failed.
        scored_with = pon if motif_outcomes["mds_z"].status == "ran" else None
        record_pon_used(
            output,
            sample_name,
            scored_with,
            resolved_pon_path,
            output_format=output_format,
            compress=compress,
        )

        # Additional aberrant k-mer analysis if PON provided
        if pon and hasattr(pon, "mds_baseline") and pon.mds_baseline:
            try:
                aberrant = pon.mds_baseline.get_aberrant_kmers(
                    result.kmer_frequencies, threshold=3.0
                )
                if aberrant:
                    logger.info(f"Found {len(aberrant)} aberrant k-mers (|z| > 3.0)")
            except Exception as e:
                logger.debug(f"Could not compute aberrant k-mers: {e}")

        logger.info(f"Motif extraction complete (MDS={result.mds_score:.4f})")

    except Exception as e:
        logger.error(f"Motif extraction failed: {e}")
        if verbose:
            import traceback

            traceback.print_exc()
        raise typer.Exit(1)
