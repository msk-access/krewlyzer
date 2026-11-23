"""
Fragment extraction from BAM files.

Extracts cfDNA fragments with configurable read filters (MAPQ, length, proper pair).
First step in the pipeline - generates .bed.gz for downstream feature extraction.

Uses Rust backend for accelerated extraction.

Output Files:
    - {sample}.bed.gz: Fragment coordinates (chrom, start, end, strand, gc%, gc_weight)
    - {sample}.correction_factors.{parquet,tsv.gz,tsv}: LOESS-based GC correction factors
    - {sample}.metadata.{tsv,tsv.gz,parquet}: The sample's record (filters, GC status)

GC Correction: Computes per-fragment GC weights for all downstream tools.
Panel Mode: GC model trained on off-target reads only for unbiased correction.

BAM Compatibility: Auto-detects duplex BAMs and warns to use --no-require-proper-pair.
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
logger = logging.getLogger("extract")

# Rust backend is required
from krewlyzer import _core

# Shared BAM utilities
from .core.bam_utils import check_bam_compatibility

# Asset resolution for --assay auto-loading
from .core.asset_resolution import resolve_exclude_regions, resolve_target_regions
from .core.gc_assets import factor_table_base, factor_table_path


def extract(
    bam_input: Path = typer.Option(
        ..., "--input", "-i", help="Input BAM file (sorted, indexed)"
    ),
    genome_reference: Path = typer.Option(
        ..., "-r", "--reference", help="Reference genome FASTA (indexed)"
    ),
    output: Path = typer.Option(..., "-o", "--output", help="Output directory"),
    # GC correction options
    genome: str = typer.Option(
        "hg19",
        "-G",
        "--genome",
        help="Genome build for bundled assets: targets, blacklist, GC reference (hg19/GRCh37/hg38/GRCh38)",
    ),
    gc_correct: bool = typer.Option(
        True, "--gc-correct/--no-gc-correct", help="Compute GC correction factors"
    ),
    # Configurable filters
    exclude_regions: Optional[Path] = typer.Option(
        None, "-x", "--exclude-regions", help="Exclude regions BED file"
    ),
    target_regions: Optional[Path] = typer.Option(
        None,
        "-T",
        "--target-regions",
        help="Target regions BED (for panel data: GC model computed on off-target reads only)",
    ),
    assay: Optional[str] = typer.Option(
        None,
        "--assay",
        "-A",
        help="Assay code (xs1/xs2) for auto-loading bundled targets",
    ),
    skip_target_regions: bool = typer.Option(
        False,
        "--skip-target-regions",
        help="Disable panel mode even when --assay has bundled targets",
    ),
    mapq: int = typer.Option(20, "--mapq", "-q", help="Minimum mapping quality"),
    minlen: int = typer.Option(65, "--minlen", help="Minimum fragment length"),
    maxlen: int = typer.Option(
        1000,
        "--maxlen",
        help="Maximum fragment length (default: 1000 for extended FSD range)",
    ),
    skip_duplicates: bool = typer.Option(
        True, "--skip-duplicates/--no-skip-duplicates", help="Skip duplicate reads"
    ),
    require_proper_pair: bool = typer.Option(
        True,
        "--require-proper-pair/--no-require-proper-pair",
        help="Require proper pairs",
    ),
    # Other options. (`--chromosomes` was removed: it was parsed and ignored, so a "chr1 only"
    # run wrote a whole-genome BED under a complete-looking record.)
    sample_name: Optional[str] = typer.Option(
        None,
        "--sample-name",
        "-s",
        help="Sample name for output files (default: derived from BAM filename)",
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
    Extract cfDNA fragments from BAM to BED.gz with full filter control.

    Filters applied (always on):
    - Skip unmapped reads
    - Skip secondary alignments
    - Skip supplementary alignments
    - Skip failed QC reads

    Filters applied (configurable):
    - Mapping quality threshold
    - Fragment length range
    - Skip duplicates
    - Require proper pairs
    - Exclude genomic regions

    Output:
    - {sample}.bed.gz: Fragment coordinates with GC content
    - {sample}.bed.gz.tbi: Tabix index
    - {sample}.metadata.{tsv,tsv.gz,parquet}: The sample's record, which extract refuses
      to replace once run-all or a single tool has written into it
    - {sample}.correction_factors.{parquet,tsv.gz,tsv}: GC correction factors (if --gc-correct)
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

    # The genome's assets: the assay's targets and the blacklist come from them whatever
    # --gc-correct says, as in run-all; the flag decides only whether factors are computed.
    # A genome AssetManager does not know is an error, as in run-all: carrying on would mean
    # no targets and no blacklist, not just no GC correction. A missing or unfetched (LFS
    # pointer) target set or blacklist is an error in the two resolvers below.
    try:
        assets = AssetManager(genome)
        logger.info(f"Genome: {assets.raw_genome} -> {assets.genome_dir}")
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)

    # Resolve target regions (auto-load from --assay if not explicit)
    try:
        resolved_target_path, target_source = resolve_target_regions(
            explicit_path=target_regions,
            assay=assay,
            skip_target_regions=skip_target_regions,
            assets=assets,
            log=logger,
        )
        if resolved_target_path and target_source == "bundled":
            logger.info(
                f"Auto-loaded target regions for assay '{assay}': {resolved_target_path.name}"
            )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)

    # The blacklist, by the rule run-all and motif use.
    try:
        exclude_regions, _ = resolve_exclude_regions(exclude_regions, assets, logger)
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)

    # Create output directory
    output.mkdir(parents=True, exist_ok=True)

    # Derive sample name (use provided or derive from BAM filename)
    if sample_name is None:
        sample_name = bam_input.stem.replace(".bam", "")

    # Before anything is written: `extract` writes the sample's record from scratch, so with other
    # settings than the folder's it would replace them, and a later rerun would pass the check.
    from .core.output_provenance import (
        ConfigurationMismatch,
        RecordInUse,
        refuse_a_different_configuration,
        refuse_to_replace_a_later_record,
    )

    try:
        refuse_a_different_configuration(
            output,
            sample_name,
            {
                "genome": genome,
                "assay": assay,
                "target_regions": resolved_target_path,
            },
            output_format=output_format,
            compress=compress,
        )
        # Then even with the same settings: a record run-all or a single tool wrote into is
        # theirs. Second, so a run with other settings is told which setting differs.
        refuse_to_replace_a_later_record(output, sample_name)
    except (ConfigurationMismatch, RecordInUse) as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)

    try:
        # Pre-check BAM compatibility with current filters
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
                    f"  krewlyzer extract {bam_input.name} -r {genome_reference.name} -o {output} [bold]{suggested}[/bold]\n"
                )

                logger.error(f"Filter mismatch detected. Re-run with: {suggested}")
                raise typer.Exit(1)
        elif compat["pass_rate"] < 0.5 and compat["issues"]:
            # Warning but continue
            console.print("\n[yellow]⚠️  Filter compatibility note:[/yellow]")
            for issue in compat["issues"]:
                console.print(f"  • {issue}")
            console.print()

        logger.info(f"Extracting fragments from {bam_input.name}")
        logger.info(
            f"Filters: mapq>={mapq}, length=[{minlen},{maxlen}], skip_dup={skip_duplicates}, proper_pair={require_proper_pair}"
        )

        if resolved_target_path:
            logger.info(
                f"Panel mode: GC model will use off-target reads only (targets: {resolved_target_path.name})"
            )

        # Use unified extract_sample() for extraction
        from .core.sample_processor import extract_sample, write_extraction_outputs

        result = extract_sample(
            input_path=bam_input,
            reference=genome_reference,
            mapq=mapq,
            minlen=minlen,
            maxlen=maxlen,
            kmer=4,
            threads=threads,
            skip_duplicates=skip_duplicates,
            require_proper_pair=require_proper_pair,
            target_regions=resolved_target_path,
            exclude_regions=exclude_regions,
            write_bed=True,
            output_dir=output,
            sample_name=sample_name,
            genome=genome,
        )

        logger.info(f"Extracted {result.fragment_count:,} fragments")

        # Write extraction outputs (BED index, metadata, GC factors)
        written = write_extraction_outputs(
            result=result,
            output_dir=output,
            genome=genome,
            assay=assay,
            compute_gc_factors=gc_correct,
            output_format=output_format,
            compress=compress,
        )

        # On-target GC correction (panel mode only) - additional to standard outputs
        if (
            resolved_target_path
            and gc_correct
            and len(result.gc_observations_ontarget) > 0
        ):
            try:
                gc_ref = assets.resolve("gc_reference")
                valid_regions = assets.resolve("valid_regions")
                factors_ontarget = factor_table_path(output, sample_name, "on_target")
                logger.info(
                    f"Computing ON-TARGET GC correction factors from {len(result.gc_observations_ontarget)} observation bins..."
                )
                n_factors_on = _core.gc.compute_and_write_gc_factors(
                    result.gc_observations_ontarget,
                    str(gc_ref),
                    str(valid_regions),
                    str(factors_ontarget),
                    output_format=output_format,
                    compress=compress,
                )
                if n_factors_on > 0:
                    logger.info(
                        f"Computed {n_factors_on} on-target correction factors: {factors_ontarget}"
                    )
                else:
                    logger.warning(
                        "No on-target GC factor could be fitted; no on-target table written"
                    )
            except Exception as e:
                logger.warning(f"On-target GC correction failed: {e}")

        # Name only what was written: `written` holds a factor table only when one exists.
        logger.info(f"Output: {result.bed_path}")
        logger.info(f"Metadata: {written['metadata']}")
        if "gc_factors" in written:
            logger.info(
                f"Correction factors: {factor_table_base(output, sample_name, 'off_target')}.*"
            )

    except Exception as e:
        logger.error(f"Fragment extraction failed: {e}")
        if verbose:
            import traceback

            traceback.print_exc()
        raise typer.Exit(1)
