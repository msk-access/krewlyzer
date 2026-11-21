"""
Per-Region MDS (Motif Diversity Score) CLI Tool

Calculates MDS at the exon/target level from BAM files.
Outputs:
- {sample}.MDS.exon.tsv: Per-exon/target MDS scores
- {sample}.MDS.gene.tsv: Gene-level aggregated MDS scores

Based on Helzer et al. (2025) methodology.
"""

from pathlib import Path
from typing import Optional
import typer
import logging

logger = logging.getLogger("krewlyzer.region_mds")

# Import asset resolution and startup banner
from .core.asset_resolution import resolve_pon_model
from .core.logging import log_startup_banner, ResolvedAsset
from .core.output_provenance import Outcome, record_stage_outcomes
from .core.output_utils import (
    read_exact_table,
    write_table,
    cleanup_intermediate_tsv,
)
from . import __version__


def _apply_exon_zscores(
    pon,
    output_exon: Path,
    output_exon_base: Path,
    output_format: str,
    compress: bool,
) -> "tuple[bool, int]":
    """Append `mds_z` to the exon table, keyed on ``(gene, name)``.

    Returns (written, n_scored): whether it wrote the final output, so the
    format-conversion block below knows not to run again on a file this already
    consumed, and how many exons got a z-score, so the PON is recorded as used
    only when it scored something.

    Split out rather than inlined: the gene block above is already deeply
    nested, and this is the same shape one level down.

    The two absences are counted separately. ``compute_zscore`` returns None
    for an exon the baseline has never seen -- a panel change, or a PON built
    for the other assay -- and NaN for one it has seen but could not measure a
    spread for. Both write NaN, correctly, but only the log distinguishes
    "rebuild the PON" from "this exon has no variance in the cohort".
    """
    import pandas as pd

    frame = read_exact_table(output_exon)
    if frame is None:
        logger.warning(f"Could not read exon output for z-score: {output_exon}")
        return False, 0
    if "gene" not in frame.columns or "name" not in frame.columns:
        logger.warning(
            f"Exon output lacks gene/name columns, cannot z-score: {output_exon}"
        )
        return False, 0

    z_scores = []
    n_absent = 0
    n_unscoreable = 0
    for _, row in frame.iterrows():
        z = pon.region_mds_exon.compute_zscore(
            str(row["gene"]), str(row["name"]), row.get("mds")
        )
        if z is None:
            n_absent += 1
            z_scores.append(float("nan"))
            continue
        if pd.isna(z):
            n_unscoreable += 1
        z_scores.append(z)

    frame["mds_z"] = z_scores
    write_table(frame, output_exon_base, output_format=output_format, compress=compress)
    cleanup_intermediate_tsv(output_exon, output_format, compress)

    n_scored = len(frame) - n_absent - n_unscoreable
    logger.info(f"Added exon z-scores for {n_scored}/{len(frame)} exons")
    if n_absent:
        logger.warning(
            f"MDS.exon: {n_absent}/{len(frame)} exons are absent from the PON's "
            "region_mds_exon baseline. Check the PON matches this panel."
        )
    if n_unscoreable:
        logger.warning(
            f"MDS.exon: {n_unscoreable}/{len(frame)} exons matched the baseline "
            "but it measured no usable spread, so the z-score is NaN rather "
            "than a value divided by a placeholder."
        )
    return True, n_scored


def _derived_sample_name(bam_input: Path) -> str:
    """The sample name when none is given: the BAM's stem, without a `.sorted` suffix."""
    name = bam_input.stem
    return name[: -len(".sorted")] if name.endswith(".sorted") else name


def region_mds_stage(
    bam_input: Path,
    reference: Path,
    output: Path,
    gene_bed: Optional[Path] = None,
    genome: str = "hg19",
    assay: Optional[str] = None,
    sample_name: Optional[str] = None,
    e1_only: bool = False,
    mapq: int = 20,
    minlen: int = 65,
    maxlen: int = 1000,
    pon_model: Optional[Path] = None,
    pon_variant: str = "all_unique",
    skip_pon: bool = False,
    threads: int = 0,
    verbose: bool = False,
    silent: bool = False,
    output_format: str = "tsv",
    compress: bool = False,
    owns_sample_record: bool = False,
):
    """Per-region MDS for one sample. Raises on failure.

    Returns ``(n_regions, n_genes, z_outcome)``: the last is what the PON z-scoring of
    `MDS.gene` did (`output_provenance.Outcome`), which carries on past its own failure.

    The body of the `region-mds` command, with plain arguments. run-all called the typer command
    as a function and passed no `pon_variant`, so a typer `OptionInfo` reached the PON resolver in
    its place. Named `region_mds_stage` because `run_region_mds` is the Rust entry point.
    ``owns_sample_record`` is set by run-all, which owns the sample's record: it skips the
    second-PON check.
    """
    from krewlyzer import _core
    from krewlyzer.assets import AssetManager

    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate inputs
    if not bam_input.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_input}")
    if not reference.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference}")

    # Validate user-provided override files
    from .core.asset_validation import validate_file, FileSchema

    if gene_bed and gene_bed.exists():
        logger.debug(f"Validating user-provided gene BED: {gene_bed}")
        validate_file(gene_bed, FileSchema.GENE_BED)

    # Resolve gene BED
    assets = AssetManager(genome)
    resolved_gene_bed = None
    gene_bed_source = "none"

    if gene_bed:
        # User-provided gene BED
        if not gene_bed.exists():
            raise FileNotFoundError(f"Gene BED not found: {gene_bed}")
        resolved_gene_bed = gene_bed
        gene_bed_source = "explicit"
        logger.debug(f"Using user-provided gene BED: {gene_bed}")
    elif assay:
        # Bundled gene BED from assay
        resolved_gene_bed = assets.get_gene_bed_for_mode(assay)
        if resolved_gene_bed is None:
            raise FileNotFoundError(
                f"Gene BED not available for assay '{assay}' and genome '{genome}'"
            )
        gene_bed_source = "bundled"
        logger.debug(f"Using bundled gene BED for assay '{assay}': {resolved_gene_bed}")
    else:
        # Default to WGS if no assay specified
        if assets.wgs_gene_bed_available:
            resolved_gene_bed = assets.wgs_gene_bed
            gene_bed_source = "bundled"
            logger.debug(f"Using WGS gene BED: {resolved_gene_bed}")
        else:
            raise FileNotFoundError(
                "No gene BED specified. Use --gene-bed, --assay, or ensure wgs.genes.bed.gz is available."
            )

    # ═══════════════════════════════════════════════════════════════════
    # ASSET RESOLUTION (PON)
    # ═══════════════════════════════════════════════════════════════════
    resolved_pon_path, pon_source = (
        resolve_pon_model(  # ValueError: --pon-model + --skip-pon
            explicit_path=pon_model,
            assay=assay,
            skip_pon=skip_pon,
            assets=assets,
            variant=pon_variant,
            log=logger,
        )
    )

    # Create output directory
    output.mkdir(parents=True, exist_ok=True)

    # Derive sample name: use CLI-provided value, fallback to BAM stem.
    # When called via run-all, wrapper.py passes the consistent sample ID;
    # standalone CLI can rely on the BAM filename fallback.
    if not sample_name:
        sample_name = _derived_sample_name(bam_input)
        logger.debug(f"Derived sample name from BAM: {sample_name}")
    else:
        logger.debug(f"Using provided sample name: {sample_name}")

    # Before any work: an unusable PON, or one other than the PON this sample's record names,
    # stops the run before a table is written.
    from krewlyzer.core.pon_integration import load_requested_pon

    requested_pon = load_requested_pon(  # ValueError (UnusablePonError) if unusable
        resolved_pon_path,
        record_dir=None if owns_sample_record else output,
        sample=sample_name,
        output_format=output_format,
        compress=compress,
    )

    output_exon_base = output / f"{sample_name}.MDS.exon"
    output_gene_base = output / f"{sample_name}.MDS.gene"
    # Rust expects .tsv paths for its initial write; Parquet added by Python after.
    # IMPORTANT: Do NOT use .with_suffix() — compound names like "MDS.exon" would
    # have ".exon" replaced by ".tsv", giving "MDS.tsv" instead of "MDS.exon.tsv".
    output_exon = output_exon_base.parent / (output_exon_base.name + ".tsv")
    output_gene = output_gene_base.parent / (output_gene_base.name + ".tsv")

    # ═══════════════════════════════════════════════════════════════════
    # STARTUP BANNER
    # ═══════════════════════════════════════════════════════════════════
    log_startup_banner(
        tool_name="region-mds",
        version=__version__,
        inputs={
            "BAM": str(bam_input.name),
            "Reference": str(reference.name),
            "Output": str(output),
        },
        config={
            "Genome": f"{assets.raw_genome} → {assets.genome_dir}",
            "Assay": assay or "None",
            "E1 Only": str(e1_only),
        },
        assets=[
            ResolvedAsset("Gene BED", resolved_gene_bed, gene_bed_source),
            ResolvedAsset("PON", resolved_pon_path, pon_source),
        ],
        logger=logger,
    )

    # Configure thread pool for Rayon parallelization
    if threads > 0:
        try:
            _core.configure_threads(threads)
            logger.debug(f"Configured {threads} threads for parallel processing")
        except RuntimeError:
            # Thread pool already configured (can only be set once per process)
            pass

    # Call Rust engine
    n_regions, n_genes = _core.region_mds.run_region_mds(
        str(bam_input),
        str(reference),
        str(resolved_gene_bed),
        str(output_exon),
        str(output_gene),
        e1_only,
        mapq,
        minlen,
        maxlen,
        silent,
    )

    logger.info(f"Region-MDS complete: {n_regions} regions, {n_genes} genes")

    # Set when the exon scorer has already written the final exon output, so
    # the format-conversion block below does not run a second time on a file
    # it already consumed.
    exon_written = False
    scored = False  # any gene or exon z-score written against the PON
    # What the PON z-scoring of MDS.gene did: `ran` only once its mds_z/mds_e1_z are written.
    # It carries on past its own failure, so this is the only place that failure is visible.
    z_outcome = Outcome("not_requested")

    # PON z-score normalization — read back Rust's TSV, append z-scores, re-write in selected format
    if resolved_pon_path and resolved_pon_path.exists():
        # Loaded through the guarded loader before the Rust pass (version floor included).
        logger.info(f"Applying PON z-score normalization: {resolved_pon_path}")
        pon = requested_pon
        try:
            if pon.region_mds is None:
                logger.warning(
                    "PON model does not have region_mds baseline - skipping z-scores"
                )
                z_outcome = Outcome("unavailable")
            else:
                # Exact, not parquet-first: Rust wrote this file moments ago
                # and a stale sibling from an earlier run would win (c92ed86).
                df = read_exact_table(output_gene)
                if df is None:
                    raise OSError(
                        f"could not read {output_gene.name} to add its z-scores"
                    )
                else:
                    # Compute z-scores for each gene
                    import pandas as pd

                    mds_z_scores = []
                    mds_e1_z_scores = []

                    for _, row in df.iterrows():
                        gene = row["gene"]
                        mds_mean = row.get("mds_mean", 0.0)
                        mds_e1 = row.get("mds_e1", 0.0)

                        # MDS mean z-score
                        z = pon.region_mds.compute_zscore(gene, mds_mean)
                        mds_z_scores.append(z if z is not None else float("nan"))

                        # E1 z-score
                        z_e1 = pon.region_mds.compute_e1_zscore(gene, mds_e1)
                        mds_e1_z_scores.append(
                            z_e1 if z_e1 is not None else float("nan")
                        )

                    # Add z-score columns and re-write in selected format
                    df["mds_z"] = mds_z_scores
                    df["mds_e1_z"] = mds_e1_z_scores
                    n_with_z = sum(
                        1 for z in mds_z_scores if z is not None and not pd.isna(z)
                    )
                    write_table(
                        df,
                        output_gene_base,
                        output_format=output_format,
                        compress=compress,
                    )
                    scored = n_with_z > 0  # set as soon as the scored table exists
                    z_outcome = Outcome("ran")
                    # Clean up Rust-produced gene TSV after write_table re-wrote it
                    cleanup_intermediate_tsv(output_gene, output_format, compress)
                    logger.info(f"Added z-scores for {n_with_z}/{len(df)} genes")

            # Exon-level z-scores, from the baseline added in this release.
            #
            # `MDS.exon` is the finest localisation krewlyzer produces and
            # shipped a raw score with nothing to compare it against. Handled
            # separately from the gene block above because the two baselines
            # can be present independently: a PON built before this release
            # has `region_mds` and no `region_mds_exon`.
            if pon.region_mds_exon is None:
                logger.info(
                    "PON has no region_mds_exon baseline; MDS.exon keeps its raw "
                    "score and gets no z-score. Rebuild the PON to add one."
                )
            else:
                exon_written, n_exon_z = _apply_exon_zscores(
                    pon, output_exon, output_exon_base, output_format, compress
                )
                scored = scored or n_exon_z > 0

        except Exception as e:
            logger.warning(f"PON z-score computation failed: {e}")
            z_outcome = Outcome.failed(e)

        # Outside the try, so a failure to record is not reported as a scoring failure, and a
        # partial score (gene table written, exon step failed) is still recorded.
        from krewlyzer.core.output_provenance import record_pon_used

        record_pon_used(
            output,
            sample_name,
            requested_pon if scored else None,
            resolved_pon_path,
            output_format,
            compress,
        )

    # Post-process exon output: convert format and/or compress.
    # Rust writes raw TSV; Python converts to Parquet/both/gzip as requested.
    # Skipped when the exon scorer above already wrote the final file.
    if not exon_written and (output_format != "tsv" or compress):
        logger.debug(
            f"region_mds: converting {output_exon.name} → format={output_format!r}, compress={compress}"
        )
        # Exact: Rust wrote this TSV in this run.
        df_exon = read_exact_table(output_exon)
        if df_exon is not None:
            write_table(
                df_exon,
                output_exon_base,
                output_format=output_format,
                compress=compress,
            )
            # Clean up Rust-produced exon TSV after write_table re-wrote it
            cleanup_intermediate_tsv(output_exon, output_format, compress)
        else:
            logger.warning(
                f"region_mds: could not read {output_exon.name} for format conversion"
            )

    # Post-process gene output: convert format and/or compress.
    # When PON z-scores are added (above), write_table already handles format
    # conversion. This block ensures conversion also happens when no PON is
    # provided, PON load fails, or PON lacks region_mds baseline.
    # Re-converting an already-converted file is safe, but only when the read
    # targets the file this run wrote: the z-score step above writes into
    # `output_gene`, and a resolving read would prefer a `.parquet` sibling
    # left by an earlier run and silently drop those z-scores.
    if output_format != "tsv" or compress:
        logger.debug(
            f"region_mds: converting {output_gene.name} → format={output_format!r}, compress={compress}"
        )
        df_gene = read_exact_table(output_gene)
        if df_gene is not None:
            write_table(
                df_gene,
                output_gene_base,
                output_format=output_format,
                compress=compress,
            )
            cleanup_intermediate_tsv(output_gene, output_format, compress)
        else:
            logger.warning(
                f"region_mds: could not read {output_gene.name} for format conversion"
            )

    return n_regions, n_genes, z_outcome


def region_mds(
    bam_input: Path = typer.Argument(..., help="Input BAM file (indexed)"),
    reference: Path = typer.Argument(..., help="Reference FASTA file"),
    output: Path = typer.Argument(..., help="Output directory"),
    gene_bed: Optional[Path] = typer.Option(
        None, "--gene-bed", "-g", help="Gene/exon BED file (panel or WGS format)"
    ),
    genome: str = typer.Option(
        "hg19", "--genome", "-G", help="Genome build (hg19/hg38)"
    ),
    assay: Optional[str] = typer.Option(
        None, "--assay", "-a", help="Assay code (xs1, xs2, wgs) for bundled gene BED"
    ),
    sample_name: Optional[str] = typer.Option(
        None,
        "--sample-name",
        "-s",
        help="Sample name for output files (default: derived from BAM filename)",
    ),
    e1_only: bool = typer.Option(
        False, "--e1-only", help="Only output E1 (first exon) results"
    ),
    mapq: int = typer.Option(20, "--mapq", help="Minimum mapping quality"),
    minlen: int = typer.Option(65, "--minlen", help="Minimum fragment length"),
    maxlen: int = typer.Option(
        1000,
        "--maxlen",
        help="Maximum fragment length (default: 1000 for consistency with run-all)",
    ),
    pon_model: Optional[Path] = typer.Option(
        None, "--pon-model", "-P", help="PON model for z-score computation"
    ),
    pon_variant: str = typer.Option(
        "all_unique",
        "--pon-variant",
        help="PON variant: 'all_unique' (default, max coverage) or 'duplex' (highest accuracy)",
    ),
    skip_pon: bool = typer.Option(
        False,
        "--skip-pon",
        help="Skip PON z-score normalization (for PON samples used as ML negatives)",
    ),
    threads: int = typer.Option(
        0, "--threads", "-t", help="Number of threads (0 = use all available cores)"
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable verbose logging"
    ),
    silent: bool = typer.Option(False, "--silent", help="Suppress progress bar"),
    output_format: str = typer.Option(
        "tsv",
        "--output-format",
        "-F",
        help="Output format: tsv | parquet | both (default: tsv)",
    ),
    compress: bool = typer.Option(
        False,
        "--compress",
        help="Gzip-compress TSV output (only applies when format is tsv or both)",
    ),
):
    """
    Calculate per-region Motif Diversity Score (MDS).

    MDS measures the diversity of 4-mer end motifs at each genomic region.
    Higher MDS indicates more diverse motif usage; lower MDS may indicate
    aberrant fragmentation patterns associated with disease.

    Examples:

        # Panel mode with bundled gene BED
        krewlyzer region-mds sample.bam ref.fa out/ --assay xs2

        # WGS mode
        krewlyzer region-mds sample.bam ref.fa out/ --assay wgs

        # Custom gene BED
        krewlyzer region-mds sample.bam ref.fa out/ --gene-bed custom.bed
    """
    try:
        _, _, z_outcome = region_mds_stage(
            bam_input,
            reference,
            output,
            gene_bed=gene_bed,
            genome=genome,
            assay=assay,
            sample_name=sample_name,
            e1_only=e1_only,
            mapq=mapq,
            minlen=minlen,
            maxlen=maxlen,
            pon_model=pon_model,
            pon_variant=pon_variant,
            skip_pon=skip_pon,
            threads=threads,
            verbose=verbose,
            silent=silent,
            output_format=output_format,
            compress=compress,
        )
    except (FileNotFoundError, ValueError) as e:
        # A missing input or an unusable PON, as before: a usage error.
        raise typer.BadParameter(str(e))
    except Exception as e:
        logger.error(f"Region-MDS failed: {e}")
        raise typer.Exit(1)
    # Into the sample's record, if it has one: rerunning this command is how a Region-MDS that
    # failed under run-all gets fixed. Success only -- a failure here exits non-zero, and may
    # not have touched the tables a good earlier run left.
    record_stage_outcomes(
        output,
        sample_name or _derived_sample_name(bam_input),
        {"region_mds": Outcome("ran"), "region_mds_z": z_outcome},
        output_format=output_format,
        compress=compress,
        owns_record=False,
    )
