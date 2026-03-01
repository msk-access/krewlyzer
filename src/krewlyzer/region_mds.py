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
from . import __version__


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
    from krewlyzer import _core
    from krewlyzer.assets import AssetManager

    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate inputs
    if not bam_input.exists():
        raise typer.BadParameter(f"BAM file not found: {bam_input}")
    if not reference.exists():
        raise typer.BadParameter(f"Reference FASTA not found: {reference}")

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
            raise typer.BadParameter(f"Gene BED not found: {gene_bed}")
        resolved_gene_bed = gene_bed
        gene_bed_source = "explicit"
        logger.debug(f"Using user-provided gene BED: {gene_bed}")
    elif assay:
        # Bundled gene BED from assay
        resolved_gene_bed = assets.get_gene_bed_for_mode(assay)
        if resolved_gene_bed is None:
            raise typer.BadParameter(
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
            raise typer.BadParameter(
                "No gene BED specified. Use --gene-bed, --assay, or ensure wgs.genes.bed.gz is available."
            )

    # ═══════════════════════════════════════════════════════════════════
    # ASSET RESOLUTION (PON)
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
    except ValueError as e:
        raise typer.BadParameter(str(e))

    # Create output directory
    output.mkdir(parents=True, exist_ok=True)

    # Derive sample name: use CLI-provided value, fallback to BAM stem.
    # When called via run-all, wrapper.py passes the consistent sample ID;
    # standalone CLI can rely on the BAM filename fallback.
    if not sample_name:
        sample_name = bam_input.stem
        if sample_name.endswith(".sorted"):
            sample_name = sample_name[:-7]
        logger.debug(f"Derived sample name from BAM: {sample_name}")
    else:
        logger.debug(f"Using provided sample name: {sample_name}")

    output_exon = output / f"{sample_name}.MDS.exon.tsv"
    output_gene = output / f"{sample_name}.MDS.gene.tsv"

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

    # PON z-score normalization
    if resolved_pon_path and resolved_pon_path.exists():
        try:
            import pandas as pd
            from krewlyzer.pon.model import PonModel

            logger.info(f"Applying PON z-score normalization: {resolved_pon_path}")
            pon = PonModel.load(resolved_pon_path)

            if pon.region_mds is None:
                logger.warning(
                    "PON model does not have region_mds baseline - skipping z-scores"
                )
            else:
                # Read gene output
                df = pd.read_csv(output_gene, sep="\t")

                # Compute z-scores for each gene
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
                    mds_e1_z_scores.append(z_e1 if z_e1 is not None else float("nan"))

                # Add z-score columns
                df["mds_z"] = mds_z_scores
                df["mds_e1_z"] = mds_e1_z_scores

                # Write updated output
                df.to_csv(output_gene, sep="\t", index=False)

                n_with_z = sum(
                    1 for z in mds_z_scores if z is not None and not pd.isna(z)
                )
                logger.info(f"Added z-scores for {n_with_z}/{len(df)} genes")

        except Exception as e:
            logger.warning(f"PON z-score computation failed: {e}")

    return n_regions, n_genes
