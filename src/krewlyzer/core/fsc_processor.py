"""
FSC (Fragment Size Coverage) processor.

Processes GC-weighted fragment counts into ML-ready features per the implementation plan:
- Raw GC-weighted counts per channel (from Rust backend)
- PoN log2 ratios: log2(channel / PoN_mean)
- Reliability scores: 1 / (PoN_variance + k)

NO z-scores - the plan specifies GC-weighting (Rust) + PoN log-ratios (Python).
"""

from pathlib import Path
from typing import Optional, Union
import pandas as pd
import numpy as np
import logging

from .output_utils import (
    read_exact_table,
    write_table,
    written_table,
    cleanup_intermediate_tsv,
)

logger = logging.getLogger("core.fsc_processor")

# 6 biologically-meaningful channels matching Rust backend
# Extended to include ultra_long for fragments 401-1000bp
CHANNELS = ["ultra_short", "core_short", "mono_nucl", "di_nucl", "long", "ultra_long"]

# Regularization constant for reliability calculation
RELIABILITY_K = 0.01

#: 100 kb bins per FSC/FSR row. Genome-wide, 50 bins make 5 Mb, arm-scale windows. A panel keeps
#: one bin per row: its targets are sparse, so a 5 Mb window would merge unrelated genes.
WGS_CONTINUE_N = 50
PANEL_CONTINUE_N = 1


def default_continue_n(panel_mode: bool) -> int:
    """Bins per FSC/FSR row when none is asked for: run-all's rule, and the single tools'."""
    return PANEL_CONTINUE_N if panel_mode else WGS_CONTINUE_N


def add_pon_channel_columns(
    frame: "pd.DataFrame", pon, ontarget: bool = False
) -> "pd.DataFrame":
    """Add ``{channel}_log2`` and ``{channel}_reliability`` in place.

    Split out of ``process_fsc`` so it can be tested without a fragment BED.
    A branch only reachable through a full pipeline run is a branch nothing
    checks, and this one silently fabricated a value for years.

    Both fall back to NaN, never to a literal:

    ``log2 = 0.0`` says "this sample sits exactly at the healthy baseline" --
    the most confident statement the column can make, asserted precisely when
    there is no baseline to compare against. Measured with ``gc_bias`` absent:
    three windows whose raw values differed fourfold all read 0.0.

    ``reliability = 1.0`` is not neutral either. With ``RELIABILITY_K`` at
    0.01 it is what a PON variance of 0.99 would produce -- a specific and
    unremarkable reading, indistinguishable from a measured one.

    The PON's ``gc_bias`` block reaches here through ``PonModel.get_mean``, so
    a model missing that block used to produce a confident zero in every window
    of every sample rather than a visibly absent column.

    ``ontarget`` selects ``gc_bias_ontarget``. Capture enrichment gives
    on-target fragments a different GC profile, which is why panel mode fits a
    separate block -- one that was built and stored by every panel PON and read
    by nothing, so on-target FSC normalised against the genome-wide curves.
    """
    for ch in CHANNELS:
        mean = pon.get_mean(ch, ontarget=ontarget) if hasattr(pon, "get_mean") else None
        var = (
            pon.get_variance(ch, ontarget=ontarget)
            if hasattr(pon, "get_variance")
            else None
        )

        if mean and mean > 0:
            frame[f"{ch}_log2"] = np.log2(frame[ch] / mean + 1e-9)
        else:
            frame[f"{ch}_log2"] = np.nan

        if var is not None:
            frame[f"{ch}_reliability"] = 1.0 / (var + RELIABILITY_K)
        else:
            frame[f"{ch}_reliability"] = np.nan
    return frame


def process_fsc(
    counts_df: pd.DataFrame,
    output_path: Path,
    windows: int = 100000,
    continue_n: int = 50,
    pon=None,
    output_format: str = "tsv",
    compress: bool = False,
    ontarget: bool = False,
) -> Path:
    """
    Process GC-weighted fragment counts into FSC ML features.

    Per the implementation plan, FSC normalization order is:
    1. GC-weight (Rust) - done upstream
    2. PoN log-ratio (Python) - done here

    Output columns:
    - chrom, start, end: Bin coordinates
    - ultra_short, core_short, mono_nucl, di_nucl, long: GC-weighted counts
    - total: All fragments
    - *_log2: log2(channel / PoN_mean) when PoN provided
    - *_reliability: 1 / (PoN_variance + k) when PoN provided

    Args:
        counts_df: DataFrame from Rust with [chrom, start, end, ultra_short,
                   core_short, mono_nucl, di_nucl, long, total, mean_gc]
        output_path: Path to write FSC.tsv output
        windows: Window size in bp (default 100000)
        continue_n: Number of bins to aggregate (default 50)
        pon: Optional PON model for log2 ratio normalization

    Returns:
        Path to the output file
    """
    logger.info(f"Processing FSC: {len(counts_df)} bins → {output_path}")

    # Validate required columns
    required_cols = ["chrom", "start", "end"] + CHANNELS + ["total"]
    missing = [c for c in required_cols if c not in counts_df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Aggregate bins into windows
    results = []
    for chrom, group in counts_df.groupby("chrom", sort=False):
        n_bins = len(group)
        n_windows = n_bins // continue_n

        if n_windows == 0:
            continue

        trunc_len = n_windows * continue_n

        # Aggregate each channel
        channel_sums = {}
        for ch in CHANNELS:
            mat = group[ch].to_numpy()[:trunc_len].reshape(n_windows, continue_n)
            channel_sums[ch] = mat.sum(axis=1)

        totals_mat = (
            group["total"].to_numpy()[:trunc_len].reshape(n_windows, continue_n)
        )
        channel_sums["total"] = totals_mat.sum(axis=1)

        # Get window coordinates from first and last bin in each window
        starts = (
            group["start"].to_numpy()[:trunc_len].reshape(n_windows, continue_n)[:, 0]
        )
        ends = group["end"].to_numpy()[:trunc_len].reshape(n_windows, continue_n)[:, -1]

        window_df = pd.DataFrame(
            {"chrom": chrom, "start": starts, "end": ends, **channel_sums}
        )
        results.append(window_df)

    if not results:
        logger.warning("No valid windows found for FSC")
        _write_empty_fsc(
            output_path, pon is not None, output_format=output_format, compress=compress
        )
        return output_path

    final_df = pd.concat(results, ignore_index=True)
    logger.info(f"Aggregated into {len(final_df)} windows")

    # Compute PoN log2 ratios and reliability if model provided
    if pon is not None:
        logger.info("Computing PoN log2 ratios...")
        add_pon_channel_columns(final_df, pon, ontarget=ontarget)

    # Write output
    _write_fsc_output(
        final_df,
        output_path,
        pon is not None,
        output_format=output_format,
        compress=compress,
    )

    logger.info(f"FSC complete: {output_path}")
    return output_path


def _write_empty_fsc(
    output_path: Path,
    with_pon: bool,
    output_format: str = "tsv",
    compress: bool = False,
):
    """Write empty FSC file with appropriate headers."""
    cols = ["chrom", "start", "end"] + CHANNELS + ["total"]
    if with_pon:
        cols += [f"{ch}_log2" for ch in CHANNELS]
        cols += [f"{ch}_reliability" for ch in CHANNELS]
    write_table(
        pd.DataFrame(columns=cols),
        output_path,
        output_format=output_format,
        compress=compress,
    )


def _write_fsc_output(
    df: pd.DataFrame,
    output_path: Path,
    with_pon: bool,
    output_format: str = "tsv",
    compress: bool = False,
):
    """Write FSC output matching the implementation plan schema."""
    # Determine columns to output
    cols = ["chrom", "start", "end"] + CHANNELS + ["total"]
    if with_pon:
        cols += [f"{ch}_log2" for ch in CHANNELS]
        cols += [f"{ch}_reliability" for ch in CHANNELS]

    # Select and write
    output_df = df[cols].copy()
    write_table(
        output_df,
        output_path,
        output_format=output_format,
        compress=compress,
        float_format="%.4f",
    )


def aggregate_by_gene(
    fragments_bed: Path,
    genes: dict,
    output_path: Path,
    pon=None,
    correction_factors: Optional[Union[str, Path]] = None,
    aggregate_by: str = "gene",
    gene_bed_path: Optional[Path] = None,
    output_format: str = "tsv",
    compress: bool = False,
) -> Path:
    """
    Aggregate fragment counts for panel-specific FSC using Rust implementation.

    This function counts fragments per gene using the high-performance Rust engine,
    applying GC correction if factors are available.

    Args:
        fragments_bed: Path to fragment BED file from extract step
        genes: Dict[gene_name, List[Region]] from gene_bed.load_gene_bed() (used only if gene_bed_path not provided)
        output_path: Path to write FSC output (base name; extension added by write_table)
        pon: Optional PON model for normalization (legacy, not used in Rust)
        correction_factors: Optional path to a six-column GC factor table (any format Rust
            loads: .parquet, .tsv.gz, .tsv). Resolve it with core.gc_assets.resolve_gc_factors.
        aggregate_by: 'gene' (sum by gene) or 'region' (one row per exon/region)
        gene_bed_path: Path to gene BED file (if provided, used directly by Rust)
        output_format: "tsv", "parquet", or "both" (default: "tsv")
        compress: Gzip TSV output (only when output_format is tsv or both)

    Returns:
        Path to output file
    """
    from krewlyzer import _core
    import tempfile  # Used for the temp gene BED

    logger.info(f"Aggregating FSC by {aggregate_by}: {len(genes)} genes")

    # Prepare gene BED for Rust
    if gene_bed_path and Path(gene_bed_path).exists():
        use_gene_bed = Path(gene_bed_path)
    else:
        # Write genes dict to temp BED for Rust.
        #
        # Eight columns, not five. The extra three carry strand and the two E1
        # flags through to the Rust aggregator; a 5-column file dropped them
        # here, which is why FSC could not do strand-aware E1 no matter what
        # the asset contained. `.` marks unknown -- a legacy gene BED or a
        # target BED parsed on the fly has no flags, and that must stay
        # distinguishable from a genuine "false".
        def _flag(value: Optional[bool]) -> str:
            return "." if value is None else ("1" if value else "0")

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            for gene, regions in genes.items():
                for region in regions:
                    name = getattr(region, "name", gene)
                    strand = getattr(region, "strand", None) or "."
                    f.write(
                        f"{region.chrom}\t{region.start}\t{region.end}\t{gene}\t"
                        f"{name}\t{strand}\t"
                        f"{_flag(getattr(region, 'is_e1', None))}\t"
                        f"{_flag(getattr(region, 'is_alt_e1', None))}\n"
                    )
            use_gene_bed = Path(f.name)

    # GC factors: a path Rust loads itself. The dict this replaced was written back out as a
    # three-column TSV the Rust loader skips row by row, so gene FSC was never corrected.
    if correction_factors is not None and not isinstance(
        correction_factors, (str, Path)
    ):
        raise TypeError(
            "correction_factors must be a path to a GC factor table, "
            f"not {type(correction_factors).__name__}"
        )
    gc_factors_path = str(correction_factors) if correction_factors else None

    # Rust always writes TSV to output_path; we convert afterward if needed.
    tsv_path = Path(
        str(output_path)
        if str(output_path).endswith(".tsv")
        else str(output_path) + ".tsv"
    )
    n_output = _core.aggregate_by_gene(
        str(fragments_bed),
        str(use_gene_bed),
        str(tsv_path),
        gc_factors_path,
        aggregate_by,
    )

    # Cleanup temp files if we created them
    if not gene_bed_path or not Path(gene_bed_path).exists():
        try:
            use_gene_bed.unlink(missing_ok=True)
        except Exception:
            pass

    logger.info(f"FSC {aggregate_by} complete: {tsv_path} ({n_output} rows)")

    # Post-Rust format conversion: read TSV, write in selected format
    base_path = Path(str(output_path).removesuffix(".tsv"))
    converted = False
    if output_format != "tsv" or compress:
        df = read_exact_table(tsv_path)
        if df is None:
            # Do not fall through quietly: the caller would then be handed a
            # path for a format that was never produced. Keep the .tsv Rust
            # wrote and say so.
            logger.warning(
                f"FSC {aggregate_by}: could not read {tsv_path.name} back for "
                f"conversion to format={output_format}, compress={compress}; "
                f"leaving the uncompressed TSV in place"
            )
        else:
            write_table(
                df,
                base_path,
                output_format=output_format,
                compress=compress,
                float_format="%.4f",
            )
            # Clean up original Rust TSV now that write_table() produced target format
            cleanup_intermediate_tsv(tsv_path, output_format, compress)
            converted = True
            logger.debug(
                f"FSC {aggregate_by}: converted to format={output_format}, compress={compress}"
            )

    # Return what was actually written. Returning the pre-conversion .tsv meant
    # callers guarding on .exists() saw a file cleanup_intermediate_tsv had just
    # deleted, so E1 generation and FSC PON z-scoring were skipped on every
    # compressed or Parquet run. parent/(name+ext) rather than with_suffix(),
    # which would eat the compound extension.
    if not converted:
        return tsv_path
    written_ext = (
        ".parquet"
        if output_format == "parquet"
        else (".tsv.gz" if compress else ".tsv")
    )
    return base_path.parent / (base_path.name + written_ext)


def apply_fsc_gene_pon(
    fsc_gene_path: Path,
    pon,
    output_path: Optional[Path] = None,
    *,
    output_format: str,
    compress: bool,
) -> Optional[Path]:
    """
    Apply FSC gene PON baseline z-score normalization.

    Reads FSC.gene.tsv (or .parquet), computes z-score for each gene's
    normalized_depth against the PON baseline, and adds a depth_zscore column.

    Args:
        fsc_gene_path: Path to FSC.gene.tsv file
        pon: PonModel with fsc_gene_baseline
        output_path: Output path (defaults to overwriting input)
        output_format: "tsv", "parquet", or "both"
        compress: Gzip TSV output when True

    Returns:
        Path to output file with depth_zscore column

    Note:
        Clinical interpretation:
        - z-score >> 0: Gene amplification
        - z-score << 0: Gene deletion
        - z-score ≈ 0: Normal copy number
    """
    if output_path is None:
        output_path = fsc_gene_path

    # The table this run wrote, exactly. Callers name the `.tsv`, but the writer
    # honours --output-format and may have written only `.parquet`. A bare
    # `.exists()` was False there, so this returned None and no `depth_zscore`
    # column was ever added under `--output-format parquet` -- the format 0.9.0
    # makes the Nextflow default. Any format that exists is no better: it can be
    # an earlier run's.
    source = written_table(
        fsc_gene_path, output_format=output_format, compress=compress
    )
    if not source.exists():
        logger.warning(f"FSC gene file not found: {source}")
        return None

    if not pon or not pon.fsc_gene_baseline:
        logger.debug("No FSC gene baseline in PON, skipping z-score")
        return fsc_gene_path

    df = read_exact_table(source)
    if df is None:
        logger.warning(f"FSC gene file could not be read: {source}")
        return fsc_gene_path

    if "gene" not in df.columns or "normalized_depth" not in df.columns:
        logger.warning(f"FSC gene file missing required columns: {fsc_gene_path}")
        return fsc_gene_path

    # Compute z-scores using PON baseline
    baseline = pon.fsc_gene_baseline
    z_scores = []
    for _, row in df.iterrows():
        gene = row["gene"]
        depth = row["normalized_depth"]
        z = baseline.compute_zscore(gene, depth)
        z_scores.append(z if z is not None else float("nan"))

    df["depth_zscore"] = z_scores

    # Write output using unified writer (TSV/Parquet/both)
    out_base = Path(output_path).with_suffix("")
    write_table(df, out_base, output_format=output_format, compress=compress)
    n_with_zscore = df["depth_zscore"].notna().sum()
    logger.debug(f"Applied FSC gene PON: {n_with_zscore}/{len(df)} genes with z-scores")

    return output_path


def apply_fsc_region_pon(
    fsc_region_path: Path,
    pon,
    output_path: Optional[Path] = None,
    *,
    output_format: str,
    compress: bool,
) -> Optional[Path]:
    """
    Apply FSC region (exon/probe) PON baseline z-score normalization.

    Reads FSC.regions.tsv (or .parquet), computes z-score for each region's
    normalized_depth against the PON baseline, adds a depth_zscore column.

    Args:
        fsc_region_path: Path to FSC.regions.tsv file
        pon: PonModel with fsc_region_baseline
        output_path: Output path (defaults to overwriting input)
        output_format: "tsv", "parquet", or "both"
        compress: Gzip TSV output when True

    Returns:
        Path to output file with depth_zscore column

    Note:
        Region IDs are formatted as "chrom:start-end".
        Clinical interpretation:
        - z-score << 0: Exon deletion
        - z-score >> 0: Exon amplification
    """
    if output_path is None:
        output_path = fsc_region_path

    # Same as apply_fsc_gene_pon above: exactly the table this run wrote.
    source = written_table(
        fsc_region_path, output_format=output_format, compress=compress
    )
    if not source.exists():
        logger.warning(f"FSC region file not found: {source}")
        return None

    if not pon or not pon.fsc_region_baseline:
        logger.debug("No FSC region baseline in PON, skipping z-score")
        return fsc_region_path

    df = read_exact_table(source)
    if df is None:
        logger.warning(f"FSC region file could not be read: {source}")
        return fsc_region_path

    required_cols = ["chrom", "start", "end", "normalized_depth"]
    if not all(c in df.columns for c in required_cols):
        logger.warning(f"FSC region file missing required columns: {fsc_region_path}")
        return fsc_region_path

    # Create region_id for lookup
    df["region_id"] = (
        df["chrom"].astype(str)
        + ":"
        + df["start"].astype(str)
        + "-"
        + df["end"].astype(str)
    )

    # Compute z-scores using PON baseline
    baseline = pon.fsc_region_baseline
    z_scores = []
    for _, row in df.iterrows():
        region_id = row["region_id"]
        depth = row["normalized_depth"]
        z = baseline.compute_zscore(region_id, depth)
        z_scores.append(z if z is not None else float("nan"))

    df["depth_zscore"] = z_scores

    # Drop temporary column, write output using unified writer
    df = df.drop(columns=["region_id"])
    out_base = Path(output_path).with_suffix("")
    write_table(df, out_base, output_format=output_format, compress=compress)
    n_with_zscore = df["depth_zscore"].notna().sum()
    logger.debug(
        f"Applied FSC region PON: {n_with_zscore}/{len(df)} regions with z-scores"
    )

    return output_path


# ============================================================================
# E1-ONLY FILTERING
# ============================================================================
#
# First exon (E1) filtering extracts only promoter-proximal regions.
# Per Helzer et al. (2025), promoters are Nucleosome Depleted Regions (NDRs)
# with distinct fragmentation patterns. E1 signal is often stronger for
# early cancer detection than whole-gene averages.
#
# E1 is identified as the first exon by genomic position per gene,
# consistent with the region-mds implementation.
# ============================================================================


def filter_fsc_to_e1(
    fsc_regions_path: Path,
    output_path: Optional[Path] = None,
    *,
    output_format: str,
    compress: bool,
) -> Optional[Path]:
    """
    Filter FSC regions to E1 (first exon) per gene.

    Extracts the first exon by genomic position for each gene from the
    full FSC.regions.tsv (or .parquet) file. This focuses on
    promoter-proximal regions where cancer fragmentation signal is strongest.

    Args:
        fsc_regions_path: Path to FSC.regions.tsv or .parquet file
        output_path: Output base path (extension added by write_table)
        output_format: "tsv", "parquet", or "both"
        compress: Gzip TSV output when True

    Returns:
        Path to E1-only output file, or None if input doesn't exist

    Scientific context:
        E1 (first exon) serves as a proxy for promoter regions which are
        Nucleosome Depleted Regions (NDRs). Differential enzymatic access
        at NDRs produces distinct fragment end motifs. Cancer signal is
        often stronger at E1 than averaging across all exons.
    """
    fsc_regions_path = Path(fsc_regions_path)

    # Same again: exactly the table this run wrote, rather than assuming the `.tsv`.
    regions_file = written_table(
        fsc_regions_path, output_format=output_format, compress=compress
    )
    if not regions_file.exists():
        logger.warning(f"FSC regions file not found: {regions_file}")
        return None

    # Default output base path: strip suffix so write_table appends the right one
    if output_path is None:
        # NB: Path.stem strips only the LAST dot-segment, so for
        # 'S1.FSC.regions.tsv.gz' it yields 'S1.FSC.regions.tsv' and the
        # '.regions' test below misses -- producing the mangled name
        # 'S1.FSC.regions.tsv.e1only.tsv.gz'. Strip compound extensions.
        from .output_utils import strip_table_extension

        stem = strip_table_extension(fsc_regions_path.name)
        if stem.endswith(".regions"):
            stem = stem.replace(".regions", ".regions.e1only")
        else:
            stem = f"{stem}.e1only"
        output_base = fsc_regions_path.parent / stem
    else:
        # Strip suffix so write_table controls the final extension
        output_base = Path(output_path).with_suffix("")

    df = read_exact_table(regions_file)
    if df is None:
        logger.warning(f"FSC regions file not found or unreadable: {regions_file}")
        return None

    # Validate required columns
    required_cols = ["gene", "chrom", "start"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        logger.warning(
            f"FSC regions file missing required columns {missing}: {fsc_regions_path}"
        )
        return None

    if df.empty:
        logger.warning(f"FSC regions file is empty: {fsc_regions_path}")
        return None

    n_input = len(df)
    n_genes = df["gene"].nunique()

    # E1 selection: prefer the precomputed flags, fall back to coordinates.
    #
    # A region qualifies if it overlaps the canonical transcript's exon 1
    # (`is_e1`) **or** another basic protein-coding transcript's first exon
    # (`is_alt_e1`). Both are genuine transcription starts, and on a targeted
    # panel restricting to the canonical one throws most of the signal away:
    # for xs1, 25 of 128 genes have a canonical-E1 tile while 40 have a tile on
    # some basic protein-coding first exon. Alternative promoters are the norm
    # -- genes carry a median of 13 distinct annotated first exons.
    #
    # Minor and non-coding isoforms are deliberately excluded upstream, in
    # `scripts/build_gene_bed.py`: their annotated starts are not evidence of a
    # live promoter, and counting them would inflate the set to 79 of 128 with
    # regions that are not promoter-proximal in any real sense.
    #
    # Genes with neither flag are **dropped**, not back-filled with their most
    # 5' region. That region is usually an internal exon, and labelling it E1
    # asserts promoter proximity that is not there -- the defect this replaces.
    flag_cols = [c for c in ("is_e1", "is_alt_e1") if c in df.columns]
    flagged = (
        df[df[flag_cols].isin(["1", 1, True]).any(axis=1)]
        if flag_cols
        else df.iloc[0:0]
    )

    if not flagged.empty:
        source = flagged
        n_kept = flagged["gene"].nunique()
        logger.info(
            "E1 filter: %d of %d gene(s) have a region on an annotated first "
            "exon; %d gene(s) have none and are omitted rather than back-filled "
            "with an internal exon",
            n_kept,
            n_genes,
            n_genes - n_kept,
        )
    else:
        # Fallback: lowest start per gene. Retained so a legacy or custom gene
        # BED still produces a table, but it is a positional guess, not an E1.
        source = df
        if flag_cols:
            logger.warning(
                "E1 filter: %s carries E1 flags but no region is marked. "
                "Falling back to the lowest-start region per gene, which "
                "ignores strand and is wrong for minus-strand genes.",
                fsc_regions_path.name,
            )
        else:
            logger.warning(
                "E1 filter: no is_e1/is_alt_e1 columns in %s, so E1 is the "
                "lowest-start region per gene. That ignores strand and picks "
                "the LAST exon of every minus-strand gene. Regenerate the gene "
                "BED with scripts/build_gene_bed.py to fix this.",
                fsc_regions_path.name,
            )

    df_sorted = source.sort_values(["gene", "chrom", "start"])
    e1_df = df_sorted.groupby("gene", as_index=False).first()

    # Preserve original column order, write output
    original_cols = [c for c in df.columns if c in e1_df.columns]
    e1_df = e1_df[original_cols]
    write_table(e1_df, output_base, output_format=output_format, compress=compress)

    # Reconstruct return path with expected extension for caller.
    # IMPORTANT: do NOT use output_base.with_suffix(ext) — Python's with_suffix()
    # replaces only the last dot-segment, so 'test.FSC.regions.e1only'.with_suffix('.tsv')
    # gives 'test.FSC.regions.tsv' (replaces '.e1only' with '.tsv').
    # Instead, build the path by appending the extension to the full stem string.
    ext = (
        ".tsv.gz"
        if compress
        else (".parquet" if output_format == "parquet" else ".tsv")
    )
    result_path = output_base.parent / (output_base.name + ext)

    # `n_genes` is the count *before* filtering; reporting it beside the output
    # row count reads as though every gene survived, which is now usually false.
    logger.info(
        "FSC E1-only filter: %d regions across %d genes → %d E1 region(s) "
        "across %d gene(s) → %s",
        n_input,
        n_genes,
        len(e1_df),
        e1_df["gene"].nunique(),
        result_path.name,
    )

    return result_path
