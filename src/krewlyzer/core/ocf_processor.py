"""
OCF (Open Chromatin Footprinting) processor.

Handles OCF post-processing:
- PON z-score normalization via Rust backend
- Output format conversion (TSV/Parquet/gzip) via write_table()

After Rust writes raw TSV, this module re-writes through write_table()
to honour --output-format and --compress flags.
"""

from pathlib import Path
from typing import Optional
import logging

from .output_utils import read_table, write_table, cleanup_intermediate_tsv

logger = logging.getLogger("core.ocf_processor")


def process_ocf_with_pon(
    ocf_path: Path,
    pon_parquet_path: Path,
    output_path: Optional[Path] = None,
    output_format: str = "tsv",
    compress: bool = False,
) -> int:
    """
    Add PON z-scores to OCF output using Rust implementation,
    then re-write through write_table() for format/compression support.

    Args:
        ocf_path: Path to sample OCF TSV (Rust-written)
        pon_parquet_path: Path to PON Parquet file
        output_path: Output path (default: overwrite input)
        output_format: One of "tsv", "parquet", or "both"
        compress: If True, gzip-compress TSV output (.tsv.gz)

    Returns:
        Number of regions with z-scores computed

    Raises:
        RuntimeError: If Rust implementation fails
    """
    from krewlyzer import _core

    if not ocf_path.exists():
        logger.warning(f"OCF file not found: {ocf_path}")
        return 0

    out_str = str(output_path) if output_path else None

    n_matched = _core.ocf.apply_pon_zscore(
        str(ocf_path), str(pon_parquet_path), out_str
    )

    logger.info(f"OCF PON z-scores: {n_matched} regions matched ({ocf_path.name})")

    # Re-write through write_table() to honour output_format and compress.
    # Rust wrote the result as plain TSV; now produce .parquet/.tsv.gz as needed.
    final_path = output_path or ocf_path
    _write_ocf_output(final_path, output_format, compress)

    return n_matched


def convert_ocf_output(
    ocf_path: Path,
    output_format: str = "tsv",
    compress: bool = False,
) -> None:
    """Convert a Rust-written OCF TSV to the requested output format.

    Used for OCF files that do NOT go through PON normalization
    (e.g., when --skip-pon is active or no PON is provided).

    Args:
        ocf_path: Path to the Rust-written OCF TSV file
        output_format: One of "tsv", "parquet", or "both"
        compress: If True, gzip-compress TSV output (.tsv.gz)
    """
    _write_ocf_output(ocf_path, output_format, compress)


def _write_ocf_output(tsv_path: Path, output_format: str, compress: bool) -> None:
    """Read Rust-written OCF TSV and re-write through write_table().

    Shared helper used by both process_ocf_with_pon() and convert_ocf_output()
    to avoid code duplication.

    Args:
        tsv_path: Path to the Rust-written TSV file
        output_format: One of "tsv", "parquet", or "both"
        compress: If True, gzip-compress TSV output
    """
    df = read_table(tsv_path)
    if df is None:
        logger.warning(f"OCF output not found for format conversion: {tsv_path}")
        return

    logger.debug(
        f"OCF format conversion: {len(df)} rows, "
        f"output_format={output_format}, compress={compress}"
    )

    write_table(df, tsv_path, output_format=output_format, compress=compress)

    # Clean up intermediate Rust TSV if write_table() produced a different format
    cleanup_intermediate_tsv(tsv_path, output_format, compress)


def apply_ocf_python_pon(
    ocf_path: Path,
    pon,
    baseline_attr: str = "ocf_baseline_ontarget",
    output_format: str = "tsv",
    compress: bool = False,
) -> int:
    """Apply PON z-scores to OCF output using Python (for on-target/off-target).

    Unlike primary OCF which uses Rust-based z-score computation, on-target and
    off-target OCF files are small (<50 rows) so Python-side computation is
    efficient and avoids Rust changes.

    Args:
        ocf_path: Path to on-target or off-target OCF TSV
        pon: Loaded PonModel with on-target/off-target OCF baseline
        baseline_attr: PON attribute name ('ocf_baseline_ontarget' or 'ocf_baseline_offtarget')
        output_format: One of "tsv", "parquet", or "both"
        compress: If True, gzip-compress TSV output

    Returns:
        Number of regions with z-scores matched
    """
    import pandas as pd

    baseline = getattr(pon, baseline_attr, None) if pon else None
    if baseline is None:
        logger.debug(
            f"No {baseline_attr} in PON, skipping OCF z-scores for {ocf_path.name}"
        )
        convert_ocf_output(ocf_path, output_format, compress)
        return 0

    df = read_table(ocf_path)
    if df is None:
        logger.warning(f"OCF file not found for Python PON: {ocf_path}")
        return 0

    # Determine column names (Rust outputs "tissue"/"OCF", PON uses "region_id"/"ocf")
    region_col = "tissue" if "tissue" in df.columns else "region_id"
    ocf_col = "OCF" if "OCF" in df.columns else "ocf"

    if region_col not in df.columns or ocf_col not in df.columns:
        logger.warning(
            f"OCF file missing required columns ({region_col}, {ocf_col}): {ocf_path.name}"
        )
        convert_ocf_output(ocf_path, output_format, compress)
        return 0

    # Two different absences, kept apart.
    #
    # `compute_zscore` returns None when the region is not in the baseline at
    # all, and NaN when it is there but the baseline measured no usable spread.
    # Both write NaN to the column -- correctly, since neither is a reading --
    # but collapsing them in the log would report "0/1 matched" for a region
    # that matched perfectly well and simply cannot be scored. That is the
    # difference between "rebuild with this assay" and "this region has no
    # variance in the cohort", and only the log can tell them apart.
    z_scores = []
    n_absent = 0
    n_unscoreable = 0
    for _, row in df.iterrows():
        z = baseline.compute_zscore(str(row[region_col]), float(row[ocf_col]))
        if z is None:
            n_absent += 1
            z_scores.append(float("nan"))
            continue
        if pd.isna(z):
            n_unscoreable += 1
        z_scores.append(z)

    df["ocf_z"] = z_scores
    n_matched = len(df) - n_absent
    n_scored = n_matched - n_unscoreable

    write_table(df, ocf_path, output_format=output_format, compress=compress)
    # Same cleanup the genome-wide path does via `_write_ocf_output`.
    #
    # Without it a `--output-format parquet` run left `{sample}.OCF.ontarget.tsv`
    # and `.offtarget.tsv` beside their Parquet, while the genome-wide table
    # correctly had none. Two files for one table, differing in age, and no way
    # for a reader to tell which is current -- the ambiguity `no_collided_columns`
    # exists to catch one level down.
    cleanup_intermediate_tsv(ocf_path, output_format, compress)
    logger.info(
        f"OCF Python PON ({baseline_attr}): {n_matched}/{len(df)} regions matched, "
        f"{n_scored} scored ({ocf_path.name})"
    )
    if n_absent:
        logger.warning(
            f"OCF Python PON ({baseline_attr}): {n_absent}/{len(df)} regions are "
            "absent from the baseline and get no z-score. Rebuild the PON for "
            "this assay if they should be there."
        )
    if n_unscoreable:
        logger.warning(
            f"OCF Python PON ({baseline_attr}): {n_unscoreable}/{len(df)} regions "
            "matched the baseline but it measured no usable spread for them, so "
            "the z-score is NaN rather than a value divided by a placeholder."
        )
    return n_matched
