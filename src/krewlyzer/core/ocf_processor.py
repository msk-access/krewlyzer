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

from .output_utils import read_table, write_table

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

    # If output_format is "parquet" only, clean up the intermediate Rust TSV
    if output_format == "parquet" and tsv_path.exists():
        tsv_path.unlink()
        logger.debug(f"Removed intermediate TSV: {tsv_path.name}")
