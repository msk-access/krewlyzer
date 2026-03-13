"""
FSD (Fragment Size Distribution) processor.

Provides FSD-specific processing including PoN log-ratio normalization
and output format conversion (TSV/Parquet/gzip).

Used by both standalone fsd.py and run-all wrapper.py.

Normalization order:
1. GC-weighting (Rust) - raw counts are GC-corrected
2. PoN log-ratio (Rust) - log2(sample / PoN_expected) via apply_pon_logratio
3. Output format conversion (Python) - write_table() handles TSV/Parquet/gzip
"""

from pathlib import Path
from typing import Optional
import logging

from .output_utils import read_table, write_table

logger = logging.getLogger("core.fsd_processor")


def process_fsd(
    fsd_raw_path: Path,
    output_path: Optional[Path] = None,
    pon_parquet_path: Optional[Path] = None,
    output_format: str = "tsv",
    compress: bool = False,
) -> Path:
    """
    Process raw FSD counts into ML-ready features.

    Converts raw GC-weighted counts from Rust into:
    - Log-ratios vs PoN (when PoN provided): log2(sample / PoN_expected)
    - PoN stability scores: 1 / (variance + k)

    After Rust processing, the result is written through write_table()
    to honour --output-format and --compress flags (TSV, Parquet, gzip).

    Args:
        fsd_raw_path: Path to raw FSD.tsv from Rust
        output_path: Output path (default: overwrite input).
            write_table() strips known extensions and appends the correct one(s).
        pon_parquet_path: Path to PON parquet for normalization
        output_format: One of "tsv", "parquet", or "both"
        compress: If True, gzip-compress TSV output (.tsv.gz)

    Returns:
        Path to processed output file (base path; actual files may have
        .tsv, .tsv.gz, or .parquet extensions depending on output_format)

    Raises:
        RuntimeError: If FSD processing fails
    """
    if not fsd_raw_path.exists():
        raise FileNotFoundError(f"FSD output not found: {fsd_raw_path}")

    output_path = output_path or fsd_raw_path

    # Apply PON log-ratio normalization via Rust
    # Rust writes the result as TSV to output_path (or in-place if same as input)
    if pon_parquet_path and pon_parquet_path.exists():
        try:
            from krewlyzer import _core

            arms_processed = _core.fsd.apply_pon_logratio(
                str(fsd_raw_path),
                str(pon_parquet_path),
                str(output_path) if output_path != fsd_raw_path else None,
            )
            if arms_processed > 0:
                logger.info(f"FSD PON: {arms_processed} arms normalized")
            else:
                logger.warning("FSD PON: no arms processed")
        except Exception as e:
            logger.error(f"FSD PON processing failed: {e}")
            raise RuntimeError(f"FSD PON processing failed: {e}")
    else:
        logger.debug(f"No PON provided for FSD: {fsd_raw_path}")

    # Re-write through write_table() to honour output_format and compress.
    # Rust always writes raw TSV; we read it back and produce the requested
    # output format(s): .tsv, .tsv.gz, .parquet, or both.
    _write_fsd_output(output_path, output_format, compress)

    return output_path


def _write_fsd_output(tsv_path: Path, output_format: str, compress: bool) -> None:
    """Read Rust-written FSD TSV and re-write through write_table().

    This ensures --output-format and --compress flags are honoured
    for FSD output, producing .parquet and/or .tsv.gz as requested.

    Args:
        tsv_path: Path to the Rust-written TSV file
        output_format: One of "tsv", "parquet", or "both"
        compress: If True, gzip-compress TSV output
    """
    df = read_table(tsv_path)
    if df is None:
        logger.warning(f"FSD output not found for format conversion: {tsv_path}")
        return

    logger.debug(
        f"FSD format conversion: {len(df)} rows, "
        f"output_format={output_format}, compress={compress}"
    )

    # write_table() strips known extensions (.tsv, .tsv.gz, .parquet) from
    # the path and appends the correct one(s), so passing the raw .tsv path
    # works correctly here.
    write_table(df, tsv_path, output_format=output_format, compress=compress)

    # If output_format is "parquet" only, clean up the intermediate Rust TSV
    # (write_table already wrote the .parquet file)
    if output_format == "parquet" and tsv_path.exists():
        tsv_path.unlink()
        logger.debug(f"Removed intermediate TSV: {tsv_path.name}")
