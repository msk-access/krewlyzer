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

from .output_utils import read_table, resolve_table_path, write_table

logger = logging.getLogger("core.fsd_processor")


def process_fsd(
    fsd_raw_path: Path,
    output_path: Optional[Path] = None,
    pon_parquet_path: Optional[Path] = None,
    output_format: str = "tsv",
    compress: bool = False,
    baseline_table: str = "fsd_baseline",
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
        baseline_table: PON baseline table name ("fsd_baseline" or "fsd_baseline_ontarget")

    Returns:
        Path to processed output file (base path; actual files may have
        .tsv, .tsv.gz, or .parquet extensions depending on output_format)

    Raises:
        RuntimeError: If FSD processing fails
    """
    # The Rust writer became format-aware, so with --output-format parquet it
    # writes {sample}.FSD.parquet and no .tsv at all -- while callers still name
    # the .tsv. `exists()` on that name was then False and the whole of this
    # function, PON normalisation included, was skipped in silence. Measured:
    # FSD came out as raw counts under `parquet`, log-ratios under `tsv`, with
    # no warning either way. 0.9.0 makes parquet the Nextflow default, which
    # would have turned that from affecting nobody into affecting every run.
    resolved = resolve_table_path(fsd_raw_path)
    if resolved is None:
        raise FileNotFoundError(
            f"FSD output not found: {fsd_raw_path} (nor .tsv.gz / .parquet)"
        )

    output_path = output_path or fsd_raw_path

    # Apply PON log-ratio normalization via Rust
    # Rust writes the result as TSV to output_path (or in-place if same as input)
    if pon_parquet_path and pon_parquet_path.exists():
        try:
            from krewlyzer import _core

            # The normaliser parses TSV line by line, so anything else has to
            # be materialised first. Kept here rather than teaching Rust every
            # container format: the normalisation is what matters, not the box
            # it arrived in.
            normalise_input = resolved
            staged: Optional[Path] = None
            if resolved.suffix != ".tsv":
                staged = output_path.with_suffix(".pon_input.tsv")
                frame = read_table(resolved)
                if frame is None:
                    raise RuntimeError(f"could not read FSD table at {resolved}")
                frame.to_csv(staged, sep="\t", index=False)
                normalise_input = staged
                logger.debug(f"staged {resolved.name} as TSV for PON normalisation")

            arms_processed = _core.fsd.apply_pon_logratio(
                str(normalise_input),
                str(pon_parquet_path),
                str(output_path) if output_path != normalise_input else None,
                baseline_table=baseline_table,
            )
            if arms_processed > 0:
                logger.info(f"FSD PON: {arms_processed} arms normalized")
            else:
                # Not a debug line: zero arms means the output is raw counts
                # wearing the same column names as log-ratios, and nothing
                # downstream can tell the difference.
                logger.warning(
                    "FSD PON: no arms matched the baseline -- output is RAW "
                    "COUNTS, not log-ratios. Check that the PON arm names match "
                    "the arms BED used for this run."
                )
        except Exception as e:
            logger.error(f"FSD PON processing failed: {e}")
            raise RuntimeError(f"FSD PON processing failed: {e}")
        finally:
            if staged is not None and staged.exists():
                staged.unlink()
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
    # Read *this* file, not whatever read_table would resolve to. read_table
    # is parquet-first: given `x.FSD.tsv` it prefers `x.FSD.parquet`, which
    # here is the raw table the single-pass Rust writer emitted before
    # normalisation. Resolving would silently discard the log-ratios we just
    # computed and write the raw counts straight back -- which is exactly what
    # it did, with "41 arms normalized" logged immediately above.
    import pandas as pd

    if not tsv_path.exists():
        logger.warning(f"FSD output not found for format conversion: {tsv_path}")
        return
    df = pd.read_csv(tsv_path, sep="\t", comment="#")

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
