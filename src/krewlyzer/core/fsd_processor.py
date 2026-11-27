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

from .output_utils import (
    cleanup_intermediate_tsv,
    read_exact_table,
    write_table,
    written_table,
)

logger = logging.getLogger("core.fsd_processor")


def process_fsd(
    fsd_raw_path: Path,
    output_path: Optional[Path] = None,
    pon_parquet_path: Optional[Path] = None,
    *,
    output_format: str,
    compress: bool,
    baseline_table: str = "fsd_baseline",
) -> Path:
    """
    Process raw FSD counts into ML-ready features.

    Converts raw GC-weighted counts from Rust into:
    - Log-ratios vs PoN (when PoN provided): log2(sample / PoN_expected)
    - PoN stability scores: 1 / (variance + k)

    The Rust writer has already written the run's formats. With a PON the normalised result is
    written over them through write_table(), honouring --output-format and --compress; without
    one there is nothing to change and they ship as Rust wrote them.

    Args:
        fsd_raw_path: Path to raw FSD.tsv from Rust
        output_path: Output path (default: overwrite input).
            write_table() strips known extensions and appends the correct one(s).
        pon_parquet_path: Path to PON parquet for normalization
        output_format: One of "tsv", "parquet", or "both"; this run's, which decides the file read
        compress: If True, gzip-compress TSV output (.tsv.gz)
        baseline_table: PON baseline table name ("fsd_baseline" or "fsd_baseline_ontarget")

    Returns:
        Path to processed output file (base path; actual files may have
        .tsv, .tsv.gz, or .parquet extensions depending on output_format)

    Raises:
        FileNotFoundError: If this run's Rust pass wrote no FSD table
        RuntimeError: If FSD processing fails
    """
    # The table this run's Rust pass wrote, exactly. The writer is format-aware: under `parquet`
    # it writes {sample}.FSD.parquet and no .tsv, while callers name the .tsv. A bare `exists()`
    # on that name skipped this whole function, PON normalisation included, under `parquet`;
    # resolving across formats instead (`.tsv` first) took the uncompressed log-ratios an earlier
    # `--compress` run had left and normalised them a second time.
    source = written_table(fsd_raw_path, output_format=output_format, compress=compress)
    if not source.exists():
        raise FileNotFoundError(
            f"FSD output not found: {source} (this run's Rust pass wrote none)"
        )

    output_path = output_path or fsd_raw_path

    if not (pon_parquet_path and pon_parquet_path.exists()):
        # Nothing to normalise, and Rust has written exactly this run's formats. Rewriting them
        # dates from when it wrote only TSV: it read a plain FSD.tsv that `parquet` and
        # compressed runs never have, warning "not found" on every one, and under --compress
        # read one an earlier run with a PON had left, shipping that run's log-ratios.
        logger.debug(f"No PON provided for FSD: {fsd_raw_path}")
        return output_path

    # Apply PON log-ratio normalization via Rust
    # Rust writes the result as TSV to output_path (or in-place if same as input)
    staged: Optional[Path] = None
    try:
        from krewlyzer import _core

        # The normaliser parses TSV line by line, so anything else has to
        # be materialised first. Kept here rather than teaching Rust every
        # container format: the normalisation is what matters, not the box
        # it arrived in.
        normalise_input = source
        if source.suffix != ".tsv":
            staged = output_path.with_suffix(".pon_input.tsv")
            frame = read_exact_table(source)
            if frame is None:
                raise RuntimeError(f"could not read FSD table at {source}")
            frame.to_csv(staged, sep="\t", index=False)
            normalise_input = staged
            logger.debug(f"staged {source.name} as TSV for PON normalisation")

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

    # The normaliser wrote plain TSV; write it over the run's formats (.tsv, .tsv.gz, .parquet).
    _write_fsd_output(output_path, output_format, compress)

    return output_path


def _write_fsd_output(tsv_path: Path, output_format: str, compress: bool) -> None:
    """Read the normalised FSD TSV and re-write it through write_table().

    This ensures --output-format and --compress flags are honoured
    for FSD output, producing .parquet and/or .tsv.gz as requested.

    Args:
        tsv_path: Path to the TSV the PON normaliser wrote
        output_format: One of "tsv", "parquet", or "both"
        compress: If True, gzip-compress TSV output
    """
    # Read *this* file, not whatever read_table would resolve to. read_table
    # is parquet-first: given `x.FSD.tsv` it prefers `x.FSD.parquet`, which
    # here is the raw table the single-pass Rust writer emitted before
    # normalisation. Resolving would silently discard the log-ratios we just
    # computed and write the raw counts straight back -- which is exactly what
    # it did, with "41 arms normalized" logged immediately above.
    df = read_exact_table(tsv_path)
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

    # The normaliser's TSV goes unless it is one of the run's outputs. It went only under
    # `parquet` before, so a compressed run left an uncompressed FSD.tsv beside its .tsv.gz,
    # and the next run into the folder read that as its own table.
    cleanup_intermediate_tsv(tsv_path, output_format, compress)
