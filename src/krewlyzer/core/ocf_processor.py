"""
OCF (Open Chromatin Footprinting) PON processor.

Adds PON-normalized z-scores to OCF output files using high-performance Rust backend.
"""

from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger("core.ocf_processor")


def process_ocf_with_pon(
    ocf_path: Path, pon_parquet_path: Path, output_path: Optional[Path] = None
) -> int:
    """
    Add PON z-scores to OCF output using Rust implementation.

    Uses high-performance Rust implementation for z-score computation.

    Args:
        ocf_path: Path to sample OCF TSV
        pon_parquet_path: Path to PON Parquet file
        output_path: Output path (default: overwrite input)

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
    return n_matched
