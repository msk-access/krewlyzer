"""
FSD (Fragment Size Distribution) processor.

Provides FSD-specific processing including PoN log-ratio normalization.
Used by both standalone fsd.py and run-all wrapper.py.

Normalization order:
1. GC-weighting (Rust) - raw counts are GC-corrected
2. PoN log-ratio (Rust) - log2(sample / PoN_expected) via apply_pon_logratio

All processing uses Rust implementation.
"""

from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger("core.fsd_processor")


def process_fsd(
    fsd_raw_path: Path,
    output_path: Optional[Path] = None,
    pon_parquet_path: Optional[Path] = None
) -> Path:
    """
    Process raw FSD counts into ML-ready features.
    
    Converts raw GC-weighted counts from Rust into:
    - Log-ratios vs PoN (when PoN provided): log2(sample / PoN_expected)
    - PoN stability scores: 1 / (variance + k)
    
    Uses Rust implementation for all processing.
    
    Args:
        fsd_raw_path: Path to raw FSD.tsv from Rust 
        output_path: Output path (default: overwrite input)
        pon_parquet_path: Path to PON parquet for normalization
        
    Returns:
        Path to processed output file
        
    Raises:
        RuntimeError: If FSD processing fails
    """
    if not fsd_raw_path.exists():
        raise FileNotFoundError(f"FSD output not found: {fsd_raw_path}")
    
    output_path = output_path or fsd_raw_path
    
    # Apply PON log-ratio normalization via Rust
    if pon_parquet_path and pon_parquet_path.exists():
        try:
            from krewlyzer import _core
            arms_processed = _core.fsd.apply_pon_logratio(
                str(fsd_raw_path),
                str(pon_parquet_path),
                str(output_path) if output_path != fsd_raw_path else None
            )
            if arms_processed > 0:
                logger.info(f"FSD PON: {arms_processed} arms normalized")
                return output_path
            else:
                logger.warning("FSD PON: no arms processed")
        except Exception as e:
            logger.error(f"FSD PON processing failed: {e}")
            raise RuntimeError(f"FSD PON processing failed: {e}")
    else:
        logger.debug(f"No PON provided for FSD: {fsd_raw_path}")
    
    return output_path
