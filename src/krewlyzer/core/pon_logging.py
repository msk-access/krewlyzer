"""
PON (Panel of Normals) logging utilities.

Provides structured logging for PON operations across all processors.
"""

import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger("pon")


def log_pon_loaded(pon, source_path: Optional[Path] = None) -> None:
    """
    Log PON model loading with component summary.
    
    Args:
        pon: Loaded PonModel instance
        source_path: Optional path to the PON file
    """
    if pon is None:
        logger.debug("No PON model loaded")
        return
    
    components = []
    if pon.gc_bias:
        components.append(f"GC-bias({len(pon.gc_bias.gc_bins)} bins)")
    if pon.fsd_baseline:
        components.append(f"FSD({len(pon.fsd_baseline.arms)} arms)")
    if pon.wps_baseline:
        components.append(f"WPS({len(pon.wps_baseline.regions)} regions)")
    if pon.ocf_baseline:
        components.append(f"OCF({len(pon.ocf_baseline.regions)} regions)")
    if pon.mds_baseline:
        components.append(f"MDS({len(pon.mds_baseline.kmer_expected)} kmers)")
    
    source = f" from {source_path.name}" if source_path else ""
    logger.info(f"PON model{source}: {pon.assay} (n={pon.n_samples})")
    if components:
        logger.info(f"  Components: {', '.join(components)}")


def log_normalization_mode(processor_name: str, has_pon: bool, mode: str = "auto") -> None:
    """
    Log the normalization mode for a processor.
    
    Args:
        processor_name: Name of the processor (FSC, FSR, FSD, WPS, OCF)
        has_pon: Whether PON model is available
        mode: Normalization mode description
    """
    if has_pon:
        logger.debug(f"{processor_name}: Using PON log-ratio normalization")
    else:
        logger.debug(f"{processor_name}: Outputting raw counts (no PON)")


def log_pon_zscore_summary(processor_name: str, n_features: int, 
                           n_normalized: int, n_missing: int = 0) -> None:
    """
    Log z-score computation summary for a processor.
    
    Args:
        processor_name: Name of the processor
        n_features: Total features processed
        n_normalized: Features with PON normalization applied
        n_missing: Features with missing PON baseline data
    """
    if n_missing > 0:
        logger.debug(
            f"{processor_name}: {n_normalized}/{n_features} normalized, "
            f"{n_missing} missing baseline"
        )
    else:
        logger.debug(f"{processor_name}: {n_normalized}/{n_features} normalized")


class PonNormalizationStats:
    """Accumulator for PON normalization statistics."""
    
    def __init__(self, processor_name: str):
        self.processor_name = processor_name
        self.n_features = 0
        self.n_normalized = 0
        self.n_missing = 0
    
    def add(self, normalized: bool, has_baseline: bool = True) -> None:
        """Record a feature's normalization status."""
        self.n_features += 1
        if normalized:
            self.n_normalized += 1
        if not has_baseline:
            self.n_missing += 1
    
    def log_summary(self) -> None:
        """Log the accumulated stats."""
        log_pon_zscore_summary(
            self.processor_name,
            self.n_features,
            self.n_normalized,
            self.n_missing
        )
