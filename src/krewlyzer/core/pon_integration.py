"""
PON (Panel of Normals) integration module.

Provides shared functionality for loading PON models and computing
background/periodicity z-scores used by the WPS processor.

Note: Tool-specific z-score functions (FSC, FSD, WPS, OCF, TFBS, ATAC) have
been migrated to Rust for performance and are called via _core.<module>.
See the following Rust entry points:
    - _core.fsd.apply_pon_logratio      → fsd_processor.py
    - _core.wps.apply_pon_zscore        → wps_processor.py
    - _core.region_entropy.apply_pon_zscore → region_entropy_processor.py
    - _core.ocf.apply_pon_zscore        → ocf_processor.py
"""

from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger("core.pon_integration")


def load_pon_model(pon_path: Path):
    """
    Load a PON model from disk.

    Args:
        pon_path: Path to the PON parquet file

    Returns:
        Loaded PonModel instance, or None if loading fails
    """
    from krewlyzer.pon.model import PonModel

    try:
        pon = PonModel.load(pon_path)
        logger.info(f"Loaded PON model: {pon.assay} (n={pon.n_samples})")
        return pon
    except Exception as e:
        logger.warning(f"Could not load PON model: {e}")
        return None


def compute_nrl_zscore(
    observed_nrl: float, pon, group_id: str = "all"
) -> Optional[float]:
    """
    Compute NRL (Nucleosome Repeat Length) z-score.

    Compares sample NRL (from Rust FFT) against healthy reference
    from WPS background baseline.

    NRL ~180-200bp in healthy plasma; deviations indicate
    altered nucleosome spacing (cancer signature).

    Args:
        observed_nrl: NRL in bp from FFT analysis (nrl_bp column)
        pon: Loaded PonModel with wps_background_baseline
        group_id: Group identifier (default: "all" for global)

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.wps_background_baseline is None:
        return None

    return pon.wps_background_baseline.compute_nrl_zscore(observed_nrl, group_id)


def compute_periodicity_zscore(
    observed_periodicity: float, pon, group_id: str = "all"
) -> Optional[float]:
    """
    Compute periodicity z-score.

    Periodicity measures the strength of nucleosome signal from FFT.
    Higher values indicate more regular nucleosome positioning.

    Args:
        observed_periodicity: Periodicity score from FFT
        pon: Loaded PonModel with wps_background_baseline
        group_id: Group identifier (default: "all" for global)

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.wps_background_baseline is None:
        return None

    stats = pon.wps_background_baseline.get_periodicity_stats(group_id)
    if stats is None:
        return None

    mean, std = stats
    if std > 0:
        return (observed_periodicity - mean) / std
    return 0.0
