"""
PON (Panel of Normals) integration module.

Provides shared functionality for loading PON models and computing z-scores
across all tools (FSC, FSR, FSD, WPS).
"""

from pathlib import Path
from typing import Optional, Tuple
import pandas as pd
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


def normalize_to_pon(df: pd.DataFrame, pon, columns: list) -> pd.DataFrame:
    """
    Divide count columns by PoN channel means.

    This normalizes raw counts to expected values from the panel of normals,
    enabling cross-sample comparison.

    Args:
        df: DataFrame with count columns
        pon: Loaded PonModel instance
        columns: List of column names to normalize

    Returns:
        DataFrame with new *_norm columns added
    """
    if pon is None:
        logger.debug("No PON provided, skipping normalization")
        return df

    for col in columns:
        mean = pon.get_mean(col) if hasattr(pon, "get_mean") else None
        if mean and mean > 0:
            df[f"{col}_norm"] = df[col] / mean
            logger.debug(f"Normalized {col} by PoN mean {mean:.2f}")
        else:
            df[f"{col}_norm"] = df[col]
            logger.debug(f"No PoN mean for {col}, using raw values")

    return df


def compute_gc_bias_correction(
    df: pd.DataFrame,
    pon,
    gc_column: str = "mean_gc",
    size_columns: Tuple[str, ...] = ("short", "intermediate", "long"),
) -> pd.DataFrame:
    """
    Apply PON-based GC bias correction to count data.

    This normalizes counts by the expected values from the PON's GC bias curve.

    Args:
        df: DataFrame with count columns
        pon: Loaded PonModel with gc_bias attribute
        gc_column: Name of column containing GC content (0-1)
        size_columns: Column names to correct

    Returns:
        DataFrame with corrected counts (modifies in place and returns)
    """
    if pon is None or pon.gc_bias is None:
        logger.debug("No PON GC bias available, skipping correction")
        return df

    logger.info("Applying PON-based GC normalization...")

    gc_arr = df[gc_column].values

    for col in size_columns:
        if col not in df.columns:
            continue

        values = df[col].values.astype(float)

        for i, gc in enumerate(gc_arr):
            expected = pon.gc_bias.get_expected(gc, col)
            if expected > 0:
                values[i] /= expected

        df[col] = values

    # Recalculate total if present
    if "total" in df.columns and all(c in df.columns for c in size_columns):
        df["total"] = sum(df[c] for c in size_columns)

    return df


def compute_fsd_zscore(
    value: float, arm: str, size_bin: Optional[str], pon
) -> Optional[float]:
    """
    Compute FSD z-score for a single (arm, size_bin) entry.

    Args:
        value: The observed FSD value
        arm: Chromosome arm (e.g., "1p", "12q")
        size_bin: Size bin identifier (e.g., "100-150")
        pon: Loaded PonModel with fsd_baseline

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.fsd_baseline is None:
        return None

    if arm not in pon.fsd_baseline.arms:
        return None

    stats = pon.fsd_baseline.get_stats(arm, size_bin)
    if stats is None:
        return None

    mean, std = stats
    if std > 0:
        return (value - mean) / std
    else:
        return 0.0


def compute_wps_zscore(value: float, gene: str, pon) -> Optional[float]:
    """
    Compute WPS z-score for a single gene entry (v1.0 scalar format).

    Args:
        value: The observed WPS value
        gene: Gene name
        pon: Loaded PonModel with wps_baseline

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.wps_baseline is None:
        return None

    stats = pon.wps_baseline.get_stats(gene)
    if stats is None:
        return None

    mean, std = stats
    if std > 0:
        return (value - mean) / std
    else:
        return 0.0


def compute_wps_z_vector(sample_vector, region_id: str, pon, column: str = "wps_nuc"):
    """
    Compute 200-element position-wise z-score vector (v2.0 format).

    Args:
        sample_vector: 200-element sample WPS vector (numpy array or list)
        region_id: Region identifier
        pon: Loaded PonModel with wps_baseline (v2.0)
        column: Vector column prefix ('wps_nuc' or 'wps_tf')

    Returns:
        200-element numpy array of z-scores or None if not computable
    """
    import numpy as np

    if pon is None or pon.wps_baseline is None:
        return None

    if pon.wps_baseline.schema_version != "2.0":
        logger.warning("compute_wps_z_vector requires v2.0 PON (vector format)")
        return None

    return pon.wps_baseline.compute_z_vector(
        region_id, np.asarray(sample_vector), column
    )


def compute_wps_panel_zscore(value: float, gene: str, pon) -> Optional[float]:
    """
    Compute WPS z-score for panel-specific anchors (v1.0 scalar format).

    Uses wps_baseline_panel instead of wps_baseline (genome-wide).
    For panel mode samples processed with panel-specific anchors.

    Args:
        value: The observed WPS value
        gene: Gene name
        pon: Loaded PonModel with wps_baseline_panel

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.wps_baseline_panel is None:
        return None

    stats = pon.wps_baseline_panel.get_stats(gene)
    if stats is None:
        return None

    mean, std = stats
    if std > 0:
        return (value - mean) / std
    return 0.0


def compute_wps_panel_z_vector(
    sample_vector, region_id: str, pon, column: str = "wps_nuc"
):
    """
    Compute 200-element position-wise z-score vector using panel baseline (v2.0 format).

    Uses wps_baseline_panel for panel-specific anchors.

    Args:
        sample_vector: 200-element sample WPS vector (numpy array or list)
        region_id: Region identifier
        pon: Loaded PonModel with wps_baseline_panel (v2.0)
        column: Vector column prefix ('wps_nuc' or 'wps_tf')

    Returns:
        200-element numpy array of z-scores or None if not computable
    """
    import numpy as np

    if pon is None or pon.wps_baseline_panel is None:
        return None

    if pon.wps_baseline_panel.schema_version != "2.0":
        logger.warning("compute_wps_panel_z_vector requires v2.0 PON (vector format)")
        return None

    return pon.wps_baseline_panel.compute_z_vector(
        region_id, np.asarray(sample_vector), column
    )


def compute_wps_shape_score(
    sample_vector, region_id: str, pon, column: str = "wps_nuc"
) -> Optional[float]:
    """
    Compute shape correlation score for cancer detection (v2.0 format).

    Returns Pearson correlation between sample WPS vector and PON mean shape.
    Healthy samples show correlation ~1.0, cancer samples show lower values.

    Args:
        sample_vector: 200-element sample WPS vector
        region_id: Region identifier
        pon: Loaded PonModel with wps_baseline (v2.0)
        column: Vector column prefix ('wps_nuc' or 'wps_tf')

    Returns:
        Correlation coefficient [-1, 1] or None if not computable

    Clinical interpretation:
        - 0.9-1.0: Healthy nucleosome positioning
        - 0.5-0.9: Mild chromatin disorganization
        - <0.5: Significant disruption (cancer signal)
    """
    import numpy as np

    if pon is None or pon.wps_baseline is None:
        return None

    if pon.wps_baseline.schema_version != "2.0":
        logger.warning("compute_wps_shape_score requires v2.0 PON (vector format)")
        return None

    return pon.wps_baseline.compute_shape_score(
        region_id, np.asarray(sample_vector), column
    )


def compute_tfbs_zscore(label: str, entropy: float, pon) -> Optional[float]:
    """
    Compute TFBS entropy z-score for a single TF label.

    Args:
        label: TF label (e.g., "CTCF", "FOXA1")
        entropy: Observed entropy value
        pon: Loaded PonModel with tfbs_baseline

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.tfbs_baseline is None:
        return None

    return pon.tfbs_baseline.get_zscore(label, entropy)


def compute_atac_zscore(label: str, entropy: float, pon) -> Optional[float]:
    """
    Compute ATAC entropy z-score for a single cancer type label.

    Args:
        label: Cancer type label (e.g., "BRCA", "LUAD")
        entropy: Observed entropy value
        pon: Loaded PonModel with atac_baseline

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.atac_baseline is None:
        return None

    return pon.atac_baseline.get_zscore(label, entropy)


def compute_tfbs_ontarget_zscore(label: str, entropy: float, pon) -> Optional[float]:
    """
    Compute TFBS entropy z-score using panel-specific (ontarget) baseline.

    Uses tfbs_baseline_ontarget which is built from panel-intersection TFBS regions.

    Args:
        label: TF label (e.g., "CTCF", "FOXA1")
        entropy: Observed entropy value
        pon: Loaded PonModel with tfbs_baseline_ontarget

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.tfbs_baseline_ontarget is None:
        return None

    return pon.tfbs_baseline_ontarget.get_zscore(label, entropy)


def compute_atac_ontarget_zscore(label: str, entropy: float, pon) -> Optional[float]:
    """
    Compute ATAC entropy z-score using panel-specific (ontarget) baseline.

    Uses atac_baseline_ontarget which is built from panel-intersection ATAC regions.

    Args:
        label: Cancer type label (e.g., "BRCA", "LUAD")
        entropy: Observed entropy value
        pon: Loaded PonModel with atac_baseline_ontarget

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.atac_baseline_ontarget is None:
        return None

    return pon.atac_baseline_ontarget.get_zscore(label, entropy)


def compute_ocf_zscore(region_id: str, ocf_value: float, pon) -> Optional[float]:
    """
    Compute OCF z-score for a single region.

    Centralized function for OCF z-score computation.
    Can be used by ocf_processor.py or directly.

    Args:
        region_id: Region identifier
        ocf_value: Observed OCF value
        pon: Loaded PonModel with ocf_baseline

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.ocf_baseline is None:
        return None

    return pon.ocf_baseline.compute_zscore(region_id, ocf_value)


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
