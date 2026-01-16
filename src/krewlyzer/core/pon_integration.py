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


def normalize_to_pon(
    df: pd.DataFrame, 
    pon, 
    columns: list
) -> pd.DataFrame:
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
        mean = pon.get_mean(col) if hasattr(pon, 'get_mean') else None
        if mean and mean > 0:
            df[f'{col}_norm'] = df[col] / mean
            logger.debug(f"Normalized {col} by PoN mean {mean:.2f}")
        else:
            df[f'{col}_norm'] = df[col]
            logger.debug(f"No PoN mean for {col}, using raw values")
    
    return df


def compute_gc_bias_correction(
    df: pd.DataFrame,
    pon,
    gc_column: str = "mean_gc",
    size_columns: Tuple[str, ...] = ("short", "intermediate", "long")
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
    value: float,
    arm: str,
    size_bin: Optional[str],
    pon
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


def compute_wps_zscore(
    value: float,
    gene: str,
    pon
) -> Optional[float]:
    """
    Compute WPS z-score for a single gene entry.
    
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
