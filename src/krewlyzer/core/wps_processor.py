"""
WPS (Windowed Protection Score) processor.

Provides WPS-specific processing including:
- PON z-score overlay (Rust-accelerated via apply_pon_zscore)
- Savitzky-Golay smoothing and FFT periodicity (Rust)

Performance:
- Rust PON z-score: 5-20x faster than Python pandas merge
- Python fallback available if Rust fails

Used by both standalone wps.py and run-all wrapper.py.
"""

from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger("core.wps_processor")


# =============================================================================
# Unified Post-Processing (shared by standalone wps and run-all)
# =============================================================================

def post_process_wps(
    wps_parquet: Path,
    wps_background_parquet: Optional[Path] = None,
    pon_baseline_parquet: Optional[Path] = None,
    pon_parquet_path: Optional[Path] = None,  # Direct path for Rust
    smooth: bool = True,
    extract_periodicity: bool = True
) -> dict:
    """
    Unified WPS post-processing pipeline.
    
    Called by both standalone `krewlyzer wps` and `run-all` for consistent output.
    Uses high-performance Rust implementation when available.
    
    Steps:
    1. Add *_smooth columns for backward compatibility (Rust already smoothed)
    2. PoN subtraction if baseline provided (adds *_delta, *_z columns)
    3. Extract periodicity score from Rust-computed nrl_bp column
    
    Args:
        wps_parquet: Path to foreground WPS Parquet (e.g., sample.WPS.parquet)
        wps_background_parquet: Path to background WPS Parquet (optional)
        pon_baseline_parquet: Path to PoN baseline Parquet (for Python fallback)
        pon_parquet_path: Direct path to PON parquet (for Rust implementation)
        smooth: Apply Savitzky-Golay smoothing (handled by Rust)
        extract_periodicity: Extract FFT periodicity from background
        
    Returns:
        dict with processing summary and metrics
    """
    result = {
        "wps_parquet": str(wps_parquet),
        "smoothed": False,
        "pon_subtracted": False,
        "periodicity_extracted": False,
        "periodicity_score": None
    }
    
    # Try Rust PON z-score first (5-20x faster)
    if pon_parquet_path and pon_parquet_path.exists() and wps_parquet.exists():
        try:
            from krewlyzer import _core
            regions_processed = _core.wps.apply_pon_zscore(
                str(wps_parquet),
                str(pon_parquet_path),
                None  # Overwrite in place
            )
            if regions_processed > 0:
                result["pon_subtracted"] = True
                logger.info(f"WPS PON (Rust): {regions_processed} regions normalized")
        except Exception as e:
            logger.debug(f"Rust WPS PON failed, will use Python fallback: {e}")
    
    # 1. Add _smooth columns for compatibility (Rust already smoothed the data)
    # Note: Rust applies Savitzky-Golay smoothing (window=11, order=3) before writing Parquet
    if smooth and wps_parquet.exists():
        try:
            logger.info(f"Smoothing foreground WPS: {wps_parquet}")
            df = pd.read_parquet(wps_parquet)
            
            for col in ["wps_nuc", "wps_tf", "prot_frac_nuc", "prot_frac_tf"]:
                if col in df.columns and f"{col}_smooth" not in df.columns:
                    # Rust already smoothed - just copy for compatibility
                    df[f"{col}_smooth"] = df[col]
            
            df.to_parquet(wps_parquet, index=False)
            result["smoothed"] = True
            logger.info(f"Added smooth columns to {wps_parquet}")
        except Exception as e:
            logger.warning(f"Failed to smooth foreground WPS: {e}")
    
    # 2. Process background WPS - periodicity is computed in Rust
    # Note: Rust handles background smoothing (window=7, order=3) and FFT periodicity
    if wps_background_parquet and wps_background_parquet.exists():
        try:
            logger.info(f"Processing background WPS: {wps_background_parquet}")
            df = pd.read_parquet(wps_background_parquet)
            
            # Periodicity and NRL are computed in Rust, stored in periodicity_score/nrl_bp columns
            if extract_periodicity and "periodicity_score" in df.columns:
                result["periodicity_extracted"] = True
                
                # Get Global_All score for the summary
                global_mask = df.get("group_id", pd.Series()) == "Global_All"
                if global_mask.any():
                    result["periodicity_score"] = float(df.loc[global_mask, "periodicity_score"].iloc[0])
                elif len(df) > 0:
                    # Use first row if no Global_All
                    result["periodicity_score"] = float(df["periodicity_score"].iloc[0])
            
            df.to_parquet(wps_background_parquet, index=False)
            logger.info(f"Processed background WPS: {wps_background_parquet}")
        except Exception as e:
            logger.warning(f"Failed to process background WPS: {e}")
    
    # 3. PoN subtraction
    if pon_baseline_parquet and Path(pon_baseline_parquet).exists() and wps_parquet.exists():
        try:
            logger.info(f"Subtracting PoN baseline: {pon_baseline_parquet}")
            subtract_pon_baseline(wps_parquet, pon_baseline_parquet, wps_parquet)
            result["pon_subtracted"] = True
        except Exception as e:
            logger.warning(f"Failed to subtract PoN baseline: {e}")
    
    return result


# =============================================================================
# PoN Subtraction for WPS Parquet
# =============================================================================

def subtract_pon_baseline(
    sample_parquet: Path,
    pon_baseline_parquet: Path,
    output_path: Optional[Path] = None,
    columns: list = None
) -> pd.DataFrame:
    """
    Subtract PoN (Panel of Normals) baseline from sample WPS profiles.
    
    For each region, computes:
    - delta_* = sample - pon_mean (raw difference)
    - z_* = (sample - pon_mean) / pon_std (z-score)
    
    Args:
        sample_parquet: Path to sample WPS Parquet
        pon_baseline_parquet: Path to PoN baseline Parquet (aggregated healthy samples)
        output_path: Output path (defaults to .pon_subtracted.parquet)
        columns: Vector columns to process (default: wps_nuc, wps_tf)
        
    Returns:
        DataFrame with PoN-subtracted profiles
    """
    if columns is None:
        columns = ["wps_nuc", "wps_tf"]
    
    logger.info(f"Loading sample: {sample_parquet}")
    sample_df = pd.read_parquet(sample_parquet)
    
    logger.info(f"Loading PoN baseline: {pon_baseline_parquet}")
    pon_df = pd.read_parquet(pon_baseline_parquet)
    
    # Build PoN lookup by region_id
    pon_lookup = {}
    region_id_col = "region_id" if "region_id" in pon_df.columns else "group_id"
    
    for _, row in pon_df.iterrows():
        region_id = row[region_id_col]
        pon_lookup[region_id] = {
            col: {
                "mean": np.array(row.get(f"{col}_mean", row.get(col, [])), dtype=np.float32),
                "std": np.array(row.get(f"{col}_std", [1.0] * len(row.get(col, []))), dtype=np.float32)
            }
            for col in columns if col in pon_df.columns or f"{col}_mean" in pon_df.columns
        }
    
    logger.info(f"PoN baseline contains {len(pon_lookup)} regions")
    
    # Process each sample row
    for col in columns:
        if col not in sample_df.columns:
            continue
        
        delta_col = []
        z_col = []
        
        for idx, row in sample_df.iterrows():
            region_id = row.get("region_id", row.get("group_id", str(idx)))
            sample_vec = np.array(row[col], dtype=np.float32)
            
            if region_id in pon_lookup and col in pon_lookup[region_id]:
                pon_mean = pon_lookup[region_id][col]["mean"]
                pon_std = pon_lookup[region_id][col]["std"]
                
                # Ensure same length
                if len(pon_mean) == len(sample_vec):
                    delta = sample_vec - pon_mean
                    # Avoid division by zero
                    pon_std_safe = np.where(pon_std > 1e-6, pon_std, 1.0)
                    z = delta / pon_std_safe
                else:
                    delta = sample_vec
                    z = np.zeros_like(sample_vec)
            else:
                delta = sample_vec
                z = np.zeros_like(sample_vec)
            
            delta_col.append(delta.tolist())
            z_col.append(z.tolist())
        
        sample_df[f"{col}_delta"] = delta_col
        sample_df[f"{col}_z"] = z_col
        logger.info(f"Added {col}_delta and {col}_z columns")
    
    if output_path is None:
        output_path = sample_parquet.with_suffix('.pon_subtracted.parquet')
    
    logger.info(f"Writing PoN-subtracted Parquet: {output_path}")
    sample_df.to_parquet(output_path, index=False)
    
    return sample_df


# =============================================================================
# Note: Smoothing and FFT periodicity are now handled by Rust for performance:
#   - Foreground WPS: Savitzky-Golay smoothing (window=11, order=3) via sci-rs
#   - Background WPS: Savitzky-Golay smoothing (window=7, order=3) via sci-rs
#   - Periodicity/NRL: FFT-based extraction via realfft, stored in nrl_bp column
# =============================================================================
