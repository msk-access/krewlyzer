"""
WPS (Windowed Protection Score) processor.

Provides WPS-specific post-processing including:
- PON z-score computation (via Rust)
- Background NRL/periodicity z-scores

All processing uses Rust implementation. No Python fallbacks.
"""

from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger("core.wps_processor")


def post_process_wps(
    wps_parquet: Path,
    wps_background_parquet: Optional[Path] = None,
    pon = None,
    extract_periodicity: bool = True
) -> dict:
    """
    Unified WPS post-processing pipeline.
    
    Called by both standalone `krewlyzer wps` and `run-all` for consistent output.
    Uses Rust implementation for PON z-score computation.
    
    Steps:
    1. Apply Rust PON z-score to foreground WPS
    2. Extract NRL/periodicity from background WPS
    3. Compute NRL/periodicity z-scores vs PON baseline
    
    Args:
        wps_parquet: Path to foreground WPS Parquet (e.g., sample.WPS.parquet)
        wps_background_parquet: Path to background WPS Parquet (optional)
        pon: Loaded PonModel for z-score computation
        extract_periodicity: Extract FFT periodicity from background
        
    Returns:
        dict with processing summary and metrics
    """
    result = {
        "wps_parquet": str(wps_parquet),
        "pon_subtracted": False,
        "periodicity_extracted": False,
        "periodicity_score": None,
        "nrl_bp": None,
        "nrl_z": None,
        "periodicity_z": None
    }
    
    # Apply Rust PON z-score to foreground WPS
    if pon is not None and wps_parquet.exists():
        try:
            from krewlyzer import _core
            # Get PON path from model if available
            pon_path = getattr(pon, '_source_path', None)
            if pon_path:
                regions_processed = _core.wps.apply_pon_zscore(
                    str(wps_parquet),
                    str(pon_path),
                    None  # Overwrite in place
                )
                if regions_processed > 0:
                    result["pon_subtracted"] = True
                    logger.info(f"WPS PON z-score: {regions_processed} regions normalized")
        except Exception as e:
            logger.error(f"WPS PON z-score failed: {e}")
            raise RuntimeError(f"WPS PON z-score computation failed: {e}")
    
    # Process background WPS - periodicity/NRL computed in Rust
    if wps_background_parquet and wps_background_parquet.exists():
        try:
            df = pd.read_parquet(wps_background_parquet)
            
            # Extract periodicity and NRL from Rust-computed columns
            if extract_periodicity and "periodicity_score" in df.columns:
                result["periodicity_extracted"] = True
                
                # Get Global_All score for the summary
                global_mask = df.get("group_id", pd.Series()) == "Global_All"
                if global_mask.any():
                    result["periodicity_score"] = float(df.loc[global_mask, "periodicity_score"].iloc[0])
                    if "nrl_bp" in df.columns:
                        result["nrl_bp"] = float(df.loc[global_mask, "nrl_bp"].iloc[0])
                elif len(df) > 0:
                    result["periodicity_score"] = float(df["periodicity_score"].iloc[0])
                    if "nrl_bp" in df.columns:
                        result["nrl_bp"] = float(df["nrl_bp"].iloc[0])
            
            # Compute NRL and periodicity z-scores if PON available
            if pon is not None and result["nrl_bp"] is not None:
                from .pon_integration import compute_nrl_zscore, compute_periodicity_zscore
                
                nrl_z = compute_nrl_zscore(result["nrl_bp"], pon)
                if nrl_z is not None:
                    result["nrl_z"] = nrl_z
                    if "nrl_z" not in df.columns:
                        df["nrl_z"] = np.nan
                    if global_mask.any():
                        df.loc[global_mask, "nrl_z"] = nrl_z
                    logger.info(f"NRL z-score: {nrl_z:.2f}")
                
                if result["periodicity_score"] is not None:
                    period_z = compute_periodicity_zscore(result["periodicity_score"], pon)
                    if period_z is not None:
                        result["periodicity_z"] = period_z
                        if "periodicity_z" not in df.columns:
                            df["periodicity_z"] = np.nan
                        if global_mask.any():
                            df.loc[global_mask, "periodicity_z"] = period_z
                        logger.info(f"Periodicity z-score: {period_z:.2f}")
            
            df.to_parquet(wps_background_parquet, index=False)
            logger.info(f"Processed background WPS: {wps_background_parquet}")
            
        except Exception as e:
            logger.error(f"Background WPS processing failed: {e}")
            raise RuntimeError(f"Background WPS processing failed: {e}")
    
    return result
