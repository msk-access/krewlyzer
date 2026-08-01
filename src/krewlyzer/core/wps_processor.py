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

from krewlyzer.pon.model import GENOME_WIDE_GROUP, zscore_or_nan

logger = logging.getLogger("core.wps_processor")


def post_process_wps(
    wps_parquet: Path,
    wps_background_parquet: Optional[Path] = None,
    pon=None,
    extract_periodicity: bool = True,
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
        "periodicity_z": None,
    }

    # WPS foreground z-scoring lives in `core/wps_pon.py` and is applied by
    # the caller, not here. What stood in this place called
    # `_core.wps.apply_pon_zscore` behind `getattr(pon, "_source_path", None)`
    # -- an attribute nothing ever set -- so it never ran, and the Rust
    # function it guarded computed a scalar z from v1.0 baseline fields that
    # no shipped PON has carried since the vector format landed.

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
                    result["periodicity_score"] = float(
                        df.loc[global_mask, "periodicity_score"].iloc[0]
                    )
                    if "nrl_bp" in df.columns:
                        result["nrl_bp"] = float(df.loc[global_mask, "nrl_bp"].iloc[0])
                elif len(df) > 0:
                    result["periodicity_score"] = float(df["periodicity_score"].iloc[0])
                    if "nrl_bp" in df.columns:
                        result["nrl_bp"] = float(df["nrl_bp"].iloc[0])

            # Score every group, not just the genome-wide one.
            #
            # The baseline holds 28 groups -- Global_All, one per chromosome,
            # and one per Alu family -- and the output has an nrl_bp for each.
            # Only Global_All was ever scored, so 27 of 28 baselines were built
            # and never used. Per-chromosome NRL drift is the point of having
            # them.
            if pon is not None and getattr(pon, "wps_background_baseline", None):
                baseline = pon.wps_background_baseline
                for column, stats in (
                    ("nrl_bp", baseline.get_nrl_stats),
                    ("periodicity_score", baseline.get_periodicity_stats),
                ):
                    if column not in df.columns:
                        continue
                    target = "nrl_z" if column == "nrl_bp" else "periodicity_z"
                    scores, absent = [], 0
                    for group_id, observed in zip(df["group_id"], df[column]):
                        pair = stats(str(group_id))
                        if pair is None:
                            absent += 1
                            scores.append(np.nan)
                            continue
                        scores.append(zscore_or_nan(float(observed), pair[0], pair[1]))
                    df[target] = scores
                    scored = int(np.isfinite(np.asarray(scores, dtype=float)).sum())
                    logger.info(f"{target}: {scored}/{len(df)} groups scored")
                    if absent:
                        logger.warning(
                            f"{target}: {absent}/{len(df)} groups are absent from "
                            "the PON's wps_background baseline. Rebuild the PON "
                            "if these groups should be there."
                        )

                # Summary values for the caller, from the genome-wide group.
                summary = df[df["group_id"] == GENOME_WIDE_GROUP]
                if not summary.empty:
                    for key, column in (
                        ("nrl_z", "nrl_z"),
                        ("periodicity_z", "periodicity_z"),
                    ):
                        if column in summary.columns:
                            value = summary[column].iloc[0]
                            result[key] = None if pd.isna(value) else float(value)

            df.to_parquet(wps_background_parquet, index=False)
            logger.info(f"Processed background WPS: {wps_background_parquet}")

        except Exception as e:
            logger.error(f"Background WPS processing failed: {e}")
            raise RuntimeError(f"Background WPS processing failed: {e}")

    return result
