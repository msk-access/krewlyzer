"""
OCF (Open Chromatin Footprinting) PON processor.

Adds PON-normalized z-scores to OCF output files.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, TYPE_CHECKING
import logging

if TYPE_CHECKING:
    from ..pon.model import OcfBaseline

logger = logging.getLogger("core.ocf_processor")


def process_ocf_with_pon(
    ocf_path: Path,
    ocf_baseline: "OcfBaseline",
    output_path: Optional[Path] = None
) -> pd.DataFrame:
    """
    Add PON z-scores to OCF output.
    
    Args:
        ocf_path: Path to sample OCF TSV
        ocf_baseline: OcfBaseline from PON model
        output_path: Output path (default: overwrite input)
    
    Returns:
        DataFrame with added ocf_z column
    """
    if not ocf_path.exists():
        logger.warning(f"OCF file not found: {ocf_path}")
        return pd.DataFrame()
    
    df = pd.read_csv(ocf_path, sep="\t")
    
    if df.empty:
        logger.warning(f"OCF file is empty: {ocf_path}")
        return df
    
    z_scores = []
    n_matched = 0
    
    for _, row in df.iterrows():
        # Try multiple column names for region ID
        region_id = row.get("region_id", row.get("name", row.get("region", "")))
        
        # Try multiple column names for OCF value
        ocf_value = row.get("ocf", row.get("ocf_score", row.get("score", 0)))
        
        stats = ocf_baseline.get_stats(str(region_id))
        if stats:
            mean, std = stats
            if std > 0:
                z = (ocf_value - mean) / std
            else:
                z = 0.0
            n_matched += 1
        else:
            z = np.nan
        z_scores.append(z)
    
    df["ocf_z"] = z_scores
    
    out = output_path or ocf_path
    df.to_csv(out, sep="\t", index=False)
    
    logger.info(f"OCF PON z-scores: {n_matched}/{len(df)} regions matched ({out.name})")
    
    return df
