"""
FSD (Fragment Size Distribution) processor.

Provides FSD-specific processing including PoN log-ratio normalization.
Used by both standalone fsd.py and run-all wrapper.py.

Normalization order:
1. GC-weighting (Rust) - raw counts are GC-corrected
2. PoN log-ratio (Python) - log2(sample / PoN_expected)
"""

from pathlib import Path
from typing import Optional, List
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger("core.fsd_processor")

# 67 bins: 65-69, 70-74, ..., 395-399
BIN_COLUMNS = [f"{s}-{s+4}" for s in range(65, 400, 5)]


def process_fsd(
    fsd_raw_path: Path,
    output_path: Optional[Path] = None,
    pon = None
) -> Path:
    """
    Process raw FSD counts into ML-ready features.
    
    Converts raw GC-weighted counts from Rust into:
    - Log-ratios vs PoN (when PoN provided): log2(sample / PoN_expected)
    - PoN stability scores: 1 / (variance + k)
    
    Args:
        fsd_raw_path: Path to raw FSD.tsv from Rust 
        output_path: Output path (default: overwrite input)
        pon: Optional PonModel with fsd_baseline
        
    Returns:
        Path to processed output file
    """
    if not fsd_raw_path.exists():
        logger.warning(f"FSD output not found: {fsd_raw_path}")
        return fsd_raw_path
    
    output_path = output_path or fsd_raw_path
    
    logger.info(f"Processing FSD: {fsd_raw_path}")
    
    # Read raw counts from Rust
    df = pd.read_csv(fsd_raw_path, sep='\t')
    
    # Validate columns
    if 'region' not in df.columns:
        logger.warning("FSD output missing 'region' column")
        return fsd_raw_path
    
    # Find bin columns present in data
    bin_cols = [c for c in BIN_COLUMNS if c in df.columns]
    if not bin_cols:
        logger.warning("No bin columns found in FSD output")
        return fsd_raw_path
    
    total_frags = df['total'].sum() if 'total' in df.columns else df[bin_cols].sum().sum()
    logger.debug(f"  {len(df)} arms, {total_frags:,.0f} total fragments")
    
    # Compute log-ratios using PoN FsdBaseline
    if pon is not None and hasattr(pon, 'fsd_baseline') and pon.fsd_baseline is not None:
        logger.info("Computing PoN log-ratios...")
        fsd_baseline = pon.fsd_baseline
        
        # Process each arm (row)
        for idx, row in df.iterrows():
            arm = row['region']
            
            # Compute log-ratio for each bin in this arm
            for bin_col in bin_cols:
                # Parse bin start (e.g., "65-69" -> 65)
                try:
                    size = int(bin_col.split('-')[0])
                except ValueError:
                    continue
                
                sample_val = row[bin_col]
                
                # Get PoN expected value for this arm/size
                pon_expected = fsd_baseline.get_expected(arm, size)
                
                if pon_expected > 0:
                    # log2(sample / PoN_expected) with pseudocount
                    df.loc[idx, f'{bin_col}_logR'] = np.log2((sample_val + 1) / (pon_expected + 1))
                else:
                    # Fallback to sample mean if PoN not available for this arm
                    sample_mean = df[bin_col].mean()
                    if sample_mean > 0:
                        df.loc[idx, f'{bin_col}_logR'] = np.log2((sample_val + 1) / (sample_mean + 1))
                    else:
                        df.loc[idx, f'{bin_col}_logR'] = 0.0
            
            # Compute per-arm stability score (inverse variance)
            # Average std across all bins for this arm
            stds = []
            for bin_col in bin_cols[:10]:  # Sample first 10 bins for speed
                try:
                    size = int(bin_col.split('-')[0])
                    std = fsd_baseline.get_std(arm, size)
                    if std > 0:
                        stds.append(std)
                except:
                    pass
            
            if stds:
                avg_var = np.mean([s**2 for s in stds])
                df.loc[idx, 'pon_stability'] = 1.0 / (avg_var + 0.01)
            else:
                df.loc[idx, 'pon_stability'] = 1.0
        
        # Log statistics
        logR_cols = [c for c in df.columns if c.endswith('_logR')]
        if logR_cols:
            logR_min = df[logR_cols].min().min()
            logR_max = df[logR_cols].max().max()
            logger.debug(f"  Log-ratio range: [{logR_min:.2f}, {logR_max:.2f}]")
    else:
        logger.debug("No PoN FSD baseline available, outputting raw counts only")
    
    # Write output
    df.to_csv(output_path, sep='\t', index=False, float_format='%.6f')
    logger.info(f"FSD processed: {output_path}")
    
    return output_path


def apply_fsd_pon(
    fsd_output: Path,
    pon,
    output_path: Optional[Path] = None
) -> Path:
    """
    Apply PON z-scores to FSD output file (legacy compatibility).
    
    Deprecated: Use process_fsd() instead for log-ratio normalization.
    """
    return process_fsd(fsd_output, output_path, pon)
