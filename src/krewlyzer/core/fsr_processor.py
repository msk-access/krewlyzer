"""
FSR (Fragment Size Ratio) processor.

Shared processing logic for aggregating fragment counts and computing ratios.
Used by both standalone fsr.py and run-all wrapper.py.
"""

from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger("core.fsr_processor")


def process_fsr(
    counts_df: pd.DataFrame,
    output_path: Path,
    windows: int = 100000,
    continue_n: int = 50,
    pon = None
) -> Path:
    """
    Process fragment counts into FSR ratios.
    
    Takes raw bin-level counts (from run_unified_pipeline) and:
    1. Aggregates into windows
    2. Computes fragment ratios for each window
    3. Optionally applies PON normalization before ratio calculation
    4. Writes output TSV
    
    Args:
        counts_df: DataFrame with columns [chrom, start, end, ultra_short, short, intermediate, long, total, mean_gc]
        output_path: Path to write FSR.tsv output
        windows: Window size in bp (default 100000)
        continue_n: Number of bins to aggregate (default 50)
        pon: Optional PON model for GC normalization
        
    Returns:
        Path to the output file
    """
    logger.info(f"Processing FSR: {len(counts_df)} bins → {output_path}")
    
    # Validate required columns
    required_cols = ['chrom', 'ultra_short', 'short', 'intermediate', 'long', 'total']
    missing = [c for c in required_cols if c not in counts_df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    # Apply PON GC correction if available
    if pon is not None and 'mean_gc' in counts_df.columns:
        from .pon_integration import compute_gc_bias_correction
        counts_df = compute_gc_bias_correction(
            counts_df.copy(), pon,
            gc_column='mean_gc',
            size_columns=('short', 'intermediate', 'long')
        )
    
    # Aggregate bins into windows and compute ratios
    results = []
    for chrom, group in counts_df.groupby('chrom', sort=False):
        n_bins = len(group)
        n_windows = n_bins // continue_n
        
        if n_windows == 0:
            continue
        
        trunc_len = n_windows * continue_n
        
        # Reshape and sum
        ultra_mat = group['ultra_short'].values[:trunc_len].reshape(n_windows, continue_n)
        short_mat = group['short'].values[:trunc_len].reshape(n_windows, continue_n)
        inter_mat = group['intermediate'].values[:trunc_len].reshape(n_windows, continue_n)
        long_mat = group['long'].values[:trunc_len].reshape(n_windows, continue_n)
        total_mat = group['total'].values[:trunc_len].reshape(n_windows, continue_n)
        
        ultra_sums = ultra_mat.sum(axis=1)
        short_sums = short_mat.sum(axis=1)
        inter_sums = inter_mat.sum(axis=1)
        long_sums = long_mat.sum(axis=1)
        total_sums = total_mat.sum(axis=1)
        
        window_starts = np.arange(n_windows) * continue_n * windows
        window_ends = (np.arange(n_windows) + 1) * continue_n * windows - 1
        
        for i in range(n_windows):
            region = f"{chrom}:{window_starts[i]}-{window_ends[i]}"
            
            u = float(ultra_sums[i])
            s = float(short_sums[i])
            m = float(inter_sums[i])
            l = float(long_sums[i])
            t = float(total_sums[i])
            
            # Calculate ratios (avoid division by zero)
            if t > 0:
                s_r = s / t
                m_r = m / t
                l_r = l / t
                u_r = u / t
            else:
                s_r = m_r = l_r = u_r = 0.0
            
            # Short/Long ratio (primary cancer biomarker)
            if l > 0:
                sl_r = s / l
            else:
                sl_r = s if s > 0 else 0.0
            
            results.append({
                'region': region,
                'ultra_short_count': int(u),
                'short_count': int(s),
                'inter_count': int(m),
                'long_count': int(l),
                'total_count': int(t),
                'short_ratio': s_r,
                'inter_ratio': m_r,
                'long_ratio': l_r,
                'short_long_ratio': sl_r,
                'ultra_short_ratio': u_r
            })
    
    if not results:
        logger.warning("No valid windows found for FSR")
        # Write empty file with header
        pd.DataFrame(columns=[
            'region', 'ultra_short_count', 'short_count', 'inter_count',
            'long_count', 'total_count', 'short_ratio', 'inter_ratio',
            'long_ratio', 'short_long_ratio', 'ultra_short_ratio'
        ]).to_csv(output_path, sep='\t', index=False)
        return output_path
    
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, sep='\t', index=False, float_format='%.6f')
    
    logger.info(f"FSR complete: {len(results_df)} windows → {output_path}")
    return output_path
