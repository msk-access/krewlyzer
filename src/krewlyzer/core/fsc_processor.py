"""
FSC (Fragment Size Counts) processor.

Shared processing logic for aggregating fragment counts and computing z-scores.
Used by both standalone fsc.py and run-all wrapper.py.
"""

from pathlib import Path
from typing import Optional, Tuple
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger("core.fsc_processor")


def process_fsc(
    counts_df: pd.DataFrame,
    output_path: Path,
    windows: int = 100000,
    continue_n: int = 50,
    pon = None
) -> Path:
    """
    Process fragment counts into FSC z-scores.
    
    Takes raw bin-level counts (from run_unified_pipeline) and:
    1. Aggregates into windows
    2. Computes z-scores for each fragment size category
    3. Optionally applies PON normalization
    4. Writes output TSV
    
    Args:
        counts_df: DataFrame with columns [chrom, start, end, short, intermediate, long, total, mean_gc]
        output_path: Path to write FSC.tsv output
        windows: Window size in bp (default 100000)
        continue_n: Number of bins to aggregate (default 50)
        pon: Optional PON model for normalization
        
    Returns:
        Path to the output file
    """
    logger.info(f"Processing FSC: {len(counts_df)} bins → {output_path}")
    
    # Validate required columns
    required_cols = ['chrom', 'start', 'end', 'short', 'intermediate', 'long', 'total']
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
    
    # Aggregate bins into windows
    results = []
    for chrom, group in counts_df.groupby('chrom', sort=False):
        n_bins = len(group)
        n_windows = n_bins // continue_n
        
        if n_windows == 0:
            continue
        
        trunc_len = n_windows * continue_n
        
        # Reshape and sum
        shorts_mat = group['short'].values[:trunc_len].reshape(n_windows, continue_n)
        inter_mat = group['intermediate'].values[:trunc_len].reshape(n_windows, continue_n)
        longs_mat = group['long'].values[:trunc_len].reshape(n_windows, continue_n)
        totals_mat = group['total'].values[:trunc_len].reshape(n_windows, continue_n)
        
        sum_shorts = shorts_mat.sum(axis=1)
        sum_inter = inter_mat.sum(axis=1)
        sum_longs = longs_mat.sum(axis=1)
        sum_totals = totals_mat.sum(axis=1)
        
        window_starts = np.arange(n_windows) * continue_n * windows
        window_ends = (np.arange(n_windows) + 1) * continue_n * windows - 1
        
        results.append(pd.DataFrame({
            'chrom': chrom,
            'start': window_starts,
            'end': window_ends,
            'short_sum': sum_shorts,
            'inter_sum': sum_inter,
            'long_sum': sum_longs,
            'total_sum': sum_totals
        }))
    
    if not results:
        logger.warning("No valid windows found for FSC")
        # Write empty file with header
        with open(output_path, 'w') as f:
            f.write("region\tshort-fragment-zscore\tintermediate-fragment-zscore\tlong-fragment-zscore\ttotal-fragment-zscore\n")
        return output_path
    
    final_df = pd.concat(results, ignore_index=True)
    logger.info(f"Aggregated into {len(final_df)} windows")
    
    # Compute z-scores
    final_df['short_z'] = _zscore(final_df['short_sum'])
    final_df['inter_z'] = _zscore(final_df['inter_sum'])
    final_df['long_z'] = _zscore(final_df['long_sum'])
    final_df['total_z'] = _zscore(final_df['total_sum'])
    
    # Write output
    with open(output_path, 'w') as f:
        f.write("region\tshort-fragment-zscore\tintermediate-fragment-zscore\tlong-fragment-zscore\ttotal-fragment-zscore\n")
        for _, row in final_df.iterrows():
            region = f"{row['chrom']}:{int(row['start'])}-{int(row['end'])}"
            f.write(f"{region}\t{row['short_z']:.6f}\t{row['inter_z']:.6f}\t{row['long_z']:.6f}\t{row['total_z']:.6f}\n")
    
    logger.info(f"FSC complete: {output_path}")
    return output_path


def _zscore(arr: np.ndarray) -> np.ndarray:
    """Compute z-score of array, handling edge cases."""
    arr = np.asarray(arr, dtype=float)
    mean = np.mean(arr)
    std = np.std(arr)
    if std == 0:
        return np.zeros_like(arr)
    return (arr - mean) / std
