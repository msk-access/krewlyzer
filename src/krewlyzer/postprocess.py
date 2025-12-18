"""
Post-processing functions for Krewlyzer.
Refactored to support both Legacy (CLI) and Unified (Single-Pass) backends.

Note: GC correction is now handled by the Rust backend (v0.3.1+).
Callers should pass already-corrected data to these functions.
"""
import pandas as pd
import numpy as np
import logging
from pathlib import Path

logger = logging.getLogger("postprocess")

def process_fsc_from_counts(
    counts_df: pd.DataFrame, 
    output_path: Path, 
    windows: int, 
    continue_n: int
):
    """
    Process raw fragment counts into FSC Z-scores.
    
    Note: GC correction should be applied BEFORE calling this function.
    Use _core.count_fragments_gc_corrected() for Rust-based LOESS correction.
    
    Args:
        counts_df: DataFrame with columns [chrom, start, end, short, intermediate, long, total]
        output_path: Path to write result
        windows: Window size (e.g. 100000)
        continue_n: Number of bins to aggregate (e.g. 50)
    """
    # Extract values (expecting pre-corrected data from Rust)
    short = counts_df['short'].values
    intermediate = counts_df['intermediate'].values
    long = counts_df['long'].values
    total = counts_df['total'].values
        
    # Create normalized DF
    df = pd.DataFrame({
        'chrom': counts_df['chrom'],
        'start': counts_df['start'],
        'end': counts_df['end'],
        'shorts': short,
        'intermediates': intermediate,
        'longs': long,
        'totals': total
    })
    
    # Aggregation
    results = []
    for chrom, group in df.groupby('chrom', sort=False):
        n_bins = len(group)
        n_windows = n_bins // continue_n
        
        if n_windows == 0:
            continue
        
        trunc_len = n_windows * continue_n
        
        shorts_mat = group['shorts'].values[:trunc_len].reshape(n_windows, continue_n)
        inter_mat = group['intermediates'].values[:trunc_len].reshape(n_windows, continue_n)
        longs_mat = group['longs'].values[:trunc_len].reshape(n_windows, continue_n)
        totals_mat = group['totals'].values[:trunc_len].reshape(n_windows, continue_n)
        
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
        logger.warning("No valid windows found for FSC.")
        # Write header only?
        return

    final_df = pd.concat(results, ignore_index=True)
    
    # Z-scores
    final_df['short_z'] = (final_df['short_sum'] - final_df['short_sum'].mean()) / final_df['short_sum'].std()
    final_df['inter_z'] = (final_df['inter_sum'] - final_df['inter_sum'].mean()) / final_df['inter_sum'].std()
    final_df['long_z'] = (final_df['long_sum'] - final_df['long_sum'].mean()) / final_df['long_sum'].std()
    final_df['total_z'] = (final_df['total_sum'] - final_df['total_sum'].mean()) / final_df['total_sum'].std()
    
    # Write
    with open(output_path, 'w') as f:
        f.write("region\tshort-fragment-zscore\titermediate-fragment-zscore\tlong-fragment-zscore\ttotal-fragment-zscore\n")
        for _, row in final_df.iterrows():
            region = f"{row['chrom']}:{int(row['start'])}-{int(row['end'])}"
            f.write(f"{region}\t{row['short_z']:.4f}\t{row['inter_z']:.4f}\t{row['long_z']:.4f}\t{row['total_z']:.4f}\n")
    
    logger.info(f"FSC processed: {output_path}")


def process_fsr_from_counts(
    counts_df: pd.DataFrame, 
    output_path: Path, 
    windows: int, 
    continue_n: int
):
    """
    Process raw fragment counts into FSR Ratios.
    
    Args:
        counts_df: DataFrame with columns [chrom, start, end, ultra_short, short, intermediate, long, total]
        output_path: Path to write result
        windows: Window size
        continue_n: Aggregation factor
    """
    
    # Aggregation
    results = []
    
    # Renaming for consistency with logic below if inputs differ
    # Expecting: ultra_short, short, intermediate, long, total
    
    for chrom, group in counts_df.groupby('chrom', sort=False):
        n_bins = len(group)
        n_windows = n_bins // continue_n
        
        if n_windows == 0:
            continue
        
        trunc_len = n_windows * continue_n
        
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
            
            u = int(ultra_sums[i])
            s = int(short_sums[i])
            m = int(inter_sums[i])
            l = int(long_sums[i])
            t = int(total_sums[i])
            
            if t > 0:
                s_r = s / t
                m_r = m / t
                l_r = l / t
                u_r = u / t
            else:
                s_r = m_r = l_r = u_r = 0.0
                
            if l > 0:
                sl_r = s / l
            else:
                sl_r = s if s > 0 else 0.0
            
            results.append({
                'region': region,
                'ultra_short_count': u,
                'short_count': s,
                'inter_count': m,
                'long_count': l,
                'total_count': t,
                'short_ratio': s_r,
                'inter_ratio': m_r,
                'long_ratio': l_r,
                'short_long_ratio': sl_r,
                'ultra_short_ratio': u_r
            })
            
    if not results:
        logger.warning("No valid windows found for FSR.")
        return

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, sep='\t', index=False, float_format='%.6f')
    
    logger.info(f"FSR processed: {output_path}")
