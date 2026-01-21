"""
FSC (Fragment Size Coverage) processor.

Processes GC-weighted fragment counts into ML-ready features per the implementation plan:
- Raw GC-weighted counts per channel (from Rust backend)
- PoN log2 ratios: log2(channel / PoN_mean)
- Reliability scores: 1 / (PoN_variance + k)

NO z-scores - the plan specifies GC-weighting (Rust) + PoN log-ratios (Python).
"""

from pathlib import Path
from typing import Optional, List
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger("core.fsc_processor")

# 5 biologically-meaningful channels matching Rust backend
CHANNELS = ['ultra_short', 'core_short', 'mono_nucl', 'di_nucl', 'long']

# Regularization constant for reliability calculation
RELIABILITY_K = 0.01


def process_fsc(
    counts_df: pd.DataFrame,
    output_path: Path,
    windows: int = 100000,
    continue_n: int = 50,
    pon = None
) -> Path:
    """
    Process GC-weighted fragment counts into FSC ML features.
    
    Per the implementation plan, FSC normalization order is:
    1. GC-weight (Rust) - done upstream
    2. PoN log-ratio (Python) - done here
    
    Output columns:
    - chrom, start, end: Bin coordinates
    - ultra_short, core_short, mono_nucl, di_nucl, long: GC-weighted counts
    - total: All fragments
    - *_log2: log2(channel / PoN_mean) when PoN provided
    - *_reliability: 1 / (PoN_variance + k) when PoN provided
    
    Args:
        counts_df: DataFrame from Rust with [chrom, start, end, ultra_short, 
                   core_short, mono_nucl, di_nucl, long, total, mean_gc]
        output_path: Path to write FSC.tsv output
        windows: Window size in bp (default 100000)
        continue_n: Number of bins to aggregate (default 50)
        pon: Optional PON model for log2 ratio normalization
        
    Returns:
        Path to the output file
    """
    logger.info(f"Processing FSC: {len(counts_df)} bins → {output_path}")
    
    # Validate required columns
    required_cols = ['chrom', 'start', 'end'] + CHANNELS + ['total']
    missing = [c for c in required_cols if c not in counts_df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    # Aggregate bins into windows
    results = []
    for chrom, group in counts_df.groupby('chrom', sort=False):
        n_bins = len(group)
        n_windows = n_bins // continue_n
        
        if n_windows == 0:
            continue
        
        trunc_len = n_windows * continue_n
        
        # Aggregate each channel
        channel_sums = {}
        for ch in CHANNELS:
            mat = group[ch].values[:trunc_len].reshape(n_windows, continue_n)
            channel_sums[ch] = mat.sum(axis=1)
        
        totals_mat = group['total'].values[:trunc_len].reshape(n_windows, continue_n)
        channel_sums['total'] = totals_mat.sum(axis=1)
        
        # Get window coordinates from first and last bin in each window
        starts = group['start'].values[:trunc_len].reshape(n_windows, continue_n)[:, 0]
        ends = group['end'].values[:trunc_len].reshape(n_windows, continue_n)[:, -1]
        
        window_df = pd.DataFrame({
            'chrom': chrom,
            'start': starts,
            'end': ends,
            **channel_sums
        })
        results.append(window_df)
    
    if not results:
        logger.warning("No valid windows found for FSC")
        _write_empty_fsc(output_path, pon is not None)
        return output_path
    
    final_df = pd.concat(results, ignore_index=True)
    logger.info(f"Aggregated into {len(final_df)} windows")
    
    # Compute PoN log2 ratios and reliability if model provided
    if pon is not None:
        logger.info("Computing PoN log2 ratios...")
        for ch in CHANNELS:
            mean = pon.get_mean(ch) if hasattr(pon, 'get_mean') else None
            var = pon.get_variance(ch) if hasattr(pon, 'get_variance') else None
            
            if mean and mean > 0:
                # log2(channel / PoN_mean)
                final_df[f'{ch}_log2'] = np.log2(final_df[ch] / mean + 1e-9)
            else:
                final_df[f'{ch}_log2'] = 0.0
            
            if var is not None:
                # reliability = 1 / (PoN_variance + k)
                final_df[f'{ch}_reliability'] = 1.0 / (var + RELIABILITY_K)
            else:
                final_df[f'{ch}_reliability'] = 1.0
    
    # Write output
    _write_fsc_output(final_df, output_path, pon is not None)
    
    logger.info(f"FSC complete: {output_path}")
    return output_path


def _write_empty_fsc(output_path: Path, with_pon: bool):
    """Write empty FSC file with appropriate headers."""
    cols = ['chrom', 'start', 'end'] + CHANNELS + ['total']
    if with_pon:
        cols += [f'{ch}_log2' for ch in CHANNELS]
        cols += [f'{ch}_reliability' for ch in CHANNELS]
    
    pd.DataFrame(columns=cols).to_csv(output_path, sep='\t', index=False)


def _write_fsc_output(df: pd.DataFrame, output_path: Path, with_pon: bool):
    """Write FSC output TSV matching the implementation plan schema."""
    # Determine columns to output
    cols = ['chrom', 'start', 'end'] + CHANNELS + ['total']
    if with_pon:
        cols += [f'{ch}_log2' for ch in CHANNELS]
        cols += [f'{ch}_reliability' for ch in CHANNELS]
    
    # Select and write
    output_df = df[cols].copy()
    output_df.to_csv(output_path, sep='\t', index=False, float_format='%.4f')


def aggregate_by_gene(
    fragments_bed: Path,
    genes: dict,
    output_path: Path,
    pon=None
) -> Path:
    """
    Aggregate fragment counts by gene for panel-specific FSC.
    
    Instead of genomic windows, this function sums fragment counts
    per gene target regions, producing gene-level FSC features.
    
    Args:
        fragments_bed: Path to fragment BED file from extract step
        genes: Dict[gene_name, List[Region]] from gene_bed.load_gene_bed()
        output_path: Path to write gene FSC output
        pon: Optional PON model for normalization
        
    Returns:
        Path to output file
        
    Output columns:
        - gene: Gene name
        - n_regions: Number of target regions for this gene
        - total_bp: Total basepairs covered
        - ultra_short, core_short, mono_nucl, di_nucl, long: Channel counts
        - total: Total fragments
        - *_ratio: Channel / total ratio
    """
    import gzip
    
    logger.info(f"Aggregating FSC by gene: {len(genes)} genes")
    
    # Build interval tree for fast region lookup
    # Structure: chrom -> [(start, end, gene), ...]
    gene_intervals = {}
    for gene, regions in genes.items():
        for region in regions:
            chrom = region.chrom
            if chrom not in gene_intervals:
                gene_intervals[chrom] = []
            gene_intervals[chrom].append((region.start, region.end, gene))
    
    # Sort intervals per chromosome
    for chrom in gene_intervals:
        gene_intervals[chrom].sort(key=lambda x: x[0])
    
    # Initialize gene counts
    gene_counts = {gene: {ch: 0 for ch in CHANNELS + ['total']} for gene in genes}
    gene_meta = {gene: {'n_regions': len(regions), 'total_bp': sum(r.end - r.start for r in regions)} 
                 for gene, regions in genes.items()}
    
    # Read fragments and count by gene
    total_fragments = 0
    assigned_fragments = 0
    
    open_func = gzip.open if str(fragments_bed).endswith('.gz') else open
    with open_func(fragments_bed, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            
            total_fragments += 1
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            
            # Determine size bin
            frag_len = end - start
            if frag_len < 65:
                continue  # Too short
            elif frag_len < 100:
                size_bin = 'ultra_short'
            elif frag_len < 150:
                size_bin = 'core_short'
            elif frag_len < 260:
                size_bin = 'mono_nucl'
            elif frag_len < 400:
                size_bin = 'di_nucl'
            else:
                size_bin = 'long'
            
            # Find overlapping gene regions
            if chrom not in gene_intervals:
                continue
            
            frag_mid = (start + end) // 2
            for reg_start, reg_end, gene in gene_intervals[chrom]:
                if reg_start > frag_mid:
                    break  # Past this fragment
                if reg_start <= frag_mid < reg_end:
                    gene_counts[gene][size_bin] += 1
                    gene_counts[gene]['total'] += 1
                    assigned_fragments += 1
                    break  # Assign to first matching gene only
    
    logger.info(f"Processed {total_fragments} fragments, {assigned_fragments} assigned to genes")
    
    # Build output DataFrame
    rows = []
    for gene in sorted(genes.keys()):
        row = {
            'gene': gene,
            'n_regions': gene_meta[gene]['n_regions'],
            'total_bp': gene_meta[gene]['total_bp'],
            **gene_counts[gene]
        }
        # Calculate channel ratios
        total = gene_counts[gene]['total']
        for ch in CHANNELS:
            row[f'{ch}_ratio'] = gene_counts[gene][ch] / total if total > 0 else 0.0
        rows.append(row)
    
    df = pd.DataFrame(rows)
    
    # Add PoN normalization if available
    if pon is not None:
        logger.info("Computing PoN log2 ratios for gene FSC...")
        for ch in CHANNELS:
            mean = pon.get_mean(f'gene_{ch}') if hasattr(pon, 'get_mean') else None
            if mean and mean > 0:
                df[f'{ch}_log2'] = np.log2(df[ch] / mean + 1e-9)
            else:
                df[f'{ch}_log2'] = 0.0
    
    # Write output
    df.to_csv(output_path, sep='\t', index=False, float_format='%.4f')
    logger.info(f"Gene FSC complete: {output_path} ({len(df)} genes)")
    
    return output_path

