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


def filter_ocf_to_panel(
    genome_ocf_path: Path,
    target_regions: Path,
    output_path: Path,
    promoter_extension: int = 2000,
) -> int:
    """
    Filter genome-wide OCF results to panel-relevant regions.
    
    Creates a focused OCF output containing only regions that overlap
    with panel targets (extended by promoter region). The genome-wide
    output is PRESERVED; this creates an ADDITIONAL panel-focused output.
    
    Args:
        genome_ocf_path: Full genome OCF.tsv from Rust pipeline
        target_regions: Panel targets BED file
        output_path: Path for filtered output
        promoter_extension: bp to extend targets upstream (default: 2kb)
        
    Returns:
        Number of filtered regions
        
    Example:
        >>> n = filter_ocf_to_panel(
        ...     Path("sample.OCF.tsv"),
        ...     Path("panel_targets.bed"),
        ...     Path("sample.OCF.panel.tsv")
        ... )
        >>> print(f"Filtered to {n} panel-relevant regions")
    """
    logger.info(f"Filtering OCF to panel targets (+{promoter_extension}bp promoter)...")
    
    if not genome_ocf_path.exists():
        logger.warning(f"Genome OCF file not found: {genome_ocf_path}")
        return 0
    
    if not target_regions.exists():
        logger.warning(f"Target regions file not found: {target_regions}")
        return 0
    
    # Load genome-wide OCF results
    ocf_df = pd.read_csv(genome_ocf_path, sep='\t')
    
    if ocf_df.empty:
        logger.warning("Genome OCF file is empty")
        return 0
    
    # Load panel targets
    targets_df = pd.read_csv(
        target_regions, sep='\t', header=None,
        usecols=[0, 1, 2], names=['chrom', 'start', 'end']
    )
    
    # Extend targets for promoter capture
    targets_df['start'] = (targets_df['start'] - promoter_extension).clip(lower=0)
    
    # Parse region from OCF (format: "chr1:1000-2000" or similar)
    # Check for various column naming conventions
    region_col = None
    for col_name in ['region', 'region_id', 'name']:
        if col_name in ocf_df.columns:
            region_col = col_name
            break
    
    if region_col is None:
        logger.warning("No region column found in OCF file")
        return 0
    
    # Parse coordinates from region string
    def parse_region(region_str):
        """Parse 'chr1:1000-2000' format."""
        try:
            if ':' in str(region_str) and '-' in str(region_str):
                chrom, coords = str(region_str).split(':')
                start, end = coords.split('-')
                return chrom, int(start), int(end)
        except (ValueError, AttributeError):
            pass
        return None, None, None
    
    parsed = ocf_df[region_col].apply(parse_region)
    ocf_df['ocf_chrom'] = [p[0] for p in parsed]
    ocf_df['ocf_start'] = [p[1] for p in parsed]
    ocf_df['ocf_end'] = [p[2] for p in parsed]
    
    # Filter to valid parsed regions
    ocf_df = ocf_df.dropna(subset=['ocf_chrom'])
    
    if ocf_df.empty:
        logger.warning("No valid regions parsed from OCF file")
        return 0
    
    # Perform intersection
    logger.debug(f"Intersecting {len(ocf_df)} OCF regions with {len(targets_df)} targets...")
    
    filtered_indices = set()
    for chrom in targets_df['chrom'].unique():
        target_chr = targets_df[targets_df['chrom'] == chrom]
        ocf_chr = ocf_df[ocf_df['ocf_chrom'] == chrom]
        
        for _, target in target_chr.iterrows():
            # Find overlapping OCF regions
            overlaps = ocf_chr[
                (ocf_chr['ocf_end'] > target['start']) & 
                (ocf_chr['ocf_start'] < target['end'])
            ]
            filtered_indices.update(overlaps.index.tolist())
    
    if not filtered_indices:
        logger.info("No OCF regions overlap with panel targets")
        # Write empty file with header
        pd.DataFrame(columns=ocf_df.columns).to_csv(output_path, sep='\t', index=False)
        return 0
    
    # Get filtered dataframe, drop temp columns
    result_df = ocf_df.loc[list(filtered_indices)].copy()
    result_df = result_df.drop(columns=['ocf_chrom', 'ocf_start', 'ocf_end'], errors='ignore')
    
    # Write output
    result_df.to_csv(output_path, sep='\t', index=False)
    
    logger.info(f"Panel OCF: {len(ocf_df)} → {len(result_df)} regions ({output_path.name})")
    
    return len(result_df)
