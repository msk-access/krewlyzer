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


def create_panel_ocf_atlas(
    ocf_atlas_path: Path,
    target_regions: Path,
    output_path: Path,
    promoter_extension: int = 2000,
) -> int:
    """
    Create a panel-filtered OCF atlas by intersecting with panel targets.
    
    Pre-filters the genome-wide OCF atlas BED file to only include 
    open chromatin regions that overlap with panel targets (extended 
    by promoter region). This filtered atlas is used as INPUT to the 
    OCF analysis, reducing noise by focusing on relevant tissue-specific 
    OCRs for the panel.
    
    For a panel with ~600 gene targets, this typically reduces the atlas
    from ~50,000 genome-wide OCRs to ~500-2000 panel-relevant OCRs.
    
    Args:
        ocf_atlas_path: Path to genome-wide OCF atlas BED (e.g., 7specificTissue.all.OC.bed.gz)
        target_regions: Panel targets BED file  
        output_path: Path for filtered atlas output (typically .bed.gz)
        promoter_extension: bp to extend targets upstream (default: 2kb)
        
    Returns:
        Number of filtered OCF regions
        
    Example:
        >>> n = create_panel_ocf_atlas(
        ...     Path("data/OpenChromatinRegion/GRCh37/7specificTissue.all.OC.bed.gz"),
        ...     Path("panel_targets.bed"),
        ...     Path("output/panel_ocf_atlas.bed.gz")
        ... )
        >>> print(f"Created panel atlas with {n} regions")
    """
    import gzip
    import tempfile
    
    logger.info(f"Creating panel OCF atlas (+{promoter_extension}bp promoter extension)...")
    
    if not ocf_atlas_path.exists():
        logger.warning(f"OCF atlas file not found: {ocf_atlas_path}")
        return 0
    
    if not target_regions.exists():
        logger.warning(f"Target regions file not found: {target_regions}")
        return 0
    
    # Load panel targets
    targets_df = pd.read_csv(
        target_regions, sep='\t', header=None, comment='#',
        usecols=[0, 1, 2], names=['chrom', 'start', 'end']
    )
    
    if targets_df.empty:
        logger.warning("Target regions file is empty")
        return 0
    
    # Normalize chromosome names (handle chr1 vs 1)
    def normalize_chrom(chrom: str) -> str:
        """Normalize chromosome to 'chr' format for consistency."""
        chrom = str(chrom)
        if not chrom.startswith('chr'):
            return f'chr{chrom}'
        return chrom
    
    targets_df['chrom'] = targets_df['chrom'].apply(normalize_chrom)
    
    # Extend targets for promoter capture
    targets_df['start'] = (targets_df['start'] - promoter_extension).clip(lower=0)
    
    # Build interval lookup by chromosome for fast intersection
    target_intervals = {}
    for chrom in targets_df['chrom'].unique():
        chrom_df = targets_df[targets_df['chrom'] == chrom]
        intervals = list(zip(chrom_df['start'], chrom_df['end']))
        target_intervals[chrom] = sorted(intervals)
    
    def overlaps_any_target(chrom: str, start: int, end: int) -> bool:
        """Check if OCR overlaps any target interval on chromosome."""
        # Normalize chromosome for comparison
        chrom_norm = normalize_chrom(chrom)
        if chrom_norm not in target_intervals:
            return False
        for t_start, t_end in target_intervals[chrom_norm]:
            if end > t_start and start < t_end:
                return True
        return False
    
    # Process OCF atlas BED (may be gzipped)
    n_total = 0
    n_kept = 0
    
    is_gzipped = str(ocf_atlas_path).endswith('.gz')
    opener = gzip.open if is_gzipped else open
    
    output_is_gzipped = str(output_path).endswith('.gz')
    out_opener = gzip.open if output_is_gzipped else open
    
    with opener(ocf_atlas_path, 'rt') as fin, out_opener(output_path, 'wt') as fout:
        for line in fin:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            n_total += 1
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            
            # Check if this OCR overlaps any panel target (+promoter)
            if overlaps_any_target(chrom, start, end):
                fout.write(line + '\n')
                n_kept += 1
    
    logger.info(f"Panel OCF atlas: {n_total:,} → {n_kept:,} regions ({output_path.name})")
    
    return n_kept
