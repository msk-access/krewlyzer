"""
Region Entropy Processor for krewlyzer.

Post-processes region entropy output from Rust:
- Applies PON Z-score normalization
- Loads/saves TSV output

Used for TFBS (transcription factor binding site) and ATAC (cancer peak) size entropy.
"""

import pandas as pd
from pathlib import Path
from typing import Optional, Dict, Tuple
import logging

logger = logging.getLogger("krewlyzer.core.region_entropy_processor")


class RegionEntropyBaseline:
    """
    Baseline for TFBS or ATAC size entropy.
    
    Stores per-label mean and std of entropy values computed from PON cohort.
    Used to convert raw entropy to Z-score for normalized comparison.
    """
    
    def __init__(self, data: Dict[str, Tuple[float, float]]):
        """
        Initialize baseline with pre-computed statistics.
        
        Args:
            data: Dict[label, (mean_entropy, std_entropy)]
        """
        self.data = data
    
    def get_zscore(self, label: str, value: float) -> float:
        """
        Compute Z-score for a given label and entropy value.
        
        Args:
            label: Region label (e.g., TF name or cancer type)
            value: Raw entropy value
            
        Returns:
            Z-score: (value - mean) / std
        """
        if label not in self.data:
            return 0.0
        mean, std = self.data[label]
        if std <= 1e-6:
            return 0.0
        return (value - mean) / std
    
    @classmethod
    def from_samples(cls, sample_data: list) -> "RegionEntropyBaseline":
        """
        Build baseline from multiple samples.
        
        Args:
            sample_data: List of dicts, each containing label -> entropy value
            
        Returns:
            RegionEntropyBaseline with mean/std per label
        """
        import numpy as np
        
        # Collect all values per label
        label_values: Dict[str, list] = {}
        for sample in sample_data:
            for label, entropy in sample.items():
                if label not in label_values:
                    label_values[label] = []
                label_values[label].append(entropy)
        
        # Compute mean/std per label
        data = {}
        for label, values in label_values.items():
            if len(values) >= 2:
                mean = float(np.mean(values))
                std = float(np.std(values, ddof=1))
                data[label] = (mean, std)
            elif len(values) == 1:
                # Single sample: use value as mean, set std=0 (no normalization)
                data[label] = (float(values[0]), 0.0)
        
        logger.info(f"Built RegionEntropyBaseline with {len(data)} labels")
        return cls(data)


def load_entropy_tsv(path: Path) -> pd.DataFrame:
    """
    Load region entropy TSV file.
    
    Columns: label, count, mean_size, entropy
    
    Args:
        path: Path to TSV file
        
    Returns:
        DataFrame with entropy data
    """
    return pd.read_csv(path, sep='\t')


def process_region_entropy(
    raw_path: Path,
    output_path: Path,
    pon_baseline: Optional[RegionEntropyBaseline] = None,
) -> Path:
    """
    Apply PON normalization to region entropy output.
    
    Adds 'z_score' column: (raw_entropy - pon_mean) / pon_std
    
    Args:
        raw_path: Path to raw entropy TSV from Rust
        output_path: Path to write normalized output
        pon_baseline: Optional PON baseline for Z-score calculation
        
    Returns:
        Path to output file
    """
    df = load_entropy_tsv(raw_path)
    
    if pon_baseline:
        logger.info(f"Applying PON normalization to {raw_path.name}...")
        z_scores = []
        for _, row in df.iterrows():
            z = pon_baseline.get_zscore(row['label'], row['entropy'])
            z_scores.append(z)
        df['z_score'] = z_scores
    else:
        logger.debug(f"No PON baseline, skipping Z-score for {raw_path.name}")
        df['z_score'] = 0.0
    
    # Write output
    df.to_csv(output_path, sep='\t', index=False, float_format='%.4f')
    logger.info(f"Wrote {len(df)} labels to {output_path.name}")
    
    return output_path


def extract_entropy_data(df: pd.DataFrame) -> Dict[str, float]:
    """
    Extract entropy values from DataFrame for PON baseline building.
    
    Args:
        df: DataFrame with columns: label, entropy
        
    Returns:
        Dict[label, entropy] for this sample
    """
    return {row['label']: row['entropy'] for _, row in df.iterrows()}


def filter_regions_by_targets(
    region_bed: Path,
    target_bed: Path,
    output_bed: Path,
) -> Optional[Path]:
    """
    Filter TFBS/ATAC regions to keep only those overlapping panel target regions.
    
    Uses bedtools-style interval intersection. Regions are kept if they overlap
    any target region by at least 1bp.
    
    Args:
        region_bed: Path to TFBS or ATAC regions BED.gz
        target_bed: Path to panel target regions BED
        output_bed: Path to write filtered regions
        
    Returns:
        Path to filtered BED file, or None if no overlaps found
    """
    import gzip
    
    # Load target regions into interval tree for fast lookup
    from collections import defaultdict
    
    # Simple interval tree using sorted lists
    targets = defaultdict(list)  # chrom -> [(start, end), ...]
    
    # Read target regions
    try:
        with open(target_bed, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.strip().split('\t')
                if len(cols) < 3:
                    continue
                chrom = cols[0].lstrip('chr')
                start = int(cols[1])
                end = int(cols[2])
                targets[chrom].append((start, end))
    except Exception as e:
        logger.warning(f"Failed to load target regions: {e}")
        return None
    
    # Sort intervals for binary search
    for chrom in targets:
        targets[chrom].sort()
    
    def overlaps_any(chrom: str, start: int, end: int) -> bool:
        """Check if region overlaps any target."""
        chrom = chrom.lstrip('chr')
        if chrom not in targets:
            return False
        for t_start, t_end in targets[chrom]:
            if start < t_end and end > t_start:  # overlap condition
                return True
        return False
    
    # Filter regions
    kept = 0
    total = 0
    
    try:
        opener = gzip.open if str(region_bed).endswith('.gz') else open
        with opener(region_bed, 'rt') as fin, open(output_bed, 'w') as fout:
            for line in fin:
                if line.startswith('#'):
                    fout.write(line)
                    continue
                total += 1
                cols = line.strip().split('\t')
                if len(cols) < 4:
                    continue
                chrom = cols[0]
                start = int(cols[1])
                end = int(cols[2])
                
                if overlaps_any(chrom, start, end):
                    fout.write(line)
                    kept += 1
    except Exception as e:
        logger.warning(f"Failed to filter regions: {e}")
        return None
    
    if kept == 0:
        logger.warning(f"No regions overlap targets (0/{total})")
        output_bed.unlink(missing_ok=True)
        return None
    
    logger.info(f"Filtered regions: {kept}/{total} overlap targets")
    return output_bed

