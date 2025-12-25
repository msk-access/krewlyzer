"""
FSD (Fragment Size Distribution) processor.

Provides FSD-specific processing including PON z-score overlay.
Used by both standalone fsd.py and run-all wrapper.py.
"""

from pathlib import Path
from typing import Optional
import pandas as pd
import logging

logger = logging.getLogger("core.fsd_processor")


def apply_fsd_pon(
    fsd_output: Path,
    pon,
    output_path: Optional[Path] = None
) -> Path:
    """
    Apply PON z-scores to FSD output file.
    
    Adds z-score columns to the existing FSD output based on the PON baseline.
    
    Args:
        fsd_output: Path to existing FSD.tsv output
        pon: Loaded PonModel with fsd_baseline
        output_path: Optional output path (defaults to overwriting fsd_output)
        
    Returns:
        Path to the output file
    """
    if pon is None or pon.fsd_baseline is None:
        logger.debug("No PON FSD baseline available, skipping z-score overlay")
        return fsd_output
    
    if not fsd_output.exists():
        logger.warning(f"FSD output not found: {fsd_output}")
        return fsd_output
    
    output_path = output_path or fsd_output
    
    logger.info(f"Applying PON z-scores to FSD: {fsd_output}")
    
    # Read existing FSD output
    df = pd.read_csv(fsd_output, sep='\t')
    
    # Expected columns: chrom_arm, size_bin, count, ratio (or similar)
    # The exact schema depends on FSD output format
    
    if 'arm' not in df.columns and 'chrom_arm' not in df.columns:
        logger.warning("FSD output missing arm column, cannot apply PON z-scores")
        return fsd_output
    
    arm_col = 'arm' if 'arm' in df.columns else 'chrom_arm'
    
    # Compute z-scores for each row
    z_scores = []
    for _, row in df.iterrows():
        arm = row[arm_col]
        size_bin = row.get('size_bin', None)
        value = row.get('ratio', row.get('count', 0))
        
        stats = pon.fsd_baseline.get_stats(arm, size_bin)
        if stats is not None:
            mean, std = stats
            if std > 0:
                z = (value - mean) / std
            else:
                z = 0.0
        else:
            z = float('nan')
        
        z_scores.append(z)
    
    df['z_score'] = z_scores
    
    # Write output
    df.to_csv(output_path, sep='\t', index=False, float_format='%.6f')
    logger.info(f"FSD PON z-scores applied: {output_path}")
    
    return output_path
