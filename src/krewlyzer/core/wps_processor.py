"""
WPS (Windowed Protection Score) processor.

Provides WPS-specific processing including PON z-score overlay.
Used by both standalone wps.py and run-all wrapper.py.
"""

from pathlib import Path
from typing import Optional
import pandas as pd
import logging

logger = logging.getLogger("core.wps_processor")


def apply_wps_pon(
    wps_output: Path,
    pon,
    output_path: Optional[Path] = None
) -> Path:
    """
    Apply PON z-scores to WPS output file.
    
    Adds z-score columns to the existing WPS output based on the PON baseline.
    
    Args:
        wps_output: Path to existing WPS.tsv output
        pon: Loaded PonModel with wps_baseline
        output_path: Optional output path (defaults to overwriting wps_output)
        
    Returns:
        Path to the output file
    """
    if pon is None or pon.wps_baseline is None:
        logger.debug("No PON WPS baseline available, skipping z-score overlay")
        return wps_output
    
    if not wps_output.exists():
        logger.warning(f"WPS output not found: {wps_output}")
        return wps_output
    
    output_path = output_path or wps_output
    
    logger.info(f"Applying PON z-scores to WPS: {wps_output}")
    
    # Read existing WPS output
    df = pd.read_csv(wps_output, sep='\t')
    
    # Expected columns: gene, wps, or similar
    # The exact schema depends on WPS output format
    
    gene_col = None
    for col in ['gene', 'Gene', 'gene_name', 'transcript']:
        if col in df.columns:
            gene_col = col
            break
    
    if gene_col is None:
        logger.warning("WPS output missing gene column, cannot apply PON z-scores")
        return wps_output
    
    value_col = None
    for col in ['wps', 'WPS', 'score', 'value']:
        if col in df.columns:
            value_col = col
            break
    
    if value_col is None:
        logger.warning("WPS output missing value column, cannot apply PON z-scores")
        return wps_output
    
    # Compute z-scores for each row
    z_scores = []
    for _, row in df.iterrows():
        gene = row[gene_col]
        value = row[value_col]
        
        stats = pon.wps_baseline.get_stats(gene)
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
    logger.info(f"WPS PON z-scores applied: {output_path}")
    
    return output_path
