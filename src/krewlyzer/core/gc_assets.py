"""
GC Correction Asset Resolution Helper.

Provides a centralized function for resolving GC correction assets
to avoid code duplication across standalone tools.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple
import logging

logger = logging.getLogger("krewlyzer.core.gc_assets")


@dataclass
class GcAssets:
    """Resolved GC correction assets."""
    gc_ref: Optional[Path] = None
    valid_regions: Optional[Path] = None
    factors_out: Optional[Path] = None
    factors_input: Optional[Path] = None
    gc_correct_enabled: bool = False


def resolve_gc_assets(
    assets,
    output_dir: Path,
    sample_name: str,
    bedgz_input: Path,
    gc_correct: bool = True,
    genome: str = "hg19"
) -> GcAssets:
    """
    Resolve GC correction assets with automatic pre-computed factors detection.
    
    This centralizes the GC asset resolution logic used by all standalone tools
    (fsc, fsd, fsr, wps, ocf, etc.)
    
    Args:
        assets: AssetManager instance
        output_dir: Output directory for correction_factors.tsv
        sample_name: Sample name for output file naming
        bedgz_input: Input BED.gz file (checked for adjacent correction_factors.tsv)
        gc_correct: Whether GC correction is requested
        genome: Genome build for logging
        
    Returns:
        GcAssets dataclass with resolved paths and enabled status
    """
    result = GcAssets()
    
    if not gc_correct:
        return result
    
    # First, check for pre-computed correction factors next to input
    potential_factors = bedgz_input.parent / f"{bedgz_input.stem.replace('.bed', '')}.correction_factors.tsv"
    if potential_factors.exists():
        result.factors_input = potential_factors
        result.gc_correct_enabled = True
        logger.info(f"Using pre-computed GC factors: {potential_factors.name}")
        return result
    
    # Resolve bundled GC assets
    try:
        result.gc_ref = assets.resolve("gc_reference")
        result.valid_regions = assets.resolve("valid_regions")
        result.factors_out = output_dir / f"{sample_name}.correction_factors.tsv"
        result.gc_correct_enabled = True
        logger.info(f"GC correction enabled ({genome})")
    except FileNotFoundError as e:
        logger.warning(f"GC correction assets not found: {e}")
        logger.warning("Proceeding without GC correction. Use --no-gc-correct to suppress.")
        result.gc_correct_enabled = False
    
    return result
