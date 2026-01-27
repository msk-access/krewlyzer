"""
Asset resolution for krewlyzer CLI tools.

Provides centralized logic for resolving bundled assets from CLI arguments:
- Target regions (for panel mode)
- PON models (for z-score normalization)

This module ensures consistent auto-loading behavior across all CLI tools,
whether invoked via `run-all` or individually.

Usage:
    from krewlyzer.core.asset_resolution import resolve_target_regions, resolve_pon_model
    
    # In any CLI tool:
    target_path, target_source = resolve_target_regions(
        explicit_path=target_regions,
        assay=assay,
        skip_target_regions=skip_target_regions,
        assets=assets
    )
"""

import logging
import gzip
from pathlib import Path
from typing import Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from krewlyzer.assets import AssetManager
    from krewlyzer.pon.model import PonModel

# Module logger
logger = logging.getLogger(__name__)


# =============================================================================
# TARGET REGIONS RESOLUTION
# =============================================================================

def resolve_target_regions(
    explicit_path: Optional[Path],
    assay: Optional[str],
    skip_target_regions: bool,
    assets: "AssetManager",
    log: Optional[logging.Logger] = None
) -> Tuple[Optional[Path], str]:
    """
    Resolve target regions from CLI arguments and bundled assets.
    
    Implements the following priority:
    1. Explicit --target-regions path (if valid) → "explicit"
    2. --skip-target-regions flag → None, "skipped"
    3. Bundled targets for --assay (if available) → "bundled"
    4. No targets available → None, "none"
    
    For --assay wgs, no bundled targets exist (by design), so it returns
    (None, "none") rather than an error.
    
    Args:
        explicit_path: Path from --target-regions CLI option, or None
        assay: Assay code from --assay option (xs1, xs2, wgs, etc.), or None
        skip_target_regions: True if --skip-target-regions was specified
        assets: AssetManager instance for loading bundled targets
        log: Optional logger for debug output (defaults to module logger)
    
    Returns:
        Tuple of (resolved_path, source) where:
        - resolved_path: Path to target BED file, or None for WGS mode
        - source: One of "explicit", "bundled", "skipped", or "none"
    
    Raises:
        ValueError: If both --target-regions and --skip-target-regions are provided
    
    Example:
        >>> path, source = resolve_target_regions(None, "xs2", False, assets)
        >>> print(f"{path.name} ({source})")
        xs2.targets.bed.gz (bundled)
    """
    log = log or logger
    
    # ==========================================================================
    # VALIDATION: Check for mutually exclusive options
    # ==========================================================================
    if explicit_path and skip_target_regions:
        raise ValueError(
            "Cannot use both --target-regions and --skip-target-regions. "
            "Use --target-regions to specify a custom targets file, or "
            "--skip-target-regions to disable panel mode entirely."
        )
    
    # ==========================================================================
    # PRIORITY 1: Explicit --target-regions path
    # ==========================================================================
    if explicit_path:
        if explicit_path.exists():
            log.debug(f"Target resolution: Using explicit path: {explicit_path}")
            return explicit_path, "explicit"
        else:
            # Path specified but doesn't exist - this is an error condition
            # Return anyway and let caller validate existence
            log.warning(f"Target resolution: Explicit path does not exist: {explicit_path}")
            return explicit_path, "explicit"
    
    # ==========================================================================
    # PRIORITY 2: --skip-target-regions flag
    # ==========================================================================
    if skip_target_regions:
        log.debug("Target resolution: --skip-target-regions specified, running in WGS mode")
        return None, "skipped"
    
    # ==========================================================================
    # PRIORITY 3: Auto-load bundled targets for assay
    # ==========================================================================
    if assay:
        try:
            bundled_path = assets.get_target_bed(assay)
            log.debug(f"Target resolution: Auto-loaded bundled targets for assay '{assay}': {bundled_path}")
            
            # Count regions for logging (optional, don't fail if error)
            try:
                n_regions = _count_bed_regions(bundled_path)
                log.debug(f"Target resolution: Bundled targets contain {n_regions} regions")
            except Exception:
                pass
            
            return bundled_path, "bundled"
            
        except FileNotFoundError:
            # No bundled targets for this assay (e.g., wgs)
            # This is expected and not an error
            log.debug(f"Target resolution: No bundled targets for assay '{assay}' (treating as WGS)")
            return None, "none"
    
    # ==========================================================================
    # DEFAULT: No targets configured
    # ==========================================================================
    log.debug("Target resolution: No target regions configured, running in WGS mode")
    return None, "none"


def _count_bed_regions(path: Path) -> int:
    """Count number of regions in a BED file (supports .gz)."""
    if str(path).endswith('.gz'):
        with gzip.open(path, 'rt') as f:
            return sum(1 for line in f if not line.startswith('#'))
    else:
        with open(path, 'r') as f:
            return sum(1 for line in f if not line.startswith('#'))


# =============================================================================
# PON MODEL RESOLUTION
# =============================================================================

def resolve_pon_model(
    explicit_path: Optional[Path],
    assay: Optional[str],
    skip_pon: bool,
    assets: "AssetManager",
    log: Optional[logging.Logger] = None
) -> Tuple[Optional[Path], str]:
    """
    Resolve PON model from CLI arguments and bundled assets.
    
    Implements the following priority:
    1. Explicit --pon-model path (if valid) → "explicit"
    2. --skip-pon flag → None, "skipped"
    3. Bundled PON for --assay (if available) → "bundled"
    4. No PON available → None, "none"
    
    Args:
        explicit_path: Path from --pon-model CLI option, or None
        assay: Assay code from --assay option (xs1, xs2, wgs, etc.), or None
        skip_pon: True if --skip-pon was specified
        assets: AssetManager instance for loading bundled PON
        log: Optional logger for debug output (defaults to module logger)
    
    Returns:
        Tuple of (resolved_path, source) where:
        - resolved_path: Path to PON parquet file, or None
        - source: One of "explicit", "bundled", "skipped", or "none"
    
    Raises:
        ValueError: If both --pon-model and --skip-pon are provided
    
    Example:
        >>> path, source = resolve_pon_model(None, "xs2", False, assets)
        >>> print(f"{path.name} ({source})")
        xs2.all_unique.pon.parquet (bundled)
    """
    log = log or logger
    
    # ==========================================================================
    # VALIDATION: Check for mutually exclusive options
    # ==========================================================================
    if explicit_path and skip_pon:
        raise ValueError(
            "Cannot use both --pon-model and --skip-pon. "
            "Use --pon-model to specify a PON for z-score normalization, or "
            "--skip-pon to output raw features without normalization."
        )
    
    # ==========================================================================
    # PRIORITY 1: Explicit --pon-model path
    # ==========================================================================
    if explicit_path:
        if explicit_path.exists():
            log.debug(f"PON resolution: Using explicit path: {explicit_path}")
            return explicit_path, "explicit"
        else:
            # Path specified but doesn't exist
            log.warning(f"PON resolution: Explicit path does not exist: {explicit_path}")
            return explicit_path, "explicit"
    
    # ==========================================================================
    # PRIORITY 2: --skip-pon flag
    # ==========================================================================
    if skip_pon:
        log.debug("PON resolution: --skip-pon specified, outputting raw features")
        return None, "skipped"
    
    # ==========================================================================
    # PRIORITY 3: Auto-load bundled PON for assay
    # ==========================================================================
    if assay:
        try:
            bundled_path = assets.get_pon(assay)
            if bundled_path and bundled_path.exists():
                log.debug(f"PON resolution: Auto-loaded bundled PON for assay '{assay}': {bundled_path}")
                return bundled_path, "bundled"
            else:
                log.debug(f"PON resolution: Bundled PON path invalid for assay '{assay}'")
                return None, "none"
            
        except FileNotFoundError:
            # No bundled PON for this assay
            log.debug(f"PON resolution: No bundled PON for assay '{assay}'")
            return None, "none"
    
    # ==========================================================================
    # DEFAULT: No PON configured
    # ==========================================================================
    log.debug("PON resolution: No PON model configured, outputting raw features")
    return None, "none"


# =============================================================================
# VALIDATION HELPERS
# =============================================================================

def validate_mutual_exclusivity(
    option_name: str,
    option_value,
    skip_flag_name: str,
    skip_flag_value: bool
) -> None:
    """
    Validate that an option and its skip flag are not both specified.
    
    Args:
        option_name: Name of the option (e.g., "--target-regions")
        option_value: Value of the option (truthy means specified)
        skip_flag_name: Name of the skip flag (e.g., "--skip-target-regions")
        skip_flag_value: Value of the skip flag
    
    Raises:
        ValueError: If both are specified
    """
    if option_value and skip_flag_value:
        raise ValueError(f"Cannot use both {option_name} and {skip_flag_name}")


def get_target_region_count(path: Optional[Path]) -> Optional[int]:
    """
    Get the number of target regions in a BED file.
    
    Args:
        path: Path to BED file (supports .gz), or None
    
    Returns:
        Number of regions, or None if path is None or error occurs
    """
    if path is None or not path.exists():
        return None
    try:
        return _count_bed_regions(path)
    except Exception:
        return None
