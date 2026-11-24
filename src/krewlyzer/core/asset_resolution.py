"""
Asset resolution for krewlyzer CLI tools.

Provides centralized logic for resolving bundled assets from CLI arguments:
- Target regions (for panel mode)
- PON models (for z-score normalization)
- Exclude regions (the blacklist)

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
    log: Optional[logging.Logger] = None,
) -> Tuple[Optional[Path], str]:
    """
    Resolve target regions from CLI arguments and bundled assets.

    Implements the following priority:
    1. Explicit --target-regions path (if valid) → "explicit"
    2. --skip-target-regions flag → None, "skipped"
    3. Bundled targets for --assay → "bundled"
    4. --assay wgs, or no --assay → None, "none"

    A named assay other than ``wgs`` whose bundled targets are missing is an
    error, not WGS. It used to be taken for WGS at debug level, so a panel
    sample was measured as WGS with ``assay`` and ``panel_mode=False`` both in
    its record and exit 0 -- a data folder that lacks the file, or holds only
    its Git LFS pointer, was enough. ``wgs`` has no bundled targets by design.
    An explicit path that does not exist is an error too.

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
        ValueError: If both --target-regions and --skip-target-regions are
            provided, if an explicit path does not exist, or if a named assay
            other than ``wgs`` has no bundled targets

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
        # One that does not exist is an error. It was warned about and returned for the
        # caller to check, and every caller took `path.exists()` as panel mode: a typo in
        # -T ran the sample as WGS, exit 0.
        if not explicit_path.exists():
            raise ValueError(f"Target regions file not found: {explicit_path}")
        _refuse_an_lfs_pointer(explicit_path, "--target-regions")
        log.debug(f"Target resolution: Using explicit path: {explicit_path}")
        return explicit_path, "explicit"

    # ==========================================================================
    # PRIORITY 2: --skip-target-regions flag
    # ==========================================================================
    if skip_target_regions:
        log.debug(
            "Target resolution: --skip-target-regions specified, running in WGS mode"
        )
        return None, "skipped"

    # ==========================================================================
    # PRIORITY 3: Auto-load bundled targets for assay
    # ==========================================================================
    if assay:
        try:
            bundled_path = assets.get_target_bed(assay)
            log.debug(
                f"Target resolution: Auto-loaded bundled targets for assay '{assay}': {bundled_path}"
            )

            # Count regions for logging (optional, don't fail if error)
            try:
                n_regions = _count_bed_regions(bundled_path)
                log.debug(
                    f"Target resolution: Bundled targets contain {n_regions} regions"
                )
            except Exception:
                pass

            return bundled_path, "bundled"

        except FileNotFoundError as e:
            if assay.lower() == "wgs":
                log.debug("Target resolution: --assay wgs has no targets, WGS mode")
                return None, "none"
            raise ValueError(
                f"No bundled target regions for assay '{assay}' ({e}). Without them a "
                f"panel sample would be measured as WGS. Pass the panel's targets with "
                f"--target-regions, run as WGS with --skip-target-regions or --assay wgs, "
                f"or point KREWLYZER_DATA_DIR at a data clone that has them."
            ) from e

    # ==========================================================================
    # DEFAULT: No targets configured
    # ==========================================================================
    log.debug("Target resolution: No target regions configured, running in WGS mode")
    return None, "none"


def _refuse_an_lfs_pointer(path: Path, option: str) -> None:
    """An explicit file that is a Git LFS pointer, as one taken from a clone that skipped
    `git lfs pull` is: it exists, and reads as a region set of nothing."""
    from krewlyzer.assets import is_lfs_pointer

    if is_lfs_pointer(path):
        raise ValueError(
            f"{option} {path} is a Git LFS pointer, not the file: run `git lfs pull` in "
            f"the clone it came from"
        )


def _count_bed_regions(path: Path) -> int:
    """Count number of regions in a BED file (supports .gz)."""
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return sum(1 for line in f if not line.startswith("#"))
    else:
        with open(path, "r") as f:
            return sum(1 for line in f if not line.startswith("#"))


# =============================================================================
# PON MODEL RESOLUTION
# =============================================================================


def resolve_pon_model(
    explicit_path: Optional[Path],
    assay: Optional[str],
    skip_pon: bool,
    assets: "AssetManager",
    variant: str = "all_unique",
    log: Optional[logging.Logger] = None,
) -> Tuple[Optional[Path], str]:
    """
    Resolve PON model from CLI arguments and bundled assets.

    Implements the following priority:
    1. Explicit --pon-model path (if valid) → "explicit"
    2. --skip-pon flag → None, "skipped"
    3. Bundled PON for --assay + --pon-variant (if available) → "bundled"
    4. No PON available → None, "none"

    Args:
        explicit_path: Path from --pon-model CLI option, or None
        assay: Assay code from --assay option (xs1, xs2, wgs, etc.), or None
        skip_pon: True if --skip-pon was specified
        assets: AssetManager instance for loading bundled PON
        variant: PON variant to use ("all_unique" or "duplex"). Defaults to "all_unique".
                 - all_unique: Maximum coverage (duplex + simplex + singletons)
                 - duplex: Highest accuracy (duplex consensus only)
        log: Optional logger for debug output (defaults to module logger)

    Returns:
        Tuple of (resolved_path, source) where:
        - resolved_path: Path to PON parquet file, or None
        - source: One of "explicit", "bundled", "skipped", or "none"

    Raises:
        ValueError: If both --pon-model and --skip-pon are provided

    Example:
        >>> path, source = resolve_pon_model(None, "xs2", False, assets, variant="duplex")
        >>> print(f"{path.name} ({source})")
        xs2.duplex.pon.parquet (bundled)
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
        # An explicit PON that is not there is an error: every caller used to skip it on
        # `.exists()`, so the run finished unscored and looked like a deliberate --skip-pon.
        # ValueError, like the other bad --pon-model arguments: every tool prints it as an error.
        raise ValueError(
            f"--pon-model {explicit_path.name}: no such file ({explicit_path}). "
            "Pass an existing PON, or --skip-pon to run without one."
        )

    # ==========================================================================
    # PRIORITY 2: --skip-pon flag
    # ==========================================================================
    if skip_pon:
        log.debug("PON resolution: --skip-pon specified, outputting raw features")
        return None, "skipped"

    # ==========================================================================
    # PRIORITY 3: Auto-load bundled PON for assay + variant
    # ==========================================================================
    if assay:
        try:
            bundled_path = assets.get_pon(assay, variant=variant)
            if bundled_path and bundled_path.exists():
                log.info(
                    f"PON resolution: Loaded bundled {variant} PON for '{assay}': {bundled_path.name}"
                )
                return bundled_path, "bundled"
            else:
                log.debug(
                    f"PON resolution: Bundled PON path invalid for assay '{assay}'"
                )
                return None, "none"

        except FileNotFoundError:
            # No bundled PON for this assay/variant
            log.debug(
                f"PON resolution: No bundled PON for assay='{assay}', variant='{variant}'"
            )
            return None, "none"

    # ==========================================================================
    # DEFAULT: No PON configured
    # ==========================================================================
    log.debug("PON resolution: No PON model configured, outputting raw features")
    return None, "none"


# =============================================================================
# EXCLUDE REGIONS (BLACKLIST) RESOLUTION
# =============================================================================


def resolve_exclude_regions(
    explicit_path: Optional[Path],
    assets: "AssetManager",
    log: Optional[logging.Logger] = None,
) -> Tuple[Optional[Path], str]:
    """The blacklist every fragment reader applies: the one given, else the genome's bundled one.

    Shared by `extract`, `run-all` and `motif`, so one sample's tables all exclude the same
    fragments. `motif` excluded none, so its motif tables counted fragments run-all's dropped.

    Returns ``(path, source)``, source ``explicit`` or ``bundled``.

    Raises:
        ValueError: an explicit file that does not exist (it used to be used as given) or is a
            directory, or no explicit file and no bundled blacklist for the genome, or only its
            Git LFS pointer (it used to be a warning, and the run kept the blacklisted fragments
            and exited 0).
    """
    log = log or logger
    if explicit_path is not None:
        if not Path(explicit_path).exists():
            raise ValueError(f"Exclude regions file not found: {explicit_path}")
        if not Path(explicit_path).is_file():
            # Its bytes are the blacklist's digest in the record (`blacklist_digest`).
            raise ValueError(
                f"Exclude regions {explicit_path} is a directory, not a BED file"
            )
        _refuse_an_lfs_pointer(Path(explicit_path), "--exclude-regions")
        log.info(f"Exclude regions: {Path(explicit_path).name} (explicit)")
        return Path(explicit_path), "explicit"
    from krewlyzer.assets import is_lfs_pointer

    bundled = assets.exclude_regions
    ways_out = (
        "Pass one with --exclude-regions (an empty file excludes nothing), or point "
        "KREWLYZER_DATA_DIR at a data clone that has it."
    )
    if bundled.exists():
        if is_lfs_pointer(bundled):
            # A pointer is read as a blacklist of nothing: every fragment kept, exit 0.
            raise ValueError(
                f"The bundled blacklist {bundled} is a Git LFS pointer, not the file: run "
                f"`git lfs pull` in the data clone. {ways_out}"
            )
        log.info(f"Exclude regions: {bundled.name} (bundled)")
        return bundled, "bundled"
    reason = assets.missing_data_reason()
    raise ValueError(
        f"No bundled blacklist at {bundled}"
        + (f" ({reason})" if reason else "")
        + f": without it blacklisted fragments would be kept. {ways_out}"
    )


# =============================================================================
# VALIDATION HELPERS
# =============================================================================


def validate_mutual_exclusivity(
    option_name: str, option_value, skip_flag_name: str, skip_flag_value: bool
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
