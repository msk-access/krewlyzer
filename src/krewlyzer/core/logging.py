"""
Standardized logging configuration for krewlyzer.

Provides consistent Rich-formatted logging across all tools and processors.
"""

import logging
from rich.console import Console
from rich.logging import RichHandler


def get_logger(name: str, level: str = "INFO") -> logging.Logger:
    """
    Get a standardized logger with Rich formatting.
    
    Args:
        name: Logger name (e.g., 'fsc', 'core.fsc_processor')
        level: Logging level (DEBUG, INFO, WARNING, ERROR)
        
    Returns:
        Configured logger instance
    """
    console = Console(stderr=True)
    
    # Configure handler if not already configured
    logger = logging.getLogger(name)
    
    if not logger.handlers:
        handler = RichHandler(
            console=console,
            show_time=True,
            show_path=False,
            rich_tracebacks=True,
        )
        handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(handler)
    
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    
    return logger


def set_verbose(logger: logging.Logger, verbose: bool = True) -> None:
    """
    Set logger to DEBUG level if verbose is True.
    
    Args:
        logger: Logger to configure
        verbose: If True, set to DEBUG; otherwise leave unchanged
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")


# =============================================================================
# STARTUP BANNER LOGGING
# =============================================================================
# Provides a consistent, informative startup banner for all CLI tools.
# Shows command, inputs, configuration, and resolved assets with sources.
# =============================================================================

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, List
import sys


@dataclass
class ResolvedAsset:
    """
    Describes a resolved asset for logging in the startup banner.
    
    Tracks both the resolved path and how it was resolved (source),
    helping users understand where configuration came from.
    
    Attributes:
        name: Display name for the asset (e.g., "PON", "Targets")
        path: Resolved path to the asset, or None if not available
        source: How the asset was resolved:
            - "explicit": User provided via CLI flag
            - "bundled": Auto-loaded from bundled data based on --assay
            - "skipped": Explicitly skipped (e.g., --skip-pon)
            - "none": Not available/configured
        detail: Optional additional info (e.g., "n=45 samples", "1852 regions")
    """
    name: str
    path: Optional[Path]
    source: str  # "explicit", "bundled", "skipped", "none"
    detail: str = ""


def log_startup_banner(
    tool_name: str,
    version: str,
    inputs: Dict[str, str],
    config: Dict[str, str],
    assets: List[ResolvedAsset],
    features: Optional[List[str]] = None,
    logger: logging.Logger = None
) -> None:
    """
    Log a comprehensive startup banner for any krewlyzer CLI tool.
    
    Provides a consistent, informative output showing:
    - Command that was run (for reproducibility)
    - Input files provided by user
    - Resolved configuration (genome, assay, mode)
    - Resolved assets with source tracking (bundled vs explicit)
    - Features enabled (for run-all only)
    
    This helps users understand what files will be used before processing starts,
    which is critical for debugging and ensuring correct configuration.
    
    Args:
        tool_name: Name of the CLI tool (e.g., "run-all", "fsc", "build-pon")
        version: krewlyzer version string
        inputs: Dict of input file names to paths (e.g., {"BAM": "sample.bam"})
        config: Dict of configuration keys to values (e.g., {"Genome": "GRCh37"})
        assets: List of ResolvedAsset objects describing resolved assets
        features: Optional list of enabled features (for run-all only)
        logger: Optional logger instance (defaults to "krewlyzer" logger)
    
    Example:
        >>> log_startup_banner(
        ...     tool_name="fsc",
        ...     version="3.0.0",
        ...     inputs={"Fragments": "sample.bed.gz"},
        ...     config={"Assay": "xs2", "Mode": "Panel"},
        ...     assets=[
        ...         ResolvedAsset("Targets", Path("xs2.targets.bed.gz"), "bundled", "1852 regions"),
        ...         ResolvedAsset("PON", Path("xs2.pon.parquet"), "bundled", "n=45"),
        ...     ]
        ... )
    """
    log = logger or logging.getLogger("krewlyzer")
    
    # === HEADER ===
    log.info("═" * 60)
    log.info(f"KREWLYZER {tool_name.upper()} v{version}")
    log.info("═" * 60)
    
    # === COMMAND LINE ===
    # Show full command for reproducibility
    log.info("")
    log.info("COMMAND:")
    log.info(f"  {' '.join(sys.argv)}")
    
    # === INPUT FILES ===
    if inputs:
        log.info("")
        log.info("INPUT:")
        max_key = max(len(k) for k in inputs.keys()) if inputs else 0
        for name, path in inputs.items():
            log.info(f"  {name:{max_key}s}  {path}")
    
    # === CONFIGURATION ===
    if config:
        log.info("")
        log.info("CONFIGURATION:")
        max_key = max(len(k) for k in config.keys()) if config else 0
        for name, value in config.items():
            log.info(f"  {name:{max_key}s}  {value}")
    
    # === RESOLVED ASSETS ===
    # Show each asset with its source (bundled, explicit, etc.)
    if assets:
        log.info("")
        log.info("RESOLVED ASSETS:")
        max_name = max(len(a.name) for a in assets) if assets else 0
        
        for asset in assets:
            # Format path
            if asset.path is None:
                path_str = "None"
            else:
                path_str = asset.path.name if hasattr(asset.path, 'name') else str(asset.path)
            
            # Format source tag
            source_tags = {
                "bundled": "(bundled)",
                "explicit": "(explicit)",
                "skipped": "(skipped)",
            }
            src_tag = source_tags.get(asset.source, "")
            
            # Format optional detail
            detail = f" [{asset.detail}]" if asset.detail else ""
            
            log.info(f"  {asset.name:{max_name}s}  {path_str} {src_tag}{detail}")
    
    # === FEATURES (run-all only) ===
    if features:
        log.info("")
        log.info("FEATURES:")
        for feat in features:
            log.info(f"  ✓ {feat}")
    
    # === FOOTER ===
    log.info("")
    log.info("═" * 60)
