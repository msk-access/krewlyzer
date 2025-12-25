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
