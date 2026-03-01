"""
Shared utility functions for krewlyzer.

Contains type resolution helpers for handling typer.OptionInfo objects
when functions decorated with typer are called directly (not via CLI).
"""

from pathlib import Path
from typing import Optional


def resolve_path(value) -> Optional[Path]:
    """Safely resolve a path value, handling typer.OptionInfo objects.

    When a typer-decorated function is called directly (not via CLI),
    Optional[Path] parameters with defaults remain as OptionInfo objects.
    This helper handles that case.

    Args:
        value: A Path, str, None, or typer.OptionInfo

    Returns:
        Path if valid, None otherwise
    """
    if value is None:
        return None
    if isinstance(value, Path):
        return value
    if isinstance(value, str):
        return Path(value)
    # If it's a typer.OptionInfo (or any other type), return None
    return None


def resolve_int(value, default: int) -> int:
    """Safely resolve an integer value, handling typer.OptionInfo objects.

    When a typer-decorated function is called directly (not via CLI),
    parameters with defaults remain as OptionInfo objects.
    This helper handles that case.

    Args:
        value: An int or typer.OptionInfo
        default: Default value if not int

    Returns:
        int value or default
    """
    if isinstance(value, int):
        return value
    return default


def resolve_bool(value, default: bool) -> bool:
    """Safely resolve a boolean value, handling typer.OptionInfo objects.

    Args:
        value: A bool or typer.OptionInfo
        default: Default value if not bool

    Returns:
        bool value or default
    """
    if isinstance(value, bool):
        return value
    return default


def resolve_str(value) -> Optional[str]:
    """Safely resolve a string value, handling typer.OptionInfo objects.

    Args:
        value: A str, None, or typer.OptionInfo

    Returns:
        str if valid, None otherwise
    """
    if value is None:
        return None
    if isinstance(value, str):
        return value
    return None
