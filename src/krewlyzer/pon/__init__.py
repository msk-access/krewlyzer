"""
Panel of Normals (PON) module for Krewlyzer.

Provides unified PON model format for:
- GC bias correction curves
- FSD baseline per chromosome arm
- WPS baseline per transcript region
- Panel configuration validation
"""

from .model import PonModel, GcBiasModel, FsdBaseline, WpsBaseline
from .validation import (
    validate_panel_config, 
    calculate_on_target_rate,
    ValidationResult,
    print_validation_warnings,
)

__all__ = [
    "PonModel", 
    "GcBiasModel", 
    "FsdBaseline", 
    "WpsBaseline",
    "validate_panel_config",
    "calculate_on_target_rate", 
    "ValidationResult",
    "print_validation_warnings",
]
