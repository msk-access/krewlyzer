"""
Core processing module for krewlyzer.

This module provides shared functionality used by both standalone CLI tools
and the run-all command, ensuring consistent behavior and output formats.

Submodules:
    - bam_utils: BAM filter compatibility checking
    - logging: Standardized logging configuration
    - pon_integration: PON model loading and z-score computation
    - fsc_processor: FSC aggregation and z-score calculation
    - fsr_processor: FSR ratio calculation
    - motif_processor: Motif file writing (EDM, BPM, MDS)
    - fsd_processor: FSD PON z-score overlay
    - wps_processor: WPS PON z-score overlay
"""

from .logging import get_logger, set_verbose
from .bam_utils import check_bam_compatibility
from . import pon_integration
from . import fsc_processor
from . import fsr_processor
from . import motif_processor
from . import fsd_processor
from . import wps_processor

__all__ = [
    'get_logger',
    'set_verbose',
    'check_bam_compatibility',
    'pon_integration',
    'fsc_processor',
    'fsr_processor',
    'motif_processor',
    'fsd_processor',
    'wps_processor',
]

