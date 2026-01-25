"""
Asset Format Validation Module
==============================

Centralized validation for all file formats used by krewlyzer.
Validates user-provided files BEFORE calling Rust code.

This module provides:
1. FileSchema enum defining all supported file types
2. validate_file() - validates a single file
3. validate_assets() - validates all files for a run
4. get_format_hint() - returns example format for error messages

Usage:
------
    from krewlyzer.core.asset_validation import validate_file, FileSchema
    
    # Validate before calling Rust
    validate_file(gene_bed_path, FileSchema.GENE_BED)
    
    # Validate all assets for a run
    validate_assets(gene_bed=path1, targets=path2, ...)

All validation happens in Python at the orchestration layer.
Rust code assumes files are already validated.
"""

import gzip
import logging
import re
from enum import Enum
from pathlib import Path
from typing import Optional, Dict, Any

logger = logging.getLogger("krewlyzer.asset_validation")


# =============================================================================
# FILE SCHEMA DEFINITIONS
# =============================================================================

class FileSchema(Enum):
    """
    Enumeration of all supported file formats.
    
    Each schema defines:
    - Minimum required columns
    - Column names and types
    - Validation rules
    - Example format for error messages
    """
    BED3 = "bed3"
    GENE_BED = "gene_bed"
    ARMS_BED = "arms_bed"
    WPS_ANCHORS = "wps_anchors"
    REGION_BED = "region_bed"          # OCF, TFBS, ATAC regions
    GC_FACTORS_TSV = "gc_factors_tsv"
    TARGETS_BED = "targets_bed"


# Schema specifications with validation rules and examples
SCHEMA_SPECS: Dict[FileSchema, Dict[str, Any]] = {
    FileSchema.BED3: {
        "min_cols": 3,
        "description": "Standard BED3 format",
        "columns": ["chrom", "start", "end"],
        "example": "chr1\t1000\t2000",
        "doc_url": "https://msk-access.github.io/krewlyzer/advanced/architecture/",
    },
    FileSchema.GENE_BED: {
        "min_cols": 4,
        "description": "Gene BED format",
        "columns": ["chrom", "start", "end", "gene", "(name)"],
        "example": "chr1\t1000\t2000\tTP53\texon_01\nchr1\t3000\t4000\tTP53\texon_02",
        "doc_url": "https://msk-access.github.io/krewlyzer/advanced/panel-mode/",
    },
    FileSchema.ARMS_BED: {
        "min_cols": 4,
        "description": "Chromosome arms BED format",
        "columns": ["chrom", "start", "end", "arm"],
        "example": "chr1\t0\t125000000\t1p\nchr1\t125000000\t249250621\t1q",
        "doc_url": "https://msk-access.github.io/krewlyzer/features/fsd/",
    },
    FileSchema.WPS_ANCHORS: {
        "min_cols": 6,
        "description": "WPS anchors BED6 format",
        "columns": ["chrom", "start", "end", "name", "score", "strand"],
        "example": "chr1\t1000\t2000\tGENE1_TSS\t0\t+",
        "doc_url": "https://msk-access.github.io/krewlyzer/features/wps/",
    },
    FileSchema.REGION_BED: {
        "min_cols": 4,
        "description": "Named region BED format (OCF/TFBS/ATAC)",
        "columns": ["chrom", "start", "end", "name"],
        "example": "chr1\t1000\t2000\tRegion_001",
        "doc_url": "https://msk-access.github.io/krewlyzer/features/ocf/",
    },
    FileSchema.GC_FACTORS_TSV: {
        "min_cols": 3,
        "description": "GC correction factors TSV",
        "columns": ["length_bin", "gc_pct", "factor"],
        "example": "Short\t40\t1.05\nShort\t41\t1.03",
        "doc_url": "https://msk-access.github.io/krewlyzer/advanced/gc-correction/",
        "has_header": True,
    },
    FileSchema.TARGETS_BED: {
        "min_cols": 4,
        "description": "Target regions BED format",
        "columns": ["chrom", "start", "end", "name"],
        "example": "chr1\t1000\t2000\tTP53_target_01",
        "doc_url": "https://msk-access.github.io/krewlyzer/advanced/panel-mode/",
    },
}


# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

def validate_file(path: Path, schema: FileSchema, max_lines: int = 10) -> int:
    """
    Validate a file against a schema. Raises ValueError on first error.
    
    Args:
        path: Path to the file to validate
        schema: Expected file schema/format
        max_lines: Maximum data lines to check (default: 10, set to 0 for all)
        
    Returns:
        Number of valid data lines found
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid (with helpful error message)
        
    Example:
        >>> validate_file(Path("genes.bed"), FileSchema.GENE_BED)
        1234  # Number of valid lines
    """
    path = Path(path)
    spec = SCHEMA_SPECS[schema]
    
    # Check file exists
    if not path.exists():
        raise FileNotFoundError(f"{spec['description']} not found: {path}")
    
    logger.info(f"Validating {spec['description']}: {path}")
    
    # Determine reader based on file extension
    open_fn = gzip.open if str(path).endswith('.gz') else open
    
    valid_lines = 0
    checked_lines = 0
    
    try:
        with open_fn(path, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                # Skip comments and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                
                # Check if this is a header line (for TSV files)
                if spec.get('has_header') and valid_lines == 0:
                    valid_lines = 0  # Don't count header
                    continue
                
                fields = line.strip().split('\t')
                
                # Validate column count
                if len(fields) < spec['min_cols']:
                    _raise_validation_error(
                        path, schema, line_num,
                        f"Expected at least {spec['min_cols']} columns, got {len(fields)}"
                    )
                
                # Validate specific columns based on schema
                _validate_columns(fields, schema, line_num, path)
                
                valid_lines += 1
                checked_lines += 1
                
                # Stop after max_lines for efficiency
                if max_lines > 0 and checked_lines >= max_lines:
                    break
    
    except gzip.BadGzipFile:
        raise ValueError(f"Invalid gzip file: {path}")
    except UnicodeDecodeError:
        raise ValueError(f"File is not valid text/UTF-8: {path}")
    
    if valid_lines == 0:
        _raise_validation_error(path, schema, 0, "No valid data lines found")
    
    logger.info(f"  ✓ Valid: {path.name} (checked {checked_lines} lines)")
    return valid_lines


def _validate_columns(fields: list, schema: FileSchema, line_num: int, path: Path):
    """
    Validate column values for a specific schema.
    
    Performs type-specific validation beyond just column count.
    """
    # Common BED validation: start < end, coordinates are numeric
    if schema in (FileSchema.BED3, FileSchema.GENE_BED, FileSchema.ARMS_BED, 
                  FileSchema.WPS_ANCHORS, FileSchema.REGION_BED, FileSchema.TARGETS_BED):
        
        # Column 2: start must be numeric
        try:
            start = int(fields[1])
            if start < 0:
                _raise_validation_error(
                    path, schema, line_num,
                    f"Start coordinate must be >= 0, got {start}"
                )
        except ValueError:
            _raise_validation_error(
                path, schema, line_num,
                f"Start coordinate (column 2) must be numeric, got '{fields[1]}'"
            )
        
        # Column 3: end must be numeric and > start
        try:
            end = int(fields[2])
            if end <= start:
                _raise_validation_error(
                    path, schema, line_num,
                    f"End ({end}) must be greater than start ({start})"
                )
        except ValueError:
            _raise_validation_error(
                path, schema, line_num,
                f"End coordinate (column 3) must be numeric, got '{fields[2]}'"
            )
    
    # Schema-specific validation
    if schema == FileSchema.GENE_BED:
        # Column 4: gene name must be non-empty
        if not fields[3].strip():
            _raise_validation_error(
                path, schema, line_num,
                "Gene name (column 4) cannot be empty"
            )
    
    elif schema == FileSchema.ARMS_BED:
        # Column 4: arm must match pattern like "1p", "1q", "22p", "Xq"
        arm_pattern = r'^(\d{1,2}|X|Y)[pq]$'
        if not re.match(arm_pattern, fields[3]):
            _raise_validation_error(
                path, schema, line_num,
                f"Arm (column 4) must match pattern like '1p', '22q', got '{fields[3]}'"
            )
    
    elif schema == FileSchema.WPS_ANCHORS:
        # Column 6: strand must be +, -, or .
        if len(fields) >= 6 and fields[5] not in ('+', '-', '.'):
            _raise_validation_error(
                path, schema, line_num,
                f"Strand (column 6) must be '+', '-', or '.', got '{fields[5]}'"
            )
    
    elif schema == FileSchema.GC_FACTORS_TSV:
        # Column 2: gc_pct must be 0-100
        try:
            gc_pct = int(fields[1])
            if not 0 <= gc_pct <= 100:
                _raise_validation_error(
                    path, schema, line_num,
                    f"GC percent (column 2) must be 0-100, got {gc_pct}"
                )
        except ValueError:
            _raise_validation_error(
                path, schema, line_num,
                f"GC percent (column 2) must be numeric, got '{fields[1]}'"
            )
        
        # Column 3: factor must be positive number
        try:
            factor = float(fields[2])
            if factor <= 0:
                _raise_validation_error(
                    path, schema, line_num,
                    f"Factor (column 3) must be positive, got {factor}"
                )
        except ValueError:
            _raise_validation_error(
                path, schema, line_num,
                f"Factor (column 3) must be numeric, got '{fields[2]}'"
            )


def _raise_validation_error(path: Path, schema: FileSchema, line_num: int, message: str):
    """
    Raise a ValueError with a helpful error message including format example.
    """
    spec = SCHEMA_SPECS[schema]
    
    error_msg = f"\n{'='*60}\n"
    error_msg += f"VALIDATION ERROR: {spec['description']}\n"
    error_msg += f"{'='*60}\n\n"
    error_msg += f"File: {path}\n"
    if line_num > 0:
        error_msg += f"Line: {line_num}\n"
    error_msg += f"Error: {message}\n\n"
    error_msg += f"Expected columns: {', '.join(spec['columns'])}\n\n"
    error_msg += f"Example format:\n"
    for line in spec['example'].split('\n'):
        error_msg += f"  {line}\n"
    error_msg += f"\nDocumentation: {spec['doc_url']}\n"
    
    raise ValueError(error_msg)


def get_format_hint(schema: FileSchema) -> str:
    """
    Get a formatted hint showing expected file format.
    
    Args:
        schema: The file schema to get hint for
        
    Returns:
        Multi-line string with format description and example
    """
    spec = SCHEMA_SPECS[schema]
    hint = f"{spec['description']}\n"
    hint += f"Required columns: {', '.join(spec['columns'])}\n\n"
    hint += "Example:\n"
    for line in spec['example'].split('\n'):
        hint += f"  {line}\n"
    return hint


# =============================================================================
# BATCH VALIDATION
# =============================================================================

def validate_assets(
    gene_bed: Optional[Path] = None,
    targets_bed: Optional[Path] = None,
    arms_bed: Optional[Path] = None,
    wps_anchors: Optional[Path] = None,
    gc_factors: Optional[Path] = None,
    ocf_regions: Optional[Path] = None,
    tfbs_regions: Optional[Path] = None,
    atac_regions: Optional[Path] = None,
) -> Dict[str, int]:
    """
    Validate multiple asset files. Fails on first error.
    
    Args:
        gene_bed: Path to gene BED file
        targets_bed: Path to targets BED file
        arms_bed: Path to chromosome arms BED file
        wps_anchors: Path to WPS anchors BED file
        gc_factors: Path to GC correction factors TSV
        ocf_regions: Path to OCF regions BED file
        tfbs_regions: Path to TFBS regions BED file
        atac_regions: Path to ATAC regions BED file
        
    Returns:
        Dictionary mapping file type to number of valid lines
        
    Raises:
        ValueError: On first validation error
    """
    logger.info("Validating input assets...")
    
    results = {}
    
    validation_map = [
        (gene_bed, FileSchema.GENE_BED, "gene_bed"),
        (targets_bed, FileSchema.TARGETS_BED, "targets_bed"),
        (arms_bed, FileSchema.ARMS_BED, "arms_bed"),
        (wps_anchors, FileSchema.WPS_ANCHORS, "wps_anchors"),
        (gc_factors, FileSchema.GC_FACTORS_TSV, "gc_factors"),
        (ocf_regions, FileSchema.REGION_BED, "ocf_regions"),
        (tfbs_regions, FileSchema.REGION_BED, "tfbs_regions"),
        (atac_regions, FileSchema.REGION_BED, "atac_regions"),
    ]
    
    for path, schema, name in validation_map:
        if path is not None:
            n_lines = validate_file(Path(path), schema)
            results[name] = n_lines
    
    logger.info(f"✓ All {len(results)} asset files validated successfully")
    return results
