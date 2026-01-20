"""
PON Validation utilities.

Validates panel configuration between PON model and sample processing,
including on-target rate calculation and config compatibility checks.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple
import logging

logger = logging.getLogger("pon.validation")


@dataclass 
class ValidationResult:
    """Result of PON/panel validation."""
    valid: bool
    warnings: List[str]
    errors: List[str]
    on_target_rate: Optional[float] = None
    
    def __bool__(self) -> bool:
        return self.valid


def validate_panel_config(
    pon_model,
    target_regions: Optional[Path] = None,
    assay: Optional[str] = None,
) -> ValidationResult:
    """
    Validate panel configuration between PON and sample processing.
    
    Checks for mismatches that could cause incorrect normalization:
    - PON built in panel mode but sample processed without --target-regions
    - PON built in WGS mode but sample processed with --target-regions  
    - Assay mismatch (PON assay != sample assay)
    
    Args:
        pon_model: Loaded PonModel instance
        target_regions: Path to target regions BED (None = WGS mode)
        assay: Sample assay code (xs1, xs2, etc.)
        
    Returns:
        ValidationResult with warnings and errors
    """
    warnings = []
    errors = []
    
    sample_panel_mode = target_regions is not None
    pon_panel_mode = getattr(pon_model, 'panel_mode', False)
    pon_assay = getattr(pon_model, 'assay', '')
    
    # Check panel mode mismatch
    if pon_panel_mode and not sample_panel_mode:
        warnings.append(
            f"PON was built in panel mode but sample is processed without --target-regions. "
            f"GC correction may be suboptimal."
        )
    
    if not pon_panel_mode and sample_panel_mode:
        warnings.append(
            f"PON was built in WGS mode but sample has --target-regions. "
            f"On-target GC correction not available from this PON."
        )
    
    # Check assay mismatch
    if assay and pon_assay and assay.lower() != pon_assay.lower():
        errors.append(
            f"Assay mismatch: sample assay '{assay}' != PON assay '{pon_assay}'. "
            f"Use a PON built for the same assay."
        )
    
    # Check for missing on-target baselines
    if sample_panel_mode and pon_panel_mode:
        if not hasattr(pon_model, 'gc_bias_ontarget') or pon_model.gc_bias_ontarget is None:
            warnings.append(
                "PON panel mode is set but gc_bias_ontarget is missing. "
                "On-target GC correction will use off-target baseline."
            )
    
    valid = len(errors) == 0
    
    return ValidationResult(
        valid=valid,
        warnings=warnings,
        errors=errors
    )


def calculate_on_target_rate(
    fragment_bed: Path,
    target_regions: Path,
    sample_size: int = 10000
) -> float:
    """
    Calculate the fraction of fragments overlapping target regions.
    
    Samples fragments for efficiency, providing an estimate of on-target rate.
    
    Args:
        fragment_bed: Path to fragment BED/BED.gz file
        target_regions: Path to target regions BED
        sample_size: Number of fragments to sample (default 10000)
        
    Returns:
        On-target rate as float (0.0 to 1.0)
    """
    import gzip
    
    # Load target regions into interval dict
    targets = {}  # chrom -> [(start, end), ...]
    
    with open(target_regions, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            if chrom not in targets:
                targets[chrom] = []
            targets[chrom].append((start, end))
    
    # Sort each chromosome's intervals
    for chrom in targets:
        targets[chrom].sort()
    
    # Sample fragments and count overlaps
    total = 0
    on_target = 0
    
    open_func = gzip.open if str(fragment_bed).endswith('.gz') else open
    with open_func(fragment_bed, 'rt') as f:
        for line in f:
            if total >= sample_size:
                break
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            
            chrom = fields[0]
            try:
                frag_start = int(fields[1])
                frag_end = int(fields[2])
            except ValueError:
                continue
            
            total += 1
            frag_mid = (frag_start + frag_end) // 2
            
            # Check if fragment overlaps any target
            if chrom in targets:
                for tgt_start, tgt_end in targets[chrom]:
                    if tgt_start > frag_mid:
                        break
                    if tgt_start <= frag_mid < tgt_end:
                        on_target += 1
                        break
    
    if total == 0:
        return 0.0
    
    rate = on_target / total
    logger.info(f"On-target rate: {rate:.2%} ({on_target}/{total} sampled fragments)")
    return rate


def print_validation_warnings(result: ValidationResult, console=None):
    """
    Print validation warnings/errors to console.
    
    Args:
        result: ValidationResult from validate_panel_config()
        console: Optional Rich Console for formatted output
    """
    if console is None:
        from rich.console import Console
        console = Console(stderr=True)
    
    for error in result.errors:
        console.print(f"[bold red]❌ ERROR:[/bold red] {error}")
    
    for warning in result.warnings:
        console.print(f"[yellow]⚠️  WARNING:[/yellow] {warning}")
    
    if result.on_target_rate is not None:
        console.print(f"[cyan]📊 On-target rate:[/cyan] {result.on_target_rate:.1%}")
