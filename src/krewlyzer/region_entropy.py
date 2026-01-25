"""
Region Entropy (TFBS/ATAC Size Entropy) calculation.

Calculates Shannon entropy of fragment size distributions at regulatory regions:
- TFBS: Transcription factor binding sites (808 factors from GTRD)
- ATAC: Cancer-specific ATAC-seq peaks (23 cancer types from TCGA)

Output Files:
    - {sample}.TFBS.tsv: Per-TF entropy scores
    - {sample}.ATAC.tsv: Per-cancer-type entropy scores
    - Panel mode: {sample}.TFBS.ontarget.tsv and {sample}.ATAC.ontarget.tsv
"""

import typer
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("region_entropy")


def region_entropy(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output files"),
    
    # Feature toggles
    tfbs: bool = typer.Option(True, "--tfbs/--no-tfbs", help="Enable TFBS entropy"),
    atac: bool = typer.Option(True, "--atac/--no-atac", help="Enable ATAC entropy"),
    
    # Region overrides
    tfbs_regions: Optional[Path] = typer.Option(None, "--tfbs-regions", help="Custom TFBS regions BED.gz"),
    atac_regions: Optional[Path] = typer.Option(None, "--atac-regions", help="Custom ATAC regions BED.gz"),
    
    # Optional settings
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    gc_factors: Optional[Path] = typer.Option(None, "--gc-factors", "-F", help="GC correction factors TSV"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for z-score computation"),
    skip_pon: bool = typer.Option(False, "--skip-pon", help="Skip PON z-score normalization (for PON samples used as ML negatives)"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED (for panel data)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
):
    """
    Calculate Region Entropy features (TFBS/ATAC size entropy) for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.TFBS.tsv and/or {sample}.ATAC.tsv
    """
    from . import _core
    from .assets import AssetManager
    from .core.region_entropy_processor import process_region_entropy
    
    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("krewlyzer.core.region_entropy_processor").setLevel(logging.DEBUG)
    
    # Input validation
    if not bedgz_input.exists():
        logger.error(f"Input file not found: {bedgz_input}")
        raise typer.Exit(1)
    
    # Validate user-provided override files
    from .core.asset_validation import validate_file, FileSchema
    if tfbs_regions and tfbs_regions.exists():
        logger.debug(f"Validating user-provided TFBS regions: {tfbs_regions}")
        validate_file(tfbs_regions, FileSchema.REGION_BED)
    if atac_regions and atac_regions.exists():
        logger.debug(f"Validating user-provided ATAC regions: {atac_regions}")
        validate_file(atac_regions, FileSchema.REGION_BED)
    if gc_factors and gc_factors.exists():
        logger.debug(f"Validating user-provided GC factors: {gc_factors}")
        validate_file(gc_factors, FileSchema.GC_FACTORS_TSV)
    if target_regions and target_regions.exists():
        logger.debug(f"Validating user-provided target regions: {target_regions}")
        validate_file(target_regions, FileSchema.BED3)
    
    # Derive sample name
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Load assets
    assets = AssetManager(genome)
    
    # Resolve GC correction factors
    gc_str = str(gc_factors) if gc_factors and gc_factors.exists() else None
    
    # Validate: -P and --skip-pon are contradictory
    if pon_model and skip_pon:
        logger.error("--pon-model (-P) and --skip-pon are contradictory")
        raise typer.Exit(1)
    
    # Use PON parquet directly for Z-score normalization (Rust implementation)
    entropy_pon_parquet = pon_model if (pon_model and pon_model.exists() and not skip_pon) else None
    
    if entropy_pon_parquet:
        logger.info(f"Using PON for z-score normalization: {pon_model.name}")
    elif skip_pon:
        logger.info("--skip-pon: skipping PON z-score normalization")
    
    try:
        # TFBS
        if tfbs:
            tfbs_path = tfbs_regions if tfbs_regions else (assets.tfbs_regions if assets.tfbs_available else None)
            if tfbs_path and Path(tfbs_path).exists():
                logger.info(f"Computing TFBS entropy...")
                out_raw = output / f"{sample_name}.TFBS.raw.tsv"
                out_final = output / f"{sample_name}.TFBS.tsv"
                
                _core.region_entropy.run_region_entropy(
                    str(bedgz_input), str(tfbs_path), str(out_raw), gc_str, not verbose
                )
                process_region_entropy(out_raw, out_final, entropy_pon_parquet, "tfbs_baseline")
                out_raw.unlink(missing_ok=True)
                logger.info(f"✅ TFBS: {out_final}")
            else:
                logger.warning("TFBS regions not available for this genome")
        
        # ATAC
        if atac:
            atac_path = atac_regions if atac_regions else (assets.atac_regions if assets.atac_available else None)
            if atac_path and Path(atac_path).exists():
                logger.info(f"Computing ATAC entropy...")
                out_raw = output / f"{sample_name}.ATAC.raw.tsv"
                out_final = output / f"{sample_name}.ATAC.tsv"
                
                _core.region_entropy.run_region_entropy(
                    str(bedgz_input), str(atac_path), str(out_raw), gc_str, not verbose
                )
                process_region_entropy(out_raw, out_final, entropy_pon_parquet, "atac_baseline")
                out_raw.unlink(missing_ok=True)
                logger.info(f"✅ ATAC: {out_final}")
            else:
                logger.warning("ATAC regions not available for this genome")
                
    except FileNotFoundError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    except RuntimeError as e:
        logger.error(f"Region entropy calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
