"""
Region Entropy (TFBS/ATAC Size Entropy) calculation.

Calculates Shannon entropy of fragment size distributions at regulatory regions:
- TFBS: Transcription factor binding sites (808 factors from GTRD)
- ATAC: Cancer-specific ATAC-seq peaks (23 cancer types from TCGA)

Output Files:
    Genome-wide (WGS-comparable, all fragments):
    - {sample}.TFBS.tsv: TFBS entropy scores
    - {sample}.ATAC.tsv: ATAC entropy scores
    
    Panel-specific (uses pre-intersected regions overlapping targets):
    - {sample}.TFBS.ontarget.tsv: TFBS entropy for panel-specific regions
    - {sample}.ATAC.ontarget.tsv: ATAC entropy for panel-specific regions
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

# Import asset resolution and startup banner
from .core.asset_resolution import resolve_target_regions, resolve_pon_model
from .core.logging import log_startup_banner, ResolvedAsset
from . import __version__


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
    skip_target_regions: bool = typer.Option(False, "--skip-target-regions", help="Disable panel mode even when --assay has bundled targets"),
    assay: Optional[str] = typer.Option(None, "--assay", "-A", help="Assay code (xs1/xs2) for auto-loading bundled assets"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0 = use all available cores)"),
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
    
    # ═══════════════════════════════════════════════════════════════════
    # ASSET RESOLUTION
    # ═══════════════════════════════════════════════════════════════════
    try:
        resolved_pon_path, pon_source = resolve_pon_model(
            explicit_path=pon_model, assay=assay, skip_pon=skip_pon,
            assets=assets, log=logger
        )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)
    
    try:
        resolved_target_path, target_source = resolve_target_regions(
            explicit_path=target_regions, assay=assay, skip_target_regions=skip_target_regions,
            assets=assets, log=logger
        )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)
    
    is_panel_mode = resolved_target_path is not None and resolved_target_path.exists()
    target_regions_str = str(resolved_target_path) if is_panel_mode else None
    
    # Resolve GC correction factors
    gc_str = str(gc_factors) if gc_factors and gc_factors.exists() else None
    
    # Use PON parquet directly for Z-score normalization (Rust implementation)
    entropy_pon_parquet = resolved_pon_path if (resolved_pon_path and resolved_pon_path.exists()) else None
    
    # Build list of enabled features
    features = []
    if tfbs:
        features.append("TFBS Entropy")
    if atac:
        features.append("ATAC Entropy")
    
    # ═══════════════════════════════════════════════════════════════════
    # STARTUP BANNER
    # ═══════════════════════════════════════════════════════════════════
    log_startup_banner(
        tool_name="region-entropy", version=__version__,
        inputs={"Fragments": str(bedgz_input.name), "Output": str(output)},
        config={"Genome": f"{assets.raw_genome} → {assets.genome_dir}", "Assay": assay or "None", "Mode": "Panel" if is_panel_mode else "WGS"},
        assets=[ResolvedAsset("PON", resolved_pon_path, pon_source), ResolvedAsset("Targets", resolved_target_path, target_source)],
        features=features,
        logger=logger
    )
    
    # Configure thread pool for Rayon parallelization
    if threads > 0:
        try:
            _core.configure_threads(threads)
            logger.debug(f"Configured {threads} threads for parallel processing")
        except RuntimeError:
            # Thread pool already configured (can only be set once per process)
            pass
    
    try:
        # TFBS
        if tfbs:
            # Always use genome-wide regions for primary output (WGS-comparable)
            tfbs_path_gw = tfbs_regions if tfbs_regions else (assets.tfbs_regions if assets.tfbs_available else None)
            
            if tfbs_path_gw and Path(tfbs_path_gw).exists():
                logger.info(f"Computing TFBS entropy (genome-wide)...")
                out_raw = output / f"{sample_name}.TFBS.raw.tsv"
                out_final = output / f"{sample_name}.TFBS.tsv"
                
                # Call Rust with optional target_regions for on/off-target split
                n_off, n_on = _core.region_entropy.run_region_entropy(
                    str(bedgz_input), str(tfbs_path_gw), str(out_raw), gc_str, 
                    target_regions_str, not verbose
                )
                
                # Process off-target output (primary)
                process_region_entropy(out_raw, out_final, entropy_pon_parquet, "tfbs_baseline")
                out_raw.unlink(missing_ok=True)
                logger.info(f"✅ TFBS: {out_final} ({n_off} TFs)")
                
            else:
                logger.warning("TFBS regions not available for this genome")
            
            # Panel-specific TFBS (using panel-intersection regions)
            # Output renamed to .ontarget.tsv for consistency with other features
            if is_panel_mode and assay:
                tfbs_path_panel = assets.get_tfbs_regions(assay)
                if tfbs_path_panel and tfbs_path_panel.exists() and tfbs_path_panel != tfbs_path_gw:
                    logger.info(f"Computing TFBS entropy (panel-specific: {assay})...")
                    out_raw_ont = output / f"{sample_name}.TFBS.ontarget.raw.tsv"
                    out_final_ont = output / f"{sample_name}.TFBS.ontarget.tsv"
                    
                    n_ont, _ = _core.region_entropy.run_region_entropy(
                        str(bedgz_input), str(tfbs_path_panel), str(out_raw_ont), gc_str, 
                        None, not verbose  # No target split for panel-specific
                    )
                    
                    process_region_entropy(out_raw_ont, out_final_ont, entropy_pon_parquet, "tfbs_baseline_ontarget")
                    out_raw_ont.unlink(missing_ok=True)
                    logger.info(f"✅ TFBS on-target: {out_final_ont} ({n_ont} TFs)")
        
        # ATAC
        if atac:
            # Always use genome-wide regions for primary output (WGS-comparable)
            atac_path_gw = atac_regions if atac_regions else (assets.atac_regions if assets.atac_available else None)
            
            if atac_path_gw and Path(atac_path_gw).exists():
                logger.info(f"Computing ATAC entropy (genome-wide)...")
                out_raw = output / f"{sample_name}.ATAC.raw.tsv"
                out_final = output / f"{sample_name}.ATAC.tsv"
                
                # Call Rust with optional target_regions for on/off-target split
                n_off, n_on = _core.region_entropy.run_region_entropy(
                    str(bedgz_input), str(atac_path_gw), str(out_raw), gc_str,
                    target_regions_str, not verbose
                )
                
                # Process off-target output (primary)
                process_region_entropy(out_raw, out_final, entropy_pon_parquet, "atac_baseline")
                out_raw.unlink(missing_ok=True)
                logger.info(f"✅ ATAC: {out_final} ({n_off} cancer types)")
                
            else:
                logger.warning("ATAC regions not available for this genome")
            
            # Panel-specific ATAC (using panel-intersection regions)
            # Output renamed to .ontarget.tsv for consistency with other features
            if is_panel_mode and assay:
                atac_path_panel = assets.get_atac_regions(assay)
                if atac_path_panel and atac_path_panel.exists() and atac_path_panel != atac_path_gw:
                    logger.info(f"Computing ATAC entropy (panel-specific: {assay})...")
                    out_raw_ont = output / f"{sample_name}.ATAC.ontarget.raw.tsv"
                    out_final_ont = output / f"{sample_name}.ATAC.ontarget.tsv"
                    
                    n_ont, _ = _core.region_entropy.run_region_entropy(
                        str(bedgz_input), str(atac_path_panel), str(out_raw_ont), gc_str, 
                        None, not verbose  # No target split for panel-specific
                    )
                    
                    process_region_entropy(out_raw_ont, out_final_ont, entropy_pon_parquet, "atac_baseline_ontarget")
                    out_raw_ont.unlink(missing_ok=True)
                    logger.info(f"✅ ATAC on-target: {out_final_ont} ({n_ont} cancer types)")
                
    except FileNotFoundError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    except RuntimeError as e:
        logger.error(f"Region entropy calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
