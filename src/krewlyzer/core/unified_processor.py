"""
Unified Feature Processor for krewlyzer.

This module provides a single entry point for running FSC, FSR, FSD, WPS, and OCF
feature extraction on BED.gz fragment files. It consolidates duplicated code 
across standalone CLI tools and wrapper.py into a single source of truth.

Architecture:
    CLI tools (fsc.py, fsd.py, etc.) → run_features() → Rust backend
    wrapper.py (run-all)             → run_features() → Rust backend

The Rust backend `run_unified_pipeline` processes all enabled features in a
single pass over the BED.gz file for maximum efficiency.

Usage:
    from krewlyzer.core.unified_processor import run_features
    
    # Single feature
    outputs = run_features(bed_path, output_dir, sample, enable_fsc=True)
    
    # All features (like run-all)
    outputs = run_features(
        bed_path, output_dir, sample,
        enable_fsc=True, enable_fsr=True, enable_fsd=True,
        enable_wps=True, enable_ocf=True
    )
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import logging
import time

from ..assets import AssetManager
from .gc_assets import resolve_gc_assets
from .fsc_processor import process_fsc
from .fsr_processor import process_fsr
from .fsd_processor import process_fsd
from .wps_processor import post_process_wps
from .pon_integration import load_pon_model
from .gene_bed import load_gene_bed

# Rust backend
from krewlyzer import _core

# Import shared resolver functions
from .utils import resolve_int as _resolve_int
from .utils import resolve_path as _resolve_path
from .utils import resolve_bool as _resolve_bool
from .utils import resolve_str as _resolve_str

logger = logging.getLogger("krewlyzer.core.unified_processor")



@dataclass
class FeatureOutputs:
    """Container for all feature output paths."""
    fsc_counts: Optional[Path] = None       # Raw FSC bin counts  
    fsc: Optional[Path] = None              # Processed FSC
    fsc_ontarget: Optional[Path] = None     # On-target FSC (panel)
    fsc_gene: Optional[Path] = None         # Gene-centric FSC
    fsr: Optional[Path] = None              # FSR ratios
    fsr_ontarget: Optional[Path] = None     # On-target FSR (panel)
    fsd: Optional[Path] = None              # FSD per arm
    fsd_ontarget: Optional[Path] = None     # On-target FSD (panel)
    wps: Optional[Path] = None              # WPS parquet
    wps_background: Optional[Path] = None   # WPS Alu stacking
    wps_panel: Optional[Path] = None        # Panel WPS (assay-specific)
    ocf: Optional[Path] = None              # OCF scores
    ocf_sync: Optional[Path] = None         # OCF sync data
    ocf_ontarget: Optional[Path] = None     # On-target OCF (panel)
    ocf_panel: Optional[Path] = None        # Panel-filtered OCF (targets + promoters)
    gc_factors: Optional[Path] = None       # GC correction factors


def run_features(
    bed_path: Path,
    output_dir: Path,
    sample_name: str,
    *,
    # Genome and asset resolution
    genome: str = "hg19",
    
    # Feature toggles
    enable_fsc: bool = False,
    enable_fsr: bool = False,
    enable_fsd: bool = False,
    enable_wps: bool = False,
    enable_ocf: bool = False,
    
    # Panel mode
    target_regions: Optional[Path] = None,
    assay: Optional[str] = None,
    
    # PON normalization
    pon_model: Optional[Path] = None,
    
    # FSC/FSR-specific options
    fsc_bins: Optional[Path] = None,
    fsc_windows: int = 100000,
    fsc_continue_n: int = 50,
    
    # FSD-specific options
    fsd_arms: Optional[Path] = None,
    
    # WPS-specific options
    wps_anchors: Optional[Path] = None,
    wps_tsv: Optional[Path] = None,          # Legacy transcript file
    wps_background: Optional[Path] = None,   # Alu BED for stacking
    wps_bait_padding: int = 50,
    wps_empty: bool = False,
    
    # OCF-specific options
    ocf_regions: Optional[Path] = None,
    
    # Shared options
    gc_correct: bool = True,
    threads: int = 0,
    verbose: bool = False,
) -> FeatureOutputs:
    """
    Run unified feature extraction pipeline.
    
    This is the single entry point for FSC/FSR/FSD/WPS/OCF feature extraction.
    Features are extracted in a single pass over the BED.gz file using the
    Rust `run_unified_pipeline` backend.
    
    Args:
        bed_path: Input .bed.gz file (from extract step)
        output_dir: Output directory for feature files
        sample_name: Sample identifier for output filenames
        genome: Genome build (hg19/GRCh37/hg38/GRCh38)
        enable_fsc: Enable Fragment Size Coverage
        enable_fsr: Enable Fragment Size Ratio (uses FSC bins)
        enable_fsd: Enable Fragment Size Distribution per arm
        enable_wps: Enable Windowed Protection Score
        enable_ocf: Enable Orientation-aware cfDNA Fragmentation
        target_regions: BED file for panel on/off-target split
        assay: Assay code (xs1/xs2) for gene-centric and panel features
        pon_model: Path to PON parquet for z-score normalization
        fsc_bins: Custom bin file for FSC/FSR (default: bundled 100kb bins)
        fsc_windows: Window aggregation size (default: 100000)
        fsc_continue_n: Bins per window (default: 50)
        fsd_arms: Custom arms file for FSD (default: bundled)
        wps_anchors: WPS anchors BED (merged TSS+CTCF)
        wps_tsv: Legacy transcript TSV for WPS
        wps_background: Alu BED for WPS background stacking
        wps_bait_padding: Bait edge padding for panel WPS (default: 50)
        wps_empty: Include regions with no coverage
        ocf_regions: Custom OCF regions file (default: bundled)
        gc_correct: Apply GC bias correction
        threads: Number of threads (0=all cores)
        verbose: Enable debug logging
        
    Returns:
        FeatureOutputs dataclass with paths to all generated files
        
    Raises:
        FileNotFoundError: If input file or required assets not found
        RuntimeError: If Rust pipeline fails
    """
    start_time = time.time()
    outputs = FeatureOutputs()
    
    # =========================================================================
    # 1. SETUP AND VALIDATION
    # =========================================================================
    
    # Configure logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("krewlyzer").setLevel(logging.DEBUG)
    
    logger.info(f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    logger.info(f"Unified Feature Processor: {sample_name}")
    logger.info(f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    
    # Input validation
    if not bed_path.exists():
        raise FileNotFoundError(f"Input BED.gz not found: {bed_path}")
    
    if not str(bed_path).endswith('.bed.gz'):
        logger.warning(f"Input may not be .bed.gz format: {bed_path}")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Resolve optional int parameters (handles typer.OptionInfo when called directly)
    resolved_threads = _resolve_int(threads, 0)
    resolved_fsc_windows = _resolve_int(fsc_windows, 100000)
    resolved_fsc_continue_n = _resolve_int(fsc_continue_n, 50)
    resolved_wps_bait_padding = _resolve_int(wps_bait_padding, 50)
    
    # Resolve optional Path parameters (handles typer.OptionInfo when called directly)
    resolved_target_regions = _resolve_path(target_regions)
    resolved_pon_model = _resolve_path(pon_model)
    resolved_fsc_bins = _resolve_path(fsc_bins)
    resolved_fsd_arms = _resolve_path(fsd_arms)
    resolved_wps_anchors = _resolve_path(wps_anchors)
    resolved_wps_tsv = _resolve_path(wps_tsv)
    resolved_wps_background = _resolve_path(wps_background)
    resolved_ocf_regions = _resolve_path(ocf_regions)
    
    # Resolve optional bool and string parameters
    resolved_wps_empty = _resolve_bool(wps_empty, False)
    resolved_assay = _resolve_str(assay)
    
    # Configure thread pool
    if resolved_threads > 0:
        try:
            _core.configure_threads(resolved_threads)
            logger.debug(f"Configured {resolved_threads} threads")
        except Exception as e:
            logger.warning(f"Could not configure threads: {e}")
    
    # Log enabled features
    features = []
    if enable_fsc: features.append("FSC")
    if enable_fsr: features.append("FSR")
    if enable_fsd: features.append("FSD")
    if enable_wps: features.append("WPS")
    if enable_ocf: features.append("OCF")
    logger.info(f"Features: {', '.join(features) or 'None enabled'}")
    
    # Panel mode detection
    is_panel_mode = resolved_target_regions is not None and resolved_target_regions.exists()
    if is_panel_mode:
        logger.info(f"Panel mode: on/off-target split enabled")
        logger.info(f"  Target regions: {resolved_target_regions.name}")
    
    if resolved_assay:
        logger.info(f"Assay: {resolved_assay} (enables gene-centric FSC, panel WPS)")
    
    # =========================================================================
    # 2. ASSET RESOLUTION
    # =========================================================================
    
    logger.debug("Resolving assets...")
    
    try:
        assets = AssetManager(genome)
        logger.info(f"Genome: {assets.raw_genome} → {assets.genome_dir}")
    except ValueError as e:
        raise FileNotFoundError(f"Invalid genome: {e}")
    
    # GC correction assets
    gc = resolve_gc_assets(
        assets, output_dir, sample_name, bed_path, gc_correct, genome
    )
    outputs.gc_factors = gc.factors_out
    
    # FSC/FSR bins (shared - only load if either enabled)
    need_fsc_bins = enable_fsc or enable_fsr
    if need_fsc_bins:
        res_bins = resolved_fsc_bins or assets.bins_100kb
        if not res_bins.exists():
            raise FileNotFoundError(f"Bin file not found: {res_bins}")
        logger.debug(f"FSC/FSR bins: {res_bins}")
    else:
        res_bins = None
    
    # FSD arms
    if enable_fsd:
        res_arms = resolved_fsd_arms or assets.arms
        if not res_arms.exists():
            raise FileNotFoundError(f"Arms file not found: {res_arms}")
        logger.debug(f"FSD arms: {res_arms}")
    else:
        res_arms = None
    
    # WPS anchors
    if enable_wps:
        if resolved_wps_anchors and resolved_wps_anchors.exists():
            res_wps = resolved_wps_anchors
        elif resolved_wps_tsv and resolved_wps_tsv.exists():
            res_wps = resolved_wps_tsv  # Legacy
        else:
            res_wps = assets.wps_anchors if hasattr(assets, 'wps_anchors') else assets.transcript_anno
        if not res_wps.exists():
            raise FileNotFoundError(f"WPS regions not found: {res_wps}")
        logger.debug(f"WPS regions: {res_wps}")
        
        # WPS background (Alu stacking)
        if resolved_wps_background and resolved_wps_background.exists():
            res_wps_bg = resolved_wps_background
        elif hasattr(assets, 'wps_background') and assets.wps_background.exists():
            res_wps_bg = assets.wps_background
        else:
            res_wps_bg = None
        
        if res_wps_bg:
            logger.debug(f"WPS background: {res_wps_bg}")
    else:
        res_wps = None
        res_wps_bg = None
    
    # OCF regions
    if enable_ocf:
        if resolved_ocf_regions and resolved_ocf_regions.exists():
            res_ocf = resolved_ocf_regions
        elif hasattr(assets, 'ocf_available') and assets.ocf_available:
            res_ocf = assets.ocf_regions
        else:
            logger.warning(f"OCF regions not available for {genome}, skipping OCF")
            enable_ocf = False
            res_ocf = None
        
        if res_ocf:
            logger.debug(f"OCF regions: {res_ocf}")
    else:
        res_ocf = None
    
    # =========================================================================
    # 3. DEFINE OUTPUT PATHS
    # =========================================================================
    
    # FSC/FSR outputs (share raw counts)
    if need_fsc_bins:
        out_fsc_counts = output_dir / f"{sample_name}.fsc_counts.tsv"
        outputs.fsc_counts = out_fsc_counts
        
        if enable_fsc:
            outputs.fsc = output_dir / f"{sample_name}.FSC.tsv"
        if enable_fsr:
            outputs.fsr = output_dir / f"{sample_name}.FSR.tsv"
    else:
        out_fsc_counts = None
    
    # FSD output
    if enable_fsd:
        outputs.fsd = output_dir / f"{sample_name}.FSD.tsv"
    
    # WPS outputs
    if enable_wps:
        outputs.wps = output_dir / f"{sample_name}.WPS.parquet"
        if res_wps_bg:
            outputs.wps_background = output_dir / f"{sample_name}.WPS_background.parquet"
    
    # OCF outputs (uses temp subdir due to Rust hardcoded names)
    if enable_ocf:
        ocf_tmp_dir = output_dir / f"{sample_name}_ocf_tmp"
        ocf_tmp_dir.mkdir(parents=True, exist_ok=True)
        outputs.ocf = output_dir / f"{sample_name}.OCF.tsv"
        outputs.ocf_sync = output_dir / f"{sample_name}.OCF.sync.tsv"
    else:
        ocf_tmp_dir = None
    
    # =========================================================================
    # 4. RUN RUST UNIFIED PIPELINE
    # =========================================================================
    
    logger.info("Running Rust unified pipeline...")
    rust_start = time.time()
    
    try:
        _core.run_unified_pipeline(
            str(bed_path),
            # GC Correction (compute)
            str(gc.gc_ref) if gc.gc_ref else None,
            str(gc.valid_regions) if gc.valid_regions else None,
            str(gc.factors_out) if gc.factors_out else None,
            # GC Correction (load pre-computed)
            str(gc.factors_input) if gc.factors_input else None,
            # FSC
            str(res_bins) if need_fsc_bins else None,
            str(out_fsc_counts) if need_fsc_bins else None,
            # WPS foreground
            str(res_wps) if enable_wps else None,
            str(outputs.wps) if enable_wps else None,
            # WPS background
            str(res_wps_bg) if res_wps_bg else None,
            str(outputs.wps_background) if res_wps_bg else None,
            resolved_wps_empty,
            # FSD
            str(res_arms) if enable_fsd else None,
            str(outputs.fsd) if enable_fsd else None,
            # OCF
            str(res_ocf) if enable_ocf else None,
            str(ocf_tmp_dir) if enable_ocf else None,
            # Target regions for on/off-target split
            str(resolved_target_regions) if is_panel_mode else None,
            resolved_wps_bait_padding,
            not verbose  # silent = not verbose
        )
    except Exception as e:
        raise RuntimeError(f"Rust pipeline failed: {e}")
    
    rust_duration = time.time() - rust_start
    logger.info(f"Rust pipeline complete in {rust_duration:.2f}s")
    
    # =========================================================================
    # 5. POST-PROCESSING (PON normalization, format conversion)
    # =========================================================================
    
    logger.info("Post-processing outputs...")
    
    # Load PON model if provided
    pon = None
    pon_parquet = None
    if resolved_pon_model:
        pon = load_pon_model(resolved_pon_model)
        pon_parquet = resolved_pon_model
        if pon:
            logger.info(f"PON loaded: {pon.assay} (n={pon.n_samples})")
    
    # Process FSC
    if enable_fsc and out_fsc_counts and out_fsc_counts.exists():
        import pandas as pd
        df_counts = pd.read_csv(out_fsc_counts, sep='\t')
        process_fsc(df_counts, outputs.fsc, resolved_fsc_windows, resolved_fsc_continue_n, pon=pon)
        logger.info(f"✓ FSC: {outputs.fsc.name}")
        
        # On-target FSC (panel mode)
        out_fsc_counts_on = output_dir / f"{sample_name}.fsc_counts.ontarget.tsv"
        if is_panel_mode and out_fsc_counts_on.exists():
            outputs.fsc_ontarget = output_dir / f"{sample_name}.FSC.ontarget.tsv"
            df_counts_on = pd.read_csv(out_fsc_counts_on, sep='\t')
            process_fsc(df_counts_on, outputs.fsc_ontarget, resolved_fsc_windows, resolved_fsc_continue_n, pon=pon)
            logger.info(f"✓ FSC on-target: {outputs.fsc_ontarget.name}")
    
    # Process FSR (same counts, different output)
    if enable_fsr and out_fsc_counts and out_fsc_counts.exists():
        import pandas as pd
        df_counts = pd.read_csv(out_fsc_counts, sep='\t')
        process_fsr(df_counts, outputs.fsr, resolved_fsc_windows, resolved_fsc_continue_n, pon=pon)
        logger.info(f"✓ FSR: {outputs.fsr.name}")
    
    # Gene-centric FSC (if assay provided)
    # Uses on-target GC correction factors for accurate copy number analysis
    if resolved_assay and enable_fsc:
        try:
            from .fsc_processor import aggregate_by_gene, load_correction_factors
            genome_map = {'hg19': 'GRCh37', 'grch37': 'GRCh37', 'hg38': 'GRCh38', 'grch38': 'GRCh38'}
            gene_genome = genome_map.get(genome.lower(), 'GRCh37')
            genes = load_gene_bed(assay=resolved_assay, genome=gene_genome)
            outputs.fsc_gene = output_dir / f"{sample_name}.FSC.gene.tsv"
            
            # Load GC correction factors (prefer on-target for panel mode)
            # On-target factors match the capture bias of gene regions
            gene_fsc_factors = None
            if is_panel_mode:
                ontarget_path = output_dir / f"{sample_name}.correction_factors.ontarget.tsv"
                if ontarget_path.exists():
                    gene_fsc_factors = load_correction_factors(ontarget_path)
                    if gene_fsc_factors:
                        logger.info(f"Gene FSC: Using on-target GC factors ({ontarget_path.name})")
                else:
                    logger.debug("Gene FSC: No on-target factors found, using raw counting")
            
            aggregate_by_gene(bed_path, genes, outputs.fsc_gene, pon=pon, 
                            correction_factors=gene_fsc_factors)
            logger.info(f"✓ FSC gene: {outputs.fsc_gene.name} ({len(genes)} genes)")
        except Exception as e:
            logger.warning(f"Gene FSC aggregation failed: {e}")
    
    # Process FSD
    if enable_fsd and outputs.fsd and outputs.fsd.exists():
        process_fsd(outputs.fsd, pon=pon, pon_parquet_path=pon_parquet)
        logger.info(f"✓ FSD: {outputs.fsd.name}")
        
        # On-target FSD
        out_fsd_on = output_dir / f"{sample_name}.FSD.ontarget.tsv"
        if is_panel_mode and out_fsd_on.exists():
            outputs.fsd_ontarget = out_fsd_on
            process_fsd(out_fsd_on, pon=pon, pon_parquet_path=pon_parquet)
            logger.info(f"✓ FSD on-target: {outputs.fsd_ontarget.name}")
    
    # Process WPS
    if enable_wps and outputs.wps and outputs.wps.exists():
        post_process_wps(
            outputs.wps,
            outputs.wps_background,
            pon_parquet_path=pon_parquet
        )
        logger.info(f"✓ WPS: {outputs.wps.name}")
    
    # Panel WPS (assay-specific)
    if resolved_assay and enable_wps:
        try:
            panel_anchors = assets.get_wps_anchors(resolved_assay)
            if panel_anchors and panel_anchors.exists():
                outputs.wps_panel = output_dir / f"{sample_name}.WPS.panel.parquet"
                logger.info(f"Running panel WPS ({resolved_assay})...")
                _core.run_unified_pipeline(
                    str(bed_path),
                    None, None, None,  # GC already computed
                    str(gc.factors_out) if gc.factors_out and gc.factors_out.exists() else str(gc.factors_input) if gc.factors_input else None,
                    None, None,  # No FSC
                    str(panel_anchors), str(outputs.wps_panel),
                    None, None, False,  # No background
                    None, None,  # No FSD
                    None, None,  # No OCF
                    str(resolved_target_regions) if is_panel_mode else None,
                    resolved_wps_bait_padding,
                    True  # silent
                )
                logger.info(f"✓ WPS panel: {outputs.wps_panel.name}")
        except Exception as e:
            logger.debug(f"Panel WPS not available: {e}")
    
    # Move OCF files from temp dir to final locations
    if enable_ocf and ocf_tmp_dir:
        import shutil
        rust_ocf = ocf_tmp_dir / "all.ocf.tsv"
        rust_sync = ocf_tmp_dir / "all.sync.tsv"
        
        if rust_ocf.exists():
            shutil.move(str(rust_ocf), str(outputs.ocf))
            logger.info(f"✓ OCF: {outputs.ocf.name}")
        if rust_sync.exists():
            shutil.move(str(rust_sync), str(outputs.ocf_sync))
        
        # On-target OCF
        if is_panel_mode:
            rust_ocf_on = ocf_tmp_dir / "all.ocf.ontarget.tsv"
            if rust_ocf_on.exists():
                outputs.ocf_ontarget = output_dir / f"{sample_name}.OCF.ontarget.tsv"
                shutil.move(str(rust_ocf_on), str(outputs.ocf_ontarget))
                logger.info(f"✓ OCF on-target: {outputs.ocf_ontarget.name}")
            
            # Panel-filtered OCF: intersect genome-wide OCF with panel targets
            # This creates a focused view for downstream analysis
            if outputs.ocf and outputs.ocf.exists() and resolved_target_regions:
                from .ocf_processor import filter_ocf_to_panel
                outputs.ocf_panel = output_dir / f"{sample_name}.OCF.panel.tsv"
                n_panel = filter_ocf_to_panel(
                    genome_ocf_path=outputs.ocf,
                    target_regions=resolved_target_regions,
                    output_path=outputs.ocf_panel,
                    promoter_extension=2000,  # 2kb upstream for promoter capture
                )
                if n_panel > 0:
                    logger.info(f"✓ OCF panel: {outputs.ocf_panel.name} ({n_panel} regions)")
        
        # Cleanup temp dir
        try:
            ocf_tmp_dir.rmdir()
        except OSError:
            pass
    
    # =========================================================================
    # 6. SUMMARY
    # =========================================================================
    
    total_duration = time.time() - start_time
    
    output_count = sum(1 for p in [
        outputs.fsc, outputs.fsr, outputs.fsd, outputs.wps, outputs.ocf,
        outputs.fsc_ontarget, outputs.fsd_ontarget, outputs.wps_panel, outputs.ocf_ontarget
    ] if p and p.exists())
    
    logger.info(f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    logger.info(f"Complete: {output_count} files generated in {total_duration:.2f}s")
    logger.info(f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    
    return outputs
