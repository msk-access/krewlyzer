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
    fsc_region: Optional[Path] = None        # Per-region FSC (exon/probe level)
    fsr: Optional[Path] = None              # FSR ratios
    fsr_ontarget: Optional[Path] = None     # On-target FSR (panel)
    fsd: Optional[Path] = None              # FSD per arm
    fsd_ontarget: Optional[Path] = None     # On-target FSD (panel)
    wps: Optional[Path] = None              # WPS parquet
    wps_background: Optional[Path] = None   # WPS Alu stacking
    wps_panel: Optional[Path] = None        # Panel WPS (assay-specific)
    ocf: Optional[Path] = None              # OCF scores (all fragments - WGS-comparable)
    ocf_sync: Optional[Path] = None         # OCF sync data (all fragments)
    ocf_ontarget: Optional[Path] = None     # On-target OCF (panel mode)
    ocf_ontarget_sync: Optional[Path] = None  # On-target OCF sync data
    ocf_offtarget: Optional[Path] = None    # Off-target OCF (panel mode)
    ocf_sync_offtarget: Optional[Path] = None  # Off-target OCF sync data
    gc_factors: Optional[Path] = None       # GC correction factors
    # Region Entropy (TFBS/ATAC size entropy)
    tfbs: Optional[Path] = None             # TFBS entropy (all fragments - genome-wide)
    tfbs_ontarget: Optional[Path] = None    # TFBS entropy (panel-specific regions)
    tfbs_sync: Optional[Path] = None        # TFBS sync (all fragments)
    tfbs_sync_ontarget: Optional[Path] = None  # TFBS sync (panel-specific regions)
    atac: Optional[Path] = None             # ATAC entropy (all fragments - genome-wide)
    atac_ontarget: Optional[Path] = None    # ATAC entropy (panel-specific regions)
    atac_sync: Optional[Path] = None        # ATAC sync (all fragments)
    atac_sync_ontarget: Optional[Path] = None  # ATAC sync (panel-specific regions)


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
    enable_tfbs: bool = False,
    enable_atac: bool = False,
    
    # Panel mode
    target_regions: Optional[Path] = None,
    assay: Optional[str] = None,
    
    # PON normalization
    pon_model: Optional[Path] = None,
    skip_pon_zscore: bool = False,  # Output raw features without z-score normalization
    
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
            res_wps = assets.wps_anchors
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
    
    # OCF regions - resolve from user input or bundled assets
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
    # Triple-output: Rust produces all/on-target/off-target with sync files
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
    pon_for_zscore = None  # Separate variable for z-score normalization
    if resolved_pon_model:
        pon = load_pon_model(resolved_pon_model)
        pon_parquet = resolved_pon_model
        if pon:
            logger.info(f"PON loaded: {pon.assay} (n={pon.n_samples})")
            # When skip_pon_zscore is True (--skip-pon mode), don't pass PON for z-scores
            if not skip_pon_zscore:
                pon_for_zscore = pon
            else:
                logger.info("  ⚠️ --skip-pon: skipping PON z-score normalization for all tools")
    
    # Process FSC
    if enable_fsc and out_fsc_counts and out_fsc_counts.exists():
        import pandas as pd
        df_counts = pd.read_csv(out_fsc_counts, sep='\t')
        if skip_pon_zscore and pon:
            logger.info("  FSC: --skip-pon active, outputting raw values (no z-scores)")
        process_fsc(df_counts, outputs.fsc, resolved_fsc_windows, resolved_fsc_continue_n, pon=pon_for_zscore)
        logger.info(f"✓ FSC: {outputs.fsc.name}")
        
        # On-target FSC (panel mode)
        out_fsc_counts_on = output_dir / f"{sample_name}.fsc_counts.ontarget.tsv"
        if is_panel_mode and out_fsc_counts_on.exists():
            outputs.fsc_ontarget = output_dir / f"{sample_name}.FSC.ontarget.tsv"
            df_counts_on = pd.read_csv(out_fsc_counts_on, sep='\t')
            process_fsc(df_counts_on, outputs.fsc_ontarget, resolved_fsc_windows, resolved_fsc_continue_n, pon=pon_for_zscore)
            logger.info(f"✓ FSC on-target: {outputs.fsc_ontarget.name}")
    
    # Process FSR (same counts, different output)
    if enable_fsr and out_fsc_counts and out_fsc_counts.exists():
        import pandas as pd
        df_counts = pd.read_csv(out_fsc_counts, sep='\t')
        if skip_pon_zscore and pon:
            logger.info("  FSR: --skip-pon active, outputting raw values (no z-scores)")
        process_fsr(df_counts, outputs.fsr, resolved_fsc_windows, resolved_fsc_continue_n, pon=pon_for_zscore)
        logger.info(f"✓ FSR: {outputs.fsr.name}")
    
    # Gene-centric FSC (if assay provided)
    # Uses on-target GC correction factors for accurate copy number analysis
    if resolved_assay and enable_fsc:
        try:
            from .fsc_processor import aggregate_by_gene, load_correction_factors, apply_fsc_gene_pon, apply_fsc_region_pon
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
            
            aggregate_by_gene(bed_path, genes, outputs.fsc_gene, pon=pon_for_zscore, 
                            correction_factors=gene_fsc_factors, aggregate_by='gene')
            
            # Apply FSC gene PON z-score normalization (unless --skip-pon)
            if pon and hasattr(pon, 'fsc_gene_baseline') and pon.fsc_gene_baseline and not skip_pon_zscore:
                apply_fsc_gene_pon(outputs.fsc_gene, pon)
                logger.info(f"✓ FSC gene: {outputs.fsc_gene.name} ({len(genes)} genes, with PON z-scores)")
            elif skip_pon_zscore and pon:
                logger.info(f"✓ FSC gene: {outputs.fsc_gene.name} ({len(genes)} genes, --skip-pon active)")
            else:
                logger.info(f"✓ FSC gene: {outputs.fsc_gene.name} ({len(genes)} genes)")
            
            # Also generate per-region FSC (exon/probe level)
            outputs.fsc_region = output_dir / f"{sample_name}.FSC.regions.tsv"
            aggregate_by_gene(bed_path, genes, outputs.fsc_region, pon=pon_for_zscore, 
                            correction_factors=gene_fsc_factors, aggregate_by='region')
            
            # Apply FSC region PON z-score normalization (unless --skip-pon)
            if pon and hasattr(pon, 'fsc_region_baseline') and pon.fsc_region_baseline and not skip_pon_zscore:
                apply_fsc_region_pon(outputs.fsc_region, pon)
                logger.info(f"✓ FSC regions: {outputs.fsc_region.name} (with PON z-scores)")
            elif skip_pon_zscore and pon:
                logger.info(f"✓ FSC regions: {outputs.fsc_region.name} (--skip-pon active)")
            else:
                logger.info(f"✓ FSC regions: {outputs.fsc_region.name}")
        except Exception as e:
            logger.warning(f"Gene FSC aggregation failed: {e}")
    
    # Process FSD
    if enable_fsd and outputs.fsd and outputs.fsd.exists():
        if skip_pon_zscore and pon:
            logger.info("  FSD: --skip-pon active, outputting raw values (no z-scores)")
        process_fsd(outputs.fsd, pon_parquet_path=pon_parquet)
        logger.info(f"✓ FSD: {outputs.fsd.name}")
        
        # On-target FSD
        out_fsd_on = output_dir / f"{sample_name}.FSD.ontarget.tsv"
        if is_panel_mode and out_fsd_on.exists():
            outputs.fsd_ontarget = out_fsd_on
            process_fsd(out_fsd_on, pon_parquet_path=pon_parquet)
            logger.info(f"✓ FSD on-target: {outputs.fsd_ontarget.name}")
    
    # Process WPS
    if enable_wps and outputs.wps and outputs.wps.exists():
        if skip_pon_zscore and pon:
            logger.info("  WPS: --skip-pon active, outputting raw values (no z-scores)")
        post_process_wps(
            outputs.wps,
            outputs.wps_background,
            pon=pon_for_zscore
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
    # Rust OCF now produces 6 files in panel mode:
    # - all.ocf.tsv, all.sync.tsv (ALL fragments)
    # - all.ocf.ontarget.tsv, all.sync.ontarget.tsv (on-target)
    # - all.ocf.offtarget.tsv, all.sync.offtarget.tsv (off-target)
    if enable_ocf and ocf_tmp_dir:
        import shutil
        
        # Primary output (ALL fragments - WGS-comparable)
        rust_ocf = ocf_tmp_dir / "all.ocf.tsv"
        rust_sync = ocf_tmp_dir / "all.sync.tsv"
        
        if rust_ocf.exists():
            shutil.move(str(rust_ocf), str(outputs.ocf))
            logger.info(f"✓ OCF: {outputs.ocf.name}")
            
            # Apply OCF PON normalization (z-scores) if PON is available and not skipped
            if pon_parquet and not skip_pon_zscore:
                from .ocf_processor import process_ocf_with_pon
                process_ocf_with_pon(outputs.ocf, pon_parquet)
                logger.info(f"  OCF PON z-scores applied ({outputs.ocf.name})")
            elif skip_pon_zscore and pon_parquet:
                logger.info("  OCF: --skip-pon active, outputting raw values (no z-scores)")
            
        if rust_sync.exists():
            shutil.move(str(rust_sync), str(outputs.ocf_sync))
        
        # On-target output (panel mode only)
        rust_ocf_ontarget = ocf_tmp_dir / "all.ocf.ontarget.tsv"
        rust_sync_ontarget = ocf_tmp_dir / "all.sync.ontarget.tsv"
        
        if rust_ocf_ontarget.exists():
            outputs.ocf_ontarget = output_dir / f"{sample_name}.OCF.ontarget.tsv"
            shutil.move(str(rust_ocf_ontarget), str(outputs.ocf_ontarget))
            logger.info(f"✓ OCF on-target: {outputs.ocf_ontarget.name}")
        if rust_sync_ontarget.exists():
            outputs.ocf_ontarget_sync = output_dir / f"{sample_name}.OCF.ontarget.sync.tsv"
            shutil.move(str(rust_sync_ontarget), str(outputs.ocf_ontarget_sync))
        
        # Off-target output (panel mode only)
        rust_ocf_offtarget = ocf_tmp_dir / "all.ocf.offtarget.tsv"
        rust_sync_offtarget = ocf_tmp_dir / "all.sync.offtarget.tsv"
        
        if rust_ocf_offtarget.exists():
            outputs.ocf_offtarget = output_dir / f"{sample_name}.OCF.offtarget.tsv"
            shutil.move(str(rust_ocf_offtarget), str(outputs.ocf_offtarget))
            logger.info(f"✓ OCF off-target: {outputs.ocf_offtarget.name}")
        if rust_sync_offtarget.exists():
            outputs.ocf_sync_offtarget = output_dir / f"{sample_name}.OCF.offtarget.sync.tsv"
            shutil.move(str(rust_sync_offtarget), str(outputs.ocf_sync_offtarget))
        
        # Cleanup temp dir
        try:
            for f in ocf_tmp_dir.glob("*"):
                f.unlink()
            ocf_tmp_dir.rmdir()
        except OSError:
            pass
        
        # Note: The panel OCF run with filtered regions is no longer needed
        # because the Rust triple-output now handles on/off-target splitting
        # automatically when target_regions is provided
    
    # =========================================================================
    # 6. TFBS/ATAC REGION ENTROPY
    # =========================================================================
    
    if enable_tfbs or enable_atac:
        from .region_entropy_processor import process_region_entropy
        # Note: load_pon_model is imported at module level (line 41)
        
        # =====================================================================
        # Load GC correction factors for region entropy
        # Off-target: WGS-like GC bias model (used for global/off-target fragments)
        # On-target: capture-specific GC bias model (used for on-target fragments)
        # =====================================================================
        gc_factors_path = output_dir / f"{sample_name}.correction_factors.tsv"
        gc_factors_ontarget_path = output_dir / f"{sample_name}.correction_factors.ontarget.tsv"
        
        gc_str = str(gc_factors_path) if gc_factors_path.exists() else None
        gc_ontarget_str = str(gc_factors_ontarget_path) if gc_factors_ontarget_path.exists() else None
        
        # Log GC factor availability for debugging
        if gc_str:
            logger.debug(f"  GC factors (off-target): {gc_factors_path.name}")
        if is_panel_mode and gc_ontarget_str:
            logger.debug(f"  GC factors (on-target): {gc_factors_ontarget_path.name}")
        elif is_panel_mode and not gc_ontarget_str:
            logger.debug("  GC factors (on-target): not available, using off-target fallback")
        
        # Use pon_parquet directly for Z-score normalization (Rust implementation)
        entropy_pon_parquet = pon_parquet if (pon_parquet and not skip_pon_zscore) else None
        
        if skip_pon_zscore and pon_parquet:
            if enable_tfbs:
                logger.info("  TFBS: --skip-pon active, outputting raw values (no z-scores)")
            if enable_atac:
                logger.info("  ATAC: --skip-pon active, outputting raw values (no z-scores)")
        
        # TFBS Entropy - Dual output architecture:
        # 1. Genome-wide: Uses all TFBS regions for WGS-comparable output
        # 2. Panel-specific: Uses pre-intersected TFBS regions that overlap with capture
        if enable_tfbs and assets.tfbs_available:
            try:
                # =========================================================
                # GENOME-WIDE TFBS (primary output - WGS-comparable)
                # =========================================================
                tfbs_regions_gw = str(assets.tfbs_regions)
                out_tfbs_raw = output_dir / f"{sample_name}.TFBS.raw.tsv"
                
                # Call Rust with genome-wide regions (no target filtering)
                n_all, _ = _core.region_entropy.run_region_entropy(
                    str(bed_path),      # Fragment BED.gz
                    tfbs_regions_gw,    # TFBS region BED (genome-wide)
                    str(out_tfbs_raw),  # Output path
                    gc_str,             # GC factors
                    None,               # No on-target GC (not splitting)
                    None,               # No target regions (no split)
                    True                # Silent mode
                )
                
                # Process output with z-score normalization
                outputs.tfbs = output_dir / f"{sample_name}.TFBS.tsv"
                process_region_entropy(out_tfbs_raw, outputs.tfbs, entropy_pon_parquet, "tfbs_baseline")
                out_tfbs_raw.unlink(missing_ok=True)
                logger.info(f"✓ TFBS entropy: {outputs.tfbs.name} ({n_all} TFs)")
                
                # Keep sync file as-is (no z-score normalization needed)
                out_tfbs_sync_raw = output_dir / f"{sample_name}.TFBS.sync.raw.tsv"
                if out_tfbs_sync_raw.exists():
                    outputs.tfbs_sync = output_dir / f"{sample_name}.TFBS.sync.tsv"
                    out_tfbs_sync_raw.rename(outputs.tfbs_sync)
                
                # =========================================================
                # PANEL-SPECIFIC TFBS (uses pre-intersected regions)
                # Only generated when assay is specified and panel regions exist
                # =========================================================
                if is_panel_mode and resolved_assay:
                    tfbs_regions_panel = assets.get_tfbs_regions(resolved_assay)
                    if tfbs_regions_panel and tfbs_regions_panel.exists() and tfbs_regions_panel != assets.tfbs_regions:
                        logger.info(f"Computing TFBS entropy (panel-specific: {resolved_assay})...")
                        out_tfbs_ont_raw = output_dir / f"{sample_name}.TFBS.ontarget.raw.tsv"
                        
                        # Call Rust with panel-specific regions
                        n_ont, _ = _core.region_entropy.run_region_entropy(
                            str(bed_path),              # Fragment BED.gz
                            str(tfbs_regions_panel),    # TFBS regions (panel-specific)
                            str(out_tfbs_ont_raw),      # Output path
                            gc_ontarget_str or gc_str,  # Use on-target GC if available
                            None,                       # No split
                            None,                       # No target filtering
                            True                        # Silent mode
                        )
                        
                        # Process output with ontarget baseline
                        outputs.tfbs_ontarget = output_dir / f"{sample_name}.TFBS.ontarget.tsv"
                        process_region_entropy(out_tfbs_ont_raw, outputs.tfbs_ontarget, 
                                             entropy_pon_parquet, "tfbs_baseline_ontarget")
                        out_tfbs_ont_raw.unlink(missing_ok=True)
                        logger.info(f"✓ TFBS on-target: {outputs.tfbs_ontarget.name} ({n_ont} TFs)")
                        
                        # Keep sync file for panel-specific
                        out_tfbs_ont_sync_raw = output_dir / f"{sample_name}.TFBS.ontarget.sync.raw.tsv"
                        if out_tfbs_ont_sync_raw.exists():
                            outputs.tfbs_sync_ontarget = output_dir / f"{sample_name}.TFBS.ontarget.sync.tsv"
                            out_tfbs_ont_sync_raw.rename(outputs.tfbs_sync_ontarget)
                    
            except Exception as e:
                logger.warning(f"TFBS entropy failed: {e}")
        
        # ATAC Entropy - Dual output architecture:
        # 1. Genome-wide: Uses all ATAC regions for WGS-comparable output
        # 2. Panel-specific: Uses pre-intersected ATAC regions that overlap with capture
        if enable_atac and assets.atac_available:
            try:
                # =========================================================
                # GENOME-WIDE ATAC (primary output - WGS-comparable)
                # =========================================================
                atac_regions_gw = str(assets.atac_regions)
                out_atac_raw = output_dir / f"{sample_name}.ATAC.raw.tsv"
                
                # Call Rust with genome-wide regions (no target filtering)
                n_all, _ = _core.region_entropy.run_region_entropy(
                    str(bed_path),      # Fragment BED.gz
                    atac_regions_gw,    # ATAC region BED (genome-wide)
                    str(out_atac_raw),  # Output path
                    gc_str,             # GC factors
                    None,               # No on-target GC (not splitting)
                    None,               # No target regions (no split)
                    True                # Silent mode
                )
                
                # Process output with z-score normalization
                outputs.atac = output_dir / f"{sample_name}.ATAC.tsv"
                process_region_entropy(out_atac_raw, outputs.atac, entropy_pon_parquet, "atac_baseline")
                out_atac_raw.unlink(missing_ok=True)
                logger.info(f"✓ ATAC entropy: {outputs.atac.name} ({n_all} cancer types)")
                
                # Keep sync file as-is (no z-score normalization needed)
                out_atac_sync_raw = output_dir / f"{sample_name}.ATAC.sync.raw.tsv"
                if out_atac_sync_raw.exists():
                    outputs.atac_sync = output_dir / f"{sample_name}.ATAC.sync.tsv"
                    out_atac_sync_raw.rename(outputs.atac_sync)
                
                # =========================================================
                # PANEL-SPECIFIC ATAC (uses pre-intersected regions)
                # Only generated when assay is specified and panel regions exist
                # =========================================================
                if is_panel_mode and resolved_assay:
                    atac_regions_panel = assets.get_atac_regions(resolved_assay)
                    if atac_regions_panel and atac_regions_panel.exists() and atac_regions_panel != assets.atac_regions:
                        logger.info(f"Computing ATAC entropy (panel-specific: {resolved_assay})...")
                        out_atac_ont_raw = output_dir / f"{sample_name}.ATAC.ontarget.raw.tsv"
                        
                        # Call Rust with panel-specific regions
                        n_ont, _ = _core.region_entropy.run_region_entropy(
                            str(bed_path),              # Fragment BED.gz
                            str(atac_regions_panel),    # ATAC regions (panel-specific)
                            str(out_atac_ont_raw),      # Output path
                            gc_ontarget_str or gc_str,  # Use on-target GC if available
                            None,                       # No split
                            None,                       # No target filtering
                            True                        # Silent mode
                        )
                        
                        # Process output with ontarget baseline
                        outputs.atac_ontarget = output_dir / f"{sample_name}.ATAC.ontarget.tsv"
                        process_region_entropy(out_atac_ont_raw, outputs.atac_ontarget, 
                                             entropy_pon_parquet, "atac_baseline_ontarget")
                        out_atac_ont_raw.unlink(missing_ok=True)
                        logger.info(f"✓ ATAC on-target: {outputs.atac_ontarget.name} ({n_ont} cancer types)")
                        
                        # Keep sync file for panel-specific
                        out_atac_ont_sync_raw = output_dir / f"{sample_name}.ATAC.ontarget.sync.raw.tsv"
                        if out_atac_ont_sync_raw.exists():
                            outputs.atac_sync_ontarget = output_dir / f"{sample_name}.ATAC.ontarget.sync.tsv"
                            out_atac_ont_sync_raw.rename(outputs.atac_sync_ontarget)
                    
            except Exception as e:
                logger.warning(f"ATAC entropy failed: {e}")
    
    # =========================================================================
    # 7. SUMMARY
    # =========================================================================
    
    total_duration = time.time() - start_time
    
    output_count = sum(1 for p in [
        outputs.fsc, outputs.fsr, outputs.fsd, outputs.wps, outputs.ocf,
        outputs.fsc_ontarget, outputs.fsd_ontarget, outputs.wps_panel, outputs.ocf_ontarget,
        outputs.tfbs, outputs.atac, outputs.tfbs_ontarget, outputs.atac_ontarget
    ] if p and p.exists())
    
    logger.info(f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    logger.info(f"Complete: {output_count} files generated in {total_duration:.2f}s")
    logger.info(f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    
    return outputs
