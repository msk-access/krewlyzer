"""
PON model building from healthy plasma samples.

This module provides the CLI and logic for building unified PON models.
"""

import typer
from pathlib import Path
from typing import Optional, List
from datetime import datetime
import logging

import numpy as np
import pandas as pd
from rich.console import Console
from rich.logging import RichHandler

from .model import PonModel, GcBiasModel, FsdBaseline, WpsBaseline

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("build-pon")

# Import core tools for processing samples
from krewlyzer import _core

# Import unified sample processor
from krewlyzer.core.sample_processor import process_sample, SampleParams, SampleOutputs


def build_pon(
    sample_list: Path = typer.Argument(..., help="Text file with paths to BAM/CRAM or BED.gz files (one per line). BAM/CRAM required for MDS baseline."),
    assay: str = typer.Option(..., "--assay", "-a", help="Assay name (e.g., msk-access-v2)"),
    reference: Path = typer.Option(..., "--reference", "-r", help="Reference FASTA file"),
    wps_anchors: Optional[Path] = typer.Option(None, "--wps-anchors", "-W", help="WPS anchors BED.gz for WPS baseline (default: bundled TSS+CTCF anchors)"),
    transcript_file: Optional[Path] = typer.Option(None, "--transcript-file", "-t", help="[LEGACY] Transcript TSV for WPS baseline (use --wps-anchors instead)"),
    bin_file: Optional[Path] = typer.Option(None, "--bin-file", "-b", help="Bin file for FSC/FSR (default: hg19_window_100kb.bed)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output PON model file (.pon.parquet)"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="BED file with target regions (panel mode - builds dual on/off-target baselines)"),
    temp_dir: Optional[Path] = typer.Option(None, "--temp-dir", help="Directory for temporary files (default: system temp)"),
    threads: int = typer.Option(4, "--threads", "-p", help="Number of threads"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    require_proper_pair: bool = typer.Option(False, "--require-proper-pair", help="Only extract properly paired reads (default: False for v1 ACCESS compatibility)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
):
    """
    Build a unified PON model from healthy plasma samples.
    
    Input can be BAM/CRAM files or pre-extracted BED.gz files:
    - BAM/CRAM: Full extraction + MDS baseline (recommended)
    - BED.gz: Faster but no MDS baseline (requires sequence data)
    
    The PON model contains:
    - GC bias curves for FSC/FSR/WPS correction
    - FSD baseline per chromosome arm
    - WPS baseline per transcript region
    - MDS baseline (BAM input only)
    
    For panel data (--target-regions):
    - Builds DUAL baselines for on-target and off-target regions
    - GC model trained on off-target only (unbiased by capture)
    - FSD/WPS include separate on-target stats
    
    Example:
        krewlyzer build-pon bams.txt --assay msk-access-v2 -r hg19.fa -o pon.parquet
        krewlyzer build-pon samples.txt -a msk-access-v2 -r hg19.fa -T targets.bed -o pon.parquet
    """
    
    # Configure logging
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Configure threads
    if threads > 0:
        try:
            _core.configure_threads(threads)
            logger.info(f"Configured {threads} threads")
        except Exception as e:
            logger.warning(f"Could not configure threads: {e}")
    
    # Validate inputs
    if not sample_list.exists():
        logger.error(f"Sample list not found: {sample_list}")
        raise typer.Exit(1)
    
    if not reference.exists():
        logger.error(f"Reference FASTA not found: {reference}")
        raise typer.Exit(1)
    
    # Read sample list
    with open(sample_list) as f:
        samples = [Path(line.strip()) for line in f if line.strip() and not line.startswith("#")]
    
    n_samples = len(samples)
    if n_samples < 1:
        logger.error("No samples found in sample list")
        raise typer.Exit(1)
    
    # Panel mode detection
    is_panel_mode = target_regions is not None and target_regions.exists()
    target_regions_str = str(target_regions) if is_panel_mode else None
    if is_panel_mode:
        logger.info(f"PANEL MODE: Building dual on/off-target baselines")
        logger.info(f"  Target regions: {target_regions}")
        logger.info(f"  GC model will use OFF-TARGET fragments only (unbiased)")
        logger.info(f"  FSD/WPS will have separate on/off-target stats")
    else:
        logger.info("WGS MODE: Building single baseline")
    
    logger.info(f"Building PON from {n_samples} samples")
    logger.info(f"Assay: {assay}")
    logger.info(f"Reference: {reference}")
    logger.info(f"Genome: {genome}")
    
    # Use AssetManager for bundled data files
    from krewlyzer.assets import AssetManager
    assets = AssetManager(genome)
    
    # Default bin file
    if bin_file is None:
        bin_file = assets.bins_100kb
    
    if not bin_file.exists():
        logger.error(f"Bin file not found: {bin_file}")
        raise typer.Exit(1)
    
    # Default WPS anchors file (prefer wps_anchors, fallback to transcript_file for legacy)
    wps_regions = wps_anchors or transcript_file
    if wps_regions is None:
        wps_regions = assets.wps_anchors if hasattr(assets, 'wps_anchors') and assets.wps_anchors.exists() else assets.transcript_anno
    
    if not wps_regions.exists():
        logger.warning(f"WPS regions file not found: {wps_regions}")
    
    import tempfile
    import shutil
    import time
    from concurrent.futures import ProcessPoolExecutor, as_completed
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UNIFIED SAMPLE PROCESSING (replaces Phase 1 + Phase 2)
    # Uses process_sample() from core/sample_processor.py for consistent processing
    # ═══════════════════════════════════════════════════════════════════════════
    
    logger.info(f"=" * 60)
    logger.info(f"SAMPLE PROCESSING ({n_samples} samples)")
    logger.info(f"=" * 60)
    
    # Setup parameters for sample processing
    sample_params = SampleParams(
        genome=genome,
        mapq=20,
        minlen=65,
        maxlen=400,
        skip_duplicates=True,
        require_proper_pair=require_proper_pair,
        kmer=4,
        threads=max(1, threads // max(1, min(4, n_samples))),  # Divide threads among samples
        bin_file=bin_file,
        wps_anchors=wps_regions,  # Use resolved WPS anchors (modern) or transcript file (legacy)
    )
    
    # Use temp directory for sample processing outputs
    temp_base = str(temp_dir) if temp_dir else None
    temp_output_dir = tempfile.mkdtemp(prefix="pon_build_", dir=temp_base)
    logger.info(f"Temporary output directory: {temp_output_dir}")
    
    # Process samples
    all_outputs: List[SampleOutputs] = []
    processing_start = time.time()
    succeeded_count = 0
    failed_count = 0
    
    try:
        for i, sample_path in enumerate(samples, 1):
            # Extract sample name from path (handle BAM, CRAM, and BED.gz)
            sample_name = sample_path.stem
            for suffix in ['.bed', '.bam', '.cram']:
                sample_name = sample_name.replace(suffix, '')
            sample_start = time.time()
            
            logger.info(f"[{i}/{n_samples}] Processing: {sample_name}")
            
            try:
                outputs = process_sample(
                    input_path=sample_path,
                    output_dir=Path(temp_output_dir) / sample_name,
                    sample_name=sample_name,
                    reference=reference,
                    params=sample_params,
                    target_regions=target_regions if is_panel_mode else None,
                    enable_fsc=True,
                    enable_fsr=False,  # Not needed for PON
                    enable_fsd=True,
                    enable_wps=True,
                    enable_ocf=True,
                    pon_mode=True,  # Skip PON normalization (we're building the PON)
                )
                
                all_outputs.append(outputs)
                elapsed = time.time() - sample_start
                logger.info(f"  ✓ Done in {elapsed:.1f}s ({outputs.fragment_count:,} frags)")
                succeeded_count += 1
                
            except Exception as e:
                logger.warning(f"  ✗ Failed: {e}")
                failed_count += 1
        
        processing_elapsed = time.time() - processing_start
        logger.info(f"-" * 60)
        logger.info(f"Processing complete: {succeeded_count} succeeded, {failed_count} failed ({processing_elapsed:.1f}s)")
        
        if len(all_outputs) < 1:
            logger.error("No samples processed successfully")
            raise typer.Exit(1)
        
        # ═══════════════════════════════════════════════════════════════════════════
        # COLLECT DATA FOR BASELINE AGGREGATION
        # ═══════════════════════════════════════════════════════════════════════════
        
        logger.info(f"=" * 60)
        logger.info(f"COLLECTING BASELINE DATA ({len(all_outputs)} samples)")
        logger.info(f"=" * 60)
        
        all_gc_data = []
        all_fsd_data = []
        all_wps_data = []
        all_ocf_data = []
        all_mds_data = []
        all_fsd_data_ontarget = []
        all_gc_data_ontarget = []
        
        for outputs in all_outputs:
            sample_name = outputs.sample_id
            feat = outputs.feature_outputs
            
            # Collect GC data from FSC output
            if feat and feat.fsc_counts and feat.fsc_counts.exists():
                try:
                    fsc_df = pd.read_csv(feat.fsc_counts, sep='\t')
                    if all(col in fsc_df.columns for col in ['ultra_short', 'core_short', 'mono_nucl', 'di_nucl', 'long', 'mean_gc']):
                        all_gc_data.append({
                            "gc": fsc_df['mean_gc'].values,
                            "short": (fsc_df['ultra_short'] + fsc_df['core_short']).values,
                            "intermediate": fsc_df['mono_nucl'].values,
                            "long": (fsc_df['di_nucl'] + fsc_df['long']).values,
                        })
                except Exception as e:
                    logger.debug(f"  Could not read GC data for {sample_name}: {e}")
            
            # Collect FSD data
            if feat and feat.fsd and feat.fsd.exists():
                try:
                    fsd_df = pd.read_csv(feat.fsd, sep='\t')
                    all_fsd_data.append(fsd_df)
                except Exception as e:
                    logger.debug(f"  Could not read FSD for {sample_name}: {e}")
            
            # Collect on-target FSD (panel mode)
            if feat and feat.fsd_ontarget and feat.fsd_ontarget.exists():
                try:
                    fsd_on_df = pd.read_csv(feat.fsd_ontarget, sep='\t')
                    all_fsd_data_ontarget.append(fsd_on_df)
                except Exception as e:
                    logger.debug(f"  Could not read on-target FSD for {sample_name}: {e}")
            
            # Collect WPS data
            if feat and feat.wps and feat.wps.exists():
                try:
                    wps_df = pd.read_parquet(feat.wps)
                    # Aggregate by region for baseline
                    if 'gene_id' in wps_df.columns:
                        region_stats = wps_df.groupby('gene_id').agg({
                            'wps_long': 'mean',
                            'wps_short': 'mean'
                        }).reset_index()
                        all_wps_data.append(region_stats)
                    elif 'region_id' in wps_df.columns:
                        all_wps_data.append(wps_df)
                except Exception as e:
                    logger.debug(f"  Could not read WPS for {sample_name}: {e}")
            
            # Collect OCF data
            if feat and feat.ocf and feat.ocf.exists():
                try:
                    ocf_df = pd.read_csv(feat.ocf, sep='\t')
                    if 'tissue' in ocf_df.columns and 'OCF' in ocf_df.columns:
                        ocf_df = ocf_df.rename(columns={'tissue': 'region_id', 'OCF': 'ocf'})
                    if 'region_id' in ocf_df.columns and 'ocf' in ocf_df.columns:
                        all_ocf_data.append(ocf_df[['region_id', 'ocf']])
                except Exception as e:
                    logger.debug(f"  Could not read OCF for {sample_name}: {e}")
            
            # Collect MDS data (from extraction)
            if outputs.mds_counts:
                all_mds_data.append({
                    "kmers": outputs.mds_counts,
                    "mds": outputs.mds_score
                })
        
        logger.info(f"  GC samples: {len(all_gc_data)}")
        logger.info(f"  FSD samples: {len(all_fsd_data)}")
        logger.info(f"  WPS samples: {len(all_wps_data)}")
        logger.info(f"  OCF samples: {len(all_ocf_data)}")
        logger.info(f"  MDS samples: {len(all_mds_data)}")
        if is_panel_mode:
            logger.info(f"  FSD on-target samples: {len(all_fsd_data_ontarget)}")
    
    finally:
        # Cleanup temp directory
        if temp_output_dir and Path(temp_output_dir).exists():
            shutil.rmtree(temp_output_dir)
            logger.debug(f"Cleaned up temp directory: {temp_output_dir}")
    
    if len(all_gc_data) < 1:
        logger.error("No samples processed successfully")
        raise typer.Exit(1)

    
    logger.info(f"=" * 60)
    logger.info(f"PHASE 3: Model Building ({len(all_gc_data)} samples)")
    logger.info(f"=" * 60)
    
    model_start = time.time()
    
    # Build GC bias model
    logger.info("  Computing GC bias curves...")
    gc_bias = _compute_gc_bias_model(all_gc_data)
    
    # Build FSD baseline
    logger.info("  Computing FSD baseline...")
    fsd_baseline = _compute_fsd_baseline(all_fsd_data)
    
    # Build WPS baseline
    logger.info("  Computing WPS baseline...")
    wps_baseline = _compute_wps_baseline(all_wps_data)
    
    # Build OCF baseline
    ocf_baseline = None
    if all_ocf_data:
        logger.info("  Computing OCF baseline...")
        ocf_baseline = _compute_ocf_baseline(all_ocf_data)
    
    # Build MDS baseline
    mds_baseline = None
    if all_mds_data:
        logger.info("  Computing MDS baseline...")
        mds_baseline = _compute_mds_baseline(all_mds_data)
    
    # Build on-target FSD baseline (panel mode only)
    fsd_baseline_ontarget = None
    if is_panel_mode and all_fsd_data_ontarget:
        logger.info("  Computing on-target FSD baseline...")
        fsd_baseline_ontarget = _compute_fsd_baseline(all_fsd_data_ontarget)
        if fsd_baseline_ontarget:
            logger.debug(f"    On-target FSD: {len(fsd_baseline_ontarget.arms)} arms")
    
    # Build on-target GC bias (panel mode only)
    # Note: gc_bias_ontarget uses the same data as gc_bias since FSC already applies
    # target region filtering. For now, set to None - full implementation requires
    # separate on-target FSC data collection which is done in extract, not build-pon.
    gc_bias_ontarget = None
    if is_panel_mode and all_gc_data_ontarget:
        logger.info("  Computing on-target GC bias curves...")
        gc_bias_ontarget = _compute_gc_bias_model(all_gc_data_ontarget)
    
    # Create PON model
    model = PonModel(
        schema_version="1.0",
        assay=assay,
        build_date=datetime.now().isoformat()[:10],
        n_samples=len(all_gc_data),
        reference=reference.name.replace(".fa", "").replace(".fasta", ""),
        panel_mode=is_panel_mode,
        target_regions_file=target_regions.name if is_panel_mode else "",
        gc_bias=gc_bias,
        fsd_baseline=fsd_baseline,
        wps_baseline=wps_baseline,
        ocf_baseline=ocf_baseline,
        mds_baseline=mds_baseline,
        fsd_baseline_ontarget=fsd_baseline_ontarget,
        gc_bias_ontarget=gc_bias_ontarget,
    )
    
    # Validate model
    errors = model.validate()
    if errors:
        for error in errors:
            logger.warning(f"Validation warning: {error}")
    
    # Save model
    logger.info(f"  Saving PON model to {output}")
    _save_pon_model(model, output)
    
    model_elapsed = time.time() - model_start
    logger.info(f"-" * 60)
    logger.info(f"Model building complete ({model_elapsed:.1f}s)")
    logger.info(f"")
    logger.info(f"✅ PON model built successfully")
    logger.info(f"   Assay: {model.assay}")
    logger.info(f"   Samples: {model.n_samples}")
    logger.info(f"   Panel mode: {is_panel_mode}")
    logger.info(f"   Output: {output}")


def _compute_gc_bias_model(all_gc_data: List[dict]) -> GcBiasModel:
    """
    Compute GC bias curves from sample data.
    
    Uses high-performance Rust implementation (3-10x faster) with Python fallback.
    For each GC bin, compute median expected coverage across samples.
    """
    # Try Rust implementation first (3-10x faster)
    try:
        from krewlyzer import _core
        result = _core.pon_builder.compute_gc_bias_model(all_gc_data)
        if result and "gc_bins" in result:
            logger.info(f"GC bias model computed (Rust): {len(result['gc_bins'])} bins")
            return GcBiasModel(
                short_expected=result["short_expected"],
                short_std=result["short_std"],
                intermediate_expected=result["intermediate_expected"],
                intermediate_std=result["intermediate_std"],
                long_expected=result["long_expected"],
                long_std=result["long_std"],
            )
    except Exception as e:
        logger.debug(f"Rust GC bias computation failed, falling back to Python: {e}")
    
    # =========================================================================
    # PYTHON FALLBACK IMPLEMENTATION
    # NOTE: This is a fallback only. Primary implementation is Rust (3-10x faster).
    # The Rust implementation is in: rust/src/pon_builder.rs::compute_gc_bias_model()
    # =========================================================================
    
    # Define GC bins
    gc_bins = np.arange(0.25, 0.75, 0.02).tolist()
    
    # Aggregate counts per GC bin across all samples
    short_by_gc = {gc: [] for gc in gc_bins}
    intermediate_by_gc = {gc: [] for gc in gc_bins}
    long_by_gc = {gc: [] for gc in gc_bins}
    
    for sample in all_gc_data:
        gc = sample["gc"]
        short = sample["short"]
        intermediate = sample["intermediate"]
        long = sample["long"]
        
        # Normalize to relative coverage (mean = 1.0)
        short_norm = short / np.nanmean(short) if np.nanmean(short) > 0 else short
        intermediate_norm = intermediate / np.nanmean(intermediate) if np.nanmean(intermediate) > 0 else intermediate
        long_norm = long / np.nanmean(long) if np.nanmean(long) > 0 else long
        
        # Bin by GC
        for i, gc_val in enumerate(gc):
            if np.isnan(gc_val) or gc_val <= 0:
                continue
            # Find nearest bin
            bin_idx = int((gc_val - 0.25) / 0.02)
            if 0 <= bin_idx < len(gc_bins):
                gc_bin = gc_bins[bin_idx]
                if not np.isnan(short_norm[i]):
                    short_by_gc[gc_bin].append(short_norm[i])
                if not np.isnan(intermediate_norm[i]):
                    intermediate_by_gc[gc_bin].append(intermediate_norm[i])
                if not np.isnan(long_norm[i]):
                    long_by_gc[gc_bin].append(long_norm[i])
    
    # Compute median and std per bin
    short_expected = []
    short_std = []
    intermediate_expected = []
    intermediate_std = []
    long_expected = []
    long_std = []
    
    for gc_bin in gc_bins:
        if short_by_gc[gc_bin]:
            short_expected.append(float(np.median(short_by_gc[gc_bin])))
            short_std.append(float(np.std(short_by_gc[gc_bin])))
        else:
            short_expected.append(1.0)
            short_std.append(0.0)
        
        if intermediate_by_gc[gc_bin]:
            intermediate_expected.append(float(np.median(intermediate_by_gc[gc_bin])))
            intermediate_std.append(float(np.std(intermediate_by_gc[gc_bin])))
        else:
            intermediate_expected.append(1.0)
            intermediate_std.append(0.0)
        
        if long_by_gc[gc_bin]:
            long_expected.append(float(np.median(long_by_gc[gc_bin])))
            long_std.append(float(np.std(long_by_gc[gc_bin])))
        else:
            long_expected.append(1.0)
            long_std.append(0.0)
    
    return GcBiasModel(
        gc_bins=gc_bins,
        short_expected=short_expected,
        short_std=short_std,
        intermediate_expected=intermediate_expected,
        intermediate_std=intermediate_std,
        long_expected=long_expected,
        long_std=long_std,
    )


def _compute_fsd_baseline(all_fsd_data: List[pd.DataFrame], fsd_paths: List[str] = None) -> Optional[FsdBaseline]:
    """
    Compute FSD baseline from sample data.
    
    Uses high-performance Rust implementation (3-10x faster) with Python fallback.
    FSD output format: DataFrame with 'region' column and size bin columns (e.g., '65-69', '70-74', ...)
    
    Aggregates size distribution proportions per chromosome arm across samples.
    """
    # Try Rust implementation first (3-10x faster)
    if fsd_paths:
        try:
            from krewlyzer import _core
            result = _core.pon_builder.compute_fsd_baseline(fsd_paths)
            if result:
                logger.info(f"FSD baseline computed (Rust): {len(result)} arms")
                arms = {}
                for arm, data in result.items():
                    arms[arm] = FsdArmBaseline(
                        size_bins=data["size_bins"],
                        expected=data["expected"],
                        std=data["std"],
                    )
                return FsdBaseline(arms=arms)
        except Exception as e:
            logger.debug(f"Rust FSD baseline computation failed, falling back to Python: {e}")
    
    # =========================================================================
    # PYTHON FALLBACK IMPLEMENTATION
    # NOTE: This is a fallback only. Primary implementation is Rust (3-10x faster).
    # The Rust implementation is in: rust/src/pon_builder.rs::compute_fsd_baseline()
    # =========================================================================
    
    if not all_fsd_data:
        return None
    
    # Filter to only DataFrames
    fsd_dfs = [df for df in all_fsd_data if isinstance(df, pd.DataFrame) and not df.empty]
    if not fsd_dfs:
        return None
    
    # Get size bin columns (all columns except 'region')
    sample_df = fsd_dfs[0]
    size_bin_cols = [col for col in sample_df.columns if col != 'region']
    
    # Convert column names to size bin integers (e.g., '65-69' -> 65)
    size_bins = []
    for col in size_bin_cols:
        try:
            size_bins.append(int(col.split('-')[0]))
        except ValueError:
            continue
    
    if not size_bins:
        return None
    
    # Aggregate per region across samples
    # First, concat all samples and compute mean/std per region
    all_data = pd.concat(fsd_dfs, ignore_index=True)
    
    # Group by region and compute stats
    arms = {}
    for region in all_data['region'].unique():
        region_data = all_data[all_data['region'] == region]
        expected = []
        std = []
        for col in size_bin_cols:
            values = region_data[col].values
            expected.append(float(np.mean(values)))
            std.append(float(np.std(values)) if len(values) > 1 else 0.0)
        arms[region] = {"expected": expected, "std": std}
    
    return FsdBaseline(size_bins=size_bins, arms=arms)


def _compute_wps_baseline(all_wps_data: List[pd.DataFrame], wps_paths: List[str] = None) -> Optional[WpsBaseline]:
    """
    Compute WPS baseline from Parquet vector format.
    
    Uses high-performance Rust implementation (3-10x faster) with Python fallback.
    For each region, computes mean and std vectors across all samples.
    This enables PoN subtraction with per-position z-scores.
    
    Expected columns: region_id, wps_nuc[200], wps_tf[200], etc.
    """
    import numpy as np
    
    # Try Rust implementation first (3-10x faster)
    if wps_paths:
        try:
            from krewlyzer import _core
            result = _core.pon_builder.compute_wps_baseline(wps_paths)
            if result:
                logger.info(f"WPS baseline computed (Rust): {len(result)} regions")
                # Convert to DataFrame for WpsBaseline
                rows = []
                for region_id, data in result.items():
                    rows.append({
                        "region_id": region_id,
                        "wps_nuc_mean": data.get("wps_long_mean", 0.0),
                        "wps_nuc_std": data.get("wps_long_std", 1.0),
                        "wps_tf_mean": data.get("wps_short_mean", 0.0),
                        "wps_tf_std": data.get("wps_short_std", 1.0),
                    })
                return WpsBaseline(regions=pd.DataFrame(rows))
        except Exception as e:
            logger.debug(f"Rust WPS baseline computation failed, falling back to Python: {e}")
    
    # =========================================================================
    # PYTHON FALLBACK IMPLEMENTATION
    # NOTE: This is a fallback only. Primary implementation is Rust (3-10x faster).
    # The Rust implementation is in: rust/src/pon_builder.rs::compute_wps_baseline()
    # =========================================================================
    
    if not all_wps_data:
        return None
    
    template = all_wps_data[0].copy()
    region_id_col = "gene_id" if "gene_id" in template.columns else ("region_id" if "region_id" in template.columns else "group_id")
    
    # Columns to aggregate - support both naming conventions
    # Rust WPS uses: wps_long, wps_short
    # Legacy naming: wps_nuc, wps_tf
    vector_cols = [c for c in ["wps_nuc", "wps_tf", "wps_long", "wps_short", "prot_frac_nuc", "prot_frac_tf"] 
                   if c in template.columns]
    
    if not vector_cols:
        logger.warning(f"No WPS vector columns found in columns: {list(template.columns)}")
        return None
    
    result_data = []
    
    for idx, row in template.iterrows():
        region_id = row[region_id_col]
        result_row = {region_id_col: region_id}
        
        # Copy non-vector metadata
        for col in ["chrom", "center", "strand", "region_type"]:
            if col in template.columns:
                result_row[col] = row[col]
        
        for col in vector_cols:
            # Collect values from all samples for this region
            values = []
            for df in all_wps_data:
                if idx < len(df) and col in df.columns:
                    val = df.iloc[idx][col]
                    if val is not None:
                        # Check if it's a scalar or array
                        if np.isscalar(val) or isinstance(val, (int, float, np.floating)):
                            values.append(float(val))
                        elif hasattr(val, '__len__') and len(val) > 0:
                            values.append(np.array(val, dtype=np.float32))
            
            if values:
                # Check if we have scalars or vectors
                if all(np.isscalar(v) or isinstance(v, (int, float, np.floating)) for v in values):
                    # Scalar values - just compute mean/std across samples
                    result_row[f"{col}_mean"] = float(np.mean(values))
                    result_row[f"{col}_std"] = float(np.std(values))
                else:
                    # Vector values - stack and compute per-bin stats
                    stacked = np.stack(values, axis=0)  # (n_samples, n_bins)
                    mean_vec = np.mean(stacked, axis=0)
                    std_vec = np.std(stacked, axis=0)
                    result_row[f"{col}_mean"] = mean_vec.tolist()
                    result_row[f"{col}_std"] = std_vec.tolist()
        
        result_data.append(result_row)
    
    result_df = pd.DataFrame(result_data)
    logger.info(f"WPS baseline: {len(result_df)} regions, {len(all_wps_data)} samples")
    
    return WpsBaseline(regions=result_df)


def _compute_ocf_baseline(all_ocf_data: List[pd.DataFrame]) -> "OcfBaseline":
    """
    Compute OCF baseline from sample data.
    
    Aggregates OCF scores per region across samples.
    """
    from .model import OcfBaseline
    
    if not all_ocf_data:
        return None
    
    # Concatenate all sample data
    combined = pd.concat(all_ocf_data, ignore_index=True)
    
    # Group by region and compute mean/std
    stats = combined.groupby("region_id").agg({
        "ocf": ["mean", "std"]
    }).reset_index()
    
    # Flatten column names
    stats.columns = ["region_id", "ocf_mean", "ocf_std"]
    stats["ocf_std"] = stats["ocf_std"].fillna(0.001)  # Handle single-sample case
    
    logger.info(f"OCF baseline: {len(stats)} regions, {len(all_ocf_data)} samples")
    
    return OcfBaseline(regions=stats)


def _compute_mds_baseline(all_mds_data: List[dict]) -> "MdsBaseline":
    """
    Compute MDS baseline from sample data.
    
    Aggregates k-mer frequencies and MDS scores across samples.
    """
    from .model import MdsBaseline
    
    if not all_mds_data:
        return None
    
    # Collect all k-mer frequencies
    all_kmers = set()
    for sample in all_mds_data:
        if sample.get("kmers"):
            all_kmers.update(sample["kmers"].keys())
    
    # Compute mean/std for each k-mer
    kmer_expected = {}
    kmer_std = {}
    
    for kmer in all_kmers:
        values = []
        for sample in all_mds_data:
            if sample.get("kmers") and kmer in sample["kmers"]:
                values.append(sample["kmers"][kmer])
        
        if values:
            kmer_expected[kmer] = np.mean(values)
            kmer_std[kmer] = np.std(values) if len(values) > 1 else 0.001
    
    # Compute MDS mean/std
    mds_values = [s["mds"] for s in all_mds_data if s.get("mds") is not None]
    mds_mean = np.mean(mds_values) if mds_values else 0.0
    mds_std = np.std(mds_values) if len(mds_values) > 1 else 1.0
    
    logger.info(f"MDS baseline: {len(kmer_expected)} k-mers, {len(all_mds_data)} samples")
    
    return MdsBaseline(
        kmer_expected=kmer_expected,
        kmer_std=kmer_std,
        mds_mean=mds_mean,
        mds_std=mds_std
    )


def _save_pon_model(model: PonModel, output: Path) -> None:
    """Save PON model to Parquet file."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    
    # Build metadata DataFrame
    metadata_df = pd.DataFrame([{
        "table": "metadata",
        "schema_version": model.schema_version,
        "assay": model.assay,
        "build_date": model.build_date,
        "n_samples": model.n_samples,
        "reference": model.reference,
        "panel_mode": model.panel_mode,
        "target_regions_file": model.target_regions_file,
    }])
    
    # Build GC bias DataFrame
    gc_bias_rows = []
    if model.gc_bias:
        for i, gc_bin in enumerate(model.gc_bias.gc_bins):
            gc_bias_rows.append({
                "table": "gc_bias",
                "gc_bin": gc_bin,
                "short_expected": model.gc_bias.short_expected[i],
                "short_std": model.gc_bias.short_std[i],
                "intermediate_expected": model.gc_bias.intermediate_expected[i],
                "intermediate_std": model.gc_bias.intermediate_std[i],
                "long_expected": model.gc_bias.long_expected[i],
                "long_std": model.gc_bias.long_std[i],
            })
    
    gc_bias_df = pd.DataFrame(gc_bias_rows) if gc_bias_rows else pd.DataFrame()
    
    # Build FSD baseline DataFrame
    fsd_rows = []
    if model.fsd_baseline:
        for arm, data in model.fsd_baseline.arms.items():
            for i, size_bin in enumerate(model.fsd_baseline.size_bins):
                if i < len(data["expected"]):
                    fsd_rows.append({
                        "table": "fsd_baseline",
                        "arm": arm,
                        "size_bin": size_bin,
                        "expected": data["expected"][i],
                        "std": data["std"][i],
                    })
    
    fsd_df = pd.DataFrame(fsd_rows) if fsd_rows else pd.DataFrame()
    
    # Build WPS baseline DataFrame
    wps_df = pd.DataFrame()
    if model.wps_baseline and model.wps_baseline.regions is not None:
        wps_df = model.wps_baseline.regions.copy()
        wps_df["table"] = "wps_baseline"
    
    # Build OCF baseline DataFrame
    ocf_df = pd.DataFrame()
    if model.ocf_baseline and model.ocf_baseline.regions is not None:
        ocf_df = model.ocf_baseline.regions.copy()
        ocf_df["table"] = "ocf_baseline"
    
    # Build MDS baseline DataFrame
    mds_df = pd.DataFrame()
    if model.mds_baseline:
        mds_df = pd.DataFrame([{
            "table": "mds_baseline",
            "mds_mean": model.mds_baseline.mds_mean,
            "mds_std": model.mds_baseline.mds_std,
            "kmer_expected": model.mds_baseline.kmer_expected,
            "kmer_std": model.mds_baseline.kmer_std,
        }])
    
    # Build on-target FSD baseline DataFrame (panel mode)
    fsd_ontarget_rows = []
    if model.fsd_baseline_ontarget:
        for arm, data in model.fsd_baseline_ontarget.arms.items():
            for i, size_bin in enumerate(model.fsd_baseline_ontarget.size_bins):
                if i < len(data["expected"]):
                    fsd_ontarget_rows.append({
                        "table": "fsd_baseline_ontarget",
                        "arm": arm,
                        "size_bin": size_bin,
                        "expected": data["expected"][i],
                        "std": data["std"][i],
                    })
    
    fsd_ontarget_df = pd.DataFrame(fsd_ontarget_rows) if fsd_ontarget_rows else pd.DataFrame()
    
    # Build on-target GC bias DataFrame (panel mode)
    gc_bias_ontarget_rows = []
    if model.gc_bias_ontarget:
        for i, gc_bin in enumerate(model.gc_bias_ontarget.gc_bins):
            gc_bias_ontarget_rows.append({
                "table": "gc_bias_ontarget",
                "gc_bin": gc_bin,
                "short_expected": model.gc_bias_ontarget.short_expected[i],
                "short_std": model.gc_bias_ontarget.short_std[i],
                "intermediate_expected": model.gc_bias_ontarget.intermediate_expected[i],
                "intermediate_std": model.gc_bias_ontarget.intermediate_std[i],
                "long_expected": model.gc_bias_ontarget.long_expected[i],
                "long_std": model.gc_bias_ontarget.long_std[i],
            })
    
    gc_bias_ontarget_df = pd.DataFrame(gc_bias_ontarget_rows) if gc_bias_ontarget_rows else pd.DataFrame()
    
    # Combine and save
    all_dfs = [metadata_df]
    if not gc_bias_df.empty:
        all_dfs.append(gc_bias_df)
    if not fsd_df.empty:
        all_dfs.append(fsd_df)
    if not wps_df.empty:
        all_dfs.append(wps_df)
    if not ocf_df.empty:
        all_dfs.append(ocf_df)
    if not mds_df.empty:
        all_dfs.append(mds_df)
    if not fsd_ontarget_df.empty:
        all_dfs.append(fsd_ontarget_df)
    if not gc_bias_ontarget_df.empty:
        all_dfs.append(gc_bias_ontarget_df)
    
    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_df.to_parquet(output, index=False)
    
    # Log summary
    n_gc = len(gc_bias_rows)
    n_fsd = len(fsd_rows)
    n_wps = len(wps_df) if not wps_df.empty else 0
    n_ocf = len(ocf_df) if not ocf_df.empty else 0
    n_mds = len(model.mds_baseline.kmer_expected) if model.mds_baseline else 0
    n_fsd_on = len(fsd_ontarget_rows)
    n_gc_on = len(gc_bias_ontarget_rows)
    logger.info(f"Saved PON model: {output}")
    if model.panel_mode:
        logger.info(f"   Panel mode: ON (targets: {model.target_regions_file})")
        if n_fsd_on > 0:
            logger.info(f"   FSD on-target: {n_fsd_on} arm×size entries")
        if n_gc_on > 0:
            logger.info(f"   GC on-target: {n_gc_on} bins")
    logger.info(f"   GC bias: {n_gc} bins")
    logger.info(f"   FSD baseline: {n_fsd} arm×size entries")
    logger.info(f"   WPS baseline: {n_wps} regions")
    logger.info(f"   OCF baseline: {n_ocf} regions")
    logger.info(f"   MDS baseline: {n_mds} k-mers")

