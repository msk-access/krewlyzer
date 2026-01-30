"""
PON model building from healthy plasma samples.

This module provides the CLI and logic for building unified PON models.
"""

import typer
from pathlib import Path
from typing import Optional, List, Dict
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

# Import asset resolution functions
from krewlyzer.core.asset_resolution import resolve_target_regions

# Import startup banner logging
from krewlyzer.core.logging import log_startup_banner, ResolvedAsset
from krewlyzer import __version__


def _process_sample_subprocess(
    sample_path: Path,
    output_dir: Path,
    sample_name: str,
    reference: Path,
    sample_params: 'SampleParams',
    target_regions: Optional[Path],
    assay: str,
    threads_per_sample: int,
    gene_bed: Optional[Path] = None,
) -> 'SampleOutputs':
    """
    Process a single sample in a subprocess.
    
    This wrapper configures the rayon thread pool for this subprocess,
    then calls process_sample. It's designed to be called via ProcessPoolExecutor.
    
    Args:
        sample_path: Path to BAM/CRAM or BED.gz file.
        output_dir: Directory for sample outputs.
        sample_name: Sample identifier.
        reference: Reference FASTA path.
        sample_params: Sample processing parameters.
        target_regions: Optional panel target regions.
        assay: Assay name for asset resolution.
        threads_per_sample: Threads to allocate to this sample's rayon pool.
        gene_bed: Optional gene BED for region-MDS.
    
    Returns:
        SampleOutputs with all extraction and feature results.
    """
    import logging
    from krewlyzer import _core
    from krewlyzer.core.sample_processor import process_sample
    
    # Configure rayon thread pool for THIS subprocess
    try:
        _core.configure_threads(threads_per_sample)
    except Exception:
        pass  # May fail if already configured, that's ok
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process sample
    outputs = process_sample(
        input_path=sample_path,
        output_dir=output_dir,
        sample_name=sample_name,
        reference=reference,
        params=sample_params,
        target_regions=target_regions,
        assay=assay,
        enable_fsc=True,
        enable_fsr=False,  # Not needed for PON
        enable_fsd=True,
        enable_wps=True,
        enable_ocf=True,
        pon_mode=True,  # Skip PON normalization (we're building the PON)
    )
    
    # Run region-mds for per-gene MDS baseline (BAM/CRAM only)
    is_bam_input = str(sample_path).endswith(('.bam', '.cram'))
    if is_bam_input and gene_bed and gene_bed.exists():
        try:
            mds_exon_output = output_dir / f"{sample_name}.MDS.exon.tsv"
            mds_gene_output = output_dir / f"{sample_name}.MDS.gene.tsv"
            
            _core.region_mds.run_region_mds(
                str(sample_path),           # BAM/CRAM path
                str(reference),             # Reference FASTA
                str(gene_bed),              # Gene BED file
                str(mds_exon_output),       # Per-exon output
                str(mds_gene_output),       # Per-gene output
                e1_only=False,
                mapq=sample_params.mapq,
                min_len=sample_params.minlen,
                max_len=400,  # Standard cfDNA range for motifs
                silent=True,
            )
        except Exception:
            pass  # Non-critical, continue without MDS
    
    return outputs


def build_pon(
    sample_list: Path = typer.Argument(..., help="Text file with paths to BAM/CRAM or BED.gz files (one per line). BAM/CRAM required for MDS baseline."),
    assay: str = typer.Option(..., "--assay", "-a", help="Assay name (e.g., msk-access-v2)"),
    reference: Path = typer.Option(..., "--reference", "-r", help="Reference FASTA file"),
    wps_anchors: Optional[Path] = typer.Option(None, "--wps-anchors", "-W", help="WPS anchors BED.gz for WPS baseline (default: bundled TSS+CTCF anchors)"),
    bin_file: Optional[Path] = typer.Option(None, "--bin-file", "-b", help="Bin file for FSC/FSR (default: hg19_window_100kb.bed)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output PON model file (.pon.parquet)"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="BED file with target regions (panel mode - builds dual on/off-target baselines)"),
    temp_dir: Optional[Path] = typer.Option(None, "--temp-dir", help="Directory for temporary files (default: system temp)"),
    threads: int = typer.Option(4, "--threads", "-p", help="Total threads (divided among parallel samples)"),
    parallel_samples: int = typer.Option(1, "--parallel-samples", "-P", help="Number of samples to process concurrently (0=auto, 1=sequential)"),
    memory_per_sample: float = typer.Option(2.0, "--memory-per-sample", help="Expected peak memory per sample in GB (for auto mode). Default: 2 GB for panels, use 16+ for WGS"),
    sample_timeout: int = typer.Option(0, "--sample-timeout", help="Max seconds per sample, 0=no timeout"),
    allow_failures: bool = typer.Option(False, "--allow-failures", help="Continue building PON even if some samples fail"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    require_proper_pair: bool = typer.Option(False, "--require-proper-pair", help="Only extract properly paired reads (default: False for v1 ACCESS compatibility)"),
    skip_target_regions: bool = typer.Option(False, "--skip-target-regions", help="Disable panel mode even when --assay has bundled targets (build as WGS)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
    validate_assets: bool = typer.Option(False, "--validate-assets", help="Validate bundled asset file formats before running"),
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
    
    # Validate user-provided override files
    from krewlyzer.core.asset_validation import validate_file, FileSchema
    if wps_anchors and wps_anchors.exists():
        logger.debug(f"Validating user-provided WPS anchors: {wps_anchors}")
        validate_file(wps_anchors, FileSchema.WPS_ANCHORS)
    if bin_file and bin_file.exists():
        logger.debug(f"Validating user-provided bin file: {bin_file}")
        validate_file(bin_file, FileSchema.BED3)
    if target_regions and target_regions.exists():
        logger.debug(f"Validating user-provided target regions: {target_regions}")
        validate_file(target_regions, FileSchema.BED3)
    
    # Read sample list
    with open(sample_list) as f:
        samples = [Path(line.strip()) for line in f if line.strip() and not line.startswith("#")]
    
    n_samples = len(samples)
    if n_samples < 1:
        logger.error("No samples found in sample list")
        raise typer.Exit(1)
    
    # Initialize AssetManager for bundled asset access
    from krewlyzer.assets import AssetManager
    assets = AssetManager(genome)
    if validate_assets:
        logger.info("Validating bundled assets (--validate-assets)...")
        assets.validate(assay=assay)
    
    # ═══════════════════════════════════════════════════════════════════
    # TARGET REGIONS RESOLUTION
    # ═══════════════════════════════════════════════════════════════════
    # Use shared resolution function for consistent behavior.
    # For build-pon, --assay is required for naming but can also auto-load targets.
    try:
        resolved_target_path, target_source = resolve_target_regions(
            explicit_path=target_regions,
            assay=assay,
            skip_target_regions=skip_target_regions,
            assets=assets,
            log=logger
        )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)
    
    # Panel mode detection (based on resolved targets)
    is_panel_mode = resolved_target_path is not None and resolved_target_path.exists()
    target_regions_str = str(resolved_target_path) if is_panel_mode else None
    
    # Count target regions for detail
    target_count_detail = ""
    if is_panel_mode:
        from krewlyzer.core.asset_resolution import get_target_region_count
        n_regions = get_target_region_count(resolved_target_path)
        if n_regions:
            target_count_detail = f"{n_regions} regions"
    
    # ═══════════════════════════════════════════════════════════════════
    # STARTUP BANNER
    # ═══════════════════════════════════════════════════════════════════
    log_startup_banner(
        tool_name="build-pon",
        version=__version__,
        inputs={
            "Sample list": str(sample_list),
            "Reference": str(reference.name),
            "Output": str(output),
            "Samples": str(n_samples),
        },
        config={
            "Genome": f"{assets.raw_genome} → {assets.genome_dir}",
            "Assay": assay,
            "Mode": "Panel" if is_panel_mode else "WGS",
        },
        assets=[
            ResolvedAsset("Targets", resolved_target_path, target_source, target_count_detail),
        ],
        logger=logger
    )
    
    # Default bin file
    if bin_file is None:
        bin_file = assets.bins_100kb
    
    if not bin_file.exists():
        logger.error(f"Bin file not found: {bin_file}")
        raise typer.Exit(1)
    
    # Default WPS anchors file
    wps_regions = wps_anchors
    if wps_regions is None:
        wps_regions = assets.wps_anchors
    
    if not wps_regions.exists():
        logger.warning(f"WPS regions file not found: {wps_regions}")
    
    # Default WPS background (Alu) file for NRL baseline
    wps_bg_file = assets.wps_background
    if wps_bg_file.exists():
        logger.debug(f"Using bundled WPS background: {wps_bg_file.name}")
    else:
        wps_bg_file = None
        logger.debug("WPS background file not found - Alu baseline will not be computed")
    
    import tempfile
    import shutil
    import time
    from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError as FuturesTimeoutError
    import multiprocessing as mp
    
    # Use 'spawn' context to avoid fork-related deadlocks with Rayon thread pools
    # on HPC systems. Fork + pre-initialized threads = deadlock.
    mp_context = mp.get_context('spawn')
    from krewlyzer.core.resource_utils import detect_resources, calculate_auto_parallel_samples
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UNIFIED SAMPLE PROCESSING (replaces Phase 1 + Phase 2)
    # Uses process_sample() from core/sample_processor.py for consistent processing
    # ═══════════════════════════════════════════════════════════════════════════
    
    # Determine parallel processing configuration
    if parallel_samples == 0:
        # Auto-detect based on resources
        detected_cpus, detected_ram = detect_resources()
        actual_parallel = calculate_auto_parallel_samples(
            threads=threads,
            memory_gb=detected_ram,
            memory_per_sample=memory_per_sample,
        )
        logger.info(f"Auto-detected resources: {detected_cpus} CPUs, {detected_ram:.1f} GB RAM")
        logger.info(f"Auto-selected {actual_parallel} parallel samples")
    else:
        actual_parallel = parallel_samples
    
    # Calculate threads per sample
    threads_per_sample = max(1, threads // max(1, actual_parallel))
    
    logger.info(f"=" * 60)
    logger.info(f"SAMPLE PROCESSING ({n_samples} samples)")
    logger.info(f"=" * 60)
    if actual_parallel > 1:
        logger.info(f"Parallel mode: {actual_parallel} concurrent samples, {threads_per_sample} threads each")
    else:
        logger.info(f"Sequential mode: {threads} threads")
    
    # Setup parameters for sample processing
    sample_params = SampleParams(
        genome=genome,
        mapq=20,
        minlen=65,
        maxlen=1000,
        skip_duplicates=True,
        require_proper_pair=require_proper_pair,
        kmer=4,
        threads=threads_per_sample,
        bin_file=bin_file,
        wps_anchors=wps_regions,  # Use resolved WPS anchors (modern) or transcript file (legacy)
        wps_background=wps_bg_file,  # Bundled Alu regions for NRL baseline
    )
    
    # Use temp directory for sample processing outputs
    temp_base = str(temp_dir) if temp_dir else None
    temp_output_dir = tempfile.mkdtemp(prefix="pon_build_", dir=temp_base)
    logger.info(f"Temporary output directory: {temp_output_dir}")
    
    # Process samples
    all_outputs: List[SampleOutputs] = []
    failed_samples: List[tuple] = []
    processing_start = time.time()

    
    try:
        # Resolve gene BED once for all samples
        gene_bed = assets.get_gene_bed_for_mode(assay)
        
        # Helper to extract sample name from path
        def get_sample_name(sample_path: Path) -> str:
            name = sample_path.stem
            for suffix in ['.bed', '.bam', '.cram']:
                name = name.replace(suffix, '')
            return name
        
        # Build sample info list
        sample_infos = [
            (sample_path, get_sample_name(sample_path))
            for sample_path in samples
        ]
        
        if actual_parallel > 1:
            # ─────────────────────────────────────────────────────────────────
            # PARALLEL PROCESSING
            # ─────────────────────────────────────────────────────────────────
            logger.info(f"Starting parallel processing with {actual_parallel} workers...")
            
            with ProcessPoolExecutor(max_workers=actual_parallel, mp_context=mp_context) as executor:
                # Submit all sample processing tasks
                futures = {
                    executor.submit(
                        _process_sample_subprocess,
                        sample_path,
                        Path(temp_output_dir) / sample_name,
                        sample_name,
                        reference,
                        sample_params,
                        resolved_target_path if is_panel_mode else None,
                        assay,
                        threads_per_sample,
                        gene_bed,
                    ): (sample_path, sample_name)
                    for sample_path, sample_name in sample_infos
                }
                
                completed = 0
                for future in as_completed(futures):
                    sample_path, sample_name = futures[future]
                    completed += 1
                    
                    try:
                        timeout_val = sample_timeout if sample_timeout > 0 else None
                        outputs = future.result(timeout=timeout_val)
                        all_outputs.append(outputs)
                        logger.info(f"[{completed}/{n_samples}] ✓ {sample_name} ({outputs.fragment_count:,} frags)")
                        
                    except FuturesTimeoutError:
                        error_msg = f"Timed out after {sample_timeout}s"
                        logger.error(f"[{completed}/{n_samples}] ✗ {sample_name}: {error_msg}")
                        failed_samples.append((sample_path, error_msg))
                        if not allow_failures:
                            raise RuntimeError(f"Sample {sample_name} timed out. Use --allow-failures to continue.")
                            
                    except Exception as e:
                        error_msg = str(e)
                        logger.error(f"[{completed}/{n_samples}] ✗ {sample_name}: {error_msg}")
                        failed_samples.append((sample_path, error_msg))
                        if not allow_failures:
                            raise RuntimeError(f"Sample {sample_name} failed: {error_msg}. Use --allow-failures to continue.")
            
            # Implicit barrier: ProcessPoolExecutor.__exit__ waits for all futures
            
        else:
            # ─────────────────────────────────────────────────────────────────
            # SEQUENTIAL PROCESSING (backward compatible)
            # ─────────────────────────────────────────────────────────────────
            for i, (sample_path, sample_name) in enumerate(sample_infos, 1):
                sample_start = time.time()
                logger.info(f"[{i}/{n_samples}] Processing: {sample_name}")
                
                try:
                    outputs = process_sample(
                        input_path=sample_path,
                        output_dir=Path(temp_output_dir) / sample_name,
                        sample_name=sample_name,
                        reference=reference,
                        params=sample_params,
                        target_regions=resolved_target_path if is_panel_mode else None,
                        assay=assay,
                        enable_fsc=True,
                        enable_fsr=False,
                        enable_fsd=True,
                        enable_wps=True,
                        enable_ocf=True,
                        pon_mode=True,
                    )
                    
                    all_outputs.append(outputs)
                    elapsed = time.time() - sample_start
                    logger.info(f"  ✓ Done in {elapsed:.1f}s ({outputs.fragment_count:,} frags)")
                    
                    # Run region-mds for per-gene MDS baseline (BAM/CRAM only)
                    is_bam_input = str(sample_path).endswith(('.bam', '.cram'))
                    if is_bam_input and gene_bed and gene_bed.exists():
                        try:
                            sample_out_dir = Path(temp_output_dir) / sample_name
                            _core.region_mds.run_region_mds(
                                str(sample_path),
                                str(reference),
                                str(gene_bed),
                                str(sample_out_dir / f"{sample_name}.MDS.exon.tsv"),
                                str(sample_out_dir / f"{sample_name}.MDS.gene.tsv"),
                                e1_only=False,
                                mapq=sample_params.mapq,
                                min_len=sample_params.minlen,
                                max_len=400,
                                silent=True,
                            )
                        except Exception as e:
                            logger.debug(f"    Could not run region-mds: {e}")
                    
                except Exception as e:
                    logger.warning(f"  ✗ Failed: {e}")
                    failed_samples.append((sample_path, str(e)))
                    if not allow_failures:
                        raise RuntimeError(f"Sample {sample_name} failed: {e}. Use --allow-failures to continue.")
        
        # ─────────────────────────────────────────────────────────────────
        # PROCESSING COMPLETE - CHECK RESULTS
        # ─────────────────────────────────────────────────────────────────
        processing_elapsed = time.time() - processing_start
        logger.info(f"-" * 60)
        logger.info(f"Processing complete: {len(all_outputs)} succeeded, {len(failed_samples)} failed ({processing_elapsed:.1f}s)")
        

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

        all_ocf_data = []
        all_mds_data = []
        all_tfbs_data = []
        all_atac_data = []
        all_fsd_data_ontarget = []
        all_gc_data_ontarget = []
        
        # On-target entropy data collectors (panel mode only)
        # These build separate baselines for on-target regions which have
        # fundamentally different depth/variance due to capture enrichment
        all_tfbs_data_ontarget = []
        all_atac_data_ontarget = []
        
        
        # Collect file paths for Rust-based baseline computation
        fsd_paths = []
        wps_paths = []
        fsd_ontarget_paths = []  # On-target FSD for panel mode
        wps_background_paths = []  # WPS Alu background for NRL baseline
        wps_panel_paths = []  # Panel-specific WPS parquets (panel mode only)
        mds_gene_paths = []  # Region MDS per-gene files
        
        # FSC gene/region data collectors (panel mode only)
        # Collect normalized_depth from FSC.gene.tsv and FSC.regions.tsv
        all_fsc_gene_data = []  # List[Dict[gene, depth]]
        all_fsc_region_data = []  # List[Dict[region_id, depth]]
        
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
            
            # Collect on-target GC data for gc_bias_ontarget (panel mode)
            # On-target FSC counts have the same format with mean_gc column
            fsc_counts_on = Path(temp_output_dir) / sample_name / f"{sample_name}.fsc_counts.ontarget.tsv"
            if fsc_counts_on.exists():
                try:
                    fsc_on_df = pd.read_csv(fsc_counts_on, sep='\t')
                    if all(col in fsc_on_df.columns for col in ['ultra_short', 'core_short', 'mono_nucl', 'di_nucl', 'long', 'mean_gc']):
                        all_gc_data_ontarget.append({
                            "gc": fsc_on_df['mean_gc'].values,
                            "short": (fsc_on_df['ultra_short'] + fsc_on_df['core_short']).values,
                            "intermediate": fsc_on_df['mono_nucl'].values,
                            "long": (fsc_on_df['di_nucl'] + fsc_on_df['long']).values,
                        })
                        logger.debug(f"  Collected on-target GC data from {sample_name}")
                except Exception as e:
                    logger.debug(f"  Could not read on-target GC data for {sample_name}: {e}")
            
            # Collect FSD data
            if feat and feat.fsd:
                if feat.fsd.exists():
                    fsd_paths.append(str(feat.fsd))  # Collect path for Rust baseline
                    try:
                        fsd_df = pd.read_csv(feat.fsd, sep='\t')
                        all_fsd_data.append(fsd_df)
                    except Exception as e:
                        logger.debug(f"  Could not read FSD for {sample_name}: {e}")
                else:
                    logger.warning(f"  FSD path set but file missing: {feat.fsd}")
            
            # Collect on-target FSD (panel mode)
            if feat and feat.fsd_ontarget and feat.fsd_ontarget.exists():
                fsd_ontarget_paths.append(str(feat.fsd_ontarget))  # Collect path for Rust baseline
                try:
                    fsd_on_df = pd.read_csv(feat.fsd_ontarget, sep='\t')
                    all_fsd_data_ontarget.append(fsd_on_df)
                except Exception as e:
                    logger.debug(f"  Could not read on-target FSD for {sample_name}: {e}")
            
            # Collect WPS paths for Rust v2.0 vector baseline
            if feat and feat.wps and feat.wps.exists():
                wps_paths.append(str(feat.wps))
            
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
            
            # Collect WPS background (Alu stacking) data
            if feat and feat.wps_background and feat.wps_background.exists():
                wps_background_paths.append(str(feat.wps_background))
            
            # Collect WPS panel parquets (panel mode only)
            if feat and feat.wps_panel and feat.wps_panel.exists():
                wps_panel_paths.append(str(feat.wps_panel))
            
            # Collect region MDS gene paths (generated during processing)
            sample_mds_gene = Path(temp_output_dir) / sample_name / f"{sample_name}.MDS.gene.tsv"
            if sample_mds_gene.exists():
                mds_gene_paths.append(str(sample_mds_gene))
            
            # Collect MDS data (from extraction)
            if outputs.mds_counts:
                all_mds_data.append({
                    "kmers": outputs.mds_counts,
                    "mds": outputs.mds_score
                })
            
            # Collect TFBS entropy data (off-target / WGS-like)
            if outputs.tfbs_data:
                all_tfbs_data.append(outputs.tfbs_data)
            
            # Collect ATAC entropy data (off-target / WGS-like)
            if outputs.atac_data:
                all_atac_data.append(outputs.atac_data)
            
            # Collect on-target TFBS entropy data (panel mode)
            # Uses panel-specific regions (pre-intersected with targets) for higher signal
            if outputs.tfbs_data_ontarget:
                all_tfbs_data_ontarget.append(outputs.tfbs_data_ontarget)
                
            # Collect on-target ATAC entropy data (panel mode)
            if outputs.atac_data_ontarget:
                all_atac_data_ontarget.append(outputs.atac_data_ontarget)
            
            # Collect FSC gene data (panel mode - normalized depth per gene)
            if feat and feat.fsc_gene and feat.fsc_gene.exists():
                try:
                    fsc_gene_df = pd.read_csv(feat.fsc_gene, sep='\t')
                    if 'gene' in fsc_gene_df.columns and 'normalized_depth' in fsc_gene_df.columns:
                        gene_depths = dict(zip(fsc_gene_df['gene'], fsc_gene_df['normalized_depth']))
                        all_fsc_gene_data.append(gene_depths)
                        logger.debug(f"  Collected FSC gene: {len(gene_depths)} genes from {sample_name}")
                except Exception as e:
                    logger.debug(f"  Could not read FSC gene for {sample_name}: {e}")
            
            # Collect FSC region data (panel mode - normalized depth per exon/probe)
            if feat and feat.fsc_region and feat.fsc_region.exists():
                try:
                    fsc_region_df = pd.read_csv(feat.fsc_region, sep='\t')
                    if all(c in fsc_region_df.columns for c in ['chrom', 'start', 'end', 'normalized_depth']):
                        fsc_region_df['region_id'] = (
                            fsc_region_df['chrom'].astype(str) + ':' + 
                            fsc_region_df['start'].astype(str) + '-' + 
                            fsc_region_df['end'].astype(str)
                        )
                        region_depths = dict(zip(fsc_region_df['region_id'], fsc_region_df['normalized_depth']))
                        all_fsc_region_data.append(region_depths)
                        logger.debug(f"  Collected FSC region: {len(region_depths)} regions from {sample_name}")
                except Exception as e:
                    logger.debug(f"  Could not read FSC region for {sample_name}: {e}")
        
        logger.info(f"  GC samples: {len(all_gc_data)}")
        logger.info(f"  FSD samples: {len(all_fsd_data)}")
        logger.info(f"  WPS sample parquets: {len(wps_paths)}")
        logger.info(f"  OCF samples: {len(all_ocf_data)}")
        logger.info(f"  MDS samples: {len(all_mds_data)}")
        logger.info(f"  TFBS samples: {len(all_tfbs_data)}")
        logger.info(f"  ATAC samples: {len(all_atac_data)}")
        if is_panel_mode:
            logger.info(f"  FSD on-target samples: {len(all_fsd_data_ontarget)}")
            logger.info(f"  FSC gene samples: {len(all_fsc_gene_data)}")
            logger.info(f"  FSC region samples: {len(all_fsc_region_data)}")
    
    except Exception as e:
        # On error, cleanup temp directory before re-raising
        if temp_output_dir and Path(temp_output_dir).exists():
            shutil.rmtree(temp_output_dir)
            logger.debug(f"Cleaned up temp directory after error: {temp_output_dir}")
        raise
    
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
    fsd_baseline = _compute_fsd_baseline(all_fsd_data, fsd_paths)
    
    # Build WPS baseline
    logger.info("  Computing WPS baseline...")
    wps_baseline = _compute_wps_baseline(wps_paths)
    
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
    
    # Build TFBS baseline (off-target / WGS-like)
    tfbs_baseline = None
    if all_tfbs_data:
        logger.info("  Computing TFBS baseline (off-target)...")
        tfbs_baseline = _compute_tfbs_baseline(all_tfbs_data)
    
    # Build ATAC baseline (off-target / WGS-like)
    atac_baseline = None
    if all_atac_data:
        logger.info("  Computing ATAC baseline (off-target)...")
        atac_baseline = _compute_atac_baseline(all_atac_data)
    
    # =========================================================================
    # Build on-target TFBS/ATAC baselines (panel mode only)
    # On-target regions have 1000x-5000x depth vs 0.5x-1x off-target
    # Using off-target baseline for on-target z-scores would artificially 
    # inflate scores due to variance mismatch
    # =========================================================================
    tfbs_baseline_ontarget = None
    if is_panel_mode and all_tfbs_data_ontarget:
        logger.info(f"  Computing on-target TFBS baseline ({len(all_tfbs_data_ontarget)} samples)...")
        tfbs_baseline_ontarget = _compute_tfbs_baseline(all_tfbs_data_ontarget)
    
    atac_baseline_ontarget = None
    if is_panel_mode and all_atac_data_ontarget:
        logger.info(f"  Computing on-target ATAC baseline ({len(all_atac_data_ontarget)} samples)...")
        atac_baseline_ontarget = _compute_atac_baseline(all_atac_data_ontarget)
    
    
    # Build on-target FSD baseline (panel mode only)
    fsd_baseline_ontarget = None
    if is_panel_mode and all_fsd_data_ontarget:
        logger.info("  Computing on-target FSD baseline...")
        fsd_baseline_ontarget = _compute_fsd_baseline(all_fsd_data_ontarget, fsd_ontarget_paths)
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
    
    # Build WPS background baseline (Alu NRL/periodicity)
    wps_background_baseline = None
    if wps_background_paths:
        logger.info("  Computing WPS background baseline (Alu NRL)...")
        wps_background_baseline = _compute_wps_background_baseline(wps_background_paths)
    
    # Build WPS panel baseline (panel mode only)
    # Uses panel-specific anchors overlapping target regions
    wps_baseline_panel = None
    if is_panel_mode and wps_panel_paths:
        logger.info(f"  Computing WPS panel baseline ({len(wps_panel_paths)} samples)...")
        wps_baseline_panel = _compute_wps_baseline(wps_panel_paths)
    
    # Build Region MDS baseline (per-gene MDS statistics)
    region_mds_baseline = None
    if mds_gene_paths:
        logger.info("  Computing Region MDS baseline...")
        try:
            region_mds_baseline = _compute_region_mds_baseline(mds_gene_paths)
        except Exception as e:
            logger.warning(f"  Region MDS baseline failed: {e}")
    
    # Compute FSC gene/region baselines (panel mode only)
    fsc_gene_baseline = None
    fsc_region_baseline = None
    if is_panel_mode and all_fsc_gene_data:
        logger.info("  Computing FSC gene baseline...")
        try:
            fsc_gene_baseline = _compute_fsc_gene_baseline(all_fsc_gene_data)
        except Exception as e:
            logger.warning(f"  FSC gene baseline failed: {e}")
    
    if is_panel_mode and all_fsc_region_data:
        logger.info("  Computing FSC region baseline...")
        try:
            fsc_region_baseline = _compute_fsc_region_baseline(all_fsc_region_data)
        except Exception as e:
            logger.warning(f"  FSC region baseline failed: {e}")
    
    # Create PON model
    model = PonModel(
        schema_version="1.0",
        assay=assay,
        build_date=datetime.now().isoformat()[:10],
        n_samples=len(all_gc_data),
        reference=reference.name.replace(".fa", "").replace(".fasta", ""),
        panel_mode=is_panel_mode,
        target_regions_file=resolved_target_path.name if is_panel_mode else "",
        gc_bias=gc_bias,
        fsd_baseline=fsd_baseline,
        wps_baseline=wps_baseline,
        wps_background_baseline=wps_background_baseline,
        wps_baseline_panel=wps_baseline_panel,
        ocf_baseline=ocf_baseline,
        mds_baseline=mds_baseline,
        region_mds=region_mds_baseline,
        tfbs_baseline=tfbs_baseline,
        atac_baseline=atac_baseline,
        # On-target baselines (panel mode - uses panel-specific regions)
        tfbs_baseline_ontarget=tfbs_baseline_ontarget,
        atac_baseline_ontarget=atac_baseline_ontarget,
        fsc_gene_baseline=fsc_gene_baseline,
        fsc_region_baseline=fsc_region_baseline,
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
    logger.info(f"=" * 60)
    logger.info(f"PON MODEL SUMMARY")
    logger.info(f"=" * 60)
    logger.info(f"  Assay: {model.assay}")
    logger.info(f"  Samples: {model.n_samples}")
    logger.info(f"  Panel mode: {is_panel_mode}")
    logger.info(f"  Output: {output}")
    logger.info(f"")
    logger.info(f"BASELINES GENERATED:")
    logger.info(f"-" * 40)
    
    # Off-target baselines
    def _baseline_status(name, baseline, detail_fn=None):
        if baseline is None:
            return f"  {name}: ❌ Not generated"
        detail = detail_fn(baseline) if detail_fn else ""
        return f"  {name}: ✅ {detail}"
    
    logger.info(_baseline_status("gc_bias", gc_bias, 
        lambda b: f"{len(b.gc_bins)} bins" if hasattr(b, 'gc_bins') else "OK"))
    logger.info(_baseline_status("fsd_baseline", fsd_baseline,
        lambda b: f"{len(b.arms)} arms" if hasattr(b, 'arms') else "OK"))
    logger.info(_baseline_status("wps_baseline", wps_baseline,
        lambda b: f"{len(b.regions)} regions" if hasattr(b, 'regions') else "OK"))
    logger.info(_baseline_status("wps_background", wps_background_baseline,
        lambda b: f"{len(b.groups)} groups" if hasattr(b, 'groups') else "OK"))
    logger.info(_baseline_status("ocf_baseline", ocf_baseline,
        lambda b: f"{len(b.regions)} regions" if hasattr(b, 'regions') else "OK"))
    logger.info(_baseline_status("mds_baseline", mds_baseline, lambda b: "OK"))
    logger.info(_baseline_status("region_mds", region_mds_baseline,
        lambda b: f"{len(b.gene_baseline)} genes" if hasattr(b, 'gene_baseline') else "OK"))
    logger.info(_baseline_status("tfbs_baseline", tfbs_baseline,
        lambda b: f"{len(b.labels)} labels" if hasattr(b, 'labels') else "OK"))
    logger.info(_baseline_status("atac_baseline", atac_baseline,
        lambda b: f"{len(b.labels)} labels" if hasattr(b, 'labels') else "OK"))
    
    # Panel-mode baselines
    if is_panel_mode:
        logger.info(f"")
        logger.info(f"ON-TARGET BASELINES (Panel Mode):")
        logger.info(f"-" * 40)
        logger.info(_baseline_status("gc_bias_ontarget", gc_bias_ontarget,
            lambda b: f"{len(b.gc_bins)} bins" if hasattr(b, 'gc_bins') else "OK"))
        logger.info(_baseline_status("fsd_baseline_ontarget", fsd_baseline_ontarget,
            lambda b: f"{len(b.arms)} arms" if hasattr(b, 'arms') else "OK"))
        # On-target Region Entropy baselines (capture-specific distributions)
        logger.info(_baseline_status("tfbs_baseline_ontarget", tfbs_baseline_ontarget,
            lambda b: f"{len(b.labels)} labels" if hasattr(b, 'labels') else "OK"))
        logger.info(_baseline_status("atac_baseline_ontarget", atac_baseline_ontarget,
            lambda b: f"{len(b.labels)} labels" if hasattr(b, 'labels') else "OK"))
        logger.info(_baseline_status("wps_baseline_panel", wps_baseline_panel,
            lambda b: f"{len(b.regions)} regions" if hasattr(b, 'regions') else "OK"))
        logger.info(_baseline_status("fsc_gene_baseline", fsc_gene_baseline,
            lambda b: f"{len(b)} genes" if b else "OK"))
        logger.info(_baseline_status("fsc_region_baseline", fsc_region_baseline,
            lambda b: f"{len(b)} regions" if b else "OK"))
        
    
    logger.info(f"")
    logger.info(f"=" * 60)
    logger.info(f"✅ PON model built successfully: {output}")
    
    # Cleanup temp directory after successful completion
    if temp_output_dir and Path(temp_output_dir).exists():
        shutil.rmtree(temp_output_dir)
        logger.debug(f"Cleaned up temp directory: {temp_output_dir}")


def _compute_gc_bias_model(all_gc_data: List[dict]) -> GcBiasModel:
    """
    Compute GC bias curves from sample data.
    
    Uses Rust implementation for high performance.
    For each GC bin, compute median expected coverage across samples.
    
    Raises:
        RuntimeError: If computation fails
    """
    from krewlyzer import _core
    
    result = _core.pon_builder.compute_gc_bias_model(all_gc_data)
    if not result or "gc_bins" not in result:
        raise RuntimeError("GC bias model computation failed: no data returned from Rust")
    
    logger.info(f"GC bias model computed: {len(result['gc_bins'])} bins")
    return GcBiasModel(
        gc_bins=result["gc_bins"],
        short_expected=result["short_expected"],
        short_std=result["short_std"],
        intermediate_expected=result["intermediate_expected"],
        intermediate_std=result["intermediate_std"],
        long_expected=result["long_expected"],
        long_std=result["long_std"],
    )


def _compute_fsd_baseline(all_fsd_data: List[pd.DataFrame], fsd_paths: List[str] = None) -> Optional[FsdBaseline]:
    """
    Compute FSD baseline from sample data.
    
    Uses Rust implementation for high performance.
    FSD output format: DataFrame with 'region' column and size bin columns.
    
    Args:
        all_fsd_data: List of FSD DataFrames (unused, kept for API compat)
        fsd_paths: List of paths to FSD TSV files
        
    Returns:
        FsdBaseline or None if no data
        
    Raises:
        RuntimeError: If computation fails
    """
    if not fsd_paths:
        logger.warning("No FSD paths provided for baseline computation")
        return None
    
    from krewlyzer import _core
    
    result = _core.pon_builder.compute_fsd_baseline(fsd_paths)
    if not result:
        raise RuntimeError("FSD baseline computation failed: no data returned from Rust")
    
    logger.info(f"FSD baseline computed: {len(result)} arms")
    arms = {}
    for arm, data in result.items():
        arms[arm] = {"expected": data["expected"], "std": data["std"]}
    
    # Get size bins from first arm
    first_arm = next(iter(result.values()))
    size_bins = first_arm.get("size_bins", list(range(65, 1000, 5)))
    
    return FsdBaseline(size_bins=size_bins, arms=arms)


def _compute_wps_baseline(wps_paths: List[str]) -> Optional[WpsBaseline]:
    """
    Compute WPS baseline from Parquet vector format (v2.0).
    
    Uses Rust implementation for high performance.
    For each region, computes element-wise mean and std vectors across all samples.
    
    Args:
        wps_paths: List of paths to WPS Parquet files
        
    Returns:
        WpsBaseline with 200-element vectors or None if no data
        
    Raises:
        RuntimeError: If computation fails
    """
    if not wps_paths:
        logger.warning("No WPS paths provided for baseline computation")
        return None
    
    from krewlyzer import _core
    
    logger.info(f"Computing WPS vector baseline from {len(wps_paths)} samples")
    result = _core.pon_builder.compute_wps_baseline(wps_paths)
    
    if not result:
        raise RuntimeError("WPS baseline computation failed: no data returned from Rust")
    
    logger.info(f"WPS vector baseline computed: {len(result)} regions")
    
    # Convert to DataFrame for WpsBaseline (v2.0 vector format)
    rows = []
    for region_id, data in result.items():
        rows.append({
            "region_id": region_id,
            "wps_nuc_mean": data.get("wps_nuc_mean", []),  # 200-element vector
            "wps_nuc_std": data.get("wps_nuc_std", []),    # 200-element vector
            "wps_tf_mean": data.get("wps_tf_mean", []),    # 200-element vector
            "wps_tf_std": data.get("wps_tf_std", []),      # 200-element vector
            "n_samples": data.get("n_samples", 0),
        })
    
    return WpsBaseline(regions=pd.DataFrame(rows), schema_version="2.0")


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


def _compute_region_mds_baseline(mds_gene_paths: List[str]) -> "RegionMdsBaseline":
    """
    Compute Region MDS baseline from sample MDS.gene.tsv files.
    
    Uses Rust implementation for high performance.
    Aggregates per-gene MDS values across all samples.
    
    Args:
        mds_gene_paths: List of paths to MDS.gene.tsv files
        
    Returns:
        RegionMdsBaseline with per-gene mean/std MDS values
        
    Raises:
        RuntimeError: If computation fails
    """
    from .model import RegionMdsBaseline
    
    if not mds_gene_paths:
        return None
    
    # Filter to existing paths
    valid_paths = [p for p in mds_gene_paths if Path(p).exists()]
    if not valid_paths:
        logger.warning("No valid MDS.gene.tsv files found for region-MDS baseline")
        return None
    
    logger.info(f"Computing region-MDS baseline from {len(valid_paths)} samples...")
    
    try:
        # Use Rust-accelerated aggregation
        result = _core.pon_builder.compute_region_mds_baseline(valid_paths)
        
        if not result:
            logger.warning("Region-MDS baseline computation returned empty result")
            return None
        
        # Convert to RegionMdsBaseline
        gene_baseline = {}
        for gene, data in result.items():
            gene_baseline[gene] = {
                "mds_mean": data.get("mds_mean", 0.0),
                "mds_std": data.get("mds_std", 1.0),
                "mds_e1_mean": data.get("mds_e1_mean", 0.0),
                "mds_e1_std": data.get("mds_e1_std", 1.0),
                "n_samples": data.get("n_samples", 0),
            }
        
        logger.info(f"Region-MDS baseline: {len(gene_baseline)} genes")
        
        return RegionMdsBaseline(gene_baseline=gene_baseline)
        
    except Exception as e:
        logger.error(f"Region-MDS baseline computation failed: {e}")
        raise RuntimeError(f"Region-MDS baseline computation failed: {e}")


def _compute_wps_background_baseline(wps_background_paths: List[str]) -> "WpsBackgroundBaseline":
    """
    Compute WPS background (Alu) baseline from sample WPS_background.parquet files.
    
    Aggregates nucleosome repeat length (NRL) and periodicity values across samples
    for Alu element stacking analysis.
    
    Args:
        wps_background_paths: List of paths to WPS_background.parquet files
        
    Returns:
        WpsBackgroundBaseline with per-group NRL/periodicity mean/std
    """
    from .model import WpsBackgroundBaseline
    
    if not wps_background_paths:
        return None
    
    # Filter to existing paths
    valid_paths = [p for p in wps_background_paths if Path(p).exists()]
    if not valid_paths:
        logger.warning("No valid WPS_background.parquet files found")
        return None
    
    logger.info(f"Computing WPS background baseline from {len(valid_paths)} samples...")
    
    try:
        # Load and aggregate WPS background data across samples
        all_groups = []
        for path in valid_paths:
            try:
                df = pd.read_parquet(path)
                # Expected columns: group_id, nrl, periodicity (or similar)
                if 'group_id' in df.columns:
                    all_groups.append(df)
            except Exception as e:
                logger.debug(f"Could not read {path}: {e}")
        
        if not all_groups:
            logger.warning("No valid WPS background data found")
            return None
        
        # Concatenate all samples
        combined = pd.concat(all_groups, ignore_index=True)
        
        # Aggregate by group_id
        group_stats = []
        for group_id in combined['group_id'].unique():
            group_data = combined[combined['group_id'] == group_id]
            
            nrl_col = 'nrl' if 'nrl' in group_data.columns else 'nucleosome_repeat_length'
            periodicity_col = 'periodicity' if 'periodicity' in group_data.columns else 'period_score'
            
            nrl_mean = group_data[nrl_col].mean() if nrl_col in group_data.columns else 167.0
            nrl_std = group_data[nrl_col].std() if nrl_col in group_data.columns else 5.0
            period_mean = group_data[periodicity_col].mean() if periodicity_col in group_data.columns else 0.0
            period_std = group_data[periodicity_col].std() if periodicity_col in group_data.columns else 1.0
            
            group_stats.append({
                'group_id': group_id,
                'nrl_mean': nrl_mean,
                'nrl_std': max(nrl_std, 0.1),  # Avoid zero std
                'periodicity_mean': period_mean,
                'periodicity_std': max(period_std, 0.01),
            })
        
        groups_df = pd.DataFrame(group_stats)
        logger.info(f"WPS background baseline: {len(groups_df)} groups")
        
        return WpsBackgroundBaseline(groups=groups_df)
        
    except Exception as e:
        logger.error(f"WPS background baseline computation failed: {e}")
        return None


def _compute_fsc_gene_baseline(all_fsc_gene_data: List[Dict]) -> "FscGeneBaseline":
    """
    Compute FSC gene depth baseline from sample FSC.gene.tsv data.
    
    Aggregates normalized_depth values per gene across all samples.
    Requires minimum 3 samples per gene for reliable statistics.
    
    Args:
        all_fsc_gene_data: List of dicts, each {gene: normalized_depth}
        
    Returns:
        FscGeneBaseline with per-gene mean/std/n_samples
        
    Note:
        Panel mode only - requires --assay parameter during build-pon.
    """
    from .model import FscGeneBaseline
    from collections import defaultdict
    
    if not all_fsc_gene_data:
        logger.debug("No FSC gene data provided for baseline computation")
        return None
    
    logger.info(f"Computing FSC gene baseline from {len(all_fsc_gene_data)} samples...")
    
    # Aggregate depths per gene across all samples
    gene_values = defaultdict(list)
    for sample_data in all_fsc_gene_data:
        for gene, depth in sample_data.items():
            if depth is not None and not np.isnan(depth):
                gene_values[gene].append(depth)
    
    # Compute mean/std for genes with sufficient samples
    MIN_SAMPLES = 3
    data = {}
    skipped = 0
    
    for gene, values in gene_values.items():
        if len(values) >= MIN_SAMPLES:
            mean_depth = float(np.mean(values))
            std_depth = float(np.std(values))
            # Ensure minimum std to avoid division by zero
            std_depth = max(std_depth, 0.001)
            data[gene] = (mean_depth, std_depth, len(values))
        else:
            skipped += 1
    
    if skipped > 0:
        logger.debug(f"FSC gene: skipped {skipped} genes with <{MIN_SAMPLES} samples")
    
    logger.info(f"FSC gene baseline: {len(data)} genes (min {MIN_SAMPLES} samples)")
    
    return FscGeneBaseline(data=data) if data else None


def _compute_fsc_region_baseline(all_fsc_region_data: List[Dict]) -> "FscRegionBaseline":
    """
    Compute FSC region (exon/probe) depth baseline from sample FSC.regions.tsv data.
    
    Aggregates normalized_depth values per region across all samples.
    Covers all exons (no filtering by variance).
    Requires minimum 3 samples per region for reliable statistics.
    
    Args:
        all_fsc_region_data: List of dicts, each {region_id: normalized_depth}
        
    Returns:
        FscRegionBaseline with per-region mean/std/n_samples
        
    Note:
        Region IDs are formatted as "chrom:start-end".
        Panel mode only - requires --assay parameter during build-pon.
    """
    from .model import FscRegionBaseline
    from collections import defaultdict
    
    if not all_fsc_region_data:
        logger.debug("No FSC region data provided for baseline computation")
        return None
    
    logger.info(f"Computing FSC region baseline from {len(all_fsc_region_data)} samples...")
    
    # Aggregate depths per region across all samples
    region_values = defaultdict(list)
    for sample_data in all_fsc_region_data:
        for region_id, depth in sample_data.items():
            if depth is not None and not np.isnan(depth):
                region_values[region_id].append(depth)
    
    # Compute mean/std for regions with sufficient samples
    MIN_SAMPLES = 3
    data = {}
    skipped = 0
    
    for region_id, values in region_values.items():
        if len(values) >= MIN_SAMPLES:
            mean_depth = float(np.mean(values))
            std_depth = float(np.std(values))
            # Ensure minimum std to avoid division by zero
            std_depth = max(std_depth, 0.001)
            data[region_id] = (mean_depth, std_depth, len(values))
        else:
            skipped += 1
    
    if skipped > 0:
        logger.debug(f"FSC region: skipped {skipped} regions with <{MIN_SAMPLES} samples")
    
    logger.info(f"FSC region baseline: {len(data)} regions (min {MIN_SAMPLES} samples)")
    
    return FscRegionBaseline(data=data) if data else None


def _compute_tfbs_baseline(all_tfbs_data: List[dict]) -> "TfbsBaseline":
    """
    Compute TFBS entropy baseline from sample data.
    
    Aggregates entropy values per TF label across samples.
    
    Args:
        all_tfbs_data: List of dicts, each {label: entropy_value}
        
    Returns:
        TfbsBaseline with per-label mean/std entropy
    """
    from .model import TfbsBaseline
    from krewlyzer.core.region_entropy_processor import RegionEntropyBaseline
    
    if not all_tfbs_data:
        return None
    
    # Build baseline using RegionEntropyBaseline.from_samples()
    baseline = RegionEntropyBaseline.from_samples(all_tfbs_data)
    
    logger.info(f"TFBS baseline: {len(baseline.data)} labels, {len(all_tfbs_data)} samples")
    
    return TfbsBaseline(baseline=baseline)


def _compute_atac_baseline(all_atac_data: List[dict]) -> "AtacBaseline":
    """
    Compute ATAC entropy baseline from sample data.
    
    Aggregates entropy values per cancer type label across samples.
    
    Args:
        all_atac_data: List of dicts, each {label: entropy_value}
        
    Returns:
        AtacBaseline with per-label mean/std entropy
    """
    from .model import AtacBaseline
    from krewlyzer.core.region_entropy_processor import RegionEntropyBaseline
    
    if not all_atac_data:
        return None
    
    # Build baseline using RegionEntropyBaseline.from_samples()
    baseline = RegionEntropyBaseline.from_samples(all_atac_data)
    
    logger.info(f"ATAC baseline: {len(baseline.data)} labels, {len(all_atac_data)} samples")
    
    return AtacBaseline(baseline=baseline)

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
    
    # Build WPS panel baseline DataFrame (panel mode only)
    wps_panel_df = pd.DataFrame()
    if model.wps_baseline_panel and model.wps_baseline_panel.regions is not None:
        wps_panel_df = model.wps_baseline_panel.regions.copy()
        wps_panel_df["table"] = "wps_baseline_panel"
    
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
    
    # Build TFBS baseline DataFrame
    tfbs_rows = []
    if model.tfbs_baseline and model.tfbs_baseline.baseline:
        for label, (mean, std) in model.tfbs_baseline.baseline.data.items():
            tfbs_rows.append({
                "table": "tfbs_baseline",
                "label": label,
                "entropy_mean": mean,
                "entropy_std": std,
            })
    tfbs_df = pd.DataFrame(tfbs_rows) if tfbs_rows else pd.DataFrame()
    
    # Build ATAC baseline DataFrame
    atac_rows = []
    if model.atac_baseline and model.atac_baseline.baseline:
        for label, (mean, std) in model.atac_baseline.baseline.data.items():
            atac_rows.append({
                "table": "atac_baseline",
                "label": label,
                "entropy_mean": mean,
                "entropy_std": std,
            })
    atac_df = pd.DataFrame(atac_rows) if atac_rows else pd.DataFrame()
    
    # =========================================================================
    # Build on-target TFBS baseline DataFrame (panel mode)
    # Separate table for on-target entropy which has different depth/variance
    # =========================================================================
    tfbs_ontarget_rows = []
    if model.tfbs_baseline_ontarget and model.tfbs_baseline_ontarget.baseline:
        for label, (mean, std) in model.tfbs_baseline_ontarget.baseline.data.items():
            tfbs_ontarget_rows.append({
                "table": "tfbs_baseline_ontarget",
                "label": label,
                "entropy_mean": mean,
                "entropy_std": std,
            })
    tfbs_ontarget_df = pd.DataFrame(tfbs_ontarget_rows) if tfbs_ontarget_rows else pd.DataFrame()
    
    # Build on-target ATAC baseline DataFrame (panel mode)
    atac_ontarget_rows = []
    if model.atac_baseline_ontarget and model.atac_baseline_ontarget.baseline:
        for label, (mean, std) in model.atac_baseline_ontarget.baseline.data.items():
            atac_ontarget_rows.append({
                "table": "atac_baseline_ontarget",
                "label": label,
                "entropy_mean": mean,
                "entropy_std": std,
            })
    atac_ontarget_df = pd.DataFrame(atac_ontarget_rows) if atac_ontarget_rows else pd.DataFrame()
    
    
    # Build FSC gene baseline DataFrame (panel mode)
    fsc_gene_rows = []
    if model.fsc_gene_baseline:
        for gene, (mean, std, n_samples) in model.fsc_gene_baseline.data.items():
            fsc_gene_rows.append({
                "table": "fsc_gene_baseline",
                "gene": gene,
                "depth_mean": mean,
                "depth_std": std,
                "n_samples": n_samples,
            })
    fsc_gene_df = pd.DataFrame(fsc_gene_rows) if fsc_gene_rows else pd.DataFrame()
    
    # Build FSC region baseline DataFrame (panel mode)
    fsc_region_rows = []
    if model.fsc_region_baseline:
        for region_id, (mean, std, n_samples) in model.fsc_region_baseline.data.items():
            fsc_region_rows.append({
                "table": "fsc_region_baseline",
                "region_id": region_id,
                "depth_mean": mean,
                "depth_std": std,
                "n_samples": n_samples,
            })
    fsc_region_df = pd.DataFrame(fsc_region_rows) if fsc_region_rows else pd.DataFrame()
    
    # Build Region MDS baseline DataFrame
    region_mds_rows = []
    if model.region_mds and model.region_mds.gene_baseline:
        for gene, data in model.region_mds.gene_baseline.items():
            region_mds_rows.append({
                "table": "region_mds",
                "gene": gene,
                "mds_mean": data.get("mds_mean", 0.0),
                "mds_std": data.get("mds_std", 1.0),
                "mds_e1_mean": data.get("mds_e1_mean"),
                "mds_e1_std": data.get("mds_e1_std"),
                "n_samples": data.get("n_samples", 0),
            })
    region_mds_df = pd.DataFrame(region_mds_rows) if region_mds_rows else pd.DataFrame()
    
    # Build WPS Background baseline DataFrame
    wps_background_df = pd.DataFrame()
    if model.wps_background_baseline and model.wps_background_baseline.groups is not None:
        wps_background_df = model.wps_background_baseline.groups.copy()
        wps_background_df["table"] = "wps_background"
    
    all_dfs = [metadata_df]
    if not gc_bias_df.empty:
        all_dfs.append(gc_bias_df)
    if not fsd_df.empty:
        all_dfs.append(fsd_df)
    if not wps_df.empty:
        all_dfs.append(wps_df)
    if not wps_panel_df.empty:
        all_dfs.append(wps_panel_df)
    if not ocf_df.empty:
        all_dfs.append(ocf_df)
    if not mds_df.empty:
        all_dfs.append(mds_df)
    if not fsd_ontarget_df.empty:
        all_dfs.append(fsd_ontarget_df)
    if not gc_bias_ontarget_df.empty:
        all_dfs.append(gc_bias_ontarget_df)
    if not tfbs_df.empty:
        all_dfs.append(tfbs_df)
    if not atac_df.empty:
        all_dfs.append(atac_df)
    # On-target entropy baselines (panel mode)
    if not tfbs_ontarget_df.empty:
        all_dfs.append(tfbs_ontarget_df)
    if not atac_ontarget_df.empty:
        all_dfs.append(atac_ontarget_df)
    # FSC gene/region baselines (panel mode)
    if not fsc_gene_df.empty:
        all_dfs.append(fsc_gene_df)
    if not fsc_region_df.empty:
        all_dfs.append(fsc_region_df)
    if not region_mds_df.empty:
        all_dfs.append(region_mds_df)
    if not wps_background_df.empty:
        all_dfs.append(wps_background_df)
    
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
    n_fsc_gene = len(fsc_gene_rows)
    n_fsc_region = len(fsc_region_rows)
    logger.info(f"Saved PON model: {output}")
    if model.panel_mode:
        logger.info(f"   Panel mode: ON (targets: {model.target_regions_file})")
        if n_fsd_on > 0:
            logger.info(f"   FSD on-target: {n_fsd_on} arm×size entries")
        if n_gc_on > 0:
            logger.info(f"   GC on-target: {n_gc_on} bins")
        if n_fsc_gene > 0:
            logger.info(f"   FSC gene: {n_fsc_gene} genes")
        if n_fsc_region > 0:
            logger.info(f"   FSC region: {n_fsc_region} regions")
    logger.info(f"   GC bias: {n_gc} bins")
    logger.info(f"   FSD baseline: {n_fsd} arm×size entries")
    logger.info(f"   WPS baseline: {n_wps} regions")
    logger.info(f"   OCF baseline: {n_ocf} regions")
    logger.info(f"   MDS baseline: {n_mds} k-mers")

