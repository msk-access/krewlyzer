"""
PON model building from healthy plasma samples.

This module provides the CLI and logic for building unified PON models.
"""

import typer
from pathlib import Path
from typing import Optional, List
from datetime import datetime
import logging
import json

import numpy as np
import pandas as pd
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn

from .model import PonModel, GcBiasModel, FsdBaseline, WpsBaseline

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("build-pon")

# Import core tools for processing samples
from krewlyzer import _core


def build_pon(
    sample_list: Path = typer.Argument(..., help="Text file with paths to BED.gz files (one per line)"),
    assay: str = typer.Option(..., "--assay", "-a", help="Assay name (e.g., msk-access-v2)"),
    reference: Path = typer.Option(..., "--reference", "-r", help="Reference FASTA file"),
    transcript_file: Optional[Path] = typer.Option(None, "--transcript-file", "-t", help="Transcript TSV for WPS baseline"),
    bin_file: Optional[Path] = typer.Option(None, "--bin-file", "-b", help="Bin file for FSC/FSR (default: hg19_window_100kb.bed)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output PON model file (.pon.parquet)"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="BED file with target regions (panel mode - builds dual on/off-target baselines)"),
    threads: int = typer.Option(4, "--threads", "-p", help="Number of threads"),
    require_proper_pair: bool = typer.Option(False, "--require-proper-pair", help="Only extract properly paired reads (default: False for v1 ACCESS compatibility)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
):
    """
    Build a unified PON model from healthy plasma samples.
    
    The PON model contains:
    - GC bias curves for FSC/FSR/WPS correction
    - FSD baseline per chromosome arm
    - WPS baseline per transcript region
    
    For panel data (--target-regions):
    - Builds DUAL baselines for on-target and off-target regions
    - GC model trained on off-target only (unbiased by capture)
    - FSD/WPS include separate on-target stats
    
    Example:
        krewlyzer build-pon samples.txt --assay msk-access-v2 -r hg19.fa -o pon.parquet
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
    
    # Default bin file
    if bin_file is None:
        pkg_dir = Path(__file__).parent.parent
        bin_file = pkg_dir / "data" / "ChromosomeBins" / "GRCh37" / "hg19_window_100kb.bed.gz"
    
    if not bin_file.exists():
        logger.error(f"Bin file not found: {bin_file}")
        raise typer.Exit(1)
    
    # Default transcript file for WPS
    if transcript_file is None:
        pkg_dir = Path(__file__).parent.parent
        transcript_file = pkg_dir / "data" / "TranscriptAnno" / "transcriptAnno-hg19-1kb.tsv"
    
    # Handle BAM vs BED.gz input - convert all to BED.gz paths
    import tempfile
    import shutil
    temp_extract_dir = None
    bed_paths = []  # Final list of BED.gz paths to process
    
    for sample_path in samples:
        if sample_path.suffix == '.bam':
            # Extract BAM to temp BED.gz
            if temp_extract_dir is None:
                temp_extract_dir = tempfile.mkdtemp(prefix="pon_extract_")
                logger.info(f"Extracting BAM files to temp directory...")
            
            try:
                bed_output_path = str(Path(temp_extract_dir) / f"{sample_path.stem}.bed.gz")
                _core.extract_motif.process_bam_parallel(
                    str(sample_path),  # bam_path
                    str(reference),    # fasta_path
                    20,                # mapq
                    65,                # min_len
                    400,               # max_len
                    4,                 # kmer
                    threads,           # threads
                    bed_output_path,   # output_bed_path
                    None,              # output_motif_prefix (skip motifs)
                    None,              # exclude_path
                    target_regions_str,  # target_regions_path - for panel mode filtering
                    True,              # skip_duplicates
                    require_proper_pair,  # from CLI argument
                    True               # silent
                )
                bed_path = Path(bed_output_path)
                if bed_path.exists():
                    bed_paths.append(bed_path)
                    logger.debug(f"Extracted: {sample_path.name} -> {bed_path.name}")
                else:
                    logger.warning(f"Extraction produced no output for {sample_path.name}")
            except Exception as e:
                logger.warning(f"Failed to extract {sample_path.name}: {e}")
        else:
            # Already BED.gz
            bed_paths.append(sample_path)
    
    if not bed_paths:
        logger.error("No valid samples after processing")
        if temp_extract_dir:
            shutil.rmtree(temp_extract_dir)
        raise typer.Exit(1)
    
    logger.info(f"Processing {len(bed_paths)} BED.gz files...")
    
    # Collect data from all samples
    all_gc_data = []  # List of (gc_values, short, intermediate, long) tuples
    all_fsd_data = []  # List of DataFrames
    all_wps_data = []  # List of (region_id, wps_long_mean, wps_short_mean) tuples
    all_ocf_data = []  # List of OCF region stats
    all_mds_data = []  # List of MDS (kmer frequencies, mds score)
    
    # On-target data (panel mode only)
    all_fsd_data_ontarget = []  # List of on-target FSD DataFrames
    all_gc_data_ontarget = []   # List of on-target GC data (from ontarget correction factors)
    
    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task = progress.add_task(f"Processing {len(bed_paths)} samples...", total=len(bed_paths))
            
            for bed_path in bed_paths:
                if not bed_path.exists():
                    logger.warning(f"Sample not found, skipping: {bed_path}")
                    progress.advance(task)
                    continue
                
                sample_name = bed_path.name.replace(".bed.gz", "")
                progress.update(task, description=f"Processing {sample_name}...")
                
                try:
                    # Get FSC data (for GC curves) - this is off-target for WGS
                    _, short, intermediate, long, _, gc_values = _core.count_fragments_by_bins(
                        str(bed_path),
                        str(bin_file)
                    )
                    all_gc_data.append({
                        "gc": np.array(gc_values),
                        "short": np.array(short),
                        "intermediate": np.array(intermediate),
                        "long": np.array(long),
                    })
                    
                    # For panel mode, look for pre-processed on-target FSC files
                    # These would be in the same directory as the BED file
                    if is_panel_mode:
                        fsc_ontarget_file = bed_path.parent / f"{sample_name}.FSC.ontarget.tsv"
                        if fsc_ontarget_file.exists():
                            try:
                                fsc_on_df = pd.read_csv(fsc_ontarget_file, sep="\t")
                                # Extract GC-binned coverage from on-target FSC
                                if "gc" in fsc_on_df.columns:
                                    all_gc_data_ontarget.append({
                                        "gc": fsc_on_df["gc"].values if "gc" in fsc_on_df.columns else np.zeros(len(fsc_on_df)),
                                        "short": fsc_on_df["short"].values if "short" in fsc_on_df.columns else np.zeros(len(fsc_on_df)),
                                        "intermediate": fsc_on_df["intermediate"].values if "intermediate" in fsc_on_df.columns else np.zeros(len(fsc_on_df)),
                                        "long": fsc_on_df["long"].values if "long" in fsc_on_df.columns else np.zeros(len(fsc_on_df)),
                                    })
                                    logger.debug(f"Collected on-target FSC for {sample_name}")
                            except Exception as fsc_on_e:
                                logger.debug(f"Could not read on-target FSC: {fsc_on_e}")
                    
                    # Get FSD data per arm
                    try:
                        with tempfile.TemporaryDirectory() as tmpdir:
                            fsd_output = Path(tmpdir) / f"{sample_name}.FSD.tsv"
                            fsd_output_ontarget = Path(tmpdir) / f"{sample_name}.FSD.ontarget.tsv"
                            pkg_dir = Path(__file__).parent.parent
                            arms_file = pkg_dir / "data" / "ChromosomeArms" / "GRCh37" / "hg19.arms.bed.gz"
                            if arms_file.exists():
                                # Calculate FSD (Rust handles target_regions internally)
                                _core.fsd.calculate_fsd(
                                    str(bed_path),
                                    str(arms_file),
                                    str(fsd_output),
                                    target_regions_str  # Pass target regions for split output
                                )
                                # Collect off-target FSD
                                if fsd_output.exists():
                                    fsd_df = pd.read_csv(fsd_output, sep="\t")
                                    all_fsd_data.append(fsd_df)
                                # Collect on-target FSD (panel mode)
                                if is_panel_mode and fsd_output_ontarget.exists():
                                    fsd_on_df = pd.read_csv(fsd_output_ontarget, sep="\t")
                                    all_fsd_data_ontarget.append(fsd_on_df)
                                    logger.debug(f"Collected on-target FSD for {sample_name}")
                    except Exception as fsd_e:
                        logger.debug(f"FSD failed for {sample_name}: {fsd_e}")
                    
                    # Get WPS data per region (if transcript_file provided)
                    if transcript_file and transcript_file.exists():
                        try:
                            with tempfile.TemporaryDirectory() as tmpdir:
                                logger.debug(f"Calling WPS for {sample_name} with transcript: {transcript_file}")
                                _core.wps.calculate_wps(
                                    str(bed_path),
                                    str(transcript_file),
                                    str(tmpdir),
                                    sample_name,
                                    False,  # empty
                                    None,   # total_fragments
                                    str(reference) if reference else None,
                                    False,  # gc_correct
                                    False   # verbose
                                )
                                wps_file = Path(tmpdir) / f"{sample_name}.WPS.tsv.gz"
                                if wps_file.exists():
                                    wps_df = pd.read_csv(wps_file, sep="\t", compression="gzip")
                                    logger.debug(f"WPS output has {len(wps_df)} rows")
                                    region_stats = wps_df.groupby("gene_id").agg({
                                        "wps_long": "mean",
                                        "wps_short": "mean"
                                    }).reset_index()
                                    all_wps_data.append(region_stats)
                                else:
                                    logger.warning(f"WPS output file not created for {sample_name}")
                        except Exception as wps_e:
                            logger.warning(f"WPS failed for {sample_name}: {wps_e}")
                    
                    # Get OCF data per region (from open chromatin regions file)
                    try:
                        pkg_dir = Path(__file__).parent.parent
                        ocr_file = pkg_dir / "data" / "OCR" / "hg19_ocr_regions.bed.gz"
                        if ocr_file.exists():
                            with tempfile.TemporaryDirectory() as tmpdir:
                                ocf_output = Path(tmpdir) / f"{sample_name}.OCF.tsv"
                                _core.ocf.calculate_ocf(
                                    str(bed_path),
                                    str(ocr_file),
                                    str(ocf_output)
                                )
                                if ocf_output.exists():
                                    ocf_df = pd.read_csv(ocf_output, sep="\t")
                                    if "region_id" in ocf_df.columns and "ocf" in ocf_df.columns:
                                        all_ocf_data.append(ocf_df[["region_id", "ocf"]])
                    except Exception as ocf_e:
                        logger.debug(f"OCF failed for {sample_name}: {ocf_e}")
                    
                    # Get MDS (Motif Diversity Score) from extract output
                    try:
                        # MDS is computed from the motif output - look for .MDS.tsv
                        mds_file = bed_path.parent / f"{sample_name}.MDS.tsv"
                        if mds_file.exists():
                            mds_df = pd.read_csv(mds_file, sep="\t")
                            if "kmer" in mds_df.columns and "frequency" in mds_df.columns:
                                kmer_freqs = dict(zip(mds_df["kmer"], mds_df["frequency"]))
                                mds_score = mds_df["mds"].iloc[0] if "mds" in mds_df.columns else None
                                all_mds_data.append({"kmers": kmer_freqs, "mds": mds_score})
                    except Exception as mds_e:
                        logger.debug(f"MDS failed for {sample_name}: {mds_e}")
                    
                except Exception as e:
                    logger.warning(f"Failed to process {sample_name}: {e}")
                
                progress.advance(task)
    
    finally:
        # Clean up temp extraction directory
        if temp_extract_dir and Path(temp_extract_dir).exists():
            shutil.rmtree(temp_extract_dir)
            logger.debug(f"Cleaned up temp directory: {temp_extract_dir}")
    
    logger.info(f"Successfully processed {len(all_gc_data)} samples")
    
    if len(all_gc_data) < 1:
        logger.error("No samples processed successfully")
        raise typer.Exit(1)
    
    # Build GC bias model
    logger.info("Computing GC bias curves...")
    gc_bias = _compute_gc_bias_model(all_gc_data)
    
    # Build FSD baseline
    logger.info("Computing FSD baseline...")
    fsd_baseline = _compute_fsd_baseline(all_fsd_data)
    
    # Build WPS baseline
    logger.info("Computing WPS baseline...")
    wps_baseline = _compute_wps_baseline(all_wps_data)
    
    # Build OCF baseline
    ocf_baseline = None
    if all_ocf_data:
        logger.info("Computing OCF baseline...")
        ocf_baseline = _compute_ocf_baseline(all_ocf_data)
    
    # Build MDS baseline
    mds_baseline = None
    if all_mds_data:
        logger.info("Computing MDS baseline...")
        mds_baseline = _compute_mds_baseline(all_mds_data)
    
    # Build on-target FSD baseline (panel mode only)
    fsd_baseline_ontarget = None
    if is_panel_mode and all_fsd_data_ontarget:
        logger.info("Computing on-target FSD baseline...")
        fsd_baseline_ontarget = _compute_fsd_baseline(all_fsd_data_ontarget)
        if fsd_baseline_ontarget:
            logger.info(f"  On-target FSD: {len(fsd_baseline_ontarget.arms)} arms")
    
    # Build on-target GC bias (panel mode only)
    # Note: gc_bias_ontarget uses the same data as gc_bias since FSC already applies
    # target region filtering. For now, set to None - full implementation requires
    # separate on-target FSC data collection which is done in extract, not build-pon.
    gc_bias_ontarget = None
    if is_panel_mode and all_gc_data_ontarget:
        logger.info("Computing on-target GC bias curves...")
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
    logger.info(f"Saving PON model to {output}")
    _save_pon_model(model, output)
    
    logger.info(f"✅ PON model built successfully")
    logger.info(f"   Assay: {model.assay}")
    logger.info(f"   Samples: {model.n_samples}")
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
    region_id_col = "region_id" if "region_id" in template.columns else "group_id"
    
    # Columns to aggregate
    vector_cols = [c for c in ["wps_nuc", "wps_tf", "prot_frac_nuc", "prot_frac_tf"] 
                   if c in template.columns]
    
    if not vector_cols:
        logger.warning("No WPS vector columns found (expected wps_nuc, wps_tf)")
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
            # Collect vectors from all samples for this region
            vectors = []
            for df in all_wps_data:
                if idx < len(df) and col in df.columns:
                    vec = df.iloc[idx][col]
                    if vec is not None and len(vec) > 0:
                        vectors.append(np.array(vec, dtype=np.float32))
            
            if vectors:
                stacked = np.stack(vectors, axis=0)  # (n_samples, n_bins)
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

