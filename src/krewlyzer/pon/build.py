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
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
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
    threads: int = typer.Option(4, "--threads", "-p", help="Number of threads"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
):
    """
    Build a unified PON model from healthy plasma samples.
    
    The PON model contains:
    - GC bias curves for FSC/FSR/WPS correction
    - FSD baseline per chromosome arm
    - WPS baseline per transcript region
    
    Example:
        krewlyzer build-pon samples.txt --assay msk-access-v2 -r hg19.fa -o pon.parquet
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
    
    logger.info(f"Building PON from {n_samples} samples")
    logger.info(f"Assay: {assay}")
    logger.info(f"Reference: {reference}")
    
    # Default bin file
    if bin_file is None:
        pkg_dir = Path(__file__).parent.parent
        bin_file = pkg_dir / "data" / "ChormosomeBins" / "hg19_window_100kb.bed"
    
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
                    True,              # skip_duplicates
                    True               # require_proper_pair
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
                    # Get FSC data (for GC curves)
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
                    
                    # Get FSD data per arm
                    try:
                        with tempfile.TemporaryDirectory() as tmpdir:
                            fsd_output = Path(tmpdir) / f"{sample_name}.FSD.tsv"
                            pkg_dir = Path(__file__).parent.parent
                            arms_file = pkg_dir / "data" / "ChormosomeArms" / "hg19.arms.bed"
                            if arms_file.exists():
                                _core.fsd.calculate_fsd(
                                    str(bed_path),
                                    str(arms_file),
                                    str(fsd_output)
                                )
                                if fsd_output.exists():
                                    fsd_df = pd.read_csv(fsd_output, sep="\t")
                                    all_fsd_data.append(fsd_df)
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
    
    # Create PON model
    model = PonModel(
        schema_version="1.0",
        assay=assay,
        build_date=datetime.now().isoformat()[:10],
        n_samples=len(all_gc_data),
        reference=reference.name.replace(".fa", "").replace(".fasta", ""),
        gc_bias=gc_bias,
        fsd_baseline=fsd_baseline,
        wps_baseline=wps_baseline,
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
    
    For each GC bin, compute median expected coverage across samples.
    """
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


def _compute_fsd_baseline(all_fsd_data: List[pd.DataFrame]) -> Optional[FsdBaseline]:
    """
    Compute FSD baseline from sample data.
    
    FSD output format: DataFrame with 'region' column and size bin columns (e.g., '65-69', '70-74', ...)
    
    Aggregates size distribution proportions per chromosome arm across samples.
    """
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


def _compute_wps_baseline(all_wps_data: List[pd.DataFrame]) -> Optional[WpsBaseline]:
    """
    Compute WPS baseline from sample data.
    
    Aggregates mean WPS per region across samples.
    """
    if not all_wps_data:
        return None
    
    # Combine all samples
    combined = pd.concat(all_wps_data, ignore_index=True)
    
    # Aggregate per region
    baseline = combined.groupby("gene_id").agg({
        "wps_long": ["mean", "std"],
        "wps_short": ["mean", "std"]
    }).reset_index()
    
    # Flatten column names
    baseline.columns = ["region_id", "wps_long_mean", "wps_long_std", "wps_short_mean", "wps_short_std"]
    
    return WpsBaseline(regions=baseline)


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
    
    # Combine and save
    all_dfs = [metadata_df]
    if not gc_bias_df.empty:
        all_dfs.append(gc_bias_df)
    if not fsd_df.empty:
        all_dfs.append(fsd_df)
    if not wps_df.empty:
        all_dfs.append(wps_df)
    
    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_df.to_parquet(output, index=False)
    
    # Log summary
    n_gc = len(gc_bias_rows)
    n_fsd = len(fsd_rows)
    n_wps = len(wps_df) if not wps_df.empty else 0
    logger.info(f"Saved PON model: {output}")
    logger.info(f"   GC bias: {n_gc} bins")
    logger.info(f"   FSD baseline: {n_fsd} arm×size entries")
    logger.info(f"   WPS baseline: {n_wps} regions")

