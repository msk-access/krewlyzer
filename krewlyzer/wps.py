import typer
from pathlib import Path
import logging
import pysam
import json

import numpy as np
from collections import defaultdict
import gzip
from rich.console import Console
from rich.logging import RichHandler

from .helpers import max_core, commonError

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("wps")

# Try to import Rust backend
try:
    import krewlyzer_core
    RUST_BACKEND_AVAILABLE = True
    logger.debug("Rust backend available for WPS")
except ImportError:
    RUST_BACKEND_AVAILABLE = False
    logger.debug("Rust backend not available")
    

import pandas as pd

def _calc_wps(
    bedgz_input: str | Path, 
    tsv_input: str | Path, 
    output_file_pattern: str, 
    empty: bool = False, 
    protect_input: int = 120, 
    min_size: int = 120, 
    max_size: int = 180
):
    """
    Calculate Windowed Protection Score (WPS) for a single .bed.gz file and transcript region file.
    Output is gzipped TSV per region.
    Optimized with vectorized operations.
    Automatically detects and uses metadata file for depth normalization if available.
    """
    try:
        # NEW: Try to load metadata for depth normalization
        total_fragments = None
        metadata_file = str(bedgz_input) + '.metadata.json'
        if Path(metadata_file).exists():
            try:
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
                    total_fragments = metadata.get('total_unique_fragments')
                    logger.info(f"Found metadata: {total_fragments:,} total fragments for WPS normalization")
            except Exception as e:
                logger.warning(f"Could not load metadata from {metadata_file}: {e}")
        else:
            logger.info("No metadata found - outputting raw WPS scores only")
        
        bedgzfile = str(bedgz_input)
        tbx = pysam.TabixFile(bedgzfile)
        protection = protect_input // 2
        
        logger.info(f"Processing {bedgz_input} with regions from {tsv_input}")
        
        with open(tsv_input, 'r') as infile:
            valid_chroms = set(map(str, list(range(1, 23)) + ["X", "Y"]))
            
            for line in infile:
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 5:
                    continue
                cid, chrom, start_str, end_str, strand = parts[:5]
                chrom = chrom.replace("chr", "")
                if chrom not in valid_chroms:
                    continue
                    
                region_start = int(float(start_str))
                region_end = int(float(end_str))
                
                if region_start < 1:
                    continue
                    
                # Region length (inclusive)
                length = region_end - region_start + 1
                
                # Arrays for the region (0-based index relative to region_start)
                cov_arr = np.zeros(length, dtype=int)
                start_arr = np.zeros(length, dtype=int)
                gcount_arr = np.zeros(length, dtype=int)
                total_arr = np.zeros(length, dtype=int)
                
                # Fetch reads
                # Pysam fetch is 0-based.
                # Region 1-based [S, E] -> 0-based [S-1, E)
                # We need reads extending 'protection' bp around the region
                fetch_start = max(0, region_start - protection - 1)
                fetch_end = region_end + protection
                
                # Ensure 'chr' prefix for fetching as per original logic
                fetch_chrom = "chr" + chrom if not chrom.startswith("chr") else chrom
                
                try:
                    rows = list(tbx.fetch(fetch_chrom, fetch_start, fetch_end, parser=pysam.asTuple()))
                except ValueError:
                    # Try without chr prefix if failed?
                    try:
                        rows = list(tbx.fetch(chrom, fetch_start, fetch_end, parser=pysam.asTuple()))
                    except ValueError:
                        rows = []
                except Exception as e:
                    logger.error(f"Error fetching region {chrom}:{region_start}-{region_end}: {e}")
                    continue
                
                if not rows:
                    if not empty:
                        continue
                
                # Process reads
                for row in rows:
                    # BED is 0-based start, 0-based exclusive end
                    rstart = int(row[1])
                    rend = int(row[2])
                    lseq = rend - rstart
                    
                    if lseq < min_size or lseq > max_size:
                        continue
                        
                    # Convert read to 1-based inclusive for easier logic with 1-based regions
                    r_start_1 = rstart + 1
                    r_end_1 = rend
                    
                    # 1. Coverage (covCount)
                    # Read spans [r_start_1, r_end_1]
                    ov_start = max(region_start, r_start_1)
                    ov_end = min(region_end, r_end_1)
                    
                    if ov_start <= ov_end:
                        idx_start = ov_start - region_start
                        idx_end = ov_end - region_start + 1
                        cov_arr[idx_start:idx_end] += 1
                        
                    # 2. Start Count (ends)
                    if region_start <= r_start_1 <= region_end:
                        start_arr[r_start_1 - region_start] += 1
                    if region_start <= r_end_1 <= region_end:
                        start_arr[r_end_1 - region_start] += 1
                        
                    # 3. WPS
                    # gcount (spanning): window [k-P, k+P] is inside read
                    # Range: [r_start_1 + P, r_end_1 - P]
                    g_start = r_start_1 + protection
                    g_end = r_end_1 - protection
                    
                    g_ov_start = max(region_start, g_start)
                    g_ov_end = min(region_end, g_end)
                    
                    if g_ov_start <= g_ov_end:
                        idx_start = g_ov_start - region_start
                        idx_end = g_ov_end - region_start + 1
                        gcount_arr[idx_start:idx_end] += 1
                        
                    # total (overlapping): window [k-P, k+P] overlaps read
                    # Range: [r_start_1 - P, r_end_1 + P]
                    t_start = r_start_1 - protection
                    t_end = r_end_1 + protection
                    
                    t_ov_start = max(region_start, t_start)
                    t_ov_end = min(region_end, t_end)
                    
                    if t_ov_start <= t_ov_end:
                        idx_start = t_ov_start - region_start
                        idx_end = t_ov_end - region_start + 1
                        total_arr[idx_start:idx_end] += 1
                
                # WPS = Spanning - (Total - Spanning)
                wps_arr = 2 * gcount_arr - total_arr
                
                # NEW: Add depth-normalized WPS if we have total fragments
                if total_fragments:
                    # Normalize by depth (per million fragments)
                    wps_norm_arr = wps_arr / (total_fragments / 1e6)
                
                # Check if we should write
                if np.sum(cov_arr) == 0 and not empty:
                    continue
                    
                # Prepare output
                filename = output_file_pattern % cid
                
                positions = np.arange(region_start, region_end + 1)
                df_data = {
                    'chrom': chrom,
                    'pos': positions,
                    'cov': cov_arr,
                    'starts': start_arr,
                    'wps': wps_arr
                }
                if total_fragments:
                    df_data['wps_norm'] = wps_norm_arr
                
                df = pd.DataFrame(df_data)
                
                if strand == "-":
                    df = df.iloc[::-1]
                
                with gzip.open(filename, 'wt') as outfile:
                    if total_fragments:
                        # NEW: Write with normalized column
                        for _, row in df.iterrows():
                            outfile.write(f"{row['chrom']}\t{int(row['pos'])}\t{int(row['cov'])}\t{int(row['starts'])}\t{int(row['wps'])}\t{row['wps_norm']}\n")
                    else:
                        # Original format
                        for _, row in df.iterrows():
                            outfile.write(f"{row['chrom']}\t{int(row['pos'])}\t{int(row['cov'])}\t{int(row['starts'])}\t{int(row['wps'])}\n")
                        
        logger.info(f"WPS calculation complete. Results written to pattern: {output_file_pattern}")
    except Exception as e:
        logger.error(f"Fatal error in _calc_wps: {e}")
        raise typer.Exit(1)


def _wps_task(bedgz_file, output_dir, tsv_input, empty, protect_input, min_size, max_size):
    """Module-level worker function for parallel WPS calculation."""
    import traceback
    try:
        output_file_pattern = str(Path(output_dir) / (Path(bedgz_file).stem.replace('.bed', '') + ".%s.WPS.tsv.gz"))
        _calc_wps(
            bedgz_input=str(bedgz_file),
            tsv_input=str(tsv_input),
            output_file_pattern=output_file_pattern,
            empty=empty,
            protect_input=protect_input,
            min_size=min_size,
            max_size=max_size
        )
        return None
    except Exception as exc:
        return traceback.format_exc()


def wps(
    bedgz_path: Path = typer.Argument(..., help="Folder containing .bed.gz files (should be the output directory from motif.py)"),
    tsv_input: Path = typer.Option(None, "--tsv-input", "-t", help="Path to transcript/region file (TSV format)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output folder for results"),
    wpstype: str = typer.Option('L', "--wpstype", "-w", help="WPS type: 'L' for long (default), 'S' for short"),
    empty: bool = typer.Option(False, "--empty", help="Keep files of empty blocks (default: False)"),
    threads: int = typer.Option(1, "--threads", "-p", help="Number of threads (default: 1)")
):
    """
    Calculate Windowed Protection Score (WPS) features for all .bed.gz files in a folder.
    """
    # Input checks
    if not bedgz_path.exists():
        logger.error(f"Input directory not found: {bedgz_path}")
        raise typer.Exit(1)
    if tsv_input and not tsv_input.exists():
        logger.error(f"Transcript region file not found: {tsv_input}")
        raise typer.Exit(1)
    try:
        output.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Could not create output directory {output}: {e}")
        raise typer.Exit(1)
    try:
        output.touch()
    except Exception as e:
        logger.error(f"Output directory {output} is not writable: {e}")
        raise typer.Exit(1)
    try:
        bedgz_files = list(Path(bedgz_path).glob("*.bed.gz"))
        if not bedgz_files:
            logger.error("No .bed.gz files found in the specified folder.")
            raise typer.Exit(1)
        if tsv_input is None:
            # Default to package data transcriptAnno-hg19-1kb.tsv
            tsv_input = Path(__file__).parent / "data" / "TranscriptAnno" / "transcriptAnno-hg19-1kb.tsv"
            logger.info(f"No tsv_input specified. Using default: {tsv_input}")
        if not tsv_input.exists():
            logger.error(f"Transcript/region file does not exist: {tsv_input}")
            raise typer.Exit(1)
        if wpstype == 'L':
            protect_input = 120
            min_size = 120
            max_size = 180
        else:
            protect_input = 16
            min_size = 35
            max_size = 80
        output.mkdir(parents=True, exist_ok=True)
        logger.info(f"Calculating WPS for {len(bedgz_files)} files...")
        if RUST_BACKEND_AVAILABLE:
            logger.info("Using Rust backend for accelerated WPS calculation")
            logger.info(f"Processing {len(bedgz_files)} files sequentially (parallelism handled internally by Rust)...")
            
            for bedgz_file in bedgz_files:
                try:
                    # Get metadata if available
                    total_fragments = None
                    try:
                        metadata_file = str(bedgz_file) + '.metadata.json'
                        if Path(metadata_file).exists():
                            with open(metadata_file, 'r') as f:
                                meta = json.load(f)
                                total_fragments = meta.get('total_unique_fragments')
                    except Exception as e:
                        logger.warning(f"Failed to load metadata: {e}")

                    bed_stem = bedgz_file.name.replace('.bed.gz', '').replace('.bed', '')
                    
                    # Call Rust
                    # calculate_wps(bedgz_path, tsv_path, output_dir, file_stem, empty, protect, min_size, max_size, total_fragments)
                    count = krewlyzer_core.calculate_wps(
                        str(bedgz_file), 
                        str(tsv_input), 
                        str(output), 
                        bed_stem,
                        empty,
                        protect_input,
                        min_size,
                        max_size,
                        total_fragments
                    )
                    logger.info(f"WPS calculated for {bedgz_file.name}: processed {count} regions")
                except Exception as e:
                    logger.error(f"Rust WPS calculation failed for {bedgz_file}: {e}")
                    # Fallback or continue? Continue for now.
        else:
            # Original Python parallel execution
            from concurrent.futures import ProcessPoolExecutor, as_completed
            from functools import partial
            
            n_procs = max_core(threads) if threads else 1
            logger.info(f"Calculating WPS for {len(bedgz_files)} files using {n_procs} processes...")
            
            worker = partial(_wps_task, output_dir=str(output), tsv_input=str(tsv_input),
                             empty=empty, protect_input=protect_input, min_size=min_size, max_size=max_size)
            
            with ProcessPoolExecutor(max_workers=n_procs) as executor:
                futures = {executor.submit(worker, str(bedgz_file)): bedgz_file for bedgz_file in bedgz_files}
                for future in as_completed(futures):
                    exc = future.result()
                    if exc:
                        logger.error(f"WPS calculation failed for {futures[future]}:\n{exc}")
        
        logger.info(f"WPS features calculated for {len(bedgz_files)} files.")
    except Exception as e:
        logger.error(f"Fatal error in wps CLI: {e}")
        raise typer.Exit(1)

