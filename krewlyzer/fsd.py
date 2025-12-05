import typer
from pathlib import Path
import logging
import pysam
import json

import numpy as np
from rich.console import Console
from rich.logging import RichHandler

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("fsd")

import pandas as pd

def _calc_fsd(
    bedgz_input: str | Path, 
    arms_file: str | Path, 
    output_file: str | Path
):
    """
    Internal: Calculate fragment size distribution (FSD) for a single .bed.gz file.
    Optimized with vectorized operations.
    Automatically detects and uses metadata file for CPM normalization if available.
    """
    try:
        # NEW: Try to load metadata for CPM normalization
        total_fragments = None
        metadata_file = str(bedgz_input) + '.metadata.json'
        if Path(metadata_file).exists():
            try:
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
                    total_fragments = metadata.get('total_unique_fragments')
                    logger.info(f"Found metadata: {total_fragments:,} total fragments")
            except Exception as e:
                logger.warning(f"Could not load metadata from {metadata_file}: {e}")
        else:
            logger.info("No metadata found - outputting region-normalized FSD only")
        
        logger.info(f"Processing {bedgz_input} with regions from {arms_file}")
        
        # Load regions
        try:
            bins_df = pd.read_csv(arms_file, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'], dtype={'chrom': str, 'start': int, 'end': int})
        except Exception as e:
            logger.error(f"Could not load regions from {arms_file}: {e}")
            raise typer.Exit(1)
            
        try:
            tbx = pysam.TabixFile(filename=bedgz_input, mode="r")
        except Exception as e:
            logger.error(f"Could not open {bedgz_input} as Tabix file: {e}")
            raise typer.Exit(1)

        results = []
        regions = []
        region_totals = []  # NEW: Track raw counts per region for CPM calculation
        
        # Define histogram bins
        hist_bins = np.arange(65, 401, 5) # 65, 70, ..., 400. 67 bins.
        
        for _, bin_row in bins_df.iterrows():
            chrom = bin_row['chrom']
            start = bin_row['start']
            end = bin_row['end']
            regions.append(f"{chrom}:{start}-{end}")
            
            try:
                rows = list(tbx.fetch(chrom, start, end, parser=pysam.asTuple()))
            except ValueError:
                rows = []
            except Exception as e:
                logger.error(f"Error fetching {chrom}:{start}-{end}: {e}")
                raise typer.Exit(1)
            
            if not rows:
                results.append(np.zeros(67))
                region_totals.append(0)  # NEW: Store zero for empty regions
                continue
                
            try:
                # Vectorized parsing
                _, starts, ends, _ = zip(*rows)
                starts = np.array(starts, dtype=int)
                ends = np.array(ends, dtype=int)
                lengths = ends - starts
                
                # Histogram
                counts, _ = np.histogram(lengths, bins=hist_bins)
                
                # Normalize (region-level)
                total = np.sum(counts)
                region_totals.append(total)  # NEW: Store for CPM calculation
                if total > 0:
                    results.append(counts / total)
                else:
                    results.append(np.zeros(67))
                    
            except Exception as e:
                logger.error(f"Error processing data in bin {chrom}:{start}-{end}: {e}")
                results.append(np.zeros(67))
                region_totals.append(0)  # NEW: Store zero for error cases
                continue

        # Write output
        try:
            with open(output_file, 'w') as f:
                # Header
                header_bins = [f"{s}-{s+4}" for s in hist_bins[:-1]]
                
                # NEW: Add CPM columns if we have total fragments
                if total_fragments:
                    header_cpm = [f"{s}-{s+4}_cpm" for s in hist_bins[:-1]]
                    f.write('region\t' + '\t'.join(header_bins) + '\t' + '\t'.join(header_cpm) + '\n')
                else:
                    f.write('region\t' + '\t'.join(header_bins) + '\n')
                
                # Data rows
                for i, region in enumerate(regions):
                    scores = results[i]  # Region-normalized [0-1]
                    row = f"{region}\t" + "\t".join(map(str, scores))
                    
                    # NEW: Add CPM columns if we have total fragments
                    if total_fragments:
                        # CPM = (region_normalized * region_total / total_fragments) * 1e6
                        region_total = region_totals[i]
                        if region_total > 0:
                            cpm_scores = (scores * region_total / total_fragments) * 1e6
                        else:
                            cpm_scores = np.zeros(67)
                        row += "\t" + "\t".join(map(str, cpm_scores))
                    
                    f.write(row + "\n")
                    
        except Exception as e:
            logger.error(f"Error writing FSD output file: {e}")
            raise typer.Exit(1)
            
        logger.info(f"FSD calculation complete. Results written to {output_file}")

    except Exception as e:
        logger.error(f"Fatal error in _calc_fsd: {e}")
        raise typer.Exit(1)

def fsd(
    bedgz_path: Path = typer.Argument(..., help="Folder containing .bed.gz files (should be the output directory from motif.py)"),
    arms_file: Path = typer.Option(..., "--arms-file", "-a", help="Path to arms/region file (BED format)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output folder for results"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads (default: 1)")
):
    """
    Calculate fragment size distribution (FSD) features for all .bed.gz files in a folder.
    The input folder should be the output directory produced by motif.py, containing the .bed.gz files.
    Output files are written to the output directory, one per .bed.gz file.
    """
    # Input checks
    if not bedgz_path.exists():
        logger.error(f"Input directory not found: {bedgz_path}")
        raise typer.Exit(1)
    if not arms_file.exists():
        logger.error(f"Arms/region file not found: {arms_file}")
        raise typer.Exit(1)
    try:
        output.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Could not create output directory {output}: {e}")
        raise typer.Exit(1)
    if not output.exists():
        logger.error(f"Output directory not found: {output}")
        raise typer.Exit(1)
    if not output.is_dir():
        logger.error(f"Output path is not a directory: {output}")
        raise typer.Exit(1)
    bedgz_files = [f for f in bedgz_path.iterdir() if f.suffixes == ['.bed', '.gz']]
    if not bedgz_files:
        logger.error("No .bed.gz files found in the specified folder.")
        raise typer.Exit(1)
    if not arms_file.exists():
        logger.error(f"Arms/region file does not exist: {arms_file}")
        raise typer.Exit(1)
    logger.info(f"Calculating FSD for {len(bedgz_files)} files...")
    from concurrent.futures import ProcessPoolExecutor, as_completed
    logger.info(f"Starting parallel FSD calculation using {threads} processes...")
    def run_fsd_file(bedgz_file):
        output_file = output / (bedgz_file.stem.replace('.bed', '') + '.FSD.txt')
        _calc_fsd(str(bedgz_file), str(arms_file), str(output_file))
        return str(output_file)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(run_fsd_file, bedgz_file): bedgz_file for bedgz_file in bedgz_files}
        for future in as_completed(futures):
            bedgz_file = futures[future]
            try:
                result = future.result()
                logger.info(f"FSD calculated: {result}")
            except Exception as exc:
                logger.error(f"FSD calculation failed for {bedgz_file}: {exc}")
    logger.info(f"FSD features calculated for {len(bedgz_files)} files.")
