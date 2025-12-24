"""
Fragment Size Ratio (FSR) calculation.

Calculates FSR features for a single sample with comprehensive counts and ratios.
Uses Rust backend for accelerated computation with GCfix correction.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

import numpy as np
import pandas as pd
from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("fsr")

# Rust backend is required
from krewlyzer import _core


def fsr(
    bedgz_input: Path = typer.Argument(..., help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for hybrid GC correction"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging")
):
    """
    Calculate Fragment Size Ratio (FSR) features for a single sample.
    
    Outputs comprehensive fragment counts and ratios for cancer biomarker analysis.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSR.tsv with columns:
        Counts: ultra_short_count, short_count, inter_count, long_count, total_count
        Ratios: short_ratio, inter_ratio, long_ratio, short_long_ratio, ultra_short_ratio
    
    Size categories (matching cfDNAFE):
        - Ultra-short: 65-100bp (TF footprints)
        - Short: 65-149bp (tumor-enriched)
        - Intermediate: 151-259bp (nucleosome dynamics)
        - Long: 261-399bp (healthy cell contribution)
        - Total: 65-399bp
    
    GC Correction:
        Uses GCfix LOESS-based correction via correction_factors.csv.
        Factors are automatically loaded from the extract step output.
        Use --no-gc-correct to disable.
    """
    from .assets import AssetManager
    
    # Configure Rust thread pool
    if threads > 0:
        try:
            _core.configure_threads(threads)
            logger.info(f"Configured {threads} threads for parallel processing")
        except Exception as e:
            logger.warning(f"Could not configure threads: {e}")
    
    # Input validation
    if not bedgz_input.exists():
        logger.error(f"Input file not found: {bedgz_input}")
        raise typer.Exit(1)
    
    if not str(bedgz_input).endswith('.bed.gz'):
        logger.error(f"Input must be a .bed.gz file: {bedgz_input}")
        raise typer.Exit(1)
    
    # Initialize Asset Manager
    try:
        assets = AssetManager(genome)
        logger.info(f"Genome: {assets.raw_genome} → {assets.genome_dir}")
    except ValueError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    
    # Resolve bin file via AssetManager
    if bin_input is None:
        try:
            bin_input = assets.resolve("bins_100kb")
            logger.info(f"Using default bin file: {bin_input}")
        except FileNotFoundError as e:
            logger.error(str(e))
            raise typer.Exit(1)
            
    if not bin_input.exists():
        logger.error(f"Bin file not found: {bin_input}")
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from input filename)
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    output_file = output / f"{sample_name}.FSR.tsv"
    fsr_counts_file = output / f"{sample_name}.fsr_counts.tsv"
    
    # Load PON model if provided
    pon = None
    if pon_model:
        from krewlyzer.pon.model import PonModel
        try:
            pon = PonModel.load(pon_model)
            logger.info(f"Loaded PON model: {pon.assay} (n={pon.n_samples})")
        except Exception as e:
            logger.warning(f"Could not load PON model: {e}")
            logger.warning("Falling back to within-sample correction only")

    try:
        logger.info(f"Processing {bedgz_input.name}")
        
        # Resolve GC correction assets
        gc_ref = None
        valid_regions = None
        factors_out = None
        
        if gc_correct:
            try:
                gc_ref = assets.resolve("gc_reference")
                valid_regions = assets.resolve("valid_regions")
                factors_out = output / f"{sample_name}.correction_factors.csv"
            except FileNotFoundError as e:
                logger.warning(f"GC correction assets not found: {e}")
                logger.warning("Proceeding without GC correction")
                gc_correct = False
        
        # Check for pre-computed correction factors (from extract step)
        factors_input = None
        if gc_correct:
            potential_factors = bedgz_input.parent / f"{bedgz_input.stem.replace('.bed', '')}.correction_factors.csv"
            if potential_factors.exists():
                factors_input = potential_factors
                logger.info(f"Using pre-computed correction factors: {factors_input}")
                # Skip computing new factors since we have pre-computed ones
                gc_ref = None
                valid_regions = None
                factors_out = None
        
        # Call Unified Pipeline (FSC bins mode - same bins as FSR)
        # FSR uses the same bin-level counts as FSC, just computes ratios instead of z-scores
        logger.info("Running fragment counting via unified pipeline...")
        _core.run_unified_pipeline(
            str(bedgz_input),
            # GC Correction (compute)
            str(gc_ref) if gc_ref else None,
            str(valid_regions) if valid_regions else None,
            str(factors_out) if factors_out else None,
            # GC Correction (load pre-computed)
            str(factors_input) if factors_input else None,
            # FSC (outputs count TSV that FSR also uses)
            str(bin_input), str(fsr_counts_file),
            # Others disabled
            None, None, False,  # WPS
            None, None,         # FSD
            None, None          # OCF
        )
        
        logger.info(f"Reading counts from {fsr_counts_file}")
        
        # Load the counts TSV
        # Format: chrom, start, end, ultra_short, short, intermediate, long, total, mean_gc
        df_counts = pd.read_csv(fsr_counts_file, sep='\t')
        
        logger.info(f"Loaded {len(df_counts)} bin counts")
        
        # Step 2: Apply PON correction if available (overlay on corrected counts)
        if pon and pon.gc_bias:
            logger.info("Applying PON-based GC normalization...")
            gc_arr = df_counts['mean_gc'].values
            short = df_counts['short'].values.astype(float)
            intermediate = df_counts['intermediate'].values.astype(float)
            long = df_counts['long'].values.astype(float)
            
            for i, gc in enumerate(gc_arr):
                short_exp = pon.gc_bias.get_expected(gc, "short")
                inter_exp = pon.gc_bias.get_expected(gc, "intermediate")
                long_exp = pon.gc_bias.get_expected(gc, "long")
                
                if short_exp > 0:
                    short[i] /= short_exp
                if inter_exp > 0:
                    intermediate[i] /= inter_exp
                if long_exp > 0:
                    long[i] /= long_exp
            
            df_counts['short'] = short
            df_counts['intermediate'] = intermediate
            df_counts['long'] = long
            df_counts['total'] = short + intermediate + long
        
        # Step 3: Aggregate bins into windows and compute ratios
        results = []
        for chrom, group in df_counts.groupby('chrom', sort=False):
            n_bins = len(group)
            n_windows = n_bins // continue_n
            
            if n_windows == 0:
                continue
            
            trunc_len = n_windows * continue_n
            
            ultra_mat = group['ultra_short'].values[:trunc_len].reshape(n_windows, continue_n)
            short_mat = group['short'].values[:trunc_len].reshape(n_windows, continue_n)
            inter_mat = group['intermediate'].values[:trunc_len].reshape(n_windows, continue_n)
            long_mat = group['long'].values[:trunc_len].reshape(n_windows, continue_n)
            total_mat = group['total'].values[:trunc_len].reshape(n_windows, continue_n)
            
            ultra_sums = ultra_mat.sum(axis=1)
            short_sums = short_mat.sum(axis=1)
            inter_sums = inter_mat.sum(axis=1)
            long_sums = long_mat.sum(axis=1)
            total_sums = total_mat.sum(axis=1)
            
            window_starts = np.arange(n_windows) * continue_n * windows
            window_ends = (np.arange(n_windows) + 1) * continue_n * windows - 1
            
            for i in range(n_windows):
                region = f"{chrom}:{window_starts[i]}-{window_ends[i]}"
                
                u = float(ultra_sums[i])
                s = float(short_sums[i])
                m = float(inter_sums[i])
                l = float(long_sums[i])
                t = float(total_sums[i])
                
                # Calculate ratios (avoid division by zero)
                if t > 0:
                    s_r = s / t
                    m_r = m / t
                    l_r = l / t
                    u_r = u / t
                else:
                    s_r = m_r = l_r = u_r = 0.0
                    
                # Short/Long ratio (primary cancer biomarker)
                if l > 0:
                    sl_r = s / l
                else:
                    sl_r = s if s > 0 else 0.0
                
                results.append({
                    'region': region,
                    'ultra_short_count': int(u),
                    'short_count': int(s),
                    'inter_count': int(m),
                    'long_count': int(l),
                    'total_count': int(t),
                    'short_ratio': s_r,
                    'inter_ratio': m_r,
                    'long_ratio': l_r,
                    'short_long_ratio': sl_r,
                    'ultra_short_ratio': u_r
                })
        
        if not results:
            logger.warning("No valid windows found for FSR.")
            raise typer.Exit(1)
        
        # Write output
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
        
        logger.info(f"FSR complete: {len(results_df)} windows → {output_file}")
        logger.info("Output columns: counts (ultra_short, short, inter, long, total) + ratios (short, inter, long, short_long, ultra_short)")

    except Exception as e:
        logger.error(f"FSR calculation failed: {e}")
        import traceback
        if verbose:
            traceback.print_exc()
        raise typer.Exit(1)
