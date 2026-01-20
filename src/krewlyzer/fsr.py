"""
Fragment Size Ratio (FSR) calculation.

Calculates FSR features for a single sample with comprehensive counts and ratios.
Uses Rust backend for accelerated computation with GC correction.

Fragment Size Bins (as defined by Rust backend):
    - ultra_short: 65-99bp (TF footprints, highly tumor-specific)
    - core_short: 100-149bp (tumor-enriched, primary biomarker)
    - mono_nucl: 150-259bp (mono-nucleosomal fragments)
    - di_nucl: 260-399bp (di-nucleosomal fragments)
    - long: 400+bp (di-nucleosomal and larger)
    
Output Columns:
    Counts: ultra_short_count, core_short_count, mono_nucl_count, di_nucl_count, long_count, total_count
    Ratios: ultra_short_ratio, core_short_ratio, mono_nucl_ratio, di_nucl_ratio, long_ratio
    Biomarker: core_short_long_ratio (primary cancer signal)
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
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("fsr")

# Rust backend is required
from krewlyzer import _core


def fsr(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for hybrid GC correction"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED (for panel data: generates on/off-target FSR)"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    format: Optional[str] = typer.Option(None, "--format", "-f", help="Output format override: tsv, parquet, json (default: tsv)"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging")
):
    """
    Calculate Fragment Size Ratio (FSR) features for a single sample.
    
    Computes fragment size distributions and ratios for cancer biomarker analysis.
    Uses Rust backend for performance with GC-corrected counts.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSR.tsv with:
        - Counts: ultra_short, core_short, mono_nucl, di_nucl, long, total
        - Ratios: per-bin fraction of total
        - Biomarker: core_short_long_ratio (tumor vs healthy signal)
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
        # Resolve GC correction assets (centralized helper)
        from .core.gc_assets import resolve_gc_assets
        gc = resolve_gc_assets(assets, output, sample_name, bedgz_input, gc_correct, genome)
        gc_ref = gc.gc_ref
        valid_regions = gc.valid_regions
        factors_out = gc.factors_out
        factors_input = gc.factors_input
        gc_correct = gc.gc_correct_enabled
        
        # Call Unified Pipeline (FSC bins mode - same bins as FSR)
        # FSR uses the same bin-level counts as FSC, just computes ratios instead of z-scores
        logger.info("Running fragment counting via unified pipeline...")
        
        is_panel_mode = target_regions and target_regions.exists()
        if is_panel_mode:
            logger.info(f"Panel mode: on/off-target split enabled (targets: {target_regions.name})")
        
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
            # WPS (disabled)
            None, None,
            # WPS Background (disabled)
            None, None, False,
            # FSD (disabled)
            None, None,
            # OCF (disabled)
            None, None,
            # Target regions for on/off-target split
            str(target_regions) if is_panel_mode else None,
            50,  # bait_padding
            False  # silent
        )
        
        logger.info(f"Reading counts from {fsr_counts_file}")
        
        # Load the counts TSV from Rust backend
        # Columns: chrom, start, end, ultra_short, core_short, mono_nucl, di_nucl, long, total, mean_gc
        df_counts = pd.read_csv(fsr_counts_file, sep='\t')
        
        logger.info(f"Loaded {len(df_counts)} bin counts")
        logger.debug(f"Columns present: {list(df_counts.columns)}")
        if 'total' in df_counts.columns:
            logger.debug(f"Total fragments: {df_counts['total'].sum():.0f}")
        
        # Validate required columns exist (Rust column names)
        # Rust outputs: ultra_short, core_short, mono_nucl, di_nucl, long, total
        required_cols = ['chrom', 'ultra_short', 'core_short', 'mono_nucl', 'di_nucl', 'long', 'total']
        missing_cols = [c for c in required_cols if c not in df_counts.columns]
        if missing_cols:
            logger.warning(f"FSR: Missing columns {missing_cols}, skipping processing")
            Path(output).mkdir(parents=True, exist_ok=True)
            pd.DataFrame().to_csv(Path(output) / f"{sample_name}.FSR.tsv", sep='\t', index=False)
            return
        
        # Check if there's any data
        if len(df_counts) == 0 or df_counts['total'].sum() == 0:
            logger.warning("FSR: No fragment data, writing empty output")
            Path(output).mkdir(parents=True, exist_ok=True)
            pd.DataFrame().to_csv(Path(output) / f"{sample_name}.FSR.tsv", sep='\t', index=False)
            return
        
        # Step 2: Apply PON correction if available (overlay on corrected counts)
        if pon and pon.gc_bias:
            logger.info("Applying PON-based GC normalization...")
            gc_arr = df_counts['mean_gc'].values
            core_short = df_counts['core_short'].values.astype(float)
            mono_nucl = df_counts['mono_nucl'].values.astype(float)
            di_nucl = df_counts['di_nucl'].values.astype(float)
            long_arr = df_counts['long'].values.astype(float)
            
            for i, gc in enumerate(gc_arr):
                core_exp = pon.gc_bias.get_expected(gc, "core_short")
                mono_exp = pon.gc_bias.get_expected(gc, "mono_nucl")
                di_exp = pon.gc_bias.get_expected(gc, "di_nucl")
                long_exp = pon.gc_bias.get_expected(gc, "long")
                
                if core_exp > 0:
                    core_short[i] /= core_exp
                if mono_exp > 0:
                    mono_nucl[i] /= mono_exp
                if di_exp > 0:
                    di_nucl[i] /= di_exp
                if long_exp > 0:
                    long_arr[i] /= long_exp
            
            df_counts['core_short'] = core_short
            df_counts['mono_nucl'] = mono_nucl
            df_counts['di_nucl'] = di_nucl
            df_counts['long'] = long_arr
            df_counts['total'] = core_short + mono_nucl + di_nucl + long_arr
        
        # Step 3: Aggregate bins into windows and compute ratios
        results = []
        for chrom, group in df_counts.groupby('chrom', sort=False):
            n_bins = len(group)
            n_windows = n_bins // continue_n
            
            if n_windows == 0:
                continue
            
            trunc_len = n_windows * continue_n
            
            ultra_mat = group['ultra_short'].values[:trunc_len].reshape(n_windows, continue_n)
            core_short_mat = group['core_short'].values[:trunc_len].reshape(n_windows, continue_n)
            mono_nucl_mat = group['mono_nucl'].values[:trunc_len].reshape(n_windows, continue_n)
            di_nucl_mat = group['di_nucl'].values[:trunc_len].reshape(n_windows, continue_n)
            long_mat = group['long'].values[:trunc_len].reshape(n_windows, continue_n)
            total_mat = group['total'].values[:trunc_len].reshape(n_windows, continue_n)
            
            ultra_sums = ultra_mat.sum(axis=1)
            core_short_sums = core_short_mat.sum(axis=1)
            mono_nucl_sums = mono_nucl_mat.sum(axis=1)
            di_nucl_sums = di_nucl_mat.sum(axis=1)
            long_sums = long_mat.sum(axis=1)
            total_sums = total_mat.sum(axis=1)
            
            window_starts = np.arange(n_windows) * continue_n * windows
            window_ends = (np.arange(n_windows) + 1) * continue_n * windows - 1
            
            for i in range(n_windows):
                region = f"{chrom}:{window_starts[i]}-{window_ends[i]}"
                
                u = float(ultra_sums[i])
                c = float(core_short_sums[i])
                m = float(mono_nucl_sums[i])
                d = float(di_nucl_sums[i])
                l = float(long_sums[i])
                t = float(total_sums[i])
                
                # Calculate ratios (avoid division by zero)
                if t > 0:
                    u_r = u / t
                    c_r = c / t
                    m_r = m / t
                    d_r = d / t
                    l_r = l / t
                else:
                    u_r = c_r = m_r = d_r = l_r = 0.0
                    
                # Core short / Long ratio (primary cancer biomarker)
                if l > 0:
                    cl_r = c / l
                else:
                    cl_r = c if c > 0 else 0.0
                
                results.append({
                    'region': region,
                    'ultra_short_count': int(u),
                    'core_short_count': int(c),
                    'mono_nucl_count': int(m),
                    'di_nucl_count': int(d),
                    'long_count': int(l),
                    'total_count': int(t),
                    'ultra_short_ratio': u_r,
                    'core_short_ratio': c_r,
                    'mono_nucl_ratio': m_r,
                    'di_nucl_ratio': d_r,
                    'long_ratio': l_r,
                    'core_short_long_ratio': cl_r,
                })
        
        if not results:
            logger.warning("No valid windows found for FSR.")
            raise typer.Exit(1)
        
        # Write output
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
        
        logger.info(f"FSR complete: {len(results_df)} windows → {output_file}")
        logger.info("Output columns: counts (ultra_short, core_short, mono_nucl, di_nucl, long, total) + ratios")

    except Exception as e:
        logger.error(f"FSR calculation failed: {e}")
        import traceback
        if verbose:
            traceback.print_exc()
        raise typer.Exit(1)
