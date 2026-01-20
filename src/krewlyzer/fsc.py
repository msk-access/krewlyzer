"""
Fragment Size Coverage (FSC) calculation.

Calculates FSC features for a single sample showing fragment coverage depth
across genomic windows, stratified by fragment size bin.

Uses Rust backend for accelerated computation with GC correction.

Fragment Size Bins (from Rust backend):
    - ultra_short: 65-99bp (TF footprints)
    - core_short: 100-149bp (tumor-enriched)
    - mono_nucl: 150-259bp (mono-nucleosomal)
    - di_nucl: 260-399bp (di-nucleosomal)
    - long: 400+bp

Output: {sample}.FSC.tsv with per-window coverage counts and z-scores.
Panel mode: Additional {sample}.FSC.ontarget.tsv for on-target regions.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

import pandas as pd
from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("fsc")

# Rust backend is required
from krewlyzer import _core


def fsc(
    bedgz_input: Path = typer.Option(..., "--input", "-i", help="Input .bed.gz file (output from extract)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory"),
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output file"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Path to bin file"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for hybrid GC correction"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED (for panel data: generates on/off-target FSC)"),
    assay: Optional[str] = typer.Option(None, "--assay", "-A", help="Assay type (xs1/xs2) for gene-centric FSC aggregation"),
    windows: int = typer.Option(100000, "--windows", "-w", help="Window size (default: 100000)"),
    continue_n: int = typer.Option(50, "--continue-n", "-c", help="Consecutive window number"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    gc_correct: bool = typer.Option(True, "--gc-correct/--no-gc-correct", help="Apply GC bias correction"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    format: Optional[str] = typer.Option(None, "--format", "-f", help="Output format override: tsv, parquet, json (default: tsv)")
):
    """
    Calculate fragment size coverage (FSC) features for a single sample.
    
    Input: .bed.gz file from extract step
    Output: {sample}.FSC.tsv file with z-scored fragment size coverage per window
    """
    from .assets import AssetManager
    from .core.fsc_processor import process_fsc
    from .core.pon_integration import load_pon_model
    
    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("core.fsc_processor").setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")
    
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
    
    logger.debug(f"Input: {bedgz_input}")
    logger.debug(f"Output directory: {output}")
    
    # Initialize Asset Manager
    try:
        assets = AssetManager(genome)
        logger.info(f"Genome: {assets.raw_genome} → {assets.genome_dir}")
    except ValueError as e:
        logger.error(str(e))
        raise typer.Exit(1)
        
    # Resolve bin file
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
    
    logger.debug(f"Bin file: {bin_input}")
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name
    if sample_name is None:
        sample_name = bedgz_input.name.replace('.bed.gz', '').replace('.bed', '')
    
    output_file = output / f"{sample_name}.FSC.tsv"
    fsc_counts_file = output / f"{sample_name}.fsc_counts.tsv"
    
    logger.debug(f"Sample name: {sample_name}")
    logger.debug(f"Output file: {output_file}")
    
    # Load PON model if provided
    pon = None
    if pon_model:
        pon = load_pon_model(pon_model)
        if pon is None:
            logger.warning("Falling back to within-sample correction only")

    try:
        logger.info(f"Processing {bedgz_input.name}")
        # Resolve GC correction assets (centralized helper)
        from .core.gc_assets import resolve_gc_assets
        gc = resolve_gc_assets(assets, output, sample_name, bedgz_input, gc_correct, genome)
        gc_ref = gc.gc_ref
        valid_regions = gc.valid_regions
        factors_out = gc.factors_out
        factors_input = gc.factors_input
        gc_correct = gc.gc_correct_enabled
        
        # Call Unified Pipeline (FSC only)
        logger.info("Running fragment counting via Rust backend...")
        
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
            # FSC
            str(bin_input), str(fsc_counts_file),
            # WPS (disabled)
            None, None,  # wps_regions, wps_output
            # WPS Background (disabled)
            None, None, False,  # wps_background
            # FSD (disabled)
            None, None,
            # OCF (disabled)
            None, None,
            # Target regions for on/off-target split
            str(target_regions) if is_panel_mode else None,
            50,  # bait_padding
            False  # silent
        )
        
        logger.info(f"Reading counts from {fsc_counts_file}")
        
        # Load the FSC counts TSV
        df_counts = pd.read_csv(fsc_counts_file, sep='\t')
        logger.debug(f"Loaded {len(df_counts)} bins, columns: {list(df_counts.columns)}")
        
        # Use shared processor for aggregation and z-score calculation
        process_fsc(
            counts_df=df_counts,
            output_path=output_file,
            windows=windows,
            continue_n=continue_n,
            pon=pon
        )
        
        # Process on-target FSC counts if panel mode
        fsc_counts_ontarget = output / f"{sample_name}.fsc_counts.ontarget.tsv"
        if is_panel_mode and fsc_counts_ontarget.exists():
            df_counts_on = pd.read_csv(fsc_counts_ontarget, sep='\t')
            output_file_on = output / f"{sample_name}.FSC.ontarget.tsv"
            process_fsc(
                counts_df=df_counts_on,
                output_path=output_file_on,
                windows=windows,
                continue_n=continue_n,
                pon=pon
            )
            logger.info(f"✅ FSC on-target: {output_file_on}")
        
        # Gene-centric FSC aggregation if assay is specified
        if assay:
            from .core.gene_bed import load_gene_bed
            from .core.fsc_processor import aggregate_by_gene
            
            # Map genome to GRCh format for gene_bed
            genome_map = {'hg19': 'GRCh37', 'grch37': 'GRCh37', 'hg38': 'GRCh38', 'grch38': 'GRCh38'}
            gene_genome = genome_map.get(genome.lower(), 'GRCh37')
            
            try:
                genes = load_gene_bed(assay=assay, genome=gene_genome)
                gene_fsc_file = output / f"{sample_name}.FSC.gene.tsv"
                
                aggregate_by_gene(
                    fragments_bed=bedgz_input,
                    genes=genes,
                    output_path=gene_fsc_file,
                    pon=pon
                )
                logger.info(f"✅ Gene FSC: {gene_fsc_file} ({len(genes)} genes)")
            except Exception as e:
                logger.warning(f"Gene FSC aggregation failed: {e}")
                if verbose:
                    import traceback
                    traceback.print_exc()
        
        logger.info(f"✅ FSC complete: {output_file}")

    except Exception as e:
        logger.error(f"FSC calculation failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
