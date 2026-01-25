"""
Run all krewlyzer feature extraction tools for a single sample.

Pipeline Overview:
    1. **Extract + Motif**: Fragment extraction, GC correction, end/breakpoint motifs
    2. **Unified Engine**: FSC, FSR, FSD, WPS, OCF (parallel via Rust backend)
    3. **Post-processing**: WPS smoothing (Rust), periodicity extraction (Rust), PoN z-scores
    4. **Optional**: UXM (bisulfite), mFSD (variants)

Key Features:
    - **BAM filter auto-detection**: Warns if duplex BAMs need --no-require-proper-pair
    - **Panel mode**: Off-target GC model, on/off-target split for all tools
    - **GC correction**: Per-fragment LOESS-based correction via Rust
    - **WPS periodicity**: FFT-based NRL extraction with deviation scoring (Rust)

Uses shared utilities from core/bam_utils.py for filter compatibility checking.
"""

import typer
from pathlib import Path
import logging
from typing import Optional

from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn


# Note: FSC/FSR processing is now handled by unified_processor, not called directly from here
from .assets import AssetManager

# Rust backend
from krewlyzer import _core

# Initialize logging globally, but level will be set in run_all
console = Console(stderr=True)
logging.basicConfig(
    level="INFO", 
    handlers=[RichHandler(console=console, show_time=True, show_path=False)], 
    format="%(message)s"
)
logger = logging.getLogger("krewlyzer")

# Import shared resolver functions
from .core.utils import resolve_path as _resolve_path
from .core.utils import resolve_int as _resolve_int
from .core.utils import resolve_str as _resolve_str

def run_all(
    bam_input: Path = typer.Option(..., "--input", "-i", help="Input BAM file (sorted, indexed)"),
    reference: Path = typer.Option(..., "--reference", "-r", help="Reference genome FASTA (indexed)"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory for all results"),
    
    # Configurable filters (exposed from filters.rs)
    mapq: int = typer.Option(20, "--mapq", "-q", help="Minimum mapping quality"),
    minlen: int = typer.Option(65, "--minlen", help="Minimum fragment length"),
    maxlen: int = typer.Option(400, "--maxlen", help="Maximum fragment length"),
    skip_duplicates: bool = typer.Option(True, "--skip-duplicates/--no-skip-duplicates", help="Skip duplicate reads"),
    require_proper_pair: bool = typer.Option(True, "--require-proper-pair/--no-require-proper-pair", help="Require proper pairs"),
    exclude_regions: Optional[Path] = typer.Option(None, "--exclude-regions", "-x", help="Exclude regions BED file"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED (for panel data: GC model from off-target reads only)"),
    assay: Optional[str] = typer.Option(None, "--assay", "-A", help="Assay code (xs1, xs2) for panel-specific assets and gene-centric FSC"),
    
    # Optional inputs for specific tools
    bisulfite_bam: Optional[Path] = typer.Option(None, "--bisulfite-bam", help="Bisulfite BAM for UXM (optional)"),
    variants: Optional[Path] = typer.Option(None, "--variants", "-v", help="VCF/MAF file for mFSD (optional)"),
    duplex: bool = typer.Option(False, "--duplex", "-D", help="Enable duplex weighting for mFSD (fgbio cD tag or Marianas read names)"),
    
    # Other options
    sample_name: Optional[str] = typer.Option(None, "--sample-name", "-s", help="Sample name for output files (default: derived from BAM filename)"),
    chromosomes: Optional[str] = typer.Option(None, "--chromosomes", help="Comma-separated chromosomes to process"),
    threads: int = typer.Option(0, "--threads", "-t", help="Number of threads (0=all cores)"),
    
    # Optional overrides
    arms_file: Optional[Path] = typer.Option(None, "--arms-file", "-a", help="Custom arms file for FSD"),
    bin_input: Optional[Path] = typer.Option(None, "--bin-input", "-b", help="Custom bin file for FSC/FSR"),
    ocr_file: Optional[Path] = typer.Option(None, "--ocr-file", help="Custom OCR file for OCF"),
    wps_file: Optional[Path] = typer.Option(None, "--wps-file", help="Custom transcript file for WPS (legacy)"),
    wps_anchors: Optional[Path] = typer.Option(None, "--wps-anchors", help="WPS anchors BED (merged TSS+CTCF) for dual-stream profiling"),
    wps_background: Optional[Path] = typer.Option(None, "--wps-background", help="WPS background Alu BED for hierarchical stacking (auto-loaded if not specified)"),
    bait_padding: int = typer.Option(50, "--bait-padding", help="Bait edge padding in bp (default 50, use 15-20 for small exon panels)"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for GC correction and z-scores"),
    skip_pon: bool = typer.Option(False, "--skip-pon", help="Skip PON z-score normalization for all tools (for PON samples used as ML negatives)"),
    
    # Output format options
    output_format: str = typer.Option("auto", "--output-format", "-F", help="Output format: auto (smart defaults), tsv, parquet, json"),
    generate_json: bool = typer.Option(False, "--generate-json", help="Generate unified sample.features.json with ALL data for ML pipelines"),
    
    # Observability
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging"),
):
    """
    Run all feature extraction tools for a single sample.
    
    Pipeline: extract → motif → [Unified Engine: FSC/FSR, FSD, WPS, OCF]
    Optional: uxm, mfsd
    """
    # ═══════════════════════════════════════════════════════════════════
    # RESOLVE TYPER PARAMETERS (handle direct function calls vs CLI)
    # ═══════════════════════════════════════════════════════════════════
    # When called directly (not via CLI), typer.Option defaults remain as
    # OptionInfo objects. Resolve them to their intended values.
    resolved_target_regions = _resolve_path(target_regions)
    resolved_exclude_regions = _resolve_path(exclude_regions)
    resolved_bisulfite_bam = _resolve_path(bisulfite_bam)
    resolved_variants = _resolve_path(variants)
    resolved_arms_file = _resolve_path(arms_file)
    resolved_ocr_file = _resolve_path(ocr_file)
    resolved_pon_model = _resolve_path(pon_model)
    resolved_bait_padding = _resolve_int(bait_padding, 50)
    resolved_assay = _resolve_str(assay)
    
    # Configure Logging
    if debug:
        logger.setLevel(logging.DEBUG)
        logging.getLogger().setLevel(logging.DEBUG)
        logger.debug("Debug logging enabled")
    else:
        logger.setLevel(logging.INFO)
        logging.getLogger().setLevel(logging.INFO)

    # Configure threads
    if threads > 0:
        try:
            _core.configure_threads(threads)
            logger.info(f"Configured {threads} threads")
        except Exception as e:
            logger.warning(f"Could not configure threads: {e}")
    
    # Input validation
    if not bam_input.exists():
        logger.error(f"BAM file not found: {bam_input}")
        raise typer.Exit(1)
    
    if not reference.exists():
        logger.error(f"Reference not found: {reference}")
        raise typer.Exit(1)
    
    # Derive sample name (use provided or derive from BAM filename)
    if sample_name is None:
        sample = bam_input.stem.replace('.bam', '')
    else:
        sample = sample_name
        
    # Initialize Asset Manager
    try:
        assets = AssetManager(genome)
        logger.info(f"Genome: {assets.raw_genome} -> {assets.genome_dir}")
    except ValueError as e:
        logger.error(str(e))
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Panel mode detection (affects FSC/FSR aggregation)
    is_panel_mode = resolved_target_regions and resolved_target_regions.exists()
    
    # === Centralized PON Model Loading ===
    # Load once here and pass to all processors (FSC, FSR, FSD, WPS, OCF, TFBS, ATAC)
    # Priority: 1. Explicit -P flag, 2. Bundled PON for assay, 3. None
    
    # Validate: -P and --skip-pon are contradictory
    if resolved_pon_model and skip_pon:
        console.print("[bold red]❌ ERROR:[/bold red] --pon-model (-P) and --skip-pon are contradictory.")
        console.print("  • Use -P to apply PON z-score normalization")
        console.print("  • Use --skip-pon with -A assay to auto-load PON but skip z-scores")
        raise typer.Exit(1)
    pon = None
    actual_pon_path = resolved_pon_model
    
    # Auto-load bundled PON when -A assay is specified and -P is not provided
    if not resolved_pon_model and resolved_assay:
        try:
            bundled_pon = assets.get_pon(resolved_assay)
            if bundled_pon and bundled_pon.exists():
                actual_pon_path = bundled_pon
                logger.info(f"Auto-loaded bundled PON for assay '{resolved_assay}'")
        except FileNotFoundError:
            logger.debug(f"No bundled PON found for assay '{resolved_assay}'")
    
    if actual_pon_path:
        from .core.pon_integration import load_pon_model
        pon = load_pon_model(actual_pon_path)
        if pon:
            logger.info(f"PON model loaded: {pon.assay} (n={pon.n_samples})")
            # Log component availability
            components = []
            if pon.gc_bias:
                components.append("GC-bias")
            if pon.fsd_baseline:
                components.append("FSD")
            if pon.wps_baseline:
                components.append("WPS")
            if pon.ocf_baseline:
                components.append("OCF")
            if pon.mds_baseline:
                components.append("MDS")
            logger.info(f"  Components: {', '.join(components) if components else 'none'}")
            
            # Panel compatibility check
            if hasattr(pon, 'check_panel_compatibility'):
                pon.check_panel_compatibility(is_panel_mode)
            
            # PON/sample panel mode validation
            resolved_assay = assay if isinstance(assay, str) else None
            from .pon.validation import validate_panel_config, print_validation_warnings
            validation = validate_panel_config(
                pon_model=pon,
                target_regions=resolved_target_regions,
                assay=resolved_assay
            )
            if validation.warnings or validation.errors:
                print_validation_warnings(validation, console)
            if not validation.valid:
                logger.error("PON validation failed. Use compatible PON or correct assay flag.")
                raise typer.Exit(1)
        else:
            logger.warning(f"Failed to load PON model: {resolved_pon_model}")
    
    # Window settings for FSC/FSR
    # - Custom bins or panel data: no aggregation (preserve gene-level resolution)
    # - WGS default: aggregate 50 bins → 5Mb windows for arm-level CNV
    if bin_input or is_panel_mode:
        fsc_windows, fsc_continue_n = 1, 1
        if bin_input:
            logger.info(f"Using custom bin file: {bin_input}")
        if is_panel_mode:
            logger.info(f"Panel mode: FSC/FSR aggregation disabled (preserving gene-level resolution)")
    else:
        fsc_windows, fsc_continue_n = 100000, 50
    
    logger.info(f"Processing sample: {sample}")
    logger.info(f"Filters: mapq>={mapq}, length=[{minlen},{maxlen}], skip_dup={skip_duplicates}, proper_pair={require_proper_pair}")
    
    # Pre-check BAM compatibility with current filters (avoid 0 fragment issue)
    from .core.bam_utils import check_bam_compatibility
    compat = check_bam_compatibility(bam_input, require_proper_pair, skip_duplicates, mapq)
    
    if compat["pass_rate"] < 0.01 and compat["total_sampled"] > 100:
        console.print("\n[bold red]⚠️  FILTER COMPATIBILITY WARNING[/bold red]\n")
        console.print(f"Only [bold]{compat['pass_rate']:.2%}[/bold] of sampled reads would pass current filters.\n")
        
        for issue in compat["issues"]:
            console.print(f"  • {issue}")
        
        if compat["suggested_flags"]:
            suggested = " ".join(compat["suggested_flags"])
            console.print(f"\n[bold yellow]Suggested command:[/bold yellow]")
            console.print(f"  krewlyzer run-all -i {bam_input.name} -r {reference.name} -o {output} [bold]{suggested}[/bold]\n")
            
            logger.error(f"Filter mismatch detected. Re-run with: {suggested}")
            raise typer.Exit(1)
    elif compat["pass_rate"] < 0.5 and compat["issues"]:
        # Warning but continue
        console.print("\n[yellow]⚠️  Filter compatibility note:[/yellow]")
        for issue in compat["issues"]:
            console.print(f"  • {issue}")
        console.print()
    
    if is_panel_mode:
        logger.info(f"Panel mode: GC model will use off-target reads only (targets: {resolved_target_regions.name})")
    
    # Determine which optional tools will run
    has_uxm = resolved_bisulfite_bam is not None and resolved_bisulfite_bam.exists()
    has_mfsd = resolved_variants is not None and resolved_variants.exists()
    
    # Multi-step progress display
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
        transient=False,  # Keep progress visible after completion
    ) as progress:
        # Create all task placeholders
        task_extract = progress.add_task("[1/5] Extract + Motif", total=100)
        task_pipeline = progress.add_task("[2/5] Unified Pipeline", total=100, start=False)
        task_postproc = progress.add_task("[3/5] Post-processing", total=100, start=False)
        task_uxm = progress.add_task(f"[4/5] UXM", total=100, start=False, visible=has_uxm)
        task_mfsd = progress.add_task(f"[5/5] mFSD", total=100, start=False, visible=has_mfsd)
        
        # ═══════════════════════════════════════════════════════════════════
        # 1. EXTRACT + MOTIF (Using unified extract_sample)
        # ═══════════════════════════════════════════════════════════════════
        from .core.sample_processor import extract_sample, write_motif_outputs, write_extraction_outputs, ExtractionResult
    
        bedgz_file = output / f"{sample}.bed.gz"
        
        # Motif Outputs (for caching check)
        edm_output = output / f"{sample}.EndMotif.tsv"
        bpm_output = output / f"{sample}.BreakPointMotif.tsv"
        mds_output = output / f"{sample}.MDS.tsv"
        
        # Caching logic: determine what needs to run
        should_run_extract = not bedgz_file.exists()
        should_run_engine = (
            not bedgz_file.exists() or 
            not edm_output.exists() or 
            not bpm_output.exists() or 
            not mds_output.exists()
        )
        
        extraction_result: ExtractionResult = None
        
        if not should_run_engine:
            progress.update(task_extract, description="[1/5] Extract + Motif  ✓ (cached)", completed=100)
        else:
            progress.update(task_extract, description="[1/5] Extract + Motif  ...")
            
            try:
                # Use unified extract_sample() for all extraction
                extraction_result = extract_sample(
                    input_path=bam_input,
                    reference=reference,
                    mapq=mapq,
                    minlen=minlen,
                    maxlen=maxlen,
                    kmer=4,
                    threads=threads,
                    skip_duplicates=skip_duplicates,
                    require_proper_pair=require_proper_pair,
                    target_regions=resolved_target_regions,
                    exclude_regions=resolved_exclude_regions or assets.exclude_regions,
                    write_bed=should_run_extract,
                    output_dir=output if should_run_extract else None,
                    sample_name=sample,
                    genome=genome,
                )
                
                # Write extraction outputs (BED index, metadata, GC factors)
                if should_run_extract and extraction_result.bed_path:
                    resolved_assay = assay if isinstance(assay, str) else None
                    write_extraction_outputs(
                        result=extraction_result,
                        output_dir=output,
                        genome=genome,
                        assay=resolved_assay,
                        compute_gc_factors=True,
                    )
                    bedgz_file = extraction_result.bed_path
                    
                    # Compute on-target GC factors for panel mode
                    # These are used by mFSD and gene-level FSC for accurate panel analysis
                    if is_panel_mode and len(extraction_result.gc_observations_ontarget) > 0:
                        try:
                            gc_ref = assets.resolve("gc_reference")
                            valid_regions = assets.resolve("valid_regions")
                            factors_ontarget = output / f"{sample}.correction_factors.ontarget.tsv"
                            logger.info(f"Computing on-target GC factors from {len(extraction_result.gc_observations_ontarget)} obs...")
                            n_factors_on = _core.gc.compute_and_write_gc_factors(
                                extraction_result.gc_observations_ontarget,
                                str(gc_ref),
                                str(valid_regions),
                                str(factors_ontarget)
                            )
                            logger.info(f"Wrote {n_factors_on} on-target GC factors: {factors_ontarget.name}")
                        except Exception as e:
                            logger.warning(f"On-target GC factor computation failed: {e}")
                
                # Write motif outputs (EDM, BPM, MDS with z-score)
                write_motif_outputs(
                    result=extraction_result,
                    output_dir=output,
                    pon=pon,
                    include_ontarget=is_panel_mode,
                )
                
                fragment_count = extraction_result.fragment_count
                progress.update(task_extract, description=f"[1/5] Extract + Motif  ✓ ({fragment_count:,} frags)", completed=100)
                    
            except Exception as e:
                progress.update(task_extract, description=f"[1/5] Extract + Motif  ✗ Error", completed=100)
                logger.error(f"Unified Extract+Motif failed: {e}")
                if debug:
                    import traceback
                    traceback.print_exc()
                raise typer.Exit(1)

        # ═══════════════════════════════════════════════════════════════════
        # 2. UNIFIED FEATURE PIPELINE (FSC, FSR, FSD, WPS, OCF)
        # Using consolidated unified_processor for all feature extraction
        # ═══════════════════════════════════════════════════════════════════
        progress.start_task(task_pipeline)
        progress.update(task_pipeline, description="[2/5] Unified Pipeline ...")
        
        # Resolve optional paths for unified processor
        resolved_wps_anchors = _resolve_path(wps_anchors)
        resolved_wps_file = _resolve_path(wps_file)
        resolved_wps_background = _resolve_path(wps_background)
        resolved_assay = assay if isinstance(assay, str) else None
        
        try:
            from .core.unified_processor import run_features
            
            # Call unified processor for FSC, FSR, FSD, WPS, OCF
            # This handles all Rust pipeline calls and Python post-processing
            feature_outputs = run_features(
                bed_path=bedgz_file,
                output_dir=output,
                sample_name=sample,
                genome=genome,
                enable_fsc=True,
                enable_fsr=True,
                enable_fsd=True,
                enable_wps=True,
                enable_ocf=assets.ocf_available,  # Skip OCF if not available for genome
                enable_tfbs=assets.tfbs_available,  # TFBS size entropy
                enable_atac=assets.atac_available,  # ATAC size entropy
                target_regions=resolved_target_regions,
                assay=resolved_assay,
                pon_model=actual_pon_path,
                skip_pon_zscore=skip_pon,  # --skip-pon: skip PON z-score normalization
                fsc_bins=bin_input,
                fsc_windows=fsc_windows,
                fsc_continue_n=fsc_continue_n,
                fsd_arms=resolved_arms_file,
                wps_anchors=resolved_wps_anchors or resolved_wps_file,
                wps_background=resolved_wps_background,
                wps_bait_padding=resolved_bait_padding,
                ocf_regions=resolved_ocr_file,
                gc_correct=True,
                threads=threads,
                verbose=debug,
            )
            
            # Report WPS periodicity if available
            if feature_outputs.wps and feature_outputs.wps.exists():
                from .core.wps_processor import post_process_wps
                wps_result = post_process_wps(
                    wps_parquet=feature_outputs.wps,
                    wps_background_parquet=feature_outputs.wps_background if feature_outputs.wps_background and feature_outputs.wps_background.exists() else None,
                    smooth=True,
                    extract_periodicity=True
                )
                if wps_result.get("periodicity_score"):
                    logger.info(f"WPS periodicity score: {wps_result['periodicity_score']:.3f}")
            
            progress.update(task_pipeline, description="[2/5] Unified Pipeline ✓", completed=100)
            
            # ═══════════════════════════════════════════════════════════════════
            # 3. POST-PROCESSING (already done by unified processor)
            # ═══════════════════════════════════════════════════════════════════
            progress.start_task(task_postproc)
            progress.update(task_postproc, description="[3/5] Post-processing ✓ (included)", completed=100)
            
        except Exception as e:
            progress.update(task_pipeline, description="[2/5] Unified Pipeline ✗ Error", completed=100)
            logger.error(f"Unified Pipeline failed: {e}")
            if debug:
                import traceback
                traceback.print_exc()
            raise typer.Exit(1)

        # ═══════════════════════════════════════════════════════════════════
        # 4. OPTIONAL TOOLS (UXM and mFSD)
        # ═══════════════════════════════════════════════════════════════════
        
        # 4a. UXM
        if has_uxm:
            progress.start_task(task_uxm)
            progress.update(task_uxm, description="[4/5] UXM ...")
            try:
                from .uxm import uxm
                uxm(resolved_bisulfite_bam, output, sample, None, 0)
                progress.update(task_uxm, description="[4/5] UXM ✓", completed=100)
            except Exception as e:
                progress.update(task_uxm, description="[4/5] UXM ✗ Error", completed=100)
                logger.warning(f"UXM failed: {e}")
        
        # 4b. mFSD
        if has_mfsd:
            progress.start_task(task_mfsd)
            progress.update(task_mfsd, description="[5/5] mFSD ...")
            try:
                from .mfsd import mfsd
                # Resolve GC correction factors using the same logic as other tools
                from .core.gc_assets import resolve_gc_assets
                gc = resolve_gc_assets(assets, output, sample, bedgz_file, gc_correct=True, genome=genome)
                
                # Panel mode: prefer on-target GC factors for mFSD
                # Variants are in captured regions (exons/promoters), which have different
                # GC bias than off-target genome-wide fragments
                if is_panel_mode:
                    # On-target factors are generated by extract.py as .ontarget.tsv
                    ontarget_factors = output / f"{sample}.correction_factors.ontarget.tsv"
                    if ontarget_factors.exists():
                        correction_factors_path = ontarget_factors
                        logger.info(f"mFSD: Using on-target GC factors ({ontarget_factors.name})")
                    else:
                        correction_factors_path = gc.factors_input or gc.factors_out
                        logger.warning(f"mFSD: On-target factors not found, using global factors")
                else:
                    correction_factors_path = gc.factors_input or gc.factors_out
                
                mfsd(
                    bam_input=bam_input,
                    input_file=resolved_variants,
                    output=output,
                    sample_name=sample,
                    reference=reference,
                    correction_factors=correction_factors_path if correction_factors_path and correction_factors_path.exists() else None,
                    mapq=mapq,
                    minlen=minlen,
                    maxlen=maxlen,
                    require_proper_pair=require_proper_pair,
                    duplex=duplex,
                    output_distributions=False,
                    verbose=debug,
                    threads=threads
                )
                progress.update(task_mfsd, description="[5/5] mFSD ✓", completed=100)
            except Exception as e:
                progress.update(task_mfsd, description="[5/5] mFSD ✗ Error", completed=100)
                logger.warning(f"mFSD failed: {e}")
    
    # ═══════════════════════════════════════════════════════════════════
    # UNIFIED JSON OUTPUT (if requested)
    # ═══════════════════════════════════════════════════════════════════
    if generate_json:
        try:
            from .core.feature_serializer import FeatureSerializer
            
            logger.info("Generating unified features JSON...")
            serializer = FeatureSerializer.from_outputs(sample, output, version="0.3.2")
            
            # Add runtime metadata
            serializer.add_metadata("bam_path", str(bam_input))
            serializer.add_metadata("genome", genome)
            serializer.add_metadata("panel_mode", is_panel_mode)
            if pon:
                serializer.add_metadata("pon_assay", pon.assay)
                serializer.add_metadata("pon_n_samples", pon.n_samples)
            
            json_path = output / sample
            serializer.save(json_path)
            logger.info(f"Unified JSON saved: {json_path}.features.json")
        except Exception as e:
            logger.warning(f"JSON generation failed: {e}")
    
    logger.info(f"✅ All feature extraction complete: {output}")

