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
import shutil
import pandas as pd

from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from .core.fsc_processor import process_fsc, aggregate_by_gene
from .core.fsr_processor import process_fsr
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

def _resolve_path(value) -> Optional[Path]:
    """Safely resolve a path value, handling typer.OptionInfo objects.
    
    When a typer-decorated function is called directly (not via CLI),
    Optional[Path] parameters with defaults remain as OptionInfo objects.
    This helper handles that case.
    """
    if value is None:
        return None
    if isinstance(value, Path):
        return value
    if isinstance(value, str):
        return Path(value)
    # If it's a typer.OptionInfo (or any other type), return None
    return None

def _resolve_int(value, default: int) -> int:
    """Safely resolve an integer value, handling typer.OptionInfo objects.
    
    When a typer-decorated function is called directly (not via CLI),
    parameters with defaults remain as OptionInfo objects.
    This helper handles that case.
    """
    if isinstance(value, int):
        return value
    # If it's a typer.OptionInfo or other non-int type, return default
    return default

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
    # Load once here and pass to all processors (FSC, FSR, FSD, WPS, OCF)
    pon = None
    if resolved_pon_model:
        from .core.pon_integration import load_pon_model
        pon = load_pon_model(resolved_pon_model)
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
        # 1. EXTRACT + MOTIF (Unified Engine)
        # ═══════════════════════════════════════════════════════════════════
    
        bedgz_file = output / f"{sample}.bed.gz"
        bed_temp = output / f"{sample}.bed.tmp"
        
        # Motif Outputs
        edm_output = output / f"{sample}.EndMotif.tsv"
        bpm_output = output / f"{sample}.BreakPointMotif.tsv"
        mds_output = output / f"{sample}.MDS.tsv"
        
        should_run_extract = not bedgz_file.exists()
        should_run_engine = (
            not bedgz_file.exists() or 
            not edm_output.exists() or 
            not bpm_output.exists() or 
            not mds_output.exists()
        )
        
        if not should_run_engine:
            progress.update(task_extract, description="[1/5] Extract + Motif  ✓ (cached)", completed=100)
        else:
            progress.update(task_extract, description="[1/5] Extract + Motif  ...")
            
            # Decide if we write BED
            bed_out_arg = str(bed_temp) if should_run_extract else None
            
            try:
                # Call Unified Engine (silent=True to hide Rust progress)
                # Returns: (count, em_off, bpm_off, gc_obs, em_on, bpm_on, gc_obs_ontarget)
                fragment_count, em_counts, bpm_counts, gc_observations, em_counts_on, bpm_counts_on, gc_observations_ontarget = _core.extract_motif.process_bam_parallel(
                    str(bam_input),
                    str(reference),
                    mapq,
                    minlen,
                    maxlen,
                    4,  # kmer
                    threads,
                    bed_out_arg,
                    "enable",  # output_motif_prefix
                    str(resolved_exclude_regions) if resolved_exclude_regions else str(assets.exclude_regions),
                    str(resolved_target_regions) if resolved_target_regions and resolved_target_regions.exists() else None,  # Panel mode: off-target GC
                    skip_duplicates,
                    require_proper_pair,
                    True  # silent=True
                )
                
                # --- Post-Process Extract (Move BGZF BED, compute GC factors) ---
                if should_run_extract and bed_temp.exists():
                    import pysam
                    
                    # Rust already writes BGZF - just move and index
                    shutil.move(str(bed_temp), str(bedgz_file))
                    pysam.tabix_index(str(bedgz_file), preset="bed", force=True)
                    
                    # Compute GC correction factors inline
                    factors_file = output / f"{sample}.correction_factors.tsv"
                    try:
                        gc_ref = assets.resolve("gc_reference")
                        valid_regions = assets.resolve("valid_regions")
                        n_factors = _core.gc.compute_and_write_gc_factors(
                            gc_observations, str(gc_ref), str(valid_regions), str(factors_file)
                        )
                    except Exception as e:
                        logger.debug(f"GC factor computation failed: {e}")
                    
                    # Write Metadata
                    import json
                    from datetime import datetime
                    meta_file = output / f"{sample}.metadata.json"
                    
                    # Resolve assay
                    resolved_assay = assay if isinstance(assay, str) else None
                    
                    metadata = {
                        "sample_id": sample,
                        "total_fragments": fragment_count,
                        "genome": genome,
                        "assay": resolved_assay,
                        "gc_correction_computed": factors_file.exists(),
                        "panel_mode": is_panel_mode,
                        "filters": { "mapq": mapq, "min_length": minlen, "max_length": maxlen },
                        "timestamp": datetime.now().isoformat()
                    }
                    
                    # Add on-target rate if available (from Rust output)
                    if is_panel_mode and hasattr(result, 'on_target_rate'):
                        metadata["on_target_rate"] = result.on_target_rate
                    
                    with open(meta_file, 'w') as f:
                        json.dump(metadata, f, indent=2)

                # --- Post-Process Motif (Write Files) ---
                from .core.motif_processor import process_motif_outputs
                
                # Off-target motifs (primary - unbiased for biomarkers)
                total_em, total_bpm, mds = process_motif_outputs(
                    em_counts=em_counts,
                    bpm_counts=bpm_counts,
                    edm_output=edm_output,
                    bpm_output=bpm_output,
                    mds_output=mds_output,
                    sample_name=sample,
                    kmer=4,
                    include_headers=True
                )
                
                # On-target motifs (for panel data comparison - PCR biased)
                if is_panel_mode and sum(em_counts_on.values()) > 0:
                    edm_on = output / f"{sample}.EndMotif.ontarget.tsv"
                    bpm_on = output / f"{sample}.BreakPointMotif.ontarget.tsv"
                    mds_on = output / f"{sample}.MDS.ontarget.tsv"
                    
                    total_em_on, total_bpm_on, mds_on_val = process_motif_outputs(
                        em_counts=em_counts_on,
                        bpm_counts=bpm_counts_on,
                        edm_output=edm_on,
                        bpm_output=bpm_on,
                        mds_output=mds_on,
                        sample_name=sample,
                        kmer=4,
                        include_headers=True
                    )
                    logger.info(f"Motif on-target: {total_em_on:,} EM, {total_bpm_on:,} BPM")
                
                # MDS z-score using PON baseline
                if pon and pon.mds_baseline and mds is not None:
                    mds_z = (mds - pon.mds_baseline.mds_mean) / max(pon.mds_baseline.mds_std, 1e-10)
                    logger.debug(f"MDS z-score: {mds_z:.3f} (raw={mds:.4f}, pon_mean={pon.mds_baseline.mds_mean:.4f})")
                    
                    # Write z-score to MDS file if it exists
                    if mds_output.exists():
                        try:
                            mds_df = pd.read_csv(mds_output, sep="\t")
                            if "mds_z" not in mds_df.columns:
                                mds_df["mds_z"] = mds_z
                                mds_df.to_csv(mds_output, sep="\t", index=False)
                        except Exception as mds_e:
                            logger.debug(f"Could not add MDS z-score: {mds_e}")
                
                progress.update(task_extract, description=f"[1/5] Extract + Motif  ✓ ({fragment_count:,} frags)", completed=100)
                    
            except Exception as e:
                progress.update(task_extract, description=f"[1/5] Extract + Motif  ✗ Error", completed=100)
                logger.error(f"Unified Extract+Motif failed: {e}")
                raise typer.Exit(1)

        # ═══════════════════════════════════════════════════════════════════
        # 2. UNIFIED SINGLE-PASS PIPELINE (FSC, FSR, FSD, WPS, OCF)
        # ═══════════════════════════════════════════════════════════════════
        progress.start_task(task_pipeline)
        progress.update(task_pipeline, description="[2/5] Unified Pipeline ...")

        # Define Resource Paths (Resolve defaults via AssetManager)
        
        # FSC/FSR Bins
        res_bin = bin_input if bin_input else assets.bins_100kb
        if not res_bin.exists(): logger.error(f"Bin file missing: {res_bin}"); raise typer.Exit(1)

        # FSD Arms
        res_arms = resolved_arms_file if resolved_arms_file else assets.arms
        if not res_arms.exists(): logger.error(f"Arms file missing: {res_arms}"); raise typer.Exit(1)

        # WPS Regions (prefer wps_anchors, fallback to wps_file, then asset defaults)
        # Resolve paths to handle typer.OptionInfo when called directly
        resolved_wps_anchors = _resolve_path(wps_anchors)
        resolved_wps_file = _resolve_path(wps_file)
        resolved_wps_background = _resolve_path(wps_background)
        
        if resolved_wps_anchors and resolved_wps_anchors.exists():
            res_wps = resolved_wps_anchors
            logger.debug(f"WPS anchors: {resolved_wps_anchors}")
        elif resolved_wps_file and resolved_wps_file.exists():
            res_wps = resolved_wps_file
            logger.debug(f"WPS transcript (legacy): {resolved_wps_file}")
        else:
            # Resolve assay for panel-specific assets
            resolved_assay = assay if isinstance(assay, str) else None
            
            # For dual WPS: always use genome-wide anchors for primary WPS
            # Panel-specific anchors are used in a second pass
            try:
                res_wps = assets.wps_anchors
                logger.debug(f"WPS anchors (genome-wide): {res_wps}")
            except (AttributeError, FileNotFoundError):
                res_wps = assets.transcript_anno
                logger.debug(f"WPS transcript (fallback): {res_wps}")
            
            # Track panel anchors for dual WPS output (will be used after primary run)
            res_wps_panel = None
            if resolved_assay:
                try:
                    res_wps_panel = assets.get_wps_anchors(resolved_assay)
                    logger.info(f"Dual WPS enabled: panel anchors ({resolved_assay}) will generate WPS.panel.parquet")
                except Exception:
                    pass
        if not res_wps.exists(): logger.error(f"WPS regions file missing: {res_wps}"); raise typer.Exit(1)

        # OCF Regions (only available for GRCh37/hg19 unless user provides custom file)
        run_ocf = True
        if resolved_ocr_file:
            res_ocf = resolved_ocr_file
        elif assets.ocf_available:
            res_ocf = assets.ocf_regions
        else:
            logger.warning(f"⚠️  OCF regions not available for {genome} (only GRCh37/hg19 currently supported)")
            logger.warning("   Skipping OCF. Provide --ocr-file to run OCF on hg38")
            run_ocf = False
            res_ocf = None
        
        if run_ocf and res_ocf and not res_ocf.exists():
            logger.error(f"OCR file missing: {res_ocf}"); raise typer.Exit(1)
        
        # GC Reference
        res_gc = assets.gc_reference
        if not res_gc.exists():
            logger.error(f"GC reference missing: {res_gc}")
            raise typer.Exit(1)
            
        res_valid_regions = assets.valid_regions
        if not res_valid_regions.exists():
            logger.error(f"Valid regions file missing: {res_valid_regions}")
            raise typer.Exit(1)

        # Define Outputs
        out_gc_factors = output / f"{sample}.correction_factors.tsv"
        out_fsc_raw = output / f"{sample}.fsc_counts.tsv"
        out_wps = output / f"{sample}.WPS.parquet"  # Foreground Parquet
        out_wps_bg = output / f"{sample}.WPS_background.parquet"  # Alu stacking
        out_fsd = output / f"{sample}.FSD.tsv"
        out_ocf_dir = output / f"{sample}_ocf_tmp"
        out_ocf_dir.mkdir(parents=True, exist_ok=True)
        
        # Resolve WPS background (explicit --wps-background takes precedence over auto-load)
        if resolved_wps_background and resolved_wps_background.exists():
            res_wps_bg = resolved_wps_background
            logger.info(f"Using explicit WPS background: {res_wps_bg}")
        elif assets.wps_background.exists():
            res_wps_bg = assets.wps_background
            logger.debug(f"Auto-loaded WPS background ({genome}): {res_wps_bg}")
        else:
            res_wps_bg = None
            logger.warning("WPS background not found - skipping background stacking")
        
        try:
            # RUN RUST PIPELINE (silent=True)
            _core.run_unified_pipeline(
                str(bedgz_file),
                str(res_gc), str(res_valid_regions), str(out_gc_factors),
                None,
                str(res_bin), str(out_fsc_raw),
                str(res_wps), str(out_wps),  # WPS foreground
                str(res_wps_bg) if res_wps_bg else None, str(out_wps_bg) if res_wps_bg else None,  # WPS background
                False,
                str(res_arms), str(out_fsd),
                str(res_ocf) if run_ocf else None, str(out_ocf_dir) if run_ocf else None,  # OCF (skip for hg38)
                str(resolved_target_regions) if resolved_target_regions and resolved_target_regions.exists() else None,
                resolved_bait_padding,  # Bait edge padding (adaptive safety applies)
                True  # silent=True
            )
            ocf_status = "OCF" if run_ocf else "OCF skipped"
            progress.update(task_pipeline, description=f"[2/5] Unified Pipeline ✓ (FSC+WPS+FSD+{ocf_status})", completed=100)
            
            # === DUAL WPS: Run panel-specific WPS if --assay was provided ===
            if 'res_wps_panel' in dir() and res_wps_panel and res_wps_panel.exists():
                out_wps_panel = output / f"{sample}.WPS.panel.parquet"
                logger.info(f"Running panel-specific WPS ({res_wps_panel.name})...")
                _core.run_unified_pipeline(
                    str(bedgz_file),
                    None, None, None,  # GC already computed
                    str(out_gc_factors) if out_gc_factors.exists() else None,  # Load pre-computed GC
                    None, None,  # No FSC needed
                    str(res_wps_panel), str(out_wps_panel),  # Panel WPS
                    None, None,  # No background for panel
                    False,
                    None, None,  # No FSD
                    None, None,  # No OCF
                    str(resolved_target_regions) if resolved_target_regions and resolved_target_regions.exists() else None,
                    resolved_bait_padding,
                    True  # silent=True
                )
                logger.info(f"Panel WPS complete: {out_wps_panel.name}")
            
            # ═══════════════════════════════════════════════════════════════════
            # 3. POST-PROCESSING (FSC/FSR, OCF, PON z-scores)
            # ═══════════════════════════════════════════════════════════════════
            progress.start_task(task_postproc)
            progress.update(task_postproc, description="[3/5] Post-processing ...")
            
            # FSC & FSR (From same raw counts)
            if out_fsc_raw.exists():
                df_counts = pd.read_csv(out_fsc_raw, sep='\t')
                
                # Use centralized PON loaded earlier (line 176+)
                
                # FSC Output (off-target - primary for biomarkers)
                final_fsc = output / f"{sample}.FSC.tsv"
                process_fsc(df_counts, final_fsc, fsc_windows, fsc_continue_n, pon=pon)
                
                # FSR Output (off-target - primary for biomarkers)
                final_fsr = output / f"{sample}.FSR.tsv"
                process_fsr(df_counts, final_fsr, fsc_windows, fsc_continue_n, pon=pon)
                
                # Process on-target files if they exist (panel mode)
                out_fsc_raw_ontarget = output / f"{sample}_raw.FSC.ontarget.tsv"
                if out_fsc_raw_ontarget.exists():
                    logger.info("Processing on-target FSC/FSR (for CNV only)")
                    df_counts_on = pd.read_csv(out_fsc_raw_ontarget, sep='\t')
                    
                    # FSC on-target (for CNV analysis only)
                    final_fsc_on = output / f"{sample}.FSC.ontarget.tsv"
                    process_fsc(df_counts_on, final_fsc_on, fsc_windows, fsc_continue_n, pon=pon)
                    
                    # FSR on-target (for comparison only - PCR biased)
                    final_fsr_on = output / f"{sample}.FSR.ontarget.tsv"
                    process_fsr(df_counts_on, final_fsr_on, fsc_windows, fsc_continue_n, pon=pon)
                
                # Gene-centric FSC (if --assay provided)
                resolved_assay = assay if isinstance(assay, str) else None
                if resolved_assay:
                    from .core.gene_bed import load_gene_bed
                    try:
                        gene_regions = load_gene_bed(assay=resolved_assay, genome=genome)
                        gene_fsc_output = output / f"{sample}.FSC.gene.tsv"
                        aggregate_by_gene(
                            df_counts=df_counts,
                            gene_regions=gene_regions,
                            output_path=gene_fsc_output,
                            pon=pon
                        )
                        logger.info(f"Gene-centric FSC: {len(gene_regions)} genes -> {gene_fsc_output.name}")
                    except Exception as e:
                        logger.warning(f"Gene-centric FSC skipped: {e}")
        
            # OCF (Move files)
            src_ocf = out_ocf_dir / "all.ocf.tsv"
            src_sync = out_ocf_dir / "all.sync.tsv"
            dst_ocf = output / f"{sample}.OCF.tsv"
            dst_sync = output / f"{sample}.OCF.sync.tsv"
            
            if src_ocf.exists(): shutil.move(str(src_ocf), str(dst_ocf))
            if src_sync.exists(): shutil.move(str(src_sync), str(dst_sync))
            
            # Apply OCF PON z-scores if available
            if dst_ocf.exists() and pon and pon.ocf_baseline:
                from .core.ocf_processor import process_ocf_with_pon
                process_ocf_with_pon(dst_ocf, pon.ocf_baseline)
            
            try:
                out_ocf_dir.rmdir()
            except:
                pass
            
            # Use centralized PON loaded earlier (line 176+)
            
            # FSD post-processing (log-ratios if PoN)
            from .core.fsd_processor import process_fsd
            if out_fsd.exists():
                process_fsd(out_fsd, out_fsd, pon=pon)
                logger.debug(f"FSD off-target: {out_fsd}")
            
            # Check for on-target FSD (panel mode)
            out_fsd_ontarget = output / f"{sample}.FSD.ontarget.tsv"
            if out_fsd_ontarget.exists():
                process_fsd(out_fsd_ontarget, out_fsd_ontarget, pon=pon)
                logger.info(f"FSD on-target: {out_fsd_ontarget}")
            
            # WPS post-processing (smoothing, PoN subtraction, FFT periodicity)
            from .core.wps_processor import post_process_wps
            out_wps_bg = output / f"{sample}.WPS_background.parquet"
            
            # Save WPS baseline from PON model if available
            pon_wps_baseline_path = None
            if pon and pon.wps_baseline and pon.wps_baseline.regions is not None:
                pon_wps_baseline_path = output / f".{sample}.pon_wps_baseline.parquet"
                pon.wps_baseline.regions.to_parquet(pon_wps_baseline_path, index=False)
                logger.debug(f"Saved PON WPS baseline: {pon_wps_baseline_path}")
            
            wps_result = post_process_wps(
                wps_parquet=out_wps,
                wps_background_parquet=out_wps_bg if out_wps_bg.exists() else None,
                pon_baseline_parquet=pon_wps_baseline_path,
                smooth=True,
                extract_periodicity=True
            )
            if wps_result.get("periodicity_score"):
                logger.info(f"WPS periodicity score: {wps_result['periodicity_score']:.3f}")
            
            progress.update(task_postproc, description="[3/5] Post-processing ✓", completed=100)
                
        except Exception as e:
            progress.update(task_pipeline, description="[2/5] Unified Pipeline ✗ Error", completed=100)
            logger.error(f"Unified Pipeline failed: {e}")
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
                mfsd(
                    bam_input=bam_input,
                    input_file=resolved_variants,
                    output=output,
                    sample_name=sample,
                    reference=reference,
                    correction_factors=out_gc_factors if out_gc_factors.exists() else None,
                    mapq=mapq,
                    minlen=minlen,
                    maxlen=maxlen,
                    require_proper_pair=require_proper_pair,
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

