"""
Run all krewlyzer feature extraction tools for a single sample.
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

from .core.fsc_processor import process_fsc
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

def run_all(
    bam_input: Path = typer.Argument(..., help="Input BAM file (sorted, indexed)"),
    reference: Path = typer.Option(..., "--reference", "-g", help="Reference genome FASTA (indexed)"),
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
    
    # Observability
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging"),
):
    """
    Run all feature extraction tools for a single sample.
    
    Pipeline: extract → motif → [Unified Engine: FSC/FSR, FSD, WPS, OCF]
    Optional: uxm, mfsd
    """
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
    is_panel_mode = target_regions and isinstance(target_regions, Path) and target_regions.exists()
    
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
    
    if is_panel_mode:
        logger.info(f"Panel mode: GC model will use off-target reads only (targets: {target_regions.name})")
    
    # Determine which optional tools will run
    has_uxm = bisulfite_bam is not None and bisulfite_bam.exists()
    has_mfsd = variants is not None and variants.exists()
    
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
                # Returns: (count, em_off, bpm_off, gc_obs, em_on, bpm_on)
                fragment_count, em_counts, bpm_counts, gc_observations, em_counts_on, bpm_counts_on = _core.extract_motif.process_bam_parallel(
                    str(bam_input),
                    str(reference),
                    mapq,
                    minlen,
                    maxlen,
                    4,  # kmer
                    threads,
                    bed_out_arg,
                    "enable",  # output_motif_prefix
                    str(exclude_regions) if exclude_regions else str(assets.exclude_regions),
                    str(target_regions) if (target_regions and isinstance(target_regions, Path) and target_regions.exists()) else None,  # Panel mode: off-target GC
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
                    factors_file = output / f"{sample}.correction_factors.csv"
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
                    metadata = {
                        "sample_id": sample,
                        "total_fragments": fragment_count,
                        "genome": genome,
                        "gc_correction_computed": factors_file.exists(),
                        "filters": { "mapq": mapq, "min_length": minlen, "max_length": maxlen },
                        "timestamp": datetime.now().isoformat()
                    }
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
                is_panel_mode = target_regions and isinstance(target_regions, Path) and target_regions.exists()
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
        res_arms = arms_file if arms_file else assets.arms
        if not res_arms.exists(): logger.error(f"Arms file missing: {res_arms}"); raise typer.Exit(1)

        # WPS Regions (prefer wps_anchors, fallback to wps_file, then asset defaults)
        if wps_anchors and wps_anchors.exists():
            res_wps = wps_anchors
            logger.debug(f"WPS anchors: {wps_anchors}")
        elif wps_file and wps_file.exists():
            res_wps = wps_file
            logger.debug(f"WPS transcript (legacy): {wps_file}")
        else:
            # Try wps_anchors asset first, fallback to transcript_anno
            try:
                res_wps = assets.wps_anchors
                logger.debug(f"WPS anchors (default): {res_wps}")
            except (AttributeError, FileNotFoundError):
                res_wps = assets.transcript_anno
                logger.debug(f"WPS transcript (fallback): {res_wps}")
        if not res_wps.exists(): logger.error(f"WPS regions file missing: {res_wps}"); raise typer.Exit(1)

        # OCF Regions
        res_ocf = ocr_file if ocr_file else assets.ocf_regions
        if not res_ocf.exists(): logger.error(f"OCR file missing: {res_ocf}"); raise typer.Exit(1)
        
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
        out_gc_factors = output / f"{sample}.correction_factors.csv"
        out_fsc_raw = output / f"{sample}.fsc_counts.tsv"
        out_wps = output / f"{sample}.WPS.parquet"  # Foreground Parquet
        out_wps_bg = output / f"{sample}.WPS_background.parquet"  # Alu stacking
        out_fsd = output / f"{sample}.FSD.tsv"
        out_ocf_dir = output / f"{sample}_ocf_tmp"
        out_ocf_dir.mkdir(parents=True, exist_ok=True)
        
        # Resolve WPS background (explicit --wps-background takes precedence over auto-load)
        if wps_background and wps_background.exists():
            res_wps_bg = wps_background
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
                str(res_ocf), str(out_ocf_dir),
                str(target_regions) if (target_regions and isinstance(target_regions, Path) and target_regions.exists()) else None,
                bait_padding,  # Bait edge padding (adaptive safety applies)
                True  # silent=True
            )
            progress.update(task_pipeline, description="[2/5] Unified Pipeline ✓ (FSC+WPS+FSD+OCF)", completed=100)
            
            # ═══════════════════════════════════════════════════════════════════
            # 3. POST-PROCESSING (FSC/FSR, OCF, PON z-scores)
            # ═══════════════════════════════════════════════════════════════════
            progress.start_task(task_postproc)
            progress.update(task_postproc, description="[3/5] Post-processing ...")
            
            # FSC & FSR (From same raw counts)
            if out_fsc_raw.exists():
                df_counts = pd.read_csv(out_fsc_raw, sep='\t')
                
                # Load PON if available
                pon = None
                if pon_model:
                    from .core.pon_integration import load_pon_model
                    pon = load_pon_model(pon_model)
                
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
        
            # OCF (Move files)
            src_ocf = out_ocf_dir / "all.ocf.tsv"
            src_sync = out_ocf_dir / "all.sync.tsv"
            dst_ocf = output / f"{sample}.OCF.tsv"
            dst_sync = output / f"{sample}.OCF.sync.tsv"
            
            if src_ocf.exists(): shutil.move(str(src_ocf), str(dst_ocf))
            if src_sync.exists(): shutil.move(str(src_sync), str(dst_sync))
            
            try:
                out_ocf_dir.rmdir()
            except:
                pass
            
            # Apply processing to FSD and WPS (PoN if provided)
            if pon_model:
                from .core.pon_integration import load_pon_model
                pon = load_pon_model(pon_model)
            else:
                pon = None
            
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
            pon_wps_baseline = pon.wps_baseline if pon else None
            wps_result = post_process_wps(
                wps_parquet=out_wps,
                wps_background_parquet=out_wps_bg if out_wps_bg.exists() else None,
                pon_baseline_parquet=None,  # TODO: extract WPS baseline path from PoN
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
                uxm(bisulfite_bam, output, sample, None, 0)
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
                    input_file=variants,
                    output=output,
                    sample_name=sample,
                    reference=reference,
                    correction_factors=out_gc_factors if out_gc_factors.exists() else None,
                    mapq=mapq,
                    output_distributions=False,
                    verbose=debug,
                    threads=threads
                )
                progress.update(task_mfsd, description="[5/5] mFSD ✓", completed=100)
            except Exception as e:
                progress.update(task_mfsd, description="[5/5] mFSD ✗ Error", completed=100)
                logger.warning(f"mFSD failed: {e}")
    
    logger.info(f"✅ All feature extraction complete: {output}")
