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
    wps_file: Optional[Path] = typer.Option(None, "--wps-file", help="Custom transcript file for WPS"),
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
    
    # Window settings for FSC/FSR
    if bin_input:
        fsc_windows, fsc_continue_n = 1, 1
        logger.info(f"Using custom bin file: {bin_input}")
    else:
        fsc_windows, fsc_continue_n = 100000, 50
    
    logger.info(f"Processing sample: {sample}")
    logger.info(f"Filters: mapq>={mapq}, length=[{minlen},{maxlen}], skip_dup={skip_duplicates}, proper_pair={require_proper_pair}")
    
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
    # Always run motif logic as part of pass if not present?
    # Or just re-run everything if any output missing?
    # Simple logic: If BED or ANY Motif file missing, run the Engine.
    should_run_engine = (
        not bedgz_file.exists() or 
        not edm_output.exists() or 
        not bpm_output.exists() or 
        not mds_output.exists()
    )
    
    if not should_run_engine:
        logger.info("Extract and Motif outputs exist. Skipping step 1 & 2.")
    else:
        logger.info("Running Unified Extract + Motif...")
        
        # Decide if we write BED
        bed_out_arg = str(bed_temp) if should_run_extract else None
        
        try:
            # Call Unified Engine
            # Returns (total_count, end_motifs, bp_motifs, gc_observations)
            fragment_count, em_counts, bpm_counts, gc_observations = _core.extract_motif.process_bam_parallel(
                str(bam_input),
                str(reference),
                mapq,
                minlen,
                maxlen,
                4, # kmer hardcoded to 4 as per tool standard
                threads,
                bed_out_arg,          # Write BED if needed
                "enable",             # Always calculate motifs
                str(exclude_regions) if exclude_regions else str(assets.exclude_regions),
                skip_duplicates,
                require_proper_pair
            )
            
            logger.info(f"Processed {fragment_count:,} fragments")
            
            # --- Post-Process Extract (Move BGZF BED, compute GC factors) ---
            if should_run_extract and bed_temp.exists():
                import shutil
                import pysam
                
                # Rust already writes BGZF - just move and index
                logger.info("Moving and indexing BED.gz...")
                shutil.move(str(bed_temp), str(bedgz_file))
                pysam.tabix_index(str(bedgz_file), preset="bed", force=True)
                
                # Compute GC correction factors inline
                factors_file = output / f"{sample}.correction_factors.csv"
                try:
                    gc_ref = assets.resolve("gc_reference")
                    valid_regions = assets.resolve("valid_regions")
                    logger.info(f"Computing GC factors from {len(gc_observations)} observation bins...")
                    n_factors = _core.gc.compute_and_write_gc_factors(
                        gc_observations, str(gc_ref), str(valid_regions), str(factors_file)
                    )
                    logger.info(f"Computed {n_factors} GC correction factors")
                except Exception as e:
                    logger.warning(f"GC factor computation failed: {e}")
                
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
            
            logger.info(f"Motif extraction complete (EM={total_em:,}, BPM={total_bpm:,}, MDS={mds:.4f})")
                
        except Exception as e:
            logger.error(f"Unified Extract+Motif failed: {e}")
            raise typer.Exit(1)


    # ═══════════════════════════════════════════════════════════════════
    # 3. UNIFIED SINGLE-PASS PIPELINE (FSC, FSR, FSD, WPS, OCF)
    # ═══════════════════════════════════════════════════════════════════
    logger.info("🚀 Running Unified Single-Pass Pipeline (FSC, FSR, FSD, WPS, OCF)...")

    # Define Resource Paths (Resolve defaults via AssetManager)
    
    # FSC/FSR Bins
    res_bin = bin_input if bin_input else assets.bins_100kb
    if not res_bin.exists(): logger.error(f"Bin file missing: {res_bin}"); raise typer.Exit(1)

    # FSD Arms
    res_arms = arms_file if arms_file else assets.arms
    if not res_arms.exists(): logger.error(f"Arms file missing: {res_arms}"); raise typer.Exit(1)

    # WPS Genes (Transcript Annotation)
    res_wps = wps_file if wps_file else assets.transcript_anno
    if not res_wps.exists(): logger.error(f"Transcript file missing: {res_wps}"); raise typer.Exit(1)

    # OCF Regions
    res_ocf = ocr_file if ocr_file else assets.ocf_regions
    if not res_ocf.exists(): logger.error(f"OCR file missing: {res_ocf}"); raise typer.Exit(1)
    
    # GC Reference (Mandatory for GCfix correction)
    res_gc = assets.gc_reference
    # We check existence in Rust or here? Here is better for CLI feedback.
    if not res_gc.exists():
        logger.error(f"GC reference missing: {res_gc}")
        logger.error("Run 'krewlyzer build-gc-reference' or ensure bundled assets are present.")
        raise typer.Exit(1)
        
    res_valid_regions = assets.valid_regions
    if not res_valid_regions.exists():
        logger.error(f"Valid regions file missing: {res_valid_regions}")
        raise typer.Exit(1)

    # Define Outputs
    out_gc_factors = output / f"{sample}.correction_factors.csv"
    out_fsc_raw = output / f"{sample}.fsc_counts.tsv"
    out_wps = output / f"{sample}.WPS.tsv.gz"
    out_fsd = output / f"{sample}.FSD.tsv"
    out_ocf_dir = output / f"{sample}_ocf_tmp" # OCF writes mult files to dir
    out_ocf_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # RUN RUST PIPELINE
        _core.run_unified_pipeline(
            str(bedgz_file),
            # GC Correction (compute fresh for run-all)
            str(res_gc), str(res_valid_regions), str(out_gc_factors),
            # Pre-computed factors (None = compute above)
            None,
            # FSC
            str(res_bin), str(out_fsc_raw),
            # WPS
            str(res_wps), str(out_wps), False,
            # FSD
            str(res_arms), str(out_fsd),
            # OCF
            str(res_ocf), str(out_ocf_dir)
        )
        logger.info("✅ Unified Pipeline execution complete.")
        
        # Post-Processing
        
        # 1. FSC & FSR (From same raw counts)
        if out_fsc_raw.exists():
            logger.info("Post-processing FSC/FSR...")
            df_counts = pd.read_csv(out_fsc_raw, sep='\t')
            
            # Load PON if available
            pon = None
            if pon_model:
                from .core.pon_integration import load_pon_model
                pon = load_pon_model(pon_model)
            
            # FSC Output (using shared processor)
            final_fsc = output / f"{sample}.FSC.tsv"
            process_fsc(df_counts, final_fsc, fsc_windows, fsc_continue_n, pon=pon)
            
            # FSR Output (using shared processor)
            final_fsr = output / f"{sample}.FSR.tsv"
            process_fsr(df_counts, final_fsr, fsc_windows, fsc_continue_n, pon=pon)
    
        # 2. OCF (Move files)
        # Rust writes 'all.ocf.csv' and 'all.sync.tsv' to out_ocf_dir
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
        
        # 3. Apply PON z-scores to FSD and WPS (if PON model provided)
        if pon_model:
            from .core.pon_integration import load_pon_model
            from .core.fsd_processor import apply_fsd_pon
            from .core.wps_processor import apply_wps_pon
            
            pon = load_pon_model(pon_model)
            
            if pon is not None:
                # Apply PON z-scores to FSD
                if out_fsd.exists():
                    apply_fsd_pon(out_fsd, pon)
                
                # Apply PON z-scores to WPS
                if out_wps.exists():
                    apply_wps_pon(out_wps, pon)
                    
                logger.info("PON z-scores applied to FSD and WPS")
            
    except Exception as e:
            logger.error(f"Unified Pipeline failed: {e}")
            raise typer.Exit(1)


    # ═══════════════════════════════════════════════════════════════════
    # 4. OPTIONAL TOOLS (Always run if requested)
    # ═══════════════════════════════════════════════════════════════════
    # 4a. UXM
    if bisulfite_bam:
        if not bisulfite_bam.exists():
            logger.warning(f"Bisulfite BAM not found: {bisulfite_bam}. Skipping UXM.")
        else:
            logger.info("Running UXM...")
            try:
                # UXM submodule functions are imported in __init__.py usually? 
                # Check import at top. import uxm from .uxm
                from .uxm import uxm
                uxm(bisulfite_bam, output, sample, None, 0)
            except Exception as e:
                logger.warning(f"UXM failed: {e}")
    else:
        logger.info("Skipping UXM (no --bisulfite-bam provided)")
    
    # 4b. mFSD
    if variants:
        if not variants.exists():
            logger.warning(f"Variants file not found: {variants}. Skipping mFSD.")
        else:
            logger.info("Running mFSD...")
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
            except Exception as e:
                logger.warning(f"mFSD failed: {e}")
    else:
        logger.info("Skipping mFSD (no --variants provided)")
    
    logger.info(f"✅ All feature extraction complete: {output}")
