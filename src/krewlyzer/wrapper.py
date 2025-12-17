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

from .postprocess import process_fsc_from_counts, process_fsr_from_counts

# Rust backend
from krewlyzer import _core

# Initialize logging globally, but level will be set in run_all
console = Console(stderr=True)
logging.basicConfig(
    level="INFO", 
    handlers=[RichHandler(console=console, show_time=False, show_path=False)], 
    format="%(message)s"
)
logger = logging.getLogger("krewlyzer")

def run_all(
    bam_input: Path = typer.Argument(..., help="Input BAM file (sorted, indexed)"),
    reference: Path = typer.Option(..., "--reference", "-g", help="Reference genome FASTA (indexed)"),
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
            # Returns (total_count, end_motifs, bp_motifs)
            fragment_count, em_counts, bpm_counts = _core.extract_motif.process_bam_parallel(
                str(bam_input),
                str(reference),
                mapq,
                minlen,
                maxlen,
                4, # kmer hardcoded to 4 as per tool standard
                threads,
                bed_out_arg,          # Write BED if needed
                "enable",             # Always calculate motifs
                str(exclude_regions) if exclude_regions else None,
                skip_duplicates,
                require_proper_pair
            )
            
            logger.info(f"Processed {fragment_count:,} fragments")
            
            # --- Post-Process Extract (Compress BED) ---
            if should_run_extract and bed_temp.exists():
                logger.info("Compressing BED...")
                import pysam
                import os
                pysam.tabix_compress(str(bed_temp), str(bedgz_file), force=True)
                pysam.tabix_index(str(bedgz_file), preset="bed", force=True)
                os.remove(str(bed_temp))
                
                # Write Metadata
                import json
                from datetime import datetime
                meta_file = output / f"{sample}.metadata.json"
                metadata = {
                    "sample_id": sample,
                    "total_fragments": fragment_count,
                    "filters": { "mapq": mapq, "min_length": minlen, "max_length": maxlen },
                    "timestamp": datetime.now().isoformat()
                }
                with open(meta_file, 'w') as f:
                    json.dump(metadata, f, indent=2)

            # --- Post-Process Motif (Write Files) ---
            # Helper imports (should be at top, but adding here locally or will add global)
            import itertools
            import numpy as np
            
            bases = ['A', 'C', 'T', 'G']
            kmer = 4
            
            # End Motif
            End_motif = {''.join(i): 0 for i in itertools.product(bases, repeat=kmer)}
            End_motif.update(em_counts)
            total_em = sum(End_motif.values())
            
            logger.info(f"Writing End Motif: {edm_output}")
            with open(edm_output, 'w') as f:
                f.write("Motif\tFrequency\n")
                for k, v in End_motif.items():
                    f.write(f"{k}\t{v/total_em if total_em else 0}\n")
            
            # Breakpoint Motif
            Breakpoint_motif = {''.join(i): 0 for i in itertools.product(bases, repeat=kmer)}
            Breakpoint_motif.update(bpm_counts)
            total_bpm = sum(Breakpoint_motif.values())
            
            logger.info(f"Writing Breakpoint Motif: {bpm_output}")
            with open(bpm_output, 'w') as f:
                f.write("Motif\tFrequency\n")
                for k, v in Breakpoint_motif.items():
                    f.write(f"{k}\t{v/total_bpm if total_bpm else 0}\n")
            
            # MDS
            logger.info(f"Writing MDS: {mds_output}")
            freq = np.array(list(End_motif.values())) / total_em if total_em else np.zeros(len(End_motif))
            mds = -np.sum(freq * np.log2(freq + 1e-12)) / np.log2(len(freq))
            with open(mds_output, 'w') as f:
                f.write("Sample\tMDS\n")
                f.write(f"{sample}\t{mds}\n")
                
        except Exception as e:
            logger.error(f"Unified Extract+Motif failed: {e}")
            raise typer.Exit(1)


    # ═══════════════════════════════════════════════════════════════════
    # 3. UNIFIED SINGLE-PASS PIPELINE (FSC, FSR, FSD, WPS, OCF)
    # ═══════════════════════════════════════════════════════════════════
    logger.info("🚀 Running Unified Single-Pass Pipeline (FSC, FSR, FSD, WPS, OCF)...")

    # Define Resource Paths (Resolve defaults if None)
    pkg_dir = Path(__file__).parent
    
    # FSC/FSR Bins
    res_bin = bin_input if bin_input else pkg_dir / "data" / "ChormosomeBins" / "hg19_window_100kb.bed"
    if not res_bin.exists(): logger.error(f"Bin file missing: {res_bin}"); raise typer.Exit(1)

    # FSD Arms
    res_arms = arms_file if arms_file else pkg_dir / "data" / "ChormosomeArms" / "hg19.arms.bed"
    if not res_arms.exists(): logger.error(f"Arms file missing: {res_arms}"); raise typer.Exit(1)

    # WPS Genes
    res_wps = wps_file if wps_file else pkg_dir / "data" / "TranscriptAnno" / "transcriptAnno-hg19-1kb.tsv"
    if not res_wps.exists(): logger.error(f"Transcript file missing: {res_wps}"); raise typer.Exit(1)

    # OCF Regions
    res_ocf = ocr_file if ocr_file else pkg_dir / "data" / "OpenChromatinRegion" / "7specificTissue.all.OC.bed"
    if not res_ocf.exists(): logger.error(f"OCR file missing: {res_ocf}"); raise typer.Exit(1)

    # Define Outputs
    out_fsc_raw = output / f"{sample}.fsc_counts.tsv"
    out_wps = output / f"{sample}.WPS.tsv.gz"
    out_fsd = output / f"{sample}.FSD.tsv"
    out_ocf_dir = output / f"{sample}_ocf_tmp" # OCF writes mult files to dir
    out_ocf_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # RUN RUST PIPELINE
        _core.run_unified_pipeline(
            str(bedgz_file),
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
            
            # FSC Output
            final_fsc = output / f"{sample}.FSC.tsv"
            process_fsc_from_counts(df_counts, final_fsc, fsc_windows, fsc_continue_n)
            
            # FSR Output
            final_fsr = output / f"{sample}.FSR.tsv"
            process_fsr_from_counts(df_counts, final_fsr, fsc_windows, fsc_continue_n)
    
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
