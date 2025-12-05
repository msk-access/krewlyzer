import typer
from pathlib import Path
import logging
from rich.console import Console
from rich.logging import RichHandler
from .motif import motif
from .fsc import fsc
from .fsr import fsr
from .fsd import fsd
from .wps import wps
from .ocf import ocf
from .uxm import uxm
from .mfsd import mfsd
from typing import Optional

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("krewlyzer-wrapper")


def run_all(
    bam_file: Path = typer.Argument(...,help="Input BAM file (sorted, indexed)"),
    reference: Path = typer.Option(..., "--reference", "-g", help="Reference genome FASTA file"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory for all results"),
    
    # Common filtering options
    chromosomes: Optional[str] = typer.Option(None, "--chromosomes", help="Comma-separated list of chromosomes (e.g., '1,2,3')"),
    mapq: int = typer.Option(20, "--mapq", "-q", help="Minimum mapping quality"),
    exclude_regions: Optional[Path] = typer.Option(None, "--exclude-regions", "-x", help="BED file of genomic regions to exclude"),
    minlen: int = typer.Option(65, "--minlen", help="Minimum fragment length"),
    maxlen: int = typer.Option(400, "--maxlen", help="Maximum fragment length"),
    
    # Optional override for FSD arms file
    arms_file: Optional[Path] = typer.Option(None, "--arms-file", "-a", help="Custom arms/regions file for FSD (default: project data/ChormosomeArms/hg19.arms.bed)"),
    
    # Existing options
    variant_input: Optional[Path] = typer.Option(None, "--variant-input", "-v", help="Input VCF/MAF file for mFSD"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of parallel processes"),
    pe_type: str = typer.Option("SE", "--type", help="Fragment type for UXM: SE or PE")
):
    """
    Run all feature extraction commands (motif, fsc, fsr, fsd, wps, ocf, uxm, mfsd) for a single BAM file.
    
    Common filtering options (chromosomes, mapq, exclude-regions, fragment lengths) are exposed.
    Advanced options can be configured by running tools individually.
    """
    # Input checks
    if not bam_file.exists() or not bam_file.is_file():
        logger.error(f"Input BAM file not found: {bam_file}")
        raise typer.Exit(1)
    if not reference.exists() or not reference.is_file():
        logger.error(f"Reference FASTA file not found: {reference}")
        raise typer.Exit(1)
    
    # Validate optional filtering files
    if exclude_regions and not exclude_regions.exists():
        logger.error(f"Exclude regions file not found: {exclude_regions}")
        raise typer.Exit(1)
    if arms_file and not arms_file.exists():
        logger.error(f"Arms file not found: {arms_file}")
        raise typer.Exit(1)
    
    # Validate fragment length ranges
    if minlen >= maxlen:
        logger.error(f"Minimum fragment length ({minlen}) must be less than maximum ({maxlen})")
        raise typer.Exit(1)
    
    try:
        output.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Could not create output directory {output}: {e}")
        raise typer.Exit(1)
    # 1. Motif extraction
    motif_output = output / "motif"
    try:
        motif(
            bam_path=bam_file,
            genome_reference=reference,      # Fixed: was reference=
            output=motif_output,
            blacklist=exclude_regions,       # Fixed: use exclude_regions parameter
            map_quality=mapq,                # Added: expose mapq
            min_length=minlen,               # Fixed: was minlen=
            max_length=maxlen,               # Fixed: was maxlen=
            kmer=3,                          # Fixed: was k=
            chromosomes=chromosomes,         # Added: pass as string (motif splits internally)
            verbose=False,                   # Changed: less noisy
            threads=threads
        )
    except Exception as e:
        logger.error(f"Motif extraction failed: {e}")
        raise typer.Exit(1)
    # 2. FSC
    fsc_output = output / "fsc"
    try:
        fsc(
            bedgz_path=motif_output,
            bin_input=None,  # Use default packaged bin file
            output=fsc_output,
            threads=threads
        )
    except Exception as e:
        logger.error(f"FSC calculation failed: {e}")
        raise typer.Exit(1)
    # 3. FSR
    fsr_output = output / "fsr"
    try:
        fsr(
            bedgz_path=motif_output,
            bin_input=None,  # Use default packaged bin file
            output=fsr_output,
            threads=threads
        )
    except Exception as e:
        logger.error(f"FSR calculation failed: {e}")
        raise typer.Exit(1)
    # 4. FSD
    fsd_output = output / "fsd"
    try:
        # Use packaged default if user doesn't provide arms file
        # Data is at project root: data/ChormosomeArms/hg19.arms.bed
        if arms_file is None:
            project_root = Path(__file__).parent.parent
            arms_file = project_root / "data" / "ChormosomeArms" / "hg19.arms.bed"
            if not arms_file.exists():
                logger.warning(f"Default arms file not found: {arms_file}. Skipping FSD.")
                raise FileNotFoundError(f"Arms file not found: {arms_file}")
        
        fsd(
            bedgz_path=motif_output,
            output=fsd_output,
            arms_file=arms_file,  # Fixed: was None
            threads=threads
        )
    except FileNotFoundError:
        logger.warning("FSD skipped due to missing arms file")
    except Exception as e:
        logger.error(f"FSD calculation failed: {e}")
        raise typer.Exit(1)
    # 5. WPS
    wps_output = output / "wps"
    try:
        wps(
            bedgz_path=motif_output,
            output=wps_output,
            threads=threads
        )
    except Exception as e:
        logger.error(f"WPS calculation failed: {e}")
        raise typer.Exit(1)
    # 6. OCF
    ocf_output = output / "ocf"
    try:
        ocf(
            bedgz_path=motif_output,
            output=ocf_output,
            threads=threads
        )
    except Exception as e:
        logger.error(f"OCF calculation failed: {e}")
        raise typer.Exit(1)
    # 7. UXM
    uxm_output = output / "uxm"
    try:
        uxm(
            bam_path=bam_file.parent,
            output=uxm_output,
            pe_type=pe_type,
            threads=threads
        )
    except Exception as e:
        logger.error(f"UXM calculation failed: {e}")
        raise typer.Exit(1)
    
    # 8. mFSD (Optional)
    if variant_input:
        if not variant_input.exists():
            logger.warning(f"Variant input file not found: {variant_input}. Skipping mFSD.")
        else:
            mfsd_output = output / "mfsd" / (bam_file.stem + ".mfsd.tsv")
            try:
                mfsd(
                    bam_path=bam_file,
                    input_file=variant_input,
                    output=mfsd_output,
                    format="auto",
                    map_quality=20
                )
            except Exception as e:
                logger.error(f"mFSD calculation failed: {e}")
                # Don't raise exit here, just log error as it's optional
    
    logger.info(f"All feature extraction complete. Results saved to {output}")
