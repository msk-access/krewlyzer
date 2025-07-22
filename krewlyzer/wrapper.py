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

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("krewlyzer-wrapper")


def run_all(
    bam_file: Path = typer.Argument(..., help="Input BAM file (sorted, indexed)"),
    reference: Path = typer.Option(..., "--reference", "-g", help="Reference genome FASTA file for motif extraction"),
    output: Path = typer.Option(..., "--output", "-o", help="Output directory for all results"),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of parallel processes for each step"),
    pe_type: str = typer.Option("SE", "--type", help="Fragment type for UXM: SE or PE (default: SE)")
):
    """
    Run all feature extraction commands (motif, fsc, fsr, fsd, wps, ocf, uxm) for a single BAM file.
    """
    # Input checks
    if not bam_file.exists() or not bam_file.is_file():
        logger.error(f"Input BAM file not found: {bam_file}")
        raise typer.Exit(1)
    if not reference.exists() or not reference.is_file():
        logger.error(f"Reference FASTA file not found: {reference}")
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
            reference=reference,
            output=motif_output,
            minlen=65,
            maxlen=400,
            k=3,
            verbose=True,
            threads=threads,
        )
    except Exception as e:
        logger.error(f"Motif extraction failed: {e}")
        raise typer.Exit(1)
    # 2. FSC
    fsc_output = output / "fsc"
    try:
        fsc(
            bedgz_path=motif_output,
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
            output=fsr_output,
            threads=threads
        )
    except Exception as e:
        logger.error(f"FSR calculation failed: {e}")
        raise typer.Exit(1)
    # 4. FSD
    fsd_output = output / "fsd"
    try:
        fsd(
            bedgz_path=motif_output,
            output=fsd_output,
            arms_file=None,
            threads=threads
        )
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
    logger.info(f"All feature extraction complete. Results saved to {output}")
