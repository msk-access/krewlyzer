"""Command-line interface for cfDNAFE"""

from typing import Optional
import typer
from pathlib import Path
from rich.console import Console
from rich.logging import RichHandler
import logging

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("krewlyzer-cli")

def set_log_level(log_level: str = typer.Option("INFO", "--log-level", help="Logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL")):
    """Set global logging level."""
    level = getattr(logging, log_level.upper(), logging.INFO)
    for handler in logging.root.handlers:
        handler.setLevel(level)
    logging.getLogger().setLevel(level)

app = typer.Typer(help="krewlyzer: A comprehensive toolkit for ctDNA fragmentomics analysis.")

from krewlyzer.motif import motif
from krewlyzer.extract import extract
from krewlyzer.fsc import fsc
from krewlyzer.fsr import fsr
from krewlyzer.fsd import fsd
from krewlyzer.wps import wps
from krewlyzer.ocf import ocf
from krewlyzer.uxm import uxm
from krewlyzer.mfsd import mfsd
from krewlyzer.region_entropy import region_entropy
from krewlyzer.wrapper import run_all
from krewlyzer.pon.build import build_pon
from krewlyzer.build_gc_reference import build_gc_reference
from krewlyzer import __version__

# -------------------------------------------------------------------
# VALIDATE command: Standalone asset validation for user-provided files
# -------------------------------------------------------------------
def validate(
    gene_bed: Path = typer.Option(None, "--gene-bed", help="Gene BED file to validate (chrom, start, end, gene, [name])"),
    targets_bed: Path = typer.Option(None, "--targets-bed", help="Targets BED file to validate"),
    arms_bed: Path = typer.Option(None, "--arms-bed", help="Chromosome arms BED file to validate (chrom, start, end, arm)"),
    wps_anchors: Path = typer.Option(None, "--wps-anchors", help="WPS anchors BED6 file to validate"),
    ocr_file: Path = typer.Option(None, "--ocr-file", help="OCF regions BED file to validate"),
    gc_factors: Path = typer.Option(None, "--gc-factors", help="GC correction factors TSV to validate"),
    bin_file: Path = typer.Option(None, "--bin-file", help="Bin BED3 file to validate (chrom, start, end)"),
):
    """
    Validate asset file formats before running analysis.
    
    Checks that user-provided override files match expected formats.
    Bundled assets (from data/) are known-good and not validated here.
    
    Examples:
        krewlyzer validate --gene-bed my_genes.bed
        krewlyzer validate --arms-bed my_arms.bed --wps-anchors my_anchors.bed
    """
    from krewlyzer.core.asset_validation import validate_file, FileSchema
    
    validated_count = 0
    
    # Validate each provided file
    if gene_bed:
        validate_file(gene_bed, FileSchema.GENE_BED)
        validated_count += 1
    if targets_bed:
        validate_file(targets_bed, FileSchema.TARGETS_BED)
        validated_count += 1
    if arms_bed:
        validate_file(arms_bed, FileSchema.ARMS_BED)
        validated_count += 1
    if wps_anchors:
        validate_file(wps_anchors, FileSchema.WPS_ANCHORS)
        validated_count += 1
    if ocr_file:
        validate_file(ocr_file, FileSchema.REGION_BED)
        validated_count += 1
    if gc_factors:
        validate_file(gc_factors, FileSchema.GC_FACTORS_TSV)
        validated_count += 1
    if bin_file:
        validate_file(bin_file, FileSchema.BED3)
        validated_count += 1
    
    if validated_count == 0:
        console.print("[yellow]No files specified. Use --help to see options.[/yellow]")
        raise typer.Exit(1)
    
    console.print(f"[green]✓ All {validated_count} files validated successfully[/green]")

app.command()(extract)
app.command()(motif)
app.command()(fsc)
app.command()(fsr)
app.command()(fsd)
app.command()(wps)
app.command()(ocf)
app.command()(uxm)
app.command()(mfsd)
app.command(name="region-entropy")(region_entropy)
app.command()(run_all)
app.command(name="build-pon")(build_pon)
app.command(name="build-gc-reference")(build_gc_reference)
app.command()(validate)

@app.callback()
def main(
    version: bool = typer.Option(False, "--version", "-v", help="Show version and exit"),
):
    if version:
        typer.echo(f"krewlyzer version: {__version__}")
        raise typer.Exit()

if __name__ == "__main__":
    app()
