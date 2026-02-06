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

app = typer.Typer(
    help="krewlyzer: A comprehensive toolkit for ctDNA fragmentomics analysis.",
    invoke_without_command=True,
)

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
from krewlyzer.region_mds import region_mds
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
    genome: str = typer.Option(None, "--genome", "-G", help="Validate ALL bundled assets for genome (hg19/hg38)"),
    assay: str = typer.Option(None, "--assay", "-A", help="Also validate assay-specific bundled assets (requires --genome)"),
):
    """
    Validate asset file formats before running analysis.
    
    Use with file options to validate user-provided override files.
    Use with --genome to validate bundled data assets.
    
    Examples:
        krewlyzer validate --gene-bed my_genes.bed
        krewlyzer validate --genome hg19
        krewlyzer validate --genome hg19 --assay xs2
    """
    from krewlyzer.core.asset_validation import validate_file, FileSchema
    
    validated_count = 0
    
    # Validate bundled assets if --genome specified
    if genome:
        from krewlyzer.assets import AssetManager
        assets = AssetManager(genome)
        results = assets.validate(assay=assay)
        validated_count += sum(1 for v in results.values() if v is True)
    
    # Validate each provided user file
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
        console.print("[yellow]No files specified. Use --genome or file options.[/yellow]")
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
app.command(name="region-mds")(region_mds)
app.command()(run_all)
app.command(name="build-pon")(build_pon)
app.command(name="build-gc-reference")(build_gc_reference)
app.command()(validate)

@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    version: bool = typer.Option(False, "--version", "-v", help="Show version and exit"),
):
    """Krewlyzer: ctDNA fragmentomics analysis toolkit."""
    if version:
        typer.echo(f"krewlyzer {__version__}")
        raise typer.Exit()
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()

if __name__ == "__main__":
    app()
