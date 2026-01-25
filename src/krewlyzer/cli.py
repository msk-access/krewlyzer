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

@app.callback()
def main(
    version: bool = typer.Option(False, "--version", "-v", help="Show version and exit"),
):
    if version:
        typer.echo(f"krewlyzer version: {__version__}")
        raise typer.Exit()

if __name__ == "__main__":
    app()
