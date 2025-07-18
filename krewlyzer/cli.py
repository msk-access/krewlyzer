"""Command-line interface for cfDNAFE"""

from typing import Optional
import typer
from pathlib import Path
from rich.console import Console
from rich.logging import RichHandler
import logging

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("krewlyzer-cli")

app = typer.Typer()

from .motif import motif
from .fsc import fsc
from .fsr import fsr
from .fsd import fsd

app.command()(motif)
app.command()(fsc)
app.command()(fsr)
app.command()(fsd)

@app.command()
def version() -> None:
    """Show version information"""
    logger.info("krewlyzer 0.1.0")

if __name__ == "__main__":
    app()
