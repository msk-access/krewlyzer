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

def set_log_level(log_level: str = typer.Option("INFO", "--log-level", help="Logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL")):
    """Set global logging level."""
    level = getattr(logging, log_level.upper(), logging.INFO)
    for handler in logging.root.handlers:
        handler.setLevel(level)
    logging.getLogger().setLevel(level)

app = typer.Typer(callback=set_log_level)

from .motif import motif
from .fsc import fsc
from .fsr import fsr
from .fsd import fsd
from .wps import wps
from .ocf import ocf
from .uxm import uxm
from .wrapper import run_all

app.command()(motif)
app.command()(fsc)
app.command()(fsr)
app.command()(fsd)
app.command()(wps)
app.command()(ocf)
app.command()(uxm)
app.command()(run_all)

@app.command()
def version() -> None:
    """Show version information"""
    logger.info("krewlyzer 0.1.0")

if __name__ == "__main__":
    app()
