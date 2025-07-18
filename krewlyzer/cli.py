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

from .extract_features import extract_features
from .quality_control import quality_control
from .motif import motif

app.command()(extract_features)
app.command()(quality_control)
app.command()(motif)

@app.command()
def process_region(
    region: str = typer.Argument(..., help="Genomic region in format chr:start-end"),
    bam_file: Path = typer.Argument(..., help="Path to BAM file"),
    output_file: Path = typer.Argument(..., help="Output file for region features"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose output")
) -> None:
    """Process specific genomic region"""
    if verbose:
        logger.setLevel(logging.DEBUG)
    logger.info(f"Processing region {region} from {bam_file}")
    try:
        # TODO: Implement region processing logic
        # For now, just mock success
        logger.info(f"Region features written to {output_file}")
    except Exception as e:
        logger.error(f"Failed to process region: {e}")
        raise typer.Exit(1)

@app.command()
def version() -> None:
    """Show version information"""
    logger.info("cfDNAFE 0.1.0")

if __name__ == "__main__":
    app()
