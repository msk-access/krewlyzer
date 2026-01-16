"""
Build GC Reference Assets CLI.

Generates pre-computed GC correction assets:
1. valid_regions.bed - Curated 100kb bins excluding problematic regions
2. ref_genome_GC.parquet - Expected fragment counts per (length, GC)

These assets are generated once per reference genome and shipped with krewlyzer.
"""

import typer
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("build_gc_reference")

# Rust backend is required
from krewlyzer import _core


def build_gc_reference(
    reference: Path = typer.Argument(..., help="Path to reference FASTA file"),
    output_dir: Path = typer.Option(..., "--output", "-o", help="Output directory for GC assets"),
    exclude_regions: Optional[Path] = typer.Option(None, "--exclude-regions", "-e", help="Path to exclude regions BED file"),
    target_regions: Optional[Path] = typer.Option(None, "--target-regions", "-T", help="Target regions BED for panel mode (generates ontarget files)"),
    bin_size: int = typer.Option(100000, "--bin-size", help="Size of each bin in bp (default: 100kb)"),
    genome_name: Optional[str] = typer.Option(None, "--genome-name", "-n", help="Genome name (default: derived from reference filename)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose/debug logging"),
):
    """
    Build GC reference assets for a reference genome.
    
    This command generates two files:
    1. valid_regions_{genome}.bed - Curated bins for GC estimation
    2. ref_genome_GC_{genome}.parquet - Expected fragment counts per (length, GC)
    
    For panel mode (--target-regions), also generates:
    3. valid_regions_{genome}.ontarget.bed - Bins intersecting target regions
    4. ref_genome_GC_{genome}.ontarget.parquet - Expected counts within targets
    
    Example:
        krewlyzer build-gc-reference hg38.fa -o data/gc/ -e exclude-regions.bed
        krewlyzer build-gc-reference hg38.fa -o data/gc/ -T msk_targets.bed
    """
    # Set log level based on verbose flag
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")
    
    # Input validation
    if not reference.exists():
        logger.error(f"Reference FASTA not found: {reference}")
        raise typer.Exit(1)
    
    # Check for .fai index
    fai_path = Path(str(reference) + ".fai")
    if not fai_path.exists():
        logger.error(f"Reference index not found: {fai_path}")
        logger.error("Run 'samtools faidx <reference>' to create the index")
        raise typer.Exit(1)
    
    # Default exclude regions from existing data
    if exclude_regions is None:
        pkg_dir = Path(__file__).parent
        # Try existing exclude regions based on genome name hint
        ref_name_lower = reference.name.lower()
        if 'hg38' in ref_name_lower or 'grch38' in ref_name_lower:
            exclude_regions = pkg_dir / "data" / "exclude-regions" / "hg38-blacklist.v2.bed.gz"
        else:
            # Default to hg19 for hg19/b37/GRCh37
            exclude_regions = pkg_dir / "data" / "exclude-regions" / "hg19-blacklist.v2.bed.gz"
        
        if exclude_regions.exists():
            logger.info(f"Using default exclude regions: {exclude_regions}")
        else:
            logger.warning(f"No exclude regions file found at {exclude_regions}")
            logger.warning("No exclude regions file provided or found. Regions may include problematic areas.")
            # Create empty exclude regions file
            exclude_regions = output_dir / "empty_exclude_regions.bed"
            output_dir.mkdir(parents=True, exist_ok=True)
            exclude_regions.write_text("# Empty exclude regions\n")
    
    if not exclude_regions.exists():
        logger.error(f"Exclude regions file not found: {exclude_regions}")
        raise typer.Exit(1)
    
    # Derive genome name
    if genome_name is None:
        genome_name = reference.stem.replace(".fa", "").replace(".fasta", "")
        # Common aliases
        if "GRCh38" in genome_name or "hg38" in genome_name:
            genome_name = "hg38"
        elif "GRCh37" in genome_name or "hg19" in genome_name or "b37" in genome_name:
            genome_name = "hg19"
        logger.info(f"Derived genome name: {genome_name}")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Output paths
    valid_regions_path = output_dir / f"valid_regions_{genome_name}.bed.gz"
    ref_gc_path = output_dir / f"ref_genome_GC_{genome_name}.parquet"
    
    try:
        # Step 1: Generate valid regions
        logger.info(f"[1/2] Generating valid regions...")
        logger.info(f"  Reference: {reference}")
        logger.info(f"  Exclude regions: {exclude_regions}")
        logger.info(f"  Bin size: {bin_size:,} bp")
        
        region_count = _core.gc.generate_valid_regions(
            str(reference),
            str(exclude_regions),
            str(valid_regions_path),
            bin_size,
        )
        
        logger.info(f"  Generated {region_count:,} valid regions: {valid_regions_path}")
        
        # Step 2: Generate ref_genome_GC
        logger.info(f"[2/2] Generating ref_genome_GC...")
        logger.info(f"  This may take several minutes for whole genome...")
        
        fragment_count = _core.gc.generate_ref_genome_gc(
            str(reference),
            str(valid_regions_path),
            str(ref_gc_path),
        )
        
        logger.info(f"  Counted {fragment_count:,} theoretical fragments: {ref_gc_path}")
        
        logger.info("GC reference assets generated successfully!")
        logger.info(f"  Valid regions: {valid_regions_path}")
        logger.info(f"  Ref genome GC: {ref_gc_path}")
        
        # Step 3: Generate on-target files if target_regions provided (panel mode)
        if target_regions:
            if not target_regions.exists():
                logger.error(f"Target regions file not found: {target_regions}")
                raise typer.Exit(1)
            
            logger.info("")
            logger.info("Panel mode: Generating on-target GC reference assets...")
            logger.info(f"  Target regions: {target_regions}")
            
            # On-target output paths
            valid_regions_ontarget_path = output_dir / f"valid_regions_{genome_name}.ontarget.bed.gz"
            ref_gc_ontarget_path = output_dir / f"ref_genome_GC_{genome_name}.ontarget.parquet"
            
            # Step 3a: Generate on-target valid regions (intersection with target regions)
            logger.info(f"[3/4] Generating on-target valid regions...")
            
            ontarget_region_count = _core.gc.generate_valid_regions_ontarget(
                str(valid_regions_path),
                str(target_regions),
                str(valid_regions_ontarget_path),
            )
            
            logger.info(f"  Generated {ontarget_region_count:,} on-target valid regions: {valid_regions_ontarget_path}")
            
            # Step 3b: Generate on-target ref_genome_GC
            logger.info(f"[4/4] Generating on-target ref_genome_GC...")
            
            ontarget_fragment_count = _core.gc.generate_ref_genome_gc(
                str(reference),
                str(valid_regions_ontarget_path),
                str(ref_gc_ontarget_path),
            )
            
            logger.info(f"  Counted {ontarget_fragment_count:,} on-target theoretical fragments: {ref_gc_ontarget_path}")
            
            logger.info("")
            logger.info("On-target GC reference assets generated successfully!")
            logger.info(f"  Valid regions (on-target): {valid_regions_ontarget_path}")
            logger.info(f"  Ref genome GC (on-target): {ref_gc_ontarget_path}")
        
    except Exception as e:
        logger.error(f"Failed to build GC reference assets: {e}")
        raise typer.Exit(1)

