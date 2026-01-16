"""
Fragment extraction from BAM files.

Extracts cfDNA fragments with configurable read filters (MAPQ, length, proper pair).
First step in the pipeline - generates .bed.gz for downstream feature extraction.

Uses Rust backend for accelerated extraction.

Output Files:
    - {sample}.bed.gz: Fragment coordinates (chrom, start, end, strand, gc%, gc_weight)
    - {sample}.correction_factors.csv: LOESS-based GC correction factors
    - {sample}.metadata.json: Processing statistics

GC Correction: Computes per-fragment GC weights for all downstream tools.
Panel Mode: GC model trained on off-target reads only for unbiased correction.

BAM Compatibility: Auto-detects duplex BAMs and warns to use --no-require-proper-pair.
"""

import typer
from pathlib import Path
from typing import Optional
import logging
import pysam
import json
import os
from datetime import datetime

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("extract")

# Rust backend is required
from krewlyzer import _core

# Shared BAM utilities
from .core.bam_utils import check_bam_compatibility


def extract(
    bam_input: Path = typer.Option(..., "--input", "-i", help="Input BAM file (sorted, indexed)"),
    genome_reference: Path = typer.Option(..., '-r', '--reference', help="Reference genome FASTA (indexed)"),
    output: Path = typer.Option(..., '-o', '--output', help="Output directory"),
    
    # GC correction options
    genome: str = typer.Option("hg19", '-G', '--genome', help="Genome build for GC correction assets (hg19/GRCh37/hg38/GRCh38)"),
    gc_correct: bool = typer.Option(True, '--gc-correct/--no-gc-correct', help="Compute GC correction factors"),
    
    # Configurable filters
    exclude_regions: Optional[Path] = typer.Option(None, '-x', '--exclude-regions', help="Exclude regions BED file"),
    target_regions: Optional[Path] = typer.Option(None, '-T', '--target-regions', help="Target regions BED (for panel data: GC model computed on off-target reads only)"),
    mapq: int = typer.Option(20, '--mapq', '-q', help="Minimum mapping quality"),
    minlen: int = typer.Option(65, '--minlen', help="Minimum fragment length"),
    maxlen: int = typer.Option(400, '--maxlen', help="Maximum fragment length"),
    skip_duplicates: bool = typer.Option(True, '--skip-duplicates/--no-skip-duplicates', help="Skip duplicate reads"),
    require_proper_pair: bool = typer.Option(True, '--require-proper-pair/--no-require-proper-pair', help="Require proper pairs"),
    
    # Other options
    chromosomes: Optional[str] = typer.Option(None, '--chromosomes', help="Comma-separated chromosomes to process"),
    sample_name: Optional[str] = typer.Option(None, '--sample-name', '-s', help="Sample name for output files (default: derived from BAM filename)"),
    verbose: bool = typer.Option(False, '--verbose', '-v', help="Enable verbose logging"),
    threads: int = typer.Option(0, '--threads', '-t', help="Number of threads (0=all cores)")
):
    """
    Extract cfDNA fragments from BAM to BED.gz with full filter control.
    
    Filters applied (always on):
    - Skip unmapped reads
    - Skip secondary alignments  
    - Skip supplementary alignments
    - Skip failed QC reads
    
    Filters applied (configurable):
    - Mapping quality threshold
    - Fragment length range
    - Skip duplicates
    - Require proper pairs
    - Exclude genomic regions
    
    Output:
    - {sample}.bed.gz: Fragment coordinates with GC content
    - {sample}.bed.gz.tbi: Tabix index
    - {sample}.metadata.json: Fragment count and metadata
    - {sample}.correction_factors.csv: GC correction factors (if --gc-correct)
    """
    from .assets import AssetManager
    
    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")
    
    # Configure Rust thread pool
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
    
    if not str(bam_input).endswith('.bam'):
        logger.error(f"Input must be a .bam file: {bam_input}")
        raise typer.Exit(1)
    
    if not genome_reference.exists():
        logger.error(f"Reference genome not found: {genome_reference}")
        raise typer.Exit(1)
    
    # Initialize Asset Manager for GC correction
    assets = None
    if gc_correct:
        try:
            assets = AssetManager(genome)
            logger.info(f"Genome: {assets.raw_genome} -> {assets.genome_dir}")
        except ValueError as e:
            logger.warning(f"AssetManager error: {e}. GC correction disabled.")
            gc_correct = False
    
    # Default exclude regions from assets
    if exclude_regions is None and assets:
        try:
            exclude_regions = assets.resolve("exclude_regions")
            logger.info(f"Using default exclude regions: {exclude_regions}")
        except FileNotFoundError:
            # Fallback to hardcoded path
            pkg_dir = Path(__file__).parent
            default_exclude = pkg_dir / "data" / "exclude-regions" / "hg19-blacklist.v2.bed.gz"
            if default_exclude.exists():
                exclude_regions = default_exclude
                logger.info(f"Using default exclude regions: {exclude_regions}")
    elif exclude_regions and not exclude_regions.exists():
        logger.error(f"Exclude regions file not found: {exclude_regions}")
        raise typer.Exit(1)
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from BAM filename)
    if sample_name is None:
        sample_name = bam_input.stem.replace('.bam', '')
    
    # Parse chromosomes
    chrom_list = chromosomes.split(',') if chromosomes else None
    
    # Output paths
    bed_temp = output / f"{sample_name}.bed"  # Temp uncompressed
    bed_output = output / f"{sample_name}.bed.gz"
    metadata_output = output / f"{sample_name}.metadata.json"
    factors_output = output / f"{sample_name}.correction_factors.csv"
    
    try:
        # Pre-check BAM compatibility with current filters
        logger.info("Checking BAM read compatibility with filters...")
        compat = check_bam_compatibility(bam_input, require_proper_pair, skip_duplicates, mapq)
        
        if compat["pass_rate"] < 0.01 and compat["total_sampled"] > 100:
            # Critical: Almost no reads would pass!
            console.print("\n[bold red]⚠️  FILTER COMPATIBILITY WARNING[/bold red]\n")
            console.print(f"Only [bold]{compat['pass_rate']:.2%}[/bold] of sampled reads would pass current filters.\n")
            
            for issue in compat["issues"]:
                console.print(f"  • {issue}")
            
            if compat["suggested_flags"]:
                suggested = " ".join(compat["suggested_flags"])
                console.print(f"\n[bold yellow]Suggested command:[/bold yellow]")
                console.print(f"  krewlyzer extract {bam_input.name} -r {genome_reference.name} -o {output} [bold]{suggested}[/bold]\n")
                
                logger.error(f"Filter mismatch detected. Re-run with: {suggested}")
                raise typer.Exit(1)
        elif compat["pass_rate"] < 0.5 and compat["issues"]:
            # Warning but continue
            console.print("\n[yellow]⚠️  Filter compatibility note:[/yellow]")
            for issue in compat["issues"]:
                console.print(f"  • {issue}")
            console.print()
        
        logger.info(f"Extracting fragments from {bam_input.name}")
        logger.info(f"Filters: mapq>={mapq}, length=[{minlen},{maxlen}], skip_dup={skip_duplicates}, proper_pair={require_proper_pair}")
        
        if target_regions:
            logger.info(f"Panel mode: GC model will use off-target reads only (targets: {target_regions.name})")
        
        # Call Unified Rust Engine (Extract Mode)
        # Returns (fragment_count, em_counts, bpm_counts, gc_observations, em_counts_on, bpm_counts_on, gc_observations_ontarget)
        fragment_count, _, _, gc_observations, _, _, gc_observations_ontarget = _core.extract_motif.process_bam_parallel(
            str(bam_input),
            str(genome_reference),
            mapq,
            minlen,
            maxlen,
            4,  # kmer
            threads,
            str(bed_temp),         # output_bed_path
            None,                  # output_motif_prefix (None = skip motif counting)
            str(exclude_regions) if exclude_regions else None,
            str(target_regions) if target_regions else None,  # For off-target GC model
            skip_duplicates,
            require_proper_pair
        )
        
        logger.info(f"Extracted {fragment_count:,} fragments")
        
        # Rust writes BGZF directly - just need to index with tabix
        # The temp file is already .bed.gz format from noodles::bgzf
        import shutil
        logger.info("Moving and indexing BED.gz...")
        shutil.move(str(bed_temp), str(bed_output))
        pysam.tabix_index(str(bed_output), preset="bed", force=True)
        
        # Compute GC correction factors inline (no second BED pass!)
        if gc_correct and assets:
            try:
                gc_ref = assets.resolve("gc_reference")
                valid_regions = assets.resolve("valid_regions")
                
                # Off-target GC correction (primary - always computed)
                logger.info(f"Computing GC correction factors from {len(gc_observations)} observation bins...")
                n_factors = _core.gc.compute_and_write_gc_factors(
                    gc_observations,
                    str(gc_ref),
                    str(valid_regions),
                    str(factors_output)
                )
                logger.info(f"Computed {n_factors} correction factors: {factors_output}")
                
                # On-target GC correction (panel mode only)
                if target_regions and len(gc_observations_ontarget) > 0:
                    factors_ontarget = output / f"{sample_name}.correction_factors.ontarget.csv"
                    logger.info(f"Computing ON-TARGET GC correction factors from {len(gc_observations_ontarget)} observation bins...")
                    n_factors_on = _core.gc.compute_and_write_gc_factors(
                        gc_observations_ontarget,
                        str(gc_ref),
                        str(valid_regions),
                        str(factors_ontarget)
                    )
                    logger.info(f"Computed {n_factors_on} on-target correction factors: {factors_ontarget}")
                
            except FileNotFoundError as e:
                logger.warning(f"GC correction assets not found: {e}")
                logger.warning("Skipping GC correction. Run 'krewlyzer build-gc-reference' to generate assets.")
                gc_correct = False
            except Exception as e:
                logger.warning(f"GC correction failed: {e}")
                gc_correct = False
        
        # Write metadata
        metadata = {
            "sample_id": sample_name,
            "total_fragments": fragment_count,
            "genome": genome,
            "gc_correction_computed": gc_correct,
            "panel_mode": target_regions is not None,
            "target_regions": str(target_regions) if target_regions else None,
            "filters": {
                "mapq": mapq,
                "min_length": minlen,
                "max_length": maxlen,
                "skip_duplicates": skip_duplicates,
                "require_proper_pair": require_proper_pair
            },
            "timestamp": datetime.now().isoformat()
        }
        with open(metadata_output, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Output: {bed_output}")
        logger.info(f"Metadata: {metadata_output}")
        if gc_correct:
            logger.info(f"Correction factors: {factors_output}")

    except Exception as e:
        logger.error(f"Fragment extraction failed: {e}")
        raise typer.Exit(1)
