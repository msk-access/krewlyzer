"""
Fragment extraction from BAM files.

Extracts cfDNA fragments with full filter control.
Uses Rust backend for accelerated extraction.
"""

import typer
from pathlib import Path
from typing import Optional
import logging
import pysam
import json
from datetime import datetime

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("extract")

# Rust backend is required
from krewlyzer import _core


def extract(
    bam_input: Path = typer.Argument(..., help="Input BAM file (sorted, indexed)"),
    genome_reference: Path = typer.Option(..., '-g', '--genome', help="Reference genome FASTA (indexed)"),
    output: Path = typer.Option(..., '-o', '--output', help="Output directory"),
    
    # Configurable filters
    exclude_regions: Optional[Path] = typer.Option(None, '-x', '--exclude-regions', help="Exclude regions BED file"),
    mapq: int = typer.Option(20, '--mapq', '-q', help="Minimum mapping quality"),
    minlen: int = typer.Option(65, '--minlen', help="Minimum fragment length"),
    maxlen: int = typer.Option(400, '--maxlen', help="Maximum fragment length"),
    skip_duplicates: bool = typer.Option(True, '--skip-duplicates/--no-skip-duplicates', help="Skip duplicate reads"),
    require_proper_pair: bool = typer.Option(True, '--require-proper-pair/--no-require-proper-pair', help="Require proper pairs"),
    
    # Other options
    chromosomes: Optional[str] = typer.Option(None, '--chromosomes', help="Comma-separated chromosomes to process"),
    sample_name: Optional[str] = typer.Option(None, '--sample-name', '-s', help="Sample name for output files (default: derived from BAM filename)"),
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
    """
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
    
    # Default exclude regions
    if exclude_regions is None:
        pkg_dir = Path(__file__).parent
        default_exclude = pkg_dir / "data" / "exclude-regions" / "hg19-blacklist.v2.bed"
        if default_exclude.exists():
            exclude_regions = default_exclude
            logger.info(f"Using default exclude regions: {exclude_regions}")
    elif not exclude_regions.exists():
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
    
    try:
        logger.info(f"Extracting fragments from {bam_input.name}")
        logger.info(f"Filters: mapq>={mapq}, length=[{minlen},{maxlen}], skip_dup={skip_duplicates}, proper_pair={require_proper_pair}")
        
        # Call Unified Rust Engine (Extract Mode)
        # Returns (fragment_count, em_counts, bpm_counts)
        # We only care about fragment_count and the side-effect (BED file)
        
        fragment_count, _, _ = _core.extract_motif.process_bam_parallel(
            str(bam_input),
            str(genome_reference),
            mapq,
            minlen,
            maxlen,
            4, # kmer (dummy, not used if motif output is None?)
               # Actually kmer=4 is default. Engine needs it?
               # process_bam_parallel signature: (bam, fasta, mapq, min, max, kmer, threads, bed, motif, exclude, skip, proper)
            threads,
            str(bed_temp),         # output_bed_path
            None,                  # output_motif_prefix (None = skip motif counting)
            str(exclude_regions) if exclude_regions else None,
            skip_duplicates,
            require_proper_pair
        )
        
        logger.info(f"Extracted {fragment_count:,} fragments")
        
        # Compress with BGZF and index with tabix
        logger.info("Compressing and indexing...")
        pysam.tabix_compress(str(bed_temp), str(bed_output), force=True)
        pysam.tabix_index(str(bed_output), preset="bed", force=True)
        
        # Remove temp file
        import os
        os.remove(str(bed_temp))
        
        # Write metadata
        metadata = {
            "sample_id": sample_name,
            "total_fragments": fragment_count,
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

    except Exception as e:
        logger.error(f"Fragment extraction failed: {e}")
        raise typer.Exit(1)

