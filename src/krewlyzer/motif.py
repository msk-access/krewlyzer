"""
Motif feature extraction from BAM files.

Extracts End Motif (EDM), Breakpoint Motif (BPM), and Motif Diversity Score (MDS)
from cfDNA fragment endpoints. These features capture tissue-of-origin signals.

Uses Rust backend for accelerated extraction.

Output Files:
    - {sample}.EndMotif.tsv: End motif 4-mer frequencies (256 columns)
    - {sample}.BreakPointMotif.tsv: Breakpoint motif frequencies
    - {sample}.MDS.tsv: Motif Diversity Score summary

Panel Mode: Generates separate on-target (.ontarget.tsv) files.

Note: For fragment extraction (BED.gz), use `krewlyzer extract` instead.
"""

import typer
from pathlib import Path
from typing import Optional
import logging
import numpy as np
import itertools

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("motif")

# Rust backend is required
from krewlyzer import _core


def motif(
    bam_input: Path = typer.Option(..., "--input", "-i", help="Input BAM file"),
    genome_reference: Path = typer.Option(..., '-r', '--reference', help="Reference genome FASTA (indexed)"),
    output: Path = typer.Option(..., '-o', '--output', help="Output directory"),
    target_regions: Optional[Path] = typer.Option(None, '-T', '--target-regions', help="Target regions BED (for panel data: generates on/off-target motifs)"),
    kmer: int = typer.Option(4, '-k', '--kmer', help="K-mer size for motif extraction"),
    chromosomes: Optional[str] = typer.Option(None, '--chromosomes', help="Comma-separated chromosomes to process"),
    sample_name: Optional[str] = typer.Option(None, '--sample-name', '-s', help="Sample name for output files (default: derived from BAM filename)"),
    require_proper_pair: bool = typer.Option(True, '--require-proper-pair/--no-require-proper-pair', help="Require proper pairs (disable for duplex/consensus BAMs)"),
    verbose: bool = typer.Option(False, '--verbose', '-v', help="Enable verbose logging"),
    threads: int = typer.Option(0, '--threads', '-t', help="Number of threads (0=all cores)"),
    format: Optional[str] = typer.Option(None, "--format", "-f", help="Output format override: tsv, parquet, json (default: tsv)")
):
    """
    Extract k-mer motif features from a BAM file.
    
    Output:
    - {sample}.EndMotif: End motif k-mer frequencies
    - {sample}.BreakPointMotif: Breakpoint motif k-mer frequencies
    - {sample}.MDS: Motif Diversity Score
    
    Note: For fragment extraction (BED.gz), use `krewlyzer extract` instead.
    """
    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("core.motif_processor").setLevel(logging.DEBUG)
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
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name (use provided or derive from BAM filename)
    if sample_name is None:
        sample_name = bam_input.stem.replace('.bam', '')
    
    # Pre-check BAM compatibility with current filters
    if require_proper_pair:
        from .core.bam_utils import check_bam_compatibility
        logger.info("Checking BAM read compatibility with filters...")
        compat = check_bam_compatibility(bam_input, require_proper_pair, True, 20)
        
        if compat["pass_rate"] < 0.01 and compat["total_sampled"] > 100:
            # Critical: Almost no reads would pass!
            console.print("\n[bold red]⚠️  FILTER COMPATIBILITY WARNING[/bold red]\n")
            console.print(f"Only [bold]{compat['pass_rate']:.2%}[/bold] of sampled reads would pass current filters.\n")
            
            for issue in compat["issues"]:
                console.print(f"  • {issue}")
            
            if compat["suggested_flags"]:
                suggested = " ".join(compat["suggested_flags"])
                console.print(f"\n[bold yellow]Suggested command:[/bold yellow]")
                console.print(f"  krewlyzer motif {bam_input.name} -r {genome_reference.name} -o {output} [bold]{suggested}[/bold]\n")
                
                logger.error(f"Filter mismatch detected. Re-run with: {suggested}")
                raise typer.Exit(1)
    
    # Output file paths
    edm_output = output / f"{sample_name}.EndMotif.tsv"
    bpm_output = output / f"{sample_name}.BreakPointMotif.tsv"
    mds_output = output / f"{sample_name}.MDS.tsv"
    
    try:
        logger.info(f"Extracting motif features from {bam_input.name}")
        
        # Initialize motif dictionaries
        bases = ['A', 'C', 'T', 'G']
        End_motif = {''.join(i): 0 for i in itertools.product(bases, repeat=kmer)}
        Breakpoint_motif = {''.join(i): 0 for i in itertools.product(bases, repeat=kmer)}
        
        # Parse chromosomes
        chroms = chromosomes.split(',') if chromosomes else None
        
        # Call Unified Rust Engine (Extract + Motif)
        # Returns: (fragment_count, em_counts, bpm_counts, gc_obs, em_counts_on, bpm_counts_on)
        
        # Call Rust extraction - returns 7 values: (count, em_off, bpm_off, gc_obs, em_on, bpm_on, gc_obs_on)
        fragment_count, em_counts, bpm_counts, _gc_obs, em_counts_on, bpm_counts_on, _gc_obs_on = _core.extract_motif.process_bam_parallel(
            str(bam_input),
            str(genome_reference),
            20,    # Default mapQ
            65,    # Default min length
            400,   # Default max length
            kmer,
            threads,
            None,  # output_bed_path
            "enable",  # output_motif_prefix (triggers counting)
            None,  # exclude_path
            str(target_regions) if target_regions else None,  # target_regions_path
            True,  # skip_duplicates
            require_proper_pair,  # proper pair filter
            False  # silent
        )
        
        End_motif.update(em_counts)
        Breakpoint_motif.update(bpm_counts)
        logger.info(f"Processed {fragment_count:,} fragments")
        
        # Write all motif outputs using shared processor
        from .core.motif_processor import process_motif_outputs
        
        total_em, total_bpm, mds = process_motif_outputs(
            em_counts=End_motif,
            bpm_counts=Breakpoint_motif,
            edm_output=edm_output,
            bpm_output=bpm_output,
            mds_output=mds_output,
            sample_name=sample_name,
            kmer=kmer,
            include_headers=True  # Consistent with run-all
        )
        
        # Process on-target motifs if target_regions was provided and we have data
        is_panel_mode = target_regions and target_regions.exists()
        if is_panel_mode and sum(em_counts_on.values()) > 0:
            edm_on = output / f"{sample_name}.EndMotif.ontarget.tsv"
            bpm_on = output / f"{sample_name}.BreakPointMotif.ontarget.tsv"
            mds_on = output / f"{sample_name}.MDS.ontarget.tsv"
            
            End_motif_on = {''.join(i): 0 for i in itertools.product(bases, repeat=kmer)}
            Breakpoint_motif_on = {''.join(i): 0 for i in itertools.product(bases, repeat=kmer)}
            End_motif_on.update(em_counts_on)
            Breakpoint_motif_on.update(bpm_counts_on)
            
            total_em_on, total_bpm_on, mds_on_val = process_motif_outputs(
                em_counts=End_motif_on,
                bpm_counts=Breakpoint_motif_on,
                edm_output=edm_on,
                bpm_output=bpm_on,
                mds_output=mds_on,
                sample_name=sample_name,
                kmer=kmer,
                include_headers=True
            )
            logger.info(f"On-target motifs: EM={total_em_on:,}, BPM={total_bpm_on:,}")
        
        logger.info(f"Motif extraction complete (EM={total_em:,}, BPM={total_bpm:,}, MDS={mds:.4f})")

    except Exception as e:
        logger.error(f"Motif extraction failed: {e}")
        raise typer.Exit(1)
