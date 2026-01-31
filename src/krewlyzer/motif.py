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

from rich.console import Console
from rich.logging import RichHandler

console = Console(stderr=True)
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console, show_time=True, show_path=False)], format="%(message)s")
logger = logging.getLogger("motif")

# Rust backend is required
from krewlyzer import _core

# Import asset resolution and startup banner
from .core.asset_resolution import resolve_target_regions, resolve_pon_model
from .core.logging import log_startup_banner, ResolvedAsset
from . import __version__


def motif(
    bam_input: Path = typer.Option(..., "--input", "-i", help="Input BAM file"),
    genome_reference: Path = typer.Option(..., '-r', '--reference', help="Reference genome FASTA (indexed)"),
    output: Path = typer.Option(..., '-o', '--output', help="Output directory"),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/GRCh37/hg38/GRCh38)"),
    target_regions: Optional[Path] = typer.Option(None, '-T', '--target-regions', help="Target regions BED (for panel data: generates on/off-target motifs)"),
    skip_target_regions: bool = typer.Option(False, "--skip-target-regions", help="Disable panel mode even when --assay has bundled targets"),
    assay: Optional[str] = typer.Option(None, "--assay", "-A", help="Assay code (xs1/xs2) for auto-loading bundled assets"),
    pon_model: Optional[Path] = typer.Option(None, '-P', '--pon-model', help="PON model for MDS z-score computation"),
    pon_variant: str = typer.Option("all_unique", "--pon-variant", help="PON variant: 'all_unique' (default, max coverage) or 'duplex' (highest accuracy)"),
    skip_pon: bool = typer.Option(False, "--skip-pon", help="Skip PON z-score normalization"),
    kmer: int = typer.Option(4, '-k', '--kmer', help="K-mer size for motif extraction"),
    chromosomes: Optional[str] = typer.Option(None, '--chromosomes', help="Comma-separated chromosomes to process"),
    sample_name: Optional[str] = typer.Option(None, '--sample-name', '-s', help="Sample name for output files (default: derived from BAM filename)"),
    require_proper_pair: bool = typer.Option(True, '--require-proper-pair/--no-require-proper-pair', help="Require proper pairs (disable for duplex/consensus BAMs)"),
    verbose: bool = typer.Option(False, '--verbose', '-v', help="Enable verbose logging"),
    threads: int = typer.Option(0, '--threads', '-t', help="Number of threads (0=all cores)"),
):
    """
    Extract k-mer motif features from a BAM file.
    
    Output:
    - {sample}.EndMotif: End motif k-mer frequencies
    - {sample}.BreakPointMotif: Breakpoint motif k-mer frequencies
    - {sample}.MDS: Motif Diversity Score
    
    Note: For fragment extraction (BED.gz), use `krewlyzer extract` instead.
    """
    from .assets import AssetManager
    
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
    
    # Initialize AssetManager
    assets = AssetManager(genome)
    
    # ═══════════════════════════════════════════════════════════════════
    # ASSET RESOLUTION
    # ═══════════════════════════════════════════════════════════════════
    try:
        resolved_pon_path, pon_source = resolve_pon_model(
            explicit_path=pon_model, assay=assay, skip_pon=skip_pon,
            assets=assets, variant=pon_variant, log=logger
        )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)
    
    try:
        resolved_target_path, target_source = resolve_target_regions(
            explicit_path=target_regions, assay=assay, skip_target_regions=skip_target_regions,
            assets=assets, log=logger
        )
    except ValueError as e:
        console.print(f"[bold red]❌ ERROR:[/bold red] {e}")
        raise typer.Exit(1)
    
    is_panel_mode = resolved_target_path is not None and resolved_target_path.exists()
    
    # ═══════════════════════════════════════════════════════════════════
    # STARTUP BANNER
    # ═══════════════════════════════════════════════════════════════════
    log_startup_banner(
        tool_name="motif", version=__version__,
        inputs={"BAM": str(bam_input.name), "Reference": str(genome_reference.name), "Output": str(output)},
        config={"Genome": f"{assets.raw_genome} → {assets.genome_dir}", "Assay": assay or "None", "Mode": "Panel" if is_panel_mode else "WGS", "K-mer": str(kmer)},
        assets=[ResolvedAsset("PON", resolved_pon_path, pon_source), ResolvedAsset("Targets", resolved_target_path, target_source)],
        logger=logger
    )
    
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
    
    try:
        logger.info(f"Extracting motif features from {bam_input.name}")
        
        # Use unified extract_sample() for extraction
        from .core.sample_processor import extract_sample, write_motif_outputs
        
        result = extract_sample(
            input_path=bam_input,
            reference=genome_reference,
            mapq=20,  # Default mapQ for motif
            minlen=65,
            maxlen=1000,
            kmer=kmer,
            threads=threads,
            skip_duplicates=True,
            require_proper_pair=require_proper_pair,
            target_regions=resolved_target_path,
            exclude_regions=None,
            write_bed=False,  # Motif command doesn't write BED
            output_dir=None,
            sample_name=sample_name,
            genome=genome,
        )
        
        logger.info(f"Processed {result.fragment_count:,} fragments")
        
        # Load PON model if resolved and not skipped
        pon = None
        if resolved_pon_path and resolved_pon_path.exists():
            try:
                from .pon.model import PonModel
                pon = PonModel.load(resolved_pon_path)
                logger.debug(f"Loaded PON model for z-score computation")
            except Exception as e:
                logger.warning(f"Could not load PON model: {e}")
        
        # Write motif outputs using shared function
        outputs = write_motif_outputs(
            result=result,
            output_dir=output,
            pon=pon,
            include_ontarget=is_panel_mode,
        )
        
        # Additional aberrant k-mer analysis if PON provided
        if pon and hasattr(pon, 'mds_baseline') and pon.mds_baseline:
            try:
                aberrant = pon.mds_baseline.get_aberrant_kmers(result.kmer_frequencies, threshold=3.0)
                if aberrant:
                    logger.info(f"Found {len(aberrant)} aberrant k-mers (|z| > 3.0)")
            except Exception as e:
                logger.debug(f"Could not compute aberrant k-mers: {e}")
        
        logger.info(f"Motif extraction complete (MDS={result.mds_score:.4f})")

    except Exception as e:
        logger.error(f"Motif extraction failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise typer.Exit(1)
