"""
Per-Region MDS (Motif Diversity Score) CLI Tool

Calculates MDS at the exon/target level from BAM files.
Outputs:
- {sample}.MDS.exon.tsv: Per-exon/target MDS scores
- {sample}.MDS.gene.tsv: Gene-level aggregated MDS scores

Based on Helzer et al. (2025) methodology.
"""

from pathlib import Path
from typing import Optional
import typer
import logging

logger = logging.getLogger("krewlyzer.region_mds")


def region_mds(
    bam_input: Path = typer.Argument(..., help="Input BAM file (indexed)"),
    reference: Path = typer.Argument(..., help="Reference FASTA file"),
    output: Path = typer.Argument(..., help="Output directory"),
    gene_bed: Optional[Path] = typer.Option(
        None, "--gene-bed", "-g", help="Gene/exon BED file (panel or WGS format)"
    ),
    genome: str = typer.Option("hg19", "--genome", "-G", help="Genome build (hg19/hg38)"),
    assay: Optional[str] = typer.Option(
        None, "--assay", "-a", help="Assay code (xs1, xs2, wgs) for bundled gene BED"
    ),
    e1_only: bool = typer.Option(False, "--e1-only", help="Only output E1 (first exon) results"),
    mapq: int = typer.Option(20, "--mapq", help="Minimum mapping quality"),
    minlen: int = typer.Option(65, "--minlen", help="Minimum fragment length"),
    maxlen: int = typer.Option(1000, "--maxlen", help="Maximum fragment length (default: 1000 for consistency with run-all)"),
    pon_model: Optional[Path] = typer.Option(None, "--pon-model", "-P", help="PON model for z-score computation"),
    skip_pon: bool = typer.Option(False, "--skip-pon", help="Skip PON z-score normalization (for PON samples used as ML negatives)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging"),
    silent: bool = typer.Option(False, "--silent", help="Suppress progress bar"),
):
    """
    Calculate per-region Motif Diversity Score (MDS).
    
    MDS measures the diversity of 4-mer end motifs at each genomic region.
    Higher MDS indicates more diverse motif usage; lower MDS may indicate
    aberrant fragmentation patterns associated with disease.
    
    Examples:
    
        # Panel mode with bundled gene BED
        krewlyzer region-mds sample.bam ref.fa out/ --assay xs2
        
        # WGS mode
        krewlyzer region-mds sample.bam ref.fa out/ --assay wgs
        
        # Custom gene BED
        krewlyzer region-mds sample.bam ref.fa out/ --gene-bed custom.bed
    """
    from krewlyzer import _core
    from krewlyzer.assets import AssetManager
    
    # Configure verbose logging
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate inputs
    if not bam_input.exists():
        raise typer.BadParameter(f"BAM file not found: {bam_input}")
    if not reference.exists():
        raise typer.BadParameter(f"Reference FASTA not found: {reference}")
    
    # Validate user-provided override files
    from .core.asset_validation import validate_file, FileSchema
    if gene_bed and gene_bed.exists():
        logger.debug(f"Validating user-provided gene BED: {gene_bed}")
        validate_file(gene_bed, FileSchema.GENE_BED)
    
    # Resolve gene BED
    assets = AssetManager(genome)
    resolved_gene_bed = None
    
    if gene_bed:
        # User-provided gene BED
        if not gene_bed.exists():
            raise typer.BadParameter(f"Gene BED not found: {gene_bed}")
        resolved_gene_bed = gene_bed
        logger.info(f"Using user-provided gene BED: {gene_bed}")
    elif assay:
        # Bundled gene BED from assay
        resolved_gene_bed = assets.get_gene_bed_for_mode(assay)
        if resolved_gene_bed is None:
            raise typer.BadParameter(f"Gene BED not available for assay '{assay}' and genome '{genome}'")
        logger.info(f"Using bundled gene BED for assay '{assay}': {resolved_gene_bed}")
    else:
        # Default to WGS if no assay specified
        if assets.wgs_gene_bed_available:
            resolved_gene_bed = assets.wgs_gene_bed
            logger.info(f"Using WGS gene BED: {resolved_gene_bed}")
        else:
            raise typer.BadParameter(
                "No gene BED specified. Use --gene-bed, --assay, or ensure wgs.genes.bed.gz is available."
            )
    
    # Create output directory
    output.mkdir(parents=True, exist_ok=True)
    
    # Derive sample name from BAM filename
    sample_name = bam_input.stem
    if sample_name.endswith(".sorted"):
        sample_name = sample_name[:-7]
    
    output_exon = output / f"{sample_name}.MDS.exon.tsv"
    output_gene = output / f"{sample_name}.MDS.gene.tsv"
    
    logger.info(f"Running region-MDS analysis")
    logger.info(f"  BAM: {bam_input}")
    logger.info(f"  Gene BED: {resolved_gene_bed}")
    logger.info(f"  Output: {output_exon}, {output_gene}")
    logger.info(f"  Filters: mapq>={mapq}, length=[{minlen},{maxlen}]")
    
    # Call Rust engine
    n_regions, n_genes = _core.region_mds.run_region_mds(
        str(bam_input),
        str(reference),
        str(resolved_gene_bed),
        str(output_exon),
        str(output_gene),
        e1_only,
        mapq,
        minlen,
        maxlen,
        silent,
    )
    
    logger.info(f"Region-MDS complete: {n_regions} regions, {n_genes} genes")
    logger.info(f"  Exon output: {output_exon}")
    logger.info(f"  Gene output: {output_gene}")
    
    # PON z-score normalization
    if pon_model and not skip_pon:
        try:
            import pandas as pd
            from krewlyzer.pon.model import PonModel
            
            logger.info(f"Applying PON z-score normalization: {pon_model}")
            pon = PonModel.load(pon_model)
            
            if pon.region_mds is None:
                logger.warning("PON model does not have region_mds baseline - skipping z-scores")
            else:
                # Read gene output
                df = pd.read_csv(output_gene, sep="\t")
                
                # Compute z-scores for each gene
                mds_z_scores = []
                mds_e1_z_scores = []
                
                for _, row in df.iterrows():
                    gene = row["gene"]
                    mds_mean = row.get("mds_mean", 0.0)
                    mds_e1 = row.get("mds_e1", 0.0)
                    
                    # MDS mean z-score
                    z = pon.region_mds.compute_zscore(gene, mds_mean)
                    mds_z_scores.append(z if z is not None else float('nan'))
                    
                    # E1 z-score
                    z_e1 = pon.region_mds.compute_e1_zscore(gene, mds_e1)
                    mds_e1_z_scores.append(z_e1 if z_e1 is not None else float('nan'))
                
                # Add z-score columns
                df["mds_z"] = mds_z_scores
                df["mds_e1_z"] = mds_e1_z_scores
                
                # Write updated output
                df.to_csv(output_gene, sep="\t", index=False)
                
                n_with_z = sum(1 for z in mds_z_scores if z is not None and not pd.isna(z))
                logger.info(f"Added z-scores for {n_with_z}/{len(df)} genes")
                
        except Exception as e:
            logger.warning(f"PON z-score computation failed: {e}")
    elif skip_pon:
        logger.debug("Skipping PON z-score normalization (--skip-pon)")
    
    return n_regions, n_genes
