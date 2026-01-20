"""
WPS Anchor Filtering for Panel Mode.

Filters genome-wide WPS anchors to panel-specific genes for
targeted WPS analysis in MSK-ACCESS and similar panels.
"""

import gzip
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

logger = logging.getLogger(__name__)


def filter_anchors_by_genes(
    anchors_bed: Path,
    gene_names: Set[str],
    output_path: Path,
    include_ctcf_near_genes: bool = True,
    ctcf_proximity_bp: int = 100000
) -> int:
    """
    Filter WPS anchors to only those related to specified genes.
    
    Keeps:
    - TSS anchors whose gene name is in gene_names
    - CTCF anchors within ctcf_proximity_bp of any gene region (if include_ctcf_near_genes)
    
    Args:
        anchors_bed: Path to genome-wide WPS anchors BED (gzipped or plain)
        gene_names: Set of gene names to keep
        output_path: Path to write filtered anchors (gzip compressed)
        include_ctcf_near_genes: Include nearby CTCF sites (default True)
        ctcf_proximity_bp: Distance from gene regions for CTCF filtering
        
    Returns:
        Number of anchors written
    """
    logger.info(f"Filtering {anchors_bed.name} to {len(gene_names)} genes")
    
    # Normalize gene names to uppercase for matching
    gene_names_upper = {g.upper() for g in gene_names}
    
    # First pass: identify TSS anchors and collect their locations for CTCF filtering
    tss_locations = []  # [(chrom, pos), ...]
    open_func = gzip.open if str(anchors_bed).endswith('.gz') else open
    
    with open_func(anchors_bed, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            
            name = fields[3]
            if name.startswith('TSS|'):
                # Format: TSS|GENENAME|ENST...
                parts = name.split('|')
                if len(parts) >= 2:
                    gene = parts[1].upper()
                    if gene in gene_names_upper:
                        chrom = fields[0]
                        pos = int(fields[1])
                        tss_locations.append((chrom, pos))
    
    logger.debug(f"Found {len(tss_locations)} TSS anchors for target genes")
    
    # Build proximity index for CTCF filtering
    tss_by_chrom = {}
    for chrom, pos in tss_locations:
        if chrom not in tss_by_chrom:
            tss_by_chrom[chrom] = []
        tss_by_chrom[chrom].append(pos)
    
    for chrom in tss_by_chrom:
        tss_by_chrom[chrom].sort()
    
    # Second pass: write filtered anchors
    count = 0
    with open_func(anchors_bed, 'rt') as f_in:
        with gzip.open(output_path, 'wt') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    f_out.write(line)
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 4:
                    continue
                
                name = fields[3]
                keep = False
                
                if name.startswith('TSS|'):
                    # TSS anchor - check gene name
                    parts = name.split('|')
                    if len(parts) >= 2:
                        gene = parts[1].upper()
                        if gene in gene_names_upper:
                            keep = True
                
                elif name.startswith('CTCF|') and include_ctcf_near_genes:
                    # CTCF anchor - check proximity to any TSS
                    chrom = fields[0]
                    pos = int(fields[1])
                    
                    if chrom in tss_by_chrom:
                        # Check if any TSS is within proximity
                        for tss_pos in tss_by_chrom[chrom]:
                            if abs(pos - tss_pos) <= ctcf_proximity_bp:
                                keep = True
                                break
                
                if keep:
                    f_out.write(line.strip() + '\n')
                    count += 1
    
    logger.info(f"Wrote {count} anchors to {output_path.name}")
    return count


def generate_panel_wps_anchors(
    genome_anchors: Path,
    genes: Dict[str, List],
    output_path: Path
) -> int:
    """
    Generate panel-specific WPS anchors from gene dictionary.
    
    Convenience wrapper around filter_anchors_by_genes.
    
    Args:
        genome_anchors: Path to genome-wide WPS anchors BED
        genes: Dict from load_gene_bed() - Dict[gene_name, List[Region]]
        output_path: Path to write filtered anchors
        
    Returns:
        Number of anchors written
    """
    gene_names = set(genes.keys())
    return filter_anchors_by_genes(
        anchors_bed=genome_anchors,
        gene_names=gene_names,
        output_path=output_path
    )


def get_bundled_wps_anchors(assay: str, genome: str = "GRCh37") -> Optional[Path]:
    """
    Get path to bundled panel-specific WPS anchors.
    
    Args:
        assay: Assay code (xs1, xs2)
        genome: Genome build (GRCh37, GRCh38)
        
    Returns:
        Path to bundled WPS anchors or None if not available
    """
    import importlib.resources
    
    try:
        with importlib.resources.files("krewlyzer.data.WpsAnchors").joinpath(genome) as genome_dir:
            anchors_path = genome_dir / f"{assay}.wps_anchors.bed.gz"
            if anchors_path.exists():
                return Path(anchors_path)
    except Exception as e:
        logger.debug(f"Bundled WPS anchors lookup failed for {assay}/{genome}: {e}")
    
    return None
