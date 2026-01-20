"""
Gene BED Parser for MSK-ACCESS Panel Data.

Parses MSK-ACCESS target BED files and groups regions by gene.
Supports V1 (xs1) and V2 (xs2) naming conventions.

V1 format: exon_GENE_exon_probe (e.g., exon_MTOR_48a.1_2)
V2 format: GENE_target_NN (e.g., MTOR_target_02)

Non-gene regions (FP_, msi_, snp_FP_, BAT-*) are filtered out.
"""

import gzip
import re
import logging
from pathlib import Path
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Dict, List, Optional, IO

logger = logging.getLogger(__name__)


@contextmanager
def _open_bed(path: Path) -> IO[str]:
    """Open a BED file, handling gzip compression transparently."""
    path = Path(path)
    if str(path).endswith('.gz'):
        f = gzip.open(path, 'rt')
    else:
        f = open(path, 'r')
    try:
        yield f
    finally:
        f.close()


@dataclass
class Region:
    """A genomic region with gene annotation."""
    chrom: str
    start: int
    end: int
    name: str
    gene: str


def detect_assay(target_bed: Path) -> Optional[str]:
    """
    Detect assay type from target BED filename.
    
    Looks for 'xs1' or 'xs2' in filename, or legacy 'v1.0'/'v2.0' patterns.
    
    Args:
        target_bed: Path to target regions BED file
        
    Returns:
        Assay code ('xs1' or 'xs2') or None if not detected
        
    Examples:
        >>> detect_assay(Path("xs2.targets.bed"))
        'xs2'
        >>> detect_assay(Path("MSK-ACCESS-v1_0-probe-A.sorted.bed"))
        'xs1'
        >>> detect_assay(Path("custom_panel.bed"))
        None
    """
    filename = target_bed.name.lower()
    
    # Direct match for standardized names
    if 'xs1' in filename:
        logger.debug(f"Detected assay 'xs1' from filename: {filename}")
        return 'xs1'
    if 'xs2' in filename:
        logger.debug(f"Detected assay 'xs2' from filename: {filename}")
        return 'xs2'
    
    # Legacy patterns
    if 'v1.0' in filename or 'v1_0' in filename or 'access-v1' in filename.replace('_', '-'):
        logger.debug(f"Detected assay 'xs1' from legacy pattern: {filename}")
        return 'xs1'
    if 'v2.0' in filename or 'v2_0' in filename or 'access-v2' in filename.replace('_', '-'):
        logger.debug(f"Detected assay 'xs2' from legacy pattern: {filename}")
        return 'xs2'
    
    logger.debug(f"Could not detect assay from filename: {filename}")
    return None


def detect_version_from_content(bed_path: Path) -> str:
    """
    Detect BED version (V1, V2, or GENE_BED) by examining content.
    
    V1 uses: exon_GENE_exon_probe pattern
    V2 uses: GENE_target_NN pattern
    GENE_BED: Pre-parsed gene BED with #chrom\\tstart\\tend\\tgene\\tname header
    
    Args:
        bed_path: Path to target BED file
        
    Returns:
        'V1', 'V2', or 'GENE_BED' based on content analysis
    """
    v1_pattern = re.compile(r'^exon_[A-Z0-9]+_\d+')
    v2_pattern = re.compile(r'^[A-Z0-9]+_target_\d+$')
    
    v1_count = 0
    v2_count = 0
    is_gene_bed = False
    
    with _open_bed(bed_path) as f:
        for i, line in enumerate(f):
            if i >= 100:  # Sample first 100 lines
                break
            
            # Check for gene BED header
            if line.startswith('#chrom'):
                if 'gene' in line.lower():
                    is_gene_bed = True
                    break
                continue
            
            if line.startswith('#') or not line.strip():
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            
            # V1 format has name in column 5 (strand in column 4)
            # V2 format has name in column 4
            # Check both columns for pattern matching
            for idx in [4, 3]:
                if len(fields) > idx:
                    name = fields[idx]
                    if v1_pattern.match(name):
                        v1_count += 1
                        break
                    elif v2_pattern.match(name):
                        v2_count += 1
                        break
    
    if is_gene_bed:
        logger.debug(f"Detected GENE_BED format (pre-parsed)")
        return 'GENE_BED'
    
    version = 'V1' if v1_count > v2_count else 'V2'
    logger.debug(f"Detected version {version} (V1: {v1_count}, V2: {v2_count})")
    return version


def _is_gene_region(name: str) -> bool:
    """
    Check if region name represents a gene target (not fingerprint/MSI).
    
    Excludes:
      - FP_* (fingerprint probes)
      - snp_FP_* (SNP fingerprint probes)
      - msi_* (microsatellite instability markers)
      - BAT-* (BAT markers)
    
    Args:
        name: Region name from BED file
        
    Returns:
        True if this is a gene region, False otherwise
    """
    name_lower = name.lower()
    
    # Exclude patterns
    exclude_prefixes = ('fp_', 'snp_fp_', 'msi_', 'bat-')
    
    for prefix in exclude_prefixes:
        if name_lower.startswith(prefix):
            return False
    
    return True


def _extract_gene_v1(name: str) -> Optional[str]:
    """
    Extract gene name from V1 format: exon_GENE_exon_probe
    
    Examples:
        exon_MTOR_48a.1_2 -> MTOR
        exon_ARID1A_1a.1_1 -> ARID1A
        exon_H3F3A_2_1 -> H3F3A
    """
    # Pattern: exon_GENE_... 
    match = re.match(r'^exon_([A-Z0-9]+)_', name, re.IGNORECASE)
    if match:
        return match.group(1).upper()
    return None


def _extract_gene_v2(name: str) -> Optional[str]:
    """
    Extract gene name from V2 format: GENE_target_NN
    
    Examples:
        MTOR_target_02 -> MTOR
        ARID1A_target_01 -> ARID1A
    """
    # Pattern: GENE_target_NN
    match = re.match(r'^([A-Z0-9]+)_target_\d+$', name, re.IGNORECASE)
    if match:
        return match.group(1).upper()
    return None


def parse_gene_bed(
    bed_path: Path, 
    version: str = "auto"
) -> Dict[str, List[Region]]:
    """
    Parse MSK-ACCESS target BED file and group regions by gene.
    
    Filters out non-gene regions (FP, MSI, BAT markers) and groups
    the remaining regions by their associated gene name.
    
    Args:
        bed_path: Path to target BED file
        version: 'V1', 'V2', or 'auto' for auto-detection
        
    Returns:
        Dictionary mapping gene names to lists of Region objects
        
    Example:
        >>> genes = parse_gene_bed(Path("xs2.targets.bed"))
        >>> genes['MTOR']
        [Region(chrom='1', start=11168235, end=11168345, ...)]
    """
    bed_path = Path(bed_path)
    
    if not bed_path.exists():
        raise FileNotFoundError(f"BED file not found: {bed_path}")
    
    # Auto-detect version if needed
    if version == "auto":
        version = detect_version_from_content(bed_path)
        logger.info(f"Auto-detected BED format: {version}")
    
    genes: Dict[str, List[Region]] = {}
    total_regions = 0
    filtered_regions = 0
    
    with _open_bed(bed_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            
            total_regions += 1
            
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            
            # Handle different formats
            if version == 'GENE_BED':
                # Pre-parsed format: chrom, start, end, gene, name
                gene = fields[3] if len(fields) > 3 else None
                name = fields[4] if len(fields) > 4 else gene or ""
            else:
                # Get region name - V1 has name in column 5, V2 in column 4
                name = ""
                if len(fields) >= 5:
                    # Check if column 5 looks like a gene name (V1 format)
                    if fields[4] and (fields[4].startswith('exon_') or 'target' in fields[3].lower()):
                        name = fields[4] if fields[4].startswith('exon_') else fields[3]
                    else:
                        name = fields[3]
                elif len(fields) >= 4:
                    name = fields[3]
                
                # Filter non-gene regions
                if not _is_gene_region(name):
                    filtered_regions += 1
                    continue
                
                # Extract gene name based on version
                extract_gene = _extract_gene_v1 if version == 'V1' else _extract_gene_v2
                gene = extract_gene(name)
                if not gene:
                    # Fallback: try other version's pattern
                    alt_extract = _extract_gene_v2 if version == 'V1' else _extract_gene_v1
                    gene = alt_extract(name)
            
            if not gene:
                logger.debug(f"Could not extract gene from: {name}")
                continue
            
            region = Region(
                chrom=chrom,
                start=start,
                end=end,
                name=name,
                gene=gene
            )
            
            if gene not in genes:
                genes[gene] = []
            genes[gene].append(region)
    
    logger.info(
        f"Parsed {len(genes)} genes from {total_regions} regions "
        f"(filtered {filtered_regions} non-gene regions)"
    )
    
    return genes


def get_gene_count(genes: Dict[str, List[Region]]) -> Dict[str, int]:
    """Get count of regions per gene."""
    return {gene: len(regions) for gene, regions in genes.items()}


def write_gene_bed(
    genes: Dict[str, List[Region]], 
    output_path: Path,
    sort: bool = True,
    compress: bool = True
) -> int:
    """
    Write gene-grouped regions to a BED file (optionally BGZF compressed).
    
    Args:
        genes: Dictionary from parse_gene_bed()
        output_path: Output path for gene BED file
        sort: If True, sort by chromosome and position
        compress: If True, write as BGZF (.bed.gz)
        
    Returns:
        Number of regions written
    """
    import gzip
    
    output_path = Path(output_path)
    
    # Collect all regions
    all_regions = []
    for gene, regions in genes.items():
        for region in regions:
            all_regions.append((region.chrom, region.start, region.end, gene, region.name))
    
    # Sort by chromosome and position
    if sort:
        # Custom chromosome sort (1, 2, ..., 22, X, Y)
        def chrom_key(r):
            chrom = r[0].replace('chr', '')
            if chrom.isdigit():
                return (0, int(chrom), r[1])
            return (1, chrom, r[1])
        
        all_regions.sort(key=chrom_key)
    
    # Build content
    lines = ["#chrom\tstart\tend\tgene\tname\n"]
    for chrom, start, end, gene, name in all_regions:
        lines.append(f"{chrom}\t{start}\t{end}\t{gene}\t{name}\n")
    content = "".join(lines)
    
    # Write output (compressed or plain)
    if compress or str(output_path).endswith('.gz'):
        with gzip.open(output_path, 'wt') as f:
            f.write(content)
    else:
        with open(output_path, 'w') as f:
            f.write(content)
    
    logger.info(f"Wrote {len(all_regions)} gene regions to {output_path}")
    return len(all_regions)


def get_bundled_gene_bed(assay: str, genome: str = "GRCh37") -> Optional[Path]:
    """
    Get path to bundled gene BED file for a known assay.
    
    Args:
        assay: Assay code ('xs1' or 'xs2')
        genome: Genome build ('GRCh37' or 'GRCh38')
        
    Returns:
        Path to bundled gene BED file, or None if not found
        
    Example:
        >>> path = get_bundled_gene_bed('xs2', 'GRCh37')
        >>> path
        PosixPath('.../data/genes/GRCh37/xs2.genes.bed.gz')
    """
    import importlib.resources
    
    assay = assay.lower()
    if assay not in ('xs1', 'xs2'):
        logger.debug(f"Unknown assay '{assay}', no bundled gene BED")
        return None
    
    # Use importlib.resources for package data
    try:
        data_pkg = f"krewlyzer.data.genes.{genome}"
        filename = f"{assay}.genes.bed.gz"
        
        # Try to get the file path
        ref = importlib.resources.files(data_pkg) / filename
        if ref.is_file():
            # Return actual path for downstream compatibility
            return Path(str(ref))
    except (ModuleNotFoundError, FileNotFoundError, TypeError):
        pass
    
    # Fallback: check relative to this file
    data_dir = Path(__file__).parent.parent / "data" / "genes" / genome
    bundled_path = data_dir / f"{assay}.genes.bed.gz"
    
    if bundled_path.exists():
        logger.debug(f"Found bundled gene BED: {bundled_path}")
        return bundled_path
    
    logger.debug(f"Bundled gene BED not found for assay={assay}, genome={genome}")
    return None


def load_gene_bed(
    target_bed: Optional[Path] = None,
    gene_bed: Optional[Path] = None,
    assay: Optional[str] = None,
    genome: str = "GRCh37"
) -> Dict[str, List[Region]]:
    """
    Load gene regions from appropriate source.
    
    Priority:
    1. Explicit gene_bed file if provided
    2. Bundled gene BED if assay is xs1/xs2
    3. Parse target_bed on-the-fly
    
    Args:
        target_bed: Path to target regions BED (for on-the-fly parsing)
        gene_bed: Explicit path to pre-parsed gene BED
        assay: Assay code ('xs1', 'xs2', or None for auto-detect)
        genome: Genome build for bundled asset lookup
        
    Returns:
        Dictionary mapping gene names to Region lists
        
    Raises:
        ValueError: If no valid source is available
    """
    # Option 1: Explicit gene BED file
    if gene_bed and Path(gene_bed).exists():
        logger.info(f"Loading gene BED from explicit path: {gene_bed}")
        return parse_gene_bed(gene_bed)
    
    # Option 2: Auto-detect or use provided assay for bundled file
    if assay is None and target_bed:
        assay = detect_assay(Path(target_bed))
    
    if assay:
        bundled = get_bundled_gene_bed(assay, genome)
        if bundled:
            logger.info(f"Loading bundled gene BED for assay {assay}: {bundled}")
            return parse_gene_bed(bundled)
    
    # Option 3: Parse target BED on-the-fly
    if target_bed and Path(target_bed).exists():
        logger.info(f"Parsing target BED on-the-fly: {target_bed}")
        return parse_gene_bed(target_bed)
    
    raise ValueError(
        "No valid gene source: provide gene_bed, target_bed, or known assay"
    )
