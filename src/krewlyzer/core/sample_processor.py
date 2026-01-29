"""
Unified sample extraction for krewlyzer.

This module provides the core extraction infrastructure used by all commands:
  - krewlyzer extract: Write BED.gz + GC factors
  - krewlyzer motif: Write motif files
  - krewlyzer run-all: Full pipeline
  - krewlyzer build-pon: Multi-sample PON building

Architecture:
    extract_sample()          - Core extraction, returns ExtractionResult
    write_motif_outputs()     - Write motif files from ExtractionResult
    write_extraction_outputs() - Write BED metadata + GC factors
    process_sample()          - High-level orchestrator (backward compatible)

The separation allows each command to use exactly what it needs while
sharing the expensive extraction logic.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, List, Any
import logging
import time
import shutil
import json
from datetime import datetime

import pandas as pd

from krewlyzer import _core
from .unified_processor import run_features, FeatureOutputs
from .motif_processor import compute_mds
from ..assets import AssetManager

logger = logging.getLogger("krewlyzer.core.sample_processor")


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class ExtractionResult:
    """
    Results from BAM/CRAM extraction.
    
    Contains all data extracted in a single pass through the input file.
    This is the contract between extraction and downstream writers/processors.
    
    Attributes:
        sample_id: Sample identifier
        input_path: Original input file path
        fragment_count: Total fragments passing filters
        bed_path: Path to BED.gz (if write_bed was True)
        
        # Off-target data (unbiased - use for biomarkers)
        em_counts: End motif k-mer counts (256 values for k=4)
        bpm_counts: Breakpoint motif k-mer counts
        gc_observations: GC content observations per bin
        
        # On-target data (PCR-affected - for comparison only in panel mode)
        em_counts_ontarget: End motif counts from on-target fragments
        bpm_counts_ontarget: Breakpoint motif counts from on-target
        gc_observations_ontarget: On-target GC observations
        
        # Derived data (computed from em_counts)
        mds_score: Motif Diversity Score (Shannon entropy)
        kmer_frequencies: Dict of k-mer -> frequency
        
        # Metadata
        is_panel_mode: Whether target_regions was provided
        extraction_time_seconds: Performance metric
        kmer: K-mer size used for motif extraction
    """
    sample_id: str
    input_path: Path
    fragment_count: int = 0
    bed_path: Optional[Path] = None
    
    # Off-target data (primary, unbiased) - dicts map k-mer -> count or (len_bin, gc) -> count
    em_counts: Dict[str, int] = field(default_factory=dict)
    bpm_counts: Dict[str, int] = field(default_factory=dict)
    gc_observations: Dict[Any, Any] = field(default_factory=dict)
    
    # On-target data (panel mode only) - dicts map k-mer -> count or (len_bin, gc) -> count
    em_counts_ontarget: Dict[str, int] = field(default_factory=dict)
    bpm_counts_ontarget: Dict[str, int] = field(default_factory=dict)
    gc_observations_ontarget: Dict[Any, Any] = field(default_factory=dict)
    
    # Derived data
    mds_score: float = 0.0
    kmer_frequencies: Dict[str, float] = field(default_factory=dict)
    
    # Metadata
    is_panel_mode: bool = False
    extraction_time_seconds: float = 0.0
    kmer: int = 4
    target_regions_path: Optional[str] = None  # For metadata output


@dataclass
class SampleParams:
    """
    Parameters for sample processing.
    
    Encapsulates all extraction parameters to simplify function signatures.
    
    Attributes:
        genome: Genome build (hg19/GRCh37/hg38/GRCh38)
        mapq: Minimum mapping quality threshold
        minlen: Minimum fragment length (bp)
        maxlen: Maximum fragment length (bp)
        skip_duplicates: Skip PCR duplicate reads
        require_proper_pair: Require proper read pairs
        kmer: K-mer size for motif extraction
        threads: Number of threads (0=auto)
        
        # Asset overrides (None = use bundled defaults)
        bin_file: Custom bin file for FSC/FSR
        arms_file: Custom chromosome arms file for FSD
        wps_anchors: Custom WPS anchor regions
        wps_background: Custom WPS background (Alu) regions
        ocf_regions: Custom OCF regions
        exclude_regions: Custom exclusion regions
    """
    genome: str = "hg19"
    
    # Extraction parameters
    mapq: int = 20
    minlen: int = 65
    maxlen: int = 1000
    skip_duplicates: bool = True
    require_proper_pair: bool = True
    kmer: int = 4
    threads: int = 4
    
    # Asset overrides (None = use bundled)
    bin_file: Optional[Path] = None
    arms_file: Optional[Path] = None
    wps_anchors: Optional[Path] = None
    wps_background: Optional[Path] = None
    ocf_regions: Optional[Path] = None
    exclude_regions: Optional[Path] = None


@dataclass
class SampleOutputs:
    """
    All outputs from processing a single sample.
    
    Contains both file paths (from feature extraction) and in-memory
    data (for PON building). This is the return type of process_sample().
    
    Note: For extraction-only use cases, use ExtractionResult directly.
    """
    sample_id: str
    bed_path: Path
    fragment_count: int = 0
    
    # Feature outputs (paths from run_features)
    feature_outputs: Optional[FeatureOutputs] = None
    
    # Motif data (in-memory for PON building)
    mds_counts: Dict[str, float] = field(default_factory=dict)
    mds_score: float = 0.0
    em_counts: List[int] = field(default_factory=list)
    bpm_counts: List[int] = field(default_factory=list)
    gc_observations: List[Any] = field(default_factory=list)
    
    # On-target motif data (panel mode)
    em_counts_ontarget: List[int] = field(default_factory=list)
    bpm_counts_ontarget: List[int] = field(default_factory=list)
    gc_observations_ontarget: List[Any] = field(default_factory=list)
    
    # Region Entropy data (in-memory for PON building)
    # Dict[label, (count, mean_size, entropy)]
    tfbs_data: Optional[Dict[str, Any]] = None  # Genome-wide (.TFBS.tsv)
    tfbs_data_ontarget: Optional[Dict[str, Any]] = None  # Panel-specific regions (.TFBS.ontarget.tsv)
    atac_data: Optional[Dict[str, Any]] = None  # Genome-wide (.ATAC.tsv)
    atac_data_ontarget: Optional[Dict[str, Any]] = None  # Panel-specific regions (.ATAC.ontarget.tsv)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Note: MDS and k-mer frequency computation now use compute_mds() from motif_processor


# =============================================================================
# CORE EXTRACTION FUNCTION
# =============================================================================

def extract_sample(
    input_path: Path,
    reference: Path,
    *,
    # Core params
    mapq: int = 20,
    minlen: int = 65,
    maxlen: int = 1000,
    kmer: int = 4,
    threads: int = 4,
    skip_duplicates: bool = True,
    require_proper_pair: bool = True,
    
    # Panel mode
    target_regions: Optional[Path] = None,
    exclude_regions: Optional[Path] = None,
    
    # Output control
    write_bed: bool = True,
    output_dir: Optional[Path] = None,
    sample_name: Optional[str] = None,
    
    # Genome for assets
    genome: str = "hg19",
) -> ExtractionResult:
    """
    Extract cfDNA fragments and motifs from BAM/CRAM in a single pass.
    
    This is the unified extraction entry point used by:
    - krewlyzer extract: Fragment extraction with GC correction
    - krewlyzer motif: Motif feature extraction
    - krewlyzer run-all: Full sample processing pipeline
    - krewlyzer build-pon: Multi-sample PON model building
    
    The function performs a single pass through the BAM file, extracting:
    - Fragment coordinates (optionally written to BED.gz)
    - End motif k-mer counts
    - Breakpoint motif k-mer counts
    - GC content observations for bias correction
    
    For panel data (when target_regions is provided), the function
    separates on-target and off-target data. Off-target GC observations
    are used for unbiased GC model training.
    
    Args:
        input_path: BAM or CRAM file (sorted, indexed)
        reference: Reference FASTA file (indexed with .fai)
        mapq: Minimum mapping quality threshold (default: 20)
        minlen: Minimum fragment length in bp (default: 65)
        maxlen: Maximum fragment length in bp (default: 400)
        kmer: K-mer size for motif extraction (default: 4)
        threads: Number of threads, 0=auto (default: 4)
        skip_duplicates: Skip PCR duplicate reads (default: True)
        require_proper_pair: Require proper read pairs (default: True)
        target_regions: BED file for panel mode (off-target GC only)
        exclude_regions: BED file for excluded regions (blacklist)
        write_bed: Whether to write BED.gz output (default: True)
        output_dir: Directory for BED.gz (required if write_bed=True)
        sample_name: Sample ID (derived from input filename if None)
        genome: Genome build for asset resolution (default: hg19)
        
    Returns:
        ExtractionResult containing all extraction data
        
    Raises:
        ValueError: If output_dir required but not provided
        FileNotFoundError: If input file doesn't exist
        RuntimeError: If Rust extraction fails
        
    Example:
        >>> result = extract_sample(
        ...     Path("sample.bam"),
        ...     Path("hg19.fa"),
        ...     target_regions=Path("targets.bed"),
        ...     write_bed=True,
        ...     output_dir=Path("./output")
        ... )
        >>> print(f"Extracted {result.fragment_count:,} fragments")
        >>> print(f"MDS: {result.mds_score:.4f}")
    """
    start_time = time.time()
    
    # Validate inputs
    input_path = Path(input_path)
    reference = Path(reference)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    if not reference.exists():
        raise FileNotFoundError(f"Reference file not found: {reference}")
    
    if write_bed and output_dir is None:
        raise ValueError("output_dir required when write_bed=True")
    
    # Derive sample name
    if sample_name is None:
        sample_name = input_path.stem.replace('.bam', '').replace('.cram', '')
    
    # Determine panel mode
    is_panel_mode = target_regions is not None and target_regions.exists()
    
    # Log extraction start
    logger.info(f"Extracting: {input_path.name}")
    logger.debug(f"  Sample: {sample_name}")
    logger.debug(f"  Params: mapq>={mapq}, len=[{minlen},{maxlen}], k={kmer}")
    logger.debug(f"  Filters: skip_dup={skip_duplicates}, proper_pair={require_proper_pair}")
    
    if is_panel_mode:
        logger.info(f"  Panel mode: off-target GC from {target_regions.name}")
    
    # Note: exclude_regions is NOT auto-resolved here.
    # The caller (wrapper.py, extract.py) should pass it explicitly.
    # This avoids issues with mismatched chromosome naming in tests.
    if exclude_regions:
        logger.debug(f"  Exclude regions: {exclude_regions.name}")
    
    # Prepare BED output path
    bed_path = None
    bed_output_arg = None
    if write_bed:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        bed_path = output_dir / f"{sample_name}.bed.gz"
        bed_output_arg = str(bed_path)
    
    # Call Rust extraction engine
    logger.debug("  Calling Rust extraction engine...")
    try:
        result = _core.extract_motif.process_bam_parallel(
            str(input_path),
            str(reference),
            mapq,
            minlen,
            maxlen,
            kmer,
            threads,
            bed_output_arg,
            "enable",  # Enable motif counting
            str(exclude_regions) if exclude_regions and exclude_regions.exists() else None,
            str(target_regions) if is_panel_mode else None,
            skip_duplicates,
            require_proper_pair,
            True  # silent=True (no Rust progress bar)
        )
    except Exception as e:
        logger.error(f"Rust extraction failed: {e}")
        raise RuntimeError(f"Extraction failed: {e}") from e
    
    # Unpack Rust results
    (fragment_count, em_counts_raw, bpm_counts_raw, gc_obs,
     em_counts_on_raw, bpm_counts_on_raw, gc_obs_on) = result
    
    # em_counts and bpm_counts from Rust are dicts (k-mer -> count)
    # Keep as dicts for compatibility with motif_processor
    em_counts = dict(em_counts_raw) if em_counts_raw else {}
    bpm_counts = dict(bpm_counts_raw) if bpm_counts_raw else {}
    # gc_observations is a dict from Rust: (length_bin, gc_pct) -> count
    # Keep as-is for Rust GC factor computation compatibility
    gc_observations = gc_obs if gc_obs else {}
    
    # On-target counts (also keep as dicts)
    em_counts_ontarget = dict(em_counts_on_raw) if em_counts_on_raw else {}
    bpm_counts_ontarget = dict(bpm_counts_on_raw) if bpm_counts_on_raw else {}
    gc_observations_ontarget = gc_obs_on if gc_obs_on else {}
    
    # Compute derived data
    kmer_frequencies = {}
    mds_score = 0.0
    if em_counts:
        # kmer_frequencies is just normalized em_counts
        total_em = sum(em_counts.values())
        if total_em > 0:
            kmer_frequencies = {k: v / total_em for k, v in em_counts.items()}
        mds_score = compute_mds(em_counts)
    
    elapsed = time.time() - start_time
    
    # Log results
    logger.info(f"  Extracted {fragment_count:,} fragments ({elapsed:.1f}s)")
    logger.debug(f"  Off-target: {sum(em_counts.values()):,} EM, {len(gc_observations):,} GC obs")
    if is_panel_mode:
        logger.debug(f"  On-target: {sum(em_counts_ontarget.values()):,} EM")
    logger.debug(f"  MDS: {mds_score:.4f}")
    
    return ExtractionResult(
        sample_id=sample_name,
        input_path=input_path,
        fragment_count=fragment_count,
        bed_path=bed_path,
        em_counts=em_counts,
        bpm_counts=bpm_counts,
        gc_observations=gc_observations,
        em_counts_ontarget=em_counts_ontarget,
        bpm_counts_ontarget=bpm_counts_ontarget,
        gc_observations_ontarget=gc_observations_ontarget,
        mds_score=mds_score,
        kmer_frequencies=kmer_frequencies,
        is_panel_mode=is_panel_mode,
        extraction_time_seconds=elapsed,
        kmer=kmer,
        target_regions_path=str(target_regions) if target_regions else None,
    )


# =============================================================================
# OUTPUT WRITING FUNCTIONS
# =============================================================================

def write_motif_outputs(
    result: ExtractionResult,
    output_dir: Path,
    *,
    pon: Optional[Any] = None,
    include_ontarget: bool = True,
) -> Dict[str, Path]:
    """
    Write motif analysis files from extraction result.
    
    Writes End Motif, Breakpoint Motif, and Motif Diversity Score files
    from the extraction data. If a PON model is provided, computes and
    includes MDS z-score for comparison against healthy baseline.
    
    Files written:
    - {sample}.EndMotif.tsv: End motif k-mer frequencies
    - {sample}.BreakPointMotif.tsv: Breakpoint motif frequencies
    - {sample}.MDS.tsv: Motif Diversity Score (with z-score if PON)
    
    For panel mode (when include_ontarget=True and on-target data exists):
    - {sample}.EndMotif.ontarget.tsv
    - {sample}.BreakPointMotif.ontarget.tsv
    - {sample}.MDS.ontarget.tsv
    
    Args:
        result: ExtractionResult from extract_sample()
        output_dir: Output directory for files
        pon: Optional PonModel for MDS z-score computation
        include_ontarget: Write on-target files for panel mode
        
    Returns:
        Dict mapping output type to file path:
        {'edm': Path, 'bpm': Path, 'mds': Path, ...}
        
    Example:
        >>> result = extract_sample(bam, ref)
        >>> outputs = write_motif_outputs(result, Path("./output"))
        >>> print(f"MDS file: {outputs['mds']}")
    """
    from .motif_processor import process_motif_outputs
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    sample = result.sample_id
    outputs = {}
    
    # Output file paths
    edm_output = output_dir / f"{sample}.EndMotif.tsv"
    bpm_output = output_dir / f"{sample}.BreakPointMotif.tsv"
    mds_output = output_dir / f"{sample}.MDS.tsv"
    
    logger.info(f"Writing motif outputs for {sample}")
    
    # Write off-target (primary) motif files
    total_em, total_bpm, mds = process_motif_outputs(
        em_counts=result.em_counts,
        bpm_counts=result.bpm_counts,
        edm_output=edm_output,
        bpm_output=bpm_output,
        mds_output=mds_output,
        sample_name=sample,
        kmer=result.kmer,
        include_headers=True
    )
    
    outputs['edm'] = edm_output
    outputs['bpm'] = bpm_output
    outputs['mds'] = mds_output
    
    logger.debug(f"  Off-target: {total_em:,} EM, {total_bpm:,} BPM, MDS={mds:.4f}")
    
    # Compute MDS z-score if PON provided
    if pon and hasattr(pon, 'mds_baseline') and pon.mds_baseline and mds is not None:
        mds_z = (mds - pon.mds_baseline.mds_mean) / max(pon.mds_baseline.mds_std, 1e-10)
        logger.debug(f"  MDS z-score: {mds_z:.3f}")
        
        # Append z-score to MDS file
        try:
            mds_df = pd.read_csv(mds_output, sep="\t")
            if "mds_z" not in mds_df.columns:
                mds_df["mds_z"] = mds_z
                mds_df.to_csv(mds_output, sep="\t", index=False)
        except Exception as e:
            logger.debug(f"  Could not add MDS z-score: {e}")
    
    # Write on-target motif files (panel mode)
    if include_ontarget and result.is_panel_mode and sum(result.em_counts_ontarget.values()) > 0:
        edm_on = output_dir / f"{sample}.EndMotif.ontarget.tsv"
        bpm_on = output_dir / f"{sample}.BreakPointMotif.ontarget.tsv"
        mds_on = output_dir / f"{sample}.MDS.ontarget.tsv"
        
        total_em_on, total_bpm_on, mds_on_val = process_motif_outputs(
            em_counts=result.em_counts_ontarget,
            bpm_counts=result.bpm_counts_ontarget,
            edm_output=edm_on,
            bpm_output=bpm_on,
            mds_output=mds_on,
            sample_name=sample,
            kmer=result.kmer,
            include_headers=True
        )
        
        outputs['edm_ontarget'] = edm_on
        outputs['bpm_ontarget'] = bpm_on
        outputs['mds_ontarget'] = mds_on
        
        logger.debug(f"  On-target: {total_em_on:,} EM, {total_bpm_on:,} BPM, MDS={mds_on_val:.4f}")
    
    # Write 1-mer End Motif (Jagged Index / C-end fraction)
    from .motif_processor import write_end_motif_1mer
    edm_1mer_output = output_dir / f"{sample}.EndMotif1mer.tsv"
    c_end_metrics = write_end_motif_1mer(
        em_counts=result.em_counts,
        output_path=edm_1mer_output,
        sample_name=sample,
        include_header=True
    )
    outputs['edm_1mer'] = edm_1mer_output
    logger.debug(f"  C-end fraction: {c_end_metrics['c_fraction']:.4f}")
    
    logger.info(f"  Wrote {len(outputs)} motif files")
    return outputs


def write_extraction_outputs(
    result: ExtractionResult,
    output_dir: Path,
    *,
    genome: str = "hg19",
    assay: Optional[str] = None,
    compute_gc_factors: bool = True,
) -> Dict[str, Path]:
    """
    Write extraction outputs (BED index, metadata, GC factors).
    
    After extract_sample() writes the BED.gz file, this function:
    - Creates tabix index for the BED.gz
    - Writes sample metadata JSON
    - Computes and writes GC correction factors
    
    Files written:
    - {sample}.bed.gz.tbi: Tabix index (for random access)
    - {sample}.metadata.json: Sample metadata and statistics
    - {sample}.correction_factors.tsv: GC correction factors
    
    Args:
        result: ExtractionResult from extract_sample()
        output_dir: Output directory
        genome: Genome build for asset resolution
        assay: Assay code for metadata
        compute_gc_factors: Whether to compute GC correction (default: True)
        
    Returns:
        Dict mapping output type to file path:
        {'tbi': Path, 'metadata': Path, 'gc_factors': Path}
        
    Raises:
        ValueError: If result.bed_path is None
    """
    if result.bed_path is None:
        raise ValueError("ExtractionResult has no bed_path - was write_bed=False?")
    
    output_dir = Path(output_dir)
    sample = result.sample_id
    outputs = {}
    
    logger.info(f"Writing extraction outputs for {sample}")
    
    # Index BED.gz with tabix
    if result.bed_path.exists():
        try:
            import pysam
            pysam.tabix_index(str(result.bed_path), preset="bed", force=True)
            tbi_path = result.bed_path.with_suffix('.gz.tbi')
            outputs['tbi'] = tbi_path
            logger.debug(f"  Created tabix index: {tbi_path.name}")
        except Exception as e:
            logger.warning(f"  Tabix indexing failed: {e}")
    
    # Write metadata JSON
    meta_path = output_dir / f"{sample}.metadata.json"
    metadata = {
        "sample_id": sample,
        "total_fragments": result.fragment_count,
        "genome": genome,
        "assay": assay,
        "panel_mode": result.is_panel_mode,
        "target_regions": result.target_regions_path,
        "extraction_time_seconds": round(result.extraction_time_seconds, 2),
        "mds_score": round(result.mds_score, 6),
        "filters": {
            "mapq": 20,
            "min_length": 65,
            "max_length": 1000,
        },
        "timestamp": datetime.now().isoformat(),
    }
    
    with open(meta_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    outputs['metadata'] = meta_path
    logger.debug(f"  Wrote metadata: {meta_path.name}")
    
    # Compute GC correction factors
    if compute_gc_factors and result.gc_observations:
        factors_path = output_dir / f"{sample}.correction_factors.tsv"
        try:
            assets = AssetManager(genome)
            gc_ref = assets.resolve("gc_reference")
            valid_regions = assets.resolve("valid_regions")
            
            n_factors = _core.gc.compute_and_write_gc_factors(
                result.gc_observations,
                str(gc_ref),
                str(valid_regions),
                str(factors_path)
            )
            
            outputs['gc_factors'] = factors_path
            # Update metadata with GC info
            metadata["gc_correction_computed"] = True
            with open(meta_path, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            logger.debug(f"  Wrote GC factors: {factors_path.name} ({n_factors} regions)")
        except Exception as e:
            logger.debug(f"  GC factor computation failed: {e}")
            metadata["gc_correction_computed"] = False
    
    logger.info(f"  Wrote {len(outputs)} extraction outputs")
    return outputs


# =============================================================================
# HIGH-LEVEL ORCHESTRATOR (BACKWARD COMPATIBLE)
# =============================================================================

def process_sample(
    input_path: Path,
    output_dir: Path,
    sample_name: str,
    *,
    reference: Optional[Path] = None,
    params: Optional[SampleParams] = None,
    target_regions: Optional[Path] = None,
    # Feature toggles
    enable_fsc: bool = True,
    enable_fsr: bool = True,
    enable_fsd: bool = True,
    enable_wps: bool = True,
    enable_ocf: bool = True,
    enable_tfbs: bool = True,
    enable_atac: bool = True,
    # PON mode (skip normalization, return raw data)
    pon_mode: bool = False,
    pon_model: Optional[Path] = None,
    assay: Optional[str] = None,
    temp_dir: Optional[Path] = None,
) -> SampleOutputs:
    """
    Process a single sample through full feature extraction pipeline.
    
    This is the high-level orchestrator that combines extraction with
    all feature extraction. It is used by:
    - wrapper.py (run-all): Single sample with all features
    - pon/build.py (build-pon): Multi-sample PON building
    
    For extraction-only use cases, use extract_sample() directly.
    
    Args:
        input_path: BAM, CRAM, or BED.gz file
        output_dir: Output directory for feature files
        sample_name: Sample identifier
        reference: Reference FASTA (required for BAM/CRAM)
        params: Extraction parameters (defaults used if None)
        target_regions: Panel target BED for on/off-target split
        enable_*: Feature toggles
        pon_mode: Skip PON normalization (for PON building)
        pon_model: PON model for normalization (ignored if pon_mode=True)
        assay: Assay code for panel-specific assets
        temp_dir: Temporary directory (unused, for compatibility)
        
    Returns:
        SampleOutputs with all results
        
    Raises:
        ValueError: If reference is required but not provided
        FileNotFoundError: If input file doesn't exist
    """
    params = params or SampleParams()
    assets = AssetManager(params.genome)
    
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    # Detect input type
    is_bam = str(input_path).lower().endswith(('.bam', '.cram'))
    
    # Initialize outputs
    outputs = SampleOutputs(sample_id=sample_name, bed_path=input_path)
    
    logger.info(f"Processing sample: {sample_name}")
    logger.debug(f"Input: {input_path} ({'BAM/CRAM' if is_bam else 'BED.gz'})")
    
    # === Phase 1: Extract from BAM if needed ===
    if is_bam:
        if reference is None:
            raise ValueError("Reference FASTA required for BAM/CRAM input")
        
        # Use new extract_sample function
        result = extract_sample(
            input_path,
            reference,
            mapq=params.mapq,
            minlen=params.minlen,
            maxlen=params.maxlen,
            kmer=params.kmer,
            threads=params.threads,
            skip_duplicates=params.skip_duplicates,
            require_proper_pair=params.require_proper_pair,
            target_regions=target_regions,
            exclude_regions=params.exclude_regions,
            write_bed=True,
            output_dir=output_dir,
            sample_name=sample_name,
            genome=params.genome,
        )
        
        # Transfer to SampleOutputs
        outputs.bed_path = result.bed_path
        outputs.fragment_count = result.fragment_count
        outputs.em_counts = result.em_counts
        outputs.bpm_counts = result.bpm_counts
        outputs.gc_observations = result.gc_observations
        outputs.em_counts_ontarget = result.em_counts_ontarget
        outputs.bpm_counts_ontarget = result.bpm_counts_ontarget
        outputs.gc_observations_ontarget = result.gc_observations_ontarget
        outputs.mds_counts = result.kmer_frequencies
        outputs.mds_score = result.mds_score
    else:
        outputs.bed_path = input_path
        logger.info(f"Using pre-extracted BED.gz")
    
    # === Phase 2: Run feature extraction ===
    logger.info(f"Running feature extraction...")
    
    feature_outputs = run_features(
        bed_path=outputs.bed_path,
        output_dir=output_dir,
        sample_name=sample_name,
        genome=params.genome,
        enable_fsc=enable_fsc,
        enable_fsr=enable_fsr,
        enable_fsd=enable_fsd,
        enable_wps=enable_wps,
        enable_ocf=enable_ocf,
        enable_tfbs=enable_tfbs,
        enable_atac=enable_atac,
        target_regions=target_regions,
        assay=assay,
        pon_model=None if pon_mode else pon_model,
        fsc_bins=params.bin_file,
        fsd_arms=params.arms_file,
        wps_anchors=params.wps_anchors,
        wps_background=params.wps_background,
        ocf_regions=params.ocf_regions,
        threads=params.threads,
    )
    
    outputs.feature_outputs = feature_outputs
    
    # === Phase 3: Extract TFBS/ATAC entropy data for PON building ===
    if pon_mode:
        from .region_entropy_processor import load_entropy_tsv, extract_entropy_data
        
        # TFBS entropy data
        if feature_outputs.tfbs and feature_outputs.tfbs.exists():
            try:
                df = load_entropy_tsv(feature_outputs.tfbs)
                outputs.tfbs_data = extract_entropy_data(df)
                logger.debug(f"  Extracted {len(outputs.tfbs_data)} TFBS labels")
            except Exception as e:
                logger.warning(f"Failed to extract TFBS data: {e}")
        
        if feature_outputs.tfbs_ontarget and feature_outputs.tfbs_ontarget.exists():
            try:
                df = load_entropy_tsv(feature_outputs.tfbs_ontarget)
                outputs.tfbs_data_ontarget = extract_entropy_data(df)
            except Exception:
                pass
        
        
        # ATAC entropy data
        if feature_outputs.atac and feature_outputs.atac.exists():
            try:
                df = load_entropy_tsv(feature_outputs.atac)
                outputs.atac_data = extract_entropy_data(df)
                logger.debug(f"  Extracted {len(outputs.atac_data)} ATAC labels")
            except Exception as e:
                logger.warning(f"Failed to extract ATAC data: {e}")
        
        if feature_outputs.atac_ontarget and feature_outputs.atac_ontarget.exists():
            try:
                df = load_entropy_tsv(feature_outputs.atac_ontarget)
                outputs.atac_data_ontarget = extract_entropy_data(df)
                logger.debug(f"  Extracted {len(outputs.atac_data_ontarget)} ATAC on-target labels")
            except Exception:
                pass
    
    logger.info(f"Sample processing complete: {sample_name}")
    return outputs
