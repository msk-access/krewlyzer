"""
Unified sample processing for krewlyzer.

Provides a single entry point for processing BAM/CRAM/BED.gz samples.
Used by:
  - wrapper.py (run-all): Single sample processing
  - pon/build.py (build-pon): Multi-sample PON building

This consolidates extraction and feature processing logic that was
previously duplicated across wrapper.py and pon/build.py.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, List
import logging
import tempfile

import pandas as pd

from krewlyzer import _core
from .unified_processor import run_features, FeatureOutputs
from ..assets import AssetManager

logger = logging.getLogger("krewlyzer.core.sample_processor")


@dataclass
class SampleParams:
    """Parameters for sample processing."""
    genome: str = "hg19"
    
    # Extraction parameters (for BAM/CRAM)
    mapq: int = 20
    minlen: int = 65
    maxlen: int = 400
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
    """All outputs from processing a single sample.
    
    Contains both file paths (from feature extraction) and in-memory
    data (for PON building).
    """
    sample_id: str
    bed_path: Path                          # Path to BED.gz
    fragment_count: int = 0                 # Total fragments
    
    # Feature outputs (paths from run_features)
    feature_outputs: Optional[FeatureOutputs] = None
    
    # Motif data (in-memory for PON building)
    mds_counts: Dict[str, float] = field(default_factory=dict)
    mds_score: float = 0.0
    em_counts: List[int] = field(default_factory=list)
    bpm_counts: List[int] = field(default_factory=list)
    gc_observations: List = field(default_factory=list)
    
    # On-target motif data (panel mode)
    em_counts_ontarget: List[int] = field(default_factory=list)
    bpm_counts_ontarget: List[int] = field(default_factory=list)
    gc_observations_ontarget: List = field(default_factory=list)


def _generate_kmers(k: int) -> List[str]:
    """Generate all possible k-mers of length k."""
    import itertools
    bases = ['A', 'C', 'G', 'T']
    return [''.join(p) for p in itertools.product(bases, repeat=k)]


def _compute_kmer_frequencies(em_counts: List[int], kmer: int) -> Dict[str, float]:
    """Convert raw k-mer counts to frequency dict."""
    kmers = _generate_kmers(kmer)
    total = sum(em_counts) if em_counts else 1
    if total == 0:
        total = 1
    return {k: c / total for k, c in zip(kmers, em_counts)}


def _compute_mds_from_counts(em_counts: List[int]) -> float:
    """Compute Motif Diversity Score from k-mer counts.
    
    MDS is the Shannon entropy of the k-mer distribution.
    """
    import math
    total = sum(em_counts) if em_counts else 0
    if total == 0:
        return 0.0
    
    entropy = 0.0
    for count in em_counts:
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)
    
    return entropy


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
    # PON mode (skip normalization, return raw data)
    pon_mode: bool = False,
    pon_model: Optional[Path] = None,
    assay: Optional[str] = None,
    temp_dir: Optional[Path] = None,
) -> SampleOutputs:
    """
    Process a single sample (BAM/CRAM/BED.gz) through feature extraction.
    
    This is the unified entry point for:
      - wrapper.py (run-all): Single sample with all features
      - pon/build.py: Multiple samples for PON building
    
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
        temp_dir: Temporary directory
        
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
        
        bed_path = output_dir / f"{sample_name}.bed.gz"
        
        # Resolve exclude regions
        exclude_path = params.exclude_regions or assets.exclude_regions
        
        logger.info(f"Extracting fragments from BAM...")
        
        # Extract + motif in one pass
        result = _core.extract_motif.process_bam_parallel(
            str(input_path),
            str(reference),
            params.mapq,
            params.minlen,
            params.maxlen,
            params.kmer,
            params.threads,
            str(bed_path),
            "enable",  # motif output
            str(exclude_path) if exclude_path and exclude_path.exists() else None,
            str(target_regions) if target_regions and target_regions.exists() else None,
            params.skip_duplicates,
            params.require_proper_pair,
            True  # silent
        )
        
        # Unpack results
        (fragment_count, em_counts, bpm_counts, gc_obs,
         em_counts_on, bpm_counts_on, gc_obs_on) = result
        
        # Store extraction results
        outputs.bed_path = bed_path
        outputs.fragment_count = fragment_count
        outputs.em_counts = list(em_counts) if em_counts else []
        outputs.bpm_counts = list(bpm_counts) if bpm_counts else []
        outputs.gc_observations = list(gc_obs) if gc_obs else []
        outputs.em_counts_ontarget = list(em_counts_on) if em_counts_on else []
        outputs.bpm_counts_ontarget = list(bpm_counts_on) if bpm_counts_on else []
        outputs.gc_observations_ontarget = list(gc_obs_on) if gc_obs_on else []
        
        # Compute MDS
        if outputs.em_counts:
            outputs.mds_counts = _compute_kmer_frequencies(outputs.em_counts, params.kmer)
            outputs.mds_score = _compute_mds_from_counts(outputs.em_counts)
        
        logger.info(f"Extracted {fragment_count:,} fragments")
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
        target_regions=target_regions,
        assay=assay,
        pon_model=None if pon_mode else pon_model,  # Skip PON in PON-build mode
        fsc_bins=params.bin_file,
        fsd_arms=params.arms_file,
        wps_anchors=params.wps_anchors,
        wps_background=params.wps_background,
        ocf_regions=params.ocf_regions,
        threads=params.threads,
    )
    
    outputs.feature_outputs = feature_outputs
    
    logger.info(f"Sample processing complete: {sample_name}")
    return outputs
