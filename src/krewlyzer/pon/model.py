"""
PON Model definitions for Krewlyzer.

This module defines the unified PON model format using Parquet storage.
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
import logging

import numpy as np
import pandas as pd

logger = logging.getLogger("pon")


@dataclass
class GcBiasModel:
    """
    GC bias correction curves for fragment types.
    
    Stores expected relative coverage as a function of GC content,
    learned from healthy plasma samples.
    """
    gc_bins: List[float]  # e.g., [0.25, 0.30, 0.35, ...]
    short_expected: List[float]  # Expected relative coverage for short fragments
    short_std: List[float]
    intermediate_expected: List[float]
    intermediate_std: List[float]
    long_expected: List[float]
    long_std: List[float]
    # WPS-specific (uses different size ranges)
    wps_long_expected: Optional[List[float]] = None  # 120-180bp
    wps_long_std: Optional[List[float]] = None
    wps_short_expected: Optional[List[float]] = None  # 35-80bp
    wps_short_std: Optional[List[float]] = None
    
    def get_expected(self, gc: float, frag_type: str) -> float:
        """
        Interpolate expected coverage for given GC content and fragment type.
        
        Args:
            gc: GC content (0.0-1.0)
            frag_type: One of "short", "intermediate", "long", "wps_long", "wps_short"
            
        Returns:
            Expected relative coverage (1.0 = no bias)
        """
        expected_map = {
            "short": self.short_expected,
            "intermediate": self.intermediate_expected,
            "long": self.long_expected,
            "wps_long": self.wps_long_expected,
            "wps_short": self.wps_short_expected,
        }
        expected = expected_map.get(frag_type)
        if expected is None:
            return 1.0
            
        # Linear interpolation
        return float(np.interp(gc, self.gc_bins, expected))


@dataclass
class FsdBaseline:
    """
    Fragment Size Distribution baseline per chromosome arm.
    
    Stores expected size proportions for healthy samples.
    """
    size_bins: List[int]  # e.g., [65, 70, 75, 80, ...]
    arms: Dict[str, Dict[str, List[float]]]  # arm -> {"expected": [...], "std": [...]}
    
    def get_expected(self, arm: str, size: int) -> float:
        """Get expected proportion for a size bin in a given arm."""
        if arm not in self.arms:
            return 0.0
        arm_data = self.arms[arm]
        return float(np.interp(size, self.size_bins, arm_data["expected"]))
    
    def get_std(self, arm: str, size: int) -> float:
        """Get standard deviation for a size bin in a given arm."""
        if arm not in self.arms:
            return 0.0
        arm_data = self.arms[arm]
        return float(np.interp(size, self.size_bins, arm_data["std"]))
    
    def get_stats(self, arm: str, size_bin: Optional[int] = None) -> Optional[tuple]:
        """
        Get (mean, std) for a chromosome arm.
        
        Used by pon_integration.compute_fsd_zscore().
        
        Args:
            arm: Chromosome arm (e.g., "1p", "12q")
            size_bin: Optional size bin - if None, returns aggregate stats
            
        Returns:
            Tuple (expected, std) or None if arm not found
        """
        if arm not in self.arms:
            return None
        arm_data = self.arms[arm]
        if size_bin is not None:
            return (self.get_expected(arm, size_bin), self.get_std(arm, size_bin))
        # Return median bin stats as aggregate
        mid_idx = len(self.size_bins) // 2
        return (arm_data["expected"][mid_idx], arm_data["std"][mid_idx])


@dataclass
class WpsBaseline:
    """
    WPS baseline per transcript region.
    
    Supports two formats:
    - Legacy (v1.0): Scalar mean/std per region (wps_long_mean, wps_short_mean)
    - Vector (v2.0): 200-bin mean/std vectors per region for ML integration
    
    Schema v2.0 enables position-specific z-score computation for cancer detection.
    """
    regions: pd.DataFrame  # Columns depend on schema version
    schema_version: str = "1.0"  # "1.0" = scalar, "2.0" = vector
    
    def get_baseline(self, region_id: str) -> Optional[Dict[str, float]]:
        """Get scalar baseline stats for a region (legacy v1.0 format)."""
        match = self.regions[self.regions["region_id"] == region_id]
        if match.empty:
            return None
        row = match.iloc[0]
        return {
            "wps_long_mean": row.get("wps_long_mean", 0),
            "wps_long_std": row.get("wps_long_std", 1),
            "wps_short_mean": row.get("wps_short_mean", 0),
            "wps_short_std": row.get("wps_short_std", 1),
        }
    
    def get_stats(self, region_id: str, wps_type: str = "wps_long") -> Optional[tuple]:
        """
        Get (mean, std) for a region.
        
        Used by pon_integration.compute_wps_zscore().
        
        Args:
            region_id: Region/gene identifier
            wps_type: Either "wps_long" or "wps_short"
            
        Returns:
            Tuple (mean, std) or None if region not found
        """
        baseline = self.get_baseline(region_id)
        if baseline is None:
            return None
        return (baseline[f"{wps_type}_mean"], baseline[f"{wps_type}_std"])
    
    # --- Vector format methods (v2.0) ---
    
    def get_baseline_vector(self, region_id: str, column: str = "wps_nuc") -> Optional[np.ndarray]:
        """
        Get 200-bin mean vector for a region (v2.0 format).
        
        Args:
            region_id: Region identifier
            column: Vector column prefix ('wps_nuc', 'wps_tf', 'prot_frac_nuc', 'prot_frac_tf')
            
        Returns:
            200-element numpy array or None if not found
        """
        if self.schema_version != "2.0":
            return None
        match = self.regions[self.regions["region_id"] == region_id]
        if match.empty:
            return None
        mean_col = f"{column}_mean"
        if mean_col not in match.columns:
            return None
        return np.array(match.iloc[0][mean_col])
    
    def get_std_vector(self, region_id: str, column: str = "wps_nuc") -> Optional[np.ndarray]:
        """Get 200-bin std vector for a region (v2.0 format)."""
        if self.schema_version != "2.0":
            return None
        match = self.regions[self.regions["region_id"] == region_id]
        if match.empty:
            return None
        std_col = f"{column}_std"
        if std_col not in match.columns:
            return None
        return np.array(match.iloc[0][std_col])
    
    def compute_z_vector(self, region_id: str, sample_vector: np.ndarray,
                          column: str = "wps_nuc") -> Optional[np.ndarray]:
        """
        Compute position-specific z-scores for sample vs PoN (v2.0 format).
        
        Args:
            region_id: Region identifier
            sample_vector: 200-element sample WPS vector
            column: Vector column prefix
            
        Returns:
            200-element z-score array or None if baseline not found
        """
        mean = self.get_baseline_vector(region_id, column)
        std = self.get_std_vector(region_id, column)
        if mean is None or std is None:
            return None
        
        # Avoid division by zero
        std_safe = np.where(std > 0, std, 1.0)
        return (sample_vector - mean) / std_safe


@dataclass
class OcfBaseline:
    """
    OCF (Open Chromatin Footprinting) baseline per region.
    
    Stores expected OCF and synchronization scores for open chromatin regions
    (promoters, enhancers) from healthy plasma samples.
    
    Used for z-score normalization of sample OCF profiles.
    """
    # DataFrame: region_id, ocf_mean, ocf_std, sync_mean, sync_std
    regions: pd.DataFrame
    
    def get_stats(self, region_id: str) -> Optional[tuple]:
        """
        Get (mean, std) for a region's OCF score.
        
        Args:
            region_id: Region identifier
            
        Returns:
            Tuple (ocf_mean, ocf_std) or None if not found
        """
        match = self.regions[self.regions["region_id"] == region_id]
        if match.empty:
            return None
        row = match.iloc[0]
        return (row["ocf_mean"], row["ocf_std"])
    
    def get_sync_stats(self, region_id: str) -> Optional[tuple]:
        """
        Get (mean, std) for a region's synchronization score.
        
        Args:
            region_id: Region identifier
            
        Returns:
            Tuple (sync_mean, sync_std) or None if not found
        """
        match = self.regions[self.regions["region_id"] == region_id]
        if match.empty:
            return None
        row = match.iloc[0]
        return (row.get("sync_mean", 0), row.get("sync_std", 1))
    
    def compute_zscore(self, region_id: str, observed_ocf: float) -> Optional[float]:
        """Compute z-score for observed OCF value."""
        stats = self.get_stats(region_id)
        if stats is None:
            return None
        mean, std = stats
        if std > 0:
            return (observed_ocf - mean) / std
        return 0.0


@dataclass
class WpsBackgroundBaseline:
    """
    WPS background (Alu) baseline for periodicity/NRL analysis.
    
    Stores nucleosome repeat length (NRL) statistics from Alu element stacking
    across healthy plasma samples. Used for:
    - NRL deviation scoring (cancer detection)
    - Periodicity z-score computation
    
    The Alu stacking method provides a robust measure of nucleosome spacing
    independent of gene expression patterns.
    """
    # DataFrame: group_id, nrl_mean, nrl_std, periodicity_mean, periodicity_std
    groups: pd.DataFrame
    
    def get_nrl_stats(self, group_id: str = "all") -> Optional[tuple]:
        """
        Get (mean, std) for nucleosome repeat length.
        
        Args:
            group_id: Group identifier (default: "all" for genome-wide)
            
        Returns:
            Tuple (nrl_mean, nrl_std) or None if not found
        """
        match = self.groups[self.groups["group_id"] == group_id]
        if match.empty:
            return None
        row = match.iloc[0]
        return (row["nrl_mean"], row["nrl_std"])
    
    def get_periodicity_stats(self, group_id: str = "all") -> Optional[tuple]:
        """
        Get (mean, std) for periodicity score.
        
        Args:
            group_id: Group identifier (default: "all" for genome-wide)
            
        Returns:
            Tuple (periodicity_mean, periodicity_std) or None if not found
        """
        match = self.groups[self.groups["group_id"] == group_id]
        if match.empty:
            return None
        row = match.iloc[0]
        return (row.get("periodicity_mean", 0), row.get("periodicity_std", 1))
    
    def compute_nrl_zscore(self, observed_nrl: float, group_id: str = "all") -> Optional[float]:
        """Compute z-score for observed NRL value."""
        stats = self.get_nrl_stats(group_id)
        if stats is None:
            return None
        mean, std = stats
        if std > 0:
            return (observed_nrl - mean) / std
        return 0.0


@dataclass
class MdsBaseline:
    """
    Motif Diversity Score baseline from healthy plasma samples.
    
    Stores expected k-mer frequencies and MDS value for healthy samples.
    Enables detection of aberrant end-motif patterns in cancer cfDNA.
    
    Note: DNA damage and tissue-of-origin affect motif patterns.
    """
    # Expected k-mer frequencies (256 4-mers)
    kmer_expected: Dict[str, float] = field(default_factory=dict)  # e.g., {"ACGT": 0.0042, ...}
    kmer_std: Dict[str, float] = field(default_factory=dict)
    
    # Global MDS baseline
    mds_mean: float = 0.0
    mds_std: float = 1.0
    
    def get_kmer_zscore(self, kmer: str, observed_freq: float) -> Optional[float]:
        """Compute z-score for a k-mer frequency."""
        if kmer not in self.kmer_expected:
            return None
        expected = self.kmer_expected[kmer]
        std = self.kmer_std.get(kmer, 0.001)
        if std > 0:
            return (observed_freq - expected) / std
        return 0.0
    
    def get_mds_zscore(self, observed_mds: float) -> float:
        """Compute z-score for MDS value."""
        if self.mds_std > 0:
            return (observed_mds - self.mds_mean) / self.mds_std
        return 0.0
    
    def get_aberrant_kmers(self, observed_freqs: Dict[str, float], 
                            threshold: float = 2.0) -> Dict[str, float]:
        """
        Find k-mers with z-scores exceeding threshold.
        
        Args:
            observed_freqs: Dict of kmer -> observed frequency
            threshold: Z-score threshold for aberrant detection
            
        Returns:
            Dict of aberrant kmers -> z-scores
        """
        aberrant = {}
        for kmer, freq in observed_freqs.items():
            zscore = self.get_kmer_zscore(kmer, freq)
            if zscore is not None and abs(zscore) > threshold:
                aberrant[kmer] = zscore
        return aberrant


@dataclass
class RegionMdsBaseline:
    """
    Per-Gene Region MDS baseline from healthy plasma samples.
    
    Stores expected MDS scores per gene/target for z-score normalization.
    Based on Helzer et al. (2025) methodology.
    
    Fields:
        gene_baseline: Dict mapping gene -> {mds_mean, mds_std, mds_e1_mean, mds_e1_std, n_samples}
    """
    # Per-gene baseline: gene -> {mds_mean, mds_std, mds_e1_mean, mds_e1_std, n_samples}
    gene_baseline: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    def get_stats(self, gene: str) -> Optional[tuple]:
        """
        Get (mean, std) for a gene's MDS score.
        
        Args:
            gene: Gene symbol
            
        Returns:
            Tuple (mds_mean, mds_std) or None if not found
        """
        if gene in self.gene_baseline:
            data = self.gene_baseline[gene]
            return (data.get("mds_mean", 0.0), data.get("mds_std", 1.0))
        return None
    
    def get_e1_stats(self, gene: str) -> Optional[tuple]:
        """
        Get (mean, std) for a gene's E1 (first exon) MDS score.
        
        Args:
            gene: Gene symbol
            
        Returns:
            Tuple (mds_e1_mean, mds_e1_std) or None if not found
        """
        if gene in self.gene_baseline:
            data = self.gene_baseline[gene]
            return (data.get("mds_e1_mean", 0.0), data.get("mds_e1_std", 1.0))
        return None
    
    def compute_zscore(self, gene: str, observed_mds: float) -> Optional[float]:
        """
        Compute z-score for observed MDS value.
        
        Args:
            gene: Gene symbol
            observed_mds: Observed MDS value
            
        Returns:
            Z-score or None if gene not in baseline
        """
        stats = self.get_stats(gene)
        if stats is None:
            return None
        mean, std = stats
        if std > 0:
            return (observed_mds - mean) / std
        return 0.0
    
    def compute_e1_zscore(self, gene: str, observed_mds_e1: float) -> Optional[float]:
        """
        Compute z-score for observed E1 MDS value.
        """
        stats = self.get_e1_stats(gene)
        if stats is None:
            return None
        mean, std = stats
        if std > 0:
            return (observed_mds_e1 - mean) / std
        return 0.0


@dataclass
class FscGeneBaseline:
    """
    Per-Gene Fragment Size Coverage (FSC) depth baseline.
    
    Stores normalized depth statistics (mean/std) across PON samples for each gene.
    Enables detection of copy number variations by comparing sample depth against
    the healthy baseline.
    
    This is panel-mode only - generated when --assay is specified during build-pon.
    
    Fields:
        data: Dict mapping gene -> (mean_depth, std_depth, n_samples)
        
    Clinical use:
        - Gene amplification: z-score >> 0
        - Gene deletion: z-score << 0
        - Normal copy number: z-score ≈ 0
        
    Note:
        Requires minimum 3 samples for reliable statistics.
    """
    # gene -> (mean_depth, std_depth, n_samples)
    data: Dict[str, tuple] = field(default_factory=dict)
    
    def get_stats(self, gene: str) -> Optional[tuple]:
        """
        Get (mean, std) for a gene's normalized depth.
        
        Args:
            gene: Gene symbol (e.g., "TP53", "EGFR")
            
        Returns:
            Tuple (mean, std) or None if gene not in baseline
        """
        if gene in self.data:
            return self.data[gene][:2]  # (mean, std), skip n_samples
        return None
    
    def compute_zscore(self, gene: str, observed_depth: float) -> Optional[float]:
        """
        Compute z-score for observed normalized depth.
        
        Args:
            gene: Gene symbol
            observed_depth: Sample's normalized depth for this gene
            
        Returns:
            Z-score or None if gene not in baseline
        """
        stats = self.get_stats(gene)
        if stats is None:
            return None
        mean, std = stats
        if std > 0:
            return (observed_depth - mean) / std
        return 0.0
    
    def __len__(self) -> int:
        """Return number of genes in baseline."""
        return len(self.data)


@dataclass
class FscRegionBaseline:
    """
    Per-Region (exon/probe) Fragment Size Coverage depth baseline.
    
    Stores normalized depth statistics (mean/std) for each exon/probe region.
    More granular than FscGeneBaseline - useful for detecting focal copy number
    changes affecting single exons.
    
    Region IDs are formatted as "chrom:start-end" (e.g., "chr17:7571720-7573008").
    
    Fields:
        data: Dict mapping region_id -> (mean_depth, std_depth, n_samples)
        
    Clinical use:
        - Exon deletion: z-score << 0
        - Exon amplification: z-score >> 0
        - Normal: z-score ≈ 0
        
    Note:
        Covers all exons (no filtering by variance).
        Requires minimum 3 samples for reliable statistics.
    """
    # region_id -> (mean_depth, std_depth, n_samples)
    data: Dict[str, tuple] = field(default_factory=dict)
    
    def get_stats(self, region_id: str) -> Optional[tuple]:
        """
        Get (mean, std) for a region's normalized depth.
        
        Args:
            region_id: Region identifier (chrom:start-end)
            
        Returns:
            Tuple (mean, std) or None if region not in baseline
        """
        if region_id in self.data:
            return self.data[region_id][:2]  # (mean, std), skip n_samples
        return None
    
    def compute_zscore(self, region_id: str, observed_depth: float) -> Optional[float]:
        """
        Compute z-score for observed normalized depth.
        
        Args:
            region_id: Region identifier (chrom:start-end)
            observed_depth: Sample's normalized depth for this region
            
        Returns:
            Z-score or None if region not in baseline
        """
        stats = self.get_stats(region_id)
        if stats is None:
            return None
        mean, std = stats
        if std > 0:
            return (observed_depth - mean) / std
        return 0.0
    
    def __len__(self) -> int:
        """Return number of regions in baseline."""
        return len(self.data)


@dataclass
class TfbsBaseline:
    """
    TFBS (Transcription Factor Binding Site) size entropy baseline.
    
    Stores per-TF mean/std entropy from healthy plasma samples.
    Used for z-score normalization of TFBS entropy values.
    """
    # Wraps RegionEntropyBaseline from region_entropy_processor
    baseline: 'RegionEntropyBaseline' = None
    
    def get_zscore(self, label: str, observed_entropy: float) -> float:
        """Compute z-score for observed entropy value."""
        if self.baseline is None:
            return 0.0
        return self.baseline.get_zscore(label, observed_entropy)
    
    def get_stats(self, label: str) -> Optional[tuple]:
        """Get (mean, std) for a TF label."""
        if self.baseline is None or label not in self.baseline.data:
            return None
        return self.baseline.data[label]


@dataclass
class AtacBaseline:
    """
    ATAC (Cancer ATAC-seq Peak) size entropy baseline.
    
    Stores per-cancer-type mean/std entropy from healthy plasma samples.
    Used for z-score normalization of ATAC entropy values.
    """
    # Wraps RegionEntropyBaseline from region_entropy_processor
    baseline: 'RegionEntropyBaseline' = None
    
    def get_zscore(self, label: str, observed_entropy: float) -> float:
        """Compute z-score for observed entropy value."""
        if self.baseline is None:
            return 0.0
        return self.baseline.get_zscore(label, observed_entropy)
    
    def get_stats(self, label: str) -> Optional[tuple]:
        """Get (mean, std) for a cancer type label."""
        if self.baseline is None or label not in self.baseline.data:
            return None
        return self.baseline.data[label]


@dataclass
class PonModel:
    """
    Unified Panel of Normals model for Krewlyzer.
    
    Contains all baselines needed for hybrid correction:
    - GC bias curves (FSC, FSR, WPS)
    - FSD baseline (per-arm size distributions)
    - WPS baseline (per-region mean/std)
    - OCF baseline (per-region open chromatin footprinting)
    - MDS baseline (k-mer frequencies and motif diversity score)
    
    For panel mode (panel_mode=True), also includes:
    - gc_bias_ontarget: GC curves from on-target fragments
    - fsd_baseline_ontarget: FSD from on-target fragments
    """
    schema_version: str = "1.0"
    assay: str = ""  # e.g., "msk-access-v2"
    build_date: str = ""
    n_samples: int = 0
    reference: str = ""  # e.g., "hg19"
    panel_mode: bool = False  # True if built with --target-regions
    target_regions_file: str = ""  # Original target regions BED file name
    
    # Off-target baselines (primary - always present)
    gc_bias: Optional[GcBiasModel] = None
    fsd_baseline: Optional[FsdBaseline] = None
    wps_baseline: Optional[WpsBaseline] = None
    wps_background_baseline: Optional[WpsBackgroundBaseline] = None  # Alu periodicity
    ocf_baseline: Optional[OcfBaseline] = None
    mds_baseline: Optional[MdsBaseline] = None
    tfbs_baseline: Optional[TfbsBaseline] = None  # TFBS size entropy
    atac_baseline: Optional[AtacBaseline] = None  # ATAC size entropy
    region_mds: Optional[RegionMdsBaseline] = None  # Per-gene MDS baseline
    fsc_gene_baseline: Optional[FscGeneBaseline] = None  # Per-gene normalized depth (panel mode)
    fsc_region_baseline: Optional[FscRegionBaseline] = None  # Per-exon normalized depth (panel mode)
    
    # On-target baselines (panel mode only)
    gc_bias_ontarget: Optional[GcBiasModel] = None
    fsd_baseline_ontarget: Optional[FsdBaseline] = None
    
    @classmethod
    def load(cls, path: Path) -> "PonModel":
        """
        Load PON model from Parquet file.
        
        Args:
            path: Path to .pon.parquet file
            
        Returns:
            Loaded PonModel instance
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"PON model not found: {path}")
        
        # Read all tables
        df_all = pd.read_parquet(path)
        
        # Split by table type
        metadata_df = df_all[df_all["table"] == "metadata"]
        gc_bias_df = df_all[df_all["table"] == "gc_bias"]
        fsd_df = df_all[df_all["table"] == "fsd_baseline"]
        wps_df = df_all[df_all["table"] == "wps_baseline"]
        
        if metadata_df.empty:
            raise ValueError("Invalid PON file: missing metadata table")
        
        # Parse metadata
        meta = metadata_df.iloc[0]
        
        # Parse GC bias model
        gc_bias = None
        if not gc_bias_df.empty:
            gc_bias = GcBiasModel(
                gc_bins=gc_bias_df["gc_bin"].tolist(),
                short_expected=gc_bias_df["short_expected"].tolist(),
                short_std=gc_bias_df["short_std"].tolist(),
                intermediate_expected=gc_bias_df["intermediate_expected"].tolist(),
                intermediate_std=gc_bias_df["intermediate_std"].tolist(),
                long_expected=gc_bias_df["long_expected"].tolist(),
                long_std=gc_bias_df["long_std"].tolist(),
            )
        
        # Parse FSD baseline
        fsd_baseline = None
        if not fsd_df.empty:
            arms_dict = {}
            for arm_name, arm_group in fsd_df.groupby("arm"):
                arm_group = arm_group.sort_values("size_bin")
                arms_dict[arm_name] = {
                    "expected": arm_group["expected"].tolist(),
                    "std": arm_group["std"].tolist(),
                }
            size_bins = sorted(fsd_df["size_bin"].unique().tolist())
            fsd_baseline = FsdBaseline(size_bins=size_bins, arms=arms_dict)
        
        # Parse WPS baseline (new format: wps_nuc_*, wps_tf_*)
        wps_baseline = None
        if not wps_df.empty:
            if "wps_nuc_mean" in wps_df.columns:
                # Standard columns: wps_nuc_mean/std, wps_tf_mean/std
                regions_df = wps_df[["region_id", "wps_nuc_mean", "wps_nuc_std", 
                                    "wps_tf_mean", "wps_tf_std"]].copy()
                # Rename to WpsBaseline expected format
                regions_df = regions_df.rename(columns={
                    "wps_nuc_mean": "wps_long_mean",
                    "wps_nuc_std": "wps_long_std",
                    "wps_tf_mean": "wps_short_mean",
                    "wps_tf_std": "wps_short_std",
                })
                wps_baseline = WpsBaseline(regions=regions_df)
            else:
                logger.warning("WPS baseline missing expected columns (wps_nuc_*), skipping")
        
        # Parse OCF baseline
        ocf_df = df_all[df_all["table"] == "ocf_baseline"]
        ocf_baseline = None
        if not ocf_df.empty:
            ocf_cols = ["region_id", "ocf_mean", "ocf_std"]
            if "sync_mean" in ocf_df.columns:
                ocf_cols.extend(["sync_mean", "sync_std"])
            regions_df = ocf_df[ocf_cols].copy()
            ocf_baseline = OcfBaseline(regions=regions_df)
        
        # Parse MDS baseline
        mds_df = df_all[df_all["table"] == "mds_baseline"]
        mds_baseline = None
        if not mds_df.empty:
            row = mds_df.iloc[0]
            # Parse k-mer frequencies from JSON columns if present
            kmer_expected = {}
            kmer_std = {}
            if "kmer_expected" in mds_df.columns:
                kmer_data = row.get("kmer_expected", {})
                if isinstance(kmer_data, dict):
                    kmer_expected = kmer_data
            if "kmer_std" in mds_df.columns:
                kmer_data = row.get("kmer_std", {})
                if isinstance(kmer_data, dict):
                    kmer_std = kmer_data
            mds_baseline = MdsBaseline(
                kmer_expected=kmer_expected,
                kmer_std=kmer_std,
                mds_mean=float(row.get("mds_mean", 0)),
                mds_std=float(row.get("mds_std", 1))
            )
        
        # Parse on-target FSD baseline (panel mode)
        fsd_on_df = df_all[df_all["table"] == "fsd_baseline_ontarget"]
        fsd_baseline_ontarget = None
        if not fsd_on_df.empty:
            size_bins = sorted(fsd_on_df["size_bin"].unique().tolist())
            arms_on = {}
            for arm in fsd_on_df["arm"].unique():
                arm_data = fsd_on_df[fsd_on_df["arm"] == arm].sort_values("size_bin")
                arms_on[arm] = {
                    "expected": arm_data["expected"].tolist(),
                    "std": arm_data["std"].tolist()
                }
            fsd_baseline_ontarget = FsdBaseline(size_bins=size_bins, arms=arms_on)
        
        # Parse on-target GC bias (panel mode)
        gc_on_df = df_all[df_all["table"] == "gc_bias_ontarget"]
        gc_bias_ontarget = None
        if not gc_on_df.empty:
            gc_on_df = gc_on_df.sort_values("gc_bin")
            gc_bias_ontarget = GcBiasModel(
                gc_bins=gc_on_df["gc_bin"].tolist(),
                short_expected=gc_on_df["short_expected"].tolist(),
                short_std=gc_on_df["short_std"].tolist(),
                intermediate_expected=gc_on_df["intermediate_expected"].tolist(),
                intermediate_std=gc_on_df["intermediate_std"].tolist(),
                long_expected=gc_on_df["long_expected"].tolist(),
                long_std=gc_on_df["long_std"].tolist()
            )
        
        # Parse TFBS baseline
        tfbs_df = df_all[df_all["table"] == "tfbs_baseline"]
        tfbs_baseline = None
        if not tfbs_df.empty:
            from krewlyzer.core.region_entropy_processor import RegionEntropyBaseline
            tfbs_data = {}
            for _, row in tfbs_df.iterrows():
                label = row["label"]
                mean = float(row["entropy_mean"])
                std = float(row["entropy_std"])
                tfbs_data[label] = (mean, std)
            tfbs_baseline = TfbsBaseline(baseline=RegionEntropyBaseline(tfbs_data))
            logger.debug(f"Loaded TFBS baseline: {len(tfbs_data)} labels")
        
        # Parse ATAC baseline
        atac_df = df_all[df_all["table"] == "atac_baseline"]
        atac_baseline = None
        if not atac_df.empty:
            from krewlyzer.core.region_entropy_processor import RegionEntropyBaseline
            atac_data = {}
            for _, row in atac_df.iterrows():
                label = row["label"]
                mean = float(row["entropy_mean"])
                std = float(row["entropy_std"])
                atac_data[label] = (mean, std)
            atac_baseline = AtacBaseline(baseline=RegionEntropyBaseline(atac_data))
            logger.debug(f"Loaded ATAC baseline: {len(atac_data)} labels")
        
        # Parse FSC gene baseline (panel mode)
        fsc_gene_df = df_all[df_all["table"] == "fsc_gene_baseline"]
        fsc_gene_baseline = None
        if not fsc_gene_df.empty:
            fsc_gene_data = {}
            for _, row in fsc_gene_df.iterrows():
                gene = row["gene"]
                mean = float(row["depth_mean"])
                std = float(row["depth_std"])
                n_samples = int(row.get("n_samples", 0))
                fsc_gene_data[gene] = (mean, std, n_samples)
            fsc_gene_baseline = FscGeneBaseline(data=fsc_gene_data)
            logger.debug(f"Loaded FSC gene baseline: {len(fsc_gene_data)} genes")
        
        # Parse FSC region baseline (panel mode)
        fsc_region_df = df_all[df_all["table"] == "fsc_region_baseline"]
        fsc_region_baseline = None
        if not fsc_region_df.empty:
            fsc_region_data = {}
            for _, row in fsc_region_df.iterrows():
                region_id = row["region_id"]
                mean = float(row["depth_mean"])
                std = float(row["depth_std"])
                n_samples = int(row.get("n_samples", 0))
                fsc_region_data[region_id] = (mean, std, n_samples)
            fsc_region_baseline = FscRegionBaseline(data=fsc_region_data)
            logger.debug(f"Loaded FSC region baseline: {len(fsc_region_data)} regions")
        
        # Build model
        model = cls(
            schema_version=str(meta.get("schema_version", "1.0")),
            assay=str(meta.get("assay", "")),
            build_date=str(meta.get("build_date", "")),
            n_samples=int(meta.get("n_samples", 0)),
            reference=str(meta.get("reference", "")),
            panel_mode=bool(meta.get("panel_mode", False)),
            target_regions_file=str(meta.get("target_regions_file", "")),
            gc_bias=gc_bias,
            fsd_baseline=fsd_baseline,
            wps_baseline=wps_baseline,
            ocf_baseline=ocf_baseline,
            mds_baseline=mds_baseline,
            fsd_baseline_ontarget=fsd_baseline_ontarget,
            gc_bias_ontarget=gc_bias_ontarget,
            tfbs_baseline=tfbs_baseline,
            atac_baseline=atac_baseline,
            fsc_gene_baseline=fsc_gene_baseline,
            fsc_region_baseline=fsc_region_baseline,
        )
        
        logger.info(f"Loaded PON model: {model.assay} (n={model.n_samples})")
        if model.panel_mode:
            logger.info(f"  Panel mode: ON (targets: {model.target_regions_file})")
            if fsd_baseline_ontarget:
                logger.info(f"  FSD on-target: {len(fsd_baseline_ontarget.arms)} arms")
            if gc_bias_ontarget:
                logger.info(f"  GC on-target: {len(gc_bias_ontarget.gc_bins)} bins")
        if gc_bias:
            logger.info(f"  GC bias: {len(gc_bias.gc_bins)} bins")
        if fsd_baseline:
            logger.info(f"  FSD baseline: {len(fsd_baseline.arms)} arms, {len(fsd_baseline.size_bins)} size bins")
        if wps_baseline:
            logger.info(f"  WPS baseline: {len(wps_baseline.regions)} regions")
        if ocf_baseline:
            logger.info(f"  OCF baseline: {len(ocf_baseline.regions)} regions")
        if mds_baseline:
            logger.info(f"  MDS baseline: {len(mds_baseline.kmer_expected)} k-mers")
        
        return model
    
    @classmethod
    def _load_from_single_table(cls, df: pd.DataFrame) -> "PonModel":
        """Load from legacy single-table format."""
        # Placeholder for legacy support
        return cls()
    
    def save(self, path: Path) -> None:
        """
        Save PON model to Parquet file.
        
        Args:
            path: Output path (should end with .pon.parquet)
        """
        path = Path(path)
        
        # Build metadata table
        metadata = pd.DataFrame([{
            "table": "metadata",
            "schema_version": self.schema_version,
            "assay": self.assay,
            "build_date": self.build_date,
            "n_samples": self.n_samples,
            "reference": self.reference,
        }])
        
        # Note: This is a simplified save method. build.py uses _save_pon_model()
        # for complete serialization including all baselines.
        
        # For now, just save metadata
        metadata.to_parquet(path, index=False)
        logger.info(f"Saved PON model: {path}")
    
    def validate(self) -> List[str]:
        """
        Validate PON model completeness.
        
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        if not self.assay:
            errors.append("Missing assay name")
        if not self.reference:
            errors.append("Missing reference genome")
        if self.n_samples < 1:
            errors.append("n_samples must be >= 1")
        if self.gc_bias is None:
            errors.append("Missing gc_bias model")
        if self.fsd_baseline is None:
            errors.append("Missing fsd_baseline")
        if self.wps_baseline is None:
            errors.append("Missing wps_baseline")
            
        return errors
    
    def check_assay_compatibility(self, sample_assay: Optional[str] = None) -> None:
        """
        Check if sample assay matches PON assay. Logs warning if mismatch.
        
        Args:
            sample_assay: Sample's assay identifier (optional)
        """
        if sample_assay and sample_assay != self.assay:
            logger.warning(
                f"PON model built for {self.assay}, sample may be from different assay ({sample_assay})"
            )
    
    def get_mean(self, channel: str) -> Optional[float]:
        """
        Get expected mean coverage for a fragment size channel.
        
        Used by FSC/FSR processors for log-ratio normalization.
        Returns the expected value at median GC (0.45) from GC bias curves.
        
        Args:
            channel: One of 'short', 'intermediate', 'long', 'ultra_short', 
                     'core_short', 'mono_nucl', 'di_nucl'
        
        Returns:
            Expected mean coverage (1.0 = no bias), or None if not available
        """
        if self.gc_bias is None:
            return None
        
        # Map FSC channels to GC bias model channels
        channel_map = {
            'ultra_short': 'short',
            'core_short': 'short', 
            'short': 'short',
            'mono_nucl': 'intermediate',
            'intermediate': 'intermediate',
            'di_nucl': 'long',
            'long': 'long',
        }
        
        gc_channel = channel_map.get(channel, channel)
        
        # Return expected at median GC (0.45)
        return self.gc_bias.get_expected(0.45, gc_channel)
    
    def get_variance(self, channel: str) -> Optional[float]:
        """
        Get variance for a fragment size channel from PoN samples.
        
        Used for reliability scoring: reliability = 1 / (variance + k)
        
        Args:
            channel: Fragment size channel name
            
        Returns:
            Variance across PoN samples, or None if not available
        """
        if self.gc_bias is None:
            return None
        
        # Map channels to GC bias std arrays
        channel_map = {
            'ultra_short': 'short', 'core_short': 'short', 'short': 'short',
            'mono_nucl': 'intermediate', 'intermediate': 'intermediate',
            'di_nucl': 'long', 'long': 'long',
        }
        gc_channel = channel_map.get(channel, channel)
        
        # Get appropriate std array
        std_map = {
            'short': self.gc_bias.short_std,
            'intermediate': self.gc_bias.intermediate_std,
            'long': self.gc_bias.long_std,
        }
        
        std_list = std_map.get(gc_channel)
        if std_list and len(std_list) > 0:
            # Return median std squared as variance
            median_idx = len(std_list) // 2
            return std_list[median_idx] ** 2
        
        return None
