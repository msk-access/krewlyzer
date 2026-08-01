"""
PON Model definitions for Krewlyzer.

This module defines the unified PON model format using Parquet storage.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional
import logging

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    # RegionEntropyBaseline is defined in region_entropy_processor; imported here
    # only for type checking to avoid a circular import at runtime. TfbsBaseline
    # and AtacBaseline reference it for their .baseline field annotations.
    from krewlyzer.core.region_entropy_processor import RegionEntropyBaseline

logger = logging.getLogger("pon")


#: What a z-score is when there is nothing to divide by.
#:
#: Every baseline class below independently ended `if std > 0: ... return 0.0`,
#: nine times over. Zero is not a cautious answer there -- it is the single
#: most confident statement the column can make ("this sample sits exactly at
#: the healthy baseline"), asserted precisely when the baseline measured no
#: spread at all, and indistinguishable from a genuine zero.
#:
#: The same reasoning that took the sigma floors out of the builder (4cd634b)
#: and `z_score = 0.0` out of region entropy: a value a reader cannot tell
#: apart from a measurement must be a measurement, or absent.
def zscore_or_nan(
    observed: Optional[float], mean: Optional[float], std: Optional[float]
) -> float:
    """``(observed - mean) / std``, or NaN when ``std`` is not usable.

    NaN propagates to an absent column value rather than a fabricated zero,
    and — unlike zero — cannot be mistaken for a reading.
    """
    if std is None or not np.isfinite(std) or std <= 0:
        return float("nan")
    if mean is None or observed is None:
        return float("nan")
    if not np.isfinite(mean) or not np.isfinite(observed):
        return float("nan")
    return (observed - mean) / std


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
        """Expected proportion for a size bin, NaN when the arm is unknown.

        Not 0.0: an arm absent from the baseline has no expectation, and zero
        is a specific -- and specifically wrong -- one.
        """
        if arm not in self.arms:
            return float("nan")
        arm_data = self.arms[arm]
        return float(np.interp(size, self.size_bins, arm_data["expected"]))

    def get_std(self, arm: str, size: int) -> float:
        """Spread for a size bin, NaN when the arm is unknown.

        Zero would be worse than useless here -- a caller dividing by it gets
        infinity rather than an absent value.
        """
        if arm not in self.arms:
            return float("nan")
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

    def get_baseline_vector(
        self, region_id: str, column: str = "wps_nuc"
    ) -> Optional[np.ndarray]:
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

    def get_std_vector(
        self, region_id: str, column: str = "wps_nuc"
    ) -> Optional[np.ndarray]:
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

    def compute_z_vector(
        self, region_id: str, sample_vector: np.ndarray, column: str = "wps_nuc"
    ) -> Optional[np.ndarray]:
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

        # NaN where sigma is not usable, never a substituted 1.0.
        #
        # Since 4cd634b the builder writes NaN for positions whose spread it
        # could not measure. `np.where(std > 0, ...)` is False for NaN, so a
        # 1.0 default would turn every one of those into a finite z -- undoing
        # the builder's honesty at the read side, which is the harder place to
        # notice it.
        with np.errstate(divide="ignore", invalid="ignore"):
            usable = np.isfinite(std) & (std > 0)
            return np.where(
                usable, (sample_vector - mean) / np.where(usable, std, 1.0), np.nan
            )

    def compute_shape_score(
        self, region_id: str, sample_vector: np.ndarray, column: str = "wps_nuc"
    ) -> Optional[float]:
        """
        Compute shape correlation score for sample vs PoN (v2.0 format).

        Returns Pearson correlation coefficient between sample vector and PON mean.
        Used for cancer detection - healthy samples have correlation ~1.0.

        Args:
            region_id: Region identifier
            sample_vector: 200-element sample WPS vector
            column: Vector column prefix ('wps_nuc' or 'wps_tf')

        Returns:
            Correlation coefficient [-1, 1] or None if baseline not found

        Clinical interpretation:
            - 0.9-1.0: Healthy nucleosome positioning
            - 0.5-0.9: Mild chromatin disorganization
            - <0.5: Significant disruption (cancer signal)
        """
        mean = self.get_baseline_vector(region_id, column)
        if mean is None or len(mean) == 0:
            return None

        if len(sample_vector) != len(mean):
            return None

        # Compute Pearson correlation
        # Handle edge cases where std is 0
        sample_std = np.std(sample_vector)
        mean_std = np.std(mean)

        if sample_std < 1e-6 or mean_std < 1e-6:
            # Undefined, not uncorrelated. Zero is a real claim about shape
            # agreement and would be indistinguishable from having measured it.
            return float("nan")

        correlation = np.corrcoef(sample_vector, mean)[0, 1]
        return float(correlation) if not np.isnan(correlation) else 0.0


#: The derived WPS shape quantities, and what each answers.
#:
#: Each is z-scored against its own mean/sigma, never derived from the
#: per-position z vector. Adjacent WPS positions have lag-1 autocorrelation
#: 0.986 -- a fragment spans ~167 bp and touches many positions at once -- so
#: an average of z across positions has none of a z-score's properties.
WPS_SHAPE_STATS = ("log_amplitude", "shape_corr_fisher")


@dataclass
class WpsShapeBaseline:
    """Per-anchor baselines for the three derived WPS shape quantities.

    Separate from :class:`WpsBaseline`, which holds the 200-element mean and
    sigma *profiles*. These are scalars per anchor, and they answer questions a
    per-position comparison cannot:

    ``log_amplitude``       is there nucleosome structure here at all
    ``shape_corr_fisher``   is it the *right* structure

    Positional displacement is measured too, but deliberately not baselined --
    see ``core/wps_pon.py``. It showed no per-sample signal and an intraclass
    correlation of 0.479 per anchor, so a z-score of it would be a plausible
    number with nothing behind it.

    All three are window-free by design. Measured on the real cohort, TSS
    anchors dip at the centre (-6.8 against -3.4 in the flanks) while CTCF
    anchors do the opposite, so any fixed centre-versus-flank definition is
    backwards for one of the two.
    """

    #: region_id -> {f"{stat}_mean": float, f"{stat}_std": float, "n_samples": int}
    regions: Dict[str, Dict[str, float]] = field(default_factory=dict)

    def compute_zscore(
        self, region_id: str, stat: str, observed: float
    ) -> Optional[float]:
        """Z for one derived quantity at one anchor.

        ``None`` when the anchor is absent from the baseline; NaN when it is
        present but the baseline measured no usable spread.
        """
        entry = self.regions.get(region_id)
        if entry is None:
            return None
        return zscore_or_nan(
            observed, entry.get(f"{stat}_mean"), entry.get(f"{stat}_std")
        )


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
        return zscore_or_nan(observed_ocf, mean, std)


#: The group id the whole-sample NRL and periodicity are keyed on.
#:
#: The default here used to be the string ``"all"``, which no PON has ever
#: contained -- the groups are ``Global_All``, ``Chr1_All`` ... ``Family_AluY``.
#: Every lookup therefore missed, ``compute_nrl_zscore`` returned None, and
#: ``nrl_z``/``periodicity_z`` never reached a single output file. Silent,
#: because a None was indistinguishable from "this PON has no baseline".
GENOME_WIDE_GROUP = "Global_All"


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

    def get_nrl_stats(self, group_id: str = GENOME_WIDE_GROUP) -> Optional[tuple]:
        """
        Get (mean, std) for nucleosome repeat length.

        Args:
            group_id: Group identifier (default: the genome-wide group)

        Returns:
            Tuple (nrl_mean, nrl_std) or None if not found
        """
        match = self.groups[self.groups["group_id"] == group_id]
        if match.empty:
            # Loud: a group id that matches nothing is a naming mistake, not a
            # PON without a baseline, and the two used to look identical.
            logger.warning(
                f"WPS background baseline has no group {group_id!r}; "
                f"available: {sorted(self.groups['group_id'].unique())[:5]}..."
            )
            return None
        row = match.iloc[0]
        return (row["nrl_mean"], row["nrl_std"])

    def get_periodicity_stats(
        self, group_id: str = GENOME_WIDE_GROUP
    ) -> Optional[tuple]:
        """
        Get (mean, std) for periodicity score.

        Args:
            group_id: Group identifier (default: the genome-wide group)

        Returns:
            Tuple (periodicity_mean, periodicity_std) or None if not found
        """
        match = self.groups[self.groups["group_id"] == group_id]
        if match.empty:
            return None
        row = match.iloc[0]
        return (row.get("periodicity_mean", 0), row.get("periodicity_std", 1))

    def compute_nrl_zscore(
        self, observed_nrl: float, group_id: str = GENOME_WIDE_GROUP
    ) -> Optional[float]:
        """Compute z-score for observed NRL value."""
        stats = self.get_nrl_stats(group_id)
        if stats is None:
            return None
        mean, std = stats
        return zscore_or_nan(observed_nrl, mean, std)


@dataclass
class MdsBaseline:
    """
    Motif Diversity Score baseline from healthy plasma samples.

    Stores expected k-mer frequencies and MDS value for healthy samples.
    Enables detection of aberrant end-motif patterns in cancer cfDNA.

    Note: DNA damage and tissue-of-origin affect motif patterns.
    """

    # Expected k-mer frequencies (256 4-mers)
    kmer_expected: Dict[str, float] = field(
        default_factory=dict
    )  # e.g., {"ACGT": 0.0042, ...}
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
        return zscore_or_nan(observed_freq, expected, std)

    def get_mds_zscore(self, observed_mds: float) -> float:
        """Compute z-score for MDS value."""
        return zscore_or_nan(observed_mds, self.mds_mean, self.mds_std)

    def get_aberrant_kmers(
        self, observed_freqs: Dict[str, float], threshold: float = 2.0
    ) -> Dict[str, float]:
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
        return zscore_or_nan(observed_mds, mean, std)

    def compute_e1_zscore(self, gene: str, observed_mds_e1: float) -> Optional[float]:
        """
        Compute z-score for observed E1 MDS value.
        """
        stats = self.get_e1_stats(gene)
        if stats is None:
            return None
        mean, std = stats
        return zscore_or_nan(observed_mds_e1, mean, std)


@dataclass
class RegionMdsExonBaseline:
    """Per-exon MDS baseline — the finest localisation krewlyzer produces.

    Keyed on ``(gene, name)``. ``name`` alone is not unique across genes, and
    coordinates would break whenever the panel BED is regenerated.

    Measured on the 0.8.3 corpus before this was built: every exon appears in
    every sample of its assay (7/7 for xs1, 19/19 for xs2) with a measurable
    spread, and under 0.25% carry fewer than 10 fragments. Exons are far
    better supported here than "per-exon" suggests, so this needs no special
    sparsity handling — only the same NaN-not-floor rule as everywhere else.
    """

    #: (gene, name) -> {mds_mean, mds_std, n_samples}
    exon_baseline: Dict[tuple, Dict[str, float]] = field(default_factory=dict)

    def get_stats(self, gene: str, name: str) -> Optional[tuple]:
        """``(mean, std)`` for one exon, or None when it is not in the baseline."""
        entry = self.exon_baseline.get((gene, name))
        if entry is None:
            return None
        return (entry.get("mds_mean"), entry.get("mds_std"))

    def compute_zscore(
        self, gene: str, name: str, observed_mds: float
    ) -> Optional[float]:
        """Z-score for one exon.

        ``None`` when the exon is absent from the baseline; NaN when it is
        present but the baseline measured no usable spread. Both write NaN to
        the output, but only the caller's log can tell them apart.
        """
        stats = self.get_stats(gene, name)
        if stats is None:
            return None
        mean, std = stats
        return zscore_or_nan(observed_mds, mean, std)


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
        return zscore_or_nan(observed_depth, mean, std)

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
        return zscore_or_nan(observed_depth, mean, std)

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

    # baseline wraps RegionEntropyBaseline from region_entropy_processor.
    # Import is under TYPE_CHECKING to avoid circular deps at runtime.
    baseline: Optional["RegionEntropyBaseline"] = None

    @property
    def labels(self) -> List[str]:
        """Return list of TF labels available in this baseline."""
        if self.baseline is None:
            return []
        return list(self.baseline.data.keys())

    def get_zscore(self, label: str, observed_entropy: float) -> float:
        """Compute z-score for observed entropy value, NaN without a baseline."""
        if self.baseline is None:
            return float("nan")
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

    # baseline wraps RegionEntropyBaseline from region_entropy_processor.
    # Import is under TYPE_CHECKING to avoid circular deps at runtime.
    baseline: Optional["RegionEntropyBaseline"] = None

    @property
    def labels(self) -> List[str]:
        """Return list of ATAC cancer-type labels available in this baseline."""
        if self.baseline is None:
            return []
        return list(self.baseline.data.keys())

    def get_zscore(self, label: str, observed_entropy: float) -> float:
        """Compute z-score for observed entropy value, NaN without a baseline."""
        if self.baseline is None:
            return float("nan")
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
    wps_shape_baseline: Optional[WpsShapeBaseline] = None
    wps_baseline_panel: Optional[WpsBaseline] = (
        None  # Panel-specific WPS (panel mode only)
    )
    ocf_baseline: Optional[OcfBaseline] = None
    ocf_baseline_ontarget: Optional[OcfBaseline] = None  # On-target OCF (panel mode)
    ocf_baseline_offtarget: Optional[OcfBaseline] = None  # Off-target OCF (panel mode)
    mds_baseline: Optional[MdsBaseline] = None
    tfbs_baseline: Optional[TfbsBaseline] = None  # TFBS size entropy
    atac_baseline: Optional[AtacBaseline] = None  # ATAC size entropy
    region_mds: Optional[RegionMdsBaseline] = None
    region_mds_exon: Optional[RegionMdsExonBaseline] = None  # Per-exon MDS baseline
    fsc_gene_baseline: Optional[FscGeneBaseline] = (
        None  # Per-gene normalized depth (panel mode)
    )
    fsc_region_baseline: Optional[FscRegionBaseline] = (
        None  # Per-exon normalized depth (panel mode)
    )

    # On-target baselines (panel mode only)
    # TFBS/ATAC on-target uses panel-specific regions (pre-intersected with targets)
    # for higher signal in target-adjacent chromatin. Other features use fragment filtering.
    gc_bias_ontarget: Optional[GcBiasModel] = None
    fsd_baseline_ontarget: Optional[FsdBaseline] = None
    mds_baseline_ontarget: Optional[MdsBaseline] = None  # On-target k-mer MDS
    tfbs_baseline_ontarget: Optional[TfbsBaseline] = None  # Panel-specific TFBS regions
    atac_baseline_ontarget: Optional[AtacBaseline] = None  # Panel-specific ATAC regions

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
                arms_dict[str(arm_name)] = {
                    "expected": arm_group["expected"].tolist(),
                    "std": arm_group["std"].tolist(),
                }
            size_bins = sorted(fsd_df["size_bin"].unique().tolist())
            fsd_baseline = FsdBaseline(size_bins=size_bins, arms=arms_dict)

        # Parse WPS baseline - detect v2.0 vector or v1.0 scalar format
        wps_baseline = None
        if not wps_df.empty:
            if "wps_nuc_mean" in wps_df.columns:
                # Check if it's vector format (list columns) or scalar format
                first_val = wps_df["wps_nuc_mean"].iloc[0] if len(wps_df) > 0 else None
                is_vector = (
                    isinstance(first_val, (list, np.ndarray))
                    if first_val is not None
                    else False
                )

                if is_vector:
                    # v2.0 vector format - keep original column names
                    cols = [
                        "region_id",
                        "wps_nuc_mean",
                        "wps_nuc_std",
                        "wps_tf_mean",
                        "wps_tf_std",
                    ]
                    if "n_samples" in wps_df.columns:
                        cols.append("n_samples")
                    regions_df = wps_df[cols].copy()
                    wps_baseline = WpsBaseline(regions=regions_df, schema_version="2.0")
                    logger.info(
                        f"Loaded WPS baseline v2.0 (vector format): {len(regions_df)} regions"
                    )
                else:
                    # v1.0 scalar format - rename to legacy column names
                    regions_df = wps_df[
                        [
                            "region_id",
                            "wps_nuc_mean",
                            "wps_nuc_std",
                            "wps_tf_mean",
                            "wps_tf_std",
                        ]
                    ].copy()
                    regions_df = regions_df.rename(
                        columns={
                            "wps_nuc_mean": "wps_long_mean",
                            "wps_nuc_std": "wps_long_std",
                            "wps_tf_mean": "wps_short_mean",
                            "wps_tf_std": "wps_short_std",
                        }
                    )
                    wps_baseline = WpsBaseline(regions=regions_df, schema_version="1.0")
            else:
                logger.warning(
                    "WPS baseline missing expected columns (wps_nuc_*), skipping"
                )

        # Parse WPS panel baseline (panel mode only) - same format as wps_baseline
        wps_panel_df = df_all[df_all["table"] == "wps_baseline_panel"]
        wps_baseline_panel = None
        if not wps_panel_df.empty:
            if "wps_nuc_mean" in wps_panel_df.columns:
                first_val = (
                    wps_panel_df["wps_nuc_mean"].iloc[0]
                    if len(wps_panel_df) > 0
                    else None
                )
                is_vector = (
                    isinstance(first_val, (list, np.ndarray))
                    if first_val is not None
                    else False
                )

                if is_vector:
                    cols = [
                        "region_id",
                        "wps_nuc_mean",
                        "wps_nuc_std",
                        "wps_tf_mean",
                        "wps_tf_std",
                    ]
                    if "n_samples" in wps_panel_df.columns:
                        cols.append("n_samples")
                    regions_df = wps_panel_df[cols].copy()
                    wps_baseline_panel = WpsBaseline(
                        regions=regions_df, schema_version="2.0"
                    )
                    logger.info(
                        f"Loaded WPS panel baseline v2.0: {len(regions_df)} regions"
                    )
                else:
                    regions_df = wps_panel_df[
                        [
                            "region_id",
                            "wps_nuc_mean",
                            "wps_nuc_std",
                            "wps_tf_mean",
                            "wps_tf_std",
                        ]
                    ].copy()
                    regions_df = regions_df.rename(
                        columns={
                            "wps_nuc_mean": "wps_long_mean",
                            "wps_nuc_std": "wps_long_std",
                            "wps_tf_mean": "wps_short_mean",
                            "wps_tf_std": "wps_short_std",
                        }
                    )
                    wps_baseline_panel = WpsBaseline(
                        regions=regions_df, schema_version="1.0"
                    )

        # Parse OCF baseline
        ocf_df = df_all[df_all["table"] == "ocf_baseline"]
        ocf_baseline = None
        if not ocf_df.empty:
            ocf_cols = ["region_id", "ocf_mean", "ocf_std"]
            if "sync_mean" in ocf_df.columns:
                ocf_cols.extend(["sync_mean", "sync_std"])
            regions_df = ocf_df[ocf_cols].copy()
            ocf_baseline = OcfBaseline(regions=regions_df)

        # Parse on-target OCF baseline (panel mode)
        ocf_on_df = df_all[df_all["table"] == "ocf_baseline_ontarget"]
        ocf_baseline_ontarget = None
        if not ocf_on_df.empty:
            ocf_on_cols = ["region_id", "ocf_mean", "ocf_std"]
            if "sync_mean" in ocf_on_df.columns:
                ocf_on_cols.extend(["sync_mean", "sync_std"])
            ocf_baseline_ontarget = OcfBaseline(regions=ocf_on_df[ocf_on_cols].copy())
            logger.debug(f"Loaded OCF on-target baseline: {len(ocf_on_df)} regions")

        # Parse off-target OCF baseline (panel mode)
        ocf_off_df = df_all[df_all["table"] == "ocf_baseline_offtarget"]
        ocf_baseline_offtarget = None
        if not ocf_off_df.empty:
            ocf_off_cols = ["region_id", "ocf_mean", "ocf_std"]
            if "sync_mean" in ocf_off_df.columns:
                ocf_off_cols.extend(["sync_mean", "sync_std"])
            ocf_baseline_offtarget = OcfBaseline(
                regions=ocf_off_df[ocf_off_cols].copy()
            )
            logger.debug(f"Loaded OCF off-target baseline: {len(ocf_off_df)} regions")

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
                mds_std=float(row.get("mds_std", 1)),
            )

        # Parse on-target MDS baseline (panel mode)
        # Uses on-target k-mer frequencies for separate panel-mode MDS baseline
        mds_on_df = df_all[df_all["table"] == "mds_baseline_ontarget"]
        mds_baseline_ontarget = None
        if not mds_on_df.empty:
            row_on = mds_on_df.iloc[0]
            kmer_expected_on = {}
            kmer_std_on = {}
            if "kmer_expected" in mds_on_df.columns:
                kmer_data = row_on.get("kmer_expected", {})
                if isinstance(kmer_data, dict):
                    kmer_expected_on = kmer_data
            if "kmer_std" in mds_on_df.columns:
                kmer_data = row_on.get("kmer_std", {})
                if isinstance(kmer_data, dict):
                    kmer_std_on = kmer_data
            mds_baseline_ontarget = MdsBaseline(
                kmer_expected=kmer_expected_on,
                kmer_std=kmer_std_on,
                mds_mean=float(row_on.get("mds_mean", 0)),
                mds_std=float(row_on.get("mds_std", 1)),
            )
            logger.debug(
                f"Loaded MDS on-target baseline: mean={mds_baseline_ontarget.mds_mean:.4f}"
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
                    "std": arm_data["std"].tolist(),
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
                long_std=gc_on_df["long_std"].tolist(),
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

        # =====================================================================
        # Parse on-target TFBS baseline (panel mode)
        # Separate baseline for on-target entropy (different depth/variance)
        # =====================================================================
        tfbs_ontarget_df = df_all[df_all["table"] == "tfbs_baseline_ontarget"]
        tfbs_baseline_ontarget = None
        if not tfbs_ontarget_df.empty:
            from krewlyzer.core.region_entropy_processor import RegionEntropyBaseline

            tfbs_ontarget_data = {}
            for _, row in tfbs_ontarget_df.iterrows():
                label = row["label"]
                mean = float(row["entropy_mean"])
                std = float(row["entropy_std"])
                tfbs_ontarget_data[label] = (mean, std)
            tfbs_baseline_ontarget = TfbsBaseline(
                baseline=RegionEntropyBaseline(tfbs_ontarget_data)
            )
            logger.debug(
                f"Loaded TFBS on-target baseline: {len(tfbs_ontarget_data)} labels"
            )

        # Parse on-target ATAC baseline (panel mode)
        atac_ontarget_df = df_all[df_all["table"] == "atac_baseline_ontarget"]
        atac_baseline_ontarget = None
        if not atac_ontarget_df.empty:
            from krewlyzer.core.region_entropy_processor import RegionEntropyBaseline

            atac_ontarget_data = {}
            for _, row in atac_ontarget_df.iterrows():
                label = row["label"]
                mean = float(row["entropy_mean"])
                std = float(row["entropy_std"])
                atac_ontarget_data[label] = (mean, std)
            atac_baseline_ontarget = AtacBaseline(
                baseline=RegionEntropyBaseline(atac_ontarget_data)
            )
            logger.debug(
                f"Loaded ATAC on-target baseline: {len(atac_ontarget_data)} labels"
            )

        # =====================================================================
        # FSC Gene/Region Baselines (panel mode)
        # =====================================================================

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

        # Parse WPS shape baseline
        wps_shape_df = df_all[df_all["table"] == "wps_shape_baseline"]
        wps_shape_baseline = None
        if not wps_shape_df.empty:
            wanted = [
                f"{stat}_{moment}"
                for stat in WPS_SHAPE_STATS
                for moment in ("mean", "std")
            ]
            shape_regions = {}
            for _, row in wps_shape_df.iterrows():
                # No defaults: a NaN sigma means the builder could not measure
                # a spread, and zscore_or_nan must see it as NaN.
                shape_regions[str(row["region_id"])] = {
                    key: float(row[key]) for key in wanted if key in row
                }
            wps_shape_baseline = WpsShapeBaseline(regions=shape_regions)
            logger.debug(f"Loaded WPS shape baseline: {len(shape_regions)} anchors")

        # Parse Region MDS exon baseline
        region_mds_exon_df = df_all[df_all["table"] == "region_mds_exon"]
        region_mds_exon = None
        if not region_mds_exon_df.empty:
            exon_baseline = {}
            for _, row in region_mds_exon_df.iterrows():
                # No `or 0.0` defaults: a NaN sigma here means the builder
                # could not measure a spread, and zscore_or_nan must see it.
                exon_baseline[(str(row["gene"]), str(row["name"]))] = {
                    "mds_mean": float(row["mds_mean"]),
                    "mds_std": float(row["mds_std"]),
                    "n_samples": int(row.get("n_samples", 0) or 0),
                }
            region_mds_exon = RegionMdsExonBaseline(exon_baseline=exon_baseline)
            logger.debug(f"Loaded Region MDS exon baseline: {len(exon_baseline)} exons")

        # Parse Region MDS baseline
        region_mds_df = df_all[df_all["table"] == "region_mds"]
        region_mds = None
        if not region_mds_df.empty:
            gene_baseline = {}
            for _, row in region_mds_df.iterrows():
                gene = row["gene"]
                gene_baseline[gene] = {
                    "mds_mean": float(row.get("mds_mean", 0.0)),
                    "mds_std": float(row.get("mds_std", 1.0)),
                    "mds_e1_mean": (
                        float(row["mds_e1_mean"])
                        if pd.notna(row.get("mds_e1_mean"))
                        else None
                    ),
                    "mds_e1_std": (
                        float(row["mds_e1_std"])
                        if pd.notna(row.get("mds_e1_std"))
                        else None
                    ),
                    "n_samples": int(row.get("n_samples", 0)),
                }
            region_mds = RegionMdsBaseline(gene_baseline=gene_baseline)  # type: ignore[arg-type]
            logger.debug(f"Loaded Region MDS baseline: {len(gene_baseline)} genes")

        # Parse WPS Background baseline
        wps_background_df = df_all[df_all["table"] == "wps_background"]
        wps_background_baseline = None
        if not wps_background_df.empty:
            groups_df = wps_background_df.drop(columns=["table"], errors="ignore")
            wps_background_baseline = WpsBackgroundBaseline(groups=groups_df)
            logger.debug(f"Loaded WPS Background baseline: {len(groups_df)} groups")

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
            wps_shape_baseline=wps_shape_baseline,
            wps_baseline_panel=wps_baseline_panel,
            ocf_baseline=ocf_baseline,
            ocf_baseline_ontarget=ocf_baseline_ontarget,
            ocf_baseline_offtarget=ocf_baseline_offtarget,
            mds_baseline=mds_baseline,
            mds_baseline_ontarget=mds_baseline_ontarget,
            fsd_baseline_ontarget=fsd_baseline_ontarget,
            gc_bias_ontarget=gc_bias_ontarget,
            tfbs_baseline=tfbs_baseline,
            atac_baseline=atac_baseline,
            # On-target entropy baselines (panel mode - panel-specific regions)
            tfbs_baseline_ontarget=tfbs_baseline_ontarget,
            atac_baseline_ontarget=atac_baseline_ontarget,
            fsc_gene_baseline=fsc_gene_baseline,
            fsc_region_baseline=fsc_region_baseline,
            region_mds=region_mds,
            region_mds_exon=region_mds_exon,
            wps_background_baseline=wps_background_baseline,
        )

        logger.info(f"Loaded PON model: {model.assay} (n={model.n_samples})")
        if model.panel_mode:
            logger.info(f"  Panel mode: ON (targets: {model.target_regions_file})")
            if fsd_baseline_ontarget:
                logger.info(f"  FSD on-target: {len(fsd_baseline_ontarget.arms)} arms")
            if gc_bias_ontarget:
                logger.info(f"  GC on-target: {len(gc_bias_ontarget.gc_bins)} bins")
            if mds_baseline_ontarget:
                logger.info(
                    f"  MDS on-target: mean={mds_baseline_ontarget.mds_mean:.4f}"
                )
            # On-target entropy baselines
            if tfbs_baseline_ontarget:
                logger.info(
                    f"  TFBS on-target: {len(tfbs_baseline_ontarget.labels)} labels"
                )
            if atac_baseline_ontarget:
                logger.info(
                    f"  ATAC on-target: {len(atac_baseline_ontarget.labels)} labels"
                )
        if gc_bias:
            logger.info(f"  GC bias: {len(gc_bias.gc_bins)} bins")
        if fsd_baseline:
            logger.info(
                f"  FSD baseline: {len(fsd_baseline.arms)} arms, {len(fsd_baseline.size_bins)} size bins"
            )
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

    # NOTE: there is deliberately no `save()` here.
    #
    # One existed and wrote only the metadata block -- a PON with no baselines
    # at all -- while `build-pon` used `_save_pon_model` in `build.py`. Two
    # writers for one format, one of them silently producing an empty model,
    # and nothing in production called it. Removed rather than completed:
    # a second serializer is a second thing to keep in step with every new
    # block, and the first one already has to be.

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
            "ultra_short": "short",
            "core_short": "short",
            "short": "short",
            "mono_nucl": "intermediate",
            "intermediate": "intermediate",
            "di_nucl": "long",
            "long": "long",
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
            "ultra_short": "short",
            "core_short": "short",
            "short": "short",
            "mono_nucl": "intermediate",
            "intermediate": "intermediate",
            "di_nucl": "long",
            "long": "long",
        }
        gc_channel = channel_map.get(channel, channel)

        # Get appropriate std array
        std_map = {
            "short": self.gc_bias.short_std,
            "intermediate": self.gc_bias.intermediate_std,
            "long": self.gc_bias.long_std,
        }

        std_list = std_map.get(gc_channel)
        if std_list and len(std_list) > 0:
            # Return median std squared as variance
            median_idx = len(std_list) // 2
            return std_list[median_idx] ** 2

        return None
