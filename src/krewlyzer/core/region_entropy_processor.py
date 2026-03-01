"""
Region Entropy Processor for krewlyzer.

Post-processes region entropy output from Rust:
- Applies PON Z-score normalization
- Loads/saves TSV output

Used for TFBS (transcription factor binding site) and ATAC (cancer peak) size entropy.
"""

import pandas as pd
from pathlib import Path
from typing import Optional, Dict, Tuple
import logging

logger = logging.getLogger("krewlyzer.core.region_entropy_processor")


class RegionEntropyBaseline:
    """
    Baseline for TFBS or ATAC size entropy.

    Stores per-label mean and std of entropy values computed from PON cohort.
    Used to convert raw entropy to Z-score for normalized comparison.
    """

    def __init__(self, data: Dict[str, Tuple[float, float]]):
        """
        Initialize baseline with pre-computed statistics.

        Args:
            data: Dict[label, (mean_entropy, std_entropy)]
        """
        self.data = data

    def get_zscore(self, label: str, value: float) -> float:
        """
        Compute Z-score for a given label and entropy value.

        Args:
            label: Region label (e.g., TF name or cancer type)
            value: Raw entropy value

        Returns:
            Z-score: (value - mean) / std
        """
        if label not in self.data:
            return 0.0
        mean, std = self.data[label]
        if std <= 1e-6:
            return 0.0
        return (value - mean) / std

    @classmethod
    def from_samples(cls, sample_data: list) -> "RegionEntropyBaseline":
        """
        Build baseline from multiple samples.

        Args:
            sample_data: List of dicts, each containing label -> entropy value

        Returns:
            RegionEntropyBaseline with mean/std per label
        """
        import numpy as np

        # Collect all values per label
        label_values: Dict[str, list] = {}
        for sample in sample_data:
            for label, entropy in sample.items():
                if label not in label_values:
                    label_values[label] = []
                label_values[label].append(entropy)

        # Compute mean/std per label
        data = {}
        for label, values in label_values.items():
            if len(values) >= 2:
                mean = float(np.mean(values))
                std = float(np.std(values, ddof=1))
                data[label] = (mean, std)
            elif len(values) == 1:
                # Single sample: use value as mean, set std=0 (no normalization)
                data[label] = (float(values[0]), 0.0)

        logger.info(f"Built RegionEntropyBaseline with {len(data)} labels")
        return cls(data)


def load_entropy_tsv(path: Path) -> pd.DataFrame:
    """
    Load region entropy TSV file.

    Columns: label, count, mean_size, entropy

    Args:
        path: Path to TSV file

    Returns:
        DataFrame with entropy data
    """
    return pd.read_csv(path, sep="\t")


def process_region_entropy(
    raw_path: Path,
    output_path: Path,
    pon_parquet_path: Optional[Path] = None,
    baseline_table: str = "tfbs_baseline",
) -> int:
    """
    Apply PON normalization to region entropy output using Rust implementation.

    Adds 'z_score' column: (raw_entropy - pon_mean) / pon_std

    Args:
        raw_path: Path to raw entropy TSV from Rust
        output_path: Path to write normalized output
        pon_parquet_path: Path to PON Parquet file
        baseline_table: "tfbs_baseline" or "atac_baseline"

    Returns:
        Number of labels with z-scores computed
    """
    from krewlyzer import _core

    if not raw_path.exists():
        logger.warning(f"Entropy file not found: {raw_path}")
        return 0

    if pon_parquet_path is None:
        # No PON - just copy file and add z_score=0 column
        df = load_entropy_tsv(raw_path)
        df["z_score"] = 0.0
        df.to_csv(output_path, sep="\t", index=False, float_format="%.4f")
        logger.debug(f"No PON baseline, wrote {len(df)} labels with z_score=0")
        return 0

    n_matched = _core.region_entropy.apply_pon_zscore(
        str(raw_path), str(pon_parquet_path), str(output_path), baseline_table
    )

    logger.info(
        f"Region entropy z-scores: {n_matched} labels matched ({output_path.name})"
    )
    return n_matched


def extract_entropy_data(df: pd.DataFrame) -> Dict[str, float]:
    """
    Extract entropy values from DataFrame for PON baseline building.

    Args:
        df: DataFrame with columns: label, entropy

    Returns:
        Dict[label, entropy] for this sample
    """
    return {row["label"]: row["entropy"] for _, row in df.iterrows()}
