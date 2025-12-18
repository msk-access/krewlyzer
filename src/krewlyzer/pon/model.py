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


@dataclass
class WpsBaseline:
    """
    WPS baseline per transcript region.
    
    Stores mean and std WPS for each region (both long and short).
    """
    regions: pd.DataFrame  # Columns: region_id, wps_long_mean, wps_long_std, wps_short_mean, wps_short_std
    
    def get_baseline(self, region_id: str) -> Optional[Dict[str, float]]:
        """Get baseline stats for a region."""
        match = self.regions[self.regions["region_id"] == region_id]
        if match.empty:
            return None
        row = match.iloc[0]
        return {
            "wps_long_mean": row["wps_long_mean"],
            "wps_long_std": row["wps_long_std"],
            "wps_short_mean": row["wps_short_mean"],
            "wps_short_std": row["wps_short_std"],
        }


@dataclass
class PonModel:
    """
    Unified Panel of Normals model for Krewlyzer.
    
    Contains all baselines needed for hybrid correction:
    - GC bias curves (FSC, FSR, WPS)
    - FSD baseline (per-arm size distributions)
    - WPS baseline (per-region mean/std)
    """
    schema_version: str = "1.0"
    assay: str = ""  # e.g., "msk-access-v2"
    build_date: str = ""
    n_samples: int = 0
    reference: str = ""  # e.g., "hg19"
    
    gc_bias: Optional[GcBiasModel] = None
    fsd_baseline: Optional[FsdBaseline] = None
    wps_baseline: Optional[WpsBaseline] = None
    
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
        
        # Read metadata table
        metadata_df = pd.read_parquet(path, filters=[("table", "==", "metadata")])
        if metadata_df.empty:
            # Try reading as single table (legacy format)
            df = pd.read_parquet(path)
            return cls._load_from_single_table(df)
        
        # Read component tables
        gc_bias_df = pd.read_parquet(path, filters=[("table", "==", "gc_bias")])
        fsd_df = pd.read_parquet(path, filters=[("table", "==", "fsd_baseline")])
        wps_df = pd.read_parquet(path, filters=[("table", "==", "wps_baseline")])
        
        # Parse metadata
        meta = metadata_df.iloc[0]
        
        # Build model
        model = cls(
            schema_version=meta.get("schema_version", "1.0"),
            assay=meta.get("assay", ""),
            build_date=meta.get("build_date", ""),
            n_samples=int(meta.get("n_samples", 0)),
            reference=meta.get("reference", ""),
        )
        
        # TODO: Parse gc_bias, fsd_baseline, wps_baseline from DataFrames
        logger.info(f"Loaded PON model: {model.assay} (n={model.n_samples})")
        
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
        
        # TODO: Build gc_bias, fsd_baseline, wps_baseline tables
        
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
