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
        
        # Parse WPS baseline
        wps_baseline = None
        if not wps_df.empty:
            regions_df = wps_df[["region_id", "wps_long_mean", "wps_long_std", 
                                "wps_short_mean", "wps_short_std"]].copy()
            wps_baseline = WpsBaseline(regions=regions_df)
        
        # Build model
        model = cls(
            schema_version=str(meta.get("schema_version", "1.0")),
            assay=str(meta.get("assay", "")),
            build_date=str(meta.get("build_date", "")),
            n_samples=int(meta.get("n_samples", 0)),
            reference=str(meta.get("reference", "")),
            gc_bias=gc_bias,
            fsd_baseline=fsd_baseline,
            wps_baseline=wps_baseline,
        )
        
        logger.info(f"Loaded PON model: {model.assay} (n={model.n_samples})")
        if gc_bias:
            logger.info(f"  GC bias: {len(gc_bias.gc_bins)} bins")
        if fsd_baseline:
            logger.info(f"  FSD baseline: {len(fsd_baseline.arms)} arms, {len(fsd_baseline.size_bins)} size bins")
        if wps_baseline:
            logger.info(f"  WPS baseline: {len(wps_baseline.regions)} regions")
        
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
