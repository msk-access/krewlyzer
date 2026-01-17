"""
Feature Serializer for unified JSON output.

Collects all krewlyzer features into a single JSON file for ML pipelines.
Includes FULL data (complete vectors, matrices) not just summary statistics.

Usage:
    serializer = FeatureSerializer("sample_001")
    serializer.add_metadata("total_fragments", 8200000)
    serializer.add_fsd(fsd_df)
    serializer.add_wps(wps_df)
    serializer.save(output_dir / "sample_001")  # -> sample_001.features.json
"""

from pathlib import Path
from typing import Any, Dict, Optional, List
import json
from datetime import datetime
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)

# Current schema version
SCHEMA_VERSION = "1.0"


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy types."""
    
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        if pd.isna(obj):
            return None
        return super().default(obj)


class FeatureSerializer:
    """
    Unified JSON serializer for krewlyzer features.
    
    Collects all feature outputs and exports as single JSON file
    with complete data for ML pipelines.
    """
    
    def __init__(self, sample_id: str, version: str = "0.3.2"):
        self.sample_id = sample_id
        self.version = version
        self.metadata: Dict[str, Any] = {}
        self.features: Dict[str, Any] = {}
        self.qc: Dict[str, Any] = {}
    
    # =========================================================================
    # Metadata
    # =========================================================================
    
    def add_metadata(self, key: str, value: Any):
        """Add metadata field."""
        self.metadata[key] = value
    
    def set_metadata(self, metadata: Dict[str, Any]):
        """Set all metadata at once."""
        self.metadata.update(metadata)
    
    # =========================================================================
    # Feature Adders - FULL DATA
    # =========================================================================
    
    def add_fsd(self, df: pd.DataFrame):
        """
        Add FULL FSD matrix (all arms, all size bins).
        
        Includes complete count matrix, not just summary stats.
        """
        if df is None or df.empty:
            return
        
        region_col = "region" if "region" in df.columns else df.columns[0]
        size_bin_cols = [c for c in df.columns if "-" in c and c != region_col]
        
        self.features["fsd"] = {
            "arms": df[region_col].tolist(),
            "size_bins": size_bin_cols,
            "counts": df[size_bin_cols].values.tolist(),
        }
        
        # Add total if present
        if "total" in df.columns:
            self.features["fsd"]["total"] = df["total"].tolist()
        
        # Add log-ratio columns if present
        logr_cols = [c for c in df.columns if "logR" in c or "log_ratio" in c]
        if logr_cols:
            self.features["fsd"]["log_ratios"] = df[logr_cols].values.tolist()
            self.features["fsd"]["log_ratio_cols"] = logr_cols
    
    def add_fsr(self, df: pd.DataFrame):
        """
        Add FULL FSR data (all regions with ratios).
        """
        if df is None or df.empty:
            return
        
        self.features["fsr"] = df.to_dict(orient="records")
    
    def add_fsc(self, df: pd.DataFrame):
        """
        Add FULL FSC data (all windows with coverage).
        """
        if df is None or df.empty:
            return
        
        self.features["fsc"] = df.to_dict(orient="records")
    
    def add_wps(self, df: pd.DataFrame):
        """
        Add FULL WPS vectors (all values per region).
        
        Includes complete wps_nuc and wps_tf vectors, not just means.
        """
        if df is None or df.empty:
            return
        
        # Determine region ID column
        region_col = None
        for col in ["region_id", "group_id", "name"]:
            if col in df.columns:
                region_col = col
                break
        
        if region_col is None:
            logger.warning("WPS DataFrame missing region ID column")
            return
        
        self.features["wps"] = {
            "regions": df[region_col].tolist(),
        }
        
        # Add all WPS-related columns
        wps_cols = ["wps_nuc", "wps_tf", "wps_nuc_smooth", "wps_tf_smooth",
                    "wps_nuc_mean", "wps_tf_mean", "wps_nuc_z", "wps_tf_z",
                    "prot_frac_nuc", "prot_frac_tf"]
        
        for col in wps_cols:
            if col in df.columns:
                self.features["wps"][col] = df[col].tolist()
        
        # Add coordinate info if present
        for col in ["chrom", "center", "start", "end"]:
            if col in df.columns:
                self.features["wps"][col] = df[col].tolist()
    
    def add_motif(
        self,
        edm_df: Optional[pd.DataFrame] = None,
        bpm_df: Optional[pd.DataFrame] = None,
        mds: Optional[float] = None,
        mds_z: Optional[float] = None
    ):
        """
        Add FULL motif frequencies (all 256 k-mers).
        """
        self.features["motif"] = {}
        
        # EDM: End motif frequencies
        if edm_df is not None and not edm_df.empty:
            # Convert to dict of kmer -> frequency
            if len(edm_df) == 1:
                self.features["motif"]["edm"] = edm_df.iloc[0].to_dict()
            else:
                self.features["motif"]["edm"] = edm_df.to_dict(orient="records")
        
        # BPM: Breakpoint motif frequencies  
        if bpm_df is not None and not bpm_df.empty:
            if len(bpm_df) == 1:
                self.features["motif"]["bpm"] = bpm_df.iloc[0].to_dict()
            else:
                self.features["motif"]["bpm"] = bpm_df.to_dict(orient="records")
        
        # MDS: Motif Diversity Score
        if mds is not None:
            self.features["motif"]["mds"] = mds
        if mds_z is not None:
            self.features["motif"]["mds_z"] = mds_z
    
    def add_ocf(self, df: pd.DataFrame):
        """
        Add FULL OCF data (all regions with scores).
        """
        if df is None or df.empty:
            return
        
        self.features["ocf"] = df.to_dict(orient="records")
    
    def add_uxm(self, df: Optional[pd.DataFrame] = None):
        """
        Add UXM (methylation) data if present.
        """
        if df is not None and not df.empty:
            self.features["uxm"] = {
                "enabled": True,
                "data": df.to_dict(orient="records"),
            }
        else:
            self.features["uxm"] = {"enabled": False}
    
    def add_mfsd(self, df: Optional[pd.DataFrame] = None):
        """
        Add mFSD (variant-centric FSD) data if present.
        """
        if df is not None and not df.empty:
            self.features["mfsd"] = {
                "enabled": True,
                "data": df.to_dict(orient="records"),
            }
        else:
            self.features["mfsd"] = {"enabled": False}
    
    # =========================================================================
    # QC Metrics
    # =========================================================================
    
    def add_qc(self, key: str, value: Any):
        """Add QC metric."""
        self.qc[key] = value
    
    def set_qc(self, qc: Dict[str, Any]):
        """Set all QC metrics at once."""
        self.qc.update(qc)
    
    # =========================================================================
    # Serialization
    # =========================================================================
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "schema_version": SCHEMA_VERSION,
            "sample_id": self.sample_id,
            "krewlyzer_version": self.version,
            "timestamp": datetime.now().isoformat(),
            "metadata": self.metadata,
            "features": self.features,
            "qc": self.qc,
        }
    
    def save(self, path: Path) -> Path:
        """
        Save to JSON file.
        
        Args:
            path: Output path (will add .features.json suffix)
        
        Returns:
            Path to saved file
        """
        # Ensure .features.json suffix
        if not str(path).endswith(".features.json"):
            path = path.with_suffix(".features.json")
        
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2, cls=NumpyEncoder)
        
        logger.info(f"Saved unified features JSON: {path}")
        return path
    
    @classmethod
    def from_outputs(
        cls,
        sample_id: str,
        output_dir: Path,
        version: str = "0.3.2"
    ) -> "FeatureSerializer":
        """
        Create FeatureSerializer by reading existing output files.
        
        Args:
            sample_id: Sample identifier
            output_dir: Directory containing krewlyzer outputs
            version: Krewlyzer version
        
        Returns:
            Populated FeatureSerializer
        """
        serializer = cls(sample_id, version)
        
        # Try to load each feature type
        fsd_path = output_dir / f"{sample_id}.FSD.tsv"
        if fsd_path.exists():
            serializer.add_fsd(pd.read_csv(fsd_path, sep="\t"))
        
        fsr_path = output_dir / f"{sample_id}.FSR.tsv"
        if fsr_path.exists():
            serializer.add_fsr(pd.read_csv(fsr_path, sep="\t"))
        
        fsc_path = output_dir / f"{sample_id}.FSC.tsv"
        if fsc_path.exists():
            serializer.add_fsc(pd.read_csv(fsc_path, sep="\t"))
        
        wps_path = output_dir / f"{sample_id}.WPS.parquet"
        if wps_path.exists():
            serializer.add_wps(pd.read_parquet(wps_path))
        
        edm_path = output_dir / f"{sample_id}.EndMotif.tsv"
        bpm_path = output_dir / f"{sample_id}.BreakPointMotif.tsv"
        mds_path = output_dir / f"{sample_id}.MDS.tsv"
        
        edm_df = pd.read_csv(edm_path, sep="\t") if edm_path.exists() else None
        bpm_df = pd.read_csv(bpm_path, sep="\t") if bpm_path.exists() else None
        mds = None
        if mds_path.exists():
            mds_df = pd.read_csv(mds_path, sep="\t")
            # Handle both 'MDS' and 'mds' column names
            mds_col = None
            for col in mds_df.columns:
                if col.lower() == "mds":
                    mds_col = col
                    break
            if mds_col:
                mds = float(mds_df[mds_col].iloc[0])
        serializer.add_motif(edm_df, bpm_df, mds)
        
        ocf_path = output_dir / f"{sample_id}.OCF.tsv"
        if ocf_path.exists():
            serializer.add_ocf(pd.read_csv(ocf_path, sep="\t"))
        
        uxm_path = output_dir / f"{sample_id}.UXM.tsv"
        if uxm_path.exists():
            serializer.add_uxm(pd.read_csv(uxm_path, sep="\t"))
        
        mfsd_path = output_dir / f"{sample_id}.mFSD.tsv"
        if mfsd_path.exists():
            serializer.add_mfsd(pd.read_csv(mfsd_path, sep="\t"))
        
        # Load metadata if present
        meta_path = output_dir / f"{sample_id}.metadata.json"
        if meta_path.exists():
            with open(meta_path) as f:
                meta = json.load(f)
                serializer.set_metadata(meta.get("metadata", meta))
        
        return serializer
