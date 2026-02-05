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
    
    def __init__(self, sample_id: str, version: str = "0.5.2"):
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
    
    def add_fsc_e1(self, df: pd.DataFrame):
        """
        Add FULL FSC E1-only data (first exon per gene).
        
        E1 (promoter-proximal) has stronger cancer signal per Helzer et al. (2025).
        """
        if df is None or df.empty:
            return
        
        self.features["fsc_region_e1"] = df.to_dict(orient="records")
    
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
        version: str = "0.5.2"
    ) -> "FeatureSerializer":
        """
        Create FeatureSerializer by reading existing output files.
        
        Args:
            sample_id: Sample identifier
            output_dir: Directory containing krewlyzer outputs
            version: Krewlyzer version
        
        Returns:
            Populated FeatureSerializer with both off-target and on-target data
        """
        serializer = cls(sample_id, version)
        output_dir = Path(output_dir)
        
        # =====================================================================
        # FSD - Fragment Size Distribution
        # =====================================================================
        fsd_path = output_dir / f"{sample_id}.FSD.tsv"
        fsd_on_path = output_dir / f"{sample_id}.FSD.ontarget.tsv"
        
        if fsd_path.exists() or fsd_on_path.exists():
            fsd_data = {}
            if fsd_path.exists():
                fsd_df = pd.read_csv(fsd_path, sep="\t")
                fsd_data["off_target"] = cls._parse_fsd(fsd_df)
            if fsd_on_path.exists():
                fsd_on_df = pd.read_csv(fsd_on_path, sep="\t")
                fsd_data["on_target"] = cls._parse_fsd(fsd_on_df)
            serializer.features["fsd"] = fsd_data
        
        # =====================================================================
        # FSR - Fragment Size Ratio
        # =====================================================================
        fsr_path = output_dir / f"{sample_id}.FSR.tsv"
        fsr_on_path = output_dir / f"{sample_id}.FSR.ontarget.tsv"
        
        if fsr_path.exists() or fsr_on_path.exists():
            fsr_data = {}
            if fsr_path.exists():
                fsr_data["off_target"] = pd.read_csv(fsr_path, sep="\t").to_dict(orient="records")
            if fsr_on_path.exists():
                fsr_data["on_target"] = pd.read_csv(fsr_on_path, sep="\t").to_dict(orient="records")
            serializer.features["fsr"] = fsr_data
        
        # =====================================================================
        # FSC - Fragment Size Coverage
        # =====================================================================
        fsc_path = output_dir / f"{sample_id}.FSC.tsv"
        fsc_on_path = output_dir / f"{sample_id}.FSC.ontarget.tsv"
        
        if fsc_path.exists() or fsc_on_path.exists():
            fsc_data = {}
            if fsc_path.exists():
                fsc_data["off_target"] = pd.read_csv(fsc_path, sep="\t").to_dict(orient="records")
            if fsc_on_path.exists():
                fsc_data["on_target"] = pd.read_csv(fsc_on_path, sep="\t").to_dict(orient="records")
            serializer.features["fsc"] = fsc_data
        
        # =====================================================================
        # FSC Gene - Gene-Centric Fragment Size Coverage (for panel mode)
        # =====================================================================
        fsc_gene_path = output_dir / f"{sample_id}.FSC.gene.tsv"
        if fsc_gene_path.exists():
            serializer.features["fsc_gene"] = pd.read_csv(fsc_gene_path, sep="\t").to_dict(orient="records")
        
        # =====================================================================
        # FSC Regions - Per-Exon/Target Fragment Size Coverage (for panel mode)
        # =====================================================================
        fsc_region_path = output_dir / f"{sample_id}.FSC.regions.tsv"
        if fsc_region_path.exists():
            serializer.features["fsc_region"] = pd.read_csv(fsc_region_path, sep="\t").to_dict(orient="records")
        
        # =====================================================================
        # FSC E1-Only - First Exon Per Gene (promoter-proximal sensitivity)
        # Per Helzer et al. (2025): E1 has stronger cancer signal than whole-gene
        # =====================================================================
        fsc_e1_path = output_dir / f"{sample_id}.FSC.regions.e1only.tsv"
        if fsc_e1_path.exists():
            serializer.features["fsc_region_e1"] = pd.read_csv(fsc_e1_path, sep="\t").to_dict(orient="records")
            logger.debug(f"  Loaded fsc_region_e1 from {fsc_e1_path.name}")
        
        # =====================================================================
        # WPS - Windowed Protection Score
        # =====================================================================
        wps_path = output_dir / f"{sample_id}.WPS.parquet"
        wps_panel_path = output_dir / f"{sample_id}.WPS.panel.parquet"
        
        if wps_path.exists():
            serializer.add_wps(pd.read_parquet(wps_path))
        
        # Panel-specific WPS (for panel mode with --assay)
        if wps_panel_path.exists():
            panel_df = pd.read_parquet(wps_panel_path)
            serializer.features["wps_panel"] = {
                "n_anchors": len(panel_df),
                "data": panel_df.to_dict(orient="records")
            }
        
        # WPS Background (Alu element stacking)
        wps_bg_path = output_dir / f"{sample_id}.WPS_background.parquet"
        if wps_bg_path.exists():
            bg_df = pd.read_parquet(wps_bg_path)
            serializer.features["wps_background"] = {
                "n_elements": len(bg_df),
                "data": bg_df.to_dict(orient="records")
            }
        
        # =====================================================================
        # Motif - EDM, BPM, MDS
        # =====================================================================
        edm_path = output_dir / f"{sample_id}.EndMotif.tsv"
        bpm_path = output_dir / f"{sample_id}.BreakPointMotif.tsv"
        mds_path = output_dir / f"{sample_id}.MDS.tsv"
        edm_on_path = output_dir / f"{sample_id}.EndMotif.ontarget.tsv"
        bpm_on_path = output_dir / f"{sample_id}.BreakPointMotif.ontarget.tsv"
        mds_on_path = output_dir / f"{sample_id}.MDS.ontarget.tsv"
        
        motif_data = {}
        
        # Off-target motifs
        if edm_path.exists():
            edm_df = pd.read_csv(edm_path, sep="\t")
            motif_data["edm"] = edm_df.iloc[0].to_dict() if len(edm_df) == 1 else edm_df.to_dict(orient="records")
        if bpm_path.exists():
            bpm_df = pd.read_csv(bpm_path, sep="\t")
            motif_data["bpm"] = bpm_df.iloc[0].to_dict() if len(bpm_df) == 1 else bpm_df.to_dict(orient="records")
        if mds_path.exists():
            mds_df = pd.read_csv(mds_path, sep="\t")
            for col in mds_df.columns:
                if col.lower() == "mds":
                    motif_data["mds"] = float(mds_df[col].iloc[0])
                    break
        
        # On-target motifs
        if edm_on_path.exists():
            edm_on_df = pd.read_csv(edm_on_path, sep="\t")
            motif_data["edm_on_target"] = edm_on_df.iloc[0].to_dict() if len(edm_on_df) == 1 else edm_on_df.to_dict(orient="records")
        if bpm_on_path.exists():
            bpm_on_df = pd.read_csv(bpm_on_path, sep="\t")
            motif_data["bpm_on_target"] = bpm_on_df.iloc[0].to_dict() if len(bpm_on_df) == 1 else bpm_on_df.to_dict(orient="records")
        if mds_on_path.exists():
            mds_on_df = pd.read_csv(mds_on_path, sep="\t")
            for col in mds_on_df.columns:
                if col.lower() == "mds":
                    motif_data["mds_on_target"] = float(mds_on_df[col].iloc[0])
                    break
        
        if motif_data:
            serializer.features["motif"] = motif_data
        
        # =====================================================================
        # OCF - Orientation-aware cfDNA Fragmentation
        # =====================================================================
        ocf_path = output_dir / f"{sample_id}.OCF.tsv"
        ocf_on_path = output_dir / f"{sample_id}.OCF.ontarget.tsv"
        
        if ocf_path.exists() or ocf_on_path.exists():
            ocf_data = {}
            if ocf_path.exists():
                ocf_data["off_target"] = pd.read_csv(ocf_path, sep="\t").to_dict(orient="records")
            if ocf_on_path.exists():
                ocf_data["on_target"] = pd.read_csv(ocf_on_path, sep="\t").to_dict(orient="records")
            serializer.features["ocf"] = ocf_data
        
        # =====================================================================
        # UXM - Methylation (optional)
        # =====================================================================
        uxm_path = output_dir / f"{sample_id}.UXM.tsv"
        if uxm_path.exists():
            serializer.add_uxm(pd.read_csv(uxm_path, sep="\t"))
        
        # =====================================================================
        # mFSD - Mutant Fragment Size Distribution (optional)
        # =====================================================================
        mfsd_path = output_dir / f"{sample_id}.mFSD.tsv"
        if mfsd_path.exists():
            mfsd_df = pd.read_csv(mfsd_path, sep="\t")
            serializer.features["mfsd"] = {
                "enabled": True,
                "variants": mfsd_df.to_dict(orient="records"),
                "n_variants": len(mfsd_df),
            }
        
        # =====================================================================
        # TFBS - Transcription Factor Binding Site Region Entropy
        # =====================================================================
        tfbs_path = output_dir / f"{sample_id}.TFBS.tsv"
        tfbs_on_path = output_dir / f"{sample_id}.TFBS.ontarget.tsv"
        
        if tfbs_path.exists() or tfbs_on_path.exists():
            tfbs_data = {}
            if tfbs_path.exists():
                tfbs_data["off_target"] = pd.read_csv(tfbs_path, sep="\t").to_dict(orient="records")
            if tfbs_on_path.exists():
                tfbs_data["on_target"] = pd.read_csv(tfbs_on_path, sep="\t").to_dict(orient="records")
            serializer.features["tfbs"] = tfbs_data
        
        # =====================================================================
        # ATAC - ATAC-seq Region Entropy
        # =====================================================================
        atac_path = output_dir / f"{sample_id}.ATAC.tsv"
        atac_on_path = output_dir / f"{sample_id}.ATAC.ontarget.tsv"
        
        if atac_path.exists() or atac_on_path.exists():
            atac_data = {}
            if atac_path.exists():
                atac_data["off_target"] = pd.read_csv(atac_path, sep="\t").to_dict(orient="records")
            if atac_on_path.exists():
                atac_data["on_target"] = pd.read_csv(atac_on_path, sep="\t").to_dict(orient="records")
            serializer.features["atac"] = atac_data
        
        # =====================================================================
        # GC Correction Factors
        # =====================================================================
        gc_factors_path = output_dir / f"{sample_id}.correction_factors.tsv"
        gc_factors_on_path = output_dir / f"{sample_id}.correction_factors.ontarget.tsv"
        
        if gc_factors_path.exists() or gc_factors_on_path.exists():
            gc_data = {}
            if gc_factors_path.exists():
                gc_data["off_target"] = pd.read_csv(gc_factors_path, sep="\t").to_dict(orient="records")
            if gc_factors_on_path.exists():
                gc_data["on_target"] = pd.read_csv(gc_factors_on_path, sep="\t").to_dict(orient="records")
            serializer.features["gc_factors"] = gc_data
        
        # =====================================================================
        # Region MDS - Per-Exon/Target Motif Diversity Score (Helzer et al.)
        # =====================================================================
        mds_exon_path = output_dir / f"{sample_id}.MDS.exon.tsv"
        mds_gene_path = output_dir / f"{sample_id}.MDS.gene.tsv"
        
        if mds_exon_path.exists() or mds_gene_path.exists():
            region_mds_data = {}
            
            if mds_exon_path.exists():
                exon_df = pd.read_csv(mds_exon_path, sep="\t")
                region_mds_data["exon"] = exon_df.to_dict(orient="records")
                region_mds_data["n_exons"] = len(exon_df)
                
                # Summary statistics
                if "mds" in exon_df.columns:
                    region_mds_data["mds_exon_mean"] = float(exon_df["mds"].mean())
                    region_mds_data["mds_exon_std"] = float(exon_df["mds"].std())
            
            if mds_gene_path.exists():
                gene_df = pd.read_csv(mds_gene_path, sep="\t")
                region_mds_data["gene"] = gene_df.to_dict(orient="records")
                region_mds_data["n_genes"] = len(gene_df)
                
                # E1 summary
                if "mds_e1" in gene_df.columns:
                    region_mds_data["mds_e1_mean"] = float(gene_df["mds_e1"].mean())
            
            serializer.features["region_mds"] = region_mds_data
        
        # =====================================================================
        # Metadata
        # =====================================================================
        meta_path = output_dir / f"{sample_id}.metadata.json"
        if meta_path.exists():
            with open(meta_path) as f:
                meta = json.load(f)
                serializer.set_metadata(meta.get("metadata", meta))
        
        return serializer
    
    @staticmethod
    def _parse_fsd(df: pd.DataFrame) -> dict:
        """Parse FSD DataFrame into dictionary."""
        region_col = "region" if "region" in df.columns else df.columns[0]
        size_bin_cols = [c for c in df.columns if "-" in c and c != region_col]
        
        result = {
            "arms": df[region_col].tolist(),
            "size_bins": size_bin_cols,
            "counts": df[size_bin_cols].values.tolist(),
        }
        
        if "total" in df.columns:
            result["total"] = df["total"].tolist()
        
        logr_cols = [c for c in df.columns if "logR" in c or "log_ratio" in c]
        if logr_cols:
            result["log_ratios"] = df[logr_cols].values.tolist()
            result["log_ratio_cols"] = logr_cols
        
        return result

