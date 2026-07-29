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
from typing import Any, Dict, Optional
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


def _resolve_output_path(output_dir: Path, sample_id: str, stem: str) -> "Path | None":
    """Locate ``{sample_id}.{stem}`` regardless of the writer's output format.

    krewlyzer emits every tabular feature as ``.tsv``, ``.tsv.gz`` (``--compress``)
    or ``.parquet`` (``--output-format parquet``). Probing only for the bare
    ``.tsv`` name silently skips the feature for any non-default format, which
    produced an EMPTY ``features.json`` for compressed and Parquet runs.

    Args:
        stem: Feature suffix without extension, e.g. ``"FSD"`` or
            ``"FSC.regions.e1only"``.

    Returns:
        The first existing candidate, or ``None`` if the feature was not written.
    """
    from .output_utils import resolve_table_path

    return resolve_table_path(output_dir / f"{sample_id}.{stem}")


class FeatureSerializer:
    """
    Unified JSON serializer for krewlyzer features.

    Collects all feature outputs and exports as single JSON file
    with complete data for ML pipelines.
    """

    def __init__(self, sample_id: str, version: str = "0.8.3"):
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
        wps_cols = [
            "wps_nuc",
            "wps_tf",
            "wps_nuc_smooth",
            "wps_tf_smooth",
            "wps_nuc_mean",
            "wps_tf_mean",
            "wps_nuc_z",
            "wps_tf_z",
            "prot_frac_nuc",
            "prot_frac_tf",
        ]

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
        mds_z: Optional[float] = None,
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
        cls, sample_id: str, output_dir: Path, version: str = "0.8.3"
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

        # Local helper: read_table returns DataFrame|None; callers here always
        # check path.exists() first, so the result is never None. This narrows
        # the type from DataFrame|None to DataFrame so mypy is satisfied.
        from .output_utils import read_table as _read_table
        import pandas as _pd

        def _rt(p: Path) -> "_pd.DataFrame":
            df = _read_table(p)
            assert df is not None, f"read_table returned None for existing file: {p}"
            return df

        # =====================================================================
        # FSD - Fragment Size Distribution
        # =====================================================================
        fsd_path = _resolve_output_path(output_dir, sample_id, "FSD")
        fsd_on_path = _resolve_output_path(output_dir, sample_id, "FSD.ontarget")

        if fsd_path is not None or fsd_on_path is not None:
            fsd_data = {}
            if fsd_path is not None:
                fsd_df = _rt(fsd_path)
                fsd_data["off_target"] = cls._parse_fsd(fsd_df)
            if fsd_on_path is not None:
                fsd_on_df = _rt(fsd_on_path)
                fsd_data["on_target"] = cls._parse_fsd(fsd_on_df)
            serializer.features["fsd"] = fsd_data

        # =====================================================================
        # FSR - Fragment Size Ratio
        # =====================================================================
        fsr_path = _resolve_output_path(output_dir, sample_id, "FSR")
        fsr_on_path = _resolve_output_path(output_dir, sample_id, "FSR.ontarget")

        if fsr_path is not None or fsr_on_path is not None:
            fsr_data = {}
            if fsr_path is not None:
                fsr_data["off_target"] = _rt(fsr_path).to_dict(orient="records")
            if fsr_on_path is not None:
                fsr_data["on_target"] = _rt(fsr_on_path).to_dict(orient="records")
            serializer.features["fsr"] = fsr_data

        # =====================================================================
        # FSC - Fragment Size Coverage
        # =====================================================================
        fsc_path = _resolve_output_path(output_dir, sample_id, "FSC")
        fsc_on_path = _resolve_output_path(output_dir, sample_id, "FSC.ontarget")

        if fsc_path is not None or fsc_on_path is not None:
            fsc_data = {}
            if fsc_path is not None:
                fsc_data["off_target"] = _rt(fsc_path).to_dict(orient="records")
            if fsc_on_path is not None:
                fsc_data["on_target"] = _rt(fsc_on_path).to_dict(orient="records")
            serializer.features["fsc"] = fsc_data

        # =====================================================================
        # FSC Gene - Gene-Centric Fragment Size Coverage (for panel mode)
        # =====================================================================
        fsc_gene_path = _resolve_output_path(output_dir, sample_id, "FSC.gene")
        if fsc_gene_path is not None:
            serializer.features["fsc_gene"] = _rt(fsc_gene_path).to_dict(
                orient="records"
            )

        # =====================================================================
        # FSC Regions - Per-Exon/Target Fragment Size Coverage (for panel mode)
        # =====================================================================
        fsc_region_path = _resolve_output_path(output_dir, sample_id, "FSC.regions")
        if fsc_region_path is not None:
            serializer.features["fsc_region"] = _rt(fsc_region_path).to_dict(
                orient="records"
            )

        # =====================================================================
        # FSC E1-Only - First Exon Per Gene (promoter-proximal sensitivity)
        # Per Helzer et al. (2025): E1 has stronger cancer signal than whole-gene
        # =====================================================================
        fsc_e1_path = _resolve_output_path(output_dir, sample_id, "FSC.regions.e1only")
        if fsc_e1_path is not None:
            serializer.features["fsc_region_e1"] = _rt(fsc_e1_path).to_dict(
                orient="records"
            )
            logger.debug(f"  Loaded fsc_region_e1 from {fsc_e1_path.name}")

        # =====================================================================
        # WPS - Windowed Protection Score
        # =====================================================================
        wps_path = _resolve_output_path(output_dir, sample_id, "WPS")
        wps_panel_path = _resolve_output_path(output_dir, sample_id, "WPS.panel")

        if wps_path is not None:
            serializer.add_wps(pd.read_parquet(wps_path))

        # Panel-specific WPS (for panel mode with --assay)
        if wps_panel_path is not None:
            panel_df = pd.read_parquet(wps_panel_path)
            serializer.features["wps_panel"] = {
                "n_anchors": len(panel_df),
                "data": panel_df.to_dict(orient="records"),
            }

        # WPS Background (Alu element stacking)
        wps_bg_path = _resolve_output_path(output_dir, sample_id, "WPS_background")
        if wps_bg_path is not None:
            bg_df = pd.read_parquet(wps_bg_path)
            serializer.features["wps_background"] = {
                "n_elements": len(bg_df),
                "data": bg_df.to_dict(orient="records"),
            }

        # =====================================================================
        # Motif - EDM, BPM, MDS
        # =====================================================================
        edm_path = _resolve_output_path(output_dir, sample_id, "EndMotif")
        bpm_path = _resolve_output_path(output_dir, sample_id, "BreakPointMotif")
        mds_path = _resolve_output_path(output_dir, sample_id, "MDS")
        edm_on_path = _resolve_output_path(output_dir, sample_id, "EndMotif.ontarget")
        bpm_on_path = _resolve_output_path(
            output_dir, sample_id, "BreakPointMotif.ontarget"
        )
        mds_on_path = _resolve_output_path(output_dir, sample_id, "MDS.ontarget")

        motif_data: Dict[str, Any] = {}

        # Off-target motifs
        if edm_path is not None:
            edm_df = _rt(edm_path)
            motif_data["edm"] = (
                edm_df.iloc[0].to_dict()
                if len(edm_df) == 1
                else edm_df.to_dict(orient="records")
            )
        if bpm_path is not None:
            bpm_df = _rt(bpm_path)
            motif_data["bpm"] = (
                bpm_df.iloc[0].to_dict()
                if len(bpm_df) == 1
                else bpm_df.to_dict(orient="records")
            )
        if mds_path is not None:
            mds_df = _rt(mds_path)
            for col in mds_df.columns:
                if col.lower() == "mds":
                    motif_data["mds"] = float(mds_df[col].iloc[0])  # type: ignore[index]
                    break
            # Also extract mds_z if PON normalization was applied
            if "mds_z" in mds_df.columns:
                val = mds_df["mds_z"].iloc[0]
                if pd.notna(val):
                    motif_data["mds_z"] = float(val)  # type: ignore[index]

        # On-target motifs
        if edm_on_path is not None:
            edm_on_df = _rt(edm_on_path)
            motif_data["edm_on_target"] = (
                edm_on_df.iloc[0].to_dict()
                if len(edm_on_df) == 1
                else edm_on_df.to_dict(orient="records")
            )
        if bpm_on_path is not None:
            bpm_on_df = _rt(bpm_on_path)
            motif_data["bpm_on_target"] = (
                bpm_on_df.iloc[0].to_dict()
                if len(bpm_on_df) == 1
                else bpm_on_df.to_dict(orient="records")
            )
        if mds_on_path is not None:
            mds_on_df = _rt(mds_on_path)
            for col in mds_on_df.columns:
                if col.lower() == "mds":
                    motif_data["mds_on_target"] = float(mds_on_df[col].iloc[0])  # type: ignore[index]
                    break
            # Also extract mds_z for on-target if PON normalization was applied
            if "mds_z" in mds_on_df.columns:
                val = mds_on_df["mds_z"].iloc[0]
                if pd.notna(val):
                    motif_data["mds_z_on_target"] = float(val)  # type: ignore[index]

        if motif_data:
            serializer.features["motif"] = motif_data

        # =====================================================================
        # EndMotif 1-mer (single-base composition at fragment ends)
        # =====================================================================
        edm1_path = _resolve_output_path(output_dir, sample_id, "EndMotif1mer")
        if edm1_path is not None:
            edm1_df = _rt(edm1_path)
            # Convert to dict: {"A": 0.197, "C": 0.345, ...}
            serializer.features.setdefault("motif", {})["edm_1mer"] = {
                row["base"]: row["fraction"] for _, row in edm1_df.iterrows()
            }
            logger.debug(f"  Loaded edm_1mer from {edm1_path.name}")

        # =====================================================================
        # OCF - Orientation-aware cfDNA Fragmentation
        # =====================================================================
        ocf_path = _resolve_output_path(output_dir, sample_id, "OCF")
        ocf_on_path = _resolve_output_path(output_dir, sample_id, "OCF.ontarget")

        if ocf_path is not None or ocf_on_path is not None:
            ocf_data = {}
            if ocf_path is not None:
                ocf_data["off_target"] = _rt(ocf_path).to_dict(orient="records")
            if ocf_on_path is not None:
                ocf_data["on_target"] = _rt(ocf_on_path).to_dict(orient="records")
            serializer.features["ocf"] = ocf_data

        # OCF off-target (panel-specific off-target scores)
        ocf_off_path = _resolve_output_path(output_dir, sample_id, "OCF.offtarget")
        if ocf_off_path is not None:
            serializer.features.setdefault("ocf", {})["offtarget"] = _rt(
                ocf_off_path
            ).to_dict(orient="records")
            logger.debug(f"  Loaded ocf.offtarget from {ocf_off_path.name}")

        # OCF sync files (positional strand-specific profiles)
        for suffix, key in [
            (".OCF.sync.tsv", "sync"),
            (".OCF.offtarget.sync.tsv", "sync_offtarget"),
            (".OCF.ontarget.sync.tsv", "sync_ontarget"),
        ]:
            sync_path = _resolve_output_path(
                output_dir, sample_id, suffix.lstrip(".").removesuffix(".tsv")
            )
            if sync_path is not None:
                serializer.features.setdefault("ocf", {})[key] = _rt(sync_path).to_dict(
                    orient="records"
                )
                logger.debug(f"  Loaded ocf.{key} from {sync_path.name}")

        # =====================================================================
        # UXM - Methylation (optional)
        # =====================================================================
        uxm_path = _resolve_output_path(output_dir, sample_id, "UXM")
        if uxm_path is not None:
            serializer.add_uxm(_rt(uxm_path))

        # =====================================================================
        # mFSD - Mutant Fragment Size Distribution (optional)
        # =====================================================================
        mfsd_path = _resolve_output_path(output_dir, sample_id, "mFSD")
        if mfsd_path is not None:
            mfsd_df = _rt(mfsd_path)
            serializer.features["mfsd"] = {
                "enabled": True,
                "variants": mfsd_df.to_dict(orient="records"),
                "n_variants": len(mfsd_df),
            }

        # =====================================================================
        # TFBS - Transcription Factor Binding Site Region Entropy
        # =====================================================================
        tfbs_path = _resolve_output_path(output_dir, sample_id, "TFBS")
        tfbs_on_path = _resolve_output_path(output_dir, sample_id, "TFBS.ontarget")

        if tfbs_path is not None or tfbs_on_path is not None:
            tfbs_data = {}
            if tfbs_path is not None:
                tfbs_data["off_target"] = _rt(tfbs_path).to_dict(orient="records")
            if tfbs_on_path is not None:
                tfbs_data["on_target"] = _rt(tfbs_on_path).to_dict(orient="records")
            serializer.features["tfbs"] = tfbs_data

        # TFBS sync files (per-TF × per-size fragment distributions)
        for suffix, key in [
            (".TFBS.sync.tsv", "sync"),
            (".TFBS.ontarget.sync.tsv", "sync_ontarget"),
        ]:
            sync_path = _resolve_output_path(
                output_dir, sample_id, suffix.lstrip(".").removesuffix(".tsv")
            )
            if sync_path is not None:
                serializer.features.setdefault("tfbs", {})[key] = _rt(
                    sync_path
                ).to_dict(orient="records")
                logger.debug(f"  Loaded tfbs.{key} from {sync_path.name}")

        # =====================================================================
        # ATAC - ATAC-seq Region Entropy
        # =====================================================================
        atac_path = _resolve_output_path(output_dir, sample_id, "ATAC")
        atac_on_path = _resolve_output_path(output_dir, sample_id, "ATAC.ontarget")

        if atac_path is not None or atac_on_path is not None:
            atac_data = {}
            if atac_path is not None:
                atac_data["off_target"] = _rt(atac_path).to_dict(orient="records")
            if atac_on_path is not None:
                atac_data["on_target"] = _rt(atac_on_path).to_dict(orient="records")
            serializer.features["atac"] = atac_data

        # ATAC sync files (per-tissue × per-size fragment distributions)
        for suffix, key in [
            (".ATAC.sync.tsv", "sync"),
            (".ATAC.ontarget.sync.tsv", "sync_ontarget"),
        ]:
            sync_path = _resolve_output_path(
                output_dir, sample_id, suffix.lstrip(".").removesuffix(".tsv")
            )
            if sync_path is not None:
                serializer.features.setdefault("atac", {})[key] = _rt(
                    sync_path
                ).to_dict(orient="records")
                logger.debug(f"  Loaded atac.{key} from {sync_path.name}")

        # =====================================================================
        # GC Correction Factors
        # =====================================================================
        gc_factors_path = _resolve_output_path(
            output_dir, sample_id, "correction_factors"
        )
        gc_factors_on_path = _resolve_output_path(
            output_dir, sample_id, "correction_factors.ontarget"
        )

        if gc_factors_path is not None or gc_factors_on_path is not None:
            gc_data = {}
            if gc_factors_path is not None:
                gc_data["off_target"] = _rt(gc_factors_path).to_dict(orient="records")
            if gc_factors_on_path is not None:
                gc_data["on_target"] = _rt(gc_factors_on_path).to_dict(orient="records")
            serializer.features["gc_factors"] = gc_data

        # =====================================================================
        # FSC Counts (raw per-bin size-class counts with GC)
        # =====================================================================
        fsc_counts_path = _resolve_output_path(output_dir, sample_id, "fsc_counts")
        if fsc_counts_path is not None:
            serializer.features["fsc_counts"] = _rt(fsc_counts_path).to_dict(
                orient="records"
            )
            logger.debug(f"  Loaded fsc_counts from {fsc_counts_path.name}")

        # =====================================================================
        # Region MDS - Per-Exon/Target Motif Diversity Score (Helzer et al.)
        # =====================================================================
        mds_exon_path = _resolve_output_path(output_dir, sample_id, "MDS.exon")
        mds_gene_path = _resolve_output_path(output_dir, sample_id, "MDS.gene")

        if mds_exon_path is not None or mds_gene_path is not None:
            region_mds_data: Dict[str, Any] = {}

            if mds_exon_path is not None:
                exon_df = _rt(mds_exon_path)
                region_mds_data["exon"] = exon_df.to_dict(orient="records")
                region_mds_data["n_exons"] = len(exon_df)

                # Summary statistics
                if "mds" in exon_df.columns:
                    region_mds_data["mds_exon_mean"] = float(exon_df["mds"].mean())  # type: ignore[index]
                    region_mds_data["mds_exon_std"] = float(exon_df["mds"].std())  # type: ignore[index]

            if mds_gene_path is not None:
                gene_df = _rt(mds_gene_path)
                region_mds_data["gene"] = gene_df.to_dict(orient="records")
                region_mds_data["n_genes"] = len(gene_df)  # type: ignore[index]

                # E1 summary
                if "mds_e1" in gene_df.columns:
                    region_mds_data["mds_e1_mean"] = float(gene_df["mds_e1"].mean())  # type: ignore[index]

            serializer.features["region_mds"] = region_mds_data

        # =====================================================================
        # Metadata — read from .metadata.tsv (Parquet-first, TSV fallback).
        # .metadata.json was removed; metadata is now consistent TSV/Parquet.
        # =====================================================================
        meta_tsv_path = _resolve_output_path(output_dir, sample_id, "metadata")
        if meta_tsv_path is not None:
            meta_df = _rt(meta_tsv_path)
            if not meta_df.empty:
                serializer.set_metadata(meta_df.iloc[0].to_dict())  # type: ignore[arg-type]

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
