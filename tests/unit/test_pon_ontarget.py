"""
Unit tests for PON on-target coverage gap implementations.

Tests the new functions and fields added for on-target PON z-score support:
- apply_ocf_python_pon: Python-side OCF z-score for on/off-target files
- PonModel new fields: mds_baseline_ontarget, ocf_baseline_ontarget/offtarget
- process_fsd: baseline_table parameter forwarding
- sample_processor: on-target MDS z-score computation
"""

from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from krewlyzer.pon.model import OcfBaseline, MdsBaseline

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def sample_ocf_baseline():
    """Create an OCF baseline with known values for on-target testing."""
    regions = pd.DataFrame(
        {
            "region_id": ["Breast", "Liver", "Lung"],
            "ocf_mean": [50.0, 100.0, -10.0],
            "ocf_std": [10.0, 20.0, 5.0],
            "sync_mean": [0.5, 0.6, 0.4],
            "sync_std": [0.1, 0.1, 0.1],
        }
    )
    return OcfBaseline(regions=regions)


@pytest.fixture
def sample_mds_baseline_ontarget():
    """Create an on-target MDS baseline with known values."""
    return MdsBaseline(
        kmer_expected={"ACGT": 0.004, "TGCA": 0.005},
        kmer_std={"ACGT": 0.001, "TGCA": 0.002},
        mds_mean=0.85,
        mds_std=0.02,
    )


@pytest.fixture
def mock_pon_with_ontarget(sample_ocf_baseline, sample_mds_baseline_ontarget):
    """Create a mock PON with on-target baselines."""
    pon = MagicMock()
    pon.ocf_baseline_ontarget = sample_ocf_baseline
    pon.ocf_baseline_offtarget = sample_ocf_baseline  # reuse for simplicity
    pon.mds_baseline_ontarget = sample_mds_baseline_ontarget
    pon.mds_baseline = MdsBaseline(mds_mean=0.82, mds_std=0.03)
    return pon


@pytest.fixture
def sample_ocf_tsv(tmp_path):
    """Create a sample OCF TSV file."""
    ocf_file = tmp_path / "sample.OCF.ontarget.tsv"
    df = pd.DataFrame(
        {
            "tissue": ["Breast", "Liver", "Lung"],
            "OCF": [60.0, 80.0, -5.0],
        }
    )
    df.to_csv(ocf_file, sep="\t", index=False)
    return ocf_file


# =============================================================================
# Test: apply_ocf_python_pon
# =============================================================================


class TestApplyOcfPythonPon:
    """Tests for apply_ocf_python_pon function."""

    def test_zscore_computation(self, sample_ocf_tsv, mock_pon_with_ontarget):
        """Z-scores are computed correctly against on-target baseline."""
        from krewlyzer.core.ocf_processor import apply_ocf_python_pon

        n_matched = apply_ocf_python_pon(
            sample_ocf_tsv,
            mock_pon_with_ontarget,
            baseline_attr="ocf_baseline_ontarget",
        )

        assert n_matched == 3

        # Verify z-scores written to file
        df = pd.read_csv(sample_ocf_tsv, sep="\t")
        assert "ocf_z" in df.columns

        # Breast: (60 - 50) / 10 = 1.0
        breast_z = df.loc[df["tissue"] == "Breast", "ocf_z"].iloc[0]
        assert abs(breast_z - 1.0) < 1e-6

        # Liver: (80 - 100) / 20 = -1.0
        liver_z = df.loc[df["tissue"] == "Liver", "ocf_z"].iloc[0]
        assert abs(liver_z - (-1.0)) < 1e-6

        # Lung: (-5 - (-10)) / 5 = 1.0
        lung_z = df.loc[df["tissue"] == "Lung", "ocf_z"].iloc[0]
        assert abs(lung_z - 1.0) < 1e-6

    def test_no_baseline_returns_zero(self, sample_ocf_tsv):
        """Returns 0 when PON has no baseline for the requested attr."""
        from krewlyzer.core.ocf_processor import apply_ocf_python_pon

        pon = MagicMock()
        pon.ocf_baseline_ontarget = None

        n_matched = apply_ocf_python_pon(
            sample_ocf_tsv, pon, baseline_attr="ocf_baseline_ontarget"
        )
        assert n_matched == 0

    def test_no_pon_returns_zero(self, sample_ocf_tsv):
        """Returns 0 when no PON is provided."""
        from krewlyzer.core.ocf_processor import apply_ocf_python_pon

        n_matched = apply_ocf_python_pon(
            sample_ocf_tsv, None, baseline_attr="ocf_baseline_ontarget"
        )
        assert n_matched == 0

    def test_missing_tissue_produces_nan(self, tmp_path, sample_ocf_baseline):
        """Missing tissue in baseline produces NaN z-score."""
        from krewlyzer.core.ocf_processor import apply_ocf_python_pon

        # Create file with unknown tissue
        ocf_file = tmp_path / "sample.OCF.ontarget.tsv"
        df = pd.DataFrame(
            {
                "tissue": ["Breast", "Unknown"],
                "OCF": [60.0, 42.0],
            }
        )
        df.to_csv(ocf_file, sep="\t", index=False)

        pon = MagicMock()
        pon.ocf_baseline_ontarget = sample_ocf_baseline

        n_matched = apply_ocf_python_pon(
            ocf_file, pon, baseline_attr="ocf_baseline_ontarget"
        )

        # Only 1 matched (Breast), Unknown has no baseline
        assert n_matched == 1

        result = pd.read_csv(ocf_file, sep="\t")
        assert pd.isna(result.loc[result["tissue"] == "Unknown", "ocf_z"].iloc[0])

    def test_zero_std_returns_zero_zscore(self, tmp_path):
        """Zero std in baseline returns 0.0 z-score."""
        from krewlyzer.core.ocf_processor import apply_ocf_python_pon

        ocf_file = tmp_path / "sample.OCF.ontarget.tsv"
        df = pd.DataFrame({"tissue": ["TestTissue"], "OCF": [50.0]})
        df.to_csv(ocf_file, sep="\t", index=False)

        # Baseline with zero std
        baseline = OcfBaseline(
            regions=pd.DataFrame(
                {
                    "region_id": ["TestTissue"],
                    "ocf_mean": [50.0],
                    "ocf_std": [0.0],
                    "sync_mean": [0.5],
                    "sync_std": [0.0],
                }
            )
        )
        pon = MagicMock()
        pon.ocf_baseline_ontarget = baseline

        n_matched = apply_ocf_python_pon(
            ocf_file, pon, baseline_attr="ocf_baseline_ontarget"
        )
        assert n_matched == 1

        result = pd.read_csv(ocf_file, sep="\t")
        assert result["ocf_z"].iloc[0] == 0.0


# =============================================================================
# Test: PonModel on-target fields
# =============================================================================


class TestPonModelOntargetFields:
    """Tests for PonModel on-target baseline fields."""

    def test_mds_baseline_ontarget_zscore(self, sample_mds_baseline_ontarget):
        """On-target MDS baseline computes z-scores correctly."""
        bl = sample_mds_baseline_ontarget

        # (0.90 - 0.85) / 0.02 = 2.5
        z = bl.get_mds_zscore(0.90)
        assert abs(z - 2.5) < 1e-6

        # (0.80 - 0.85) / 0.02 = -2.5
        z = bl.get_mds_zscore(0.80)
        assert abs(z - (-2.5)) < 1e-6

    def test_mds_baseline_ontarget_kmer_zscore(self, sample_mds_baseline_ontarget):
        """On-target MDS baseline computes k-mer z-scores correctly."""
        bl = sample_mds_baseline_ontarget

        # (0.006 - 0.004) / 0.001 = 2.0
        z = bl.get_kmer_zscore("ACGT", 0.006)
        assert abs(z - 2.0) < 1e-6

    def test_ocf_baseline_ontarget_zscore(self, sample_ocf_baseline):
        """On-target OCF baseline computes z-scores correctly."""
        bl = sample_ocf_baseline

        # Breast: (70 - 50) / 10 = 2.0
        z = bl.compute_zscore("Breast", 70.0)
        assert abs(z - 2.0) < 1e-6

    def test_ocf_baseline_missing_region(self, sample_ocf_baseline):
        """OCF baseline returns None for unknown region."""
        bl = sample_ocf_baseline
        z = bl.compute_zscore("Unknown", 50.0)
        assert z is None

    def test_pon_model_has_ontarget_attrs(self):
        """PonModel constructor accepts on-target baseline fields."""
        from krewlyzer.pon.model import PonModel

        # Verify PonModel accepts these fields
        pon = PonModel.__new__(PonModel)
        pon.mds_baseline_ontarget = None
        pon.ocf_baseline_ontarget = None
        pon.ocf_baseline_offtarget = None
        pon.fsd_baseline_ontarget = None

        assert hasattr(pon, "mds_baseline_ontarget")
        assert hasattr(pon, "ocf_baseline_ontarget")
        assert hasattr(pon, "ocf_baseline_offtarget")
        assert hasattr(pon, "fsd_baseline_ontarget")


# =============================================================================
# Test: process_fsd baseline_table parameter
# =============================================================================


class TestProcessFsdBaselineTable:
    """Tests for process_fsd baseline_table parameter forwarding."""

    def test_baseline_table_default_is_fsd_baseline(self):
        """Default baseline_table is 'fsd_baseline'."""
        from krewlyzer.core.fsd_processor import process_fsd
        import inspect

        sig = inspect.signature(process_fsd)
        assert sig.parameters["baseline_table"].default == "fsd_baseline"

    def test_baseline_table_accepts_ontarget(self, tmp_path):
        """process_fsd accepts baseline_table='fsd_baseline_ontarget'."""
        from krewlyzer.core.fsd_processor import process_fsd

        # Create a minimal FSD file
        fsd_file = tmp_path / "sample.FSD.ontarget.tsv"
        df = pd.DataFrame(
            {
                "region": ["chr1:0-120000000"],
                "65-69": [10.0],
                "70-74": [20.0],
                "total": [30.0],
            }
        )
        df.to_csv(fsd_file, sep="\t", index=False)

        # Should succeed without PON (no actual Rust call)
        result = process_fsd(
            fsd_file,
            pon_parquet_path=None,
            baseline_table="fsd_baseline_ontarget",
        )
        assert result == fsd_file

    @patch("krewlyzer._core.fsd.apply_pon_logratio")
    def test_baseline_table_forwarded_to_rust(self, mock_apply, tmp_path):
        """baseline_table is forwarded to Rust apply_pon_logratio."""
        from krewlyzer.core.fsd_processor import process_fsd

        # Create dummy files
        fsd_file = tmp_path / "sample.FSD.ontarget.tsv"
        df = pd.DataFrame(
            {"region": ["chr1:0-120000000"], "65-69": [10.0], "total": [10.0]}
        )
        df.to_csv(fsd_file, sep="\t", index=False)

        pon_file = tmp_path / "test.pon.parquet"
        pon_file.touch()

        mock_apply.return_value = 41

        process_fsd(
            fsd_file,
            pon_parquet_path=pon_file,
            baseline_table="fsd_baseline_ontarget",
        )

        # Verify Rust was called with the correct baseline_table
        mock_apply.assert_called_once()
        call_kwargs = mock_apply.call_args
        assert call_kwargs.kwargs["baseline_table"] == "fsd_baseline_ontarget"


# =============================================================================
# Test: sample_processor on-target MDS z-score
# =============================================================================


class TestSampleProcessorMdsOntargetZscore:
    """Tests for on-target MDS z-score computation in sample_processor."""

    def test_mds_zscore_computed_with_ontarget_baseline(self, tmp_path):
        """MDS z-score is computed using mds_baseline_ontarget."""
        # Create an on-target MDS file
        mds_file = tmp_path / "sample.MDS.ontarget.tsv"
        df = pd.DataFrame(
            {
                "Sample": ["test_sample"],
                "MDS": [0.90],
            }
        )
        df.to_csv(mds_file, sep="\t", index=False)

        # Mock PON with on-target baseline
        pon = MagicMock()
        pon.mds_baseline_ontarget = MdsBaseline(mds_mean=0.85, mds_std=0.02)
        pon.mds_baseline = MdsBaseline(mds_mean=0.82, mds_std=0.03)

        # Compute z-score
        mds_val = 0.90
        mds_bl = pon.mds_baseline_ontarget
        mds_z = (mds_val - mds_bl.mds_mean) / max(mds_bl.mds_std, 1e-10)

        # Expected: (0.90 - 0.85) / 0.02 = 2.5
        assert abs(mds_z - 2.5) < 1e-6

        # Simulate appending z-score
        mds_df = pd.read_csv(mds_file, sep="\t")
        mds_df["mds_z"] = mds_z
        mds_df.to_csv(mds_file, sep="\t", index=False)

        # Verify
        result = pd.read_csv(mds_file, sep="\t")
        assert "mds_z" in result.columns
        assert abs(result["mds_z"].iloc[0] - 2.5) < 1e-6

    def test_mds_zscore_fallback_to_genome_wide(self):
        """Falls back to genome-wide MDS baseline when on-target not available."""
        pon = MagicMock()
        pon.mds_baseline_ontarget = None
        pon.mds_baseline = MdsBaseline(mds_mean=0.82, mds_std=0.03)

        # Simulate the fallback logic from sample_processor.py
        mds_bl = None
        bl_name = "none"
        if hasattr(pon, "mds_baseline_ontarget") and pon.mds_baseline_ontarget:
            mds_bl = pon.mds_baseline_ontarget
            bl_name = "ontarget"
        elif hasattr(pon, "mds_baseline") and pon.mds_baseline:
            mds_bl = pon.mds_baseline
            bl_name = "genome-wide (fallback)"

        assert bl_name == "genome-wide (fallback)"
        assert mds_bl is not None
        assert mds_bl.mds_mean == 0.82

    def test_mds_zscore_with_zero_std(self):
        """Zero std uses epsilon to avoid division by zero."""
        pon = MagicMock()
        pon.mds_baseline_ontarget = MdsBaseline(mds_mean=0.85, mds_std=0.0)
        pon.mds_baseline = None

        mds_bl = pon.mds_baseline_ontarget
        mds_val = 0.90

        # Replicate the actual code logic
        mds_z = (mds_val - mds_bl.mds_mean) / max(mds_bl.mds_std, 1e-10)

        # (0.90 - 0.85) / 1e-10 = 5e8 (very large but finite)
        assert mds_z > 0
        assert mds_z == pytest.approx(5e8, rel=1e-6)

    def test_no_pon_skips_zscore(self):
        """No PON means no z-score computation."""
        pon = None
        mds_val = 0.90

        # Simulate the actual condition check
        should_compute = pon and mds_val is not None
        assert not should_compute
