"""
Unit tests for PON model classes.

Tests the new baseline classes and methods added in the PON enhancement.
"""

import pytest
import numpy as np
import pandas as pd
from krewlyzer.pon.model import (
    PonModel,
    GcBiasModel,
    FsdBaseline,
    WpsBaseline,
    OcfBaseline,
    MdsBaseline,
    WpsBackgroundBaseline,
    RegionMdsBaseline,
    FscGeneBaseline,
    FscRegionBaseline,
)


class TestFsdBaseline:
    """Tests for FsdBaseline class."""

    @pytest.fixture
    def sample_fsd_baseline(self):
        return FsdBaseline(
            size_bins=[65, 70, 75, 80, 85],
            arms={
                "1p": {
                    "expected": [0.1, 0.15, 0.2, 0.15, 0.1],
                    "std": [0.01, 0.02, 0.03, 0.02, 0.01],
                },
                "1q": {
                    "expected": [0.12, 0.18, 0.22, 0.18, 0.12],
                    "std": [0.02, 0.03, 0.04, 0.03, 0.02],
                },
            },
        )

    def test_get_expected(self, sample_fsd_baseline):
        """Test get_expected returns correct value."""
        expected = sample_fsd_baseline.get_expected("1p", 75)
        assert expected == 0.2

    def test_get_std(self, sample_fsd_baseline):
        """Test get_std returns correct value."""
        std = sample_fsd_baseline.get_std("1p", 75)
        assert std == 0.03

    def test_get_stats(self, sample_fsd_baseline):
        """Test get_stats returns (expected, std) tuple."""
        stats = sample_fsd_baseline.get_stats("1p", 75)
        assert stats == (0.2, 0.03)

    def test_get_stats_missing_arm(self, sample_fsd_baseline):
        """Test get_stats returns None for missing arm."""
        stats = sample_fsd_baseline.get_stats("unknown", 75)
        assert stats is None


class TestWpsBaseline:
    """Tests for WpsBaseline class."""

    @pytest.fixture
    def sample_wps_v1(self):
        df = pd.DataFrame(
            {
                "region_id": ["gene1", "gene2"],
                "wps_long_mean": [100.0, 120.0],
                "wps_long_std": [10.0, 12.0],
                "wps_short_mean": [50.0, 55.0],
                "wps_short_std": [5.0, 6.0],
            }
        )
        return WpsBaseline(regions=df, schema_version="1.0")

    def test_get_baseline(self, sample_wps_v1):
        """Test get_baseline returns dict for region."""
        baseline = sample_wps_v1.get_baseline("gene1")
        assert baseline is not None
        assert baseline["wps_long_mean"] == 100.0

    def test_get_stats(self, sample_wps_v1):
        """Test get_stats returns (mean, std) tuple."""
        stats = sample_wps_v1.get_stats("gene1", "wps_long")
        assert stats == (100.0, 10.0)

    def test_v1_vector_returns_none(self, sample_wps_v1):
        """Test v1.0 schema returns None for vector methods."""
        assert sample_wps_v1.get_baseline_vector("gene1") is None


class TestWpsBaselineV2:
    """Tests for WpsBaseline v2.0 vector format."""

    @pytest.fixture
    def sample_wps_v2(self):
        """Create v2.0 baseline with 200-element vectors."""
        # Create simple 200-element vectors for testing
        mean_vec = [float(i) for i in range(200)]  # 0, 1, 2, ..., 199
        std_vec = [0.1] * 200  # Constant std for easy z-score calculation

        df = pd.DataFrame(
            {
                "region_id": ["gene1", "gene2"],
                "wps_nuc_mean": [mean_vec, [x + 10 for x in mean_vec]],
                "wps_nuc_std": [std_vec, std_vec],
                "wps_tf_mean": [
                    [x * 0.5 for x in mean_vec],
                    [x * 0.5 for x in mean_vec],
                ],
                "wps_tf_std": [std_vec, std_vec],
                "n_samples": [5, 5],
            }
        )
        return WpsBaseline(regions=df, schema_version="2.0")

    def test_schema_version(self, sample_wps_v2):
        """Test schema_version is correctly set."""
        assert sample_wps_v2.schema_version == "2.0"

    def test_get_baseline_vector(self, sample_wps_v2):
        """Test get_baseline_vector returns 200-element array."""
        vec = sample_wps_v2.get_baseline_vector("gene1", "wps_nuc")
        assert vec is not None
        assert len(vec) == 200
        assert vec[0] == 0.0
        assert vec[100] == 100.0

    def test_get_std_vector(self, sample_wps_v2):
        """Test get_std_vector returns 200-element array."""
        std = sample_wps_v2.get_std_vector("gene1", "wps_nuc")
        assert std is not None
        assert len(std) == 200
        assert all(s == 0.1 for s in std)

    def test_get_baseline_vector_missing(self, sample_wps_v2):
        """Test get_baseline_vector returns None for missing region."""
        vec = sample_wps_v2.get_baseline_vector("unknown", "wps_nuc")
        assert vec is None

    def test_compute_z_vector(self, sample_wps_v2):
        """Test compute_z_vector returns position-wise z-scores."""
        # Sample vector with +1 offset from mean
        sample = np.array([float(i) + 1 for i in range(200)])  # 1, 2, 3, ..., 200
        z_vec = sample_wps_v2.compute_z_vector("gene1", sample, "wps_nuc")

        assert z_vec is not None
        assert len(z_vec) == 200
        # z = (sample - mean) / std = 1 / 0.1 = 10 for all positions
        assert all(pytest.approx(z, abs=0.01) == 10.0 for z in z_vec)

    def test_compute_z_vector_missing(self, sample_wps_v2):
        """Test compute_z_vector returns None for missing region."""
        sample = np.array([1.0] * 200)
        z_vec = sample_wps_v2.compute_z_vector("unknown", sample, "wps_nuc")
        assert z_vec is None

    def test_compute_shape_score_perfect(self, sample_wps_v2):
        """Test compute_shape_score returns 1.0 for identical shape."""
        # Use exact same shape as PON mean (should have correlation = 1.0)
        sample = np.array([float(i) for i in range(200)])
        score = sample_wps_v2.compute_shape_score("gene1", sample, "wps_nuc")

        assert score is not None
        assert pytest.approx(score, abs=0.01) == 1.0

    def test_compute_shape_score_offset(self, sample_wps_v2):
        """Test compute_shape_score ~1.0 for offset but same shape."""
        # Offset doesn't affect correlation - still same shape
        sample = np.array([float(i) + 100 for i in range(200)])
        score = sample_wps_v2.compute_shape_score("gene1", sample, "wps_nuc")

        assert score is not None
        # Correlation should still be very high (same trend, just shifted)
        assert pytest.approx(score, abs=0.01) == 1.0

    def test_compute_shape_score_inverted(self, sample_wps_v2):
        """Test compute_shape_score returns -1.0 for inverted shape."""
        # Inverted shape (negatively correlated)
        sample = np.array([float(199 - i) for i in range(200)])
        score = sample_wps_v2.compute_shape_score("gene1", sample, "wps_nuc")

        assert score is not None
        assert pytest.approx(score, abs=0.01) == -1.0

    def test_compute_shape_score_missing(self, sample_wps_v2):
        """Test compute_shape_score returns None for missing region."""
        sample = np.array([1.0] * 200)
        score = sample_wps_v2.compute_shape_score("unknown", sample, "wps_nuc")
        assert score is None


class TestOcfBaseline:
    """Tests for OcfBaseline class."""

    @pytest.fixture
    def sample_ocf_baseline(self):
        df = pd.DataFrame(
            {
                "region_id": ["promoter_1", "enhancer_1"],
                "ocf_mean": [0.5, 0.3],
                "ocf_std": [0.1, 0.05],
            }
        )
        return OcfBaseline(regions=df)

    def test_get_stats(self, sample_ocf_baseline):
        """Test get_stats returns (mean, std) tuple."""
        stats = sample_ocf_baseline.get_stats("promoter_1")
        assert stats[0] == 0.5
        assert stats[1] == 0.1

    def test_compute_zscore(self, sample_ocf_baseline):
        """Test z-score computation."""
        zscore = sample_ocf_baseline.compute_zscore("promoter_1", 0.7)
        assert pytest.approx(zscore) == 2.0  # (0.7 - 0.5) / 0.1 = 2.0


class TestMdsBaseline:
    """Tests for MdsBaseline class."""

    @pytest.fixture
    def sample_mds_baseline(self):
        return MdsBaseline(
            kmer_expected={"ACGT": 0.01, "TGCA": 0.02},
            kmer_std={"ACGT": 0.002, "TGCA": 0.004},
            mds_mean=4.5,
            mds_std=0.5,
        )

    def test_get_kmer_zscore(self, sample_mds_baseline):
        """Test k-mer z-score computation."""
        zscore = sample_mds_baseline.get_kmer_zscore("ACGT", 0.014)
        assert pytest.approx(zscore) == 2.0  # (0.014 - 0.01) / 0.002 = 2.0

    def test_get_mds_zscore(self, sample_mds_baseline):
        """Test MDS z-score computation."""
        zscore = sample_mds_baseline.get_mds_zscore(5.0)
        assert zscore == 1.0  # (5.0 - 4.5) / 0.5 = 1.0

    def test_get_aberrant_kmers(self, sample_mds_baseline):
        """Test aberrant k-mer detection."""
        observed = {"ACGT": 0.016, "TGCA": 0.02}  # ACGT has z=3, TGCA has z=0
        aberrant = sample_mds_baseline.get_aberrant_kmers(observed, threshold=2.0)
        assert "ACGT" in aberrant
        assert "TGCA" not in aberrant


class TestWpsBackgroundBaseline:
    """Tests for WpsBackgroundBaseline class."""

    @pytest.fixture
    def sample_wps_background(self):
        df = pd.DataFrame(
            {
                "group_id": ["all", "chr1"],
                "nrl_mean": [196.0, 195.0],
                "nrl_std": [3.0, 4.0],
                "periodicity_mean": [0.85, 0.82],
                "periodicity_std": [0.05, 0.06],
            }
        )
        return WpsBackgroundBaseline(groups=df)

    def test_get_nrl_stats(self, sample_wps_background):
        """Test get_nrl_stats returns (mean, std) tuple."""
        stats = sample_wps_background.get_nrl_stats("all")
        assert stats[0] == 196.0
        assert stats[1] == 3.0

    def test_get_nrl_stats_missing(self, sample_wps_background):
        """Test get_nrl_stats returns None for missing group."""
        stats = sample_wps_background.get_nrl_stats("chr99")
        assert stats is None

    def test_get_periodicity_stats(self, sample_wps_background):
        """Test get_periodicity_stats returns (mean, std) tuple."""
        stats = sample_wps_background.get_periodicity_stats("all")
        assert stats[0] == 0.85
        assert stats[1] == 0.05

    def test_compute_nrl_zscore(self, sample_wps_background):
        """Test NRL z-score computation."""
        # observed=202, mean=196, std=3 -> z=(202-196)/3=2.0
        zscore = sample_wps_background.compute_nrl_zscore(202.0, "all")
        assert pytest.approx(zscore) == 2.0

    def test_compute_nrl_zscore_missing(self, sample_wps_background):
        """Test NRL z-score returns None for missing group."""
        zscore = sample_wps_background.compute_nrl_zscore(200.0, "chr99")
        assert zscore is None


class TestRegionMdsBaseline:
    """Tests for RegionMdsBaseline class."""

    @pytest.fixture
    def sample_region_mds(self):
        return RegionMdsBaseline(
            gene_baseline={
                "TP53": {
                    "mds_mean": 4.5,
                    "mds_std": 0.3,
                    "mds_e1_mean": 4.2,
                    "mds_e1_std": 0.4,
                    "n_samples": 20,
                },
                "EGFR": {
                    "mds_mean": 4.8,
                    "mds_std": 0.5,
                    "mds_e1_mean": 4.6,
                    "mds_e1_std": 0.6,
                    "n_samples": 20,
                },
            }
        )

    def test_get_stats(self, sample_region_mds):
        """Test get_stats returns (mean, std) tuple."""
        stats = sample_region_mds.get_stats("TP53")
        assert stats[0] == 4.5
        assert stats[1] == 0.3

    def test_get_stats_missing(self, sample_region_mds):
        """Test get_stats returns None for missing gene."""
        stats = sample_region_mds.get_stats("UNKNOWN_GENE")
        assert stats is None

    def test_get_e1_stats(self, sample_region_mds):
        """Test get_e1_stats returns (mean, std) tuple."""
        stats = sample_region_mds.get_e1_stats("TP53")
        assert stats[0] == 4.2
        assert stats[1] == 0.4

    def test_compute_zscore(self, sample_region_mds):
        """Test z-score computation."""
        # observed=5.1, mean=4.5, std=0.3 -> z=(5.1-4.5)/0.3=2.0
        zscore = sample_region_mds.compute_zscore("TP53", 5.1)
        assert pytest.approx(zscore) == 2.0

    def test_compute_zscore_missing(self, sample_region_mds):
        """Test z-score returns None for missing gene."""
        zscore = sample_region_mds.compute_zscore("UNKNOWN_GENE", 5.0)
        assert zscore is None

    def test_compute_e1_zscore(self, sample_region_mds):
        """Test E1 z-score computation."""
        # observed=5.0, mean=4.2, std=0.4 -> z=(5.0-4.2)/0.4=2.0
        zscore = sample_region_mds.compute_e1_zscore("TP53", 5.0)
        assert pytest.approx(zscore) == 2.0


class TestFscGeneBaseline:
    """Tests for FscGeneBaseline class."""

    @pytest.fixture
    def sample_fsc_gene(self):
        return FscGeneBaseline(
            data={
                "TP53": (1.0, 0.1, 5),  # mean=1.0, std=0.1, n_samples=5
                "EGFR": (0.8, 0.15, 4),  # mean=0.8, std=0.15, n_samples=4
                "KRAS": (1.2, 0.05, 6),  # mean=1.2, std=0.05, n_samples=6
            }
        )

    def test_get_stats(self, sample_fsc_gene):
        """Test get_stats returns (mean, std) tuple."""
        stats = sample_fsc_gene.get_stats("TP53")
        assert stats == (1.0, 0.1)

    def test_get_stats_missing(self, sample_fsc_gene):
        """Test get_stats returns None for missing gene."""
        assert sample_fsc_gene.get_stats("UNKNOWN") is None

    def test_compute_zscore(self, sample_fsc_gene):
        """Test z-score computation."""
        # observed=1.2, mean=1.0, std=0.1 -> z=(1.2-1.0)/0.1=2.0
        zscore = sample_fsc_gene.compute_zscore("TP53", 1.2)
        assert pytest.approx(zscore) == 2.0

    def test_compute_zscore_missing(self, sample_fsc_gene):
        """Test z-score returns None for missing gene."""
        zscore = sample_fsc_gene.compute_zscore("UNKNOWN", 1.0)
        assert zscore is None

    def test_len(self, sample_fsc_gene):
        """Test __len__ returns gene count."""
        assert len(sample_fsc_gene) == 3


class TestFscRegionBaseline:
    """Tests for FscRegionBaseline class."""

    @pytest.fixture
    def sample_fsc_region(self):
        return FscRegionBaseline(
            data={
                "chr17:7571720-7573008": (1.0, 0.1, 5),
                "chr7:55086725-55086925": (0.9, 0.12, 4),
                "chr12:25398284-25398490": (1.1, 0.08, 6),
            }
        )

    def test_get_stats(self, sample_fsc_region):
        """Test get_stats returns (mean, std) tuple."""
        stats = sample_fsc_region.get_stats("chr17:7571720-7573008")
        assert stats == (1.0, 0.1)

    def test_get_stats_missing(self, sample_fsc_region):
        """Test get_stats returns None for missing region."""
        assert sample_fsc_region.get_stats("chr1:1-100") is None

    def test_compute_zscore(self, sample_fsc_region):
        """Test z-score computation."""
        # observed=1.2, mean=1.0, std=0.1 -> z=(1.2-1.0)/0.1=2.0
        zscore = sample_fsc_region.compute_zscore("chr17:7571720-7573008", 1.2)
        assert pytest.approx(zscore) == 2.0

    def test_compute_zscore_missing(self, sample_fsc_region):
        """Test z-score returns None for missing region."""
        zscore = sample_fsc_region.compute_zscore("chr1:1-100", 1.0)
        assert zscore is None

    def test_len(self, sample_fsc_region):
        """Test __len__ returns region count."""
        assert len(sample_fsc_region) == 3


class TestPonModel:
    """Tests for PonModel class."""

    @pytest.fixture
    def sample_gc_bias(self):
        return GcBiasModel(
            gc_bins=[0.3, 0.4, 0.45, 0.5, 0.6],
            short_expected=[0.8, 0.9, 1.0, 0.95, 0.85],
            short_std=[0.05, 0.06, 0.07, 0.06, 0.05],
            intermediate_expected=[0.85, 0.92, 1.0, 0.97, 0.88],
            intermediate_std=[0.04, 0.05, 0.06, 0.05, 0.04],
            long_expected=[0.9, 0.95, 1.0, 0.98, 0.92],
            long_std=[0.03, 0.04, 0.05, 0.04, 0.03],
        )

    @pytest.fixture
    def sample_pon_model(self, sample_gc_bias):
        return PonModel(gc_bias=sample_gc_bias, assay="test", n_samples=10)

    def test_get_mean(self, sample_pon_model):
        """Test get_mean returns expected value at median GC."""
        mean = sample_pon_model.get_mean("short")
        assert mean == 1.0  # Expected at GC=0.45

    def test_get_mean_channel_mapping(self, sample_pon_model):
        """Test get_mean maps FSC channels correctly."""
        # ultra_short should map to short
        mean = sample_pon_model.get_mean("ultra_short")
        assert mean == 1.0

    def test_get_variance(self, sample_pon_model):
        """Test get_variance returns squared std."""
        var = sample_pon_model.get_variance("short")
        # Median std is 0.07, variance = 0.07^2 = 0.0049
        assert pytest.approx(var, rel=0.01) == 0.0049

    def test_get_mean_no_gc_bias(self):
        """Test get_mean returns None without GC bias."""
        pon = PonModel()
        assert pon.get_mean("short") is None


class TestPonModelPanelMode:
    """Tests for PonModel panel mode fields."""

    def test_panel_mode_default(self):
        """Test panel_mode defaults to False."""
        pon = PonModel()
        assert not pon.panel_mode
        assert pon.target_regions_file == ""

    def test_panel_mode_set(self):
        """Test panel_mode can be set."""
        pon = PonModel(
            panel_mode=True, target_regions_file="targets.bed", assay="msk-access"
        )
        assert pon.panel_mode
        assert pon.target_regions_file == "targets.bed"


class TestBuildComputeFunctions:
    """Tests for build.py compute functions."""

    def test_compute_ocf_baseline(self):
        """Test _compute_ocf_baseline aggregates correctly."""
        import pandas as pd
        from krewlyzer.pon.build import _compute_ocf_baseline

        data = [
            pd.DataFrame({"region_id": ["R1", "R2"], "ocf": [0.5, 0.3]}),
            pd.DataFrame({"region_id": ["R1", "R2"], "ocf": [0.7, 0.5]}),
        ]
        result = _compute_ocf_baseline(data)

        assert result is not None
        assert len(result.regions) == 2
        # R1 mean should be (0.5 + 0.7) / 2 = 0.6
        r1 = result.regions[result.regions["region_id"] == "R1"].iloc[0]
        assert pytest.approx(r1["ocf_mean"], rel=0.01) == 0.6

    def test_compute_mds_baseline(self):
        """Test _compute_mds_baseline aggregates correctly."""
        from krewlyzer.pon.build import _compute_mds_baseline

        data = [
            {"kmers": {"ACGT": 0.01, "TGCA": 0.02}, "mds": 4.0},
            {"kmers": {"ACGT": 0.02, "TGCA": 0.03}, "mds": 5.0},
        ]
        result = _compute_mds_baseline(data)

        assert result is not None
        assert len(result.kmer_expected) == 2
        assert pytest.approx(result.kmer_expected["ACGT"], rel=0.01) == 0.015
        assert pytest.approx(result.mds_mean, rel=0.01) == 4.5
