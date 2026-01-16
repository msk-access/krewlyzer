"""
Unit tests for PON model classes.

Tests the new baseline classes and methods added in the PON enhancement.
"""

import pytest
import numpy as np
import pandas as pd
from krewlyzer.pon.model import (
    PonModel, GcBiasModel, FsdBaseline, WpsBaseline,
    OcfBaseline, MdsBaseline
)


class TestFsdBaseline:
    """Tests for FsdBaseline class."""
    
    @pytest.fixture
    def sample_fsd_baseline(self):
        return FsdBaseline(
            size_bins=[65, 70, 75, 80, 85],
            arms={
                '1p': {'expected': [0.1, 0.15, 0.2, 0.15, 0.1], 'std': [0.01, 0.02, 0.03, 0.02, 0.01]},
                '1q': {'expected': [0.12, 0.18, 0.22, 0.18, 0.12], 'std': [0.02, 0.03, 0.04, 0.03, 0.02]},
            }
        )
    
    def test_get_expected(self, sample_fsd_baseline):
        """Test get_expected returns correct value."""
        expected = sample_fsd_baseline.get_expected('1p', 75)
        assert expected == 0.2
    
    def test_get_std(self, sample_fsd_baseline):
        """Test get_std returns correct value."""
        std = sample_fsd_baseline.get_std('1p', 75)
        assert std == 0.03
    
    def test_get_stats(self, sample_fsd_baseline):
        """Test get_stats returns (expected, std) tuple."""
        stats = sample_fsd_baseline.get_stats('1p', 75)
        assert stats == (0.2, 0.03)
    
    def test_get_stats_missing_arm(self, sample_fsd_baseline):
        """Test get_stats returns None for missing arm."""
        stats = sample_fsd_baseline.get_stats('unknown', 75)
        assert stats is None


class TestWpsBaseline:
    """Tests for WpsBaseline class."""
    
    @pytest.fixture
    def sample_wps_v1(self):
        df = pd.DataFrame({
            'region_id': ['gene1', 'gene2'],
            'wps_long_mean': [100.0, 120.0],
            'wps_long_std': [10.0, 12.0],
            'wps_short_mean': [50.0, 55.0],
            'wps_short_std': [5.0, 6.0]
        })
        return WpsBaseline(regions=df, schema_version='1.0')
    
    def test_get_baseline(self, sample_wps_v1):
        """Test get_baseline returns dict for region."""
        baseline = sample_wps_v1.get_baseline('gene1')
        assert baseline is not None
        assert baseline['wps_long_mean'] == 100.0
    
    def test_get_stats(self, sample_wps_v1):
        """Test get_stats returns (mean, std) tuple."""
        stats = sample_wps_v1.get_stats('gene1', 'wps_long')
        assert stats == (100.0, 10.0)
    
    def test_v1_vector_returns_none(self, sample_wps_v1):
        """Test v1.0 schema returns None for vector methods."""
        assert sample_wps_v1.get_baseline_vector('gene1') is None


class TestOcfBaseline:
    """Tests for OcfBaseline class."""
    
    @pytest.fixture
    def sample_ocf_baseline(self):
        df = pd.DataFrame({
            'region_id': ['promoter_1', 'enhancer_1'],
            'ocf_mean': [0.5, 0.3],
            'ocf_std': [0.1, 0.05]
        })
        return OcfBaseline(regions=df)
    
    def test_get_stats(self, sample_ocf_baseline):
        """Test get_stats returns (mean, std) tuple."""
        stats = sample_ocf_baseline.get_stats('promoter_1')
        assert stats[0] == 0.5
        assert stats[1] == 0.1
    
    def test_compute_zscore(self, sample_ocf_baseline):
        """Test z-score computation."""
        zscore = sample_ocf_baseline.compute_zscore('promoter_1', 0.7)
        assert pytest.approx(zscore) == 2.0  # (0.7 - 0.5) / 0.1 = 2.0


class TestMdsBaseline:
    """Tests for MdsBaseline class."""
    
    @pytest.fixture
    def sample_mds_baseline(self):
        return MdsBaseline(
            kmer_expected={'ACGT': 0.01, 'TGCA': 0.02},
            kmer_std={'ACGT': 0.002, 'TGCA': 0.004},
            mds_mean=4.5,
            mds_std=0.5
        )
    
    def test_get_kmer_zscore(self, sample_mds_baseline):
        """Test k-mer z-score computation."""
        zscore = sample_mds_baseline.get_kmer_zscore('ACGT', 0.014)
        assert pytest.approx(zscore) == 2.0  # (0.014 - 0.01) / 0.002 = 2.0
    
    def test_get_mds_zscore(self, sample_mds_baseline):
        """Test MDS z-score computation."""
        zscore = sample_mds_baseline.get_mds_zscore(5.0)
        assert zscore == 1.0  # (5.0 - 4.5) / 0.5 = 1.0
    
    def test_get_aberrant_kmers(self, sample_mds_baseline):
        """Test aberrant k-mer detection."""
        observed = {'ACGT': 0.016, 'TGCA': 0.02}  # ACGT has z=3, TGCA has z=0
        aberrant = sample_mds_baseline.get_aberrant_kmers(observed, threshold=2.0)
        assert 'ACGT' in aberrant
        assert 'TGCA' not in aberrant


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
            long_std=[0.03, 0.04, 0.05, 0.04, 0.03]
        )
    
    @pytest.fixture
    def sample_pon_model(self, sample_gc_bias):
        return PonModel(
            gc_bias=sample_gc_bias,
            assay='test',
            n_samples=10
        )
    
    def test_get_mean(self, sample_pon_model):
        """Test get_mean returns expected value at median GC."""
        mean = sample_pon_model.get_mean('short')
        assert mean == 1.0  # Expected at GC=0.45
    
    def test_get_mean_channel_mapping(self, sample_pon_model):
        """Test get_mean maps FSC channels correctly."""
        # ultra_short should map to short
        mean = sample_pon_model.get_mean('ultra_short')
        assert mean == 1.0
    
    def test_get_variance(self, sample_pon_model):
        """Test get_variance returns squared std."""
        var = sample_pon_model.get_variance('short')
        # Median std is 0.07, variance = 0.07^2 = 0.0049
        assert pytest.approx(var, rel=0.01) == 0.0049
    
    def test_get_mean_no_gc_bias(self):
        """Test get_mean returns None without GC bias."""
        pon = PonModel()
        assert pon.get_mean('short') is None


class TestPonModelPanelMode:
    """Tests for PonModel panel mode fields."""
    
    def test_panel_mode_default(self):
        """Test panel_mode defaults to False."""
        pon = PonModel()
        assert pon.panel_mode == False
        assert pon.target_regions_file == ""
    
    def test_panel_mode_set(self):
        """Test panel_mode can be set."""
        pon = PonModel(
            panel_mode=True,
            target_regions_file="targets.bed",
            assay="msk-access"
        )
        assert pon.panel_mode == True
        assert pon.target_regions_file == "targets.bed"


class TestBuildComputeFunctions:
    """Tests for build.py compute functions."""
    
    def test_compute_ocf_baseline(self):
        """Test _compute_ocf_baseline aggregates correctly."""
        import pandas as pd
        from krewlyzer.pon.build import _compute_ocf_baseline
        
        data = [
            pd.DataFrame({'region_id': ['R1', 'R2'], 'ocf': [0.5, 0.3]}),
            pd.DataFrame({'region_id': ['R1', 'R2'], 'ocf': [0.7, 0.5]}),
        ]
        result = _compute_ocf_baseline(data)
        
        assert result is not None
        assert len(result.regions) == 2
        # R1 mean should be (0.5 + 0.7) / 2 = 0.6
        r1 = result.regions[result.regions['region_id'] == 'R1'].iloc[0]
        assert pytest.approx(r1['ocf_mean'], rel=0.01) == 0.6
    
    def test_compute_mds_baseline(self):
        """Test _compute_mds_baseline aggregates correctly."""
        from krewlyzer.pon.build import _compute_mds_baseline
        
        data = [
            {'kmers': {'ACGT': 0.01, 'TGCA': 0.02}, 'mds': 4.0},
            {'kmers': {'ACGT': 0.02, 'TGCA': 0.03}, 'mds': 5.0},
        ]
        result = _compute_mds_baseline(data)
        
        assert result is not None
        assert len(result.kmer_expected) == 2
        assert pytest.approx(result.kmer_expected['ACGT'], rel=0.01) == 0.015
        assert pytest.approx(result.mds_mean, rel=0.01) == 4.5
