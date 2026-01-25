"""
Unit tests for PON validation module.

Tests panel configuration validation and on-target rate calculation.
"""

import pytest
from pathlib import Path
import tempfile
import gzip

from krewlyzer.pon.validation import (
    validate_panel_config,
    calculate_on_target_rate,
    ValidationResult,
)


# =============================================================================
# Mock PON Model for testing
# =============================================================================

class MockPonModel:
    """Mock PON model for testing validation."""
    
    def __init__(self, panel_mode=False, assay="", gc_bias_ontarget=None):
        self.panel_mode = panel_mode
        self.assay = assay
        self.gc_bias_ontarget = gc_bias_ontarget


# =============================================================================
# validate_panel_config tests
# =============================================================================

class TestValidatePanelConfig:
    """Tests for PON/panel configuration validation."""
    
    def test_wgs_pon_wgs_sample_valid(self):
        """WGS PON + WGS sample = valid, no warnings."""
        pon = MockPonModel(panel_mode=False)
        result = validate_panel_config(pon, target_regions=None)
        
        assert result.valid
        assert len(result.warnings) == 0
        assert len(result.errors) == 0
    
    def test_panel_pon_panel_sample_valid(self, tmp_path):
        """Panel PON + panel sample = valid."""
        pon = MockPonModel(panel_mode=True, gc_bias_ontarget="mock")
        target = tmp_path / "targets.bed"
        target.touch()
        
        result = validate_panel_config(pon, target_regions=target)
        
        assert result.valid
        assert len(result.errors) == 0
    
    def test_panel_pon_wgs_sample_warns(self):
        """Panel PON + WGS sample = warning."""
        pon = MockPonModel(panel_mode=True)
        result = validate_panel_config(pon, target_regions=None)
        
        assert result.valid  # Warning, not error
        assert len(result.warnings) > 0
        assert "without --target-regions" in result.warnings[0]
    
    def test_wgs_pon_panel_sample_warns(self, tmp_path):
        """WGS PON + panel sample = warning."""
        pon = MockPonModel(panel_mode=False)
        target = tmp_path / "targets.bed"
        target.touch()
        
        result = validate_panel_config(pon, target_regions=target)
        
        assert result.valid  # Warning, not error
        assert len(result.warnings) > 0
        assert "WGS mode" in result.warnings[0]
    
    def test_assay_mismatch_is_error(self):
        """Assay mismatch = error."""
        pon = MockPonModel(panel_mode=True, assay="xs1")
        result = validate_panel_config(pon, target_regions=None, assay="xs2")
        
        assert not result.valid
        assert len(result.errors) > 0
        assert "Assay mismatch" in result.errors[0]
    
    def test_missing_ontarget_gc_warns(self, tmp_path):
        """Panel PON without gc_bias_ontarget = warning."""
        pon = MockPonModel(panel_mode=True, gc_bias_ontarget=None)
        target = tmp_path / "targets.bed"
        target.touch()
        
        result = validate_panel_config(pon, target_regions=target)
        
        assert result.valid
        assert any("gc_bias_ontarget" in w for w in result.warnings)


# =============================================================================
# calculate_on_target_rate tests
# =============================================================================

class TestCalculateOnTargetRate:
    """Tests for on-target rate calculation."""
    
    def test_all_on_target(self, tmp_path):
        """All fragments on target = 100%."""
        # Create target regions
        targets = tmp_path / "targets.bed"
        targets.write_text("chr1\t0\t1000\tgene1\n")
        
        # Create fragments all within target
        frags = tmp_path / "frags.bed"
        frags.write_text("chr1\t100\t250\nchr1\t200\t350\nchr1\t300\t450\n")
        
        rate = calculate_on_target_rate(frags, targets, sample_size=100)
        assert rate == 1.0
    
    def test_no_on_target(self, tmp_path):
        """No fragments on target = 0%."""
        targets = tmp_path / "targets.bed"
        targets.write_text("chr1\t0\t100\tgene1\n")
        
        frags = tmp_path / "frags.bed"
        frags.write_text("chr1\t500\t650\nchr1\t600\t750\n")
        
        rate = calculate_on_target_rate(frags, targets, sample_size=100)
        assert rate == 0.0
    
    def test_partial_on_target(self, tmp_path):
        """Some fragments on target."""
        targets = tmp_path / "targets.bed"
        targets.write_text("chr1\t0\t200\tgene1\n")
        
        frags = tmp_path / "frags.bed"
        frags.write_text("chr1\t50\t150\nchr1\t500\t650\n")  # 1 on, 1 off
        
        rate = calculate_on_target_rate(frags, targets, sample_size=100)
        assert rate == 0.5
    
    def test_gzipped_fragments(self, tmp_path):
        """Works with gzipped fragment file."""
        targets = tmp_path / "targets.bed"
        targets.write_text("chr1\t0\t500\tgene1\n")
        
        frags = tmp_path / "frags.bed.gz"
        with gzip.open(frags, 'wt') as f:
            f.write("chr1\t100\t250\nchr1\t200\t350\n")
        
        rate = calculate_on_target_rate(frags, targets, sample_size=100)
        assert rate == 1.0


# =============================================================================
# --skip-pon flag tests
# =============================================================================

class TestSkipPonValidation:
    """Tests for --skip-pon flag handling and mutual exclusivity."""
    
    def test_skip_pon_flag_valid_alone(self):
        """--skip-pon without -P is valid."""
        from krewlyzer.pon.validation import validate_skip_pon
        
        result = validate_skip_pon(skip_pon=True, pon_model=None)
        assert result.valid
        assert len(result.errors) == 0
    
    def test_pon_model_valid_alone(self, tmp_path):
        """-P without --skip-pon is valid."""
        from krewlyzer.pon.validation import validate_skip_pon
        
        pon_path = tmp_path / "test.pon.parquet"
        pon_path.touch()
        
        result = validate_skip_pon(skip_pon=False, pon_model=pon_path)
        assert result.valid
        assert len(result.errors) == 0
    
    def test_skip_pon_and_pon_model_mutually_exclusive(self, tmp_path):
        """-P and --skip-pon together should raise error."""
        from krewlyzer.pon.validation import validate_skip_pon
        
        pon_path = tmp_path / "test.pon.parquet"
        pon_path.touch()
        
        result = validate_skip_pon(skip_pon=True, pon_model=pon_path)
        assert not result.valid
        assert len(result.errors) > 0
        assert "mutually exclusive" in result.errors[0].lower()
    
    def test_neither_skip_pon_nor_pon_model_valid(self):
        """Neither -P nor --skip-pon is valid (normal mode)."""
        from krewlyzer.pon.validation import validate_skip_pon
        
        result = validate_skip_pon(skip_pon=False, pon_model=None)
        assert result.valid
        assert len(result.errors) == 0
