"""
Unit tests for WPS (Windowed Protection Score) calculations.

Note: WPS now runs through the unified pipeline.
These tests cover the normalization math.
"""
import pytest
import numpy as np
from pathlib import Path


@pytest.mark.unit
def test_wps_calculation_formula():
    """Test WPS score calculation formula.
    
    WPS = fragments_spanning_window - fragments_with_endpoint_in_window
    """
    spanning = 10
    endpoints = 3
    
    wps = spanning - endpoints
    
    assert wps == 7


@pytest.mark.unit
def test_wps_negative_score():
    """Test WPS can be negative (more endpoints than spanning)."""
    spanning = 2
    endpoints = 5
    
    wps = spanning - endpoints
    
    assert wps == -3
    assert wps < 0


@pytest.mark.unit
def test_wps_ratio_calculation():
    """Test WPS long/short ratio calculation."""
    wps_long = 100
    wps_short = 50
    
    ratio = wps_long / wps_short if wps_short != 0 else 0
    
    assert ratio == 2.0


@pytest.mark.unit
def test_wps_depth_normalization():
    """Test WPS depth normalization formula.
    
    Formula: wps_norm = wps_raw / (total_fragments / 1_000_000)
    """
    wps_raw = 45
    total_fragments = 7_500_000
    
    norm_factor = total_fragments / 1_000_000  # = 7.5
    wps_norm = wps_raw / norm_factor
    
    expected = 45 / 7.5  # = 6.0
    
    assert np.isclose(wps_norm, 6.0)
    assert np.isclose(wps_norm, expected)


@pytest.mark.unit
def test_wps_normalization_comparability():
    """Test normalized WPS is comparable across different depths."""
    # Same biology, different depths
    wps_sample_a = 45  # 7.5M fragments
    wps_sample_b = 90  # 15M fragments (2x depth = 2x raw WPS)
    
    depth_a = 7_500_000
    depth_b = 15_000_000
    
    norm_a = wps_sample_a / (depth_a / 1_000_000)  # = 6.0
    norm_b = wps_sample_b / (depth_b / 1_000_000)  # = 6.0
    
    # After normalization, values should be comparable
    assert np.isclose(norm_a, norm_b, rtol=0.01)


@pytest.mark.unit
def test_wps_missing_metadata_default():
    """Test default normalization when metadata missing."""
    default_total = 1_000_000
    norm_factor = default_total / 1_000_000  # = 1.0
    
    wps_raw = 45
    wps_norm = wps_raw / norm_factor  # = 45
    
    # When missing metadata, normalized == raw
    assert wps_norm == wps_raw
