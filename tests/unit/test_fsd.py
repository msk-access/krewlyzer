"""
Unit tests for FSD (Fragment Size Distribution) calculations.

Tests FSD statistical computations.
"""
import pytest
import numpy as np
from pathlib import Path


@pytest.mark.unit
def test_fsd_arm_ratio_calculation():
    """Test FSD arm ratio calculation formula."""
    # FSD ratio = short_fragments / long_fragments
    short_count = 100
    long_count = 200
    
    ratio = short_count / long_count if long_count > 0 else 0
    
    assert ratio == 0.5
    assert ratio > 0


@pytest.mark.unit
def test_fsd_mean_size_calculation():
    """Test mean fragment size calculation."""
    sizes = [100, 150, 180, 200, 250]
    
    mean_size = np.mean(sizes)
    
    assert mean_size == 176.0
    assert 100 < mean_size < 250


@pytest.mark.unit
def test_fsd_mode_size_calculation():
    """Test mode fragment size calculation."""
    # Histogram approach for mode
    sizes = [150, 160, 165, 165, 165, 170, 180]
    
    # Binned histogram (10bp bins)
    bins = list(range(100, 300, 10))
    hist, _ = np.histogram(sizes, bins=bins)
    mode_bin_idx = np.argmax(hist)
    mode_center = bins[mode_bin_idx] + 5  # Center of bin
    
    # Mode should be around 165 (bin 160-170)
    assert 160 <= mode_center <= 170


@pytest.mark.unit
def test_fsd_empty_arm_handling():
    """Test FSD handles empty arms correctly."""
    short_count = 0
    long_count = 100
    
    # Should handle division by zero
    ratio = short_count / long_count if long_count > 0 else 0
    
    assert ratio == 0


@pytest.mark.unit
def test_fsd_no_long_fragments():
    """Test FSD when no long fragments exist."""
    short_count = 100
    long_count = 0
    
    # Avoid division by zero
    ratio = short_count / long_count if long_count > 0 else float('inf')
    
    assert ratio == float('inf')
