"""
Unit tests for WPS (Windowed Protection Score) calculations.

Tests cover:
- WPS formula math
- Dual-stream weighted fragment classification
- Savitzky-Golay smoothing
- FFT periodicity extraction
"""
import pytest
import numpy as np
from pathlib import Path


# ============================================================================
# WPS Formula Tests
# ============================================================================

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


# ============================================================================
# Dual-Stream Weighted Fragment Classification Tests
# ============================================================================

def _nuc_weight(length: int) -> float:
    """Python reference implementation of nucleosome weight calculation.
    
    Mirrors WpsConfig.nuc_weight() in Rust:
    - 1.0 for primary range [160, 175]
    - 0.5 for secondary range [120, 159] ∪ [176, 180]
    - 0.0 otherwise
    """
    if 160 <= length <= 175:
        return 1.0
    elif 120 <= length <= 180:
        return 0.5
    else:
        return 0.0


def _tf_weight(length: int) -> float:
    """Python reference implementation of TF weight calculation.
    
    Mirrors WpsConfig.tf_weight() in Rust:
    - 1.0 for [35, 80]
    - 0.0 otherwise
    """
    if 35 <= length <= 80:
        return 1.0
    else:
        return 0.0


@pytest.mark.unit
def test_wps_nuc_weight_primary_range():
    """Test nucleosome weight for primary range [160, 175]."""
    assert _nuc_weight(160) == 1.0
    assert _nuc_weight(167) == 1.0
    assert _nuc_weight(175) == 1.0


@pytest.mark.unit
def test_wps_nuc_weight_secondary_range():
    """Test nucleosome weight for secondary ranges [120,159] and [176,180]."""
    # Left secondary [120, 159]
    assert _nuc_weight(120) == 0.5
    assert _nuc_weight(140) == 0.5
    assert _nuc_weight(159) == 0.5
    
    # Right secondary [176, 180]
    assert _nuc_weight(176) == 0.5
    assert _nuc_weight(180) == 0.5


@pytest.mark.unit
def test_wps_nuc_weight_outside_range():
    """Test nucleosome weight is 0 outside valid ranges."""
    assert _nuc_weight(100) == 0.0
    assert _nuc_weight(119) == 0.0
    assert _nuc_weight(181) == 0.0
    assert _nuc_weight(200) == 0.0


@pytest.mark.unit
def test_wps_tf_weight_valid_range():
    """Test TF weight for valid range [35, 80]."""
    assert _tf_weight(35) == 1.0
    assert _tf_weight(55) == 1.0
    assert _tf_weight(80) == 1.0


@pytest.mark.unit
def test_wps_tf_weight_outside_range():
    """Test TF weight is 0 outside valid range."""
    assert _tf_weight(34) == 0.0
    assert _tf_weight(81) == 0.0
    assert _tf_weight(100) == 0.0


@pytest.mark.unit
def test_wps_fragment_dual_stream_independent():
    """Verify nucleosome and TF weights are independent."""
    # Short fragments: TF only
    assert _nuc_weight(50) == 0.0
    assert _tf_weight(50) == 1.0
    
    # Long fragments: Nuc only
    assert _nuc_weight(167) == 1.0
    assert _tf_weight(167) == 0.0
    
    # Gap between ranges: neither stream
    assert _nuc_weight(100) == 0.0
    assert _tf_weight(100) == 0.0


# ============================================================================
# Note: Savitzky-Golay smoothing and FFT periodicity tests removed
# These functions are now handled by Rust (sci-rs and realfft crates)
# See rust/src/wps.rs for implementation
# ============================================================================


# ============================================================================
# Adaptive Bait Padding Tests
# ============================================================================

def _adaptive_trim(bait_start: int, bait_end: int, user_trim: int) -> int:
    """Python reference implementation of adaptive bait padding.
    
    Mirrors the Rust check_position() logic:
    effective_trim = min(user_trim, bait_length / 4)
    """
    bait_length = bait_end - bait_start
    max_safe_trim = bait_length // 4  # Never trim more than 25% per side
    return min(user_trim, max_safe_trim)


@pytest.mark.unit
def test_adaptive_trim_large_target():
    """Large targets get full user-specified trim."""
    # 1000bp target, user wants 50bp trim
    effective = _adaptive_trim(0, 1000, 50)
    assert effective == 50  # 50 < 250 (1000/4), so full trim


@pytest.mark.unit
def test_adaptive_trim_small_exon():
    """Small exons get reduced trim to preserve data."""
    # 120bp exon, user wants 50bp trim
    effective = _adaptive_trim(0, 120, 50)
    assert effective == 30  # 120/4 = 30 < 50, so reduced


@pytest.mark.unit
def test_adaptive_trim_very_small_exon():
    """Very small exons get minimal trim."""
    # 80bp exon, user wants 50bp trim
    effective = _adaptive_trim(0, 80, 50)
    assert effective == 20  # 80/4 = 20 < 50


@pytest.mark.unit
def test_adaptive_trim_user_zero():
    """User can disable trimming entirely."""
    effective = _adaptive_trim(0, 1000, 0)
    assert effective == 0  # User wants no trim


@pytest.mark.unit
def test_adaptive_trim_preserves_at_least_50_percent():
    """Adaptive trim never masks more than 50% of target."""
    for length in [50, 100, 200, 500, 1000]:
        effective = _adaptive_trim(0, length, 100)
        # With max_safe_trim = length/4, trimming both sides
        # means we keep at least 50% of the region
        masked_pct = (2 * effective) / length
        assert masked_pct <= 0.5, f"Masked {masked_pct*100}% for {length}bp target"
