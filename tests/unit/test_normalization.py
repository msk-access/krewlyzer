"""
Unit tests for depth normalization calculations.

Covers normalization formulas used across tools.
"""

import pytest
import numpy as np
import json


@pytest.mark.unit
def test_per_million_normalization():
    """Test standard per-million normalization."""
    raw_value = 1000
    total_fragments = 5_000_000

    normalized = raw_value / (total_fragments / 1_000_000)

    assert normalized == 200.0


@pytest.mark.unit
def test_normalization_preserves_relative_values():
    """Test normalization preserves relative ordering."""
    values = [100, 200, 300]
    total = 2_000_000

    factor = total / 1_000_000
    normalized = [v / factor for v in values]

    # Relative ordering preserved
    assert normalized[0] < normalized[1] < normalized[2]


@pytest.mark.unit
def test_gc_correction_weight():
    """Test GC correction weight calculation."""
    # GC weight = expected_count / observed_count
    expected = 100
    observed = 80

    weight = expected / observed if observed > 0 else 1.0

    assert weight == 1.25
    assert weight > 1.0  # Under-represented GC bin


@pytest.mark.unit
def test_gc_correction_underrepresented():
    """Test GC correction for underrepresented bins."""
    expected = 100
    observed = 200

    weight = expected / observed if observed > 0 else 1.0

    assert weight == 0.5
    assert weight < 1.0  # Over-represented GC bin


@pytest.mark.unit
def test_metadata_json_format():
    """Test metadata.json structure for normalization."""
    metadata = {
        "sample_id": "test_sample",
        "total_fragments": 7500000,
        "filters": {"mapq": 20, "min_length": 65, "max_length": 400},
        "timestamp": "2025-12-17T09:00:00",
    }

    # Verify required fields
    assert "total_fragments" in metadata
    assert isinstance(metadata["total_fragments"], int)
    assert metadata["total_fragments"] > 0


@pytest.mark.unit
def test_metadata_file_flow(tmp_path):
    """Test metadata file read for normalization."""
    sample_name = "test_sample"
    metadata = {
        "sample_id": sample_name,
        "total_fragments": 5_000_000,
        "filters": {"mapq": 20},
        "timestamp": "2025-12-17T09:00:00",
    }

    metadata_file = tmp_path / f"{sample_name}.metadata.json"
    with open(metadata_file, "w") as f:
        json.dump(metadata, f, indent=2)

    # Load and use for normalization
    with open(metadata_file) as f:
        loaded = json.load(f)

    assert loaded["total_fragments"] == 5_000_000

    # Calculate normalization
    raw = 30
    factor = loaded["total_fragments"] / 1_000_000
    normalized = raw / factor

    assert np.isclose(normalized, 6.0)
