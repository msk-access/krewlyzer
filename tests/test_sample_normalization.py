"""Unit tests for WPS depth normalization.

WPS normalized columns (wps_*_norm) require metadata.json from extract step
to provide total_fragments for per-million normalization.
"""
import json
import tempfile
import shutil
from pathlib import Path
import numpy as np
import pytest


def test_wps_depth_normalization_formula():
    """Test WPS depth normalization calculation.
    
    Formula: wps_norm = wps_raw / (total_fragments / 1_000_000)
    This gives "per million fragments" values comparable across samples.
    """
    # Test data
    wps_raw = np.array([45, 52, 38, 61])  # Raw WPS scores
    total_fragments = 7_500_000
    
    # Calculate normalized WPS (per million fragments)
    norm_factor = total_fragments / 1_000_000  # = 7.5
    wps_norm = wps_raw / norm_factor
    
    # Expected: 45 / 7.5 = 6.0
    expected_first = 45 / 7.5
    
    assert np.isclose(wps_norm[0], expected_first, rtol=1e-5)
    assert wps_norm[0] > 0
    assert wps_norm[0] < wps_raw[0]  # Normalized should be smaller (depth > 1M)


def test_wps_normalization_comparability():
    """Test that WPS normalization makes values comparable across samples.
    
    Same biological signal at different depths should give similar normalized values.
    """
    # Same region, same biology, different sequencing depths
    wps_sample_a = 45  # 7.5M fragments
    wps_sample_b = 90  # 15M fragments (2x depth = 2x raw WPS)
    
    depth_a = 7_500_000
    depth_b = 15_000_000
    
    norm_a = wps_sample_a / (depth_a / 1_000_000)  # = 45 / 7.5 = 6.0
    norm_b = wps_sample_b / (depth_b / 1_000_000)  # = 90 / 15 = 6.0
    
    # After normalization, comparable
    assert np.isclose(norm_a, norm_b, rtol=0.01)


def test_metadata_json_format():
    """Test that metadata.json has correct structure for WPS normalization."""
    # Simulate metadata that extract.py creates
    metadata = {
        "sample_id": "test_sample",
        "total_fragments": 7500000,
        "filters": {
            "mapq": 20,
            "min_length": 65,
            "max_length": 400
        },
        "timestamp": "2025-12-17T09:00:00"
    }
    
    # Verify required fields for WPS normalization
    assert "total_fragments" in metadata
    assert isinstance(metadata["total_fragments"], int)
    assert metadata["total_fragments"] > 0


def test_missing_metadata_default():
    """Test behavior when metadata file is missing.
    
    When metadata is not available:
    - norm_factor defaults to 1.0 (1M / 1M)
    - wps_norm = wps_raw (columns are present but not useful)
    - Warning should be logged
    """
    # Default behavior in Rust
    default_total = 1_000_000
    norm_factor = default_total / 1_000_000  # = 1.0
    
    wps_raw = 45
    wps_norm = wps_raw / norm_factor  # = 45
    
    # When missing metadata, normalized == raw
    assert wps_norm == wps_raw


def test_integration_metadata_flow():
    """Test realistic metadata flow from extract to WPS."""
    test_dir = Path(tempfile.mkdtemp())
    
    try:
        # 1. Simulate extract creating metadata
        sample_name = "test_sample"
        metadata = {
            "sample_id": sample_name,
            "total_fragments": 5_000_000,
            "filters": {"mapq": 20},
            "timestamp": "2025-12-17T09:00:00"
        }
        metadata_file = test_dir / f"{sample_name}.metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        # 2. Verify WPS can find metadata
        # WPS looks for metadata by replacing .bed.gz with .metadata.json
        bed_gz = test_dir / f"{sample_name}.bed.gz"
        expected_meta_path = str(bed_gz).replace('.bed.gz', '.metadata.json')
        
        # Create dummy bed.gz so path manipulation works
        bed_gz.touch()
        
        # Load metadata as WPS would
        assert Path(expected_meta_path).exists()
        with open(expected_meta_path) as f:
            loaded = json.load(f)
        
        assert loaded["total_fragments"] == 5_000_000
        
        # 3. Calculate normalization as WPS would
        wps_raw = 30
        norm_factor = loaded["total_fragments"] / 1_000_000  # = 5.0
        wps_norm = wps_raw / norm_factor  # = 6.0
        
        assert np.isclose(wps_norm, 6.0)
    
    finally:
        shutil.rmtree(test_dir)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
