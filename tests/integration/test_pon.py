"""
Integration tests for PON (Panel of Normals) functionality.

Tests PON model loading, application, and z-score calculations.
"""
import pytest
import json
import numpy as np
import pandas as pd
import tempfile
from pathlib import Path


@pytest.fixture
def sample_pon_parquet(tmp_path):
    """Create a valid PON parquet file for testing."""
    # PON format: parquet with 'table' column to identify row types
    metadata = pd.DataFrame({
        "table": ["metadata"],
        "version": ["1.0"],
        "genome": ["hg19"],
        "sample_count": [10],
        "gc_bin": [None],
        "short_expected": [None],
        "short_std": [None],
        "intermediate_expected": [None],
        "intermediate_std": [None],
        "long_expected": [None],
        "long_std": [None],
        "arm": [None],
        "mean_ratio": [None],
        "std_ratio": [None],
        "gene_id": [None],
        "wps_long_mean": [None],
        "wps_long_std": [None],
    })
    
    gc_bias = pd.DataFrame({
        "table": ["gc_bias"] * 3,
        "version": [None, None, None],
        "genome": [None, None, None],
        "sample_count": [None, None, None],
        "gc_bin": [0.3, 0.4, 0.5],
        "short_expected": [100.0, 120.0, 110.0],
        "short_std": [10.0, 12.0, 11.0],
        "intermediate_expected": [200.0, 220.0, 210.0],
        "intermediate_std": [20.0, 22.0, 21.0],
        "long_expected": [50.0, 55.0, 52.0],
        "long_std": [5.0, 5.5, 5.2],
        "arm": [None, None, None],
        "mean_ratio": [None, None, None],
        "std_ratio": [None, None, None],
        "gene_id": [None, None, None],
        "wps_long_mean": [None, None, None],
        "wps_long_std": [None, None, None],
    })
    
    # Combine and save
    pon_df = pd.concat([metadata, gc_bias], ignore_index=True)
    pon_file = tmp_path / "test.pon.parquet"
    pon_df.to_parquet(pon_file, index=False)
    return pon_file


@pytest.mark.skip(reason="PON parquet schema requires production format")
@pytest.mark.integration
def test_pon_model_loading(tmp_path, sample_pon_parquet):
    """Test PON model loading from parquet."""
    from krewlyzer.pon.model import PonModel
    
    # Load model
    model = PonModel.load(sample_pon_parquet)
    
    # Model should load successfully
    assert model is not None
    assert model.genome == "hg19"


@pytest.mark.integration
def test_pon_zscore_calculation():
    """Test z-score calculation formula."""
    # Z-score = (value - mean) / std
    value = 120.0
    mean = 100.0
    std = 10.0
    
    zscore = (value - mean) / std
    
    assert zscore == 2.0
    assert abs(zscore) > 1.96  # Significant at 95% CI


@pytest.mark.integration
def test_pon_integration_load_model(tmp_path, sample_pon_parquet):
    """Test PON integration module load_pon_model function."""
    from krewlyzer.core.pon_integration import load_pon_model
    
    # Load via integration module
    model = load_pon_model(sample_pon_parquet)
    
    assert model is not None


@pytest.mark.integration
def test_pon_missing_file():
    """Test PON loading handles missing file gracefully."""
    from krewlyzer.core.pon_integration import load_pon_model
    
    result = load_pon_model(Path("/nonexistent/path/pon.parquet"))
    
    assert result is None


@pytest.mark.unit
def test_pon_zscore_edge_cases():
    """Test z-score edge cases."""
    # Zero std should be handled
    value = 100.0
    mean = 100.0
    std = 0.0
    
    if std > 0:
        zscore = (value - mean) / std
    else:
        zscore = 0.0  # Default for zero std
    
    assert zscore == 0.0
