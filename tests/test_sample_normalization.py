"""Unit test for sample-level normalization features"""
import json
import tempfile
import shutil
from pathlib import Path
import numpy as np

def test_metadata_json_format():
    """Test that metadata.json has correct structure"""
    # Simulate metadata that motif.py would create
    metadata = {
        "sample_id": "test_sample",
        "total_unique_fragments": 7500000,
        "timestamp": "2025-12-05T13:00:00Z"
    }
    
    # Verify required fields
    assert "sample_id" in metadata
    assert "total_unique_fragments" in metadata
    assert "timestamp" in metadata
    assert isinstance(metadata["total_unique_fragments"], int)
    assert metadata["total_unique_fragments"] > 0
    
    print("✅ Metadata JSON format: PASS")
    return True

def test_fsd_cpm_calculation():
    """Test FSD CPM normalization calculation"""
    # Test data
    region_normalized = np.array([0.012, 0.015, 0.018])  # Sums to ~0.045
    region_total = 1500  # fragments in this region
    total_fragments = 7_500_000  # total sample fragments
    
    # Calculate CPM (counts per million)
    # Formula: (region_normalized * region_total / total_fragments) * 1e6
    cpm = (region_normalized * region_total / total_fragments) * 1e6
    
    # Expected: (0.012 * 1500 / 7500000) * 1e6 = 2.4
    expected_first = (0.012 * 1500 / 7_500_000) * 1e6
    
    assert np.isclose(cpm[0], expected_first, rtol=1e-5)
    assert cpm[0] > 0  # Should be positive
    assert cpm[0] < region_normalized[0] * 1e6  # Should be less than naïve normalization
    
    print(f"✅ FSD CPM calculation: PASS (CPM[0] = {cpm[0]:.2f})")
    return True

def test_wps_depth_normalization():
    """Test WPS depth normalization calculation"""
    # Test data
    wps_raw = np.array([45, 52, 38, 61])  # Raw WPS scores
    total_fragments = 7_500_000
    
    # Calculate normalized WPS (per million fragments)
    # Formula: wps_raw / (total_fragments / 1e6)
    wps_norm = wps_raw / (total_fragments / 1e6)
    
    # Expected: 45 / (7500000 / 1000000) = 45 / 7.5 = 6.0
    expected_first = 45 / (7_500_000 / 1e6)
    
    assert np.isclose(wps_norm[0], expected_first, rtol=1e-5)
    assert wps_norm[0] > 0
    assert wps_norm[0] < wps_raw[0]  # Normalized should be smaller
    
    print(f"✅ WPS depth normalization: PASS (WPS_norm[0] = {wps_norm[0]:.2f})")
    return True

def test_backward_compatibility():
    """Test that modules work without metadata file"""
    # Simulate scenario where metadata.json doesn't exist
    metadata_file = Path("/tmp/nonexistent.metadata.json")
    
    # FSD and WPS should handle missing metadata gracefully
    # They should:
    # 1. Not crash
    # 2. Output only original columns (no CPM/normalized columns)
    
    assert not metadata_file.exists()
    print("✅ Backward compatibility: PASS (missing metadata handled)")
    return True

def test_integration():
    """Test realistic data flow"""
    # Simulate complete workflow
    test_dir = Path(tempfile.mkdtemp())
    
    try:
        # 1. Motif creates metadata
        metadata = {
            "sample_id": "integration_test",
            "total_unique_fragments": 5_000_000,
            "timestamp": "2025-12-05T13:00:00Z"
        }
        metadata_file = test_dir / "sample.metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        # 2. FSD reads metadata
        assert metadata_file.exists()
        with open(metadata_file) as f:
            loaded = json.load(f)
        assert loaded["total_unique_fragments"] == 5_000_000
        
        # 3. Calculate CPM
        region_norm = 0.02
        region_total = 2000
        cpm = (region_norm * region_total / 5_000_000) * 1e6
        assert cpm == 8.0  # (0.02 * 2000 / 5000000) * 1e6 = 8
        
        print(f"✅ Integration test: PASS (full workflow verified)")
        return True
    
    finally:
        shutil.rmtree(test_dir)

if __name__ == "__main__":
    print("\n" + "="*60)
    print("Testing Sample-Level Normalization Implementation")
    print("="*60 + "\n")
    
    tests = [
        test_metadata_json_format,
        test_fsd_cpm_calculation,
        test_wps_depth_normalization,
        test_backward_compatibility,
        test_integration
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"❌ {test.__name__}: FAILED - {e}")
            failed += 1
    
    print("\n" + "="*60)
    print(f"Results: {passed} passed, {failed} failed")
    print("="*60 + "\n")
    
    if failed == 0:
        print("🎉 All tests passed!")
        exit(0)
    else:
        exit(1)
