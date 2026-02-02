"""
Unit tests for AssetManager assay-aware methods.

Tests the new panel-specific asset resolution methods.
"""

import pytest
from pathlib import Path

from krewlyzer.assets import AssetManager


# =============================================================================
# Assay-aware asset resolution tests
# =============================================================================

class TestAssetManagerAssayMethods:
    """Tests for assay-aware asset resolution."""
    
    @pytest.fixture
    def manager(self):
        """Create an AssetManager for GRCh37."""
        return AssetManager("hg19")
    
    def test_list_available_assays(self, manager):
        """List available assays."""
        assays = manager.list_available_assays()
        assert "xs1" in assays
        assert "xs2" in assays
    
    def test_get_gene_bed_xs1(self, manager):
        """Get xs1 gene BED."""
        path = manager.get_gene_bed("xs1")
        assert path.exists()
        assert "xs1" in str(path)
        assert path.suffix == ".gz"
    
    def test_get_gene_bed_xs2(self, manager):
        """Get xs2 gene BED."""
        path = manager.get_gene_bed("xs2")
        assert path.exists()
        assert "xs2" in str(path)
    
    def test_get_gene_bed_unknown_raises(self, manager):
        """Unknown assay raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            manager.get_gene_bed("unknown")
    
    def test_get_target_bed_xs1(self, manager):
        """Get xs1 target BED."""
        path = manager.get_target_bed("xs1")
        assert path.exists()
        assert "xs1" in str(path)
    
    def test_get_target_bed_xs2(self, manager):
        """Get xs2 target BED."""
        path = manager.get_target_bed("xs2")
        assert path.exists()
    
    def test_get_wps_anchors_panel(self, manager):
        """Get panel-specific WPS anchors."""
        path = manager.get_wps_anchors("xs2")
        assert path.exists()
        assert "xs2" in str(path)
    
    def test_get_wps_anchors_genome_wide(self, manager):
        """Get genome-wide WPS anchors when no assay specified."""
        path = manager.get_wps_anchors()
        assert path.exists()
        assert "hg19" in str(path)
    
    def test_get_pon_xs1(self, manager):
        """Get xs1 PON (may be legacy naming)."""
        try:
            path = manager.get_pon("xs1")
            assert path.exists()
        except FileNotFoundError:
            pytest.skip("PON for xs1 not bundled")
    
    def test_get_pon_xs2(self, manager):
        """Get xs2 PON (may be legacy naming)."""
        try:
            path = manager.get_pon("xs2")
            assert path.exists()
        except FileNotFoundError:
            pytest.skip("PON for xs2 not bundled")
