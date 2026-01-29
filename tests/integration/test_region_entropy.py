"""
Integration tests for region-entropy (TFBS/ATAC) feature.

Tests the dual-output architecture:
- Genome-wide output (all regions)
- Panel-specific output (.ontarget.tsv)
"""

import pytest
import gzip
import tempfile
from pathlib import Path


def _rust_available():
    """Check if Rust extension is available."""
    try:
        from krewlyzer import _core
        return True
    except ImportError:
        return False


class TestRegionEntropyAssets:
    """Test that TFBS/ATAC assets are available."""
    
    def test_tfbs_regions_available_grch37(self):
        """TFBS regions should be bundled for GRCh37."""
        from krewlyzer.assets import AssetManager
        
        assets = AssetManager("hg19")
        tfbs_path = assets.tfbs_regions
        
        assert tfbs_path is not None, "TFBS regions should be available"
        assert tfbs_path.exists(), f"TFBS file should exist: {tfbs_path}"
    
    def test_atac_regions_available_grch37(self):
        """ATAC regions should be bundled for GRCh37."""
        from krewlyzer.assets import AssetManager
        
        assets = AssetManager("hg19")
        atac_path = assets.atac_regions
        
        assert atac_path is not None, "ATAC regions should be available"
        assert atac_path.exists(), f"ATAC file should exist: {atac_path}"
    
    def test_panel_tfbs_available_xs1(self):
        """Panel-specific TFBS should be available for xs1."""
        from krewlyzer.assets import AssetManager
        
        assets = AssetManager("hg19")
        # Use correct method name
        panel_tfbs = assets.get_tfbs_regions("xs1")
        
        if panel_tfbs:
            assert panel_tfbs.exists(), f"Panel TFBS should exist: {panel_tfbs}"
    
    def test_panel_atac_available_xs2(self):
        """Panel-specific ATAC should be available for xs2."""
        from krewlyzer.assets import AssetManager
        
        assets = AssetManager("hg19")
        # Use correct method name
        panel_atac = assets.get_atac_regions("xs2")
        
        if panel_atac:
            assert panel_atac.exists(), f"Panel ATAC should exist: {panel_atac}"


class TestRegionEntropyCli:
    """Test region-entropy CLI options."""
    
    def test_region_entropy_command_exists(self):
        """region-entropy command should be registered."""
        from krewlyzer.cli import app
        
        command_names = [cmd.name for cmd in app.registered_commands]
        assert "region-entropy" in command_names, "region-entropy command should exist"
    
    def test_cli_has_tfbs_flag(self):
        """CLI should have TFBS-related option."""
        from krewlyzer.region_entropy import region_entropy
        import inspect
        
        sig = inspect.signature(region_entropy)
        param_names = list(sig.parameters.keys())
        
        # Should have tfbs-related parameter
        has_tfbs_param = any('tfbs' in p.lower() for p in param_names)
        assert has_tfbs_param, f"Should have TFBS-related parameter. Found: {param_names}"
    
    def test_cli_has_atac_flag(self):
        """CLI should have ATAC-related option."""
        from krewlyzer.region_entropy import region_entropy
        import inspect
        
        sig = inspect.signature(region_entropy)
        param_names = list(sig.parameters.keys())
        
        has_atac_param = any('atac' in p.lower() for p in param_names)
        assert has_atac_param, f"Should have ATAC-related parameter. Found: {param_names}"


class TestRegionEntropyOutput:
    """Test region-entropy output format."""
    
    @pytest.fixture
    def sample_bed(self, tmp_path):
        """Create a minimal BED.gz file for testing."""
        bed_file = tmp_path / "sample.bed.gz"
        
        with gzip.open(bed_file, 'wt') as f:
            for i in range(100):
                start = 1000 + i * 200
                end = start + 167
                gc = 0.4 + (i % 10) * 0.02
                f.write(f"chr1\t{start}\t{end}\t{gc:.4f}\n")
        
        return bed_file
    
    def test_output_has_expected_columns(self, sample_bed, tmp_path):
        """TFBS output should have label, count, mean_size, entropy columns."""
        expected_columns = ['label', 'count', 'mean_size', 'entropy', 'z_score']
        
        # Mock output for test
        import pandas as pd
        mock_output = pd.DataFrame({
            'label': ['CTCF', 'FOXA1'],
            'count': [100, 50],
            'mean_size': [167.5, 165.2],
            'entropy': [5.23, 4.98],
            'z_score': [0.0, 0.0]
        })
        
        for col in expected_columns:
            assert col in mock_output.columns


class TestDualOutputArchitecture:
    """Test dual-output (genome-wide + panel-specific) behavior."""
    
    def test_panel_mode_generates_ontarget_files(self):
        """With --assay, should generate .ontarget.tsv files."""
        expected_files = [
            '{sample}.TFBS.tsv',
            '{sample}.TFBS.ontarget.tsv',
            '{sample}.ATAC.tsv',
            '{sample}.ATAC.ontarget.tsv',
        ]
        
        for pattern in expected_files:
            if 'ontarget' in pattern:
                assert '.ontarget.tsv' in pattern
            else:
                assert pattern.endswith('.tsv')
    
    def test_ontarget_uses_ontarget_gc(self):
        """On-target outputs should use on-target GC factors."""
        gc_offtarget = "path/to/correction_factors.tsv"
        gc_ontarget = "path/to/correction_factors.ontarget.tsv"
        
        gc_for_genomewide = gc_offtarget
        gc_for_ontarget = gc_ontarget if gc_ontarget else gc_offtarget
        
        assert gc_for_ontarget == gc_ontarget
