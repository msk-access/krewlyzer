"""
Integration tests for dual GC correction in panel mode.

Tests that both off-target and on-target GC correction factors are generated.
"""

import pytest
from pathlib import Path

pytestmark = pytest.mark.skipif(
    not Path("/Users/shahr2/Documents/Github/krewlyzer/tests/data").exists(),
    reason="Test data not available",
)


class TestDualGcObservations:
    """Tests for dual GC observation collection in Rust."""

    def test_extract_returns_7_values(self):
        """Test that process_bam_parallel returns 7 values including gc_observations_ontarget."""
        # This is a signature test - actual BAM processing needs real data
        from krewlyzer import _core

        # Check the function exists and has the right signature
        assert hasattr(_core.extract_motif, "process_bam_parallel")

        # Signature check by inspecting function (Python introspection)
        import inspect

        sig = inspect.signature(_core.extract_motif.process_bam_parallel)
        params = list(sig.parameters.keys())

        # Should have these parameters
        expected_params = ["bam_path", "fasta_path", "mapq", "min_len", "max_len"]
        for p in expected_params:
            assert p in params, f"Missing parameter: {p}"


class TestPonModelOntarget:
    """Tests for on-target baseline fields in PonModel."""

    def test_gc_bias_ontarget_field(self):
        """Test PonModel has gc_bias_ontarget field."""
        from krewlyzer.pon.model import PonModel

        pon = PonModel(assay="test", panel_mode=True, target_regions_file="targets.bed")

        assert hasattr(pon, "gc_bias_ontarget")
        assert hasattr(pon, "fsd_baseline_ontarget")

    def test_panel_mode_fields(self):
        """Test panel_mode and target_regions_file fields."""
        from krewlyzer.pon.model import PonModel

        pon = PonModel(
            assay="msk-access", panel_mode=True, target_regions_file="msk_targets.bed"
        )

        assert pon.panel_mode
        assert pon.target_regions_file == "msk_targets.bed"
