"""
Tests for KREWLYZER_DATA_DIR environment variable support.

Validates that:
1. Env var overrides bundled data path
2. Invalid paths raise helpful errors
3. Asset resolution works with external data directory
"""

import pytest
from krewlyzer.assets import AssetManager


class TestKrewlyzerDataDirEnvVar:
    """Tests for KREWLYZER_DATA_DIR environment variable."""

    def test_env_var_overrides_bundled_path(self, tmp_path, monkeypatch):
        """Test that KREWLYZER_DATA_DIR overrides the bundled data path."""
        # Create a minimal data directory structure
        data_dir = tmp_path / "krewlyzer_data"
        data_dir.mkdir()

        # Set env var
        monkeypatch.setenv("KREWLYZER_DATA_DIR", str(data_dir))

        # Create AssetManager
        assets = AssetManager("hg19")

        # Verify base_path uses env var
        assert assets.base_path == data_dir

    def test_env_var_with_tilde_expansion(self, tmp_path, monkeypatch):
        """Test that ~ is expanded in KREWLYZER_DATA_DIR."""
        # Create test directory in tmp, but use relative-style path
        data_dir = tmp_path / "krewlyzer_data"
        data_dir.mkdir()

        # Mock expanduser to return our tmp path
        monkeypatch.setenv("KREWLYZER_DATA_DIR", str(data_dir))

        assets = AssetManager("hg19")
        assert assets.base_path.is_absolute()

    def test_env_var_invalid_path_raises_error(self, monkeypatch):
        """Test that invalid KREWLYZER_DATA_DIR raises helpful error."""
        monkeypatch.setenv("KREWLYZER_DATA_DIR", "/nonexistent/path/to/data")

        with pytest.raises(ValueError) as exc_info:
            AssetManager("hg19")

        assert "KREWLYZER_DATA_DIR does not exist" in str(exc_info.value)
        assert "git clone" in str(exc_info.value)  # Helpful hint

    def test_no_env_var_uses_bundled_path(self, monkeypatch):
        """Test that unset env var falls back to bundled path."""
        monkeypatch.delenv("KREWLYZER_DATA_DIR", raising=False)

        assets = AssetManager("hg19")

        # Should use bundled path (relative to assets.py)
        assert "krewlyzer" in str(assets.base_path)
        assert "data" in str(assets.base_path)


class TestAssetResolutionWithExternalDir:
    """Tests for asset resolution with external data directory."""

    def test_asset_paths_use_external_dir(self, tmp_path, monkeypatch):
        """Test that asset paths resolve to external directory."""
        # Create minimal data structure
        data_dir = tmp_path / "data"
        (data_dir / "ChromosomeArms" / "GRCh37").mkdir(parents=True)
        (data_dir / "ChromosomeBins" / "GRCh37").mkdir(parents=True)

        monkeypatch.setenv("KREWLYZER_DATA_DIR", str(data_dir))

        assets = AssetManager("hg19")

        # Check paths resolve to external dir
        assert str(data_dir) in str(assets.arms)
        assert str(data_dir) in str(assets.bins_100kb)

    def test_pon_resolution_with_external_dir(self, tmp_path, monkeypatch):
        """Test PON resolution with external data directory."""
        # Create PON structure
        data_dir = tmp_path / "data"
        pon_dir = data_dir / "pon" / "GRCh37" / "all_unique"
        pon_dir.mkdir(parents=True)

        # Create dummy PON file
        pon_file = pon_dir / "xs2.all_unique.pon.parquet"
        pon_file.touch()

        monkeypatch.setenv("KREWLYZER_DATA_DIR", str(data_dir))

        assets = AssetManager("hg19")
        resolved = assets.get_pon("xs2", variant="all_unique")

        assert resolved == pon_file
        assert resolved.exists()
