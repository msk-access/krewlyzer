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


class TestTheDataFolderIsChecked:
    """`bundled_data_dir` checked the folder only when KREWLYZER_DATA_DIR named it.

    The default, the package's own `data/`, was returned unchecked, and the wheel has none
    (PyPI's 100 MB limit). The AssetManager still builds without it -- a command given every
    input explicitly needs nothing bundled -- but a bundled lookup that fails says the folder
    itself is missing, rather than only that one file is.
    """

    def _package_at(self, root, monkeypatch):
        import krewlyzer.assets as assets_module

        package = root / "site-packages" / "krewlyzer"
        package.mkdir(parents=True)
        monkeypatch.setattr(assets_module, "__file__", str(package / "assets.py"))
        monkeypatch.delenv("KREWLYZER_DATA_DIR", raising=False)
        return package

    def test_a_package_without_its_data_folder_still_builds(
        self, tmp_path, monkeypatch
    ):
        package = self._package_at(tmp_path, monkeypatch)
        assets = AssetManager("hg19")
        assert assets.base_path == package / "data"
        assert "pip install" in assets.missing_data_reason()

    def test_a_package_with_its_data_folder_has_no_reason(self, tmp_path, monkeypatch):
        package = self._package_at(tmp_path, monkeypatch)
        (package / "data").mkdir()
        assets = AssetManager("hg19")
        assert assets.base_path == package / "data"
        assert assets.missing_data_reason() is None

    def test_the_resolvers_name_the_missing_folder(self, tmp_path, monkeypatch):
        from krewlyzer.core.asset_resolution import (
            resolve_exclude_regions,
            resolve_target_regions,
        )

        package = self._package_at(tmp_path, monkeypatch)
        assets = AssetManager("hg19")
        for call in (
            lambda: resolve_target_regions(None, "xs1", False, assets),
            lambda: resolve_exclude_regions(None, assets),
        ):
            with pytest.raises(ValueError) as exc_info:
                call()
            assert f"no data folder at {package / 'data'}" in str(exc_info.value)
            assert "pip install" in str(exc_info.value)

    def test_env_var_naming_a_file_is_an_error(self, tmp_path, monkeypatch):
        a_file = tmp_path / "data.tar"
        a_file.write_text("")
        monkeypatch.setenv("KREWLYZER_DATA_DIR", str(a_file))
        with pytest.raises(ValueError, match="KREWLYZER_DATA_DIR"):
            AssetManager("hg19")


class TestLfsPointersAreNotAssets:
    """A clone that skipped `git lfs pull` has every file, as text pointers. An existence
    check passes them: the targets became a panel of nothing, the blacklist excluded nothing.
    """

    def _pointer(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            "version https://git-lfs.github.com/spec/v1\noid sha256:"
            + "0" * 64
            + "\nsize 1024\n"
        )
        return path

    def test_targets_that_are_a_pointer_are_not_found(self, tmp_path, monkeypatch):
        monkeypatch.setenv("KREWLYZER_DATA_DIR", str(tmp_path))
        self._pointer(tmp_path / "targets" / "GRCh37" / "xs1.targets.bed.gz")
        with pytest.raises(FileNotFoundError, match="git lfs pull"):
            AssetManager("hg19").get_target_bed("xs1")

    def test_a_blacklist_that_is_a_pointer_is_an_error(self, tmp_path, monkeypatch):
        from krewlyzer.core.asset_resolution import resolve_exclude_regions

        monkeypatch.setenv("KREWLYZER_DATA_DIR", str(tmp_path))
        self._pointer(
            tmp_path / "exclude-regions" / "GRCh37" / "hg19-blacklist.v2.bed.gz"
        )
        with pytest.raises(ValueError, match="git lfs pull"):
            resolve_exclude_regions(None, AssetManager("hg19"))

    def test_real_files_are_not_pointers(self, tmp_path):
        from krewlyzer.assets import is_lfs_pointer

        real = tmp_path / "x.bed"
        real.write_text("chr1\t0\t10\n")
        assert not is_lfs_pointer(real)
        assert is_lfs_pointer(self._pointer(tmp_path / "p.bed.gz"))
        assert not is_lfs_pointer(tmp_path / "absent.bed")


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


class TestGeneBedRespectsExternalDir:
    """`get_bundled_gene_bed` honours KREWLYZER_DATA_DIR, like AssetManager.

    It did not, and nothing noticed for a long time. Both of its lookups --
    `importlib.resources` and a fallback relative to its own `__file__` -- look
    inside the installed package, and the wheel excludes `data/` to stay under
    PyPI's 100 MB limit. So the layout the README prescribes for pip users, a
    separate data clone plus KREWLYZER_DATA_DIR, resolved AssetManager's assets
    correctly and returned `None` here. `None` is not an error: the caller just
    produced no gene-level output.

    The tests that would have caught it existed and were skipped, because
    `conftest.DATA_AVAILABLE` asked the same dataless installed package.
    """

    def _make_gene_bed(self, tmp_path, assay="xs1", genome="GRCh37"):
        genes = tmp_path / "data" / "genes" / genome
        genes.mkdir(parents=True)
        path = genes / f"{assay}.genes.bed.gz"
        path.write_bytes(b"\x1f\x8b")  # gzip magic; contents are not read here
        return path

    def test_gene_bed_resolves_into_the_external_dir(self, tmp_path, monkeypatch):
        from krewlyzer.core.gene_bed import get_bundled_gene_bed

        expected = self._make_gene_bed(tmp_path)
        monkeypatch.setenv("KREWLYZER_DATA_DIR", str(tmp_path / "data"))

        assert get_bundled_gene_bed("xs1", "GRCh37") == expected

    def test_external_dir_wins_over_the_bundled_copy(self, tmp_path, monkeypatch):
        """The env var is an override, so it must win even when both exist.

        Without this the previous behaviour would still pass the test above on
        a checkout that happens to carry bundled data.
        """
        from krewlyzer.core.gene_bed import get_bundled_gene_bed

        expected = self._make_gene_bed(tmp_path)
        monkeypatch.setenv("KREWLYZER_DATA_DIR", str(tmp_path / "data"))

        resolved = get_bundled_gene_bed("xs1", "GRCh37")
        assert resolved == expected
        assert "site-packages" not in str(resolved)
        assert str(tmp_path) in str(resolved)

    def test_absent_from_the_external_dir_is_still_none(self, tmp_path, monkeypatch):
        """An override that lacks the file must not silently fall back.

        Falling through to the bundled copy would make the override a
        suggestion, and would hide a misconfigured data clone behind whatever
        the package happened to ship.
        """
        from krewlyzer.core.gene_bed import get_bundled_gene_bed

        self._make_gene_bed(tmp_path, assay="xs1")
        monkeypatch.setenv("KREWLYZER_DATA_DIR", str(tmp_path / "data"))

        assert get_bundled_gene_bed("xs2", "GRCh37") is None
