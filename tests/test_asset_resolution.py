"""
Unit tests for core/asset_resolution.py

Tests the centralized asset resolution functions:
- resolve_target_regions()
- resolve_pon_model()
- resolve_exclude_regions()
- get_target_region_count()
"""

import pytest
from pathlib import Path
from unittest.mock import Mock
import logging

from krewlyzer.core.asset_resolution import (
    resolve_target_regions,
    resolve_pon_model,
    resolve_exclude_regions,
    get_target_region_count,
)


@pytest.fixture
def mock_assets():
    """Create a mock AssetManager."""
    assets = Mock()
    assets.get_target_bed = Mock()
    assets.get_pon = Mock()
    assets.missing_data_reason = Mock(return_value=None)  # the data folder is there
    return assets


def _lfs_pointer(path):
    path.write_text(
        "version https://git-lfs.github.com/spec/v1\noid sha256:"
        + "0" * 64
        + "\nsize 9\n"
    )
    return path


@pytest.fixture
def mock_logger():
    """Create a mock logger."""
    return Mock(spec=logging.Logger)


# ═══════════════════════════════════════════════════════════════════════════════
# resolve_target_regions tests
# ═══════════════════════════════════════════════════════════════════════════════


class TestResolveTargetRegions:
    """Tests for resolve_target_regions function."""

    def test_explicit_path_takes_priority(self, mock_assets, mock_logger, tmp_path):
        """Explicit --target-regions path should override everything."""
        explicit_file = tmp_path / "explicit.bed"
        explicit_file.write_text("chr1\t100\t200\n")

        path, source = resolve_target_regions(
            explicit_path=explicit_file,
            assay="xs2",  # Should be ignored
            skip_target_regions=False,
            assets=mock_assets,
            log=mock_logger,
        )

        assert path == explicit_file
        assert source == "explicit"
        mock_assets.get_target_bed.assert_not_called()

    def test_an_explicit_path_that_does_not_exist_is_an_error(
        self, mock_assets, mock_logger, tmp_path
    ):
        """It was warned about and returned, and every caller then ran WGS (`exists()`)."""
        with pytest.raises(ValueError, match="typo.bed"):
            resolve_target_regions(
                tmp_path / "typo.bed", "xs1", False, mock_assets, mock_logger
            )

    def test_an_explicit_path_that_is_an_lfs_pointer_is_an_error(
        self, mock_assets, mock_logger, tmp_path
    ):
        """`-T` pointed into a clone that skipped `git lfs pull`: a panel of nothing."""
        with pytest.raises(ValueError, match="git lfs pull"):
            resolve_target_regions(
                _lfs_pointer(tmp_path / "xs1.targets.bed.gz"),
                "xs1",
                False,
                mock_assets,
                mock_logger,
            )

    def test_skip_flag_returns_none(self, mock_assets, mock_logger):
        """--skip-target-regions should return (None, 'skipped')."""
        path, source = resolve_target_regions(
            explicit_path=None,
            assay="xs2",  # Should be ignored
            skip_target_regions=True,
            assets=mock_assets,
            log=mock_logger,
        )

        assert path is None
        assert source == "skipped"
        mock_assets.get_target_bed.assert_not_called()

    def test_bundled_asset_from_assay(self, mock_assets, mock_logger, tmp_path):
        """--assay should load bundled targets."""
        bundled_file = tmp_path / "xs2.targets.bed.gz"
        bundled_file.write_text("")
        mock_assets.get_target_bed.return_value = bundled_file

        path, source = resolve_target_regions(
            explicit_path=None,
            assay="xs2",
            skip_target_regions=False,
            assets=mock_assets,
            log=mock_logger,
        )

        assert path == bundled_file
        assert source == "bundled"
        mock_assets.get_target_bed.assert_called_once_with("xs2")

    def test_wgs_assay_returns_none(self, mock_assets, mock_logger):
        """--assay wgs should return (None, 'none') since WGS has no targets."""
        mock_assets.get_target_bed.side_effect = FileNotFoundError("No targets for wgs")

        path, source = resolve_target_regions(
            explicit_path=None,
            assay="wgs",
            skip_target_regions=False,
            assets=mock_assets,
            log=mock_logger,
        )

        assert path is None
        assert source == "none"

    @pytest.mark.parametrize("assay", ["WGS", "Wgs"])
    def test_wgs_is_wgs_in_any_case(self, mock_assets, mock_logger, assay):
        mock_assets.get_target_bed.side_effect = FileNotFoundError("No targets")
        assert resolve_target_regions(None, assay, False, mock_assets, mock_logger) == (
            None,
            "none",
        )

    @pytest.mark.parametrize("assay", ["xs1", "xs2", "XS1", "mypanel"])
    def test_a_named_assay_without_bundled_targets_is_an_error(
        self, mock_assets, mock_logger, assay
    ):
        """It was taken for WGS, at debug level: a panel sample measured as WGS, exit 0."""
        mock_assets.get_target_bed.side_effect = FileNotFoundError(
            "tried /data/targets/GRCh37/xs1.targets.bed.gz"
        )
        with pytest.raises(ValueError) as exc:
            resolve_target_regions(None, assay, False, mock_assets, mock_logger)
        message = str(exc.value)
        assert f"'{assay}'" in message
        assert "xs1.targets.bed.gz" in message  # where it looked
        for way_out in (
            "--target-regions",
            "--skip-target-regions",
            "KREWLYZER_DATA_DIR",
        ):
            assert way_out in message

    def test_no_assay_no_skip_returns_none(self, mock_assets, mock_logger):
        """No assay and no skip should return (None, 'none')."""
        path, source = resolve_target_regions(
            explicit_path=None,
            assay=None,
            skip_target_regions=False,
            assets=mock_assets,
            log=mock_logger,
        )

        assert path is None
        assert source == "none"

    def test_explicit_and_skip_raises_error(self, mock_assets, mock_logger, tmp_path):
        """--target-regions and --skip-target-regions together should raise ValueError."""
        explicit_file = tmp_path / "explicit.bed"
        explicit_file.write_text("chr1\t100\t200\n")

        with pytest.raises(ValueError, match="Cannot use both"):
            resolve_target_regions(
                explicit_path=explicit_file,
                assay=None,
                skip_target_regions=True,
                assets=mock_assets,
                log=mock_logger,
            )


# ═══════════════════════════════════════════════════════════════════════════════
# resolve_pon_model tests
# ═══════════════════════════════════════════════════════════════════════════════


class TestResolvePonModel:
    """Tests for resolve_pon_model function."""

    def test_explicit_path_takes_priority(self, mock_assets, mock_logger, tmp_path):
        """Explicit --pon-model path should override everything."""
        explicit_file = tmp_path / "custom.pon.parquet"
        explicit_file.write_bytes(b"")

        path, source = resolve_pon_model(
            explicit_path=explicit_file,
            assay="xs2",  # Should be ignored
            skip_pon=False,
            assets=mock_assets,
            log=mock_logger,
        )

        assert path == explicit_file
        assert source == "explicit"
        mock_assets.get_pon.assert_not_called()

    def test_an_explicit_pon_that_does_not_exist_is_an_error(
        self, mock_assets, mock_logger, tmp_path
    ):
        """It used to warn and be skipped by every caller, giving an unscored run."""
        # ValueError, like every other bad --pon-model argument: each tool prints it as an
        # ERROR line instead of a traceback.
        with pytest.raises(ValueError, match="nope.pon.parquet"):
            resolve_pon_model(
                explicit_path=tmp_path / "nope.pon.parquet",
                assay="xs2",
                skip_pon=False,
                assets=mock_assets,
                log=mock_logger,
            )

    def test_skip_flag_returns_none(self, mock_assets, mock_logger):
        """--skip-pon should return (None, 'skipped')."""
        path, source = resolve_pon_model(
            explicit_path=None,
            assay="xs2",  # Should be ignored
            skip_pon=True,
            assets=mock_assets,
            log=mock_logger,
        )

        assert path is None
        assert source == "skipped"
        mock_assets.get_pon.assert_not_called()

    def test_bundled_asset_from_assay(self, mock_assets, mock_logger, tmp_path):
        """--assay should load bundled PON."""
        bundled_file = tmp_path / "xs2.all_unique.pon.parquet"
        bundled_file.write_bytes(b"")
        mock_assets.get_pon.return_value = bundled_file

        path, source = resolve_pon_model(
            explicit_path=None,
            assay="xs2",
            skip_pon=False,
            assets=mock_assets,
            log=mock_logger,
        )

        assert path == bundled_file
        assert source == "bundled"
        mock_assets.get_pon.assert_called_once_with("xs2", variant="all_unique")

    def test_no_bundled_pon_returns_none(self, mock_assets, mock_logger):
        """Assay without bundled PON should return (None, 'none')."""
        mock_assets.get_pon.side_effect = FileNotFoundError("No PON for assay")

        path, source = resolve_pon_model(
            explicit_path=None,
            assay="custom",
            skip_pon=False,
            assets=mock_assets,
            log=mock_logger,
        )

        assert path is None
        assert source == "none"

    def test_no_assay_no_skip_returns_none(self, mock_assets, mock_logger):
        """No assay and no skip should return (None, 'none')."""
        path, source = resolve_pon_model(
            explicit_path=None,
            assay=None,
            skip_pon=False,
            assets=mock_assets,
            log=mock_logger,
        )

        assert path is None
        assert source == "none"

    def test_explicit_and_skip_raises_error(self, mock_assets, mock_logger, tmp_path):
        """--pon-model and --skip-pon together should raise ValueError."""
        explicit_file = tmp_path / "custom.pon.parquet"
        explicit_file.write_bytes(b"")

        with pytest.raises(ValueError, match="Cannot use both"):
            resolve_pon_model(
                explicit_path=explicit_file,
                assay=None,
                skip_pon=True,
                assets=mock_assets,
                log=mock_logger,
            )

    def test_variant_default_is_all_unique(self, mock_assets, mock_logger, tmp_path):
        """Default variant should be 'all_unique' when not specified."""
        bundled_file = tmp_path / "xs2.all_unique.pon.parquet"
        bundled_file.write_bytes(b"")
        mock_assets.get_pon.return_value = bundled_file

        path, source = resolve_pon_model(
            explicit_path=None,
            assay="xs2",
            skip_pon=False,
            assets=mock_assets,
            log=mock_logger,
            # variant not specified, should default to "all_unique"
        )

        assert path == bundled_file
        mock_assets.get_pon.assert_called_once_with("xs2", variant="all_unique")

    def test_variant_duplex_passed_to_assets(self, mock_assets, mock_logger, tmp_path):
        """--pon-variant duplex should be passed to assets.get_pon()."""
        bundled_file = tmp_path / "xs2.duplex.pon.parquet"
        bundled_file.write_bytes(b"")
        mock_assets.get_pon.return_value = bundled_file

        path, source = resolve_pon_model(
            explicit_path=None,
            assay="xs2",
            skip_pon=False,
            assets=mock_assets,
            variant="duplex",
            log=mock_logger,
        )

        assert path == bundled_file
        assert source == "bundled"
        mock_assets.get_pon.assert_called_once_with("xs2", variant="duplex")

    def test_explicit_path_ignores_variant(self, mock_assets, mock_logger, tmp_path):
        """Explicit --pon-model should ignore --pon-variant."""
        explicit_file = tmp_path / "custom.pon.parquet"
        explicit_file.write_bytes(b"")

        path, source = resolve_pon_model(
            explicit_path=explicit_file,
            assay="xs2",
            skip_pon=False,
            assets=mock_assets,
            variant="duplex",  # Should be ignored
            log=mock_logger,
        )

        assert path == explicit_file
        assert source == "explicit"
        mock_assets.get_pon.assert_not_called()


# ═══════════════════════════════════════════════════════════════════════════════
# get_target_region_count tests
# ═══════════════════════════════════════════════════════════════════════════════


class TestGetTargetRegionCount:
    """Tests for get_target_region_count function."""

    def test_counts_bed_lines(self, tmp_path):
        """Should count non-empty lines in BED file."""
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("chr1\t100\t200\nchr1\t300\t400\nchr2\t500\t600\n")

        count = get_target_region_count(bed_file)
        assert count == 3

    def test_returns_none_for_missing_file(self):
        """Should return None if file doesn't exist."""
        count = get_target_region_count(Path("/nonexistent/file.bed"))
        assert count is None

    def test_returns_none_for_none_path(self):
        """Should return None if path is None."""
        count = get_target_region_count(None)
        assert count is None

    def test_handles_gzip_files(self, tmp_path):
        """Should handle .bed.gz files."""
        import gzip

        bed_file = tmp_path / "test.bed.gz"
        with gzip.open(bed_file, "wt") as f:
            f.write("chr1\t100\t200\nchr1\t300\t400\n")

        count = get_target_region_count(bed_file)
        assert count == 2


# ═══════════════════════════════════════════════════════════════════════════════
# resolve_exclude_regions tests
# ═══════════════════════════════════════════════════════════════════════════════


class TestResolveExcludeRegions:
    """One blacklist rule for extract, run-all and motif. `motif` applied none, so its tables
    counted blacklisted fragments that run-all's excluded."""

    def test_explicit_path_takes_priority(self, mock_assets, mock_logger, tmp_path):
        bed = tmp_path / "mine.bed"
        bed.write_text("chr1\t0\t10\n")
        mock_assets.exclude_regions = tmp_path / "bundled.bed.gz"
        assert resolve_exclude_regions(bed, mock_assets, mock_logger) == (
            bed,
            "explicit",
        )

    def test_an_explicit_file_that_does_not_exist_is_an_error(
        self, mock_assets, mock_logger, tmp_path
    ):
        with pytest.raises(ValueError, match="nope.bed"):
            resolve_exclude_regions(tmp_path / "nope.bed", mock_assets, mock_logger)

    def test_an_explicit_pipe_is_an_error(self, mock_assets, mock_logger, tmp_path):
        """A pipe, as process substitution gives, would be read once for the record's digest
        and the filter would get nothing. Refused here, before anything reads it: a run that
        read it would wait for a writer."""
        import os

        pipe = tmp_path / "regions.pipe"
        os.mkfifo(pipe)
        with pytest.raises(ValueError, match="is not a regular file"):
            resolve_exclude_regions(pipe, mock_assets, mock_logger)

    def test_an_explicit_device_is_an_error(self, mock_assets, mock_logger):
        """`/dev/zero` was read without end for the record's digest, and a terminal waited for
        input (seventh review)."""
        with pytest.raises(ValueError, match="is not a regular file"):
            resolve_exclude_regions(Path("/dev/zero"), mock_assets, mock_logger)

    def test_dev_null_is_taken_as_excluding_nothing(self, mock_assets, mock_logger):
        """The one device taken: it reads as an empty file."""
        import os

        assert resolve_exclude_regions(Path(os.devnull), mock_assets, mock_logger) == (
            Path(os.devnull),
            "explicit",
        )

    def test_the_genomes_bundled_blacklist_by_default(
        self, mock_assets, mock_logger, tmp_path
    ):
        bundled = tmp_path / "hg19-blacklist.v2.bed.gz"
        bundled.write_bytes(b"")
        mock_assets.exclude_regions = bundled
        assert resolve_exclude_regions(None, mock_assets, mock_logger) == (
            bundled,
            "bundled",
        )

    def test_no_bundled_blacklist_is_an_error(self, mock_assets, mock_logger, tmp_path):
        """It was a warning, and the run went on keeping blacklisted fragments, exit 0."""
        mock_assets.exclude_regions = tmp_path / "absent.bed.gz"
        with pytest.raises(ValueError) as exc:
            resolve_exclude_regions(None, mock_assets, mock_logger)
        message = str(exc.value)
        assert "absent.bed.gz" in message
        for way_out in ("--exclude-regions", "KREWLYZER_DATA_DIR"):
            assert way_out in message

    def test_an_explicit_file_that_is_an_lfs_pointer_is_an_error(
        self, mock_assets, mock_logger, tmp_path
    ):
        """`-x` pointed into a clone that skipped `git lfs pull`: it excluded nothing."""
        with pytest.raises(ValueError, match="git lfs pull"):
            resolve_exclude_regions(
                _lfs_pointer(tmp_path / "hg19-blacklist.v2.bed.gz"),
                mock_assets,
                mock_logger,
            )

    def test_an_explicit_file_needs_no_bundled_blacklist(
        self, mock_assets, mock_logger, tmp_path
    ):
        empty = tmp_path / "nothing.bed"
        empty.write_text("")
        mock_assets.exclude_regions = tmp_path / "absent.bed.gz"
        assert resolve_exclude_regions(empty, mock_assets, mock_logger) == (
            empty,
            "explicit",
        )
