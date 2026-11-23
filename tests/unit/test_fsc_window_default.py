"""`krewlyzer fsc`/`fsr` choose the window by mode, as run-all does.

Run-all aggregated a panel's FSC/FSR one 100 kb bin per row and WGS 50 bins (5 Mb) per row; the
single tools always defaulted to 50. So a panel sample's single-tool tables held 1/50th of
run-all's rows, and tool-level Nextflow mode published those (invariant #6).
"""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

from krewlyzer.cli import app
from krewlyzer.core.fsc_processor import (
    PANEL_CONTINUE_N,
    WGS_CONTINUE_N,
    default_continue_n,
)

pytestmark = pytest.mark.unit


class _NoOutputs:
    def __getattr__(self, name):
        return None


@pytest.fixture
def given(monkeypatch):
    """The `fsc_continue_n` each call passed to `run_features`."""
    seen = []

    def run_features(**kwargs):
        seen.append(kwargs["fsc_continue_n"])
        return _NoOutputs()

    monkeypatch.setattr("krewlyzer.core.unified_processor.run_features", run_features)
    return seen


def _invoke(tool, tmp_path, *extra):
    bed = tmp_path / "S.bed.gz"
    bed.write_bytes(b"")
    bins = tmp_path / "bins.bed"
    bins.write_text("chr1\t0\t100000\n")
    result = CliRunner().invoke(
        app,
        [
            tool,
            "-i",
            str(bed),
            "-o",
            str(tmp_path / "out"),
            "-s",
            "S",
            "-b",
            str(bins),
            "--skip-pon",
            *extra,
        ],
    )
    assert result.exit_code == 0, result.output


def test_the_rule_is_run_alls():
    assert (default_continue_n(True), default_continue_n(False)) == (1, 50)
    assert (PANEL_CONTINUE_N, WGS_CONTINUE_N) == (1, 50)


@pytest.mark.parametrize("tool", ["fsc", "fsr"])
def test_a_panel_defaults_to_one_bin_per_row_and_wgs_to_fifty(tool, tmp_path, given):
    targets = tmp_path / "targets.bed"
    targets.write_text("chr1\t100\t200\n")
    _invoke(tool, tmp_path, "-T", str(targets))
    _invoke(tool, tmp_path, "--skip-target-regions")
    assert given == [1, 50]


@pytest.mark.parametrize("tool", ["fsc", "fsr"])
def test_an_explicit_window_count_is_kept(tool, tmp_path, given):
    targets = tmp_path / "targets.bed"
    targets.write_text("chr1\t100\t200\n")
    _invoke(tool, tmp_path, "-T", str(targets), "--continue-n", "7")
    assert given == [7]
