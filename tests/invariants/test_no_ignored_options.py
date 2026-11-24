"""An option a command parses and then ignores is refused, not accepted.

`run-all --chromosomes chr1` and `motif --chromosomes chr1` parsed the value and never read it,
so each wrote a whole-genome sample, and `run-all` wrote it under `{sample}.metadata.parquet`,
downstream's completion marker (invariant #2). Passing the value through was not the fix:
nothing records a chromosome subset, so a one-chromosome sample under a completion marker reads
downstream as a whole one. The option is gone; using it is a usage error that writes nothing.
"""

from __future__ import annotations

import pytest
from rich.text import Text
from typer.testing import CliRunner

from conftest import requires_data
from krewlyzer.cli import app

from .runall import SAMPLE, write_arms
from .synth import BASELINE, make_sample

pytestmark = [pytest.mark.invariant, requires_data]


@pytest.fixture(scope="module")
def sample(tmp_path_factory):
    return make_sample(tmp_path_factory.mktemp("ignored_options_in"), BASELINE)


def _run_all_args(sample, out):
    """The settings `synth_run_all` uses, so the only thing a test adds is what it is about."""
    arms = write_arms(sample.directory / "arms.validated.bed", sample)
    return [
        "run-all",
        "-i",
        sample.bam,
        "-r",
        sample.reference,
        "-o",
        out,
        "-s",
        SAMPLE,
        "-b",
        sample.assets.bins,
        "-a",
        arms,
        "--ocr-file",
        sample.assets.regions,
        "--wps-anchors",
        sample.assets.anchors,
        "--wps-background",
        sample.assets.alu,
        "-x",
        sample.assets.blacklist,
        "--output-format",
        "parquet",
        "--skip-pon",
        "--skip-target-regions",
        "--no-tfbs",
        "--no-atac",
        "-t",
        "1",
    ]


def _motif_args(sample, out):
    return [
        "motif",
        "-i",
        sample.bam,
        "-r",
        sample.reference,
        "-o",
        out,
        "-s",
        SAMPLE,
        "-x",
        sample.assets.blacklist,
        "--skip-pon",
        "--output-format",
        "parquet",
    ]


@pytest.mark.parametrize(
    "build", [_run_all_args, _motif_args], ids=["run-all", "motif"]
)
def test_chromosomes_is_refused_rather_than_ignored(build, sample, tmp_path):
    out = tmp_path / "out"
    args = build(sample, out) + ["--chromosomes", "chr1"]
    result = CliRunner().invoke(app, [str(a) for a in args])
    # Typer styles its error box when GITHUB_ACTIONS is set, which splits the option name
    # with escape codes; the message is compared as text. Rich, not Click: Click is not a
    # krewlyzer dependency, and Typer 0.29 no longer installs it.
    message = Text.from_ansi(result.output).plain
    assert result.exit_code == 2, f"not refused as a usage error: {message}"
    assert "No such option" in message
    assert "--chromosomes" in message
    assert (
        not out.exists()
    ), f"a refused run wrote {sorted(p.name for p in out.iterdir())}"
