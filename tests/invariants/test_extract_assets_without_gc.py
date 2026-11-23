"""`extract --no-gc-correct` still resolves the assay's targets and the bundled blacklist.

`extract` built its AssetManager only under `--gc-correct`, and target and blacklist
resolution both run through it. So `extract --no-gc-correct -A xs1` ran in WGS mode with no
targets and no blacklist, exit 0, while run-all (which always builds the AssetManager) ran
the same sample in panel mode with the blacklist applied. Single tool and run-all must match
(AGENTS.md invariant #6); GC correction is the only thing `--no-gc-correct` turns off.

The resolved paths are read where both entry points hand them to the one extraction routine,
`core.sample_processor.extract_sample`, and the panel mode from the sample's metadata row.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

import krewlyzer.core.sample_processor as sample_processor
from conftest import requires_data
from krewlyzer.assets import AssetManager
from krewlyzer.cli import app

from .runall import SAMPLE, record, synth_run_all
from .synth import BASELINE, make_sample

pytestmark = [pytest.mark.invariant, requires_data]

ASSAY = "xs1"


@pytest.fixture(scope="module")
def sample(tmp_path_factory):
    return make_sample(tmp_path_factory.mktemp("extract_assets_in"), BASELINE)


@pytest.fixture
def extraction_inputs(monkeypatch):
    """The target and exclude regions each call of `extract_sample` received."""
    seen = []
    real = sample_processor.extract_sample

    def spy(*args, **kwargs):
        seen.append(
            {
                "target_regions": kwargs.get("target_regions"),
                "exclude_regions": kwargs.get("exclude_regions"),
            }
        )
        return real(*args, **kwargs)

    monkeypatch.setattr(sample_processor, "extract_sample", spy)
    return seen


def _extract(sample, out: Path, *extra):
    args = [
        "extract",
        "-i",
        sample.bam,
        "--reference",
        sample.reference,
        "--output",
        out,
        "--sample-name",
        SAMPLE,
        "--output-format",
        "parquet",
        "--threads",
        1,
        *extra,
    ]
    return CliRunner().invoke(app, [str(a) for a in args])


@pytest.mark.parametrize("gc_flag", ["--no-gc-correct", "--gc-correct"])
def test_extract_applies_the_assays_targets_and_blacklist(
    sample, tmp_path, extraction_inputs, gc_flag
):
    out = tmp_path / "out"
    result = _extract(sample, out, "-A", ASSAY, gc_flag)
    assert result.exit_code == 0, result.output

    assets = AssetManager("hg19")
    assert extraction_inputs == [
        {
            "target_regions": assets.get_target_bed(ASSAY),
            "exclude_regions": assets.exclude_regions,
        }
    ], f"extract {gc_flag} did not resolve the assay's targets and the blacklist"

    row = record(out)
    assert bool(row["panel_mode"]) is True
    assert row["target_regions"] == str(assets.get_target_bed(ASSAY))


def test_extract_without_gc_matches_run_all(sample, tmp_path, extraction_inputs):
    """The same sample and assay: both entry points extract with the same regions.

    Each resolves its own default blacklist here: the harness gives run-all the synthetic
    one, and this compares the resolution, not the harness."""
    code = synth_run_all(
        sample,
        tmp_path / "run_all",
        assay=ASSAY,
        skip_target_regions=False,
        exclude_regions=None,
    )
    assert code is None
    result = _extract(sample, tmp_path / "extract", "-A", ASSAY, "--no-gc-correct")
    assert result.exit_code == 0, result.output

    run_all_inputs, extract_inputs = extraction_inputs
    assert extract_inputs == run_all_inputs
    assert all(v is not None for v in extract_inputs.values())


@pytest.mark.parametrize("gc_flag", ["--no-gc-correct", "--gc-correct"])
def test_extract_refuses_an_unknown_genome(sample, tmp_path, gc_flag):
    """run-all exits on a genome it has no assets for. `extract` disabled GC correction and
    carried on -- with no targets and no blacklist either, which nothing on screen said.
    """
    out = tmp_path / "out"
    result = _extract(sample, out, "-A", ASSAY, "-G", "hg99", gc_flag)
    assert result.exit_code == 1, result.output
    assert "Unsupported genome" in result.output and "Traceback" not in result.output
    assert not (out / f"{SAMPLE}.bed.gz").exists()
