"""A GC factor table must never outlive the BED it was fitted on.

The resolver finds a sample's table in any format, Parquet first. Nothing removed the tables an
earlier run left in the same directory, so a rerun reused the old run's weights whenever it wrote
no table of its own (no fittable length bin, no observations, a failed fit) or wrote a different
format (an old `.parquet` beside a new `.tsv` wins). The metadata then recorded the correction as
applied. Review of PR #89.
"""

from pathlib import Path

import pytest

from krewlyzer.core.gc_assets import FACTOR_KINDS, POPULATIONS, factor_table_base
from krewlyzer.core.sample_processor import extract_sample
from krewlyzer.core.unified_processor import run_features

from .synth import make_sample

pytestmark = [pytest.mark.invariant, pytest.mark.rust]


def _plant_stale_tables(directory: Path, sample: str, populations=POPULATIONS):
    directory.mkdir(parents=True, exist_ok=True)
    planted = []
    for population in populations:
        base = factor_table_base(directory, sample, population)
        for _, ext in FACTOR_KINDS:
            path = Path(f"{base}{ext}")
            path.write_text("an earlier run's table\n")
            planted.append(path)
    return planted


def test_a_fresh_extraction_removes_the_earlier_runs_factor_tables(tmp_path):
    sample = make_sample(tmp_path / "s")
    out = tmp_path / "out"
    planted = _plant_stale_tables(out, "s")
    extract_sample(
        sample.bam, sample.reference, output_dir=out, sample_name="s", threads=1
    )
    assert (out / "s.bed.gz").exists()
    assert [p.name for p in planted if p.exists()] == []


def test_an_extraction_that_writes_no_bed_leaves_the_tables_alone(tmp_path):
    """`krewlyzer motif` extracts without writing a BED: the tables still belong to the BED."""
    sample = make_sample(tmp_path / "s")
    out = tmp_path / "out"
    planted = _plant_stale_tables(out, "s")
    extract_sample(sample.bam, sample.reference, write_bed=False, sample_name="s")
    assert all(p.exists() for p in planted)


def test_refitting_from_a_bed_removes_the_earlier_table_in_the_output_directory(
    tmp_path,
):
    """With no table beside the BED, run_features refits into the output directory; an older
    table there, in any format, must not survive to be picked up by the next lookup."""
    sample = make_sample(tmp_path / "s")
    bed_dir = tmp_path / "bed"
    extract_sample(
        sample.bam, sample.reference, output_dir=bed_dir, sample_name="s", threads=1
    )
    out = tmp_path / "features"
    planted = _plant_stale_tables(out, "s", populations=("off_target",))
    run_features(
        bed_path=bed_dir / "s.bed.gz",
        output_dir=out,
        sample_name="s",
        fsc_bins=sample.assets.bins,
        enable_fsc=True,
        gc_correct=True,
        threads=1,
    )
    assert [
        p.name
        for p in planted
        if p.exists() and p.read_bytes().startswith(b"an earlier")
    ] == []
