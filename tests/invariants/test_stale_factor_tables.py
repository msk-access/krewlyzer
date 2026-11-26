"""A GC factor table must never outlive the BED it was fitted on.

The resolver finds a sample's table in any format, Parquet first. Nothing removed the tables an
earlier run left in the same directory, so a rerun reused the old run's weights whenever it wrote
no table of its own (no fittable length bin, no observations, a failed fit) or wrote a different
format (an old `.parquet` beside a new `.tsv` wins). The metadata then recorded the correction as
applied. Review of PR #89.

Only extraction writes a table, so only extraction removes one. A table in a feature step's output
folder is never read (`test_factor_table_writes.py`), so that step leaves it alone.
"""

from pathlib import Path

import pytest

from krewlyzer.core.gc_assets import FACTOR_KINDS, POPULATIONS, factor_table_base
from krewlyzer.core.sample_processor import extract_sample

from .synth import make_sample

pytestmark = [pytest.mark.invariant, pytest.mark.rust]


def _plant_stale_tables(directory: Path, sample: str):
    directory.mkdir(parents=True, exist_ok=True)
    planted = []
    for population in POPULATIONS:
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
