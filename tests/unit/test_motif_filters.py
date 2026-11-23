"""`krewlyzer motif` reads with the run's filters and blacklist, everywhere it reads.

It took no filter options and fixed mapq 20, duplicates skipped and no blacklist, where run-all
applies the run's. When the options came, its BAM pre-check kept the fixed values, so its
warning described filters the run was not using.
"""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

from conftest import requires_data
from krewlyzer.cli import app

pytestmark = pytest.mark.unit


class _Stop(Exception):
    pass


@pytest.fixture
def calls(monkeypatch):
    seen = {}

    def check(bam, require_proper_pair, skip_duplicates, mapq):
        seen["check"] = (skip_duplicates, mapq)
        return {"pass_rate": 1.0, "total_sampled": 0}

    def extract(**kwargs):
        seen["extract"] = kwargs
        raise _Stop

    monkeypatch.setattr("krewlyzer.core.bam_utils.check_bam_compatibility", check)
    monkeypatch.setattr("krewlyzer.core.sample_processor.extract_sample", extract)
    return seen


# Asserts the genome's bundled blacklist is applied: needs the data package, not its LFS
# pointers (a pointer used to pass as the blacklist, and this test with it).
@requires_data
def test_the_run_filters_reach_the_check_and_the_extraction(tmp_path, calls):
    bam, ref = tmp_path / "S.bam", tmp_path / "ref.fa"
    bam.write_bytes(b"")
    ref.write_text(">chr1\nACGT\n")
    CliRunner().invoke(
        app,
        [
            "motif",
            "-i",
            str(bam),
            "-r",
            str(ref),
            "-o",
            str(tmp_path / "out"),
            "-s",
            "S",
            "--skip-pon",
            "--skip-target-regions",
            "--mapq",
            "30",
            "--minlen",
            "80",
            "--maxlen",
            "500",
            "--no-skip-duplicates",
        ],
    )
    assert calls["check"] == (False, 30)
    got = calls["extract"]
    assert (got["mapq"], got["minlen"], got["maxlen"], got["skip_duplicates"]) == (
        30,
        80,
        500,
        False,
    )
    assert got["exclude_regions"] is not None, "the bundled blacklist was not applied"
