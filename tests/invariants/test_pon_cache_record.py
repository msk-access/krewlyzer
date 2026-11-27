"""A kept `build-pon` cache carries each donor's record, as run-all's and `extract`'s folders do.

`build-pon --keep-sample-outputs` ends by printing the `--from-outputs` command that rebuilds the
PON from the cache without re-reading a BAM, and that rebuild refused every donor as written before
the GC-lookup fix: it reads each donor's GC provenance from its record (`{sample}.metadata`), and
the PON build wrote none. run-all and `extract` write it
through `write_extraction_outputs`; the build's per-donor path (`process_sample`) wrote the factor
tables itself and skipped the record, so `run_features` had nothing to stamp. In a cache that
already held an `extract` record it stamped that one instead: another extraction's word on whether
a factor table was written, beside the table this build had just fitted.
"""

from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest

from conftest import requires_data

from .pon_donor import build_donor, dense_panel, extract, factor_tables
from .runall import SAMPLE, write_arms

pytestmark = [pytest.mark.invariant, requires_data]

#: Record columns that differ between two extractions of the same input by design.
RUN_SPECIFIC = {"timestamp", "extraction_time_seconds"}


@pytest.fixture(scope="module")
def panel(tmp_path_factory):
    """The dense panel donor, with two target sets: targets on chr1 only, so extraction fits both
    factor tables, and the dense ones, which leave too few off-target fragments to fit one. And the
    arms file FSD reads: FSD on turns on the PON's GC gate."""
    root = tmp_path_factory.mktemp("pon_cache_record")
    donor, unfittable = dense_panel(root)
    fitted = root / "chr1.bed"
    fitted.write_text("chr1\t0\t300000\n")
    return SimpleNamespace(
        donor=donor,
        fitted=fitted,
        unfittable=unfittable,
        arms=write_arms(donor.directory / "arms.validated.bed", donor),
    )


def _record(folder):
    path = folder / f"{SAMPLE}.metadata.tsv"
    assert path.exists(), f"no record in {folder.name}"
    return pd.read_csv(path, sep="\t").iloc[0]


def test_the_pon_build_writes_the_record_extract_writes(panel, tmp_path):
    """One donor through `extract` and through the build's per-donor path, with the same settings:
    every column `extract` records, the same, and the same factor tables beside it."""
    extracted, built = tmp_path / "extract", tmp_path / "build"
    result = extract(panel.donor, extracted, panel.fitted, output_format="tsv")
    assert result.exit_code == 0, result.output
    build_donor(panel.donor, built, panel.fitted)

    want, got = _record(extracted), _record(built)
    assert want[
        "gc_correction_computed"
    ], "extract fitted no table: the comparison is moot"
    differ = {
        column: (want[column], got.get(column))
        for column in want.index
        if column not in RUN_SPECIFIC and str(want[column]) != str(got.get(column))
    }
    assert not differ, f"the build's record differs from extract's: {differ}"
    populations = ("off_target", "on_target")
    assert factor_tables(extracted, *populations) == factor_tables(built, *populations)
    assert len(factor_tables(built, *populations)) == 2, "a population fitted no table"


def test_from_outputs_gc_check_accepts_the_donor_the_build_accepted(panel, tmp_path):
    """The build accepted the donor on GC grounds; `--from-outputs`' GC check of its cache
    (`gc_reason`) must too. It found no record and refused the donor as written before the
    GC-lookup fix."""
    from krewlyzer.pon import from_outputs

    cache = tmp_path / SAMPLE
    build_donor(panel.donor, cache, panel.fitted, arms=panel.arms)  # raises if refused

    assert from_outputs.gc_reason(cache, SAMPLE) is None


def test_from_outputs_gc_check_refuses_the_donor_for_the_build_s_reason(
    panel, tmp_path
):
    """A donor the build refuses on GC: `--from-outputs`' GC check of the cache it leaves refuses
    it for that reason, not as a build from before the GC-lookup fix."""
    from krewlyzer.pon import from_outputs

    cache = tmp_path / SAMPLE
    with pytest.raises(RuntimeError, match="no off-target GC factor table"):
        build_donor(panel.donor, cache, panel.unfittable, arms=panel.arms)

    reason = from_outputs.gc_reason(cache, SAMPLE)
    assert reason is not None and reason.startswith(
        "not GC-corrected, off-target"
    ), reason


def test_an_earlier_record_in_the_folder_is_replaced(panel, tmp_path):
    """The per-donor path into a folder holding an `extract --no-gc-correct` record of the donor.
    It fits both tables, and that record, saying none was written, disowned them: the GC gate
    refused the donor for having no off-target table. The record must be this run's. (`build-pon`
    itself refuses a used `--keep-sample-outputs` directory; the per-donor path must not depend on
    that.)"""
    cache = tmp_path / SAMPLE
    result = extract(
        panel.donor, cache, panel.fitted, "--no-gc-correct", output_format="tsv"
    )
    assert result.exit_code == 0, result.output
    assert not _record(cache)["gc_correction_computed"]

    build_donor(panel.donor, cache, panel.fitted, arms=panel.arms)  # raises if refused

    assert _record(cache)["gc_correction_computed"]
