"""Normalising FSD twice must give the same answer as normalising once.

`apply_pon_logratio` reads a table, appends log-ratio columns and writes it
back. It found its size bins by taking any header whose text before the first
`-` parsed as an integer -- so `65-69_logR`, a column it had written itself,
came back as bin 65 on the next run. Each pass then appended a fresh set of
log-ratios *and* log-ratios of the previous log-ratios.

Measured on a real 67-bin sample: 69 columns raw, 137 after one pass, 273 after
two, with every `_logR` name repeated. `_write_fsd_output` reads that back with
pandas, which renames the collisions to `_logR.1` / `_logR.2` and writes those
to disk -- which is what the shipped FSD tables contained.

The right answer was in the file three times with no way to say which was
current, and `read_csv` returns the first, which is the oldest.

Nothing in the pipeline marks the work as done, so a rerun into the same
directory, a Nextflow retry or a resumed job all trigger it. That is the case
this file pins.
"""

from __future__ import annotations

import csv

import numpy as np
import pandas as pd
import pytest

BINS = [f"{s}-{s + 4}" for s in range(65, 105, 5)]
#: 9q is deliberately absent from the baseline below.
ARMS = ("1p", "1q", "9q")


@pytest.fixture
def fsd_and_pon(tmp_path):
    # The sample must differ from the baseline, bin by bin, or every log-ratio
    # is log2(1) = 0 and the test cannot tell a bin-matched lookup from a
    # collapsed one. It used to be identical to the baseline and passed only
    # because `size_bin` failed to parse, so all 8 bins were scored against the
    # last one and the ratios varied for the wrong reason.
    raw = pd.DataFrame(
        [
            {
                "region": a,
                **{b: float((i + 1) * 10 + 5 * i) for i, b in enumerate(BINS)},
            }
            for a in ARMS
        ]
    )
    fsd = tmp_path / "s.FSD.tsv"
    raw.to_csv(fsd, sep="\t", index=False)

    rows = [
        {
            "table": "fsd_baseline",
            "arm": arm,
            "size_bin": int(b.split("-")[0]),
            "expected": float((i + 1) * 10),
            "std": 1.0 + i,
        }
        for arm in ("1p", "1q")
        for i, b in enumerate(BINS)
    ]
    pon = tmp_path / "p.pon.parquet"
    pd.concat(
        [
            pd.DataFrame([{"table": "metadata", "schema_version": "1.0"}]),
            pd.DataFrame(rows),
        ],
        ignore_index=True,
    ).to_parquet(pon)
    return fsd, pon


def _normalise(fsd, pon):
    from krewlyzer import _core

    return _core.fsd.apply_pon_logratio(
        str(fsd), str(pon), None, baseline_table="fsd_baseline"
    )


def _header(path):
    with open(path) as handle:
        return next(csv.reader(handle, delimiter="\t"))


def test_normalising_twice_changes_nothing(fsd_and_pon):
    """The property that was missing. Byte-identical, not merely same-shaped."""
    fsd, pon = fsd_and_pon
    _normalise(fsd, pon)
    once = fsd.read_text()
    _normalise(fsd, pon)
    twice = fsd.read_text()
    _normalise(fsd, pon)
    thrice = fsd.read_text()
    assert once == twice == thrice


def test_the_column_count_does_not_grow(fsd_and_pon):
    fsd, pon = fsd_and_pon
    raw_cols = len(_header(fsd))
    counts = []
    for _ in range(3):
        _normalise(fsd, pon)
        counts.append(len(_header(fsd)))
    # raw + one logR per bin + pon_stability, and it stays there
    assert counts == [raw_cols + len(BINS) + 1] * 3, counts


def test_no_duplicated_or_suffixed_column_names(fsd_and_pon):
    """`_logR.1` in a shipped file is the fingerprint of this bug."""
    fsd, pon = fsd_and_pon
    _normalise(fsd, pon)
    _normalise(fsd, pon)
    header = _header(fsd)
    assert len(header) == len(set(header)), "duplicate field names written to disk"
    assert not [h for h in header if h.rstrip("0123456789").endswith(".")]


def test_a_derived_column_is_never_treated_as_a_size_bin(fsd_and_pon):
    """The mechanism, pinned directly: no log-ratio of a log-ratio."""
    fsd, pon = fsd_and_pon
    _normalise(fsd, pon)
    _normalise(fsd, pon)
    header = _header(fsd)
    assert not [h for h in header if h.count("_logR") > 1]
    assert sum(1 for h in header if h.endswith("_logR")) == len(BINS)
    assert header.count("pon_stability") == 1


def test_an_arm_absent_from_the_baseline_gets_no_log_ratio(fsd_and_pon):
    """Not 0.0, which says the sample sits exactly at the healthy baseline."""
    fsd, pon = fsd_and_pon
    _normalise(fsd, pon)
    frame = pd.read_csv(fsd, sep="\t")
    absent = frame[frame.region == "9q"]
    assert np.isnan(absent["65-69_logR"].iloc[0]), "fabricated a log-ratio"
    assert np.isnan(absent["pon_stability"].iloc[0]), "fabricated a stability"


def test_an_arm_in_the_baseline_still_gets_real_numbers(fsd_and_pon):
    """So the fix cannot be 'return NaN always'."""
    fsd, pon = fsd_and_pon
    _normalise(fsd, pon)
    frame = pd.read_csv(fsd, sep="\t")
    present = frame[frame.region == "1p"]
    values = [present[f"{b}_logR"].iloc[0] for b in BINS]
    assert all(np.isfinite(v) for v in values)
    assert len(set(values)) > 1, "log-ratios must vary across bins"

    # Each bin against *its own* baseline, not the last one.
    #
    # `size_bin` is a double in every PON -- the long-format table carries a
    # null there for every other block's rows -- and the reader called
    # `get_int`, which errors on a Double rather than coercing. Every size_bin
    # became 0, so `size >= size_bins.last()` held for every size and all 67
    # bins were scored against the final row. Measured on a shipped PON and a
    # real sample: 41/41 arms matched the last-bin baseline exactly.
    pseudocount = 1.0
    for i, b in enumerate(BINS):
        sample = float((i + 1) * 10 + 5 * i)
        expected = float((i + 1) * 10)
        want = np.log2((sample + pseudocount) / (expected + pseudocount))
        assert values[i] == pytest.approx(want, rel=1e-6), (
            f"bin {b} scored against the wrong baseline: {values[i]:.4f} "
            f"instead of {want:.4f}"
        )
    assert np.isfinite(present["pon_stability"].iloc[0])


def test_process_fsd_run_twice_is_also_stable(fsd_and_pon):
    """The Python entry point, not just the Rust one.

    `process_fsd` reads back and rewrites through `write_table`, which is where
    the `.1` names were introduced; running it twice exercises that half too.
    """
    from krewlyzer.core.fsd_processor import process_fsd

    fsd, pon = fsd_and_pon
    process_fsd(fsd, pon_parquet_path=pon, output_format="tsv")
    once = fsd.read_text()
    process_fsd(fsd, pon_parquet_path=pon, output_format="tsv")
    assert fsd.read_text() == once
