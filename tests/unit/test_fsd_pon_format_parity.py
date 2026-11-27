"""FSD must be PON-normalised whatever `--output-format` says.

It was not. Under `--output-format tsv` the bins came out as log2 ratios;
under `parquet` they came out as raw counts, same column names, no warning
either way. 0.9.0 makes parquet the Nextflow default, which would have turned
that from affecting nobody into affecting every pipeline run.

Two independent bugs stacked, and each alone was enough:

1. The caller guarded on ``outputs.fsd.exists()`` where ``outputs.fsd`` names
   the ``.tsv``. The Rust writer honours ``--output-format``, so under parquet
   no ``.tsv`` was ever written and the whole post-processing block — PON
   normalisation included — was skipped.

2. ``_write_fsd_output`` read back through ``read_table``, which is
   *parquet-first*: given ``x.FSD.tsv`` it prefers ``x.FSD.parquet``. That
   file is the raw table the single-pass writer emitted *before*
   normalisation, so the log-ratios were computed, logged as
   "41 arms normalized", and then overwritten with the raw counts.

A third came with ``--compress``. The normaliser writes an uncompressed ``FSD.tsv``, which was
removed only under ``parquet``, and ``process_fsd`` looked for its input TSV-first. So the next
run into the folder normalised the previous run's log-ratios a second time instead of the table
its own Rust pass had just written. The inputs below are what the Rust writer leaves for each
format, since it became format-aware: a ``.tsv`` planted under ``parquet`` is a state no run
produces.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from krewlyzer.core.fsd_processor import process_fsd
from krewlyzer.core.output_utils import read_exact_table, table_paths

FIXTURE_PON = (
    Path(__file__).resolve().parents[2] / "tests/data/fixtures/test.pon.parquet"
)

#: Every (--output-format, --compress) a run can ask for; compress changes no Parquet file.
RUNS = [
    ("tsv", False),
    ("tsv", True),
    ("parquet", False),
    ("both", False),
    ("both", True),
]


def _rust_writes(directory: Path, frame: pd.DataFrame, fmt: str, compress: bool):
    """What the format-aware Rust FSD writer leaves for one run: exactly the run's formats."""
    for path in table_paths(directory / "s.FSD", output_format=fmt, compress=compress):
        if path.suffix == ".parquet":
            frame.to_parquet(path, index=False)
        else:
            frame.to_csv(path, sep="\t", index=False)  # gzip inferred from `.gz`


def _shipped(directory: Path, fmt: str, compress: bool) -> pd.DataFrame:
    """The table the run ships: its Parquet when it writes one (downstream reads Parquet)."""
    frame = read_exact_table(
        table_paths(directory / "s.FSD", output_format=fmt, compress=compress)[0]
    )
    assert frame is not None, f"nothing shipped for {fmt}, compress={compress}"
    return frame


def _fsd_files(directory: Path) -> set:
    return {p.name for p in directory.glob("s.FSD*")}


def _fsd_frame() -> pd.DataFrame:
    """Arms and bins the fixture PON has a baseline for."""
    pon = pd.read_parquet(FIXTURE_PON)
    arms = pon[pon.table == "fsd_baseline"].arm.dropna().unique()[:4]
    bins = [f"{s}-{s + 4}" for s in range(65, 400, 5)]
    rng = np.random.default_rng(0)
    data = {"region": list(arms)}
    for i, b in enumerate(bins):
        data[b] = rng.integers(50, 500, size=len(arms)).astype(float) + i
    frame = pd.DataFrame(data)
    frame["total"] = frame[bins].sum(axis=1)
    return frame


def _bins(frame: pd.DataFrame) -> np.ndarray:
    cols = [c for c in frame.columns if "-" in c]
    return frame.sort_values("region")[cols].to_numpy(dtype=float)


def _was_transformed(before: pd.DataFrame, after: pd.DataFrame) -> bool:
    """Did the values change at all?

    Deliberately not "are any negative": a log-ratio is only negative when the
    sample sits below the baseline, so a synthetic frame with counts above
    every PON expectation is entirely positive *and* correctly normalised.
    Sign would have made this test pass for the wrong reason on real data and
    fail for the wrong reason here.
    """
    a, b = _bins(before), _bins(after)
    return a.shape != b.shape or not np.allclose(a, b, equal_nan=True)


@pytest.mark.parametrize("output_format,compress", RUNS)
def test_the_pon_is_applied_under_every_output_format(
    tmp_path, output_format, compress
):
    """The parity that was broken. `both` is included because it exercises
    the same read-back path as `parquet` while still leaving a TSV behind."""
    _rust_writes(tmp_path, _fsd_frame(), output_format, compress)

    process_fsd(
        tmp_path / "s.FSD.tsv",
        pon_parquet_path=FIXTURE_PON,
        output_format=output_format,
        compress=compress,
    )

    assert _was_transformed(
        _fsd_frame(), _shipped(tmp_path, output_format, compress)
    ), (
        f"--output-format {output_format} (compress={compress}) returned the input "
        "unchanged. The PON was silently not applied."
    )


@pytest.mark.parametrize("output_format,compress", RUNS)
@pytest.mark.parametrize("pon", [FIXTURE_PON, None], ids=["pon", "no-pon"])
def test_a_run_leaves_exactly_its_own_formats(tmp_path, output_format, compress, pon):
    """Under `--compress` the normaliser's uncompressed `FSD.tsv` stayed beside the `.tsv.gz`,
    and the next run into the folder read it as its own input."""
    _rust_writes(tmp_path, _fsd_frame(), output_format, compress)
    process_fsd(
        tmp_path / "s.FSD.tsv",
        pon_parquet_path=pon,
        output_format=output_format,
        compress=compress,
    )
    own = table_paths(
        tmp_path / "s.FSD", output_format=output_format, compress=compress
    )
    assert _fsd_files(tmp_path) == {p.name for p in own}


@pytest.mark.parametrize("output_format,compress", RUNS)
def test_a_rerun_normalises_this_run_s_table(tmp_path, output_format, compress):
    """Run 2's Rust pass writes counts seven times run 1's; run 2 must ship what a fresh
    folder ships for those counts, not run 1's log-ratios normalised again."""
    reused, fresh = tmp_path / "reused", tmp_path / "fresh"
    reused.mkdir()
    fresh.mkdir()
    bins = [c for c in _fsd_frame().columns if "-" in c or c == "total"]
    seven = _fsd_frame()
    seven[bins] = seven[bins] * 7.0
    kw = dict(
        pon_parquet_path=FIXTURE_PON, output_format=output_format, compress=compress
    )

    _rust_writes(reused, _fsd_frame(), output_format, compress)
    process_fsd(reused / "s.FSD.tsv", **kw)
    _rust_writes(reused, seven, output_format, compress)
    process_fsd(reused / "s.FSD.tsv", **kw)

    _rust_writes(fresh, seven, output_format, compress)
    process_fsd(fresh / "s.FSD.tsv", **kw)

    pd.testing.assert_frame_equal(
        _shipped(reused, output_format, compress),
        _shipped(fresh, output_format, compress),
    )


def test_a_parquet_only_input_is_still_normalised(tmp_path):
    """The caller names the .tsv; under parquet only the .parquet exists.

    Guarding on `.exists()` against the .tsv name skipped the whole block.
    """
    _fsd_frame().to_parquet(tmp_path / "s.FSD.parquet")
    process_fsd(
        tmp_path / "s.FSD.tsv",
        pon_parquet_path=FIXTURE_PON,
        output_format="parquet",
        compress=False,
    )
    assert _was_transformed(_fsd_frame(), pd.read_parquet(tmp_path / "s.FSD.parquet"))


def test_a_stale_raw_parquet_never_wins_over_the_normalised_result(tmp_path):
    """Bug 2 in isolation.

    A raw `.parquet` sits beside the `.tsv` the normaliser writes. The
    parquet-first reader preferred it, so the result was the input.
    """
    frame = _fsd_frame()
    frame.to_csv(tmp_path / "s.FSD.tsv", sep="\t", index=False)
    frame.to_parquet(tmp_path / "s.FSD.parquet")  # the stale raw sibling

    process_fsd(
        tmp_path / "s.FSD.tsv",
        pon_parquet_path=FIXTURE_PON,
        output_format="parquet",
        compress=False,
    )
    assert _was_transformed(_fsd_frame(), pd.read_parquet(tmp_path / "s.FSD.parquet"))


@pytest.mark.parametrize("output_format,compress", RUNS)
def test_without_a_pon_rust_s_tables_ship_as_written(
    tmp_path, output_format, compress, caplog
):
    """No PON must not mean a silent transform -- nor a false warning. Under `parquet` and
    `--compress` the run has no plain `FSD.tsv`, and every such run without a PON logged
    "FSD output not found for format conversion" about a table that was there."""
    _rust_writes(tmp_path, _fsd_frame(), output_format, compress)
    with caplog.at_level(logging.WARNING):
        process_fsd(
            tmp_path / "s.FSD.tsv",
            pon_parquet_path=None,
            output_format=output_format,
            compress=compress,
        )
    assert not _was_transformed(
        _fsd_frame(), _shipped(tmp_path, output_format, compress)
    )
    assert not [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]


@pytest.mark.parametrize("output_format", ["tsv", "both"])
def test_without_a_pon_an_earlier_run_s_tsv_is_never_shipped(tmp_path, output_format):
    """A compressed run with a PON used to leave its normalised `FSD.tsv` behind. A rerun
    without one (`--skip-pon`) read that file back and wrote it into its own `.tsv.gz`: the
    earlier run's log-ratios, shipped as this run's raw counts."""
    earlier = _fsd_frame()
    earlier["65-69_logR"] = 9.0
    earlier.to_csv(tmp_path / "s.FSD.tsv", sep="\t", index=False)
    _rust_writes(tmp_path, _fsd_frame(), output_format, True)

    process_fsd(
        tmp_path / "s.FSD.tsv",
        pon_parquet_path=None,
        output_format=output_format,
        compress=True,
    )
    shipped = read_exact_table(tmp_path / "s.FSD.tsv.gz")
    pd.testing.assert_frame_equal(shipped, _fsd_frame())
