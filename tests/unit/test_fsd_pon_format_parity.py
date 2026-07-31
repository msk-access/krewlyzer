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
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from krewlyzer.core.fsd_processor import process_fsd

FIXTURE_PON = (
    Path(__file__).resolve().parents[2] / "tests/data/fixtures/test.pon.parquet"
)


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


@pytest.mark.parametrize("output_format", ["tsv", "parquet", "both"])
def test_the_pon_is_applied_under_every_output_format(tmp_path, output_format):
    """The parity that was broken. `both` is included because it exercises
    the same read-back path as `parquet` while still leaving a TSV behind."""
    base = tmp_path / "s.FSD.tsv"
    _fsd_frame().to_csv(base, sep="\t", index=False)

    process_fsd(base, pon_parquet_path=FIXTURE_PON, output_format=output_format)

    written = tmp_path / ("s.FSD.parquet" if output_format != "tsv" else "s.FSD.tsv")
    assert written.exists(), f"nothing written for {output_format}"
    frame = (
        pd.read_parquet(written)
        if written.suffix == ".parquet"
        else pd.read_csv(written, sep="\t", comment="#")
    )
    assert _was_transformed(_fsd_frame(), frame), (
        f"--output-format {output_format} returned the input unchanged. "
        "The PON was silently not applied."
    )


def test_a_parquet_only_input_is_still_normalised(tmp_path):
    """The caller names the .tsv; under parquet only the .parquet exists.

    Guarding on `.exists()` against the .tsv name skipped the whole block.
    """
    _fsd_frame().to_parquet(tmp_path / "s.FSD.parquet")
    process_fsd(
        tmp_path / "s.FSD.tsv", pon_parquet_path=FIXTURE_PON, output_format="parquet"
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
        tmp_path / "s.FSD.tsv", pon_parquet_path=FIXTURE_PON, output_format="parquet"
    )
    assert _was_transformed(_fsd_frame(), pd.read_parquet(tmp_path / "s.FSD.parquet"))


def test_without_a_pon_the_counts_are_left_alone(tmp_path):
    """The other direction: no PON must not mean a silent transform."""
    base = tmp_path / "s.FSD.tsv"
    _fsd_frame().to_csv(base, sep="\t", index=False)
    process_fsd(base, pon_parquet_path=None, output_format="parquet")
    assert not _was_transformed(
        _fsd_frame(), pd.read_parquet(tmp_path / "s.FSD.parquet")
    )
