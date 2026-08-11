"""The product must say what produced it.

`{sample}.metadata.parquet` is the consumer's completion marker, and it
recorded the genome, the assay and the filters — but not the code version or
the PON. The version was in `{sample}.features.json`, which invariant #2 says
downstream never reads.

That became load-bearing in 0.9.0, which changed what several columns *mean*:
every FSD `_logR` had been scored against the wrong size bin, `pon_stability`
was wrong by a median of 4709%, and WPS gained seven columns. All of those
produced output that was present, finite and plausible, so given a results
directory there was no way to answer *was this made before or after the fix?*
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from krewlyzer.core.output_provenance import (
    PROVENANCE_COLUMNS,
    provenance_fields,
    stamp_metadata,
)

pytestmark = pytest.mark.unit


class _Pon:
    """Only the fields the stamp reads."""

    cohort_digest = "1e7158a60f3e3ca7"
    krewlyzer_version = "0.9.0"


def _marker(directory: Path, sample: str = "S1") -> Path:
    path = directory / f"{sample}.metadata.parquet"
    pd.DataFrame(
        [{"sample_id": sample, "total_fragments": 1000, "genome": "hg19"}]
    ).to_parquet(path)
    return path


def test_a_scored_sample_records_what_it_was_scored_against(tmp_path):
    marker = _marker(tmp_path)
    assert stamp_metadata(
        tmp_path,
        "S1",
        pon=_Pon(),
        pon_path=Path("/data/pon/xs1.all_unique.pon.parquet"),
        output_format="parquet",
    )

    row = pd.read_parquet(marker).iloc[0]
    assert bool(row["pon_applied"]) is True
    assert row["pon_model"] == "xs1.all_unique.pon.parquet"
    assert row["pon_cohort_digest"] == "1e7158a60f3e3ca7"
    assert row["pon_krewlyzer_version"] == "0.9.0"
    assert row["krewlyzer_version"], "the code version must be recorded"


def test_the_pon_path_is_never_recorded_in_full(tmp_path):
    """A shipped table must not carry the operator's directory layout.

    Sample directories are named for patients (invariant #4), and a PON path
    typically sits under a home or scratch directory. The basename answers
    *which model*; the rest answers *whose machine*.
    """
    # A deep absolute path, but deliberately *not* under `/Users` or `/home`:
    # `scripts/check_phi_guard.sh` refuses either in any tracked file, and it
    # is right to — a test that asserts we never publish a home path should not
    # itself contain one. The behaviour under test is basename extraction from
    # a long absolute path, which this exercises identically.
    _marker(tmp_path)
    stamp_metadata(
        tmp_path,
        "S1",
        pon=_Pon(),
        pon_path=Path("/mnt/scratch/analyst-a/private/xs1.duplex.pon.parquet"),
        output_format="parquet",
    )
    recorded = str(
        pd.read_parquet(tmp_path / "S1.metadata.parquet")["pon_model"].iloc[0]
    )
    assert recorded == "xs1.duplex.pon.parquet"
    assert "/" not in recorded and "scratch" not in recorded


def test_an_unscored_run_says_so_rather_than_saying_nothing(tmp_path):
    """`--skip-pon`, no PON, and a refused PON are all legitimately unscored.

    All three produce a directory with no z-scores, and a reader counting
    columns cannot tell them from a run whose scoring silently failed.
    """
    _marker(tmp_path)
    stamp_metadata(tmp_path, "S1", output_format="parquet")
    row = pd.read_parquet(tmp_path / "S1.metadata.parquet").iloc[0]
    assert bool(row["pon_applied"]) is False
    assert row["pon_model"] == ""
    assert row["krewlyzer_version"], "still records which build wrote it"


def test_a_refused_pon_is_not_an_applied_one():
    """The version guard's case: a path, but nothing loaded from it.

    `unified_processor` clears `pon_parquet` when the load is refused, but the
    fields function must not depend on the caller remembering to.
    """
    assert provenance_fields(None, Path("x.pon.parquet"))["pon_applied"] is False
    assert provenance_fields(_Pon(), None)["pon_applied"] is False
    assert provenance_fields(_Pon(), Path("x.pon.parquet"))["pon_applied"] is True


def test_no_metadata_table_is_not_a_failure(tmp_path):
    """A single-feature CLI run writes one table and no completion marker."""
    assert stamp_metadata(tmp_path, "S1") is False


def test_stamping_twice_does_not_accumulate(tmp_path):
    """Re-running a sample into the same directory must not duplicate columns.

    Frame-collision duplicates (`pon_applied.1`) are exactly what
    `no_collided_columns` exists to catch downstream, and a re-run into an
    existing directory is routine — a Nextflow retry, a resumed job.
    """
    _marker(tmp_path)
    for _ in range(3):
        stamp_metadata(
            tmp_path,
            "S1",
            pon=_Pon(),
            pon_path=Path("p.parquet"),
            output_format="parquet",
        )
    columns = list(pd.read_parquet(tmp_path / "S1.metadata.parquet").columns)
    assert len(columns) == len(set(columns)), columns
    for name in PROVENANCE_COLUMNS:
        assert columns.count(name) == 1


def test_the_original_columns_survive(tmp_path):
    """The stamp adds; it must not rewrite the marker's existing content."""
    _marker(tmp_path)
    before = pd.read_parquet(tmp_path / "S1.metadata.parquet").iloc[0].to_dict()
    stamp_metadata(tmp_path, "S1", output_format="parquet")
    after = pd.read_parquet(tmp_path / "S1.metadata.parquet").iloc[0]
    for key, value in before.items():
        assert after[key] == value, f"{key} changed"


def test_run_features_stamps_the_marker():
    """The caller must actually do it, in the one place that knows the answer.

    Written as a source assertion rather than a full pipeline run: this pins
    that `run_features` — not `wrapper.py` — owns the call, which is what gives
    a single-feature CLI run the same provenance as `run-all` (invariant #6).
    """
    import inspect

    from krewlyzer.core import unified_processor

    source = inspect.getsource(unified_processor)
    assert "stamp_metadata(" in source, "run_features no longer records provenance"
    assert "pon=pon_for_zscore" in source, (
        "must record the PON that was *used for scoring*; `pon` alone is still "
        "set under --skip-pon and would claim a scoring that did not happen"
    )
