"""What each optional stage did, recorded in the product.

run-all carries on when an optional stage fails, by design: one stage must not cost a sample its
other features. The outcome was known only to the log, so "not requested" (no variants file, no
assay) and "requested and failed" looked identical to a reader counting tables -- the
`pon_applied` problem again. `{sample}.metadata` now says which, one writer for every stage.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from krewlyzer.core.output_provenance import Outcome, record_stage_outcomes
from krewlyzer.validate.contract import (
    CONTRACT,
    NOT_CONSUMED,
    STAGE_STATUSES,
    STAGE_TABLES,
    STAGES,
    Kind,
    Vary,
)

pytestmark = pytest.mark.unit

#: Shaped like a real sample directory, which is named for the patient (invariant #4).
PATIENT_DIR = "/mnt/results/P-0000000/C-000000-L000"


def _marker(directory: Path, fmt: str = "parquet", sample: str = "S1") -> None:
    frame = pd.DataFrame([{"sample_id": sample, "total_fragments": 1000}])
    if fmt in ("parquet", "both"):
        frame.to_parquet(directory / f"{sample}.metadata.parquet")
    if fmt in ("tsv", "both"):
        frame.to_csv(directory / f"{sample}.metadata.tsv", sep="\t", index=False)


def _row(directory: Path, fmt: str = "parquet", sample: str = "S1") -> dict:
    if fmt == "parquet":
        return (
            pd.read_parquet(directory / f"{sample}.metadata.parquet").iloc[0].to_dict()
        )
    frame = pd.read_csv(
        directory / f"{sample}.metadata.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    return frame.iloc[0].to_dict()


@pytest.mark.parametrize("fmt", ["parquet", "tsv"])
def test_each_stage_records_its_status_and_error_class(tmp_path, fmt):
    _marker(tmp_path, fmt)
    assert record_stage_outcomes(
        tmp_path,
        "S1",
        {
            "mfsd": Outcome.failed(RuntimeError("boom")),
            "uxm": Outcome("not_requested"),
            "region_mds": Outcome("ran"),
        },
        output_format=fmt,
        compress=False,
    )
    row = _row(tmp_path, fmt)
    assert row["mfsd_status"] == "failed" and row["mfsd_error"] == "RuntimeError"
    assert row["uxm_status"] == "not_requested" and row["uxm_error"] == ""
    assert row["region_mds_status"] == "ran" and row["region_mds_error"] == ""
    assert row["total_fragments"] in (1000, "1000"), "existing columns survive"


def test_two_outcomes_give_two_records(tmp_path):
    """Invariant #1, for the record itself: a failed run and a clean one must not read alike."""
    _marker(tmp_path)
    record_stage_outcomes(
        tmp_path,
        "S1",
        {"mfsd": Outcome.failed(OSError("x"))},
        output_format="parquet",
        compress=False,
    )
    failed = _row(tmp_path)
    record_stage_outcomes(
        tmp_path,
        "S1",
        {"mfsd": Outcome("ran")},
        output_format="parquet",
        compress=False,
    )
    ran = _row(tmp_path)
    assert (failed["mfsd_status"], failed["mfsd_error"]) == ("failed", "OSError")
    assert (ran["mfsd_status"], ran["mfsd_error"]) == ("ran", "")


def test_a_recorded_error_is_the_class_and_never_the_message():
    """Messages carry input paths, and sample directories are named for the patient."""
    outcome = Outcome.failed(
        FileNotFoundError(f"BAM file not found: {PATIENT_DIR}/sample.bam")
    )
    assert outcome.error_class == "FileNotFoundError"
    assert "/" not in outcome.error_class and "P-0000000" not in outcome.error_class
    with pytest.raises(ValueError, match="class name"):
        Outcome("failed", f"{PATIENT_DIR}/sample.bam")


def test_an_outcome_is_one_of_the_declared_statuses():
    for status in STAGE_STATUSES:
        if status != "failed":
            assert Outcome(status).status == status
    with pytest.raises(ValueError):
        Outcome("reused")
    with pytest.raises(ValueError, match="failed"):
        Outcome("failed")  # a failure without a class says nothing
    with pytest.raises(ValueError):
        Outcome("ran", "RuntimeError")  # only a failure carries a class


def test_an_unknown_stage_is_refused(tmp_path):
    _marker(tmp_path)
    with pytest.raises(ValueError, match="unknown stage"):
        record_stage_outcomes(
            tmp_path,
            "S1",
            {"mfds": Outcome("ran")},
            output_format="parquet",
            compress=False,
        )


def test_a_single_tool_never_writes_not_requested(tmp_path):
    """A single tool knows what it ran, not what run-all was asked for: `krewlyzer fsc` without
    `--assay` must not record gene FSC as not requested over run-all's `ran`."""
    _marker(tmp_path)
    record_stage_outcomes(
        tmp_path,
        "S1",
        {"gene_fsc": Outcome("ran"), "tfbs": Outcome("ran")},
        output_format="parquet",
        compress=False,
    )
    record_stage_outcomes(
        tmp_path,
        "S1",
        {"gene_fsc": Outcome("not_requested"), "tfbs": Outcome.failed(OSError("x"))},
        output_format="parquet",
        owns_record=False,
        compress=False,
    )
    row = _row(tmp_path)
    assert row["gene_fsc_status"] == "ran"
    # What it did attempt, it records: it rewrote those tables, so the last writer wins.
    assert (row["tfbs_status"], row["tfbs_error"]) == ("failed", "OSError")


def test_no_metadata_table_is_not_a_failure(tmp_path):
    """A single-feature run into an empty directory writes no completion marker."""
    assert (
        record_stage_outcomes(
            tmp_path, "S1", {"uxm": Outcome("ran")}, output_format="tsv", compress=False
        )
        is False
    )


def test_a_run_reads_the_record_of_its_own_format(tmp_path):
    """`both` reads the Parquet copy and writes both: resolving `.tsv` first and then reading
    Parquet first wrote a stale row back (the stamp's defect, fixed for this writer too).
    """
    pd.DataFrame([{"sample_id": "S1", "total_fragments": 1}]).to_csv(
        tmp_path / "S1.metadata.tsv", sep="\t", index=False
    )
    pd.DataFrame([{"sample_id": "S1", "total_fragments": 2}]).to_parquet(
        tmp_path / "S1.metadata.parquet"
    )
    record_stage_outcomes(
        tmp_path, "S1", {"uxm": Outcome("ran")}, output_format="both", compress=False
    )
    assert _row(tmp_path, "tsv")["total_fragments"] == "2"
    assert _row(tmp_path, "parquet")["uxm_status"] == "ran"


def test_every_stage_table_is_declared_somewhere():
    """A stage's tables decide its severity: ERROR if any is consumed. A suffix that is in
    neither list would silently make a failure look harmless."""
    known = {rule.suffix for rule in CONTRACT} | set(NOT_CONSUMED)
    for stage, suffixes in STAGE_TABLES.items():
        assert suffixes, stage
        for suffix in suffixes:
            assert suffix in known, f"{stage}: {suffix} is in neither list"


def test_the_status_columns_are_in_the_contract():
    metadata = next(r for r in CONTRACT if r.suffix == ".metadata.parquet")
    columns = {c.name: c for c in metadata.columns}
    for stage in STAGES:
        for name in (f"{stage}_status", f"{stage}_error"):
            rule = columns[name]
            assert rule.kind is Kind.STRING and rule.vary is Vary.NEVER
            assert not rule.required, "directories from before this release lack them"


@pytest.mark.parametrize(
    "record_fmt, tool_fmt, tool_compress",
    [("parquet", "tsv", False), ("tsv", "parquet", False), ("tsv", "tsv", True)],
)
def test_a_single_tool_never_records_into_a_format_it_did_not_write(
    tmp_path, record_fmt, tool_fmt, tool_compress
):
    """The record must describe the tables beside it. A single tool at `-F tsv` after a Parquet
    run-all rewrote `.TFBS.tsv`, not the `.TFBS.parquet` downstream reads; recording `ran` in
    `.metadata.parquet` would hide that stale table from the gate. Nor may it create a copy:
    a `.metadata.parquet` in a TSV directory is a completion marker (invariant #2)."""
    _marker(tmp_path, record_fmt)
    before = sorted(p.name for p in tmp_path.iterdir())
    assert not record_stage_outcomes(
        tmp_path,
        "S1",
        {"uxm": Outcome("ran")},
        output_format=tool_fmt,
        compress=tool_compress,
        owns_record=False,
    )
    assert sorted(p.name for p in tmp_path.iterdir()) == before
    assert "uxm_status" not in _row(tmp_path, record_fmt)


@pytest.mark.parametrize(
    "tool_fmt, updated, untouched",
    [("tsv", "tsv", "parquet"), ("parquet", "parquet", "tsv")],
)
def test_a_single_tool_updates_only_the_copies_in_its_own_format(
    tmp_path, tool_fmt, updated, untouched
):
    _marker(tmp_path, "both")
    record_stage_outcomes(
        tmp_path,
        "S1",
        {"uxm": Outcome("ran")},
        output_format=tool_fmt,
        owns_record=False,
        compress=False,
    )
    assert _row(tmp_path, updated)["uxm_status"] == "ran"
    assert "uxm_status" not in _row(tmp_path, untouched)


def test_a_single_tool_writing_both_updates_both(tmp_path):
    _marker(tmp_path, "both")
    record_stage_outcomes(
        tmp_path,
        "S1",
        {"uxm": Outcome("ran")},
        output_format="both",
        owns_record=False,
        compress=False,
    )
    assert _row(tmp_path, "parquet")["uxm_status"] == "ran"
    assert _row(tmp_path, "tsv")["uxm_status"] == "ran"
