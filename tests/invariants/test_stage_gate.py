"""validate-output reads what each optional stage did.

Before, a failed Region-MDS surfaced as "MDS.gene is absent" with no reason, a failed mFSD not at
all (its table is not consumed), and a failed MDS z-score passed outright, since `mds_z` is a
column the contract never required. Each directory here is well formed except for one recorded
outcome, and the gate must say what that outcome means -- and say nothing new about a directory
that records none.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from krewlyzer.validate.contract import CONTRACT, STAGE_TABLES, STAGES
from krewlyzer.validate.findings import Severity
from krewlyzer.validate.gate import check_sample

from .synth_outputs import write_sample

pytestmark = pytest.mark.invariant


def _sample(tmp_path: Path, **statuses: str) -> Path:
    """A well-formed sample whose record says ``stage=status`` (``status:Class`` for failed)."""
    directory = write_sample(tmp_path, "S0", 0)
    marker = directory / "S0.metadata.parquet"
    frame = pd.read_parquet(marker)
    for stage, value in statuses.items():
        status, _, error = value.partition(":")
        frame[f"{stage}_status"] = status
        frame[f"{stage}_error"] = error
    frame.to_parquet(marker, index=False)
    return directory


def _drop(directory: Path, *suffixes: str) -> None:
    for suffix in suffixes:
        (directory / f"S0{suffix}").unlink(missing_ok=True)


def _add(directory: Path, suffix: str) -> None:
    pd.DataFrame({"x": [1.0, 2.0]}).to_parquet(directory / f"S0{suffix}")


def _drop_column(directory: Path, suffix: str, column: str) -> None:
    path = directory / f"S0{suffix}"
    pd.read_parquet(path).drop(columns=[column]).to_parquet(path, index=False)


def _findings(directory: Path):
    return check_sample("S0", directory)[0]


def _ids(findings, severity=None):
    return {f.id for f in findings if severity is None or f.severity is severity}


def test_a_failed_consumed_stage_is_one_error_that_says_why(tmp_path):
    """It replaces the two "absent" findings, which gave no reason, and keeps their severity."""
    failed = _sample(tmp_path / "a", region_mds="failed:RuntimeError")
    ran = _sample(tmp_path / "b", region_mds="ran")
    for directory in (failed, ran):
        _drop(directory, *STAGE_TABLES["region_mds"])

    reported = _findings(failed)
    stage = [f for f in reported if f.id == "STAGE.FAILED"]
    assert len(stage) == 1 and stage[0].severity is Severity.ERROR
    assert "region_mds" in stage[0].message and "RuntimeError" in stage[0].message
    assert not {"MDS.gene.ABSENT", "MDS.exon.ABSENT"} & _ids(reported)

    # Two inputs, two outputs: the same absent tables under a stage that says it ran.
    assert {"MDS.gene.ABSENT", "MDS.exon.ABSENT"} <= _ids(
        _findings(ran), Severity.ERROR
    )
    assert "STAGE.FAILED" not in _ids(_findings(ran))


def test_a_failed_stage_nothing_downstream_reads_is_a_warning(tmp_path):
    directory = _sample(tmp_path, mfsd="failed:FileNotFoundError")
    reported = _findings(directory)
    stage = [f for f in reported if f.id == "STAGE.FAILED"]
    assert [f.severity for f in stage] == [Severity.WARN]
    assert "mFSD.ABSENT" not in _ids(reported)
    assert not _ids(reported, Severity.ERROR)


def test_a_table_beside_a_failed_stage_is_an_error(tmp_path):
    """A deliberate new ERROR: the table is a previous run's, or a partial one, and downstream
    would read it as this run's."""
    directory = _sample(
        tmp_path, region_mds="failed:RuntimeError"
    )  # tables left in place
    stale = [f for f in _findings(directory) if f.id == "STAGE.STALE_TABLE"]
    assert {f.table for f in stale} == set(STAGE_TABLES["region_mds"])
    assert all(f.severity is Severity.ERROR for f in stale)


def test_a_stale_table_nothing_reads_is_a_warning(tmp_path):
    directory = _sample(tmp_path, mfsd="failed:OSError")
    _add(directory, ".mFSD.parquet")
    stale = [f for f in _findings(directory) if f.id == "STAGE.STALE_TABLE"]
    assert [f.severity for f in stale] == [Severity.WARN]


def test_a_stage_nobody_asked_for_and_nothing_reads_gives_no_finding(tmp_path):
    """No `--variants` means no mFSD, correctly; it was a WARN on every such run."""
    directory = _sample(tmp_path, uxm="not_requested", mfsd="no_variants")
    tables = {".UXM.parquet", ".mFSD.parquet"}
    assert not [f for f in _findings(directory) if f.table in tables]
    old = write_sample(tmp_path / "old", "S0", 0)
    assert {"UXM.ABSENT", "mFSD.ABSENT"} <= _ids(_findings(old), Severity.WARN)


def test_a_consumed_stage_not_requested_is_still_absent(tmp_path):
    """Downstream reads `MDS.gene` whether or not `--assay` was given; making the contract
    conditional on the assay is its own change."""
    directory = _sample(tmp_path, region_mds="not_requested")
    _drop(directory, *STAGE_TABLES["region_mds"])
    errors = _ids(_findings(directory), Severity.ERROR)
    assert {"MDS.gene.ABSENT", "MDS.exon.ABSENT"} <= errors


def test_a_table_beside_a_stage_not_requested_is_a_warning(tmp_path):
    """Left from an earlier configuration of the same directory."""
    directory = _sample(tmp_path, uxm="not_requested")
    _add(directory, ".UXM.parquet")
    found = [f for f in _findings(directory) if f.table == ".UXM.parquet"]
    assert [(f.id, f.severity) for f in found] == [
        ("STAGE.UNEXPECTED_TABLE", Severity.WARN)
    ]


def test_a_stage_whose_input_was_missing_says_so(tmp_path):
    directory = _sample(tmp_path, uxm="unavailable")
    found = [f for f in _findings(directory) if f.table == ".UXM.parquet"]
    assert [(f.id, f.severity) for f in found] == [("STAGE.UNAVAILABLE", Severity.WARN)]


def test_a_directory_with_no_statuses_behaves_as_before(tmp_path):
    """Every directory written before this release. The one change is `MDS.gene`'s z columns,
    which are no longer required without a record that its z-scoring ran."""
    directory = write_sample(tmp_path, "S0", 0)
    reported = _findings(directory)
    assert not _ids(reported, Severity.ERROR)
    assert not {i for i in _ids(reported) if i.startswith("STAGE.")}
    assert {"UXM.ABSENT", "mFSD.ABSENT"} <= _ids(reported, Severity.WARN)

    _drop_column(directory, ".MDS.gene.parquet", "mds_z")
    _drop_column(directory, ".MDS.parquet", "mds_z")
    assert not _ids(_findings(directory), Severity.ERROR)


@pytest.mark.parametrize(
    "suffix, stage",
    [
        (".MDS.parquet", "mds_z"),
        (".MDS.ontarget.parquet", "mds_z_ontarget"),
        (".MDS.gene.parquet", "region_mds_z"),
    ],
)
def test_a_z_score_is_required_only_where_its_stage_ran(tmp_path, suffix, stage):
    """`mds_z` vanishing passed every check. Now it is required where the record says it was
    written -- and only there: a BED-built PON has no MDS baseline, and requiring it on
    `pon_applied` alone would be a new false ERROR."""
    results = {}
    for status in ("ran", "unavailable", "failed:TypeError"):
        directory = _sample(tmp_path / status.split(":")[0], **{stage: status})
        _drop_column(directory, suffix, "mds_z")
        results[status] = _findings(directory)

    family = suffix.lstrip(".").replace(".parquet", "")
    assert f"{family}.MISSING_COLUMN.mds_z" in _ids(results["ran"], Severity.ERROR)
    assert not _ids(results["unavailable"], Severity.ERROR)
    failed = results["failed:TypeError"]
    assert f"{family}.MISSING_COLUMN.mds_z" not in _ids(failed)
    assert [f.severity for f in failed if f.id == "STAGE.FAILED"] == [Severity.ERROR]


def test_every_column_stage_is_named_by_a_column_rule():
    """`mds_z` and friends add columns, not tables; the columns say which stage writes them."""
    named = {c.requires_stage for r in CONTRACT for c in r.columns if c.requires_stage}
    assert named <= set(STAGES)
    assert named | set(STAGE_TABLES) == set(STAGES)


@pytest.mark.parametrize("status", ["not_requested", "unavailable"])
def test_a_consumed_table_beside_a_stage_that_did_not_run_is_an_error(tmp_path, status):
    """`MDS.gene` from an earlier `--assay` run, beside a run without one: downstream reads it
    as this run's, the same harm as a table beside a failed stage."""
    directory = _sample(tmp_path, region_mds=status)  # tables left in place
    found = [f for f in _findings(directory) if f.id == "STAGE.UNEXPECTED_TABLE"]
    assert {f.table for f in found} == set(STAGE_TABLES["region_mds"])
    assert all(f.severity is Severity.ERROR for f in found)


@pytest.mark.parametrize(
    "suffix, stage",
    [
        (".MDS.parquet", "mds_z"),
        (".MDS.ontarget.parquet", "mds_z_ontarget"),
        (".MDS.gene.parquet", "region_mds_z"),
    ],
)
def test_a_table_rewritten_unscored_in_a_scored_sample_is_an_error(
    tmp_path, suffix, stage
):
    """`motif --skip-pon` into a folder run-all scored: `MDS` loses `mds_z`, the record says
    `mds_z=not_requested` (true), and `pon_applied` stays True (a single tool may only raise it).
    Downstream reads a scored sample with an unscored table -- flagged, not hidden. In an unscored
    sample, `not_requested` is simply the truth and needs no column."""
    family = suffix.lstrip(".").replace(".parquet", "")
    found = {}
    for applied in (True, False):
        directory = _sample(tmp_path / str(applied), **{stage: "not_requested"})
        marker = directory / "S0.metadata.parquet"
        frame = pd.read_parquet(marker)
        frame["pon_applied"] = applied
        frame.to_parquet(marker, index=False)
        _drop_column(directory, suffix, "mds_z")
        found[applied] = {f.id: f for f in _findings(directory)}

    missing = f"{family}.MISSING_COLUMN.mds_z"
    assert missing in found[True] and found[True][missing].severity is Severity.ERROR
    assert "not_requested" in found[True][missing].message
    assert missing not in found[False]
