"""`validate-output` gated the results. Nothing gated the reference.

Every PON defect this release fixes sat in four shipped files, visible in the
file and invisible to every check:

- `wps_background` held a hardcoded `167.0 / 5.0 / 0.0 / 1.0`, identical across
  all 28 groups and all four models, from cohorts of 21 and 47 samples
- six σ floors turned "no spread measured" into a divisor
- four blocks were built, shipped and read by nothing
- no model recorded what it was built from, so none could be reproduced

The load-bearing assertion is the same invariant the output gate enforces one
level up: **a baseline that cannot vary with its cohort was not fitted to one.**
A single check that σ differs between groups would have caught `wps_background`
in March.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from krewlyzer.validate.findings import Category, Severity
from krewlyzer.validate.pon_gate import MIN_SAMPLES, check_pon, describe, exit_code

pytestmark = pytest.mark.unit

_SHIPPED = sorted(
    (Path(__file__).resolve().parents[2] / "src/krewlyzer/data/pon/GRCh37").glob(
        "*/*.pon.parquet"
    )
)


def _write(tmp_path: Path, *frames: pd.DataFrame) -> Path:
    path = tmp_path / "t.pon.parquet"
    pd.concat(frames, ignore_index=True).to_parquet(path)
    return path


def _metadata(**overrides) -> pd.DataFrame:
    row = {
        "table": "metadata",
        "schema_version": "1.0",
        "assay": "xs2",
        "n_samples": 21.0,
        "krewlyzer_version": "0.9.0",
        "cohort_digest": "deadbeefcafe0000",
        "cohort_label": "healthy-donors",
    }
    row.update(overrides)
    return pd.DataFrame([row])


def _background(nrl_std) -> pd.DataFrame:
    groups = [f"Chr{i}_All" for i in range(1, 6)]
    return pd.DataFrame(
        {
            "table": "wps_background",
            "group_id": groups,
            "nrl_mean": [190.0] * len(groups),
            "nrl_std": (
                [nrl_std] * len(groups) if not isinstance(nrl_std, list) else nrl_std
            ),
        }
    )


# ---------------------------------------------------------------------------
# the acceptance criterion from the plan
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _SHIPPED, reason="bundled PONs not available (git lfs pull)")
@pytest.mark.parametrize("path", _SHIPPED, ids=lambda p: p.parent.name + "/" + p.stem)
def test_the_gate_rejects_the_currently_shipped_models(path):
    """A gate that passes the models we know are wrong is not a gate.

    All four carry the fabricated `wps_background` and no provenance. When the
    0.9.0 rebuild lands this test flips to asserting they pass — and that flip
    is the acceptance record for the rebuild.
    """
    findings = check_pon(path)
    assert exit_code(findings) == 1, "expected a contract violation"
    ids = {f.id for f in findings}
    assert "PON.BLOCK_DEGENERATE" in ids, "the fabricated wps_background was not caught"
    assert "PON.NO_VERSION" in ids, "the missing provenance was not caught"


@pytest.mark.skipif(not _SHIPPED, reason="bundled PONs not available")
def test_the_fabricated_sigma_is_named_by_value(path=_SHIPPED[0] if _SHIPPED else None):
    """The message has to be actionable without opening the file."""
    degenerate = [f for f in check_pon(path) if f.id == "PON.BLOCK_DEGENERATE"]
    assert degenerate
    assert any(
        "5.0" in f.message for f in degenerate
    ), "the message should quote the repeated value"


# ---------------------------------------------------------------------------
# the individual checks
# ---------------------------------------------------------------------------


def test_a_fitted_baseline_passes(tmp_path):
    """The other direction: real variation must not trip the gate."""
    path = _write(tmp_path, _metadata(), _background([4.0, 5.1, 6.2, 4.8, 5.5]))
    assert exit_code(check_pon(path)) == 0, [f.message for f in check_pon(path)]


def test_one_sigma_repeated_everywhere_is_an_error(tmp_path):
    """The `wps_background` signature, in isolation."""
    path = _write(tmp_path, _metadata(), _background(5.0))
    findings = [f for f in check_pon(path) if f.id == "PON.BLOCK_DEGENERATE"]
    assert findings and findings[0].severity is Severity.ERROR
    assert findings[0].category is Category.DEGENERACY


def test_a_single_entry_block_is_not_called_degenerate(tmp_path):
    """One value cannot be constant *across* anything.

    Guards against the check firing on blocks that legitimately hold one row,
    which would train people to ignore it.
    """
    single = pd.DataFrame(
        {
            "table": ["wps_background"],
            "group_id": ["Global_All"],
            "nrl_mean": [190.0],
            "nrl_std": [5.0],
        }
    )
    path = _write(tmp_path, _metadata(), single)
    assert not [f for f in check_pon(path) if f.id == "PON.BLOCK_DEGENERATE"]


def test_a_nonpositive_sigma_is_an_error(tmp_path):
    """A z divided by zero is infinite, not conservative."""
    path = _write(tmp_path, _metadata(), _background([4.0, 0.0, 6.0, 5.0, 5.5]))
    assert any(f.id == "PON.NONPOSITIVE_SIGMA" for f in check_pon(path))


def test_a_missing_version_is_an_error(tmp_path):
    """0.9.0 changes what every feature means, so the version is load-bearing."""
    path = _write(
        tmp_path,
        _metadata(krewlyzer_version=""),
        _background([4.0, 5.0, 6.0, 4.5, 5.5]),
    )
    findings = check_pon(path)
    assert any(f.id == "PON.NO_VERSION" for f in findings)
    assert exit_code(findings) == 1


def test_a_missing_cohort_digest_warns_but_does_not_gate(tmp_path):
    """Not reproducible is bad; not comparable is worse. Only one blocks."""
    path = _write(
        tmp_path, _metadata(cohort_digest=""), _background([4.0, 5.0, 6.0, 4.5, 5.5])
    )
    findings = check_pon(path)
    assert any(
        f.id == "PON.NO_COHORT" and f.severity is Severity.WARN for f in findings
    )
    assert exit_code(findings) == 0


def test_too_few_samples_is_an_error(tmp_path):
    path = _write(
        tmp_path, _metadata(n_samples=2.0), _background([4.0, 5.0, 6.0, 4.5, 5.5])
    )
    assert any(f.id == "PON.TOO_FEW_SAMPLES" for f in check_pon(path))
    assert MIN_SAMPLES == 3


def test_entries_backed_by_too_few_samples_are_an_error(tmp_path):
    block = pd.DataFrame(
        {
            "table": "region_mds_exon",
            "gene": ["A", "B", "C"],
            "mds_mean": [0.9, 0.9, 0.9],
            "mds_std": [0.01, 0.02, 0.03],
            "n_samples": [19, 19, 1],
        }
    )
    path = _write(tmp_path, _metadata(), block)
    thin = [f for f in check_pon(path) if f.id == "PON.THIN_ENTRIES"]
    assert thin and thin[0].evidence["thin"] == 1


# ---------------------------------------------------------------------------
# structural
# ---------------------------------------------------------------------------


def test_an_unreadable_file_exits_2_not_1(tmp_path):
    """Structural failures are retryable; contract violations are not."""
    path = tmp_path / "broken.pon.parquet"
    path.write_bytes(b"not a parquet")
    findings = check_pon(path)
    assert exit_code(findings) == 2


def test_something_that_is_not_a_pon_exits_2(tmp_path):
    path = tmp_path / "other.parquet"
    pd.DataFrame({"a": [1]}).to_parquet(path)
    assert exit_code(check_pon(path)) == 2


def test_describe_reports_provenance_without_identifiers(tmp_path):
    path = _write(tmp_path, _metadata(), _background([4.0, 5.0, 6.0, 4.5, 5.5]))
    line = describe(path)
    assert "0.9.0" in line and "deadbeefcafe0000" in line
    assert "healthy-donors" in line
