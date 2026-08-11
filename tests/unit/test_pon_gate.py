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
    """A PON complete enough to reach the check under test.

    The core blocks come from `pon_fixtures`: since the gate started checking a
    packing list, a file holding only metadata plus one block fails for being
    incomplete, which would mask whatever the test was actually about.
    """
    from . import pon_fixtures

    meta = [f for f in frames if (f["table"] == "metadata").all()]
    rest = [f for f in frames if not (f["table"] == "metadata").all()]
    supplied = {str(f["table"].iloc[0]) for f in rest}
    filler = pon_fixtures.core_blocks(exclude=tuple(supplied))

    path = tmp_path / "t.pon.parquet"
    pd.concat(meta + filler + rest, ignore_index=True).to_parquet(path)
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
        # Fragment BEDs, so region_mds and mds_baseline are legitimately
        # absent and warn rather than fail.
        "input_kind": "bed",
        "panel_mode": False,
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
@pytest.mark.xfail(
    strict=True,
    reason=(
        "The models in the tree predate `breakpoint_motif_baseline`, added "
        "later in 0.9.0 when BreakPointMotif was found to be scored against "
        "the EndMotif baseline. They are not wrong, only incomplete, and the "
        "re-aggregation that adds the block is the next release step. "
        "strict=True on purpose: when the rebuilt models land this test starts "
        "passing, the strict marker turns that into a failure, and whoever "
        "sees it removes the marker. It cannot be forgotten."
    ),
)
def test_the_shipped_models_pass_their_own_gate(path):
    """This is the acceptance record for the 0.9.0 rebuild.

    It used to assert the opposite. The four models shipped before 0.9.0 each
    carried a hardcoded `167.0 / 5.0` in `wps_background` -- byte-identical
    across two different cohorts -- and no provenance at all, and this test
    asserted the gate caught both. Its docstring said it would flip when the
    rebuild landed, and this is that flip.

    The four now in the tree were built from run-all output via
    `--from-outputs`, and pass with zero findings: no fabricated baseline, no
    non-positive sigma, every required block present, and a cohort digest that
    identifies what each was built from.
    """
    findings = check_pon(path)
    assert exit_code(findings) == 0, [
        f"{f.id} {f.table}.{f.column}: {f.message}" for f in findings
    ]


@pytest.mark.skipif(not _SHIPPED, reason="bundled PONs not available")
@pytest.mark.parametrize("path", _SHIPPED, ids=lambda p: p.parent.name + "/" + p.stem)
def test_the_shipped_models_record_what_they_were_built_from(path):
    """None of the pre-0.9.0 models could be reproduced or told apart."""
    import pandas as pd

    meta = pd.read_parquet(path)
    meta = meta[meta["table"] == "metadata"].iloc[0]
    assert str(meta.get("cohort_digest", "")).strip(), "no cohort digest"
    assert float(meta["n_samples"]) >= MIN_SAMPLES
    assert str(meta.get("input_kind", "")).strip(), "no input_kind"


@pytest.mark.skipif(len(_SHIPPED) < 2, reason="need two bundled PONs to compare")
def test_no_two_shipped_models_share_a_baseline():
    """Invariant #1, applied to the shipped artifacts themselves.

    The defect that started this release was four models carrying a
    byte-identical `wps_background`. A single assertion that two of them differ
    would have caught it in March.
    """
    import pandas as pd

    seen = {}
    for path in _SHIPPED:
        table = pd.read_parquet(path)
        block = table[table["table"] == "wps_background"].set_index("group_id")
        seen[path.stem] = pd.to_numeric(block["nrl_mean"], errors="coerce")

    names = list(seen)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            common = seen[a].index.intersection(seen[b].index)
            assert (
                not seen[a].loc[common].equals(seen[b].loc[common])
            ), f"{a} and {b} carry an identical wps_background"


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


# ---------------------------------------------------------------------------
# the packing list
# ---------------------------------------------------------------------------

#: A sigma that varies across groups, i.e. one that was actually fitted.
_FITTED = [4.0, 5.1, 6.2, 4.8, 5.5]


def test_a_block_that_vanished_is_noticed(tmp_path):
    """Every other check reads the blocks present and skips the rest.

    Demonstrated on a real model before this existed: deleting `region_mds`, or
    `fsc_gene_baseline`, or almost every baseline at once, produced an
    identical finding list each time. A block absent from the file is invisible
    to a check that iterates the file.
    """
    from . import pon_fixtures

    path = pon_fixtures.complete(
        tmp_path, _background(_FITTED), exclude=("ocf_baseline",)
    )
    findings = check_pon(path)
    missing = [f for f in findings if f.id == "PON.BLOCK_MISSING"]
    assert missing, "a missing core block went unnoticed"
    assert "ocf_baseline" in missing[0].message
    assert exit_code(findings) == 1


def test_a_complete_model_is_not_accused_of_missing_anything(tmp_path):
    """The other direction, so the check cannot be 'always complain'."""
    from . import pon_fixtures

    path = pon_fixtures.complete(tmp_path, _background(_FITTED))
    ids = {f.id for f in check_pon(path)}
    assert "PON.BLOCK_MISSING" not in ids
    assert exit_code(check_pon(path)) == 0


def test_bam_only_blocks_warn_for_a_bed_cohort_and_fail_for_a_bam_one(tmp_path):
    """`input_kind` is what makes this decidable.

    `mds_baseline` and `region_mds` need a BAM, so a fragment-BED cohort lacks
    them legitimately. Without the field the gate could only warn either way --
    which is how `region_mds` vanishing from a BAM build stayed a warning.
    """
    from . import pon_fixtures

    for name in ("a", "b"):
        (tmp_path / name).mkdir(parents=True, exist_ok=True)
    bed = pon_fixtures.complete(tmp_path / "a", _background(_FITTED), input_kind="bed")
    bam = pon_fixtures.complete(tmp_path / "b", _background(_FITTED), input_kind="bam")

    bed_finding = [f for f in check_pon(bed) if f.id == "PON.BLOCK_MISSING_BAM_ONLY"]
    bam_finding = [f for f in check_pon(bam) if f.id == "PON.BLOCK_MISSING_BAM_ONLY"]
    assert bed_finding and bed_finding[0].severity is Severity.WARN
    assert bam_finding and bam_finding[0].severity is Severity.ERROR
    assert exit_code(check_pon(bed)) == 0, "a BED cohort must not be failed for this"
    assert exit_code(check_pon(bam)) == 1


def test_panel_mode_requires_the_on_target_blocks(tmp_path):
    """A panel PON without them scores on-target regions against nothing."""
    from . import pon_fixtures

    path = pon_fixtures.complete(
        tmp_path, _background(_FITTED), panel_mode=True, input_kind="bed"
    )
    missing = [f for f in check_pon(path) if f.id == "PON.BLOCK_MISSING"]
    assert missing
    assert "gc_bias_ontarget" in missing[0].message
    assert "wps_baseline_panel" in missing[0].message


# ---------------------------------------------------------------------------
# the blocks whose sigma is shaped differently
# ---------------------------------------------------------------------------


def test_a_degenerate_gc_bias_is_caught(tmp_path):
    """It was unchecked because its sigma sits in per-bin columns.

    `gc_bias` is the worst one to leave unchecked: it reaches every FSC
    log-ratio and every FSR normalisation through `PonModel.get_mean`.
    """
    from . import pon_fixtures

    flat = pd.DataFrame(
        {
            "table": "gc_bias",
            "gc_bin": [0.30, 0.35, 0.40, 0.45, 0.50],
            "short_expected": [1.0] * 5,
            "short_std": [0.05] * 5,  # the signature: one value, every bin
            "intermediate_expected": [1.0] * 5,
            "intermediate_std": [0.06, 0.07, 0.08, 0.09, 0.10],
            "long_expected": [1.0] * 5,
            "long_std": [0.07, 0.08, 0.09, 0.10, 0.11],
        }
    )
    path = pon_fixtures.complete(
        tmp_path, _background(_FITTED), flat, exclude=("gc_bias",)
    )
    hits = [
        f
        for f in check_pon(path)
        if f.id == "PON.BLOCK_DEGENERATE" and f.table == "gc_bias"
    ]
    assert hits, "a constant gc_bias sigma was not caught"
    assert hits[0].column == "short_std"


def test_a_nonpositive_fsd_sigma_is_caught(tmp_path):
    """FSD supplies the log-ratios; a zero sigma there is an infinite z."""
    from . import pon_fixtures

    bad = pd.DataFrame(
        {
            "table": "fsd_baseline",
            "arm": ["1p"] * 5,
            "size_bin": [65, 70, 75, 80, 85],
            "expected": [0.02, 0.03, 0.04, 0.05, 0.06],
            "std": [0.004, 0.0, 0.006, 0.007, 0.008],
        }
    )
    path = pon_fixtures.complete(
        tmp_path, _background(_FITTED), bad, exclude=("fsd_baseline",)
    )
    hits = [f for f in check_pon(path) if f.id == "PON.NONPOSITIVE_SIGMA"]
    assert any(f.table == "fsd_baseline" for f in hits)
