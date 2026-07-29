"""The gate must fail on bad input.

On real output the gate currently only ever reports the same four findings, so
without deliberately broken fixtures an assertion that it "passes" would be
indistinguishable from `return True`. Everything here builds a directory that
is wrong in one specific way and checks the gate says so.
"""

from __future__ import annotations

import pandas as pd
import pytest

from krewlyzer.validate import EXIT_PASS, EXIT_STRUCTURAL, EXIT_VIOLATION, run
from krewlyzer.validate.contract import CONTRACT, ColumnRule, Kind, Vary
from krewlyzer.validate.findings import Category, Severity

from .synth_outputs import write_cohort

pytestmark = pytest.mark.invariant


def _ids(result, category=None):
    return {f.id for f in result.findings if category is None or f.category is category}


def test_wellformed_cohort_passes(tmp_path):
    write_cohort(tmp_path, n_samples=3)
    result = run(tmp_path)
    assert result.exit_code == EXIT_PASS, [
        f.message for f in result.findings if f.severity is Severity.ERROR
    ]


def test_missing_completion_marker_is_an_error(tmp_path):
    write_cohort(tmp_path, n_samples=3)
    (tmp_path / "S1" / "S1.metadata.parquet").unlink()

    result = run(tmp_path)

    assert result.exit_code == EXIT_VIOLATION
    marker = [f for f in result.findings if f.category is Category.COMPLETION]
    assert (
        marker
    ), "a sample without .metadata.parquet is dropped from the cohort silently"
    assert marker[0].samples == ["S1"]


def test_constant_column_is_a_degeneracy_error(tmp_path):
    """The WPS_background failure mode, reproduced deliberately."""
    write_cohort(tmp_path, n_samples=3)
    for sample in ("S0", "S1", "S2"):
        path = tmp_path / sample / f"{sample}.WPS_background.parquet"
        df = pd.read_parquet(path)
        df["nrl_bp"] = 150.0
        df.to_parquet(path, index=False)

    result = run(tmp_path)

    assert result.exit_code == EXIT_VIOLATION
    assert "WPS_background.DEGENERACY.nrl_bp" in _ids(result, Category.DEGENERACY)


def test_single_sample_skips_rather_than_passes_degeneracy(tmp_path):
    """One sample cannot demonstrate variation; SKIP is the honest verdict."""
    write_cohort(tmp_path, n_samples=1)

    result = run(tmp_path, min_samples=3)

    skipped = [
        f
        for f in result.findings
        if f.category is Category.DEGENERACY and f.severity is Severity.SKIP
    ]
    assert skipped, "cross-sample degeneracy must not silently pass at n=1"
    assert not [
        f
        for f in result.findings
        if f.category is Category.DEGENERACY and f.severity is Severity.ERROR
    ], "a single sample cannot prove a column is constant either"


def test_min_samples_cannot_be_lowered_below_two(tmp_path):
    """--min-samples 1 must not turn every column into a false positive.

    With a single sample every signature is trivially identical, so an
    unclamped floor would report the entire cohort as degenerate.
    """
    write_cohort(tmp_path, n_samples=1)

    result = run(tmp_path, min_samples=1)

    assert not [
        f
        for f in result.findings
        if f.category is Category.DEGENERACY and f.severity is Severity.ERROR
    ], "one sample cannot establish that anything is constant"


def test_unreadable_parquet_is_structural_not_a_violation(tmp_path):
    """Exit 2 means 'not comparable' so a caller can retry instead of paging."""
    write_cohort(tmp_path, n_samples=3)
    (tmp_path / "S1" / "S1.FSC.parquet").write_bytes(b"not a parquet file")

    result = run(tmp_path)

    assert result.exit_code == EXIT_STRUCTURAL
    assert any(f.category is Category.STRUCTURAL for f in result.findings)


def test_empty_directory_is_structural(tmp_path):
    result = run(tmp_path)
    assert result.exit_code == EXIT_STRUCTURAL
    assert "INPUT.NO_SAMPLES" in _ids(result)


def test_stray_numeric_column_in_fsd_is_caught(tmp_path):
    """Consumers treat every unexpected numeric FSD column as a size bucket."""
    write_cohort(tmp_path, n_samples=3)
    path = tmp_path / "S0" / "S0.FSD.parquet"
    df = pd.read_parquet(path)
    df["mean_gc"] = 0.41
    df.to_parquet(path, index=False)

    result = run(tmp_path)

    assert result.exit_code == EXIT_VIOLATION
    assert "FSD.FSD_ONLY_SIZE_BINS" in _ids(result, Category.DOMAIN)


def test_bare_chromosome_names_are_caught(tmp_path):
    write_cohort(tmp_path, n_samples=3)
    path = tmp_path / "S0" / "S0.FSC.parquet"
    df = pd.read_parquet(path)
    df["chrom"] = [c.removeprefix("chr") for c in df["chrom"]]
    df.to_parquet(path, index=False)

    result = run(tmp_path)

    assert result.exit_code == EXIT_VIOLATION
    assert "FSC.CHR_PREFIXED" in _ids(result, Category.DOMAIN)


def test_missing_consumed_table_is_an_error(tmp_path):
    write_cohort(tmp_path, n_samples=3)
    (tmp_path / "S2" / "S2.MDS.gene.parquet").unlink()

    result = run(tmp_path)

    assert result.exit_code == EXIT_VIOLATION
    assert any(
        f.category is Category.MISSING
        and f.severity is Severity.ERROR
        and f.table == ".MDS.gene.parquet"
        for f in result.findings
    )


def test_declaring_a_column_constant_requires_a_reason():
    """Silencing a degeneracy finding must cost a written justification."""
    with pytest.raises(ValueError, match="constant_reason"):
        ColumnRule("nrl_bp", Kind.NUMERIC, Vary.NEVER)

    ok = ColumnRule("group_id", Kind.STRING, Vary.NEVER, constant_reason="a key")
    assert ok.vary is Vary.NEVER


def test_contract_covers_the_consumer_surface():
    """Guards against a table quietly leaving the contract."""
    suffixes = {r.suffix for r in CONTRACT}
    for required in (
        ".metadata.parquet",
        ".WPS_background.parquet",
        ".MDS.gene.parquet",
        ".FSC.gene.parquet",
        ".EndMotif.parquet",
        ".OCF.ontarget.parquet",
    ):
        assert required in suffixes, f"{required} dropped out of the contract"
    assert len(suffixes) == len(CONTRACT), "duplicate suffix in CONTRACT"
