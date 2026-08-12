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


def test_a_skip_pon_cohort_passes_without_the_pon_columns(tmp_path):
    """`--skip-pon` is a legitimate run, not a broken one.

    A run with no PON, one whose PON the version guard refused, and a
    deliberate `--skip-pon` all produce output with no z-scores. Requiring the
    PON-derived columns unconditionally would fail all three.
    """
    write_cohort(tmp_path, n_samples=3, pon_applied=False)
    result = run(tmp_path)
    assert result.exit_code == EXIT_PASS, [
        f.message for f in result.findings if f.severity is Severity.ERROR
    ]


def test_pon_scoring_that_vanished_is_caught(tmp_path):
    """The gap this closes.

    Twice in 0.9.0 the WPS scoring was written somewhere nothing reads -- once
    to a `.WPS.tsv`, once to `{sample}.parquet` after `with_suffix` ate the
    `.WPS`. Both times a raw `.WPS.parquet` shipped and passed every check,
    because the contract declared none of the columns PON scoring produces.

    Now a sample that says it was scored must show the evidence.
    """
    write_cohort(tmp_path, n_samples=3)
    for sample in ("S0", "S1", "S2"):
        path = tmp_path / sample / f"{sample}.WPS.parquet"
        frame = pd.read_parquet(path)
        frame.drop(
            columns=[c for c in frame.columns if c.startswith("wps_nuc_z")]
        ).to_parquet(path, index=False)

    result = run(tmp_path)
    assert result.exit_code != EXIT_PASS, "scoring vanished and the gate passed"
    ids = _ids(result)
    assert "WPS.MISSING_COLUMN.wps_nuc_z" in ids, ids
    message = next(f.message for f in result.findings if "wps_nuc_z" in f.message)
    assert "pon_applied=True" in message, (
        "the finding must say why the column was expected, or a reader cannot "
        "tell it from a --skip-pon run"
    )


def test_an_unscored_sample_claiming_otherwise_is_caught(tmp_path):
    """The marker and the tables must agree.

    `pon_applied=True` with no z-scores anywhere is the shape of a run whose
    scoring step failed after the metadata was stamped.
    """
    write_cohort(tmp_path, n_samples=2, pon_applied=False)
    for sample in ("S0", "S1"):
        marker = tmp_path / sample / f"{sample}.metadata.parquet"
        frame = pd.read_parquet(marker)
        frame["pon_applied"] = True
        frame.to_parquet(marker, index=False)

    result = run(tmp_path)
    assert result.exit_code != EXIT_PASS
    ids = _ids(result)
    assert any(i.startswith("WPS.MISSING_COLUMN.") for i in ids), ids
    assert any(i.startswith("FSD.MISSING_COLUMN.") for i in ids), ids


def test_a_directory_with_no_provenance_is_reported(tmp_path):
    """A cohort you cannot identify is a cohort you cannot bless.

    0.9.0 changed what several columns mean, so "which build wrote this?" is
    not optional metadata. A directory from an older build says nothing, and
    the gate says so rather than assuming.
    """
    write_cohort(tmp_path, n_samples=2)
    for sample in ("S0", "S1"):
        marker = tmp_path / sample / f"{sample}.metadata.parquet"
        frame = pd.read_parquet(marker)
        frame.drop(columns=["krewlyzer_version"]).to_parquet(marker, index=False)

    result = run(tmp_path)
    assert result.exit_code != EXIT_PASS
    assert "metadata.MISSING_COLUMN.krewlyzer_version" in _ids(result)


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


def test_scatter_gather_matches_a_single_pass(tmp_path):
    """Per-sample fingerprints plus a gather must equal the monolithic run.

    This is the property the Nextflow split depends on: the scatter step never
    sees another sample, so degeneracy has to survive the round trip through
    fingerprints intact.
    """
    from krewlyzer.validate.gate import Fingerprint, check_sample, evaluate_cohort

    write_cohort(tmp_path, n_samples=3)
    for sample in ("S0", "S1", "S2"):
        path = tmp_path / sample / f"{sample}.WPS_background.parquet"
        df = pd.read_parquet(path)
        df["nrl_bp"] = 150.0
        df.to_parquet(path, index=False)

    fps = []
    for sample in ("S0", "S1", "S2"):
        sample_findings, fp = check_sample(sample, tmp_path / sample)
        assert not [f for f in sample_findings if f.severity is Severity.ERROR], (
            "a constant column is invisible to a single sample -- that is why "
            "the cohort step exists"
        )
        out = tmp_path / f"{sample}.fingerprint.json"
        fp.save(out)
        fps.append(Fingerprint.load(out))

    cohort = evaluate_cohort(fps)

    assert any(
        f.id == "WPS_background.DEGENERACY.nrl_bp" and f.severity is Severity.ERROR
        for f in cohort
    ), "the gather step must catch what no single sample can"


def test_run_all_validation_artifacts(tmp_path):
    """Pins the composition run-all uses, so a signature change can't break it.

    run-all writes both artifacts on every Parquet run and never changes its
    exit code without --strict-validation. The fingerprint in particular has to
    land: it is what makes validate-cohort affordable later.
    """
    from krewlyzer.validate.gate import Result, check_sample
    from krewlyzer.validate.report import to_json

    from .synth_outputs import write_sample

    write_sample(tmp_path, "S0", 0)
    sample_dir = tmp_path / "S0"

    findings, fingerprint = check_sample("S0", sample_dir)
    fingerprint.save(sample_dir / "S0.fingerprint.json")
    to_json(
        Result(findings=findings, samples=["S0"], fingerprints=[fingerprint]),
        sample_dir / "S0.validation.json",
        "0.8.3",
    )

    assert (sample_dir / "S0.fingerprint.json").exists()
    assert (sample_dir / "S0.validation.json").exists()
    assert not [f for f in findings if f.severity is Severity.ERROR]


def test_fingerprint_version_mismatch_is_rejected(tmp_path):
    """A stale fingerprint must not be silently compared against a fresh one."""
    from krewlyzer.validate.gate import Fingerprint

    path = tmp_path / "old.json"
    path.write_text('{"fingerprint_version": "0", "sample": "S0", "observations": {}}')

    with pytest.raises(ValueError, match="fingerprint version"):
        Fingerprint.load(path)


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
