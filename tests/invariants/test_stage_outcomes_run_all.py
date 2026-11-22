"""Every optional stage records what it did, in the real run-all.

run-all carries on past a failed optional stage on purpose. Until now the only trace of a
failure was a log line, and some paths left none: a `--variants` file that does not exist was
skipped silently, a missing `--mfsd-bam` silently ran mFSD on the main BAM with duplex off, and
panel WPS logged its failures at debug.

Each stage is run four ways on the same synthetic sample -- ran, failed on a real broken input,
unavailable, not requested -- and the record must tell them apart (invariant #1: two inputs, two
outputs). The broken input is what a clone that skipped `git lfs pull` has: a pointer file where
a gzipped asset should be.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import krewlyzer.core.gene_bed as gene_bed_module
from conftest import requires_data
from krewlyzer.assets import AssetManager

from .runall import lfs_pointer, outcomes, record, synth_run_all, write_vcf
from .synth import BASELINE, make_sample

pytestmark = [pytest.mark.invariant, requires_data]

ASSET_STAGES = ("region_mds", "gene_fsc", "panel_wps", "tfbs", "atac")


@pytest.fixture(scope="module")
def sample(tmp_path_factory):
    return make_sample(tmp_path_factory.mktemp("stage_in"), BASELINE)


def _run(sample, out, patches=(), **overrides):
    with pytest.MonkeyPatch.context() as mp:
        for target, name, value in patches:
            mp.setattr(target, name, value)
        code = synth_run_all(sample, out, **overrides)
    assert (
        code is None
    ), f"run-all exited {code}: an optional stage must not stop the run"
    return outcomes(out)


def _everything(sample, out, **overrides):
    """Every optional stage asked for: an assay, a bisulfite BAM, variants, TFBS and ATAC."""
    options = dict(
        assay="xs1",
        bisulfite_bam=sample.bam,
        variants=write_vcf(out.parent / f"{out.name}.vcf"),
        no_tfbs=False,
        no_atac=False,
    )
    options.update(overrides)
    return options


@pytest.fixture(scope="module")
def ran(sample, tmp_path_factory):
    out = tmp_path_factory.mktemp("ran") / "out"
    return _run(sample, out, **_everything(sample, out))


@pytest.fixture(scope="module")
def failed(sample, tmp_path_factory):
    root = tmp_path_factory.mktemp("failed")
    out = root / "out"
    pointer = lfs_pointer(root / "asset.bed.gz")
    # UXM reads a .bam only. (A file of junk named .bam does not fail it: the Rust reader
    # writes a table regardless, which is its own defect and not this one.)
    corrupt = root / "bisulfite.cram"
    corrupt.write_bytes(b"not a BAM file")
    patches = (
        # Region-MDS reads the gene BED the asset manager names; gene FSC its own lookup.
        (AssetManager, "get_gene_bed_for_mode", lambda self, assay=None: pointer),
        (gene_bed_module, "get_bundled_gene_bed", lambda *a, **k: pointer),
        (AssetManager, "get_wps_anchors", lambda self, assay=None: pointer),
        (AssetManager, "tfbs_regions", property(lambda self: pointer)),
        (AssetManager, "atac_regions", property(lambda self: pointer)),
    )
    result = _run(
        sample,
        out,
        patches,
        **_everything(
            sample,
            out,
            bisulfite_bam=corrupt,
            mfsd_bam=root / "absent.duplex.bam",
        ),
    )
    return result, out


@pytest.fixture(scope="module")
def unavailable(sample, tmp_path_factory):
    root = tmp_path_factory.mktemp("unavailable")
    out = root / "out"
    absent = root / "absent.bed.gz"
    patches = (
        (AssetManager, "get_gene_bed_for_mode", lambda self, assay=None: None),
        (AssetManager, "get_wps_anchors", lambda self, assay=None: absent),
        (AssetManager, "tfbs_regions", property(lambda self: absent)),
        (AssetManager, "atac_regions", property(lambda self: absent)),
    )
    return _run(
        sample,
        out,
        patches,
        **_everything(
            sample,
            out,
            bisulfite_bam=root / "absent.bisulfite.bam",
            variants=root / "absent.vcf",
        ),
    )


@pytest.fixture(scope="module")
def not_requested_run(sample, tmp_path_factory):
    out = tmp_path_factory.mktemp("not_requested") / "out"
    return _run(sample, out), out


@pytest.fixture(scope="module")
def not_requested(not_requested_run):
    return not_requested_run[0]


@pytest.mark.parametrize("stage", ("uxm", "mfsd", *ASSET_STAGES))
def test_a_stage_that_runs_records_ran_and_one_that_breaks_records_failed(
    stage, ran, failed
):
    failed_outcomes, _ = failed
    assert ran[stage] == ("ran", ""), stage
    status, error = failed_outcomes[stage]
    assert status == "failed", (stage, failed_outcomes[stage])
    assert error and error.isidentifier(), (stage, error)


def test_the_recorded_classes_are_the_real_ones(failed):
    """The classes, not `Exit`: the stage bodies raise, only the commands convert."""
    failed_outcomes, _ = failed
    assert failed_outcomes["mfsd"] == ("failed", "FileNotFoundError")
    assert failed_outcomes["uxm"] == ("failed", "ValueError")
    assert failed_outcomes["gene_fsc"] == ("failed", "BadGzipFile")
    assert all(failed_outcomes[s][1] != "Exit" for s in failed_outcomes)


def test_a_missing_mfsd_bam_fails_mfsd_instead_of_using_the_main_bam(failed):
    """It used to run on the primary BAM with duplex weighting off, and say nothing."""
    _, out = failed
    assert not (out / "S.mFSD.parquet").exists()


def test_a_stage_whose_input_is_missing_records_unavailable(unavailable):
    for stage in ("uxm", "mfsd", "region_mds", "panel_wps", "tfbs", "atac"):
        assert unavailable[stage] == ("unavailable", ""), stage
    assert unavailable["gene_fsc"] == ("ran", ""), "its gene BED was not withheld"


def test_a_stage_nobody_asked_for_records_not_requested(not_requested):
    for stage in ("uxm", "mfsd", *ASSET_STAGES, "region_mds_z"):
        assert not_requested[stage] == ("not_requested", ""), stage


def test_without_a_pon_neither_mds_z_score_was_requested(ran, not_requested):
    for result in (ran, not_requested):
        assert result["mds_z"] == ("not_requested", "")
        assert result["mds_z_ontarget"] == ("not_requested", "")
        assert result["region_mds_z"] == ("not_requested", "")


def test_a_variants_file_with_no_variants_records_no_variants(sample, tmp_path):
    out = tmp_path / "out"
    empty = write_vcf(tmp_path / "empty.vcf", n_variants=0)
    assert _run(sample, out, variants=empty)["mfsd"] == ("no_variants", "")


def test_a_recorded_error_never_contains_a_path(failed):
    """Sample directories are named for the patient (invariant #4)."""
    failed_outcomes, out = failed
    row = record(out)
    for stage in failed_outcomes:
        value = str(row[f"{stage}_error"])
        assert "/" not in value and str(out.parent) not in value, (stage, value)


def test_a_rerun_after_a_failure_records_ran(sample, tmp_path):
    out = tmp_path / "out"
    corrupt = tmp_path / "bisulfite.cram"
    corrupt.write_bytes(b"not a BAM file")
    assert _run(sample, out, bisulfite_bam=corrupt)["uxm"][0] == "failed"
    assert _run(sample, out, bisulfite_bam=sample.bam)["uxm"] == ("ran", "")


def test_a_single_tool_never_writes_not_requested_over_run_all(sample, tmp_path):
    """`krewlyzer fsc`, rerun with run-all's settings, did not attempt panel WPS: it must leave
    run-all's `ran` alone rather than record `not_requested`. (Rerun without the assay, it is
    refused before any work: `test_configuration_guard_runs.py`.)"""
    from krewlyzer.core.unified_processor import run_features

    out = tmp_path / "out"
    before = _run(sample, out, assay="xs1")
    assert before["gene_fsc"] == ("ran", "") and before["panel_wps"] == ("ran", "")
    run_features(
        bed_path=out / "S.bed.gz",
        output_dir=out,
        sample_name="S",
        assay="xs1",
        enable_fsc=True,
        fsc_bins=sample.assets.bins,
        gc_correct=False,
        output_format="parquet",
        threads=1,
        compress=False,
    )
    assert outcomes(out) == before


def test_a_single_region_mds_after_a_failed_run_all_records_ran(sample, tmp_path):
    """The single tool records its own outcome: rerunning it is how a failure gets fixed."""
    from typer.testing import CliRunner

    from krewlyzer.cli import app

    out = tmp_path / "out"
    pointer = lfs_pointer(tmp_path / "genes.bed.gz")
    patch = ((AssetManager, "get_gene_bed_for_mode", lambda self, assay=None: pointer),)
    assert _run(sample, out, patch, assay="xs1")["region_mds"][0] == "failed"

    result = CliRunner().invoke(
        app,
        [
            "region-mds",
            str(sample.bam),
            str(sample.reference),
            str(out),
            "--assay",
            "xs1",
            "--sample-name",
            "S",
            "--skip-pon",
            "--output-format",
            "parquet",
        ],
    )
    assert result.exit_code == 0, result.output
    assert outcomes(out)["region_mds"] == ("ran", "")


def _pon_without_mds_baseline(path: Path) -> Path:
    """A PON the loader accepts, with no MDS baseline, as a BED-built PON has none."""
    import pandas as pd

    pd.DataFrame(
        [
            {
                "table": "metadata",
                "schema_version": "1.0",
                "assay": "xs1",
                "n_samples": 21.0,
                "krewlyzer_version": "0.9.0",
            }
        ]
    ).to_parquet(path)
    return path


def test_a_pon_without_an_mds_baseline_gives_mds_z_unavailable(sample, tmp_path):
    """No MDS baseline is not a failed z-score."""
    pon = _pon_without_mds_baseline(tmp_path / "no_mds.pon.parquet")
    result = _run(sample, tmp_path / "out", pon_model=pon, skip_pon=False)
    assert result["mds_z"] == ("unavailable", "")
    assert result["region_mds_z"] == ("not_requested", ""), "no assay, so no Region-MDS"


def test_a_cached_motif_step_leaves_the_mds_z_record_alone(sample, tmp_path):
    """A rerun that reuses the motif tables (same PON) did not z-score them; their record
    stands, failures included.

    TSV, because the motif cache check looks for `.tsv` names only (noted in the plan, not
    fixed here). The earlier failure is written into the record by hand: it is what a rewrite
    from the reusing run would lose.
    """
    import pandas as pd

    out = tmp_path / "out"
    pon = _pon_without_mds_baseline(tmp_path / "no_mds.pon.parquet")
    _run_tsv(sample, out, pon_model=pon, skip_pon=False)
    path = out / "S.metadata.tsv"
    row = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    row["mds_z_status"], row["mds_z_error"] = "failed", "OSError"
    row.to_csv(path, sep="\t", index=False)
    _run_tsv(sample, out, pon_model=pon, skip_pon=False)
    assert _outcomes_tsv(out)["mds_z"] == ("failed", "OSError")
    assert _outcomes_tsv(out)["uxm"] == ("not_requested", ""), "the rest is rewritten"


def _run_tsv(sample, out, **overrides):
    assert synth_run_all(sample, out, output_format="tsv", **overrides) is None


def _outcomes_tsv(out: Path) -> dict:
    import pandas as pd

    row = pd.read_csv(
        out / "S.metadata.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    return {
        c[: -len("_status")]: (
            row[c].iloc[0],
            row[c.replace("_status", "_error")].iloc[0],
        )
        for c in row.columns
        if c.endswith("_status")
    }


# ---------------------------------------------------------------------------
# What run-all's own validation.json makes of the record
# ---------------------------------------------------------------------------


def _validation(out: Path) -> dict:
    import json

    return json.loads((out / "S.validation.json").read_text())


def _stage_findings(out: Path, finding_id: str) -> list:
    return [f for f in _validation(out)["findings"] if f["id"] == finding_id]


def test_run_all_reports_each_failed_stage_with_its_reason(failed):
    """At ERROR when the stage's tables are consumed, else WARN: under `--strict-validation`
    this is exactly today's outcome, with the reason instead of "absent"."""
    failed_outcomes, out = failed
    reported = {
        stage: f["severity"]
        for f in _stage_findings(out, "STAGE.FAILED")
        for stage in failed_outcomes
        if f["message"].startswith(f"{stage} ")
    }
    assert reported == {
        "uxm": "warn",
        "mfsd": "warn",
        "region_mds": "error",
        "gene_fsc": "error",
        "panel_wps": "error",
        "tfbs": "error",
        "atac": "error",
    }
    absent = {
        f["id"] for f in _validation(out)["findings"] if f["id"].endswith(".ABSENT")
    }
    assert not absent & {
        "MDS.gene.ABSENT",
        "FSC.gene.ABSENT",
        "TFBS.ABSENT",
        "mFSD.ABSENT",
    }


def test_a_recorded_error_never_reaches_validation_json_as_a_path(failed):
    """Checked where it is published, not only where it is written (invariant #4)."""
    _, out = failed
    text = (out / "S.validation.json").read_text()
    assert str(out.parent) not in text
    for finding in _stage_findings(out, "STAGE.FAILED"):
        assert "/" not in finding["message"], finding["message"]


def test_a_table_left_beside_a_failed_rerun_is_an_error(sample, tmp_path):
    """The rerun failed, so `MDS.gene` is the first run's: downstream would read it as this one."""
    out = tmp_path / "out"
    _run(sample, out, assay="xs1")
    pointer = lfs_pointer(tmp_path / "genes.bed.gz")
    patch = ((AssetManager, "get_gene_bed_for_mode", lambda self, assay=None: pointer),)
    _run(sample, out, patch, assay="xs1")
    stale = {f["table"] for f in _stage_findings(out, "STAGE.STALE_TABLE")}
    assert {".MDS.gene.parquet", ".MDS.exon.parquet"} <= stale
    assert all(
        f["severity"] == "error" for f in _stage_findings(out, "STAGE.STALE_TABLE")
    )


def test_no_mds_baseline_is_not_a_missing_z_score(sample, tmp_path):
    out = tmp_path / "out"
    pon = _pon_without_mds_baseline(tmp_path / "no_mds.pon.parquet")
    _run(sample, out, pon_model=pon, skip_pon=False)
    about_mds_z = [
        f
        for f in _validation(out)["findings"]
        if f["column"] == "mds_z" or "mds_z" in f["message"]
    ]
    assert not about_mds_z, about_mds_z


def test_uxm_and_mfsd_not_requested_give_no_finding(not_requested_run):
    _, out = not_requested_run
    tables = {".UXM.parquet", ".mFSD.parquet"}
    assert not [f for f in _validation(out)["findings"] if f["table"] in tables]


def test_a_single_region_entropy_after_a_failed_run_all_records_ran(sample, tmp_path):
    """`region-entropy` writes the TFBS and ATAC tables without `run_features`; rerunning it
    must clear run-all's `failed`, or the fixed tables read as stale."""
    from typer.testing import CliRunner

    from krewlyzer.cli import app

    out = tmp_path / "out"
    pointer = lfs_pointer(tmp_path / "regions.bed.gz")
    patch = ((AssetManager, "tfbs_regions", property(lambda self: pointer)),)
    assert _run(sample, out, patch, no_tfbs=False)["tfbs"][0] == "failed"

    result = CliRunner().invoke(
        app,
        [
            "region-entropy",
            "-i",
            str(out / "S.bed.gz"),
            "-o",
            str(out),
            "-s",
            "S",
            "--no-atac",
            "--skip-pon",
            "--no-gc-correct",
            "--output-format",
            "parquet",
        ],
    )
    assert result.exit_code == 0, result.output
    assert outcomes(out)["tfbs"] == ("ran", "")
    assert outcomes(out)["atac"] == ("not_requested", ""), "not touched: it was not run"


def test_a_single_tool_in_another_format_leaves_a_stale_parquet_table_flagged(
    sample, tmp_path
):
    """region-entropy at its default `-F tsv` rewrote `.TFBS.tsv`, and the `.TFBS.parquet`
    downstream reads stayed the first run's. It is now refused before writing (the folder's
    `table_formats`), and the gate must keep saying the Parquet table is stale."""
    from typer.testing import CliRunner

    from krewlyzer.cli import app
    from krewlyzer.validate.gate import check_sample

    out = tmp_path / "out"
    _run(sample, out, no_tfbs=False)  # TFBS.parquet written
    pointer = lfs_pointer(tmp_path / "regions.bed.gz")
    patch = ((AssetManager, "tfbs_regions", property(lambda self: pointer)),)
    assert _run(sample, out, patch, no_tfbs=False)["tfbs"][0] == "failed"

    result = CliRunner().invoke(
        app,
        [
            "region-entropy",
            "-i",
            str(out / "S.bed.gz"),
            "-o",
            str(out),
            "-s",
            "S",
            "--no-atac",
            "--skip-pon",
            "--no-gc-correct",
        ],
    )
    assert result.exit_code == 1, result.output
    assert "table_formats" in result.output
    assert not (out / "S.TFBS.tsv").exists()
    assert outcomes(out)["tfbs"][0] == "failed"
    stale = [f for f in check_sample("S", out)[0] if f.id == "STAGE.STALE_TABLE"]
    assert ".TFBS.parquet" in {f.table for f in stale}


def _mds_pon(path: Path) -> Path:
    """A PON the loader accepts that can score the whole-sample MDS."""
    import pandas as pd

    pd.DataFrame(
        [
            {
                "table": "metadata",
                "schema_version": "1.0",
                "assay": "xs1",
                "n_samples": 21.0,
                "krewlyzer_version": "0.9.0",
            },
            {"table": "mds_baseline", "mds_mean": 0.9, "mds_std": 0.05},
        ]
    ).to_parquet(path)
    return path


def test_cached_motif_tables_made_without_the_pon_now_in_use_are_redone(
    sample, tmp_path
):
    """The motif cache skipped the step whatever the PON: a run without one, then with one,
    shipped `pon_applied=True` with an unscored MDS."""
    import pandas as pd

    out = tmp_path / "out"
    _run_tsv(sample, out)
    assert _outcomes_tsv(out)["mds_z"] == ("not_requested", "")
    _run_tsv(
        sample, out, pon_model=_mds_pon(tmp_path / "mds.pon.parquet"), skip_pon=False
    )
    assert _outcomes_tsv(out)["mds_z"] == ("ran", "")
    assert "mds_z" in pd.read_csv(out / "S.MDS.tsv", sep="\t", comment="#").columns


def test_cached_motif_tables_are_redone_when_the_record_names_the_pon_but_mds_was_not_scored(
    sample, tmp_path
):
    """A single tool can raise `pon_applied` and name the PON (it scored its own tables), while
    the motif tables from the earlier run were never scored. A later run-all with that PON must
    not reuse them (review of #88 PR B)."""
    import pandas as pd

    from krewlyzer.core.output_provenance import record_pon_used
    from krewlyzer.core.pon_integration import load_requested_pon

    out = tmp_path / "out"
    _run_tsv(sample, out)
    pon_path = _mds_pon(tmp_path / "mds.pon.parquet")
    record_pon_used(
        out,
        "S",
        load_requested_pon(pon_path, output_format="tsv", compress=False),
        pon_path,
        output_format="tsv",
        compress=False,
    )
    assert _outcomes_tsv(out)["mds_z"] == ("not_requested", "")

    _run_tsv(sample, out, pon_model=pon_path, skip_pon=False)
    assert _outcomes_tsv(out)["mds_z"] == ("ran", "")
    assert "mds_z" in pd.read_csv(out / "S.MDS.tsv", sep="\t", comment="#").columns


# -- the motif cache, by format (Greptile on #98 review, F5) ---------------------------------
#
# The cache looked for `.tsv` names only, so a Parquet run-all into a folder a TSV run used
# skipped the motif step and wrote no Parquet motif tables. That sequence is now refused before
# anything is written (`test_a_tsv_run_all_then_a_parquet_run_all_is_refused`), so the cache only
# ever sees its own format.

_MOTIF_TABLES = ("EndMotif", "EndMotif1mer", "BreakPointMotif", "MDS")


def test_a_parquet_run_reuses_its_own_motif_tables(sample, tmp_path):
    """The other half: under Parquet the cache never fired at all."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out, output_format="parquet") is None
    before = {t: (out / f"S.{t}.parquet").stat().st_mtime_ns for t in _MOTIF_TABLES}
    assert synth_run_all(sample, out, output_format="parquet") is None
    after = {t: (out / f"S.{t}.parquet").stat().st_mtime_ns for t in _MOTIF_TABLES}
    assert after == before, "the motif step ran again with nothing changed"


def test_a_parquet_run_with_a_pon_does_not_reuse_unscored_motif_tables(
    sample, tmp_path
):
    """Guard: a format-aware cache must still check the PON."""
    import pandas as pd

    out = tmp_path / "out"
    assert synth_run_all(sample, out, output_format="parquet") is None
    pon_path = _mds_pon(tmp_path / "mds.pon.parquet")
    assert (
        synth_run_all(
            sample, out, output_format="parquet", pon_model=pon_path, skip_pon=False
        )
        is None
    )
    assert "mds_z" in pd.read_parquet(out / "S.MDS.parquet").columns


# -- a rewritten table owns its z-score stage (Greptile on #98, F1) ---------------------------


def _mds_z_findings(out: Path) -> list:
    """The gate run now, on the folder as it is: `S.validation.json` is what run-all wrote, before
    any single tool ran after it."""
    from krewlyzer.validate.gate import check_sample

    return [
        f.to_dict()
        for f in check_sample("S", out)[0]
        if f.column in ("mds_z", "mds_e1_z") or "mds_z" in f.message
    ]


def test_motif_without_a_pon_after_a_scored_run_records_its_mds_z_as_not_requested(
    sample, tmp_path
):
    """`motif --skip-pon` rewrites `MDS` with no `mds_z`. Its `not_requested` was dropped as
    "not attempted", so run-all's `mds_z=ran` stayed and the gate required the missing column.
    """
    out = tmp_path / "out"
    pon = _mds_pon(tmp_path / "mds.pon.parquet")
    assert synth_run_all(sample, out, pon_model=pon, skip_pon=False) is None
    assert outcomes(out)["mds_z"] == ("ran", "")

    from typer.testing import CliRunner

    from krewlyzer.cli import app

    result = CliRunner().invoke(
        app,
        [
            "motif",
            "-i",
            str(sample.bam),
            "--reference",
            str(sample.reference),
            "-o",
            str(out),
            "--sample-name",
            "S",
            "--skip-pon",
            "--skip-target-regions",
            "--output-format",
            "parquet",
        ],
    )
    assert result.exit_code == 0, result.output
    assert outcomes(out)["mds_z"] == ("not_requested", "")
    # The record is now true, and the gate still flags the table: the sample stays recorded as
    # PON-scored, and downstream would read its MDS so (review of the branch).
    flagged = [f for f in _mds_z_findings(out) if f["table"] == ".MDS.parquet"]
    assert flagged and all("not_requested" in f["message"] for f in flagged), flagged


def test_region_mds_without_a_pon_records_its_z_scores_as_not_requested(
    sample, tmp_path
):
    """The same for `region-mds --skip-pon` and `MDS.gene`."""
    from typer.testing import CliRunner

    from krewlyzer.cli import app
    from krewlyzer.core.output_provenance import Outcome, record_stage_outcomes

    out = tmp_path / "out"
    _run(sample, out, assay="xs1")
    # As a run-all scored against a PON with a Region-MDS baseline leaves it.
    record_stage_outcomes(
        out,
        "S",
        {"region_mds_z": Outcome("ran")},
        output_format="parquet",
        compress=False,
    )
    result = CliRunner().invoke(
        app,
        [
            "region-mds",
            str(sample.bam),
            str(sample.reference),
            str(out),
            "--assay",
            "xs1",
            "--sample-name",
            "S",
            "--skip-pon",
            "--output-format",
            "parquet",
        ],
    )
    assert result.exit_code == 0, result.output
    assert outcomes(out)["region_mds_z"] == ("not_requested", "")


def test_a_single_tool_still_leaves_a_table_stage_it_did_not_run_alone(
    sample, tmp_path
):
    """Guard: only a stage whose outputs a tool rewrote is recorded. `fsc` without an assay did
    not run gene FSC, whatever the record says about it."""
    from krewlyzer.core.output_provenance import Outcome, record_stage_outcomes
    from krewlyzer.core.unified_processor import run_features

    out = tmp_path / "out"
    _run(sample, out)
    record_stage_outcomes(
        out,
        "S",
        {"gene_fsc": Outcome.failed(RuntimeError("x"))},
        output_format="parquet",
        compress=False,
    )
    run_features(
        bed_path=out / "S.bed.gz",
        output_dir=out,
        sample_name="S",
        enable_fsc=True,
        fsc_bins=sample.assets.bins,
        gc_correct=False,
        threads=1,
        output_format="parquet",
        compress=False,
    )
    assert outcomes(out)["gene_fsc"] == ("failed", "RuntimeError")
