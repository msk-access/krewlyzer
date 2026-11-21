"""A rerun with other settings is refused before it writes anything (Greptile on #98, F2).

The leftover Greptile found needs a second run into a folder an earlier run used, with settings
that write fewer tables. Each entry point that writes a settings-dependent set of tables checks
the folder's record first: run-all, `extract` (which rewrites the record), the `run_features`
tools, `motif`, `region-entropy` and `region-mds`. The comparison rules are pinned in
`tests/unit/test_configuration_guard.py`; these pin who checks, and that nothing is written.
"""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

from conftest import requires_data
from krewlyzer.cli import app

from .runall import synth_run_all
from .synth import BASELINE, make_sample

pytestmark = [pytest.mark.invariant, requires_data]


@pytest.fixture(scope="module")
def sample(tmp_path_factory):
    return make_sample(tmp_path_factory.mktemp("config_in"), BASELINE)


def _files(out):
    return {p.name: p.stat().st_mtime_ns for p in sorted(out.iterdir()) if p.is_file()}


def _cli(*args):
    return CliRunner().invoke(app, [str(a) for a in args])


def test_run_all_without_the_assay_into_an_assay_folder_is_refused(sample, tmp_path):
    out = tmp_path / "out"
    assert synth_run_all(sample, out, assay="xs1") is None
    before = _files(out)
    assert synth_run_all(sample, out) == 1
    assert _files(out) == before, "a refused run wrote into the folder"


def test_run_all_with_the_same_settings_is_allowed(sample, tmp_path):
    out = tmp_path / "out"
    assert synth_run_all(sample, out, assay="xs1") is None
    assert synth_run_all(sample, out, assay="xs1") is None


def test_extract_without_the_assay_into_an_assay_folder_is_refused(sample, tmp_path):
    """`extract` writes the record from scratch: unchecked, it would replace run-all's assay with
    its own, and a later run-all without the assay would then pass."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out, assay="xs1") is None
    before = _files(out)
    result = _cli(
        "extract",
        "-i",
        sample.bam,
        "--reference",
        sample.reference,
        "--output",
        out,
        "--sample-name",
        "S",
        "--skip-target-regions",
    )
    assert result.exit_code == 1, result.output
    assert "assay" in result.output and "Traceback" not in result.output
    assert _files(out) == before


def test_region_entropy_without_the_assay_into_an_assay_folder_is_refused(
    sample, tmp_path
):
    out = tmp_path / "out"
    assert synth_run_all(sample, out, assay="xs1") is None
    before = _files(out)
    result = _cli(
        "region-entropy",
        "-i",
        out / "S.bed.gz",
        "-o",
        out,
        "-s",
        "S",
        "--skip-pon",
        "--no-gc-correct",
        "--output-format",
        "parquet",
    )
    assert result.exit_code == 1, result.output
    assert "assay" in result.output and "Traceback" not in result.output
    assert _files(out) == before


def test_fsc_without_the_e1_setting_of_the_run_all_is_refused(sample, tmp_path):
    out = tmp_path / "out"
    assert synth_run_all(sample, out, assay="xs1", disable_e1_aggregation=True) is None
    before = _files(out)
    result = _cli(
        "fsc",
        "-i",
        out / "S.bed.gz",
        "-o",
        out,
        "-s",
        "S",
        "--assay",
        "xs1",
        "--skip-target-regions",
        "--skip-pon",
        "-b",
        sample.assets.bins,
        "--output-format",
        "parquet",
    )
    assert result.exit_code == 1, result.output
    assert "fsc_e1_aggregation" in result.output and "Traceback" not in result.output
    assert _files(out) == before


def test_a_tool_that_writes_no_gene_fsc_ignores_the_e1_setting(sample, tmp_path):
    """Guard: E1 decides a table only gene FSC writes."""
    from krewlyzer.core.unified_processor import run_features

    out = tmp_path / "out"
    assert synth_run_all(sample, out, assay="xs1", disable_e1_aggregation=True) is None
    run_features(
        bed_path=out / "S.bed.gz",
        output_dir=out,
        sample_name="S",
        assay="xs1",
        enable_fsd=True,
        fsd_arms=sample.directory / "arms.validated.bed",
        gc_correct=False,
        threads=1,
        output_format="parquet",
        compress=False,
    )


def test_uxm_is_never_refused(sample, tmp_path):
    """Guard: its tables do not depend on assay or targets."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out, assay="xs1") is None
    result = _cli("uxm", "-i", sample.bam, "-o", out, "--sample-name", "S")
    assert "fresh output folder" not in result.output


def test_the_in_process_pon_build_is_never_refused(sample, tmp_path):
    """Guard: it owns the record it writes, as run-all does."""
    from krewlyzer.core.unified_processor import run_features

    out = tmp_path / "out"
    assert synth_run_all(sample, out, assay="xs1") is None
    run_features(
        bed_path=out / "S.bed.gz",
        output_dir=out,
        sample_name="S",
        owns_sample_record=True,
        enable_fsc=True,
        fsc_bins=sample.assets.bins,
        gc_correct=False,
        threads=1,
        output_format="parquet",
        compress=False,
    )


def test_extract_keeps_the_settings_it_does_not_decide(sample, tmp_path):
    """`extract` rewrites the record from scratch, and does not take E1. Dropping run-all's
    `fsc_e1_aggregation` let a later run-all with the other setting pass the check, and the first
    run's E1 table read as current (review of the branch)."""
    import pandas as pd

    out = tmp_path / "out"
    assert synth_run_all(sample, out, assay="xs1") is None
    result = _cli(
        "extract",
        "-i",
        sample.bam,
        "--reference",
        sample.reference,
        "--output",
        out,
        "--sample-name",
        "S",
        "--assay",
        "xs1",
        "--skip-target-regions",
        "--output-format",
        "parquet",
    )
    assert result.exit_code == 0, result.output
    assert bool(
        pd.read_parquet(out / "S.metadata.parquet")["fsc_e1_aggregation"].iloc[0]
    )
    assert synth_run_all(sample, out, assay="xs1", disable_e1_aggregation=True) == 1


def test_uxm_with_another_genome_is_refused(sample, tmp_path):
    """UXM chooses its marker file by genome: an hg38 run into an hg19 folder would compute UXM
    over the wrong markers and record it as `ran` (review of the branch)."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    before = _files(out)
    result = _cli(
        "uxm", "-i", sample.bam, "-o", out, "--sample-name", "S", "-G", "hg38"
    )
    assert result.exit_code == 1, result.output
    assert "genome" in result.output and "Traceback" not in result.output
    assert _files(out) == before


def test_e1_aggregation_does_not_matter_without_gene_fsc(sample, tmp_path):
    """Guard: without an assay no gene FSC is written, so E1 decides no table. run-all compared it
    anyway and refused a WGS rerun that differed only in it; `run_features` did not (review of the
    branch)."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    assert synth_run_all(sample, out, disable_e1_aggregation=True) is None
