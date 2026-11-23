"""A rerun with other settings is refused before it writes anything (Greptile on #98, F2).

The leftover Greptile found needs a second run into a folder an earlier run used, with settings
that write fewer tables. Each entry point that writes a settings-dependent set of tables checks
the folder's record first: run-all, `extract` (which rewrites the record), the `run_features`
tools, `motif`, `region-entropy`, `region-mds`, `uxm` and `mfsd`. The output format is one such
setting for every one of them. The comparison rules are pinned in
`tests/unit/test_configuration_guard.py`; these pin who checks, and that nothing is written.
"""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

from conftest import requires_data
from krewlyzer.cli import app

from .runall import record, synth_run_all, write_vcf
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


def test_uxm_is_not_refused_for_the_assay_or_targets(sample, tmp_path):
    """Guard: its tables do not depend on assay or targets. At the folder's own format: any tool
    is refused another one (`test_a_tool_without_the_format_into_a_parquet_folder_is_refused`).
    """
    out = tmp_path / "out"
    assert synth_run_all(sample, out, assay="xs1") is None
    result = _cli(
        "uxm", "-i", sample.bam, "-o", out, "--sample-name", "S", "-F", "parquet"
    )
    assert result.exit_code == 0, result.output


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
    """`extract` does not take E1. When it rewrote run-all's record from scratch, dropping
    `fsc_e1_aggregation` let a later run-all with the other setting pass the check, and the first
    run's E1 table read as current (review of the branch). It now refuses run-all's record
    outright (`test_extract_record.py`), which keeps the setting; this pins that the later run-all
    is still refused."""
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
    assert result.exit_code == 1, result.output
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


# The output format. Every CLI defaults to `-F tsv` and Nextflow writes Parquet, so a tool run by
# hand without `-F` into a Nextflow folder rewrote its tables as TSV and left the Parquet ones,
# which downstream reads, from the earlier run.


def _fsc(out, sample, *format_args):
    return _cli(
        "fsc",
        "-i",
        out / "S.bed.gz",
        "-o",
        out,
        "-s",
        "S",
        "--skip-target-regions",
        "--skip-pon",
        "-b",
        sample.assets.bins,
        *format_args,
    )


def test_a_tool_without_the_format_into_a_parquet_folder_is_refused(sample, tmp_path):
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    before = _files(out)
    result = _fsc(out, sample)
    assert result.exit_code == 1, result.output
    assert "table_formats" in result.output and "Traceback" not in result.output
    assert "--output-format parquet" in " ".join(result.output.split())
    assert _files(out) == before, "a refused run wrote into the folder"


def test_a_tool_at_the_folder_format_is_allowed(sample, tmp_path):
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    result = _fsc(out, sample, "--output-format", "parquet")
    assert result.exit_code == 0, result.output


def test_a_both_folder_allows_a_both_tool(sample, tmp_path):
    """Guard: `both` writes two records, which are not by themselves a mixed folder."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out, output_format="both") is None
    result = _fsc(out, sample, "--output-format", "both")
    assert result.exit_code == 0, result.output


@pytest.mark.parametrize(
    "fmt",
    [["both"], ["both", "--compress"], ["parquet"]],
    ids=["both", "both-compress", "parquet"],
)
def test_single_tools_after_run_all_in_its_format_are_allowed_and_not_mixed(
    sample, tmp_path, fmt
):
    """A single tool records its stage in each copy of the record, with the run's formats.

    Recording each copy's own format instead (the `.parquet` copy of a `both` record says
    `parquet`) passed every other test. Yet after one `mfsd -F both`, every later `both` run was
    refused as mixed, and `validate-output` reported a clean folder `FORMAT.MIXED`.
    """
    from krewlyzer.validate.gate import check_sample

    out = tmp_path / "out"
    vcf = write_vcf(tmp_path / "v.vcf")
    output_format, compress = fmt[0], "--compress" in fmt
    flags = ["--output-format", *fmt]
    runs = dict(variants=vcf, bisulfite_bam=sample.bam)
    assert (
        synth_run_all(
            sample, out, output_format=output_format, compress=compress, **runs
        )
        is None
    )
    bed = out / "S.bed.gz"
    for args in (
        [
            "mfsd",
            "-i",
            sample.bam,
            "-V",
            vcf,
            "-g",
            sample.reference,
            "--sample-name",
            "S",
        ],
        ["uxm", "-i", sample.bam, "--sample-name", "S"],
        ["fsd", "-i", bed, "-s", "S", "--skip-target-regions", "--skip-pon"],
        [
            "mfsd",
            "-i",
            sample.bam,
            "-V",
            vcf,
            "-g",
            sample.reference,
            "--sample-name",
            "S",
        ],
    ):
        result = _cli(*args, "-o", out, *flags)
        assert result.exit_code == 0, f"{args[0]} refused:\n{result.output[-1500:]}"
    assert (
        synth_run_all(
            sample, out, output_format=output_format, compress=compress, **runs
        )
        is None
    )
    findings, _ = check_sample("S", out)
    assert not [f.id for f in findings if f.id.startswith("FORMAT")]


def test_a_tsv_run_all_then_a_parquet_run_all_is_refused(sample, tmp_path):
    out = tmp_path / "out"
    assert synth_run_all(sample, out, output_format="tsv") is None
    before = _files(out)
    assert synth_run_all(sample, out) == 1
    assert _files(out) == before


def test_mfsd_without_the_format_into_a_parquet_folder_is_refused(sample, tmp_path):
    """`mfsd` was no configuration check site: its tables depend on no assay or targets. Their
    format still decides which copy downstream reads."""
    out = tmp_path / "out"
    assert synth_run_all(sample, out) is None
    vcf = write_vcf(tmp_path / "variants.vcf")
    before = _files(out)
    result = _cli(
        "mfsd",
        "-i",
        sample.bam,
        "-V",
        vcf,
        "-o",
        out,
        "-s",
        "S",
        "-g",
        sample.reference,
    )
    assert result.exit_code == 1, result.output
    assert "table_formats" in result.output and "Traceback" not in result.output
    assert _files(out) == before


def test_the_record_says_the_format_it_was_written_in(sample, tmp_path):
    """Invariant #1: two formats, two recorded values."""
    assert synth_run_all(sample, tmp_path / "p") is None
    assert synth_run_all(sample, tmp_path / "b", output_format="both") is None
    assert record(tmp_path / "p")["table_formats"] == ".parquet"
    assert record(tmp_path / "b")["table_formats"] == ".parquet,.tsv"


def test_extract_records_its_own_format(sample, tmp_path):
    """`extract` writes the record from scratch and takes `-F`: it decides this one itself."""
    import pandas as pd

    out = tmp_path / "out"
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
        "--output-format",
        "tsv",
        "--compress",
    )
    assert result.exit_code == 0, result.output
    row = pd.read_csv(out / "S.metadata.tsv.gz", sep="\t").iloc[0]
    assert row["table_formats"] == ".tsv.gz"
