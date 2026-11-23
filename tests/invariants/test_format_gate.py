"""validate-output sees a sample folder holding tables in two formats (`FORMAT.MIXED`).

The configuration guard refuses a rerun in another format before it writes, but not under
Nextflow: each task starts in an empty folder, and `publishDir` then copies a rerun's files over an
`--outdir` that may hold an earlier format. The published folder is the only place that is
visible. A table in an extension outside the record's `table_formats` means a later run wrote here
in another format, so the Parquet downstream reads may be the earlier run's.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from conftest import requires_data
from krewlyzer.validate.findings import Severity
from krewlyzer.validate.gate import check_sample

from .synth_outputs import write_sample

pytestmark = pytest.mark.invariant


def _sample(tmp_path: Path, table_formats=None) -> Path:
    """A well-formed Parquet sample, its record naming ``table_formats`` (None: an older record)."""
    directory = write_sample(tmp_path, "S0", 0)
    if table_formats is not None:
        marker = directory / "S0.metadata.parquet"
        frame = pd.read_parquet(marker)
        frame["table_formats"] = table_formats
        frame.to_parquet(marker, index=False)
    return directory


def _add_tsv(directory: Path, name: str) -> None:
    pd.DataFrame({"x": [1.0, 2.0]}).to_csv(directory / name, sep="\t", index=False)


def _format_findings(directory: Path):
    return [f for f in check_sample("S0", directory)[0] if f.id.startswith("FORMAT.")]


def test_a_tsv_table_beside_a_parquet_record_is_an_error(tmp_path):
    directory = _sample(tmp_path, ".parquet")
    _add_tsv(directory, "S0.FSC.tsv")
    (finding,) = _format_findings(directory)
    assert finding.id == "FORMAT.MIXED" and finding.severity is Severity.ERROR
    assert "S0.FSC.tsv" in finding.evidence["files"]


def test_a_second_record_is_mixed(tmp_path):
    directory = _sample(tmp_path, ".parquet")
    _add_tsv(directory, "S0.metadata.tsv.gz")
    assert [f.id for f in _format_findings(directory)] == ["FORMAT.MIXED"]


def test_a_both_folder_is_not_mixed(tmp_path):
    """Guard: `both` writes every table, and the record, twice."""
    directory = _sample(tmp_path, ".parquet,.tsv")
    _add_tsv(directory, "S0.FSC.tsv")
    _add_tsv(directory, "S0.metadata.tsv")
    assert _format_findings(directory) == []


def test_a_file_that_is_no_table_is_not_mixed(tmp_path):
    """Guard: only tables the contract names count; BED, JSON and logs are not tables."""
    directory = _sample(tmp_path, ".parquet")
    (directory / "S0.bed.gz").write_bytes(b"")
    (directory / "S0.metadata.json").write_text("{}")
    _add_tsv(directory, "S0.fsc_counts.tsv")
    assert _format_findings(directory) == []


def test_an_older_record_with_tables_in_two_formats_is_not_judged(tmp_path):
    """Without the field a `both` folder and a Parquet-then-TSV folder look the same: say so,
    rather than pass it or fail it."""
    directory = _sample(tmp_path)
    _add_tsv(directory, "S0.FSC.tsv")
    (finding,) = _format_findings(directory)
    assert finding.id == "FORMAT.UNJUDGED" and finding.severity is Severity.SKIP


def test_an_older_record_in_one_format_gives_no_format_finding(tmp_path):
    """Guard: nothing to judge, so nothing to say about every 0.9.0 folder."""
    assert _format_findings(_sample(tmp_path)) == []


# A folder a real run wrote. The fixtures above are written by hand; run-all's own leftovers
# are what the gate meets in practice (review of the branch: mFSD and UXM left Rust's raw `.tsv`
# beside the `.parquet` they were converted to, and every Parquet sample running either failed).


@pytest.fixture(scope="module")
def synth(tmp_path_factory):
    from .synth import BASELINE, make_sample

    return make_sample(tmp_path_factory.mktemp("format_gate_in"), BASELINE)


@requires_data
@pytest.mark.parametrize(
    "output_format,compress", [("parquet", False), ("both", True), ("tsv", True)]
)
def test_a_run_all_with_every_stage_leaves_one_format(
    synth, tmp_path, output_format, compress
):
    from krewlyzer.core.output_utils import table_extensions

    from .runall import synth_run_all, write_vcf

    out = tmp_path / "out"
    assert (
        synth_run_all(
            synth,
            out,
            bisulfite_bam=synth.bam,
            variants=write_vcf(tmp_path / "v.vcf"),
            output_format=output_format,
            compress=compress,
        )
        is None
    )
    wanted = table_extensions(output_format=output_format, compress=compress)
    for table in ("mFSD", "UXM"):
        written = [
            e
            for e in (".parquet", ".tsv", ".tsv.gz")
            if (out / f"S.{table}{e}").exists()
        ]
        assert sorted(written) == sorted(wanted), (table, written)
    if output_format != "tsv":  # the gate reads Parquet folders
        assert _format_findings_of(out, "S") == []


def _format_findings_of(directory: Path, sample: str):
    return [f for f in check_sample(sample, directory)[0] if f.id.startswith("FORMAT.")]
