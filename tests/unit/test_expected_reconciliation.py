"""A sample that produced nothing is invisible, and there are two ways to vanish.

Everything else in the gate validates the samples it *finds*. Two failure modes
leave nothing to find, and both leave the cohort quietly smaller than intended:

* no output directory -- the job never ran, or died before creating one;
* a directory holding no Parquet -- ``discover_samples`` requires at least one
  ``{sample}.*.parquet``, so a job killed between mkdir and the first write is
  skipped exactly like a sample nobody submitted.

Only the third case, Parquet present but no completion marker, was ever caught.
At a handful of samples the other two are obvious; across 16,000 they are not,
and the consumer reads a short cohort without knowing it was short.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from krewlyzer.validate.cli import read_expected
from krewlyzer.validate.gate import reconcile_expected

pytestmark = pytest.mark.unit

# Deliberately NOT identifier-shaped. Only the two documented placeholders are
# allowlisted in .gitleaks.toml, so a set of *distinct* ids in the real patient
# format -- which is what this test first used, varying the final digit -- is
# precisely what the PHI guard exists to stop, and the pre-commit hook rejected
# them. Reconciliation is plain string matching, so a realistic shape buys this
# test nothing and costs it a guard violation.
SAMPLES = [f"sample_{i:02d}" for i in range(1, 6)]


def _ids(findings, finding_id):
    for f in findings:
        if f.id == finding_id:
            return f.samples
    return None


def test_a_sample_with_no_directory_is_reported(tmp_path: Path):
    (tmp_path / SAMPLES[0]).mkdir()
    findings = reconcile_expected(tmp_path, SAMPLES[:2], [SAMPLES[0]])
    assert _ids(findings, "EXPECTED.NO_OUTPUT_DIRECTORY") == [SAMPLES[1]]


def test_a_directory_with_no_tables_is_reported_separately(tmp_path: Path):
    """Distinct from "never ran" because the remedies differ.

    An empty directory means the job started and died, so its logs exist and
    the failure is diagnosable. No directory at all usually means the job was
    never submitted. Collapsing them would send someone hunting for logs that
    were never written.
    """
    (tmp_path / SAMPLES[0]).mkdir()  # exists, discovery found nothing in it
    findings = reconcile_expected(tmp_path, SAMPLES[:1], [])
    assert _ids(findings, "EXPECTED.DIRECTORY_HAS_NO_TABLES") == [SAMPLES[0]]
    assert _ids(findings, "EXPECTED.NO_OUTPUT_DIRECTORY") is None


def test_an_unexpected_sample_warns_rather_than_fails(tmp_path: Path):
    """A WARN, because it does not corrupt anything -- but it must be visible.

    The usual cause is an identifier that does not round-trip between the
    samplesheet and the pipeline, and that makes every other count here wrong.
    """
    findings = reconcile_expected(tmp_path, SAMPLES[:1], [SAMPLES[0], "sample_99"])
    unexpected = [f for f in findings if f.id == "EXPECTED.UNEXPECTED_SAMPLE"]
    assert unexpected and unexpected[0].severity.value == "warn"
    assert unexpected[0].samples == ["sample_99"]


def test_a_complete_cohort_reports_nothing(tmp_path: Path):
    assert reconcile_expected(tmp_path, SAMPLES, list(SAMPLES)) == []


def test_duplicate_ids_do_not_inflate_the_denominator(tmp_path: Path):
    """A samplesheet listing a sample twice must not report it missing twice."""
    findings = reconcile_expected(tmp_path, [SAMPLES[0], SAMPLES[0]], [])
    missing = [f for f in findings if f.id == "EXPECTED.NO_OUTPUT_DIRECTORY"]
    assert missing[0].samples == [SAMPLES[0]]
    assert missing[0].evidence["expected"] == 1


def test_reads_the_sample_column_of_a_samplesheet(tmp_path: Path):
    sheet = tmp_path / "samplesheet.csv"
    sheet.write_text(
        "sample,bam,mfsd_bam,assay\n"
        + "".join(f"{s},/x/{s}.bam,/x/{s}.d.bam,xs1\n" for s in SAMPLES)
    )
    assert read_expected(sheet) == SAMPLES


def test_reads_a_plain_list_ignoring_comments(tmp_path: Path):
    listing = tmp_path / "ids.txt"
    listing.write_text("# cohort\n" + "\n".join(SAMPLES) + "\n\n")
    assert read_expected(listing) == SAMPLES


def test_a_csv_without_a_sample_column_is_an_error(tmp_path: Path):
    """Never guess at the first column.

    Reconciling against the wrong field reports every sample missing, which
    reads as a catastrophic run failure rather than a malformed input -- the
    most expensive possible way to be wrong here.
    """
    sheet = tmp_path / "wrong.csv"
    sheet.write_text("sample_id,bam\nsample_01,/x/a.bam\n")
    with pytest.raises(Exception, match="no 'sample' column"):
        read_expected(sheet)
