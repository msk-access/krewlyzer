"""`krewlyzer describe-output` — shape, content, and what it must never print.

The report exists to answer "what are these files?" without opening them, and
is meant to be hosted. That second part is why the redaction tests matter more
than the formatting ones.

Sample directories in this project are named for the patient, and several
tables carry the sample id as a *column value*. The first version of this
report was generated from deliberately renamed files and still leaked a real
identifier through the `Sample` and `sample_id` example values — the filenames
were clean and the contents were not.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

from krewlyzer.validate.describe import (
    REDACTED,
    _column_facts,
    describe_sample,
    human_bytes,
    is_identifier_column,
    render_markdown,
)

SAMPLE = "P-0000000-T01-XS1"

#: A *different* placeholder, written into `metadata.sample_id`.
#:
#: Three constraints meet here. It must differ from the directory name, or the
#: redaction test passes merely because the two coincide. It must be
#: identifier-shaped, or the leak regex below matches nothing and the test is
#: vacuous. And it must use the sanctioned `P-0000000` prefix -- a realistic
#: fake was the first thing tried, and the PHI pre-commit hook rejected it,
#: correctly: an invented-but-plausible id is indistinguishable from a real one
#: to everyone downstream of the person who invented it.
OTHER_SAMPLE = "P-0000000-T09-XS2"

#: Shape of an MSK-ACCESS identifier, excluding the all-zero placeholder.
#:
#: Written as a negative lookahead rather than by filtering matches afterwards,
#: so the pattern itself encodes "a real one", and this file never has to
#: contain an example of one to test against.
_IDENTIFIER = re.compile(r"P-(?!0000000)0\d{6}-[A-Z]\d{2}|C-[0-9A-Z]{6}-[A-Z]\d{3}")


@pytest.fixture
def sample_dir(tmp_path: Path) -> Path:
    """A minimal but realistic output directory.

    `Sample` and `sample_id` deliberately carry a *different* placeholder from
    the directory name, so a test asserting redaction cannot pass merely
    because the two happen to match.
    """
    d = tmp_path / SAMPLE
    d.mkdir()
    pd.DataFrame(
        {
            "region": ["chr1:0-100000", "chr2:0-100000"],
            "total": [1234.5, 2345.6],
            "65-69": [1.0, 2.0],
        }
    ).to_parquet(d / f"{SAMPLE}.FSD.parquet")
    pd.DataFrame(
        {
            "Motif": ["AAAA", "CCCC"],
            "Frequency": [0.25, 0.75],
        }
    ).to_parquet(d / f"{SAMPLE}.EndMotif.parquet")
    pd.DataFrame(
        {
            "sample_id": [OTHER_SAMPLE],
            "total_fragments": [8_200_000],
        }
    ).to_parquet(d / f"{SAMPLE}.metadata.parquet")
    return d


def test_it_finds_the_tables_that_are_there(sample_dir):
    report = describe_sample(sample_dir)
    families = {t.family for t in report.present}
    assert {"FSD", "EndMotif", "metadata"} <= families


def test_absent_tables_are_reported_not_skipped(sample_dir):
    """A missing table is the interesting case; silence would hide it."""
    report = describe_sample(sample_dir)
    assert report.missing, "nothing reported missing from a 3-table directory"
    assert "FSC" in {t.family for t in report.missing}


def test_shape_is_measured_from_the_file(sample_dir):
    report = describe_sample(sample_dir)
    fsd = next(t for t in report.present if t.family == "FSD")
    assert (fsd.n_rows, fsd.n_cols) == (2, 3)
    assert fsd.size_bytes and fsd.size_bytes > 0


def test_numeric_ranges_are_reported(sample_dir):
    report = describe_sample(sample_dir)
    fsd = next(t for t in report.present if t.family == "FSD")
    total = next(c for c in fsd.columns if c.name == "total")
    assert total.minimum == pytest.approx(1234.5)
    assert total.maximum == pytest.approx(2345.6)


def test_consumption_status_comes_from_the_contract(sample_dir):
    """Not a hand-kept list: it follows `NOT_CONSUMED` automatically."""
    report = describe_sample(sample_dir)
    assert next(t for t in report.present if t.family == "FSD").consumed is True


# ---------------------------------------------------------------------------
# Redaction -- the tests that decide whether this is safe to host
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name", ["sample", "Sample", "sample_id", "SAMPLE_ID", "patient", "subject", "id"]
)
def test_identifier_columns_are_recognised(name):
    assert is_identifier_column(name)


@pytest.mark.parametrize("name", ["total", "Frequency", "mds_e1", "region", "chrom"])
def test_measurement_columns_are_not_redacted(name):
    """Over-redaction would gut the report; the set is deliberately narrow."""
    assert not is_identifier_column(name)


def test_an_identifier_in_a_column_value_never_reaches_the_report(sample_dir):
    """The exact leak. Filenames were clean; the contents were not.

    `metadata.parquet` carries a sample id as a value, and printing it as an
    "example" put a real patient identifier into a document intended for
    hosting.
    """
    text = render_markdown(describe_sample(sample_dir))
    assert OTHER_SAMPLE not in text, (
        f"{OTHER_SAMPLE} reached the report from metadata.sample_id -- the "
        "exact leak this redaction exists to prevent"
    )
    leaked = _IDENTIFIER.findall(text)
    assert not leaked, f"identifier(s) leaked into the report: {sorted(set(leaked))}"
    assert REDACTED in text, "the identifier column should be shown as redacted"


def test_the_column_is_still_described_after_redaction(sample_dir):
    """Redacting the value must not hide that the column exists.

    Knowing `metadata` carries a `sample_id` is the useful fact; which id it is
    is the part that must not leave the machine.
    """
    report = describe_sample(sample_dir)
    meta = next(t for t in report.present if t.family == "metadata")
    col = next(c for c in meta.columns if c.name == "sample_id")
    assert col.example == REDACTED
    assert col.n_distinct == 1 and col.n_null == 0


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def test_the_report_names_the_completion_marker_when_absent(tmp_path):
    """A sample without it is dropped from the cohort silently, so the report
    is the last place it can be said."""
    d = tmp_path / SAMPLE
    d.mkdir()
    pd.DataFrame({"Motif": ["AAAA"], "Frequency": [1.0]}).to_parquet(
        d / f"{SAMPLE}.EndMotif.parquet"
    )
    text = render_markdown(describe_sample(d))
    assert "metadata.parquet" in text and "silently" in text


def test_an_empty_directory_produces_a_report_rather_than_an_error(tmp_path):
    d = tmp_path / SAMPLE
    d.mkdir()
    report = describe_sample(d)
    assert report.present == []
    assert "0 of" in render_markdown(report)


@pytest.mark.parametrize(
    "n,expected", [(None, "—"), (0, "—"), (512, "512 B"), (2048, "2.0 KB")]
)
def test_human_bytes(n, expected):
    assert human_bytes(n) == expected


def test_a_numeric_identifier_column_is_redacted_too():
    """Redacting `example` alone was not enough.

    The renderer prefers a min…max range over the example whenever one exists,
    so a numeric identifier — an accession, an integer patient key — was printed
    in full by the very branch meant to hide it. A range over one distinct value
    *is* the value.
    """
    facts = _column_facts(pd.Series([88123456, 88123456]), "patient_id", {})
    assert facts.example == REDACTED
    assert facts.minimum is None and facts.maximum is None


def test_redaction_does_not_cost_the_column_its_shape():
    """Knowing a column holds an identifier is the useful fact; which one is not."""
    facts = _column_facts(pd.Series([1, 2, 3]), "sample_id", {})
    assert facts.dtype.startswith("int")
    assert facts.n_distinct == 3
    assert facts.n_null == 0


def test_a_numeric_measurement_still_reports_its_range():
    """The guard is on the name, not on being numeric — ranges are the point."""
    facts = _column_facts(pd.Series([0.1, 0.9]), "short_long_log2", {})
    assert facts.minimum == pytest.approx(0.1)
    assert facts.maximum == pytest.approx(0.9)
