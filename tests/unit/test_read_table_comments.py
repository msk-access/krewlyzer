"""Regression tests for metadata-footer handling in read_table()."""

import pytest

from krewlyzer.core.output_utils import read_table


def test_read_table_ignores_hash_footer_lines(tmp_path):
    """'#'-prefixed footers must not be parsed as data rows.

    EndMotif1mer.tsv appends its summary metrics as '#' comment lines AFTER
    the data. Without comment='#', pandas turned each into a data row, and
    FeatureSerializer then emitted JSON keys '# c_fraction', '# entropy',
    '# c_bias' and '# sample' with NaN values inside features.motif.edm_1mer.
    """
    path = tmp_path / "sample.EndMotif1mer.tsv"
    path.write_text(
        "base\tcount\tfraction\n"
        "A\t661428\t0.188214\n"
        "C\t1261449\t0.358953\n"
        "G\t849433\t0.241712\n"
        "T\t741931\t0.211121\n"
        "# c_fraction\t0.358953\n"
        "# entropy\t1.952998\n"
        "# c_bias\t0.108953\n"
        "# sample\tsample\n"
    )

    df = read_table(path)

    assert df is not None
    assert len(df) == 4, f"expected only the 4 base rows, got:\n{df}"
    assert sorted(df["base"].tolist()) == ["A", "C", "G", "T"]
    assert not any(str(b).startswith("#") for b in df["base"]), (
        "comment footer leaked into the data frame"
    )
    # The numeric column must stay numeric (junk rows previously forced object dtype)
    assert df["fraction"].sum() > 0.99


def test_read_table_comment_can_be_overridden(tmp_path):
    """Callers may still opt out of comment handling explicitly."""
    path = tmp_path / "sample.tsv"
    path.write_text("a\tb\n1\t2\n# note\t3\n")

    df = read_table(path, comment=None)

    assert df is not None
    assert len(df) == 2, "comment=None should keep the '#' row"


def test_read_table_recovers_gzip_with_plain_text_footer(tmp_path):
    """Files written by krewlyzer <= 0.8.3 with --compress must stay readable.

    write_table() gzipped the table, then write_1mer_motifs() appended the
    '#' metadata footer with a PLAIN open(path, "a"). The result is a valid
    gzip member followed by raw bytes, which gzip.decompress() and pandas both
    reject outright -- so the whole EndMotif1mer table was unreadable, not
    merely polluted. Every production run using --compress is in this state.
    """
    import pandas as pd

    path = tmp_path / "sample.EndMotif1mer.tsv.gz"
    frame = pd.DataFrame(
        [
            {"base": "A", "count": 661428, "fraction": 0.188214},
            {"base": "C", "count": 1261449, "fraction": 0.358953},
            {"base": "G", "count": 849433, "fraction": 0.241712},
            {"base": "T", "count": 741931, "fraction": 0.211121},
        ]
    )
    frame.to_csv(path, sep="\t", index=False, compression="gzip")
    with open(path, "a") as fh:  # the historical bug, verbatim
        fh.write("# c_fraction\t0.358953\n")
        fh.write("# entropy\t1.952998\n")
        fh.write("# c_bias\t0.108953\n")
        fh.write("# sample\tsample\n")

    # Sanity: this really is unreadable by the normal path.
    with pytest.raises(Exception):
        pd.read_csv(path, sep="\t", compression="gzip", comment="#")

    df = read_table(path)

    assert df is not None, "read_table must recover the gzip member"
    assert len(df) == 4
    assert sorted(df["base"].tolist()) == ["A", "C", "G", "T"]
    assert not any(str(b).startswith("#") for b in df["base"])


def test_write_1mer_footer_is_written_through_gzip(tmp_path):
    """New writes must append the footer THROUGH the gzip codec."""
    import gzip as _gzip

    from krewlyzer.core.motif_processor import write_end_motif_1mer

    em_counts = {"AAAA": 10, "CCCC": 20, "GGGG": 10, "TTTT": 10}
    base = tmp_path / "sample.EndMotif1mer"
    write_end_motif_1mer(em_counts, base, output_format="tsv", compress=True,
                         sample_name="sample")

    path = tmp_path / "sample.EndMotif1mer.tsv.gz"
    assert path.exists()

    # The whole file must decompress as ONE valid gzip stream.
    text = _gzip.decompress(path.read_bytes()).decode()
    assert "# c_fraction" in text, "footer must be inside the gzip member"

    df = read_table(path)
    assert df is not None and len(df) == 4
