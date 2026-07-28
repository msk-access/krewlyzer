"""Regression tests for metadata-footer handling in read_table()."""

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
