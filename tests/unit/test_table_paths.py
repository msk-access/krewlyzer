"""`table_paths`: the files a run of a given format writes, exactly. One statement of it, for the
record readers, the stage-outcome writer, the motif tables and the run-all motif cache, which each
spelled it out before."""

from pathlib import Path

import pytest

from krewlyzer.core.output_utils import table_paths

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    "fmt,compress,names",
    [
        ("parquet", False, ["S.EndMotif.parquet"]),
        ("parquet", True, ["S.EndMotif.parquet"]),
        ("tsv", False, ["S.EndMotif.tsv"]),
        ("tsv", True, ["S.EndMotif.tsv.gz"]),
        ("both", False, ["S.EndMotif.parquet", "S.EndMotif.tsv"]),
        ("both", True, ["S.EndMotif.parquet", "S.EndMotif.tsv.gz"]),
    ],
)
def test_the_paths_a_format_writes(fmt, compress, names):
    paths = table_paths(Path("out/S.EndMotif"), output_format=fmt, compress=compress)
    assert [p.name for p in paths] == names
    assert all(p.parent == Path("out") for p in paths)


def test_an_unknown_format_is_refused():
    with pytest.raises(ValueError, match="unknown output format"):
        table_paths(Path("S.MDS"), output_format="csv", compress=False)
