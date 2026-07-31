"""E1 selection in `filter_fsc_to_e1`.

It used to sort by ``start`` and take the first row per gene: no strand input
existed, so it returned the *last* exon of every minus-strand gene, and it
emitted a row for every gene whether or not any of its regions were near a
transcription start.

`FSC.regions` now carries `strand`, `is_e1` and `is_alt_e1`, resolved at build
time from a GENCODE transcript. Selection is on the flags, and a gene with
neither is **omitted** rather than back-filled with an internal exon.

Both fallbacks are covered too. They exist so a legacy or custom gene BED still
produces a table, and they are the paths a real run will silently take if the
asset is stale — so they need to behave predictably and say what they did.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import pandas as pd

from krewlyzer.core.fsc_processor import filter_fsc_to_e1

#: Column order as written by `write_region_output` in rust/src/fsc.rs.
_BASE = ["chrom", "start", "end", "gene", "region_name", "region_bp"]
_FLAGS = ["strand", "is_e1", "is_alt_e1"]
_COUNTS = ["total", "normalized_depth"]


def _regions(rows: List[dict], with_flags: bool = True) -> pd.DataFrame:
    cols = _BASE + (_FLAGS if with_flags else []) + _COUNTS
    return pd.DataFrame(rows)[cols]


def _row(
    gene: str,
    start: int,
    name: str,
    strand: str = "+",
    is_e1: str = "0",
    is_alt_e1: str = "0",
) -> dict:
    return {
        "chrom": "1",
        "start": start,
        "end": start + 100,
        "gene": gene,
        "region_name": name,
        "region_bp": 100,
        "strand": strand,
        "is_e1": is_e1,
        "is_alt_e1": is_alt_e1,
        "total": 50.0,
        "normalized_depth": 1.0,
    }


def _run(tmp_path: Path, df: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Run the filter and read back what it wrote.

    ``output_path`` is deliberately not passed: the function derives the name,
    and supplying one containing a dot trips the ``Path.with_suffix`` trap its
    own comments warn about (``.e1only`` is stripped as if it were a suffix).
    Letting it derive the name is also what the pipeline does.
    """
    src = tmp_path / "S1.FSC.regions.tsv"
    df.to_csv(src, sep="\t", index=False)
    out = filter_fsc_to_e1(src)
    if out is None:
        return None
    written = next(p for p in tmp_path.glob("*e1only*") if p.suffix == ".tsv")
    return pd.read_csv(written, sep="\t")


def test_a_canonical_e1_region_is_selected(tmp_path):
    df = _regions(
        [
            _row("TP53", 100, "a"),
            _row("TP53", 900, "b", is_e1="1"),
        ]
    )
    out = _run(tmp_path, df)
    assert list(out["region_name"]) == [
        "b"
    ], "selection fell back to the lowest start instead of using is_e1"


def test_an_alternative_promoter_region_also_qualifies(tmp_path):
    """The decision this implements: `is_e1 OR is_alt_e1`.

    Restricting to the canonical transcript discards most of the panel — 25 of
    128 xs1 genes against 40 — because alternative promoters are the norm.
    """
    df = _regions(
        [
            _row("MYC", 100, "a"),
            _row("MYC", 900, "b", is_alt_e1="1"),
        ]
    )
    out = _run(tmp_path, df)
    assert list(out["region_name"]) == ["b"]


def test_a_gene_with_neither_flag_is_omitted(tmp_path):
    """Not back-filled. Its most 5' region is an internal exon, and calling
    that E1 asserts promoter proximity that is not there."""
    df = _regions(
        [
            _row("TP53", 100, "a", is_e1="1"),
            _row("NOFLAG", 100, "x"),
            _row("NOFLAG", 200, "y"),
        ]
    )
    out = _run(tmp_path, df)
    assert set(out["gene"]) == {"TP53"}, (
        "a gene with no annotated first exon must be dropped, not represented "
        "by an arbitrary internal exon"
    )


def test_minus_strand_selection_does_not_depend_on_coordinate_order(tmp_path):
    """The original defect. The flagged region is at the highest coordinate."""
    df = _regions(
        [
            _row("BRCA1", 100, "last", strand="-"),
            _row("BRCA1", 500, "mid", strand="-"),
            _row("BRCA1", 900, "first", strand="-", is_e1="1"),
        ]
    )
    out = _run(tmp_path, df)
    assert list(out["region_name"]) == ["first"]


def test_one_row_per_gene_when_several_regions_are_flagged(tmp_path):
    """A wide tile can overlap more than one annotated first exon."""
    df = _regions(
        [
            _row("EGFR", 100, "a", is_e1="1"),
            _row("EGFR", 300, "b", is_alt_e1="1"),
        ]
    )
    out = _run(tmp_path, df)
    assert len(out) == 1 and out["gene"].iloc[0] == "EGFR"


def test_flag_columns_survive_into_the_output(tmp_path):
    """A consumer must be able to tell canonical from alternative."""
    df = _regions([_row("TP53", 100, "a", is_alt_e1="1")])
    out = _run(tmp_path, df)
    for col in _FLAGS:
        assert col in out.columns, f"{col} was dropped from the E1 output"
    assert str(out["is_alt_e1"].iloc[0]) == "1"


# ---------------------------------------------------------------------------
# Fallbacks -- the paths a stale asset takes
# ---------------------------------------------------------------------------


def test_without_flag_columns_it_falls_back_and_warns(tmp_path, caplog):
    """A legacy gene BED still yields a table, loudly labelled as a guess."""
    df = _regions([_row("TP53", 900, "b"), _row("TP53", 100, "a")], with_flags=False)
    with caplog.at_level("WARNING"):
        out = _run(tmp_path, df)
    assert list(out["region_name"]) == ["a"], "fallback is lowest start"
    assert any(
        "ignores strand" in r.getMessage() for r in caplog.records
    ), "the fallback must say it is not strand-aware"


def test_flags_present_but_none_set_falls_back_and_warns(tmp_path, caplog):
    """Distinct from the column-absent case, and worth its own message.

    An all-zero flag column means the asset was built but matched nothing,
    which is a different problem from an asset that predates the columns.
    """
    df = _regions([_row("TP53", 100, "a"), _row("TP53", 900, "b")])
    with caplog.at_level("WARNING"):
        out = _run(tmp_path, df)
    assert list(out["region_name"]) == ["a"]
    assert any(
        "no region is marked" in r.getMessage() for r in caplog.records
    ), "an all-zero flag column needs its own warning"


def test_a_missing_input_returns_none_rather_than_raising(tmp_path):
    assert filter_fsc_to_e1(tmp_path / "absent.FSC.regions.tsv") is None
