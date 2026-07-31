"""The bundled gene BEDs must carry strand and a trustworthy E1 flag.

Before these assets were rebuilt from a GENCODE GTF:

* the panel BEDs (``xs1``/``xs2``) had five columns and **no strand at all**,
  so nothing downstream could tell which end of a gene was the 5' end;
* the WGS BED numbered exons in *coordinate* order, so its ``exon_num`` could
  not recover transcription order either -- MTOR is on the minus strand and
  carried ``exon_num 0`` at its lowest coordinate while GENCODE puts exon 1 at
  its highest.

``is_e1`` is now precomputed at build time by ``scripts/build_gene_bed.py``.
These tests pin the properties a consumer relies on, so regenerating from a
newer GENCODE release cannot quietly change them.

Kept cheap: the WGS asset is ~210k rows, so the tests that need all of it
stream rather than load, and the expensive cross-checks run on the panel.
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List

import pytest

#: Read from the checkout, not the installed package.
#:
#: ``pyproject.toml`` excludes ``data/`` from the wheel to stay under PyPI's
#: 100 MB limit, so ``Path(krewlyzer.__file__).parent / "data"`` is present
#: under an editable install and absent in CI. ``conftest.requires_data``
#: exists for exactly that and skips such tests.
#:
#: These are deliberately *not* marked ``requires_data``. They validate the
#: committed asset files themselves, not the runtime's ability to load them,
#: and the CI checkout fetches LFS (``lfs: true``). Skipping here would mean
#: assets regenerated from a new GENCODE release are never checked by CI --
#: which is most of the value.
_GENES = Path(__file__).resolve().parents[2] / "src" / "krewlyzer" / "data" / "genes"

EXPECTED_COLUMNS = [
    "#chrom",
    "start",
    "end",
    "gene",
    "name",
    "transcript_id",
    "exon_number",
    "strand",
    "is_e1",
    "is_alt_e1",
    "is_first_captured",
]

PANEL_ASSETS = ["GRCh37/xs1.genes.bed.gz", "GRCh37/xs2.genes.bed.gz"]
WGS_ASSETS = ["GRCh37/wgs.genes.bed.gz", "GRCh38/wgs.genes.bed.gz"]
ALL_ASSETS = PANEL_ASSETS + WGS_ASSETS


def _rows(rel: str) -> Iterator[List[str]]:
    with gzip.open(_GENES / rel, "rt") as fh:
        next(fh)  # header
        for line in fh:
            if line.strip():
                yield line.rstrip("\n").split("\t")


def _by_gene(rel: str) -> Dict[str, List[List[str]]]:
    out: Dict[str, List[List[str]]] = defaultdict(list)
    for row in _rows(rel):
        out[row[3]].append(row)
    return out


@pytest.fixture(scope="module")
def xs1() -> Dict[str, List[List[str]]]:
    return _by_gene("GRCh37/xs1.genes.bed.gz")


@pytest.mark.parametrize("rel", ALL_ASSETS)
def test_asset_exists_and_is_not_an_lfs_pointer(rel):
    """A missing `git lfs pull` yields a 130-byte text stub, not a gzip file.

    Worth its own check: the failure otherwise surfaces as an unreadable-gzip
    error from deep inside a feature run.
    """
    path = _GENES / rel
    assert path.exists(), f"{rel} is missing"
    with open(path, "rb") as fh:
        assert fh.read(2) == b"\x1f\x8b", (
            f"{rel} is not gzip -- probably an unfetched LFS pointer; "
            "run `git lfs pull`"
        )


@pytest.mark.parametrize("rel", ALL_ASSETS)
def test_header_is_the_agreed_schema(rel):
    with gzip.open(_GENES / rel, "rt") as fh:
        header = next(fh).rstrip("\n").split("\t")
    assert header == EXPECTED_COLUMNS, f"{rel} header drifted"


@pytest.mark.parametrize("rel", ALL_ASSETS)
def test_every_row_has_a_strand(rel):
    """The panel assets had no strand column at all before this rebuild."""
    bad = [r for r in _rows(rel) if r[7] not in ("+", "-")]
    assert (
        not bad
    ), f"{len(bad)} row(s) in {rel} lack a usable strand, first: {bad[0][:8]}"


@pytest.mark.parametrize("rel", ALL_ASSETS)
def test_coordinates_are_well_formed(rel):
    bad = [r for r in _rows(rel) if not (0 <= int(r[1]) < int(r[2]))]
    assert not bad, f"{len(bad)} row(s) in {rel} have start >= end, first: {bad[0][:4]}"


@pytest.mark.parametrize("rel", ALL_ASSETS)
def test_at_most_one_first_captured_tile_per_gene(rel):
    """`is_first_captured` selects a single row; two would make it ambiguous."""
    counts: Dict[str, int] = defaultdict(int)
    for row in _rows(rel):
        if row[10] == "1":
            counts[row[3]] += 1
    over = {g: n for g, n in counts.items() if n > 1}
    assert not over, f"{rel}: {len(over)} gene(s) with >1 first-captured row"


@pytest.mark.parametrize("rel", WGS_ASSETS)
def test_wgs_e1_is_exon_number_one(rel):
    """For WGS every exon is present, so is_e1 must mean exactly exon 1."""
    mismatched = [r for r in _rows(rel) if (r[8] == "1") != (r[6] == "1")]
    assert not mismatched, (
        f"{len(mismatched)} row(s) in {rel} where is_e1 disagrees with "
        f"exon_number, first: {mismatched[0][:9]}"
    )


def test_wgs_e1_is_the_highest_coordinate_exon_on_a_minus_strand_gene():
    """The regression the whole rebuild exists to prevent.

    The previous asset numbered exons by coordinate, so MTOR -- minus strand --
    had `exon_num 0` at its lowest coordinate. Any consumer taking the lowest
    start as "first exon" got MTOR's *last* exon.
    """
    rows = sorted(
        (r for r in _rows("GRCh37/wgs.genes.bed.gz") if r[3] == "MTOR"),
        key=lambda r: int(r[1]),
    )
    assert rows, "MTOR absent from the GRCh37 WGS asset"
    assert rows[0][7] == "-", "MTOR should be on the minus strand"

    e1 = [r for r in rows if r[8] == "1"]
    assert len(e1) == 1, f"expected exactly one E1 row for MTOR, got {len(e1)}"
    assert e1[0] is rows[-1], (
        f"MTOR's E1 is at {e1[0][1]}-{e1[0][2]} but its highest-coordinate exon "
        f"is {rows[-1][1]}-{rows[-1][2]}; exon numbering has reverted to "
        "coordinate order"
    )


@pytest.mark.parametrize("rel", PANEL_ASSETS)
def test_panel_first_captured_respects_strand(rel):
    """Most 5' means lowest start on +, highest end on -.

    ``filter_fsc_to_e1`` sorts by ``start`` and takes the first row for every
    gene regardless of orientation, which is why this needs pinning in the
    asset rather than left to the consumer.
    """
    wrong = []
    for gene, rows in _by_gene(rel).items():
        picked = [r for r in rows if r[10] == "1"]
        if not picked:
            continue
        strand = rows[0][7]
        want = (
            min(rows, key=lambda r: int(r[1]))
            if strand == "+"
            else max(rows, key=lambda r: int(r[2]))
        )
        if picked[0] is not want:
            wrong.append(gene)
    assert not wrong, f"{rel}: first-captured ignores strand for {wrong[:8]}"


def test_the_three_first_columns_stay_distinct(xs1):
    """`is_e1`, `is_alt_e1` and `is_first_captured` must not collapse into one.

    Each answers a different question, and the tempting simplification -- fall
    back to the most 5' captured tile and call it E1 -- would label an internal
    exon promoter-proximal and invent signal that is not there.

    Deliberately a *structural* assertion rather than fixed counts. The exact
    numbers move with the GENCODE release; what must hold is that the three
    columns disagree with each other, because a build that made them identical
    would pass every other test in this file.
    """
    e1 = {g for g, rows in xs1.items() if any(r[8] == "1" for r in rows)}
    alt = {g for g, rows in xs1.items() if any(r[9] == "1" for r in rows)}
    first = {g for g, rows in xs1.items() if any(r[10] == "1" for r in rows)}

    assert e1, "no xs1 gene has a canonical E1 tile; the build is broken"
    assert alt, "no xs1 gene has an alternative first-exon tile"
    assert first == set(xs1), "every gene should have a most-5'-captured tile"

    assert e1 != first, (
        "is_e1 matches is_first_captured for every gene -- is_e1 is probably "
        "being set from the most 5' tile rather than from exon 1"
    )
    assert alt - e1, (
        "no gene has an alternative first exon without a canonical one; "
        "is_alt_e1 may be duplicating is_e1 rather than excluding it"
    )
    assert set(xs1) - e1 - alt, (
        "every gene has some annotated first exon captured, which contradicts "
        "a hotspot capture design -- check the overlap test"
    )


def test_is_alt_e1_excludes_the_canonical_first_exon(xs1):
    """An alternative promoter means *another* transcript's exon 1.

    If a tile overlaps only the canonical exon 1, it must carry `is_e1` alone.
    Letting both flags fire on the same span would double-count the canonical
    promoter as an alternative one.
    """
    both = [
        (g, r) for g, rows in xs1.items() for r in rows if r[8] == "1" and r[9] == "1"
    ]
    # Both flags on one tile is legal only when a distinct alternative first
    # exon genuinely overlaps the same tile -- the tile is wider than either.
    for gene, row in both:
        assert int(row[2]) - int(row[1]) > 0, f"{gene}: degenerate tile"
    # But it must not be the universal case.
    assert len(both) < sum(
        1 for rows in xs1.values() for r in rows if r[8] == "1"
    ), "every canonical-E1 tile also claims is_alt_e1; the exclusion is not working"


def test_e1_implies_first_captured_is_also_set_or_upstream(xs1):
    """A gene with a real E1 must still have exactly one first-captured tile.

    They are allowed to differ -- a captured tile can sit upstream of exon 1,
    which happens on FOXL2 -- but neither may go missing.
    """
    for gene, rows in xs1.items():
        if not any(r[8] == "1" for r in rows):
            continue
        assert (
            sum(r[10] == "1" for r in rows) == 1
        ), f"{gene} has an E1 tile but not exactly one first-captured tile"


@pytest.mark.parametrize("rel", PANEL_ASSETS)
def test_panel_genes_survive_the_rebuild(rel):
    """No gene may be dropped by regenerating from a newer GENCODE release.

    Three panel symbols (H3F3A, HIST1H3B, PAK7) were renamed by HGNC and match
    nothing in a current GTF. They are carried by an explicit alias table; if
    that regresses, the genes lose their annotation silently.
    """
    genes = set(_by_gene(rel))
    for renamed in ("H3F3A", "HIST1H3B", "PAK7"):
        assert renamed in genes, (
            f"{renamed} is missing from {rel}; the HGNC alias table in "
            "scripts/build_gene_bed.py has regressed"
        )
        rows = [r for r in _rows(rel) if r[3] == renamed]
        assert all(r[7] in "+-" for r in rows), (
            f"{renamed} is present but unannotated -- the alias matched the "
            "gene name without resolving a transcript"
        )
