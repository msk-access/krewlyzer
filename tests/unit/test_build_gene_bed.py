"""Transcript selection and override handling in ``scripts/build_gene_bed.py``.

The asset tests check the *output*; these check the decisions that produced it,
including the failure paths a successful build never exercises.

An override that silently reverts to MANE would produce an asset disagreeing
with the file the operator wrote, with nothing to indicate it — the failure mode
this repo keeps finding. Every rejection below is asserted rather than assumed.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import List

import pytest

_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "build_gene_bed.py"


def _load():
    """Import the script by path; ``scripts/`` is not a package."""
    spec = importlib.util.spec_from_file_location("build_gene_bed", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["build_gene_bed"] = module
    spec.loader.exec_module(module)
    return module


bgb = _load()


def _tx(tx_id: str, gene: str = "TP53", **flags) -> "bgb.Transcript":
    t = bgb.Transcript(
        transcript_id=tx_id,
        gene=gene,
        chrom="chr17",
        strand="-",
        is_mane=flags.get("mane", False),
        is_canonical=flags.get("canonical", False),
        is_basic=flags.get("basic", False),
        is_protein_coding=flags.get("coding", False),
        cds_length=flags.get("cds", 0),
    )
    t.exons = flags.get("exons", [(1, 100, 200), (2, 300, 400)])
    return t


# ---------------------------------------------------------------------------
# Override file parsing
# ---------------------------------------------------------------------------


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "ov.tsv"
    p.write_text(text)
    return p


def test_override_file_reads_pairs_and_skips_comments(tmp_path):
    path = _write(tmp_path, "# header\nTP53\tENST1\n\nMTOR\tENST2\n")
    assert bgb.read_transcript_overrides(path) == {"TP53": "ENST1", "MTOR": "ENST2"}


@pytest.mark.parametrize(
    "text,fragment",
    [
        ("TP53 ENST1\n", "expected 'gene<TAB>transcript_id'"),
        ("TP53\n", "expected 'gene<TAB>transcript_id'"),
        ("TP53\tENST1\textra\n", "expected 'gene<TAB>transcript_id'"),
        ("TP53\t\n", "expected 'gene<TAB>transcript_id'"),
        ("TP53\tENST1\nTP53\tENST2\n", "listed twice"),
    ],
)
def test_malformed_override_lines_are_rejected(tmp_path, text, fragment):
    """Every bad line is fatal.

    Skipping one would drop that gene back to MANE without saying so, which is
    exactly what the override file exists to prevent.
    """
    with pytest.raises(bgb.BuildError, match=fragment):
        bgb.read_transcript_overrides(_write(tmp_path, text))


def test_a_gene_repeated_with_the_same_transcript_is_allowed(tmp_path):
    """Harmless duplication should not fail a build; only a conflict should."""
    path = _write(tmp_path, "TP53\tENST1\nTP53\tENST1\n")
    assert bgb.read_transcript_overrides(path) == {"TP53": "ENST1"}


# ---------------------------------------------------------------------------
# Canonical transcript policy
# ---------------------------------------------------------------------------


def test_mane_wins_over_everything_else():
    txs = [
        _tx("ENST_long", cds=9999, basic=True, coding=True),
        _tx("ENST_canon", canonical=True, basic=True, coding=True),
        _tx("ENST_mane", mane=True),
    ]
    assert bgb.pick_canonical({"TP53": txs})["TP53"].transcript_id == "ENST_mane"


def test_basic_protein_coding_outranks_ensembl_canonical():
    """The tier the user's policy inserts above Ensembl canonical."""
    txs = [
        _tx("ENST_canon_noncoding", canonical=True),
        _tx("ENST_basic_coding", basic=True, coding=True),
    ]
    assert (
        bgb.pick_canonical({"TP53": txs})["TP53"].transcript_id == "ENST_basic_coding"
    )


def test_longest_cds_is_the_last_resort():
    txs = [_tx("ENST_a", cds=100), _tx("ENST_b", cds=500)]
    assert bgb.pick_canonical({"TP53": txs})["TP53"].transcript_id == "ENST_b"


def test_an_override_beats_mane():
    """The whole point: an explicit choice outranks the default."""
    txs = [_tx("ENST_mane", mane=True), _tx("ENST_other", basic=True, coding=True)]
    chosen = bgb.pick_canonical({"TP53": txs}, {"TP53": "ENST_other"})
    assert chosen["TP53"].transcript_id == "ENST_other"


def test_an_unknown_override_transcript_is_fatal():
    """Never a silent fallback to MANE."""
    txs = [_tx("ENST_mane", mane=True)]
    with pytest.raises(bgb.BuildError, match="not in the GTF under that gene"):
        bgb.pick_canonical({"TP53": txs}, {"TP53": "ENST_nonexistent"})


def test_an_override_naming_another_genes_transcript_is_fatal():
    """Scoped per gene, so a copy-paste slip cannot silently mis-annotate."""
    genes = {
        "TP53": [_tx("ENST_tp53", mane=True)],
        "MTOR": [_tx("ENST_mtor", gene="MTOR", mane=True)],
    }
    with pytest.raises(bgb.BuildError, match="not in the GTF under that gene"):
        bgb.pick_canonical(genes, {"TP53": "ENST_mtor"})


def test_an_unversioned_override_matches_a_versioned_transcript():
    """Override files must not have to track GENCODE version suffixes.

    lift37 ids carry two (``ENST00000361445.9_9``); requiring an exact match
    would make every override file build-specific.
    """
    txs = [_tx("ENST00000269305.9_9", mane=True), _tx("ENST00000420246.6_5")]
    chosen = bgb.pick_canonical({"TP53": txs}, {"TP53": "ENST00000420246"})
    assert chosen["TP53"].transcript_id == "ENST00000420246.6_5"


def test_an_ambiguous_unversioned_override_is_fatal():
    """Two transcripts sharing a base id must not be resolved by guessing."""
    txs = [_tx("ENST00000269305.9_9"), _tx("ENST00000269305.10_1")]
    with pytest.raises(bgb.BuildError, match="specify the full versioned id"):
        bgb.pick_canonical({"TP53": txs}, {"TP53": "ENST00000269305"})


def test_an_override_for_an_absent_gene_warns_but_does_not_fail(caplog):
    """One override file may serve several assays, so this is not an error.

    It must still be said: otherwise a typo looks exactly like a working
    override.
    """
    txs = [_tx("ENST_mane", mane=True)]
    with caplog.at_level("WARNING"):
        bgb.pick_canonical({"TP53": txs}, {"NOTAGENE": "ENST_x"})
    assert any("had no effect" in r.getMessage() for r in caplog.records)


# ---------------------------------------------------------------------------
# Alternative first exons
# ---------------------------------------------------------------------------


def test_alternative_first_exons_exclude_the_canonical_one():
    canonical = _tx("ENST_mane", mane=True, exons=[(1, 100, 200), (2, 300, 400)])
    other = _tx("ENST_alt", basic=True, coding=True, exons=[(1, 50, 90)])
    same_start = _tx(
        "ENST_same", basic=True, coding=True, exons=[(1, 100, 200), (2, 500, 600)]
    )
    alt = bgb.alternative_first_exons(
        {"TP53": [canonical, other, same_start]}, {"TP53": canonical}
    )
    assert alt["TP53"] == [(50, 90)], (
        "a transcript whose exon 1 coincides with the canonical one is not an "
        "alternative promoter"
    )


def test_alternative_first_exons_ignore_non_basic_and_non_coding():
    """Minor isoforms' annotated starts are not evidence of a live promoter."""
    canonical = _tx("ENST_mane", mane=True, exons=[(1, 100, 200)])
    noncoding = _tx("ENST_nc", basic=True, coding=False, exons=[(1, 10, 20)])
    nonbasic = _tx("ENST_nb", basic=False, coding=True, exons=[(1, 30, 40)])
    good = _tx("ENST_ok", basic=True, coding=True, exons=[(1, 50, 60)])
    alt = bgb.alternative_first_exons(
        {"TP53": [canonical, noncoding, nonbasic, good]}, {"TP53": canonical}
    )
    assert alt["TP53"] == [(50, 60)]


def test_a_gene_with_only_a_canonical_transcript_has_no_alternatives():
    canonical = _tx("ENST_mane", mane=True, exons=[(1, 100, 200)])
    assert bgb.alternative_first_exons({"TP53": [canonical]}, {"TP53": canonical}) == {}


# ---------------------------------------------------------------------------
# Policy is documented where it is implemented
# ---------------------------------------------------------------------------


def test_the_policy_constant_matches_the_sort_key_arity():
    """`CANONICAL_POLICY` is prose; `sort_key` is what runs.

    They drift silently unless something ties them together. The override tier
    short-circuits rather than scoring, so the key carries the remaining four.
    """
    key: tuple = _tx("ENST", mane=True).sort_key()
    assert len(bgb.CANONICAL_POLICY) == len(key) + 1, (
        "CANONICAL_POLICY lists a different number of tiers than sort_key "
        "implements (plus the override tier, which short-circuits)"
    )


def test_header_lists_every_column_the_builders_emit():
    """A row shorter or longer than the header writes a malformed BED."""
    canonical = _tx("ENST_mane", mane=True, exons=[(1, 100, 200)])
    rows, _ = bgb.build_wgs_rows({"TP53": canonical}, prefix=True)
    n_header = len(bgb.HEADER.split("\t"))
    assert rows, "no rows produced"
    for row in rows:
        assert (
            len(row) == n_header
        ), f"row has {len(row)} fields, header declares {n_header}"


def test_panel_and_wgs_rows_have_the_same_width():
    """One schema for both assay families, so consumers never branch."""
    canonical = _tx("ENST_mane", mane=True, exons=[(1, 100, 200)])
    wgs, _ = bgb.build_wgs_rows({"TP53": canonical}, prefix=False)
    panel: List[List] = bgb.build_panel_rows(
        [("chr17", 100, 200, "TP53", "tile_1")],
        {"TP53": canonical},
        {"TP53": "TP53"},
        prefix=False,
    )[0]
    assert len(wgs[0]) == len(panel[0]) == len(bgb.HEADER.split("\t"))
