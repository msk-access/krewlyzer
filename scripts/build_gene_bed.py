#!/usr/bin/env python3
"""Build the bundled gene BED assets from a GENCODE GTF.

Why this exists
---------------
E1 -- the first exon of a gene -- is where krewlyzer looks for promoter-proximal
(nucleosome-depleted) fragmentation signal. Deriving it at runtime from the
coordinates alone gets it wrong for every gene on the minus strand, because
there the first exon is at the *highest* coordinate, not the lowest. Two
independent bugs followed from that:

* ``fsc_processor.filter_fsc_to_e1`` sorts by ``start`` and takes the first row
  per gene, with no strand input available to it at all.
* the bundled WGS BED numbers exons in *coordinate* order, so its ``exon_num``
  column cannot be used to recover the real order either. Verified: MTOR is on
  the minus strand and carries ``exon_num 0`` at 11,167,541, its lowest
  coordinate, while GENCODE places exon 1 at 11,322,502.

The panel BEDs (``xs1``/``xs2``) have no strand column whatsoever -- five
columns, and the rows are capture-probe tiles (``exon_MTOR_48a.1_2``) rather
than exons. Nothing in the existing assets can answer "which exon is first".

A GTF can. GENCODE's ``exon_number`` attribute is transcription-ordered, and
``tag "MANE_Select"`` identifies the clinically agreed transcript. This script
resolves the question once, at build time, and writes the answer into the asset
as an explicit ``is_e1`` column -- so the runtime reads a fact instead of
re-deriving it, and the bug class disappears rather than this instance of it.

Canonical transcript policy
---------------------------
See :data:`CANONICAL_POLICY`. An explicit ``--transcript-overrides`` entry wins,
then MANE Select, then GENCODE basic protein-coding, then Ensembl canonical,
then longest CDS. Decided here rather than per-run, and recorded in
``.agents/rules/architecture.md``.

Three notions of "first"
------------------------
Genes carry a median of 13 distinct annotated first exons (max 61 across these
panels), so "the first exon" is not one thing and a single boolean cannot say
what a tile is. Three columns, because they answer different questions:

``is_e1``
    Overlaps exon 1 of the resolved canonical transcript. Strict, reproducible,
    and what an override changes.
``is_alt_e1``
    Overlaps exon 1 of some *other* basic protein-coding transcript -- an
    alternative promoter. On xs1, 25 genes have a canonical E1 tile while 79
    have a tile on some annotated first exon; that gap is this column.
``is_first_captured``
    The most 5' captured tile in transcription order. Always exists, but an
    internal exon is not a promoter proxy and must not be read as E1.

A fallback chain cannot substitute for ``is_alt_e1``: every xs1 and xs2 gene has
a MANE transcript, so the lower tiers never fire and widening them changes
nothing.

Genome builds
-------------
All current MSK-ACCESS data is hg19/GRCh37, so that is the build that matters.
Ensembl's own GRCh37 line is frozen at release 87 (2016) and predates MANE
entirely; GENCODE's lift37 is a modern release mapped back, and does carry MANE
tags. Verified before use rather than assumed -- see ``_require_mane``.

Which GENCODE file
------------------
Use the **comprehensive** annotation (``gencode.vNN[lift37].annotation.gtf.gz``),
not ``.basic.``.

For the panel assays it makes no difference -- building ``xs1`` from both gives
the same 128 MANE transcripts, the same 25 genes with a true E1, and the same
1,725 rows, because MANE Select transcripts are in both files. It matters for
WGS: ``basic`` is a curated subset, and the ~5% of genes with no MANE tag fall
through to the Ensembl-canonical and longest-CDS branches, which need the full
transcript set to choose from sensibly.

The assets in this repo were built from ``gencode.v47lift37`` (GRCh37) and
``gencode.v50`` (GRCh38), comprehensive in both cases.

Usage
-----
    python scripts/build_gene_bed.py \\
        --gtf ~/Downloads/gencode.v47lift37.annotation.gtf.gz \\
        --genome GRCh37 --assay xs1 \\
        --probes src/krewlyzer/data/genes/GRCh37/xs1.genes.bed.gz \\
        --output src/krewlyzer/data/genes/GRCh37/xs1.genes.bed.gz

    python scripts/build_gene_bed.py \\
        --gtf ~/Downloads/gencode.v47lift37.annotation.gtf.gz \\
        --genome GRCh37 --assay wgs \\
        --output src/krewlyzer/data/genes/GRCh37/wgs.genes.bed.gz
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

#: Emitted as the first line of every generated asset.
#:
#: Three separate "first" columns -- see the module docstring. They are not
#: alternatives to each other and a tile can carry more than one: on xs1, 25
#: genes have a canonical-E1 tile, 15 more have only an alternative-promoter
#: tile, and 88 have neither.
#:
#: The E1 rationale (Helzer 2025) is promoter-proximal nucleosome depletion, so
#: labelling a merely-most-5' internal exon ``is_e1`` would state something
#: false. Each column is emitted separately and left for the model to weigh.
#:
#: The first five columns are unchanged from the previous GENE_BED format, so
#: readers indexing fields[3] (gene) and fields[4] (name) keep working against
#: a regenerated asset. Columns six onward are additive.
HEADER = (
    "#chrom\tstart\tend\tgene\tname\ttranscript_id\texon_number\tstrand"
    "\tis_e1\tis_alt_e1\tis_first_captured"
)

_ATTR = {
    key: re.compile(rf'{key} "([^"]+)"')
    for key in ("gene_name", "transcript_id", "gene_type", "transcript_type")
}
_EXON_NUMBER = re.compile(r"exon_number (\d+)")

#: HGNC renamed these after the panel BEDs were authored. Without the mapping
#: they simply fail to match and the gene silently loses its E1 annotation --
#: the exact failure mode this codebase keeps producing. Kept explicit and
#: asserted (see --allow-missing-genes) rather than resolved by fuzzy matching.
LEGACY_GENE_ALIASES: Dict[str, str] = {
    "H3F3A": "H3-3A",
    "HIST1H3B": "H3C2",
    "PAK7": "PAK5",
}


#: How one transcript per gene is chosen. Ordered; first match wins.
#:
#: The override tier exists because a capture panel is designed around
#: particular transcripts, and imposing MANE on it annotates a structure the
#: assay was not built for.
#:
#: Note what a *fallback* chain can and cannot do. Every gene in xs1 and xs2
#: has a MANE Select transcript, so tiers 3-5 never fire for those panels and
#: the chain resolves to MANE for every gene without an override. Widening the
#: fallback would not change that; capturing alternative promoters needs a
#: separate signal, which is what ``is_alt_e1`` is for.
CANONICAL_POLICY = (
    "1. --transcript-overrides entry for the gene",
    "2. MANE Select",
    "3. GENCODE basic AND protein_coding",
    "4. Ensembl canonical",
    "5. longest CDS",
)


class BuildError(RuntimeError):
    """A condition that must stop the build rather than degrade the asset."""


@dataclass
class Transcript:
    transcript_id: str
    gene: str
    chrom: str
    strand: str
    is_mane: bool = False
    is_canonical: bool = False
    is_basic: bool = False
    is_protein_coding: bool = False
    cds_length: int = 0
    #: (exon_number, start0, end) with exon_number straight from the GTF.
    exons: List[Tuple[int, int, int]] = field(default_factory=list)

    def sort_key(self) -> Tuple[int, int, int, int]:
        """Higher is better. Mirrors :data:`CANONICAL_POLICY` exactly.

        The override tier is not represented here -- it short-circuits in
        :func:`pick_canonical` rather than competing on score, so an override
        wins even against a MANE transcript.

        There is deliberately no separate ``is_protein_coding`` tier: a
        non-coding transcript accumulates no CDS records, so its
        ``cds_length`` is 0 and the final tier already ranks coding above
        non-coding. Adding one would be a second expression of the same rule.
        """
        return (
            int(self.is_mane),
            int(self.is_basic and self.is_protein_coding),
            int(self.is_canonical),
            self.cds_length,
        )

    def first_exon(self) -> Tuple[int, int, int]:
        """The exon GENCODE numbers 1 -- transcription order, not coordinate."""
        return min(self.exons, key=lambda e: e[0])


def normalise_chrom(chrom: str, prefix: bool) -> str:
    """Apply or strip a ``chr`` prefix.

    The runtime strips ``chr`` on both sides before matching
    (``fsc.rs`` calls ``trim_start_matches("chr")``), so this is consistency
    rather than correctness -- the panel assets are bare and the WGS assets are
    prefixed today, which makes them needlessly hard to diff against each other.
    """
    bare = chrom[3:] if chrom.startswith("chr") else chrom
    return f"chr{bare}" if prefix else bare


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def parse_gtf(
    gtf: Path,
    wanted: Optional[Set[str]] = None,
    protein_coding_only: bool = False,
) -> Dict[str, List[Transcript]]:
    """Collect exon and CDS records per transcript, keyed by gene symbol.

    ``wanted`` restricts to an explicit gene set -- the panel case, where the
    probe BED already decides the scope.

    ``protein_coding_only`` is what the WGS asset uses. Without it GENCODE
    yields 78,000 genes against the existing asset's 19,990, the difference
    being lncRNA and pseudogenes. That is a fourfold change in the size of
    every WGS gene table, so it is an explicit choice rather than a side
    effect of switching annotation source.
    """
    by_id: Dict[str, Transcript] = {}
    n_lines = 0

    with _open(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            n_lines += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature, attrs = fields[2], fields[8]
            if feature not in ("transcript", "exon", "CDS"):
                continue

            gene_match = _ATTR["gene_name"].search(attrs)
            if gene_match is None:
                continue
            gene = gene_match.group(1)
            if wanted is not None and gene not in wanted:
                continue
            if protein_coding_only:
                gtype = _ATTR["gene_type"].search(attrs)
                if gtype is None or gtype.group(1) != "protein_coding":
                    continue

            tx_match = _ATTR["transcript_id"].search(attrs)
            if tx_match is None:
                continue
            tx_id = tx_match.group(1)

            tx = by_id.get(tx_id)
            if tx is None:
                tx = by_id[tx_id] = Transcript(
                    transcript_id=tx_id,
                    gene=gene,
                    chrom=fields[0],
                    strand=fields[6],
                )

            if feature == "transcript":
                tx.is_mane = 'tag "MANE_Select"' in attrs
                tx.is_canonical = 'tag "Ensembl_canonical"' in attrs
                tx.is_basic = 'tag "basic"' in attrs
                ttype = _ATTR["transcript_type"].search(attrs)
                tx.is_protein_coding = (
                    ttype is not None and ttype.group(1) == "protein_coding"
                )
            elif feature == "CDS":
                tx.cds_length += int(fields[4]) - int(fields[3]) + 1
            else:  # exon
                num = _EXON_NUMBER.search(attrs)
                if num is None:
                    # Every GENCODE exon carries one. Its absence means the
                    # input is not the format this script was written for, and
                    # guessing an order is how E1 got wrong in the first place.
                    raise BuildError(
                        f"exon without exon_number in {gtf.name}: {attrs[:120]}"
                    )
                tx.exons.append((int(num.group(1)), int(fields[3]) - 1, int(fields[4])))

    logger.info(
        "parsed %s: %d records -> %d transcripts", gtf.name, n_lines, len(by_id)
    )

    by_gene: Dict[str, List[Transcript]] = defaultdict(list)
    for tx in by_id.values():
        if tx.exons:
            by_gene[tx.gene].append(tx)
    return dict(by_gene)


def _require_mane(by_gene: Dict[str, List[Transcript]], gtf_name: str) -> None:
    """Fail if the GTF carries no MANE tags at all.

    Ensembl's frozen GRCh37 release 87 has none. Building from it would fall
    straight through to the longest-CDS fallback for every gene and produce an
    asset that looks fine and encodes a different transcript set.
    """
    if any(tx.is_mane for txs in by_gene.values() for tx in txs):
        return
    raise BuildError(
        f"{gtf_name} contains no MANE_Select tags. The canonical-transcript "
        "policy needs them; a GRCh37 build must use GENCODE's lift37, not "
        "Ensembl's frozen release 87."
    )


def read_transcript_overrides(path: Path) -> Dict[str, str]:
    """Read a two-column ``gene<TAB>transcript_id`` TSV.

    Exists because a panel is designed around particular transcripts, and
    imposing MANE on it silently annotates a different gene structure than the
    assay was built for. An override says "for this gene, that one".

    Every malformed line is an error. A skipped row here would drop the gene
    back to MANE without saying so, which is the failure this file exists to
    prevent.
    """
    overrides: Dict[str, str] = {}
    with _open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2 or not all(parts):
                raise BuildError(
                    f"{path}:{lineno}: expected 'gene<TAB>transcript_id', got "
                    f"{line!r}"
                )
            gene, tx_id = parts
            if gene in overrides and overrides[gene] != tx_id:
                raise BuildError(
                    f"{path}:{lineno}: {gene} is listed twice with different "
                    f"transcripts ({overrides[gene]!r} and {tx_id!r})"
                )
            overrides[gene] = tx_id
    logger.info("read %d transcript override(s) from %s", len(overrides), path.name)
    return overrides


def _match_override(txs: List[Transcript], wanted: str) -> Optional[Transcript]:
    """Find the requested transcript, tolerating a missing version suffix.

    GENCODE ids carry a version (``ENST00000361445.9``) and lift37 appends a
    further suffix (``_9``). Requiring an exact match would make every override
    file build-specific, so a bare ``ENST00000361445`` matches too -- but only
    if it is unambiguous.
    """
    exact = [t for t in txs if t.transcript_id == wanted]
    if exact:
        return exact[0]
    base = wanted.split(".")[0]
    loose = [t for t in txs if t.transcript_id.split(".")[0] == base]
    if len(loose) == 1:
        return loose[0]
    if len(loose) > 1:
        raise BuildError(
            f"{wanted!r} matches {len(loose)} transcripts "
            f"({', '.join(sorted(t.transcript_id for t in loose))}); "
            "specify the full versioned id"
        )
    return None


def pick_canonical(
    by_gene: Dict[str, List[Transcript]],
    overrides: Optional[Dict[str, str]] = None,
) -> Dict[str, Transcript]:
    """Resolve one transcript per gene, following :data:`CANONICAL_POLICY`."""
    overrides = overrides or {}
    chosen: Dict[str, Transcript] = {}
    counts = {
        "override": 0,
        "mane": 0,
        "basic_protein_coding": 0,
        "canonical": 0,
        "longest_cds": 0,
    }

    for gene, txs in by_gene.items():
        wanted = overrides.get(gene)
        if wanted is not None:
            match = _match_override(txs, wanted)
            if match is None:
                # Hard error, never a fallback: an override that silently
                # reverts to MANE gives an asset that disagrees with the file
                # the operator wrote, with nothing to indicate it.
                raise BuildError(
                    f"override transcript {wanted!r} for {gene} is not in the "
                    f"GTF under that gene. Present for {gene}: "
                    f"{', '.join(sorted(t.transcript_id for t in txs)[:6])}"
                    f"{' ...' if len(txs) > 6 else ''}"
                )
            chosen[gene] = match
            counts["override"] += 1
            continue

        best = max(txs, key=Transcript.sort_key)
        chosen[gene] = best
        if best.is_mane:
            counts["mane"] += 1
        elif best.is_basic and best.is_protein_coding:
            counts["basic_protein_coding"] += 1
        elif best.is_canonical:
            counts["canonical"] += 1
        else:
            counts["longest_cds"] += 1

    logger.info(
        "canonical transcripts: %d override, %d MANE, %d basic protein-coding, "
        "%d Ensembl-canonical, %d longest-CDS",
        counts["override"],
        counts["mane"],
        counts["basic_protein_coding"],
        counts["canonical"],
        counts["longest_cds"],
    )
    unused = sorted(set(overrides) - set(by_gene))
    if unused:
        # Not fatal -- one override file may serve several assays -- but silence
        # would let a typo look like a working override.
        logger.warning(
            "%d override(s) name a gene absent from this build and had no "
            "effect: %s%s",
            len(unused),
            ", ".join(unused[:10]),
            ", ..." if len(unused) > 10 else "",
        )
    return chosen


def resolve_aliases(
    wanted: Set[str], available: Set[str]
) -> Tuple[Set[str], Dict[str, str]]:
    """Map panel symbols the GTF no longer uses onto its current ones.

    Returns the lookup set and a GTF-symbol -> panel-symbol map, so output
    keeps the name the panel and every downstream table already use.
    """
    lookup, back = set(), {}
    for gene in wanted:
        if gene in available:
            lookup.add(gene)
            back[gene] = gene
            continue
        alias = LEGACY_GENE_ALIASES.get(gene)
        if alias and alias in available:
            lookup.add(alias)
            back[alias] = gene
            logger.info("alias: panel %s -> GTF %s", gene, alias)
    return lookup, back


def read_probe_bed(path: Path) -> List[Tuple[str, int, int, str, str]]:
    """Read a panel probe BED: chrom, start, end, gene, name."""
    rows = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            name = f[4] if len(f) > 4 else f[3]
            rows.append((f[0], int(f[1]), int(f[2]), f[3], name))
    logger.info("read %d probe tiles from %s", len(rows), path.name)
    return rows


def alternative_first_exons(
    by_gene: Dict[str, List[Transcript]],
    canonical: Dict[str, Transcript],
) -> Dict[str, List[Tuple[int, int]]]:
    """First exons of basic protein-coding transcripts *other* than the canonical.

    Genes carry a median of 13 distinct annotated first exons (max 61 across
    these panels), because alternative promoters are the norm rather than the
    exception. A tile can therefore sit on a genuine transcription start while
    missing the canonical transcript's exon 1 entirely -- 79 of 128 xs1 genes
    have a tile on *some* first exon against 25 on the MANE one.

    Restricted to basic protein-coding transcripts: the full set includes
    lowly-supported and non-coding isoforms whose annotated start is not
    evidence of an active promoter, and counting those would inflate the signal
    rather than measure it.
    """
    out: Dict[str, List[Tuple[int, int]]] = {}
    for gene, txs in by_gene.items():
        chosen = canonical.get(gene)
        spans = {
            (s, e)
            for tx in txs
            if tx.is_basic
            and tx.is_protein_coding
            and (chosen is None or tx.transcript_id != chosen.transcript_id)
            for n, s, e in tx.exons
            if n == 1
        }
        if chosen is not None:
            # Same span as the canonical exon 1 is is_e1, not an alternative.
            _, cs, ce = chosen.first_exon()
            spans.discard((cs, ce))
        if spans:
            out[gene] = sorted(spans)
    return out


def build_panel_rows(
    probes: Sequence[Tuple[str, int, int, str, str]],
    canonical: Dict[str, Transcript],
    alias_back: Dict[str, str],
    prefix: bool,
    alt_e1: Optional[Dict[str, List[Tuple[int, int]]]] = None,
) -> Tuple[List[List[str]], Set[str]]:
    """Keep the probe tiles; annotate strand and the three notions of "first".

    The tiles are the unit FSC and MDS count over, so they stay the rows. What
    the GTF adds is the strand (absent from the panel BEDs entirely), the
    transcript they belong to, and:

    * ``is_e1`` -- overlaps the canonical transcript's exon 1
    * ``is_alt_e1`` -- overlaps another basic protein-coding transcript's exon 1
    * ``is_first_captured`` -- most 5' tile for the gene, in transcription order

    Three columns rather than one because they answer different questions and
    conflating them is what made the previous E1 table unreadable.
    """
    by_panel_gene = {
        alias_back[g]: tx for g, tx in canonical.items() if g in alias_back
    }
    alt_by_panel_gene = {
        alias_back[g]: spans for g, spans in (alt_e1 or {}).items() if g in alias_back
    }

    per_gene: Dict[str, List[int]] = defaultdict(list)
    staged: List[List] = []

    for chrom, start, end, gene, name in probes:
        tx = by_panel_gene.get(gene)
        if tx is None:
            # No canonical transcript. Emit the tile -- dropping it would
            # change what FSC counts -- but leave every derived field blank
            # rather than guessing.
            staged.append(
                [
                    normalise_chrom(chrom, prefix),
                    start,
                    end,
                    gene,
                    name,
                    ".",
                    ".",
                    ".",
                    "0",
                    "0",
                    "0",
                ]
            )
            continue
        e1_num, e1_start, e1_end = tx.first_exon()
        overlaps = start < e1_end and e1_start < end
        alt = any(
            start < ae and as_ < end for as_, ae in alt_by_panel_gene.get(gene, ())
        )
        staged.append(
            [
                normalise_chrom(chrom, prefix),
                start,
                end,
                gene,
                name,
                tx.transcript_id,
                e1_num if overlaps else ".",
                tx.strand,
                "1" if overlaps else "0",
                "1" if alt else "0",
                "0",  # is_first_captured, filled in below
            ]
        )
        per_gene[gene].append(len(staged) - 1)

    # Most 5' captured tile, in transcription order: lowest start on the plus
    # strand, highest end on the minus strand. This is the strand awareness
    # `filter_fsc_to_e1` never had -- it sorted by start and took the first row
    # for every gene regardless of orientation.
    e1_genes = set()
    for gene, idxs in per_gene.items():
        strand = staged[idxs[0]][7]
        pick = (
            min(idxs, key=lambda i: staged[i][1])
            if strand == "+"
            else max(idxs, key=lambda i: staged[i][2])
        )
        staged[pick][10] = "1"
        if any(staged[i][8] == "1" for i in idxs):
            e1_genes.add(gene)

    return staged, e1_genes


def build_wgs_rows(
    canonical: Dict[str, Transcript],
    prefix: bool,
    alt_e1: Optional[Dict[str, List[Tuple[int, int]]]] = None,
) -> Tuple[List[List[str]], Set[str]]:
    """Emit every exon of each canonical transcript, transcription-numbered.

    Every exon of the transcript is present, so exon 1 always exists and
    ``is_first_captured`` is identical to ``is_e1``. Both are still written, so
    the panel and WGS assets share one schema and a consumer never has to
    branch on which it loaded.

    ``is_alt_e1`` is meaningful here too: an exon of the canonical transcript
    can coincide with the first exon of a different basic protein-coding
    transcript, which is the alternative-promoter case.
    """
    alt_e1 = alt_e1 or {}
    rows, e1_genes = [], set()
    for gene, tx in canonical.items():
        spans = alt_e1.get(gene, ())
        for num, start, end in sorted(tx.exons, key=lambda e: e[1]):
            is_e1 = num == 1
            if is_e1:
                e1_genes.add(gene)
            flag = "1" if is_e1 else "0"
            alt = any(start < ae and as_ < end for as_, ae in spans)
            rows.append(
                [
                    normalise_chrom(tx.chrom, prefix),
                    start,
                    end,
                    gene,
                    f"exon_{gene}_{num}",
                    tx.transcript_id,
                    num,
                    tx.strand,
                    flag,
                    "1" if alt else "0",
                    flag,
                ]
            )
    return rows, e1_genes


def _chrom_sort_key(chrom: str) -> Tuple[int, str]:
    bare = chrom[3:] if chrom.startswith("chr") else chrom
    try:
        return (int(bare), "")
    except ValueError:
        return (100, bare)


def write_bed(rows: List[List[str]], output: Path) -> None:
    rows = sorted(
        rows, key=lambda r: (_chrom_sort_key(str(r[0])), int(r[1]), int(r[2]))
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if output.suffix == ".gz" else open
    with opener(output, "wt") as fh:
        fh.write(HEADER + "\n")
        for row in rows:
            fh.write("\t".join(str(v) for v in row) + "\n")
    logger.info("wrote %d rows to %s", len(rows), output)


def main(argv: Optional[Iterable[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--gtf", type=Path, required=True, help="GENCODE GTF (.gtf/.gtf.gz)"
    )
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--genome", choices=("GRCh37", "GRCh38"), required=True)
    ap.add_argument(
        "--assay",
        required=True,
        help="xs1/xs2 (panel, needs --probes) or wgs (all canonical exons)",
    )
    ap.add_argument(
        "--probes",
        type=Path,
        help="Existing panel probe BED whose tiles are kept as the rows",
    )
    ap.add_argument(
        "--chr-prefix",
        action="store_true",
        help="Emit chr1 rather than 1 (the runtime strips it either way)",
    )
    ap.add_argument(
        "--all-gene-types",
        action="store_true",
        help="WGS only: keep lncRNA and pseudogenes too (~78k genes rather "
        "than ~19k). Quadruples every WGS gene table; off by default so the "
        "annotation-source change does not silently redefine the asset's scope",
    )
    ap.add_argument(
        "--transcript-overrides",
        type=Path,
        help="Two-column TSV (gene<TAB>transcript_id) naming the transcript to "
        "use for specific genes, taking precedence over MANE. A transcript "
        "that is not in the GTF under that gene is a hard error",
    )
    ap.add_argument(
        "--allow-missing-genes",
        action="store_true",
        help="Downgrade unmatched panel genes from an error to a warning",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    if not args.gtf.exists():
        logger.error("GTF not found: %s", args.gtf)
        return 2
    if args.transcript_overrides and not args.transcript_overrides.exists():
        logger.error("override file not found: %s", args.transcript_overrides)
        return 2

    try:
        overrides = (
            read_transcript_overrides(args.transcript_overrides)
            if args.transcript_overrides
            else {}
        )

        if args.assay == "wgs":
            by_gene = parse_gtf(args.gtf, protein_coding_only=not args.all_gene_types)
            _require_mane(by_gene, args.gtf.name)
            canonical = pick_canonical(by_gene, overrides)
            alt_e1 = alternative_first_exons(by_gene, canonical)
            rows, e1_genes = build_wgs_rows(canonical, args.chr_prefix, alt_e1)
            expected = set(canonical)
        else:
            if not args.probes or not args.probes.exists():
                logger.error("--probes is required for a panel assay")
                return 2
            probes = read_probe_bed(args.probes)
            wanted = {p[3] for p in probes}
            # Two passes: the alias map needs to know which symbols the GTF has.
            everything = parse_gtf(args.gtf)
            _require_mane(everything, args.gtf.name)
            lookup, alias_back = resolve_aliases(wanted, set(everything))
            # Overrides are written against panel symbols; translate to the
            # GTF's current ones so a renamed gene can still be overridden.
            gtf_overrides = {
                LEGACY_GENE_ALIASES.get(g, g): tx for g, tx in overrides.items()
            }
            subset = {g: everything[g] for g in lookup}
            canonical = pick_canonical(subset, gtf_overrides)
            alt_e1 = alternative_first_exons(subset, canonical)

            unmatched = sorted(wanted - set(alias_back.values()))
            if unmatched:
                msg = (
                    f"{len(unmatched)} panel gene(s) absent from {args.gtf.name}: "
                    f"{', '.join(unmatched)}. Add them to LEGACY_GENE_ALIASES if "
                    "they were renamed; otherwise they get no E1 annotation."
                )
                if args.allow_missing_genes:
                    logger.warning(msg)
                else:
                    raise BuildError(msg)

            rows, e1_genes = build_panel_rows(
                probes, canonical, alias_back, args.chr_prefix, alt_e1
            )
            expected = wanted
    except BuildError as exc:
        logger.error("%s", exc)
        return 1

    alt_genes = {r[3] for r in rows if r[9] == "1"}
    without_any = sorted(expected - e1_genes - alt_genes)
    only_alt = sorted(alt_genes - e1_genes)

    if only_alt:
        # Not a warning: an alternative promoter is a real transcription start,
        # and on a hotspot panel it is often the only one captured. Reported so
        # the split between the two columns is visible at build time rather
        # than discovered during modelling.
        logger.info(
            "%d gene(s) have no tile on the canonical exon 1 but do have one on "
            "another basic protein-coding transcript's first exon; these carry "
            "is_alt_e1 without is_e1",
            len(only_alt),
        )
    if without_any:
        # These genuinely have no captured transcription start. A gene here is
        # absent from every E1-only table downstream rather than erroring, so
        # say it loudly.
        logger.warning(
            "%d of %d gene(s) have no tile on any annotated first exon -- they "
            "will be absent from E1-only outputs (%s%s)",
            len(without_any),
            len(expected),
            ", ".join(without_any[:10]),
            ", ..." if len(without_any) > 10 else "",
        )
        logger.warning(
            "  is_first_captured still marks their most 5' captured tile, but "
            "an internal exon is not a promoter proxy -- do not read it as E1."
        )

    # Partitioned deliberately: canonical-E1 and alternative-promoter genes
    # overlap (a tile can be both), so reporting the raw counts side by side
    # gives three numbers that do not add up to the total.
    logger.info(
        "%s/%s: %d rows, %d genes = %d with a canonical E1 tile + %d with only "
        "an alternative first exon + %d with neither (%d carry both flags)",
        args.genome,
        args.assay,
        len(rows),
        len(expected),
        len(e1_genes),
        len(only_alt),
        len(without_any),
        len(alt_genes & e1_genes),
    )
    write_bed(rows, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
