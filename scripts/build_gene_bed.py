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
MANE Select, then Ensembl canonical, then longest CDS. Decided here rather than
per-run, and recorded in ``.agents/rules/architecture.md``.

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
#: ``is_e1`` and ``is_first_captured`` are deliberately separate columns.
#: On a targeted panel they are usually *not* the same tile: MSK-ACCESS tiles
#: coding hotspot exons, and for 103 of xs1's 128 genes the canonical
#: transcript's exon 1 is not captured at all. AKT1's is 15 kb past the
#: panel's most 5' tile.
#:
#: The E1 rationale (Helzer 2025) is promoter-proximal nucleosome depletion.
#: An internal exon that merely happens to be the most 5' one captured is not
#: a promoter proxy, so labelling it ``is_e1`` would state something false.
#: Both are emitted and left for the model to weigh.
#: The first five columns are unchanged from the previous GENE_BED format,
#: so readers indexing fields[3] (gene) and fields[4] (name) keep working
#: against a regenerated asset. Columns six onward are additive.
HEADER = (
    "#chrom\tstart\tend\tgene\tname\ttranscript_id\texon_number\tstrand"
    "\tis_e1\tis_first_captured"
)

_ATTR = {
    key: re.compile(rf'{key} "([^"]+)"')
    for key in ("gene_name", "transcript_id", "gene_type")
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
    cds_length: int = 0
    #: (exon_number, start0, end) with exon_number straight from the GTF.
    exons: List[Tuple[int, int, int]] = field(default_factory=list)

    def sort_key(self) -> Tuple[int, int, int]:
        """Higher is better. Mirrors the documented policy exactly."""
        return (int(self.is_mane), int(self.is_canonical), self.cds_length)

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


def pick_canonical(by_gene: Dict[str, List[Transcript]]) -> Dict[str, Transcript]:
    chosen = {}
    counts = {"mane": 0, "canonical": 0, "longest_cds": 0}
    for gene, txs in by_gene.items():
        best = max(txs, key=Transcript.sort_key)
        chosen[gene] = best
        counts[
            (
                "mane"
                if best.is_mane
                else "canonical" if best.is_canonical else "longest_cds"
            )
        ] += 1
    logger.info(
        "canonical transcripts: %d MANE, %d Ensembl-canonical, %d longest-CDS",
        counts["mane"],
        counts["canonical"],
        counts["longest_cds"],
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


def build_panel_rows(
    probes: Sequence[Tuple[str, int, int, str, str]],
    canonical: Dict[str, Transcript],
    alias_back: Dict[str, str],
    prefix: bool,
) -> Tuple[List[List[str]], Set[str]]:
    """Keep the probe tiles; annotate strand, true E1, and most-5' captured.

    The tiles are the unit FSC and MDS count over, so they stay the rows. What
    the GTF adds is the strand (absent from the panel BEDs entirely), the
    transcript they belong to, and the two distinct notions of "first".
    """
    by_panel_gene = {
        alias_back[g]: tx for g, tx in canonical.items() if g in alias_back
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
                ]
            )
            continue
        e1_num, e1_start, e1_end = tx.first_exon()
        overlaps = start < e1_end and e1_start < end
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
        staged[pick][9] = "1"
        if any(staged[i][8] == "1" for i in idxs):
            e1_genes.add(gene)

    return staged, e1_genes


def build_wgs_rows(
    canonical: Dict[str, Transcript], prefix: bool
) -> Tuple[List[List[str]], Set[str]]:
    """Emit every exon of each canonical transcript, transcription-numbered.

    For WGS every exon of the transcript is present, so exon 1 is always
    captured and ``is_first_captured`` is identical to ``is_e1``. It is still
    written, so the two assay families share one schema and a consumer does not
    have to branch on which asset it loaded.
    """
    rows, e1_genes = [], set()
    for gene, tx in canonical.items():
        for num, start, end in sorted(tx.exons, key=lambda e: e[1]):
            is_e1 = num == 1
            if is_e1:
                e1_genes.add(gene)
            flag = "1" if is_e1 else "0"
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
        "--allow-missing-genes",
        action="store_true",
        help="Downgrade unmatched panel genes from an error to a warning",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    if not args.gtf.exists():
        logger.error("GTF not found: %s", args.gtf)
        return 2

    try:
        if args.assay == "wgs":
            by_gene = parse_gtf(args.gtf, protein_coding_only=not args.all_gene_types)
            _require_mane(by_gene, args.gtf.name)
            canonical = pick_canonical(by_gene)
            rows, e1_genes = build_wgs_rows(canonical, args.chr_prefix)
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
            canonical = pick_canonical({g: everything[g] for g in lookup})

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
                probes, canonical, alias_back, args.chr_prefix
            )
            expected = wanted
    except BuildError as exc:
        logger.error("%s", exc)
        return 1

    without_e1 = sorted(expected - e1_genes)
    if without_e1:
        # Loud, because a gene with no E1 row simply vanishes from every
        # E1-only table downstream instead of erroring there. On a hotspot
        # panel this is the normal case, not a build failure: the canonical
        # exon 1 is usually 5'UTR and outside the capture design.
        logger.warning(
            "%d of %d gene(s) have no tile overlapping the canonical exon 1 "
            "-- they will be absent from E1-only outputs (%s%s)",
            len(without_e1),
            len(expected),
            ", ".join(without_e1[:10]),
            ", ..." if len(without_e1) > 10 else "",
        )
        logger.warning(
            "  is_first_captured still marks their most 5' captured tile, but "
            "an internal exon is not a promoter proxy -- do not read it as E1."
        )

    logger.info(
        "%s/%s: %d rows, %d genes, %d with a true E1 tile",
        args.genome,
        args.assay,
        len(rows),
        len(expected),
        len(e1_genes),
    )
    write_bed(rows, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
