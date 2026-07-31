# E1 is `is_e1` OR `is_alt_e1`

**Decided 2026-07-31.** Applies to `FSC.regions.e1only` and to `region-mds`.

## The question

"First exon" is not one thing. Genes carry a **median of 13 distinct annotated
first exons** (max 61 across these panels), because alternative promoters are
the norm. Which of them counts as E1 decides what a promoter-proximal feature
actually measures.

## The numbers that drove it

Measured against `gencode.v47lift37` and the bundled panel BEDs:

| a tile overlaps… | xs1 / 128 | xs2 / 146 |
|---|---:|---:|
| the canonical (MANE) transcript's exon 1 | 25 | 33 |
| **either that or another basic protein-coding transcript's exon 1** | **40** | **48** |
| any transcript's exon 1, including minor isoforms | 79 | 90 |

MSK-ACCESS tiles coding hotspot exons, so the canonical exon 1 is frequently
5'UTR and outside the capture design — AKT1's sits 15 kb past the panel's most
5' tile.

## Why not the alternatives

**MANE only (25/33).** Discards genes the panel captures at a real alternative
promoter, which is most of the gap. Reproducible, but throws away a third of
the available signal for no biological reason.

**Any transcript (79/90).** Includes lowly-supported and non-coding isoforms
whose annotated start is not evidence of a live promoter. That inflates the
count rather than measuring anything.

**Most 5' captured tile (128/146).** Always available, and wrong: for the genes
with no captured first exon it is an internal exon, and labelling it E1 asserts
promoter proximity that is not there. This was the pre-0.9.0 behaviour.

## Consequences

Genes with neither flag are **omitted** from E1 outputs, not back-filled. On a
real xs2 sample `FSC.regions.e1only` went from 146 rows to 48, and every
remaining row is on an annotated first exon (24 were, before).

## What would reopen this

Evidence that alternative-promoter regions behave differently from canonical
ones in the fragmentation signal. Both flags are emitted separately precisely
so that can be tested rather than assumed — if `is_alt_e1` regions turn out to
carry no signal, narrowing to `is_e1` is a one-line change.

The canonical transcript itself is configurable per gene via
`--transcript-overrides` (see `scripts/build_gene_bed.py`), because a panel
designed around specific clinical transcripts should not have MANE imposed on
it.
