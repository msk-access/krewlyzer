---
description: What downstream consumers rely on, and how to change it safely
alwaysApply: false
---

# Output contract

The immediate consumer is `kreview`. It reads **Parquet only** — its own rules
forbid `features.json` — and resolves each sample as
`{results_dir}/{sample_id}/{sample_id}{suffix}`.

Two properties make this contract unusually unforgiving:

**`.metadata.parquet` is a completion marker.** A sample lacking it is dropped
from the cohort silently. Not warned about, not errored on — the sample simply
never enters the analysis.

**Every reader swallows exceptions** and yields an empty feature dict. A
malformed table degrades to missing columns rather than a crash. So a defect
travels all the way to a model fit without anyone seeing it, and nothing
downstream can tell us we shipped something wrong.

That is why the contract is asserted here: `src/krewlyzer/validate/contract.py`
declares it, and `krewlyzer validate-output` checks it.

## Changing the contract

- **Adding a table or renaming a suffix** needs an entry in `contract.py` *and*
  a matching declaration downstream. A suffix nothing declares is a suffix
  nothing reads.
- **Adding a numeric column** requires a `must_vary` setting. If it is
  legitimately constant, `Vary.NEVER` demands a written `constant_reason` —
  silencing a degeneracy finding should cost a sentence the next person can
  judge.
- **Changing a value's semantics** (a boundary, a denominator, a units change)
  is a release-note event. Precedent: the region-MDS E1 entry, which carries an
  explicit "values are not comparable" warning.

## Things that look like bugs and are not

- **WPS is Parquet-only by design.** Do not "fix" it to honour
  `--output-format`.
- **`ultra_long` (401-1000bp) is deliberate**, added for necrosis and fetal
  cfDNA detection. The 1000bp ceiling is load-bearing.
- **Channels are non-overlapping on purpose** — overlapping ones would make the
  ML features collinear.

## Known divergences

Recorded in `src/krewlyzer/validate/claims.py` under `KNOWN_DIVERGENCES`, which
asserts they are *still* present so that fixing one fails the test and forces
the entry out rather than leaving a stale note. Two today: the GC length-bin
ceiling, and `mono_nucl` meaning 150-220bp genome-wide but 150-259bp per gene.
