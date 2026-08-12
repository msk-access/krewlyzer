# Working conventions

How changes are expected to be made here. Not style rules — those are in
`.agents/rules/development.md` — but the habits that catch this codebase's
characteristic failures.

## Review before and after each step

The after-step is the one that matters. Several fixes in this repository
introduced a fresh defect that only a second look caught:

- correcting reverse-strand fragment coordinates left the BED **unsorted**,
  which breaks tabix and takes every downstream feature with it;
- a diagnostic warning added for "no silent failures" **deadlocked the
  pipeline**, because it logged per fragment from inside the rayon loop and
  pyo3-log routes into Python under the GIL;
- adding columns to the gene BED assets silently **disabled strand handling**
  in `region-mds`, because format detection keys on column count.

None of these were in the change being reviewed. All were introduced by it.

## No second copy of anything

Every large defect found here was a value that existed twice and drifted:

| duplicated | outcome |
|---|---|
| FSC size bands, in the genome and gene paths | a column named `di_nucl` meant two different ranges |
| the version string, in three files | maintained by four hand-written regexes |
| the MDS scale, in two sections of one page | the page contradicted itself for a year |
| the whole `docs/` tree, in `krewlyzer_all_docs.md` | four documents missing |

The fix is always to delete the copy, not to correct it. Prefer a shared
function, a generated file, or a build-time asset over an edit in two places.

## No silent failures

A metric that cannot vary with the data is worse than a missing one — it is
present, plausible, and gets modelled. Related: a value that means "we could
not determine this" must not be spelled the same way as a real measurement.
`mds_e1` used `0.0` for both "E1 had no fragments" and "there is no E1", on a
scale where lower means more abnormal; it now uses `NaN` for the second.

Where a flag can be unknown, spell it `.`, not `0`.

## A test proves nothing until it has failed

Revert the fix, watch the test fail, restore the fix. Three tests in this
repository passed against broken code before this was enforced. See
`.agents/rules/validation-gates.md` for the examples and the mechanics.
