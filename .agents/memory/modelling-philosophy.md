# What to emit, and what not to decide here

krewlyzer runs upstream of kreview's model fitting. Deciding a signal is
uninteresting *here* removes it from consideration permanently, and kreview
cannot tell the difference between "absent because uninformative" and "absent
because broken" — every reader there swallows exceptions and yields an empty
feature dict.

So: **measure everything, and let the modelling step decide what is relevant.**

## What that means in practice

**Prefer an extra column to a collapsed one.** When a value is a boundary hit
rather than a measurement, add a flag beside it instead of dropping or clamping
the row — `nrl_at_band_limit` exists because `nrl_bp = 250.0` could mean either
"the repeat length is 250bp" or "no peak was found", and a consumer had no way
to tell.

**When two notions could share a name, emit both under distinct names.**
`is_e1` (the canonical transcript's exon 1) and `is_alt_e1` (an alternative
promoter) are separate columns for this reason. Collapsing them would have
forced a choice between discarding two thirds of the captured signal and
labelling internal exons as promoter-proximal.

## The limit

This is not a licence to emit anything. **Never mislabel.** Widening a
definition and annotating it is right; renaming something into a claim it
cannot support is not.

`is_first_captured` exists for every gene, and it was tempting to call it E1 so
that the E1 table stayed full. It is usually an internal exon. Calling it E1
would have asserted promoter proximity that is not there — inventing signal
rather than measuring it.

The test: if a downstream model treated this column exactly as its name
suggests, would it be misled? If yes, rename it or add the flag that
disambiguates.
