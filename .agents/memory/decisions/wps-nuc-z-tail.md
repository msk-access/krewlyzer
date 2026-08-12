# `wps_nuc_z` has no calibrated tail, and shrinkage will not give it one

**Decided 2026-08-11.** Applies to the per-position `wps_nuc_z` / `wps_tf_z`
vectors, and to the `WPS.PLAUSIBLE_Z_SCORES` finding that survives 0.9.0.

## The question

0.9.0 removed a real divisor defect: σ values that were floating-point residue
(~1e-17) rather than measured spread, which produced z up to 6.1 × 10¹⁸. After
that fix the gate still reported implausible z on WPS, and 9.1% of positions
exceeded |z| = 3 where a normal distribution expects 0.27%.

Three explanations were open: a second divisor defect, an estimator problem
that variance shrinkage would fix, or the data. Getting this wrong in either
direction is expensive — shipping a known-bad divisor, or adding a statistical
method to the product to silence an alarm that was telling the truth.

## The numbers that drove it

Measured on `xs1.all_unique` against one XS1 plasma sample, 5.0 M positions:

| | observed | N(0,1) expects |
|---|---:|---:|
| median \|z\| | 0.447 | 0.674 |
| \|z\| > 2 | 12.79% | 4.55% |
| \|z\| > 3 | 9.13% | 0.27% |
| \|z\| > 5 | 4.92% | 0.0001% |
| \|z\| > 10 | 1.4% | ~0 |
| \|z\| > 100 | 4,808 | 0 |

A core that is *too narrow* alongside a tail that is *too fat* is the mixture
signature shrinkage is supposed to correct. It is also **not** the residue-σ
defect: only **14** of those 4,808 positions have σ < 1e-4.

## Why not the alternatives

**Variance shrinkage (empirical Bayes, d₀ = 5).** The textbook remedy, and it
was costed against the cohort rather than assumed. It fixes only the extreme
tail and leaves the actual problem untouched:

| | current | eBayes d₀=5 | target |
|---|---:|---:|---:|
| \|z\| > 100 | 4,808 | **211** ✓ | 0 |
| \|z\| > 5 | 4.92% | 4.47% | 0.0001% |
| \|z\| > 2 | 12.79% | 12.49% | 4.55% |
| median \|z\| | 0.447 | 0.425 ✗ | 0.674 |

The median moves *away* from 0.674. **If this were a variance-estimation
problem, shrinkage would fix it. It doesn't, so it isn't.** It would have been
cheap — one builder function, a d₀ fit, and it rides the re-aggregation — which
is exactly why it needed measuring rather than adopting.

**Anchor-level offset.** Removing each anchor's own median z makes it worse:
|z| > 2 rises from 12.8% to 13.6%.

**Depth.** Correlation between anchor offset and `local_depth` is +0.022.

**Whole bad anchors.** 64.8% of anchors have no position beyond |z| = 10, and
none have more than half — so it is a position-level subset, not a set of
anchors to exclude.

## What it actually is

Each anchor's baseline is fitted from a **median of 8 donors**. The xs1 cohort
is 47, but genome-wide coverage is patchy and an anchor is kept once ≥ 3 donors
cover it. At n = 8, `(x − μ)/σ` is not a z-score — even a t₇ expects ~2% beyond
|z| = 3, against the 9.1% observed — and WPS at a position is a difference of
coverage-dependent weighted counts, which is genuinely heavy-tailed.

## Consequences

**The values are unchanged.** They are what the data supports.

The `PLAUSIBLE_Z_SCORES` finding on WPS stays, and its message — *"a z that
large is a near-zero sigma, not a measurement"* — is demonstrably false for
LIST columns. The check needs its claim corrected, not its threshold raised.

The numbers above are **not** in `validate/claims.py`. That registry exists to
catch doc-versus-code drift in constants; these are measurements against one
PON and one sample, and pinning them there would assert a stability they do not
have. They are dated in prose instead, so they are falsifiable.

For a calibrated per-anchor statistic use `wps_shape_corr_z` and
`wps_log_amplitude_z`, which are scalars fitted against baselines of themselves
and are unaffected by any of this.

## What would reopen this

**More donors per anchor.** That is the only remedy that addresses the cause,
and it is a cohort-design question rather than a code one — worth knowing
before anyone builds a model on `wps_nuc_z`.

These figures are one sample against one PON, and a healthy-donor run gave a
different rate (0.83% beyond |z| = 100 against 0.095%), confounded by both
version and cohort. Re-measure before quoting them; do not treat the table
above as a constant.
