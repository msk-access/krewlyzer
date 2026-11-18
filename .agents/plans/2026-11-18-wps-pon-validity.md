> **SUPERSEDED (18 Nov 2026)** by `2026-11-18-integrated-plan.md`, which folds this plan into
> the fragment-model release (`2026-11-18-fragment-model-spec.md`). Phase 1 (#83) is done and
> carried over; the #80/#77 tasks (7-12) move to the integrated plan's Step 6 unchanged in
> substance. Kept for its review findings and decision rules, which the integrated plan cites.

# WPS coordinates, amplitude depth-independence, and PON scorability — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the WPS coordinate convention (#83), remove the depth dependence of `wps_log_amplitude` (#80), and stop WPS z-scores being computed against PON baselines that cannot serve as a yardstick (#77) — in one release, with one per-sample recompute and one PON rebuild.

**Architecture:** Phase 1 changes only per-sample WPS computation (Rust, local). Phase 2 is the single per-sample recompute on that code, from retained fragment BEDs, followed by measurements that fix the two numbers Phase 3 needs. Phase 3 changes only PON building and scoring, so it consumes Phase 2's outputs without another per-sample pass over donors. Phase 4 rebuilds the four PONs and releases.

**Tech Stack:** Rust (`rust/src/wps.rs`, `rust/src/pon_builder.rs`, pyo3), Python (`src/krewlyzer/core/wps_pon.py`, `src/krewlyzer/validate/`), pytest, cargo test, maturin.

**Spec:** issues [#83](https://github.com/msk-access/krewlyzer/issues/83), [#80](https://github.com/msk-access/krewlyzer/issues/80), [#77](https://github.com/msk-access/krewlyzer/issues/77), amended by **Review findings** below. Read the issues; this plan argues from them.

## Decisions to confirm before Task 1

Each has the recommended option the plan assumes. A different answer changes the tasks named.

| | decision | recommended | alternative | changes |
|---|---|---|---|---|
| D1 | WPS protection window | **symmetric, 2p+1 bases `[x−p, x+p]`** — the only convention under which a fragment's contribution is mirror-symmetric, which matters because minus-strand profiles are reversed before modelling | 2p bases, accepting a half-base offset between strands | Task 2 |
| D2 | Fix binned strand reflection in this release | **yes** — same rebuild, no extra cost | defer, record in #83 | Task 3 |
| D3 | #80 fix | **depth-normalise if it passes M3's gate; otherwise document the dependence** | condition the PON baseline on depth | Task 10 |
| D4 | #77 mechanism | **withhold z with a label, keyed on donor support** (the issue's agreed requirements 1–3) | predictive variance `sqrt(σ² + σ_noise²)`: no exclusions, but needs a calibrated noise model in WPS units | Tasks 7, 9 |
| D5 | Name of the changed amplitude feature | **keep `wps_log_amplitude`**; record the meaning change in the contract and CHANGELOG, tell kreview | rename — but kreview reads by name and a missing column becomes an empty feature silently (invariant #2) | Task 10 |
| D6 | Branching | **one integration branch `feature/pon-validity`, merged to `develop` once, with the rebuilt PONs** | merge Phase 1 alone — leaves `develop` scoring new-semantics samples against old-semantics PONs | all |

## Review findings that amend the issues

1. **#83 — the index is 0-based; the WPS terms were written as if it were 1-based.** `parse_regions_with_window` builds `[center − 1000, center + 1000]` inclusive from a 0-based anchor base (`len = end − start + 1 = 2001`), so the anchor is at index 1000. `f_start = fragment.start + 1` places every fragment one base right. Under D1 the per-stream deltas are: coverage loses the fragment's first base (confirmed: 2200 not 2220 in `test_local_depth_is_mean_coverage_over_the_window`); spanning is shifted +1 with correct width — **not** "one base narrower" as #83 says; both endpoint windows are one base narrow and the end-side one is shifted. Update #83 in Task 2.
2. **#83 — binned strand reflection is about the wrong point.** Bins cover indices 0–1999 (index 2000 is never binned) and `binned_*.reverse()` reflects about 999.5, not the anchor at 1000. The legacy per-base path reverses 2001 rows and is exact.
3. **#77 — the rebuild is not "minutes, no recompute" any more.** #79 changed `local_depth` and #83 changes WPS values, so donor outputs must be recomputed. If fragment BEDs are retained (the #80 measurement used them), that is a BED pass, not a BAM pass.
4. **#77 — the |μ| cut conflated two populations.** WPS is a raw weighted count difference (binned `wps_nuc` is the mean of raw WPS over the bin, not depth-normalised), so |μ| ≈ 0 both where there is no coverage and in well-covered transition zones where spanning ≈ ending. Keying on donor coverage may withhold far less than the third of positions the |μ| cut discarded. M2 tests this.
5. **#80 — three implementations.** `pon_builder.rs::log_amplitude` (f32, `(1.0 + range).ln()`), `wps.rs::wps_log_amplitude` (f64, `ln_1p`), and the Python oracle in `tests/unit/test_rust_python_equivalence.py`. In f32, `1.0 + x` is exactly 1.0 below about 6e-8, so near-zero-coverage anchors can get a donor baseline of exactly 0 while the sample side gets a small positive value — in the regime #77 is about. `pon_builder.rs:374-377` still asserts the debunked `+0.512 → −0.036`.

## Global Constraints

- Rebuild the extension with `maturin develop --release` **from the project root** before any `pytest`; otherwise the suite tests the previous binary.
- Gates per change: `cargo test --manifest-path rust/Cargo.toml --lib`; `cargo clippy --manifest-path rust/Cargo.toml --tests -- -D warnings -A clippy::too_many_arguments -A clippy::type_complexity`; `pytest tests/ -q`; `bash scripts/lint.sh`; `python scripts/check_output_format.py` for output-writing changes; `pytest tests/unit/test_claims_registry.py` for constants and docs.
- **One per-sample recompute.** Phase 3 must not change any value computed in the per-sample pass (`wps_nuc`, `wps_tf`, `prot_frac_*`, `local_depth`, `capture_mask`). A task that would must say so, and the plan then needs a second pass.
- `MIN_PON_VERSION` moves off `(0, 9, 0)`; all four PONs are **rebuilt, not re-stamped** (release guide Phase 2.8).
- New or changed columns: listed in `src/krewlyzer/validate/contract.py` with a `Vary` rule (invariant #1), in `docs/reference/output-files.md`, and identical between the single `wps` tool and `run-all` (invariant #6). Exclusions are labelled in their own column, never a bare NaN (invariant #3). Numeric constants a doc quotes live in `src/krewlyzer/validate/claims.py` (invariant #5).
- No patient identifiers, cluster paths or operator-specific settings in any committed file; use `P-0000000` / `C-000000-L000`. HPC steps follow the operator's own rules, which are not in this repository: every read and copy is approved first.
- Measurements are reported as counts and summary statistics, never per-sample records.

## Review Focus

1. **Anchors within 1 kb of a chromosome start.** `center.saturating_sub(window_bp)` truncates the window, so the anchor is not at index 1000; bins and strand reflection both shift. Expected: the anchor's own base lands in the same bin as everywhere else, and positions before the chromosome start are NaN, not 0, in both orientations. Test owned by Task 3.
2. **Fragments straddling a window edge.** `add_range` clamps; a clamp must clip, never shift. Expected: a fragment starting 50 bases before the window contributes exactly its in-window bases. Test owned by Task 2.
3. **Fragments of length 2p and 2p+1.** Expected: length 2p spans nothing; length 2p+1 spans exactly one position, the midpoint. Test owned by Task 2.
4. **A sample with zero coverage at a well-supported anchor.** `(0 − μ)/σ` is a measurement of a failed capture, not biology. Expected: withheld and labelled, not a large negative z. Test owned by Task 9.
5. **σ NaN (below `SIGMA_FLOOR`) at an anchor that passes the support criterion.** Expected: withheld and labelled, so "z is null" always means "absent or withheld", never an unexplained gap. Test owned by Task 9.

---

## Phase 1 — #83, local

### Task 1: Extract fragment accumulation into a testable function

**Files:**
- Modify: `rust/src/wps.rs` — the per-fragment block inside the region loop (~1221-1272)

**Interfaces:**
- Produces: `fn accumulate_fragment(acc: &mut RegionAccumulator, region_start: u64, frag_start: u64, frag_end: u64, nuc_weight: f64, tf_weight: f64, nuc_window: i64, tf_window: i64)`. `frag_start`/`frag_end` are the fragment's 0-based half-open coordinates exactly as read from the BED; a weight of `0.0` means that stream ignores the fragment. Behaviour identical to today.

- [ ] **Step 1:** Move the block into `accumulate_fragment` unchanged and call it from the loop.
- [ ] **Step 2:** Run `cargo test --manifest-path rust/Cargo.toml --lib` — expect 121 passed, unchanged.
- [ ] **Step 3:** Rebuild, then `pytest tests/integration/test_wps_cli.py -q` — expect `test_local_depth_is_mean_coverage_over_the_window` still passing at `2200/2001`. That fixed value is what proves the refactor changed nothing.
- [ ] **Step 4:** Commit `refactor: extract WPS fragment accumulation into a function`.

### Task 2: One stated coordinate convention, asserted per stream

**Files:**
- Modify: `rust/src/wps.rs` — `accumulate_fragment`, new `#[cfg(test)] mod coordinate_convention_tests`
- Modify: `tests/integration/test_wps_cli.py` — the `2200/2001` expectation and its derivation comment
- Modify: `docs/features/core/wps.md` — state the convention; `src/krewlyzer/validate/claims.py` — `wps.nuc_window = 60`, `wps.tf_window = 8`
- Modify: `CHANGELOG.md` `[Unreleased] / Fixed`

**Interfaces:**
- Consumes: `accumulate_fragment` (Task 1)
- Produces: the convention, written once in the doc comment of `accumulate_fragment`: index `i` is 0-based genomic position `region_start + i`; a fragment occupies `[frag_start, frag_end)`; the protection window at `x` is `[x − p, x + p]`. In region-relative half-open terms, with `a = frag_start − region_start`, `z = frag_end − region_start`: coverage `[a, z)`; spanning `[a + p, z − p)`; start endpoint `[a − p, a + p + 1)`; end endpoint `[z − 1 − p, z + p)`.

- [ ] **Step 1: Write the failing tests** (one fragment, `a = 100`, `z = 267`, nucleosome weight 1.0, `p = 60`; read each array's running sum):
  - `a_fragment_covers_exactly_its_own_bases`: coverage is 1.0 on 100..=266, 0.0 at 99 and 267, total 167.0
  - `spanning_is_where_the_whole_window_is_covered`: spanning is 1.0 on 160..=206, 0.0 at 159 and 207
  - `each_endpoint_window_is_centred_on_its_endpoint`: start window 1.0 on 40..=160; end window 1.0 on 206..=326
  - `a_fragments_contribution_is_mirror_symmetric`: for every array, `value[m + d] == value[m − d]` about `m = (a + z − 1) / 2` — the convention-free property
  - Review Focus 2 — `a_fragment_straddling_the_window_edge_is_clipped_not_shifted`: `frag_start = region_start − 50`, length 167 → coverage on `0..=116` exactly
  - Review Focus 3 — `length_2p_spans_nothing_and_2p_plus_1_spans_its_midpoint`
- [ ] **Step 2:** Run `cargo test ... coordinate_convention_tests` — expect all six to fail, the first with coverage `0.0` at index 100.
- [ ] **Step 3:** Change the five `add_range` calls in `accumulate_fragment` to the intervals in Interfaces. No other change.
- [ ] **Step 4:** Run the six tests — expect PASS; full `cargo test --lib` — expect PASS.
- [ ] **Step 5:** Change the integration expectation to `2220 / 2001` and rewrite its derivation comment: each fragment contributes `length` bases; remove the `(length − 1)` explanation. Rebuild, run `pytest tests/integration/test_wps_cli.py -q` — expect PASS.
- [ ] **Step 6:** Audit, without fixing: `grep -nE "start \+ 1|\.start \+ 1" rust/src/*.rs` and any region-relative coordinate arithmetic in `ocf.rs`, `region_entropy.rs`, `fsc.rs`. Record each finding as a comment on #83, with whether it uses the same convention.
- [ ] **Step 7:** State the window as "±60 bp (121 bases)" and "±8 bp (17 bases)" in `wps.md`; add both claims; run `pytest tests/unit/test_claims_registry.py`.
- [ ] **Step 8:** Edit #83 with finding 1 (the convention, and the corrected per-stream deltas).
- [ ] **Step 9:** Commit `fix: one stated WPS coordinate convention, asserted per stream`.

### Task 3: Reflect minus-strand profiles about the anchor (D2)

**Files:**
- Modify: `rust/src/wps.rs` — binned output path (~920-990)

**Interfaces:**
- Produces: `fn anchor_index(region: &Region) -> usize` = `region.center − region.start`. Bins are defined by **offset from the anchor base**, not by array index: bin `k` covers offsets `[10·(k − 100), 10·(k − 100) + 10)`. In a full window (`anchor_index == 1000`) this is exactly today's array indices 0–1999; in a window truncated by a chromosome start it keeps the anchor in bin 100, and offsets with no array position are NaN. Minus-strand per-base arrays are reflected about the anchor (`j ↦ 2 · anchor_index − j`) **before** binning, and the bin `reverse()` calls are removed.

- [ ] **Step 1: Write the failing tests:**
  - `opposite_strands_of_a_mirror_profile_bin_identically`: a plus-strand region with an asymmetric per-base profile and a minus-strand region carrying its mirror image produce equal 200-bin vectors for `wps_nuc`, `wps_tf`, `prot_frac_nuc`, `prot_frac_tf`, `capture_mask`
  - Review Focus 1 — `an_anchor_near_a_chromosome_start_keeps_its_bin`: `center = 500` → the anchor's base lands in bin 100 on both strands; bins entirely before position 0 are NaN
- [ ] **Step 2:** Run — expect the first to fail by one position, the second to fail on bin placement.
- [ ] **Step 3:** Implement the reflection as specified.
- [ ] **Step 4:** Run the two tests and `cargo test --lib` — expect PASS.
- [ ] **Step 5:** Rebuild; `pytest tests/ -q`; `bash scripts/lint.sh`; `python scripts/check_output_format.py`.
- [ ] **Step 6:** Commit `fix: reflect minus-strand WPS about the anchor, not the window midpoint`.

---

## Phase 2 — one per-sample pass, then measure (HPC, operator-approved)

### Task 4: Confirm the pass can start from fragment BEDs

**Files:** none in the repository. A local manifest only.

- [ ] **Step 1:** With the operator's approval, list for each of the four PONs whether every donor's fragment BED is retained, and the same for a patient subset of 20 per assay (the selection used for #77's measurement 5).
- [ ] **Step 2:** Confirm which entry point recomputes features from a fragment BED on the `feature/pon-validity` build, and that it writes the `--from-outputs` cache that `build-pon` accepts.
- [ ] **Step 3: Stop condition.** Any missing BED means that sample needs a BAM pass. Report the count to the user and wait; do not start one.

### Task 5: Recompute, build provisional PONs

- [ ] **Step 1:** Recompute per-sample features for all donors of all four PONs and the patient subset, on the Phase 1 code.
- [ ] **Step 2:** For each output: `krewlyzer validate-output` — no new error class; `local_depth ≥ 0` on every row.
- [ ] **Step 3:** `build-pon --from-outputs` for each PON, current scoring semantics. These are provisional and are never shipped.
- [ ] **Step 4:** Donor invariant (#77 measurement 6) on each provisional PON: observed/required `RMS(z)` within 5% and zero donor `|z| > 100`. A failure means build and score disagree on the new code: stop.

### Task 6: Measure #83's effect size (M1)

- [ ] **Step 1:** Per donor, old vs new `wps_nuc` vectors: report the median `|Δ| / σ` and the share of anchors whose maximum bin moves. Informational; it says how much this release changes WPS.

### Task 7: Fix #77's support threshold (M2)

**Interfaces:**
- Produces: `S*`, a constant named `WPS_MIN_DONOR_SUPPORT`, used by Task 9.

- [ ] **Step 1:** Per anchor, donor support `S` = median over donors of `local_depth`. Bin anchors by `S` decade.
- [ ] **Step 2:** For each cell (anchor × bin), the cross-sample SD of `wps_nuc_z` over the patient subset, scored against the provisional PON. Per decade and assay report: `p99` of that SD over the decade's cells; the share of cells with any `|z| > 100`; the share of genome-wide anchors. `plateau` = the same `p99` taken over all cells in the top three decades together. `inflation(decade) = p99 / plateau`.
- [ ] **Step 3: Decision rule.** `S*` = the lower edge of the lowest decade such that **it and every higher decade** have `inflation ≤ 2`, in **both** assays. The tail statistic, not the median: the defect lives in the tail, and #77's own tables show central tendencies understating it. "Every higher decade" because a single low decade can pass by chance — the earlier xs1 run showed 0.30 next to 1.69.
- [ ] **Step 4: Granularity check.** Share of `|z| > 100` cells lying in anchors with `S < S*`. If ≥ 95%, a per-anchor criterion suffices. If not, stop: the criterion needs per-bin support, which needs a new per-sample column — a second per-sample pass — and goes back to the user.
- [ ] **Step 5: Affordability check.** Share of anchors withheld at `S*`, next to the |μ| < 0.1 cut on the same anchors. If `S*` withholds as much or more, stop and report: finding 4's hypothesis failed and D4 returns to the user.

### Task 8: Choose #80's form (M3)

- [ ] **Step 1:** Per anchor, compute `c0 = ln(1 + range)` (today) and `c1 = ln(1 + range / local_depth)`, where `range = max(wps_nuc) − min(wps_nuc)`.
- [ ] **Step 2: Gate:** `|Spearman(c, local_depth)| ≤ 0.1` within donors and pooled with patients, in both assays, and nonzero cross-sample variance.
- [ ] **Step 3:** Report donors-vs-patients AUC for `c0` and `c1`, not gated: if donors and patients differ in depth, `c0`'s separation is partly that confound, so a drop is expected.
- [ ] **Step 4: Decision rule.** `c1` passes → Task 10 implements it. `c1` fails → Task 10 documents the dependence instead (D3).
- [ ] **Step 5:** Write `.agents/memory/decisions/2026-11-XX-wps-pon-validity.md` following the existing records there: `S*` with the M2 table, the M3 table and choice, M1, M4. Commit it.

---

## Phase 3 — PON building and scoring only (local)

### Task 9: Withhold and label z where the baseline cannot serve (#77)

**Files:**
- Modify: `rust/src/pon_builder.rs` — store per-anchor donor support in the WPS baseline
- Modify: `rust/src/wps.rs` — `apply_pon_zscore_inner`
- Modify: `src/krewlyzer/validate/contract.py`, `src/krewlyzer/validate/claims.py`, `docs/reference/output-files.md`, `docs/features/core/wps.md`, `CHANGELOG.md`
- Test: `rust/src/wps.rs` `#[cfg(test)] mod scorability_tests`; `tests/invariants/test_antidegeneracy.py`

**Interfaces:**
- Consumes: `S*` (Task 7)
- Produces: PON column `wps_baseline.donor_support: f32` (median donor `local_depth`); const `WPS_MIN_DONOR_SUPPORT: f32 = S*`, passed to `apply_pon_zscore_inner` as a `min_donor_support: f32` parameter so a test can set it to `0.0` and recover the old behaviour; output column `wps_z_withheld: bool`, one per anchor. Every other feature keeps scoring against `SIGMA_FLOOR` exactly as today — this criterion is WPS's alone (#77 requirement 1). True when donor support < `S*`, the sample's `local_depth == 0`, or σ is NaN; then every `*_z` vector for that anchor is null. Null z with `wps_z_withheld == false` means the anchor is absent from the PON.

- [ ] **Step 1: Write the failing tests:**
  - `a_baseline_no_donor_covered_is_not_a_yardstick`: a synthetic PON where every donor has zero coverage at one anchor, and a sample with trace coverage there → `wps_z_withheld`, z null. This reproduces #77's `|z| = 9.09e3`: also assert that with `min_donor_support = 0.0` the same input gives `|z| > 100`, so the test is shown to bite.
  - Review Focus 4 — `a_sample_with_no_coverage_is_withheld_not_scored`
  - Review Focus 5 — `an_unmeasurable_sigma_is_withheld`
  - `a_supported_anchor_scores_exactly_as_before`: `withheld == false`, z equals `(x − μ) / σ` to 1e-6
- [ ] **Step 2:** Run — expect the first three to fail.
- [ ] **Step 3:** Implement as specified.
- [ ] **Step 4:** Run, plus `cargo test --lib` — expect PASS.
- [ ] **Step 5:** Register `wps_z_withheld` in the contract with a `Vary` rule; add the antidegeneracy case; add `wps.min_donor_support` to claims; make `PLAUSIBLE_Z_SCORES` ignore withheld anchors, with a test.
- [ ] **Step 6:** Rebuild; full gates; confirm `wps` and `run-all` write the identical column.
- [ ] **Step 7:** Commit `fix: withhold WPS z where no donor covered the anchor`.

### Task 10: One log-amplitude, depth-free if M3 allows (#80)

**Files:**
- Modify: `rust/src/pon_builder.rs` (remove `log_amplitude`, use the shared one; delete the stale comment at 374-377), `rust/src/wps.rs` (`wps_log_amplitude`), `tests/unit/test_rust_python_equivalence.py`, `docs/reference/output-files.md`, `CHANGELOG.md`

**Interfaces:**
- Produces: one `pub(crate) fn wps_log_amplitude(profile: &[f64], local_depth: f64) -> f64` used by both builder and scorer. Under `c1`: `ln_1p(range / local_depth)`, NaN when `local_depth ≤ 0` or fewer than two finite points. Under the documented fallback: today's formula, with `local_depth` unused.

- [ ] **Step 1: Write the failing tests:**
  - `builder_and_scorer_compute_the_same_amplitude`: bit-identical on the same profile, including `range = 1e-8` — fails today on f32 vs f64
  - Under `c1` — `amplitude_is_invariant_to_uniform_depth_scaling`: profile × k with `local_depth` × k gives the same value for k ∈ {0.1, 10, 1000}. This is the unit-level form of "no rank correlation with depth".
  - Under the fallback — rename `the_log_form_compresses_the_depth_dependence` to state the dependence it pins, and say so in its comment.
- [ ] **Step 2:** Run — expect failures.
- [ ] **Step 3:** Implement; update the Python oracle to the same definition.
- [ ] **Step 4:** Run; rebuild; full gates.
- [ ] **Step 5:** Record the meaning change of `wps_log_amplitude` in `output-files.md` and `.agents/rules/output-contract.md` (D5).
- [ ] **Step 6:** Commit `fix: one wps_log_amplitude for builder and scorer`.

### Task 11: The donor invariant as a permanent check (#77 requirement 6)

**Files:**
- Modify: the `validate-pon` command; Test: `tests/unit/test_pon_donor_invariant.py`

**Interfaces:**
- Produces: `krewlyzer validate-pon PON --donor-outputs DIR`. For each vector block, score the donors against the PON and report observed/required `RMS(z)` with required `sqrt((n − 1) / n)`. Fails when the ratio is outside `[0.90, 1.02]` or any donor `|z| > 100`.

- [ ] **Step 1: Write the failing test** `rms_z_of_the_donors_is_sqrt_n_minus_1_over_n`: build a PON from 12 synthetic donors, score the same donors, assert `RMS(z) == sqrt(11/12)` to 1e-6.
- [ ] **Step 2:** Run — expect failure (no `--donor-outputs`).
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run; full gates.
- [ ] **Step 5:** Commit `feat: validate-pon asserts the donor RMS(z) invariant`.

### Task 12: Assert the premise (#77 requirements 4 and 5)

**Files:**
- Test: `tests/invariants/test_z_spread_commensurate.py`; `tests/real_data/` (runs locally with `KREWLYZER_TEST_CORPUS`)

- [ ] **Step 1:** Synthetic invariant `non_withheld_z_spread_is_commensurate`: a cohort with mixed donor support → cross-sample SD of every non-withheld z ≤ 50.
- [ ] **Step 2:** Real-data test: over non-withheld cells, p99 SD(z)/plateau ≤ 2 per assay — the same quantity as Task 7's rule, so the criterion and its check cannot drift apart.
- [ ] **Step 3:** Run both; commit `test: assert the scorability premise, not just its constant`.

---

## Phase 4 — rebuild and release (HPC, operator-approved)

### Task 13: Rebuild, validate, release

- [ ] **Step 1:** Bump `MIN_PON_VERSION`; `build-pon --from-outputs` for all four PONs from Task 5's donor outputs, on the final code. No second donor pass: Phase 3 changed only building and scoring.
- [ ] **Step 2:** `validate-pon --donor-outputs` on each — pass.
- [ ] **Step 3:** Re-score the patient subset (a BED recompute of 40 patients, or a re-score entry point if Task 4 found one); `validate-output` on each. Expected: xs1 `PLAUSIBLE_Z_SCORES` findings close to 0. Report xs2's separately — the issue's asymmetry is out of scope.
- [ ] **Step 4:** Follow `docs/development/release-guide.md` for the data change: bundled PONs, data-image pin, `full-assets.yml`.
- [ ] **Step 5:** Tell kreview: `wps_z_withheld` is new; `wps_log_amplitude` changed meaning (if `c1`).
- [ ] **Step 6:** Merge `feature/pon-validity` to `develop` (D6).

### Task 14: Close out

- [ ] **Step 1:** Comment the measured results on #83, #80 and #77, then close them.
- [ ] **Step 2:** Open a new issue for the xs1/xs2 asymmetry: it needs external healthy samples scored against the xs2 PON, which this plan does not have.
