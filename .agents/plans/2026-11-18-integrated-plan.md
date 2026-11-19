# Integrated plan — fragment model + PON validity (#83, #80, #77)

**Revision 2, 18 Nov 2026.** Revision 1 was reviewed adversarially by a fresh reviewer; its
findings were checked against the code and are folded in below. Where revision 1 was wrong, this
says so rather than quietly changing it.

**Sources, in order of authority**

1. `2026-11-18-fragment-model-spec.md` — the maintainer's plan, **verbatim**. Authoritative for
   everything it covers, except where §1 below records a maintainer decision that changes it.
   This document maps the spec; it does not restate it.
2. `2026-11-18-wps-pon-validity.md` — superseded. Its #80/#77 decision rules (Tasks 7–14) are
   cited from there, with the corrections in §4 below.
3. This document — the decisions taken since, where the PON-validity work lands, the
   interactions between the two, and the measurement design.

Branch: `feature/fragment-model`, off `develop`, per spec §3. Each step is a PR into it.

---

## 1. Decisions (18 Nov 2026)

| Topic | Decision | Changes the spec? |
|---|---|---|
| Scope | One release: everything in the spec, plus #83, #80, #77. No kreview runs until it ships | — |
| PONs | All four rebuilt: xs1/xs2 × All Unique/Duplex | Spec §3 Step 6 named All Unique only |
| `fragments.parquet` | **Off by default.** `--fragment-export` on `extract` and `run-all` writes it in the same BAM pass as the BED. No separate command: `extract --fragment-export` already is a BAM pass without features. Refuses to start if free space at the output is below the estimate (bytes per fragment from Step 0 × expected fragments). The maintainer confirms space before any cohort-wide export | Spec §2e: "the cohort configuration turns it on" |
| Measurement harness | Test samples are extracted **once**, with `--fragment-export --emit-bed`. Each definition or policy is then evaluated as a projection: parquet → BED under that policy → feature engine. **Measurement only**; production runs still start from the BAM | Exception to spec §0 "every run starts from the BAM", for measurement |
| Re-score | New command: recompute z-scores for an existing per-sample output against a named PON, from the raw values already stored, with no BAM | New |
| Molecule length | **In both**: a fifth BED column and the parquet. A four-column BED still loads, falling back to `end − start` and logging that it did. Which length each feature uses stays spec §2d row 5's measured decision, applied to FSD and WPS stream assignment alike | Spec §2e: BED stays four columns |
| Refused PON | **Fails loudly.** Today `load_pon_model` logs and returns `None`, and the run finishes without z columns | New |
| Test set | Samples from two cohorts — one sequenced at IGO (1,579 samples), and krewlyzer output from clinical ACCESS BAMs mirrored for research — stratified by source × assay | New |
| Final run | 16,000+ samples from clinical ACCESS runs | Spec §5 sized storage for 16,000 |
| BAM access | SLURM jobs read clinical BAMs in place, in parallel. They run as the operator | — |
| Cluster work | **Everything runs as a SLURM job** — env builds, probes, inventories, measurements. Nothing runs on the login node | — |
| Healthy donors | Sequenced on the **clinical** side, and used for research too. Every cohort goes through the **same processing workflow**; the only batch difference is the sequencing site. The final 16,000-sample run therefore shares the PON donors' batch; the IGO test patients do not (§4) | — |

## 2. Every spec section, and where it lands

| Spec § | Content | Step |
|---|---|---|
| §0 | Maintainer decisions | Global constraints, as amended by §1 |
| §1 F1 | Mates never joined; global motifs and region MDS count different ends | 1, 3, 4 |
| §1 F2 | Fragment record is a formatted string | 1 |
| §1 F3 | GC 0.0 doubles as "unknown" | 3 (row 1) |
| §1 F4 | Soft clips ignored at fragment ends | 1, 3 (row 7) |
| §1 F5 | GC length bins capped at 400 bp vs 188 declared | 3 (row 2) |
| §1 F6 | Whole-genome buffering; linear scan over exclude regions | 2 |
| §1 F7 | Strand-family counts summed in mFSD | 1, 4 |
| §1 F8 | Off-target = "no target overlap"; near-target pooled | 0, 3 (row 9) |
| §1 F9 | BED is an intermediate; tabix never queried | 5 |
| §1 F10 | Run-all cache dead and unsafe | 5 |
| §1 F11 | Molecule length = reference span | 1, 3 (row 5), 5 (fifth column) |
| §1 "Not verified" | Target BED probe vs exon coordinates; `output_provenance.py` unread | 0 |
| §2a | `fragments.parquet` schema, groups A–F | 1 (record), 5 (writer) |
| §2b | `MatePairer` | 1 |
| §2c | `FilterPolicy` for all four BAM readers | 1, 4 |
| §2d | Value-changing policies, 9 rows | 3 (rows 1, 2, 4–9), 4 (row 3: region MDS reads the BAM until migrated) |
| §2e | Outputs, BED lifecycle, cache removal, Nextflow, metadata counts | 5, as amended by §1 |
| §2f | Acceptance rule | Every value-changing commit; test 1 not used for judgment calls (§5) |
| §2g | Region classes; penumbra map | 5 (mechanism), 6e (asset) |
| §2h | Provenance; GC-bias summary | 5 |
| §2i | Runtime and memory budget | 0 (baseline), 2 and 5 (report) |
| §3 | Steps 0–6 | §3 below |
| §4 | Pre-merge tests | Gate before merge; additions in §5 |
| §5 | Consequences | Stated now, in §8 — the spec says state them before starting |
| §6 | Open maintainer items | §7 |
| §7 | Downstream work | Not this release |
| §8 | Guardrails | Global constraints |

## 3. Step order

Each step's task-level plan is written when the step starts, against the spec section it implements.

- **Step 0a — #83.** Done: `70a7c35`, `71aa424`, `62ad4ec`. First PR into the branch.
- **Step 0 — Baseline, oracle, truth set, test set** (spec §3 Step 0, plus):
  - **Test-set manifest.** Random within each stratum (source × assay), fixed seed, frozen
    before any measurement. Each stratum split into a measure half and a held-out half of at
    least 20 patients each: up to 160 patients. Recorded as a manifest outside the repository.
  - **A few Duplex BAMs per assay** in the tag survey and the Step 1 diff. Duplex BAMs otherwise
    first appear in 6b, too late to catch a Duplex-specific defect.
  - `MQ` and `MC` added to the tag-presence survey.
  - **Runtime baseline:** per-sample `run-all` wall time and peak RSS at production threads. This
    sizes every later pass; revision 1 stated none.
  - The **extraction** oracle is identical from `develop` or this branch: #83 changed no extractor
    code. The **feature** baseline for Step 3's measurements comes from the branch after 0a,
    so that #83's own WPS change is not charged to the first policy.
- **Step 0b — #77 pilot**, before any final per-sample pass. Run superseded Task 7 Steps 1–5 on a
  donor and patient subset, with `local_depth` recomputed from the fragment BEDs: the retained
  0.9.x outputs predate #79, so their `local_depth` is the broken column. Settles per-anchor vs
  per-bin support and the affordability stop while a per-bin column can still be added to Step 5.
- **Step 1 — `fragment.rs`, `MatePairer`, `FilterPolicy`** (spec §3 Step 1).
- **Step 2 — Projections onto `FragmentRecord`** (spec §3 Step 2).
- **Step 3 — Value-changing policies** (spec §3 Step 3), measured with the harness (§1): one
  extraction per test sample, each policy a projection of it.
  - Each policy against **the same base** (its main effect), and **all together** (the
    interaction is the difference). A commit sequence alone measures each policy on top of
    whatever came before, and a rejected policy would leave every later measurement on the wrong base.
  - Still one commit per policy, each carrying its evidence.
- **Step 4 — Migrate `region_mds.rs`, `uxm.rs`, `mfsd.rs`** (spec §3 Step 4). First a **neutral**
  migration commit per module that must reproduce today's values exactly, then **one commit per
  definition change**. Region MDS alone changes five things: duplicates, proper pair, MAPQ rule,
  coordinate source, both ends.
- **Step 5 — Outputs** (spec §3 Step 5), plus: the `--fragment-export` flag and the free-space
  check; the re-score command; the BED fifth column with four-column fallback; refused PON fails
  loudly; Nextflow tool-level mode passes `--emit-bed`, because it feeds `KREWLYZER_EXTRACT.out.bed`
  to every tool and an empty channel would run nothing and still succeed.
- **Step 6 — Release work** (spec §3 Step 6):
  - **6a. Versioning, in this order:** bump `__version__`, then `MIN_PON_VERSION`, then build and
    `stamp-pon`. Raising the floor first makes every new PON refused, because `pon/build.py`
    stamps the current `__version__`.
  - **6b. One per-sample pass** on the final code, with `--emit-bed`: the normals (All Unique and
    Duplex) and the test-set patients. The normals also with `--fragment-export`, as penumbra
    input, once space is confirmed.
  - **6c. #77 and #80.** Build **provisional** PONs from 6b's donor outputs; **re-score** the test
    patients against them. Measure with the corrected rules (§4), implement, re-score again.
  - **6d. Final PONs.** Rebuild all four; `validate-pon` with the donor invariant (superseded
    Task 11); premise assertions (Task 12); re-score the held-out patients against the final PONs.
    **Acceptance:** xs1 `PLAUSIBLE_Z_SCORES` findings near zero on the held-out half (superseded
    Task 13). xs2 is reported separately.
  - **6e. Penumbra map** from the normals' parquet, reported before freezing (spec §2g), or a
    recorded null result.
  - **6f. Cohort run**, 16,000+ clinical ACCESS samples, with the frozen assets. **Acceptance:**
    every completion marker present; `validate-cohort` clean; `validate-output` with no new
    error class; the xs1 `PLAUSIBLE_Z_SCORES` rate at or below 6d's.
  - **6g. Close out:** results on #77, #80, #83; a new issue for the xs1/xs2 asymmetry (superseded
    Task 14).

## 4. Corrections to the #77/#80 decision rules (superseded Tasks 7–8)

- **#80's gate measured the wrong direction.** A within-donor Spearman across anchors reflects
  capture efficiency and biology. Depth matters across samples at each anchor, which is where z is
  computed. New gate: for each anchor, the Spearman across donors of amplitude against sample
  depth, summarised as the median over anchors, compared against its noise floor (standard error
  ≈ 1/√(n−1): 0.22 at n = 21). Watch for spurious ratio correlation when the candidate divides by
  `local_depth`.
- **Batch: decide on the clinical stratum.** The PON donors are clinical-side. The IGO test
  patients are not, so their z-scores carry a batch difference that the final 16,000-sample run
  will not. `S*` and #80's form are decided on the clinical strata (the clinical research mirror);
  the IGO strata are a robustness check, where a failure
  may be batch rather than definition.
- **The xs2 asymmetry may be batch.** #77 found xs2 badly behaved even at well-covered positions
  (median |z| 6.5–12.7), but measured it on the IGO cohort against clinical donors. Batch is a third
  explanation, beside "21 donors too homogeneous" and genuine patient deviation, and the stratified
  test set tests it at no extra cost: if clinical xs2 patients behave like xs1, it is batch.
  Report it in 6c; it shapes the issue opened in 6g.
- **`S*` needs power and an honest held-out check.** At least 20 patients per stratum half; a
  minimum cell count per decade; bootstrap intervals on each decade's inflation. A p99 rises with
  sample size, so a small set biases `S*` low. `S*` ships only if the held-out half meets the same
  rule in every stratum.

## 5. Additions to the spec's test list (§4)

- **Done (#83):** `coordinate_convention_tests` (six), the Alu background test, `binning_tests` (two).
- **BED projection contract.** Synthetic BAM → `extract` → `wps`, checked against the known
  fragment coordinates: 0-based half-open `[start, end)`. Revision 1 cited the `2220/2001`
  integration pin here, but that test reads a hand-built BED and cannot fail on extractor changes.
- **Re-score parity.** Re-scoring an output against PON *P* gives the same z-scores as a fresh run
  with *P*.
- **Export parity.** `--fragment-export` through `extract` and `run-all` gives identical files; a
  run refuses to start when free space is below the estimate.
- **Fifth column.** Five-column and four-column BEDs both load; the four-column path logs its
  fallback.
- **Refused PON fails loudly**, with a test that a run given an old PON exits non-zero.
- **Nextflow tool-level mode** produces features with `--emit-bed`.
- **§2f test 3b — feature impact**, reported for each policy: per feature, median
  `|Δ| / between-donor SD` over cells that pass `S*` (or with a floored SD), the share of samples
  above 0.1 SD, and the Spearman of sample order before and after. **Informational, not a gate.**
  Restricting to supported cells matters: at uncovered cells the between-donor SD is about 1e-6,
  so the ratio would reproduce #77's own defect. Not used for mFSD (normals carry no variants) or
  UXM (no methylation libraries); those take simulated-truth evidence only.
- **§2f test 1 is not used for judgment calls.** For those, the simulated truth is the policy by
  construction — molecule length is "true" because the simulator defines it so. They rest on
  tests 2 and 3 and the maintainer's decision.
- **#77/#80:** scorability (Task 9), amplitude unification and depth-scaling invariance (Task 10),
  donor invariant (Task 11), premise assertion (Task 12), held-out `S*` (§4).

## 6. Interactions neither document states

| # | Interaction | Resolution | Step |
|---|---|---|---|
| I1 | WPS reads `[start, end)` from the BED under the #83 convention | Asserted by the synthetic-BAM contract test (§5) | 1 |
| I2 | WPS assigns fragments to the nucleosome or TF stream by length | Molecule length in both (§1); the same choice applies to FSD and WPS | 3, 5 |
| I3 | Unresolved-GC fragments leave GC-weighted features, WPS and `local_depth` included | #77 measured after Step 3 | 6c |
| I4 | **Revision 1 was wrong.** It said uncapping GC length bins moves WPS. WPS weights only 35–180 bp fragments; uncapping moves only fragments of 395 bp and longer | Point row 2's evidence at FSC `ultra_long` and FSD | 3 |
| I5 | Spec §5 says state the consequences before starting | Stated in §8 | — |
| I6 | One PON version bump for everything | Order fixed in 6a; `feat/pon-version-guard` is merged into `develop` | 6a |
| I8 | Every run starts from the BAM | Jobs read clinical BAMs in place, in parallel (§1) | 0, 6b, 6f |
| I9 | `WpsBackgroundConsumer` reads fragments through the same BED | Fixed with #83; covered by I1 | 0a, 1 |
| I10 | Definitions measured by fragments gained and lost only | Extended by test 3b (§5) | 3, 4 |
| I11 | Patient z-scores are written inside the per-sample pass (`wps_pon` rewrites `{sample}.WPS.parquet`) | Re-score command; provisional PONs in 6c | 5, 6c |
| I14 | Four PONs ship; an un-rebuilt one is refused, and today that drops z columns **silently** | All four rebuilt; refusal fails loudly | 5, 6d |

## 7. Decision log

**Open, from the spec (§6):** each stringency or judgment policy on its §2f evidence; the storage
allocation before any cohort-wide export; the penumbra map's thresholds, or a null result.

**Open, new:**

None beyond the spec's.

**Answered (18 Nov):** the healthy donors were sequenced on the clinical side; the second test
cohort is krewlyzer output from clinical ACCESS BAMs mirrored for research, so also clinical-side.
See §4 for what that changes.

**Decided (18 Nov), from Step 0: identical-coordinate fragments stay distinct rows** — no
collapse and no flag column. Closes the spec's §2d row "Identical-coordinate fragments". Measured on
the 228 test-set All Unique BAMs, comparing UMIs among fragments that share exact coordinates, against
a chance control built from UMIs at other coordinates in the same BAM:

| | XS1 (Marianas, read-name UMI) | XS2 (fgbio, `RX`) |
|---|---|---|
| UMI pairs within one mismatch, on-target (chance) | 6.1% (1.9%) | 1.0% (1.1%) |
| identical-UMI pairs, off-target | 0.2% | 3.3% |
| fragments a collapse would remove: exact / within one mismatch | 0.01% / 2.6% | 0.20% / 0.34% |
| share of those that are singletons (one mismatch) | 76% | 63% |

Duplex BAMs sit at the chance level in both assays (1.0–1.1%), so mFSD, which reads them, is
unaffected. The XS1 residual is real but cannot be separated from chance matches: roughly a third of a
one-mismatch collapse's on-target removals would be distinct molecules. The gain does not justify
re-doing the upstream collapse, and without a collapse a flag column has no consumer. Residual per BAM
is 0.9–6.0% (XS1) and 0.13–1.0% (XS2); the output documentation states it.

**Decided (18 Nov), from Step 0: the §2i budget is confirmed at the spec's starting values** —
wall time ≤ 1.5× baseline and peak RSS ≤ baseline — measured on four fixed benchmark samples, the
median-memory sample of each assay × sequencing site, always at a 64 GB limit: the reported peak
grows with the limit it is given (file cache), so a comparison at different limits measures the
limit. Memory is a guard against regression, not an optimisation target.

| benchmark sample | wall time | peak RSS at 64 GB |
|---|---|---|
| XS1 clinical | 35 min | 23.1 GB |
| XS1 IGO | 33 min | 27.1 GB |
| XS2 clinical | 49 min | 27.7 GB |
| XS2 IGO | 53 min | 35.0 GB |

At production's 32 GB first attempt the XS2 IGO median sat on the limit (31.9 GB) yet completed.
Extraction is 2% (XS1) and 5% (XS2) of `run-all` wall time, so mate pairing moves little of it; the
largest single stage is writing the TFBS region-entropy QC table, 17–22 min, through an unbuffered
writer. If that is fixed first, the wall-time baseline is re-measured before Step 2. Run-to-run,
141 of 160 tables are bit-identical and FSD/MDS differ at ~1e-16, so Step 1's diff against the
oracle uses a relative tolerance of 1e-12, not equality.

**Settled:** everything in §1.

## 8. Consequences, stated before starting (spec §5)

- Every feature value can move. The measurement harness quantifies by how much, per policy.
- kreview is paused. When it resumes it needs a pre-registration note for the re-baseline, and
  declarations for the new outputs: `fragments.parquet` (not read by its loader), the by-class
  tables, and the `wps_z_withheld` column if #77 adds it.
- The BED no longer appears in a finished run by default: release notes give `--emit-bed`.
- `fragments.parquet` carries read-derived bases at every fragment end, so it is as sensitive as
  the BAM: it stays on HPC, and no tool or CI job prints its rows.

## 9. Global constraints

The spec's §0 (as amended by §1) and §8 apply in full. Also:

- **All cluster work runs as SLURM jobs**; nothing runs on the login node beyond submitting jobs
  and reading their logs. Env rebuilds included.
- Rebuild the extension with `maturin develop --release` from the project root before `pytest`.
- Gates per change: `cargo test --lib`; `cargo clippy --tests` at CI's flags; `pytest tests/ -q`;
  `bash scripts/lint.sh`; `check_output_format.py` for output-writing changes; the claims registry.
- No patient identifiers, cluster paths or operator-specific settings in committed files.
- Measurements are reported as aggregates, never per-sample records.
- After 6b, a change to anything computed in the per-sample **extraction** needs another per-sample
  pass and must say so. Scoring changes do not: that is what re-score is for.
