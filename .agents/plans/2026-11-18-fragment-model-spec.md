# krewlyzer — one fragment model at the source: review and implementation plan

> Hand-off brief for Claude Code. Written against `msk-access/krewlyzer`
> **develop @ 0abb037** (18 Nov 2026). The extraction path is identical to
> main @ 504815b. Read-only code review; nothing was built or run. Read `AGENTS.md`,
> `.agents/rules/*.md` and `.agents/memory/modelling-philosophy.md` before touching
> anything. Every repo invariant applies.
>
> The defects below share one cause, so the plan fixes the cause. Nothing is deferred out
> of scope; work is ordered by dependency and ships as **one release**.
>
> Companion doc for the downstream side: `claude/kreview_followups_plan.md`.
>
> Confidence tags: **[Certain]** = read in the code at the cited line. **[Likely]** =
> strong inference. **[Guessing]** = unmeasured; the baseline step measures it.

---

## 0. Maintainer decisions (18 Nov 2026)

| Topic | Decision |
|---|---|
| Flow | **BAM in; feature tables and `fragments.parquet` out.** Every run starts from the BAM. |
| Input BAM | The cohort runs on the **All Unique** BAM. One export per sample. |
| `fragments.parquet` | The persisted fragment-level output. `{sample}.fragments.parquet`, declared in the output contract like every other table. **All fragments**; no subsampling, no off-target-only mode. |
| Row unit | **One row per fragment (molecule).** Not per reference position and not binned (§2e). |
| Schema coverage | **Capture everything the BAM pass can give per fragment** (§2a). |
| `bed.gz` | Temporary working file. The feature calculations keep reading it as today. Deleted at the end of a successful run unless `--emit-bed` is given. |
| Reuse / cache | None. The existence-only cache in `run-all` is removed (§1, finding 10). |
| Target labelling | Every fragment carries its relation to the bundled XS1/XS2 target intervals. |
| Value changes | Allowed when they move toward accuracy, as defined in §2f. |
| Measure before decide | **Required.** No value-changing policy is merged on argument alone (§2f). |
| Data location | Stays on HPC under the same access as the BAMs. Revisit only if the file, or a model trained on it, is to leave HPC. |
| Out of scope | **Methylation:** the cohort has no methylation libraries, so no UXM feature work; `uxm.rs` is touched only to share the filter policy. **ACCESS-HEME:** excluded for now. |

## 1. Root cause: krewlyzer has no single definition of a fragment

Four modules read the BAM, and each rebuilds "a fragment" from single read records with
its own filters. [Certain]

| Module | Duplicates | Proper pair | Which ends give motifs | Fragment coordinates |
|---|---|---|---|---|
| `extract_motif.rs` (BED, global motifs, GC observations) | Configurable | Configurable | Every passing read, so both ends when both mates pass | R1 only, from `pos` and TLEN sign |
| `region_mds.rs:519-565` | Always skipped | Not checked | **R1 only** — one end per fragment | R1 only, same TLEN logic, re-implemented |
| `uxm.rs:123-126` | Always skipped | Own rule | n/a | Own rule |
| `mfsd.rs:1398-1411` | Always skipped | Own rule | n/a | Own rule; only module that reads family size |

Every defect found in the review is a symptom of this:

1. **Mates are never joined.** A fragment is inferred from R1 plus TLEN. The mate's
   MAPQ, alignment end and soft clips are never seen. Global motifs and per-region MDS
   count different ends of the same molecules. [Certain]
2. **The fragment record is a formatted string**, `chrom, start, end, gc`
   (`extract_motif.rs:331`). Everything else available in the loop is dropped. [Certain]
3. **GC = 0.0 doubles as "unknown".** A failed fetch, an absent contig or an all-N insert
   writes 0.0 (`motif_utils.rs:119`), and `engine.rs:147` parses an unreadable value with
   `unwrap_or(0.0)`. Those fragments enter GC bin 0 of the correction model. Invariant #3.
   [Certain]
4. **Soft clips are ignored.** The motif is `seq[0..k]`; the boundary is the aligned
   position. On a clipped read the two disagree. [Certain unhandled; Guessing on prevalence]
5. **GC length bins are capped at 400 bp at accumulation** while the model declares 188
   bins to 999 bp (`validate/claims.py:339`, known divergence). [Certain]
6. **Whole-genome buffering** of fragment strings before writing; **linear scan** over all
   regions per fragment (`extract_motif.rs:87`) although `coitrees` is already a
   dependency. [Certain]
7. **Strand-family information is collapsed or lost.** Only `mfsd.rs:253-295` reads family
   size, and for Marianas read names it **sums** `posCount + negCount`. In an All Unique
   BAM that sum hides whether a molecule is duplex, simplex or a singleton. [Certain]
8. **"Off-target" means "does not overlap a target interval".** Only a targets BED ships;
   the WPS `BaitMask` is built from the same file. A fragment just outside a target is
   pooled with fragments megabases away when the "unbiased" off-target GC model is fit.
   [Certain for the definition; Guessing on near-target enrichment]
9. **The BED is an intermediate, not a product.** The tabix index is created in
   `core/sample_processor.py:733` and nothing in `rust/src` or `src/krewlyzer` queries it.
   kreview reads Parquet only. [Certain]
10. **The run-all cache is dead and unsafe.** `wrapper.py:608-621` skips extraction only
    if `{sample}.bed.gz`, `EndMotif.tsv`, `BreakPointMotif.tsv` and `MDS.tsv` exist. The
    names are hard-coded to `.tsv`; Nextflow has defaulted to Parquet since 0.9.0
    (`nextflow.config:72`), so the check never matches and every run starts from the BAM.
    When it does match (TSV output, same directory), it checks existence only, so changed
    filters or targets would silently reuse a stale BED. [Certain for the code; Likely
    for the consequence]
11. **Molecule length is taken as the reference span.** An insertion or deletion inside
    the reads makes the molecule longer or shorter than `end - start`; the CIGAR is not
    consulted. [Certain]

What already exists and needs no new work: per-fragment GC fraction, reference-derived
over the insert with Ns excluded. [Certain]

Not verified in this review:
- The bundled target BEDs are Git LFS objects and were not readable, so whether their
  intervals are probe footprints or target (exon) coordinates is unknown. Step 0 measures it.
- `core/output_provenance.py` exists and was not read. §2h builds on whatever it already
  records; read it first.

## 2. Target design

One module, `rust/src/fragment.rs`, owns the definition. Everything else is a projection.

```
BAM ──► MatePairer ──► FragmentRecord ──► FilterPolicy ──► projections
                                                         ├─ {sample}.fragments.parquet (persisted, all fragments)
                                                         ├─ {sample}.bed.gz            (temporary; feeds the feature engine)
                                                         ├─ GC observations            (off/on target)
                                                         ├─ global end + breakpoint motifs (optionally by region class)
                                                         └─ region MDS
```

### 2a. `FragmentRecord` and the `fragments.parquet` schema

**Inclusion rule.** A field is stored when it is **lost once the BAM pass ends**. A value
that can be recomputed later from coordinates plus the reference or a versioned annotation
is not stored, with two stated exceptions (group F). This keeps the file complete for
what cannot be recovered, without multiplying storage across every fragment or freezing
versioned assets into the data.

**A. Geometry**

| Field | Type | Notes |
|---|---|---|
| `chrom` | dictionary | `chr` prefix stripped |
| `start`, `end` | u32 | Outermost aligned bases of the two mates |
| `left_span`, `right_span` | u16 | Reference span of the left and right mate; mate overlap and inner gap follow from these and `end - start` |
| `ins_bp_left`, `del_bp_left`, `ins_bp_right`, `del_bp_right` | u8 | Inserted / deleted bases from each mate's CIGAR |
| `tlen_reported` | i32 | Aligner's TLEN, for audit against `end - start` |
| `mate_seen` | bool | False when the mate was absent or filtered; coordinates then come from TLEN |

**B. Orientation and pairing**

| Field | Type | Notes |
|---|---|---|
| `r1_is_left` | bool | Which mate is leftmost |
| `left_reverse`, `right_reverse` | bool | Strand of each mate; pair orientation (FR/RF/FF/RR) is derived |
| `proper_pair` | bool | The aligner's flag, recorded not required |
| `has_supplementary` | bool | Either mate carries an `SA` tag |

**C. Alignment quality**

| Field | Type | Notes |
|---|---|---|
| `mapq_left`, `mapq_right` | u8, nullable | Each mate's own MAPQ |
| `nm_left`, `nm_right` | u8, nullable | Edit distance (`NM`) |
| `as_left`, `as_right` | i16, nullable | Alignment score (`AS`) |
| `mean_bq_left`, `mean_bq_right` | u8 | Mean base quality of each mate |

**D. Fragment ends**

| Field | Type | Notes |
|---|---|---|
| `clip_end_left`, `clip_end_right` | u8 | Soft clip at each **fragment end** (the 5' clip of each mate) |
| `clip_inner_left`, `clip_inner_right` | u8 | Soft clip at each mate's 3' end (read-through / adapter side) |
| `hard_clipped` | bool | Any hard clip on either mate |
| `read_kmer_left`, `read_kmer_right` | u16-coded 8-mer, nullable | **Read-derived** bases from the first aligned base inward at each fragment end; null if any base is N. The 4-mer motif is its prefix |
| `end_mismatch_left`, `end_mismatch_right` | u8, nullable | Mismatches between that 8-mer and the reference |

**E. Molecule (consensus) information**

| Field | Type | Notes |
|---|---|---|
| `family_pos`, `family_neg` | u16, nullable | Raw reads supporting each strand: Marianas `posCount` / `negCount`, or fgbio `aD` / `bD` |
| `consensus_depth_max`, `consensus_depth_min` | u16, nullable | fgbio `cD` / `cM` |
| `consensus_error` | f32, nullable | fgbio `cE` |
| `name_start_left`, `name_start_right` | u32, nullable | Start coordinates encoded in Marianas read names; an independent check on the paired coordinates |

Duplex is `family_pos > 0 && family_neg > 0`; simplex is one strand with ≥2 reads;
singleton is a total of 1. The counts are stored; the class is derived. File metadata
records `bam_type = all_unique` and **`family_source`** (`marianas_name` or `fgbio_tags`),
because the two encodings may not count the same thing (Step 0, comparability check).

**F. Derived, stored by exception**

| Field | Type | Why stored although derivable |
|---|---|---|
| `gc` | f32, nullable | Already computed and used by the run's GC model; null when unresolved, never 0.0 |
| `target_overlap_bp` | u16, nullable | Records the target relation the run actually used |
| `target_distance` | u32, nullable | Distance to the nearest target interval; 0 when overlapping |

`on_target` is `target_overlap_bp > 0`; "near-target" is a consumer-chosen threshold on
`target_distance`. The targets file name and checksum go in file metadata.

**Deliberately not stored** (recompute from coordinates + reference or annotation):
insert sequence; reference k-mers at and beyond the boundaries (breakpoint context);
CpG or dinucleotide content; GC correction weight and any GC residual; overlap with TFBS,
ATAC, gene, exon, bin or arm annotations; region class (§2g); mappability; any
PON-relative value. **Not obtainable** from this assay: jagged-end length, methylation.
**Considered and excluded:** the UMI string (no analysis use identified; poorly
compressible). The UMI is still **read once in Step 0** for the residual-duplicate check.

**Tag-dependent columns.** Groups C and E depend on tags and read-name formats that may
differ between XS1 (Marianas) and XS2 (fgbio). Step 0 reports presence per assay. A
column that would be null for every fragment of every sample is **dropped, not shipped**
(invariant #1); a column present for one assay only is kept and declared in the contract
with its reason.

### 2b. `MatePairer`

- Pair by read name inside each chunk with a hash map. Fetch each chunk with
  `max_len` padding on both sides; a fragment is **owned by the chunk containing its
  start**, so every fragment is emitted exactly once.
- Evict a pending read once the cursor passes `pos + max_len`; an evicted read becomes a
  record with `mate_seen = false`, flagged, not dropped silently.
- Pending-map size is bounded by depth × `max_len`, not by the genome.
- Consensus BAMs are not proper-paired by flag (`docs/features/core/extract.md`), so
  pairing must not depend on the proper-pair bit.

### 2c. `FilterPolicy` — one struct, used by all four BAM readers

- Unmapped, secondary, supplementary, QC-fail: always excluded.
- Duplicates and proper-pair: one setting each, honoured everywhere.
- MAPQ: one rule for both mates (default in §2d).
- Length: one `min_len` / `max_len`.
- `uxm.rs` and `mfsd.rs` stay locus-driven but take their filters, fragment-coordinate
  function and family parsing from the shared modules. No copies.

### 2d. Policies that change values

| Topic | Policy | Kind | Decided by |
|---|---|---|---|
| Unresolved GC | Excluded from GC observations and GC-weighted features; counted in metadata; BED column 4 written as `NA`; `engine.rs` parses to `Option` | Correctness fix | §2f tests 1–2 |
| GC length bin | `LengthBin::from_len` is the only binning function; the cap pattern is deleted | Correctness fix | §2f tests 1–2; the claims test forces the divergence entry out |
| Region MDS ends | Both ends, same as global motifs | Correctness fix | §2f tests 1–2 |
| Fragment coordinates | From both mates, not R1 plus TLEN | Correctness fix | §2f tests 1–2 |
| Fragment length for size features | **Candidate:** molecule length (reference span corrected for indels) instead of reference span | Judgment call | §2f measurement; not merged until reported |
| Mate MAPQ | **Candidate default:** both mates must pass | Stringency choice | §2f measurement; not merged until reported |
| Soft-clipped end | **Candidate default:** excluded from aggregate motif counts, tallied separately | Judgment call | §2f measurement; not merged until reported |
| Identical-coordinate fragments | Kept as distinct rows unless Step 0's UMI check shows residual duplicates; then the policy is decided on that evidence | Open | Step 0 measurement |
| Off-target GC model | Unchanged unless the Step 0 GC-by-distance curve shows near-target enrichment | Open | Step 0 measurement |

For the candidate policies the per-fragment fields are in `fragments.parquet` either way,
so the aggregate policy can be revisited without re-extracting.

### 2e. Outputs

#### `{sample}.fragments.parquet`

Every fragment that passes `FilterPolicy`, with the §2a columns. Written by `extract` and
`run-all` under `--fragment-export`; the cohort configuration turns it on. It is an
**output only**: the feature engine does not read it.

**Row unit: one row per fragment.** The file is not per reference position and not
binned. Each row carries its own `start` and `end`, so the only interval is the fragment
itself.

| Case | Handling |
|---|---|
| Two fragments covering the same position | Two rows. Nothing is merged, clipped or split |
| Two fragments with identical coordinates | Two rows, when they are distinct records in the BAM. No coordinate-based deduplication is applied here (see the open policy in §2d) |
| The two mates of one fragment overlapping each other | One row. The overlap is recoverable from `left_span`, `right_span` and `end - start`, and is never counted twice |
| A fragment crossing a 10 Mb processing chunk | One row, owned by the chunk containing its `start` (§2b) |
| Position-level or bin-level quantities (coverage, WPS, bin counts) | Not stored. They are computed from fragment rows; the existing feature tables remain the interval-level products |

**Row order (deterministic).** Sorted by `chrom` in reference-header order, then `start`,
then `end`, then `left_reverse`, then `r1_is_left`, then `read_kmer_left`,
`read_kmer_right`, then `family_pos`, `family_neg`. Two runs on the same BAM with the same
settings must produce identical row order, whatever the thread count. If two records
still tie on every key, their order must be made stable by the writer (for example by
input order within the owning chunk) and the determinism test must cover it.

**Row groups (physical storage only).** A row group never spans two chromosomes. Within a
chromosome, groups are cut by **row count**, not by genomic window: on a capture panel
most fragments fall on a small part of the genome, so fixed windows would give a few
enormous groups and many near-empty ones. The row cap is a named constant (starting
value 1,000,000) set from Step 0's fragment counts and registered in
`validate/claims.py`. Dictionary-encode `chrom`; write min/max statistics for `start`,
`end` and `target_distance`.

**Reading a region.** Rows are filed under `start`, so a reader that wants every fragment
overlapping `[a, b)` must scan from `a - max_len` and keep rows with `end > a` and
`start < b`. `max_len` is in file metadata. State this in the output documentation.

**File metadata.** Genome build, assay, `bam_type`, `family_source`, targets file and
checksum, filter settings including `max_len`, k-mer length, row-group cap, sort key,
krewlyzer version.

#### `{sample}.bed.gz`

Same four columns and sort order as today, a projection of `FragmentRecord`, streamed per
chromosome. The feature engine reads it exactly as now.

- Deleted at the end of the run unless `--emit-bed` is given. Deletion happens only after
  every enabled feature has finished, including any retry, and `.metadata.parquet` has
  been written. On any failure the BED is kept and its path logged.
- The tabix index is built only with `--emit-bed`.
- Single tools still take `-i sample.bed.gz`; produce one with `extract --emit-bed`.

#### Removed

The existence-only cache in `wrapper.py:608-630`. Every run starts from the BAM, and the
progress line no longer reports "cached".

#### Other changes that follow

- Nextflow `extract` and `runall` modules: `*.bed.gz` and its index become
  `optional: true`; `*.fragments.parquet` is emitted; resource labels are revised from
  the Step 0 and Step 2 measurements (§2i).
- `pon/from_outputs.py` documents its completion order starting at `bed.gz`; update the
  note. The completion marker itself (`.MDS.gene`) is unaffected.
- `metadata`: new counts — `n_mate_unseen`, `n_gc_unresolved`, `n_clipped_end_left/right`,
  `n_mapq_mate_fail`, `n_on_target`, `n_off_target`, `n_duplex`, `n_simplex`,
  `n_singleton`, `n_with_indel`, `n_identical_coordinates`, `bam_type`, `family_source`,
  pairing settings, `bed_emitted`, and the GC-bias summary of §2h.

### 2f. Acceptance rule — measure before decide (required)

"More accurate" has no ground truth in real plasma, so it is defined operationally. Every
value-changing commit must carry the evidence below in its PR description. A commit
without it is not merged.

1. **Simulated truth.** A generator in `tests/invariants/` writes synthetic BAMs from
   fragments with known coordinates, ends, soft clips, indels, mate qualities,
   strand-family counts and target relation. The change must move the output **closer to
   the simulated truth**, reported as error before and after.
2. **Internal consistency on real data** (aggregates only), before and after:
   - `tlen_reported` agrees with `end - start`;
   - paired coordinates agree with `name_start_left/right` where read names carry them;
   - left-end and right-end motif distributions are symmetric;
   - global motif total equals the sum over regions plus the remainder;
   - fragment counts are conserved across chunk boundaries.
3. **Cost.** Fragments gained or lost by the change, overall and by on/off-target and by
   molecule class.
4. **Never downstream performance.** kreview AUC or sensitivity is not evidence for an
   extraction policy. Choosing a policy because the endpoint improved tunes extraction on
   the endpoint and invalidates the pre-registered numbers.

Decision rule: a correctness fix merges when tests 1–2 improve or hold. A stringency or
judgment policy merges only when test 1 improves **and** the cost in test 3 is reported
and accepted by the maintainer; otherwise the current behaviour stays and the
per-fragment fields carry the information.

### 2g. Region-class stratification and the penumbra map

Two separate pieces. The mechanism is generic; the penumbra map is one asset that uses it.

**Mechanism — `--region-classes <BED>`.** A BED whose fourth column is a class label.
When given, the global end-motif / MDS projection and the fragment-size distribution are
emitted **additionally** per class, as long-format tables with a `region_class` column
(`{sample}.EndMotif.byclass`, `{sample}.FSD.byclass`). A fragment end is assigned by the
position of that end; a fragment for FSD by its midpoint. Fragments outside every class
go to an explicit `unclassified` row, never dropped. Existing tables are unchanged. The
class BED name and checksum are recorded in metadata.

These are **within-sample contrasts** of motif and size, not coverage, so they do not
measure copy number.

**Asset — the penumbra map.** A static, versioned class BED per assay, built like the GC
reference:

- Input: `fragments.parquet` from the **healthy-normal set only** (the same normals used
  for the PON rebuild in Step 6), never from patient samples.
- Statistic: per-region unique-molecule yield (duplex + simplex + singleton counts from
  `family_pos` / `family_neg`), corrected for GC and for probe performance, aggregated
  across normals; regions are classed by residual depletion.
- Region definition, thresholds and the number of classes are **not fixed here**. They
  are set from the normal-pool distribution and reported before the asset is frozen
  (measure before decide). Thresholds from low-pass WGS studies are not to be lifted.
- A **per-sample yield or penumbra score is not emitted as a feature.** In a patient
  sample, per-region yield tracks amplifications and deletions, and copy number stays out
  of the model.
- A null result is acceptable and is reported: if no residual depletion survives GC and
  probe correction, no map ships and the mechanism remains available for other class BEDs.

### 2h. Provenance and GC-bias summary in the outputs

- **Provenance in every output table.** Read `core/output_provenance.py` first and extend
  what it already writes so that every feature table and `metadata` carry: assay, genome
  build, `bam_type`, `family_source`, krewlyzer version, fragment-model version, targets
  checksum, PON version, and region-class BED checksum when used. One writer; no table
  hand-rolls its own.
- **Per-sample GC-bias summary** in `metadata`, computed from the correction factors the
  run already fits, separately for the off-target and on-target models: the 5th, 50th and
  95th percentile of the fragment-weighted correction factor, and the ratio of the mean
  factor in a high-GC band to a low-GC band for mononucleosome-length fragments. Band
  edges are named constants in `validate/claims.py`. These are QC and nuisance variables
  for downstream audits, not model features. Each needs a `must_vary` setting.

### 2i. Runtime and memory budget

Mate pairing adds a hash map and a padded fetch; the Parquet writer adds a second output.

- Step 0 records wall time and peak RSS for `run-all` on a real-size All Unique sample at
  the production thread count, and the share spent in extraction.
- **Budget (starting values, to be confirmed by the maintainer after Step 0):** wall
  time ≤ 1.5× baseline and peak RSS ≤ baseline, since whole-genome buffering is removed.
  The values live in one place and are asserted by a local benchmark script.
- Steps 2 and 5 each report the same two numbers. Exceeding the budget blocks the merge
  until the cause is reported.
- Nextflow resource labels for `extract` and `runall` are set from the measured figures.

## 3. Work order (dependency order; one release, e.g. 0.10.0)

Each step is a PR into a long-lived `feature/fragment-model` branch off `develop`. The
branch merges once, so the cohort is rerun once and PONs are rebuilt once.

**Step 0 — Baseline, oracle and truth set.**
- Freeze the current extractor's outputs on the fixtures and, with permission, the
  real-data corpus, as an equivalence oracle (precedent:
  `tests/unit/test_rust_python_equivalence.py`).
- Build the §2f simulated-truth generator.
- Measure on All Unique BAMs, aggregates only: fragments per sample and per chromosome;
  on/off-target split; duplex / simplex / singleton shares; `gc == 0.0` rate; soft-clip
  rate at each end and inner; indel rate; mate-MAPQ failure rate; presence of `NM`, `AS`,
  `SA`, `aD/bD/cD/cM/cE` and Marianas read names, **for XS1 and XS2 separately**; bytes
  per fragment for the full §2a schema and per column group; BED size per sample; wall
  time, extraction share and peak RSS (§2i).
- **Residual-duplicate check.** Among fragments with identical coordinates, compare UMIs
  (read name or `RX`): report the share with identical UMIs, the share within one
  mismatch, and the share clearly distinct, on-target and off-target. Aggregates only;
  the UMI is not written anywhere.
- **Family comparability, XS1 vs XS2.** Document from the upstream tools what Marianas
  `posCount` / `negCount` and fgbio `aD` / `bD` each count. Compare the distributions of
  per-strand counts and the duplex / simplex / singleton shares between assays at matched
  depth. If the definitions differ, say so in the output documentation and keep
  `family_source` as a required stratifier; do not harmonise silently.
- Measure the bundled target BEDs: interval count and length distribution for XS1 and
  XS2, and the GC curve and fragment density by `target_distance` band.
- From the per-chromosome counts, confirm or revise the row-group cap.

**Step 1 — `fragment.rs` + `MatePairer` + `FilterPolicy`.** No output changes yet: build
records, project the old BED from them, and diff against the oracle. Expected: identical
for well-formed FR pairs without clips; every difference classified and counted.

**Step 2 — Move all projections onto `FragmentRecord`.** BED, GC observations, global
motifs. Interval lookup replaces the linear scan. Per-chromosome streaming replaces
whole-genome buffering. Report runtime and memory against §2i.

**Step 3 — Value-changing policies (§2d), one commit each,** each with the §2f evidence
and a CHANGELOG "values are not comparable" note.

**Step 4 — Migrate `region_mds.rs`, `uxm.rs`, `mfsd.rs`** to the shared policy,
coordinate function and `molecule.rs`. Delete the re-implementations. `mfsd`'s duplex
weight keeps its current formula but reads `family_pos + family_neg` from the shared
parser.

**Step 5 — Outputs (§2e, §2g mechanism, §2h).** `fragments.parquet` writer with the
deterministic sort and row-count row groups; `--emit-bed` and end-of-run BED deletion;
removal of the cache; `--region-classes` and the two by-class tables; provenance and
GC-bias summary; Nextflow output declarations; contract entries with a `must_vary`
setting per column; anti-degeneracy test; `validate-output` coverage; docs, including
the region-reading rule. Report runtime and memory against §2i.

**Step 6 — Release work.** PON format/version bump so an old PON cannot score new
features (see `feat/pon-version-guard`); run the healthy-normal set; rebuild XS1/XS2 All
Unique PONs; **build the penumbra map from the same normals (§2g) and report its
distribution before freezing it**; rerun the cohort with the frozen assets; hand over to
kreview (§5).

## 4. Tests that must exist before merge

- The §2f truth-set generator, with cases for: R1 rightmost; mate in the next chunk;
  mate missing; soft clip at a fragment end and at an inner end; hard clip; insertion and
  deletion inside a mate; overlapping and non-overlapping mates; mate MAPQ below
  threshold; contig absent from FASTA; all-N insert; N inside an end 8-mer; mismatch
  inside an end 8-mer; fragment ≥400 bp; supplementary alignment present; fgbio-tagged
  record; Marianas-named duplex, simplex and singleton records; fragment overlapping a
  target, abutting one, and far from any; two fragments overlapping each other; two
  distinct fragments with identical coordinates.
- **Each §2a column is asserted by at least one case**, by value.
- Fragment emitted exactly once across chunk boundaries (extends
  `tests/invariants/test_conservation.py`).
- Row count of `fragments.parquet` equals the line count of the BED for the same run.
- **Determinism:** the same BAM run with 1 thread and with many threads gives identical
  row order and identical row-group boundaries, including fully tied records.
- **Row groups:** no group spans two chromosomes; no group exceeds the cap; `start` is
  non-decreasing within a chromosome.
- **Region read:** a reference reader applying the §2e rule returns exactly the fragments
  overlapping a region, including one that starts before it.
- BED lifecycle: absent after a successful run without `--emit-bed`; present with it;
  present after a failed feature.
- Global motif total equals the sum of per-fragment ends counted under the active policy;
  the 4-mer motif equals the prefix of `read_kmer_*`.
- **Region classes:** by-class motif counts and by-class FSD counts each sum to the
  global table, including the `unclassified` row; without `--region-classes` no by-class
  table is written.
- **Provenance:** every output table of a run carries the same provenance block.
- Single-tool vs `run-all` parity for `--fragment-export`, `--emit-bed` and
  `--region-classes` (invariant #6).
- Anti-degeneracy: two different BAMs give two different `fragments.parquet` tables, and
  no shipped column is constant across them without a written reason.

## 5. Consequences outside this repo (state these before starting)

- **Every feature value can move**, because fragment membership and motif ends change.
  [Likely] Small for consensus BAMs. Step 1's diff measures it.
- **kreview's frozen numbers are re-baselined.** Under `ANALYSIS_PLAN.md`'s between-runs
  rule this is a behaviour change needing a pre-registration note before the next
  production run. Batch it with the patient-grouped CV change so kreview reruns once.
- **kreview must know the new suffixes.** The output contract here requires a matching
  declaration downstream for any new table. `fragments.parquet` and the two by-class
  tables are declared in kreview as known outputs; `fragments.parquet` is marked
  **not read** by the tabular loader. See `claude/kreview_followups_plan.md`.
- **Storage is a commitment.** Every fragment, roughly 35 columns, one BAM type per
  sample. [Guessing] 5–25 TB for 16,000 samples; deleting the BED offsets part of it.
  Step 0's per-column-group figures replace the guess, and the HPC allocation is
  confirmed against them before Step 6.
- **The BED no longer appears in a finished run by default.** A visible change for
  PyPI and pipeline users; it goes in the release notes with the `--emit-bed` instruction.

## 6. Open items for the maintainer

1. Accept or reject each stringency/judgment policy once its §2f evidence is reported.
2. Confirm the HPC storage allocation once Step 0 reports bytes per fragment. If it is
   over budget, group C (alignment quality) is the first candidate to trim.
3. Confirm the runtime and memory budget (§2i) once Step 0 reports the baseline.
4. Approve the penumbra map's region definition and thresholds once the normal-pool
   distribution is reported (§2g), or accept a null result.

## 7. Downstream work this unblocks (separate repos or assets)

- D5 dilution mixing at fragment level, valid for motif features as well as coordinate
  ones (design in the kreview follow-ups plan).
- Fragment-encoder (LEAF-style) experiment, reading `fragments.parquet` with its own
  loader and reconstructing insert sequence from coordinates.

## 8. Guardrails (from AGENTS.md — re-read before each PR)

- `pytest tests/ -q`, `bash scripts/lint.sh`; `cargo test` (macOS) and `cargo clippy` on
  Rust changes; `maturin develop --release` before any Python test.
- `python scripts/check_output_format.py` on output-writing changes;
  `pytest tests/unit/test_claims_registry.py` on constants.
- Any constant a doc quotes goes in `validate/claims.py`.
- Placeholders `P-0000000`, `C-000000-L000` only. Never commit data artifacts.
- Conventional commits; docs and CHANGELOG in the same commit as the change.
- Ask before: pushing to shared branches, any cluster submission, releases.
