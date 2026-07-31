# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added
- **`scripts/build_gene_bed.py`, and gene BED assets rebuilt from GENCODE.**
  The bundled gene BEDs could not answer "which exon is first". The panel
  assets (`xs1`/`xs2`) had five columns and **no strand at all**; the WGS asset
  numbered exons in *coordinate* order, so MTOR — minus strand — carried
  `exon_num 0` at its lowest coordinate while its real exon 1 is at its
  highest. Anything deriving E1 from those files got the wrong end of every
  minus-strand gene.

  GENCODE's `exon_number` is transcription-ordered and `MANE_Select` names the
  agreed transcript, so the question is resolved once at build time and written
  into the asset as an explicit `is_e1` column. Canonical-transcript policy is
  `--transcript-overrides` → MANE Select → GENCODE basic protein-coding →
  Ensembl canonical → longest CDS; the build **fails** if the GTF carries no
  MANE tags at all, which is what silently happens with Ensembl's frozen
  GRCh37 release 87.

  New columns are additive — the first five are unchanged and the existing
  parser produces identical genes and coordinates: `transcript_id`,
  `exon_number`, `strand`, `is_e1`, `is_alt_e1`, `is_first_captured`.

  Three panel symbols (`H3F3A`, `HIST1H3B`, `PAK7`) were renamed by HGNC and
  match nothing in a current GTF. They are carried by an explicit alias table
  and asserted, because the failure mode is a gene silently losing its
  annotation rather than an error.

  **"First exon" is not one thing.** Genes carry a median of 13 distinct
  annotated first exons, because alternative promoters are the norm, so a
  single boolean cannot describe what a tile is. Three columns:

  | tile overlaps… | column | xs1 / 128 | xs2 / 146 |
  |---|---|---:|---:|
  | canonical transcript's exon 1 | `is_e1` | 25 | 33 |
  | another basic protein-coding transcript's exon 1 | `is_alt_e1` | +15 | +15 |
  | most 5′ captured tile (always exists) | `is_first_captured` | 128 | 146 |

  MSK-ACCESS tiles coding hotspot exons, so AKT1's canonical exon 1 sits 15 kb
  past the panel's most 5′ tile — but many genes are captured at an
  *alternative* promoter instead, which a MANE-only view misses entirely. 40
  `xs1` genes have a tile on some basic protein-coding first exon; 88 have
  none. `is_first_captured` exists for every gene but is frequently an internal
  exon, which is not a promoter proxy and must not be read as E1.

  Which transcript is canonical is configurable per gene via
  `--transcript-overrides` (a `gene<TAB>transcript_id` TSV), because a panel
  designed around specific clinical transcripts should not have MANE imposed on
  it. A transcript that is absent from the GTF, or belongs to a different gene,
  is a **hard error** — a silent fall back to MANE would produce an asset that
  disagrees with the file the operator wrote.

  This also bounds the existing `filter_fsc_to_e1` defect: of the 25 `xs1`
  genes with a canonical E1 tile, its lowest-start pick is correct for 18 and
  wrong for 7 (6 minus-strand). The larger issue is the 88 with no captured
  first exon at all, for which it emits a row labelled E1 that is an arbitrary
  internal exon.

- **`krewlyzer validate-output`** — checks a finished output directory against
  the contract its consumers rely on. Three layers, in increasing order of what
  a schema alone can catch: every consumed table present and shaped correctly;
  domain invariants (frequencies sum to 1, chromosomes `chr`-prefixed, the six
  FSC channels partition `total`, FSD carries only size-bin columns); and
  **anti-degeneracy** — a metric that is identical across every sample is an
  error, not a pass.

  That last layer is the point. Run against a real 0.8.3 cohort the schema
  passes completely, while `WPS_background.nrl_bp` ≡ 150.0,
  `nrl_deviation_bp` ≡ 40.0, `periodicity_score` ≡ 0.3333 and
  `adjusted_score` ≡ 0.0 in every sample: four of the five columns that table
  contributes carry no information. A schema-only gate certifies that directory
  as good.

  Exit codes are `0` satisfied, `1` contract violation, `2` structural (missing
  directory, unreadable Parquet) so a workflow can retry on 2 and escalate on 1.
  Cross-sample degeneracy below two samples is reported SKIP, never PASS.
  Declaring a column legitimately constant requires a written
  `constant_reason`, so silencing a finding costs a justification.

  `--json-report` emits stable finding ids for trend tracking.
  `scripts/validate_output.py` runs it from a checkout.

- **`krewlyzer validate-cohort`** plus the `KREWLYZER_VALIDATE_SAMPLE` and
  `KREWLYZER_VALIDATE_COHORT` Nextflow modules — the gate split into a scatter
  and a gather so it scales to a real cohort.

  A sample directory is ~1.5 GB (`WPS.parquet` alone is ~120 MB), so re-reading
  tens of thousands of them is not viable. `validate-output --fingerprint-out`
  reduces each sample to a ~20 KB fingerprint — a hash and two counts per
  column — and `validate-cohort` compares those. Reading is projected to the
  declared columns and bounded by `TableRule.scan_rows`, which roughly halves
  per-sample time and, more importantly, makes memory independent of cohort
  size.

  The split is not just an optimisation: degeneracy is inherently cross-sample.
  On a real cohort every sample passes `validate-output` individually while
  `validate-cohort` fails, because no single sample can distinguish "this
  metric is a constant" from "this is its value here".

- **`run-all` now writes `{sample}.validation.json` and
  `{sample}.fingerprint.json`** on Parquet runs, and `--strict-validation`
  makes a contract violation fail the run.

  Emitting is on by default because the fingerprint is a cheap byproduct here —
  the tables are already written — and it is what makes `validate-cohort`
  affordable later; leaving it opt-in would in practice disable the only check
  that catches a constant metric. Failing is opt-in because a contract rule
  that turns out too strict should not take down a cohort. Skipped entirely for
  `tsv`-only runs, since the contract describes the Parquet surface downstream
  reads. A checker that throws is caught and logged: it must never lose a
  completed run.

- **`WPS_background.nrl_at_band_limit`** — marks a right-censored NRL estimate,
  where the spectral peak sat on the edge of the 140-250bp search band so
  `nrl_bp` is that bound rather than a measurement.

  On real plasma this is **21% of Alu groups for XS1 and 43% for XS2**, and
  those groups are indistinguishable from interior ones by `periodicity_score`
  or fragment support — so a consumer reading `nrl_bp = 250.0` had no way to
  tell "the repeat length is 250bp" from "no nucleosomal peak was found". The
  same *present, plausible, wrong* failure mode the NRL fix itself addressed,
  one level down.

  Not triggered by a long period: 400bp and 2000bp synthetic signals both
  resolve to 225-235bp inside the band. Only the *absence* of periodicity pins
  the edge, so the flag means "no peak found", not "period too long".

  Additive and optional, so existing consumers are unaffected; the contract
  gate accepts directories written before it existed.

- **Every feature tool now accepts `--output-format` and `--compress`.** Eight
  of eleven had neither, so `krewlyzer fsc` could only write TSV while
  `run-all --output-format parquet` wrote Parquet -- and Parquet is all the
  downstream consumer reads. Anyone running a single tool got output nothing
  could load, with no error.

  The computation was never the problem: standalone `wps` and `run-all` produce
  byte-identical `WPS_background` on the same input, because both route through
  the same single-pass `run_features()`. Only the CLI layer hardcoded the
  format defaults, and the underlying processors already accepted both
  parameters.

  `scripts/check_output_format.py` could not catch this -- it verifies that
  internal call sites *forward* the parameters, not that the CLI ever lets
  anyone set them. A new test asserts every feature command exposes both.

  The eleven per-tool Nextflow modules now pass `params.output_format` and
  `params.compress_tsv` through as well; previously only `runall` did.

### Fixed
- **The Nextflow pipeline discarded its own validation artifacts, and never ran
  the cohort check.** Three compounding gaps, so the scatter/gather validation
  added earlier was inert end to end:

  1. `run-all` writes `{sample}.validation.json` and `{sample}.fingerprint.json`
     on Parquet runs, but the `runall` module declared neither as an output —
     both were produced and left to die in the work directory.
  2. `KREWLYZER_VALIDATE_COHORT` existed as a module and was **included by
     nothing and called by nothing**. The gather half of the check had no
     inputs and never ran.
  3. `--strict-validation` was not exposed as a pipeline parameter at all.

  The workflow now collects the fingerprints and runs the cohort step, and
  emits `cohort_report`. Tool-level mode (`--use_runall false`) produces no
  fingerprints, so the cohort step is skipped rather than run on a partial set,
  which would report degeneracy that is an artefact of the missing samples.

  The pipeline also **warns at start** when `output_format` is tsv-only: kreview
  reads Parquet only and a tsv-only run skips validation, so the default
  produces a cohort that is invisible downstream — silently, because every
  reader there swallows exceptions and yields an empty feature dict. The
  default is unchanged; the warning and the documentation are new.

  Three parameters were also read by modules but **declared nowhere**:
  `validate_min_samples`, `silent` and `assay`. Nextflow returns `null` for an
  undeclared parameter rather than failing, so the `?: default` swallowed it —
  the parameters were invisible to `--help` and the docs, and would become hard
  errors the moment `nextflow.enable.strict` is enabled. All three are now
  declared and documented, and a test cross-checks every `params.*` read in a
  `.nf` file against the config.

  `docs/nextflow/outputs.md` gains the validation artifacts and a table of the
  six output families whose values change in 0.9.0; `parameters.md` gains
  `--strict_validation`, `--validate_min_samples`, `--gc_correct`,
  `--queue_size` and the Parquet caveat.

- **`FSC.regions.e1only` selected by coordinate, and emitted a row for every
  gene whether or not it had one.** `filter_fsc_to_e1` sorted by `start` and
  took the first row per gene. It had no strand to work with — the temp BED
  handed to the Rust aggregator was five columns, so the asset's strand and E1
  flags were discarded at that boundary regardless of what the asset contained.

  On a real xs2 sample, **24 of 146 rows were the canonical exon 1**; the rest
  were internal exons labelled E1.

  `FSC.regions` now carries `strand`, `is_e1` and `is_alt_e1` end to end
  (Python `Region` → 8-column temp BED → Rust `GeneRegion` → output), and
  selection is on the flags: a region qualifies if it overlaps the canonical
  transcript's exon 1 **or** another basic protein-coding transcript's first
  exon. Both are real transcription starts, and requiring the canonical one
  discards most of a panel — 25 of 128 xs1 genes against 40.

  **Genes with neither flag are now omitted rather than back-filled.** The
  table drops from 146 rows to 48 on that sample, and every remaining row is
  on an annotated first exon.

  A legacy or custom gene BED still works: the flags are absent, and selection
  falls back to lowest-start with a warning that it is not strand-aware. `.`
  rather than `0` marks an unknown flag, so "the source could not say" stays
  distinguishable from "no".

  **Output-contract impact.** `FSC.regions` gains three columns;
  `FSC.regions.e1only` changes row count and membership. Neither is read by
  kreview today (`e1only` is in `NOT_CONSUMED`), so nothing downstream breaks,
  but any local analysis keyed on the old 146-row shape will.

- **`MDS.gene` row order was non-deterministic.** The writer iterated a
  `HashMap`, and Rust randomises hash iteration per process, so two runs on
  identical input produced byte-different files. A comment above the loop
  claimed "stable, reproducible output ordering" — it sorted regions *within* a
  gene, not the genes themselves. `fsc.rs` already sorted its genes; this did
  not. Verified on a real sample: identical SHA-256 across two runs after the
  fix.

- **`region-mds` E1 was never strand-aware on panel data, and the new assets
  would have silently disabled it for WGS too.** Two compounding problems.

  The strand fix in `identify_e1_regions()` reads `RegionInfo.strand`, but
  `parse_gene_bed()` hard-coded `'+'` for the panel format — the 5-column panel
  BEDs have no strand column to read. So `mds_e1` reported the **last** exon
  for every minus-strand `xs1`/`xs2` gene, which the docs claimed was fixed.

  Worse, format detection keys on **column count** (5 → panel, 8 → WGS,
  anything else → warn and treat as panel). Regenerating the assets took them
  to 11 columns, so the WGS BED would have fallen into the panel branch and
  lost its strand as well — turning a panel-only defect into a universal one.
  The existing E1 tests build `RegionInfo` by hand and never call
  `parse_gene_bed`, so none of them would have failed.

  `region-mds` now recognises the annotated format, reads strand from it, and
  consumes the precomputed `is_e1` instead of re-deriving it — the build-time
  flag comes from a GENCODE `exon_number` and is simply better than the
  coordinate heuristic. Legacy formats still parse; the 5-column one now logs
  a warning that E1 will not be strand-aware rather than quietly producing a
  plausible number.

  `mds_e1` now distinguishes three states instead of two: a value, `0.0` for
  "E1 exists but had no fragments", and **`NaN` for "this gene has no E1 at
  all"**. The last previously collapsed into `0.0` — the worst available
  choice, since MDS lives in `[0, 1]` and lower means more abnormal, so a
  fabricated `0.0` read as maximal tumour signal. It affects 88 of 128 `xs1`
  genes, and every gene when a strandless legacy BED is supplied.

  A legacy 5-column gene BED still parses and still produces per-exon MDS, but
  `region-mds` now **refuses to derive E1 from it** — without strand the
  heuristic returns the last exon for every minus-strand gene. `mds_e1` is
  `NaN` for that input, with a warning naming the fix.

  **Output-contract impact.** `MDS.gene.mds_e1` and `mds_e1_z` change for every
  minus-strand gene on panel data, and for panel genes whose canonical exon 1
  is not the most 5' captured region. Not comparable across the 0.9.0 boundary.

- **`region-mds.md` documented the wrong MDS scale, and contradicted itself.**
  The Formulas section showed an unnormalised Shannon entropy and a "~6.0 to
  ~8.0" range — raw bits, which the tool has never emitted — while the clinical
  table further down the same page quoted ~0.95–1.0. `motif_utils.rs` divides by
  `log2(256) = 8`, so MDS has always been in `[0, 1]`.

  Anyone who built a threshold from the formula section is out by a factor of 8.
  The divisor is now pinned in `validate/claims.py`, so the doc and the code
  cannot drift apart again silently.

- **Fragment coordinates were wrong whenever R1 was the rightmost mate.** The
  BED writer computed `pos() + |tlen|`, correct only when R1 is leftmost. For
  the other orientation the interval was shifted right by roughly
  `tlen - read_length` — measured at 105bp on a real read.

  > **Output change:** `{sample}.bed.gz` moves for every affected fragment, and
  > with it every positional feature — OCF end phasing, WPS fragment centres,
  > TFBS/ATAC and gene/exon overlap, and the GC value in BED column 4. Fragment
  > *lengths* are unchanged, so FSD/FSR/FSC size distributions are unaffected.
  > **Values from before this fix are not comparable**; PON baselines built on
  > uncollapsed input should be rebuilt.

  Measured incidence: **~48% of R1 reads on an uncollapsed BAM**, versus 0.4%
  after collapsing, 0.8% simplex and 0.0% duplex. Consensus callers normalise
  R1 to forward — Marianas (XS1) and fgbio (XS2) both — which is why this
  survived: the primary MSK-ACCESS path barely triggers it, while WGS and
  uncollapsed input run at the full rate.

  The pre-Rust implementation was strand-aware; the 2025-12 rewrite dropped the
  reverse branch while writing a correct sign-aware version forty lines below,
  in the same function. No commit message, doc, comment or test ever mentioned
  it. The same pattern was later copied into `region_mds.rs`, fixed here too.

  Branches on the TLEN sign rather than `is_reverse()`: htslib gives the
  rightmost segment a negative TLEN, and real data contains records flagged
  forward whose mate lies to their left, where the two disagree.

- **The fragment BED could be emitted unsorted**, which makes tabix indexing
  fail and takes every downstream feature with it. A consequence of the fix
  above: a fragment whose R1 is rightmost begins before the read that produced
  it, so read order stopped being coordinate order. The writer now sorts each
  chromosome before writing. Caught by the invariant suite, not by review.
- **Gene- and region-level FSC used different size bands than the genome-level
  table, so identically-named columns meant different things.** The genome-bin
  counter split at `<=100 / <=149 / <=220 / <=260 / <=400` with a sixth
  `ultra_long` channel. The gene/region aggregator split at
  `<100 / <150 / <260 / <400` with no sixth channel. `mono_nucl` therefore meant
  150–220 bp in `FSC.parquet` and 150–259 bp in `FSC.gene.parquet`; a column
  labelled `di_nucl` held the genome table's `long` range; and gene `long` held
  what the genome table called `ultra_long`. A consumer reading both tables into
  one feature matrix was combining channels that do not describe the same thing.

  The bands were never a deliberate gene-level choice. Three places already
  documented the *correct* (genome) bands for these tables: the column table in
  `docs/features/core/fsc.md`, the column table in
  `docs/reference/output-files.md`, and a line in the latter describing
  "146 genes × 6 channels" against an implementation that emitted five. Only the
  code disagreed, so this is a correction, not a redefinition.

  Both paths now call a single `SizeChannel::of`, and the bounds exist in
  exactly one place. Duplicating them is what allowed the drift, so the fix is
  the deduplication rather than a second edit to the same numbers. Gene and
  region tables gain `ultra_long` and `ultra_long_ratio`.

  A fragment that passes the length filter but matches no channel is now
  counted and reported once at write time. It is *not* logged where it is
  detected: that point is inside the rayon loop, and routing a per-fragment
  warning through pyo3-log into Python deadlocks under the GIL — observed while
  testing this change, where it hung the run outright rather than reporting
  anything.

  **Output-contract impact.** `FSC.gene.*` and `FSC.regions.*` change value for
  every fragment in a moved band, and gain two columns. The six ratios sum to 1;
  the five that existed before sum to `1 - ultra_long_ratio`, and the contract
  gate now asserts both so a downstream renormalisation cannot quietly use the
  wrong base. Values from before this fix are not comparable. `FSC.parquet` and
  `FSC.ontarget.parquet` are unchanged — the genome path was always correct.

- **`--generate-json` silently dropped most features for compressed and Parquet
  runs.** Every probe in `FeatureSerializer.from_outputs()` was
  `Path(f"{sample}.FOO.tsv").exists()`, and that gate ran *before* `read_table()`
  (which does understand `.tsv.gz` / `.parquet`). Reconstructing a real
  MSK-ACCESS xs2 output directory (written with both `--output-format both` and
  `--compress`, so every table exists as `.parquet` **and** `.tsv.gz`, but never
  as bare `.tsv`) shows the payload going from **5 to 16 feature families**:

  | | features captured |
  |---|---|
  | before | `mfsd`, `ocf`, `wps`, `wps_background`, `wps_panel` |
  | after | the above plus `fsd`, `fsr`, `fsc`, `fsc_gene`, `fsc_region`, `fsc_counts`, `motif`, `tfbs`, `atac`, `gc_factors`, `region_mds` |

  The three WPS entries survived only because they are probed as `.parquet`;
  `mfsd` and `ocf` survived only because a separate cleanup defect leaves stray
  uncompressed `.tsv` copies of those two outputs behind. Every fragmentomics
  feature — the entire size, motif and regulatory signal — was absent from the
  ML payload.

  > **Note:** this means tidying up those stray `.tsv` leftovers *without* this
  > fix would have made `features.json` strictly worse.

- **E1-only FSC was never generated for compressed or Parquet runs, and was
  misnamed when it was.** Two compounding bugs: `unified_processor` hard-coded
  `outputs.fsc_region = ...FSC.regions.tsv`, so the `outputs.fsc_region.exists()`
  guard was False whenever the real file was `.tsv.gz` / `.parquet` and E1
  generation was skipped entirely; and `filter_fsc_to_e1()` derived its output
  name via `Path.stem`, which strips only the last dot-segment, so a
  `.tsv.gz` input produced `S1.FSC.regions.tsv.e1only.tsv.gz` instead of
  `S1.FSC.regions.e1only.tsv.gz`. Added `strip_table_extension()` for compound
  suffixes.

- **EndMotif1mer was unreadable whenever `--compress` was used.**
  `write_end_motif_1mer()` gzipped the table via `write_table()` and then
  appended the `# c_fraction` / `# entropy` / `# c_bias` / `# sample` footer
  with a plain `open(path, "a")`. The result is a valid gzip member followed by
  raw bytes, which `gzip.decompress()` and pandas both reject with
  `BadGzipFile: Not a gzipped file (b'# ')` — so the **entire table** was lost,
  not merely polluted with footer rows. Every run using `--compress` (i.e. any
  pipeline emitting `.tsv.gz`) is affected. The footer is now written through
  the gzip codec, and `read_table()` gained a recovery path that decompresses
  the first gzip member and skips trailing bytes, so files already produced by
  <= 0.8.3 remain readable.


- **WPS: nucleosome repeat length (NRL) was a constant, not a measurement.**
  The Alu background profile covered only the ~300bp Alu body in 30 bins, so
  the DFT period grid was `300 / i` and the *only* index falling inside the
  nucleosomal search band was `i = 2` (150bp). Every sample therefore produced
  `nrl_bp = 150.0`, `periodicity_score = 0.3333` and `adjusted_score = 0.0`
  regardless of its data, and the documented "healthy NRL ~190bp, quality >
  0.7" was unreachable. The profile is now 2000bp (850bp flank + 300bp body +
  850bp flank) binned at 200 x 10bp, the DFT is zero-padded 8x, the search band
  is 140-250bp, and the deviation tolerance is the 50bp already documented
  (the code used 20bp, which forced `adjusted_score` to 0 given the pinned
  150bp NRL).

  > **Output change:** `{sample}.WPS_background.parquet` `nrl_bp`,
  > `nrl_deviation_bp`, `periodicity_score` and `adjusted_score` become
  > data-dependent. Earlier values are constants and must be discarded, not
  > compared.

- **UXM: the `X` (mixed-methylation) class was unreachable.** The CLI passed
  `methy_threshold = unmethy_threshold = 0.5`. Because the backend evaluates
  `ratio >= methy_threshold` first, every fragment collapsed into `M` or `U`
  and the published `X` column was identically `0.0` for every region of every
  sample. Thresholds are now `METHY_THRESHOLD = 0.75` / `UNMETHY_THRESHOLD =
  0.25`, matching the documented Loyfer et al. (2022) definition.

  > **Output change:** `{sample}.UXM.tsv` `X` becomes non-zero and `U`/`M`
  > shrink correspondingly. Models trained on the previous columns must be
  > refit.

- **region-MDS: E1 (first exon) selection ignored strand.** `identify_e1_regions`
  always chose the lowest start coordinate, so for every **minus-strand gene**
  the reported `mds_e1` was the gene's *last* exon — roughly half of a typical
  panel. E1 is now selected by transcription order (lowest start on `+`,
  highest on `-`). `write_gene_output` no longer re-derives E1 by coordinate;
  it reads the strand-aware `is_e1` flag, so `mds_e1` also stops silently
  falling through to the next covered exon when E1 has no fragments.

  > **Output change:** `{sample}.MDS.gene.tsv` `mds_e1` / `mds_e1_z` change for
  > all minus-strand genes. Earlier values are not comparable.

- **region-entropy crashed when the PON lacked a matching baseline.** The Rust
  `apply_pon_zscore` returns early *without writing an output file* if the PON
  has no `tfbs_baseline` / `atac_baseline` table. The Python caller then ran
  `load_entropy_tsv()` on the missing file, tripping its `assert df is not
  None`, and the raw Rust output was unlinked afterwards — losing the entropy
  results entirely. It now degrades to the documented no-PON behaviour
  (`z_score = 0`) with a warning.

- **Metadata footers parsed as data (`read_table`).** `{sample}.EndMotif1mer.tsv`
  appends `# c_fraction`, `# entropy`, `# c_bias` and `# sample` after the data
  rows, but `read_table()` called `pd.read_csv` without `comment="#"`. Those
  lines became data rows and `FeatureSerializer` propagated them into
  `features.motif.edm_1mer` as junk keys with NaN values. `read_table()` now
  defaults to `comment="#"` (callers may override).

  > **Output change:** `features.motif.edm_1mer` loses the four `"# ..."` keys.

- **The tool-level Nextflow path could not run at all.**
  `modules/local/krewlyzer/extract/main.nf` declared `path("*.json")` as a
  required output, but metadata moved from `.metadata.json` to
  `.metadata.tsv`/`.parquet` in **0.7.0**. The tool has not written a `.json`
  since, so the process failed on a missing required output. Unnoticed because
  `use_runall = true` is the default. Now emits the three real metadata
  variants as optional outputs, and the stub writes TSV to match.

- **Rust test suite did not compile.** `src/output_utils.rs` used
  `tempfile::tempdir` with no `[dev-dependencies]` section in `rust/Cargo.toml`
  and referenced `Arc` without importing it, so `cargo test` failed to build
  and no Rust unit test was running in CI. Also corrected the stale
  `gc_reference` bin-index test, which asserted a 68-bin ceiling although
  `get_length_bin_index` spans 60..999 across 188 bins.

- **`--generate-json` crashed on any output directory without a metadata
  table.** The metadata probe was the one call site that handed its resolver
  result straight to `read_table()` without an `is not None` guard, raising
  `AttributeError: 'NoneType' object has no attribute 'suffix'`. The
  format-matrix test always wrote a metadata table, so it never reached the
  branch; it now covers the absent case across tsv/tsv.gz/parquet.

- **E1-only FSC and FSC PON z-scoring were still skipped on compressed and
  Parquet runs.** `aggregate_by_gene()` returned the pre-conversion `.tsv` path
  even after `cleanup_intermediate_tsv()` had deleted it, so every caller
  guarding on `.exists()` saw a missing file. It now returns the path it
  actually wrote, which fixes the E1 guard, `apply_fsc_region_pon` and
  `apply_fsc_gene_pon` together; `outputs.fsc_gene` was additionally never
  resolved at all. The two "✓ … (with PON z-scores)" log lines no longer claim
  success when z-scoring bailed.

- Removed `src/krewlyzer.libs`, a broken symlink into a build sandbox
  (`/usr/local/lib/python3.11/dist-packages/`) committed by mistake.

- Restored `black` and `ruff` cleanliness. Both pass on `develop` but were
  broken by the changes above this line, and CI's lint job gates on them.

### Documentation

- Corrected feature documentation against the implementation: the inverted
  "higher MDS = tumor" summary in `motif.md`; mFSD LLR sigmas (`human()` is
  167/30 and 145/35, not 167/35 and 145/25), the canine/ssdna presets, their
  unreachability, and the `MIN_FOR_KS = 2` floor; the non-existent 20-500bp
  region-entropy window and its unusable absolute entropy bands; OCF's 1bp
  resolution and per-10,000 per-label normalization; the FSR channel table
  (FSR and FSC share one counter) and the fact that `total` spans 65-1000bp
  while only five channels are returned, so `channel / total` ratios do not
  sum to 1; and that WPS `z_vector` / `shape_score` / `z_amplitude` are
  documented but never emitted, while `apply_pon_zscore` silently yields
  `z = (0 - mean) / std`.

- Corrected three doc statements the pass above missed: `fsr.py`'s module
  docstring still listed the retired 65-99/100-149/150-259/260-399/400+ channel
  table (FSR reads FSC's per-bin counts, so it uses FSC's channels), and now
  also records that `short_frac` / `long_frac` divide by a `total` spanning
  65-1000bp; `estimate_periodicity`'s doc comment still claimed a 140-200bp
  search band after the NRL fix widened it to 140-250bp; and `mfsd.md` still
  asserted a KS floor of ≥5 in its tip and worked example, contradicting
  `MIN_FOR_KS = 2` and its own column reference.

  > Note: those retired FSR boundaries are not fictional — `GeneResult` in
  > `fsc.rs` still uses exactly that split, so `mono_nucl` means 150-259bp in
  > `FSC.gene`/`FSC.regions` but 150-220bp in `FSC`. Tracked separately.

## [0.8.3] - 2026-04-28

### Fixed
- **GIL Deadlock in Parallel Modules**: Wrapped all Rayon `par_iter()` calls in
  `py.allow_threads()` across `mfsd.rs`, `uxm.rs`, `extract_motif.rs`, and
  `region_mds.rs`. Previously, `pyo3-log` attempted to acquire the GIL from
  Rayon worker threads while the main thread held it, causing indefinite hangs
  on multi-threaded runs (observed as 16-hour wall time on IRIS HPC).
- **mFSD Silent Error Swallowing**: Replaced `Err(_) => continue` in BAM record
  iterator with a logged error breaker (max 1000 consecutive errors). Previously,
  corrupt BAM regions could cause infinite silent loops.
- **FILTER_MAF Substring Match**: Changed `sample_id in tsb` to exact match
  (`sample_id == tsb`) preventing `T01` from matching `T01-XS1`.
- **FILTER_MAF Comment Lines**: Stripped `#` comment lines from filtered MAF output
  in both multi-sample and single-sample modes.
- **mFSD 0-Variant Guard**: Added early exit at both Python (`wrapper.py`) and Rust
  layers when MAF has 0 data lines. Produces header-only TSV instead of attempting
  BAM access.
- **mFSD GC Correction Fallback**: When reference FASTA is unavailable or GC lookup
  fails for a region, GC correction is now skipped (weight=1.0) instead of silently
  using a hardcoded 50% GC.

### Added
- **mFSD BAM I/O Diagnostics**: Per-variant timing, BAM open/fetch latency logging,
  record counts, and slow-variant warnings (>30s). All diagnostics use structured
  `log` macros (visible with `--verbose`).
- **mFSD Header Constant**: Extracted 46-column TSV header into `MFSD_HEADER` module
  constant, eliminating duplication between normal output and 0-variant early-exit.

### Changed
- **Minimum Python**: `requires-python` raised from `>=3.8` to `>=3.10` (Python
  3.8 and 3.9 are EOL).
- **Minimum Rust**: Added `rust-version = "1.87"` MSRV to `Cargo.toml`.
- **Cargo Dependencies**: Updated 81 semver-compatible transitive dependencies
  (rayon 1.11→1.12, anyhow 1.0.100→1.0.102, flate2 1.1.5→1.1.9, etc.).
- **Diagnostic Logging**: Converted all `eprintln!` breadcrumbs inside parallel
  closures to structured `debug!`/`warn!` macros for proper timestamp/level
  formatting through Python's logging framework.

### Documentation
- **Developer Guide**: Added "PyO3 + Rayon GIL deadlock" Known Gotcha with correct
  and incorrect code examples.
- **Developer Guide**: Updated contributing checklist test count (244→357).

### Tests
- Added `test_mfsd_zero_variants` — verifies 0-variant input produces header-only TSV.
- Added `test_mfsd_maf_with_comment_lines` — verifies MAFs with `#` comment headers
  are parsed correctly.
- Total: 357 passed, 4 skipped.

## [0.8.2] - 2026-03-26

### Fixed
- **Region MDS Output Collision**: Fixed `Path.with_suffix()` compound-extension bug that
  caused `MDS.exon.tsv` and `MDS.gene.tsv` to both resolve to `MDS.tsv`, silently
  destroying per-exon and per-gene data. Affects standalone `region-mds`, `motif`, and
  `run-all` commands.
- **Motif Tracking Paths**: Fixed 5 additional `with_suffix()` bugs in motif output
  tracking (EndMotif, BreakPointMotif, MDS, MDS.ontarget, EndMotif1mer) that caused
  path collisions and silent PON z-score normalization failures.
- **Gene Format Conversion**: Region MDS gene output now correctly undergoes format
  conversion (Parquet/gzip) when no PON is provided.
- **MDS Z-score Logging**: Z-score append failures upgraded from `debug` to `warning`
  level to surface issues in production runs.

### Added
- **Compound Extension Tests**: 57 new parametrized tests (`test_compound_extension.py`)
  covering all 13 compound base names across TSV, Parquet, both formats, and roundtrip.
- **Developer Guide**: "Known Gotchas" section documenting the `Path.with_suffix()` anti-
  pattern with safe alternative and complete table of affected names.

### Changed
- **Nextflow Outputs**: Added missing Parquet emit declarations for TFBS sync (2), ATAC
  sync (2), ATAC/TFBS ontarget sync (2), and fsc_counts (4 — TSV+Parquet, off/on-target).
- **SLURM Script**: Tuned head process memory (16G→32G), queue size (100→200), added
  `--output_format both --compress_tsv true --verbose true` for 14K+ sample production runs.

## [0.8.1] - 2026-03-24

### Fixed
- **WPS Panel GC Correction**: Panel WPS now uses on-target GC correction factors instead
  of off-target. Panel anchors overlap capture regions, making on-target correction more
  accurate. Falls back to off-target factors when on-target unavailable.
- **Rust WPS Background**: Fixed `coitrees` metadata type mismatch in WPS background consumer
  (`wps.rs:1684`) — uses `.clone()` for cross-platform compatibility (returns `&usize` on
  macOS but `usize` on CI Linux).

### Changed
- **Dead Code Cleanup**: Removed 11 unused Python z-score functions from `pon_integration.py`
  (all replaced by Rust equivalents). Module reduced from 448 to 98 lines, 14 to 3 functions.
  Remaining: `load_pon_model`, `compute_nrl_zscore`, `compute_periodicity_zscore`.
- **GC Factor Resolution**: `gc_str`/`gc_ontarget_str` path resolution hoisted to a shared
  section in `unified_processor.py`, eliminating duplication between panel WPS and TFBS/ATAC.

### Added
- **build-pon Logging**: OCF on-target/off-target baseline status now logged during PON build.
- **sbatch Script**: `scripts/build_pon_unfiltered.sh` for building PON from high-coverage
  unfiltered BAMs on SLURM clusters.

### Data
- **PON Models Updated**: Rebuilt xs1/xs2 PON models for both `all_unique` and `duplex`
  variants with krewlyzer 0.8.x.

## [0.8.0] - 2026-03-17

### Added
- **On-target PON Z-scores**: Panel mode now computes on-target/off-target PON baselines
  for MDS, OCF, and FSD features, providing clinically-scoped z-scores in panel assays.
- **FSR On-target Output**: FSR now emits a separate `.FSR.ontarget.tsv` file in panel mode.
- **FSR Real Genomic Coordinates**: Region labels now reflect true genomic window coordinates
  instead of internal indices.

### Fixed
- **Output Format / Compress**: All output files now correctly respect `--output-format` and
  `--compress` flags (gzip compression, path handling, GC correction factors loading).
- **Rust GC Correction**: Added missing `PathBuf` import; fixed path handling for correction
  factor files.
- **Rust Clippy**: Replaced manual string strip with `strip_suffix` (`manual_strip` lint).

### Documentation
- Corrected stale FSR column names in `concepts.md` and `json-output.md`.
- Updated FSR on-target output docs and window/panel-mode descriptions.
- Documented `_core.pyi` stub maintenance requirements.

### CI
- **GitHub Actions Node.js 24**: Bumped `actions/checkout` v4→v5, `actions/setup-python` v5→v6,
  `actions/cache` v4→v5 to address Node.js 20 deprecation (enforced June 2, 2026).

## [0.7.0] - 2026-03-02

### Added
- **Configurable Output Formats**: `--output-format tsv|parquet|both` (default: `tsv`) controls
  all tabular feature outputs. `--compress` gzip-compresses TSV outputs (`.tsv.gz`).
  WPS outputs remain always-Parquet regardless of setting.

### Fixed
- **build-pon Intermediate Files**: Explicitly force `output_format="tsv"` and `compress=False`
  in all `process_sample()` calls within `build-pon` (both parallel and sequential paths).
  Prevents silent failure if default output format changes — intermediate files are internal
  scratch consumed by `pd.read_csv(sep="\t")`.
- **Feature Serializer**: Include `mds_z` in JSON output for the `from_outputs()` code path.
- **OCF Base File**: `OCF.tsv` = all reads (on + off combined), not off-target.
  `OCF.offtarget.tsv` is the true panel off-target score. Corrected in docs and code comments.
- **Rust wps.rs**: Remove erroneous `*` dereference on `node.metadata` (E0614 — `usize` is Copy).
- **Rust gc_correction.rs**: Prefix unused `valid_regions_path` param with `_` to silence
  compiler warning; parameter retained for API symmetry.
- **MkDocs Snippets**: Fix `--8<-- "CHANGELOG.md"` / `--8<-- "CONTRIBUTING.md"` broken includes
  by changing `pymdownx.snippets.base_path` from scalar `docs` to list `['.', 'docs']` so
  repo-root files resolve without path traversal blocked by CI deploy sandbox.

### Documentation
- **Post-0.6.0 Docs Audit** (12 files, 7 issue categories):
  - Fixed broken `--8<-- "../CHANGELOG.md"` / `--8<-- "../CONTRIBUTING.md"` snippet includes
  - Corrected outdated "no global `--output-format` flag" note in `cli/run-all.md`
  - Updated `metadata.json` → `metadata.tsv` across 7 files (8 references total)
  - Added WPS always-Parquet exception note to `reference/output-files.md`
  - Added build-pon intermediate TSV format note to `guides/building-pon.md`
  - Updated test count (248 → 244 + 4 skipped) and added CI lint steps to `developer-guide.md`
  - Added `--output_format` and `--compress` parameters to `nextflow/parameters.md`
- **Output Format Options section**: New section in `reference/output-files.md` documenting
  `--output-format`, `--compress`, WPS always-Parquet exception, and `--generate-json`
- **OCF Variant Clarification**: Added 3-variant table and note block in `reference/output-files.md`
  explaining `OCF.tsv` (all reads) vs `OCF.ontarget.tsv` vs `OCF.offtarget.tsv`
- **docs/index.md**: Replaced `:latest` Docker tag with explicit `:0.7.0` per release policy

### CI
- **Lint Job**: Added parallel `lint` job (`ruff · black · mypy · cargo clippy -- -D warnings`)
  running alongside tests on all push/PR events

## [0.6.0] - 2026-02-28

### Added
- **mFSD Base Quality Filtering**: `--min-baseq` / `-Q` (default 20) gates variant evidence by base quality
- **mFSD GC Correction**: Rust-native LOESS GC bias correction for variant fragment size distributions
- **mFSD Duplex Weighting**: Proper consensus fragment handling via `--duplex`
- **Region MDS `--sample-name`**: Consistent output naming without post-hoc rename
- **Feature Serializer**: Auto-load `fsc_counts`, `region_mds`, `uxm` in `from_outputs()`
- **IRIS Batch Submission**: `scripts/run_krewlyzer_iris.sh` for SLURM/IRIS cluster runs with `--generate_json`
- **nf-core Institutional Configs**: `custom_config_base` param and IRIS profile support
- **Versioned Documentation**: Implemented `mike` for dev/stable doc versions
- **Nextflow mfsd Module**: Full standalone params (`--reference`, `--correction-factors`, `--mapq`, `--minlen`, `--maxlen`, `--min-baseq`, `--duplex`, `--no-skip-duplicates`)
- **Nextflow runall**: `fsc_counts.tsv` output declaration, `--min-baseq` wired
- **mFSD Integration Tests**: 161 lines of new test coverage

### Fixed
- **mFSD MAF Parsing**: Header-based column lookup (fixes column-index mismatch with different MAF flavors)
- **Nextflow FILTER_MAF**: Complete overhaul — eliminated join operator blocking, replaced regex with substring match, fixed SyntaxError in versions.yml, dynamic maxForks for SLURM
- **Nextflow Workflow Streaming**: Fixed RUNALL blocking from `remainder:true`, `failOnMismatch`, channel round-robin; used `multiMap` for proper routing
- **Nextflow RUNALL Outputs**: Added 14 output declarations, fixed BreakPointMotif casing, explicit publishDir
- **Region MDS Nextflow**: `--sample-name` replaces `mv` workaround
- **Nextflow Config**: Executor queueSize placement, `-qs` CLI flag, global publishDir removal
- **WPS CLI Tests**: Fixed `--input` flag (was positional arg) — recovered 2 skipped tests
- **Pandas FutureWarning**: Fixed `pd.concat` with all-NA columns in PON test fixture

### Changed
- **Code Quality**: Black formatted 71 files, ruff fixed 129 lint errors, cargo clippy applied
- **Ruff Config**: Added `[tool.ruff]` to `pyproject.toml` with documented E402/F821 ignores
- **Agent Config**: `.agent/` → `.agent/rules/` with `alwaysApply` frontmatter

### Documentation
- **45-item Audit**: Corrections across 25 doc files including `.csv`→`.tsv` (7 files), `.WPS.tsv.gz`→`.parquet` (3 files), phantom `--output-format` removed, Docker versions→`X.Y.Z`, parameters.md 12→28, outputs.md 14→41, JSON schema corrected, developer guide Rust table 10→19, architecture pipeline signature updated
- **PDF Embedding**: Fixed rendering with mkdocs-pdf plugin

## [0.5.3] - 2026-02-06

### Documentation
- **Admonition Formatting**: Converted blockquotes to proper MkDocs admonitions across 10+ docs files for consistent styling
- **Table Rendering**: Fixed tables not rendering by adding required blank lines before table headers
- **LaTeX Formulas**: Fixed underscore escaping in math blocks (mFSD, FSR, WPS formulas)
- **Abbreviations**: Added glossary with auto-append for consistent acronym tooltips (cfDNA, WPS, etc.)
- **Docs CI**: Added PR validation trigger with `mkdocs build --strict` to catch issues before deployment

## [0.5.2] - 2026-02-05

### Added
- **Dual BAM Support (mFSD)**: New `--mfsd-bam` parameter in `run-all` to use a dedicated duplex BAM for mFSD while other features use the main all_unique BAM
  - Auto-duplex detection from filename patterns (`-duplex`, `_duplex`, `.duplex`)
  - Startup banner shows mFSD BAM path when specified
  - Nextflow: Added `mfsd_bam` column to samplesheet schema

### Fixed
- **Correction Factors Loading**: Fixed delimiter detection for `.correction_factors.tsv` files (CSV format with TSV extension)

### Documentation
- **Release Guide**: Added version format note (use `0.5.2` not `v0.5.2` in code)
- **mFSD docs**: Updated with dual BAM workflow examples
- **Samplesheet docs**: Added `mfsd_bam` column documentation
- **Footer**: Added site footer with author attribution and Antigravity acknowledgment
- **Theme**: Changed color scheme to blood-red (#ef5552) to match "krew" (blood) branding

## [0.5.1] - 2026-02-04

### Fixed
- **Docker Image**: Data folder now bundled in container (was missing in 0.5.0)
- **CI Tests**: Tests requiring bundled data now skip in PyPI installs via `@requires_data` marker

### Changed
- **Docker Tags**: Use versioned tags only (e.g., `:0.5.1`); no `:latest` tag published

### Documentation
- **Installation**: Added Singularity/Apptainer section for HPC clusters
- **Installation**: Clarified Docker uses release tags, not `:latest`
- **Nextflow Examples**: Added IRIS HPC cluster profile example (`-profile iris`)
- **Testing Guide**: Added Data Availability section explaining PyPI vs source
- **Structure**: Moved Changelog to Development section; removed duplicate files
- **Mermaid Diagrams**: Upgraded to `mkdocs-mermaid2` plugin for theme-aware rendering
- **Removed**: Full-site PDF export (unreliable rendering)

## [0.5.0] - 2026-02-02

### Added

#### PON Framework
- **`build-pon` command**: Generate PON models from cohort samples with FSD/WPS/OCF/MDS baselines
- **Bundled PON assets**: Pre-computed xs1/xs2 PON models for `all_unique` and `duplex` variants
- **`--pon-variant` flag**: Select between `all_unique` (max coverage) or `duplex` (highest accuracy) PON models
- **Duplex tag warning**: mFSD warns when `--duplex` is used but no cD tags found in BAM

#### Panel Mode (MSK-ACCESS)
- **Assay-aware asset resolution**: Auto-load targets, anchors, PON via `-A/--assay xs1|xs2`
- **On/Off-target splitting**: All tools produce separate `.ontarget.tsv` outputs
- **Gene-level FSC**: New `{sample}.FSC.gene.tsv` and `{sample}.FSC.regions.e1only.tsv`
- **Bait edge masking**: WPS `--bait-padding` to remove capture artifacts
- **`--skip-target-regions`**: Force WGS mode for panel assays

#### Duplex Sequencing
- **`--duplex` flag (mFSD)**: Enable family size (cD tag) weighting for duplex BAMs
- **LLR scoring**: Log-likelihood ratio classification for cross-species support
- **GC-weighted mFSD**: 5 new GC-corrected mean fragment size columns

#### Region-Based Analysis
- **Region Entropy**: TFBS/ATAC dual-output architecture with per-cluster entropy
- **Region MDS**: Per-gene Motif Diversity Score with E1-only filtering
- **Rust backend**: New `region_entropy.rs` for high-performance calculation

#### Feature Enhancements
- **Jagged Index**: 1-mer End Motif analysis with C-end fraction
- **WPS v2.0**: Hierarchical stacking, extended anchors, panel-specific normalization
- **Output formats**: `--format` flag (tsv/parquet/json) for all tools
- **`--generate-json`**: Unified JSON export with all features

#### Infrastructure
- **Git LFS**: Large files (.gz, .parquet, .bed) tracked via git-lfs
- **BGZF compression**: All BED outputs use block-gzip format
- **GC references**: Pre-computed Parquet format for GRCh37/GRCh38
- **Unified processor**: Single-pass Rust pipeline for Extract→Motif→FSC→FSD→WPS→OCF

### Changed

#### Nextflow Pipeline
- **nf-core compliance**: Refactored to shared INPUT_CHECK subworkflow
- **KREWLYZER_RUNALL**: Unified module with 30+ output channels
- **`pon_variant` parameter**: New pipeline parameter (default: `all_unique`)

#### CLI
- **9 tools with `--pon-variant`**: run-all, fsc, fsd, fsr, wps, ocf, motif, region-entropy, region-mds
- **Centralized asset resolution**: `resolve_pon_model()` and `resolve_target_regions()`
- **Filter flags**: `--mapq`, `--minlen`, `--maxlen`, `--skip-duplicates`, `--require-proper-pair`
- **Parallel processing**: `--threads` option across all tools

#### Data & Assets
- **Terminology**: Renamed "blacklist" to "exclude_regions" for inclusive language
- **Directory structure**: `data/{type}/{genome}/{variant}/` for organized assets
- **Memory optimization**: Parallel sample processing with spawn context in build-pon

### Fixed
- mFSD: Filter discordant reads with extreme TLEN values
- mFSD: Fix verbose mode hanging by moving debug logging outside parallel loop
- Proper pair detection for legacy BAMs without proper pair flags
- GC correction for gene-level FSC in panel mode
- CIGAR handling improvements for INDELs and complex variants
- BAM extraction for v1 ACCESS and bgzip output format

### Documentation
- **11 docs updated**: PON Variant Selection across pon.md, pipeline.md, usage.md, feature docs
- **Release Guide**: New `docs/advanced/release-guide.md` for version release process
- **Math rendering**: LaTeX formulas in all feature documentation
- **Glossary**: New terms and definitions

### Tests
- **28 new test files**: Unit, integration, and e2e coverage
- **4 PON variant tests**: test_asset_resolution.py
- **Real data tests**: test_real_data.py for end-to-end validation

## [0.3.2] - 2025-12-18

### Fixed
- **CI Build**: Removed `gfortran` and `scikit-misc` - GC correction now fully in Rust
- **FSC GC Correction**: Added `--gc-correct/--no-gc-correct` and `--verbose` flags to FSC CLI
- **Python LOESS Removed**: Removed `gc_correct()` from `helpers.py` and `postprocess.py`

### Changed
- **Dockerfile**: Removed `gfortran` dependency (no longer needed)
- **GitHub Actions**: Removed `gfortran` from CI workflows
- **Nextflow Modules**: Updated all container versions to `0.3.2`

## [0.3.1] - 2025-12-18

### Added
- **Rust LOESS GC Correction**: New `rust/src/gc_correction.rs` module using the `lowess` crate
  - Per-fragment-type correction (short, intermediate, long)
  - Configurable LOESS parameters (fraction, iterations, delta)
- **FSR GC Correction**: `--gc-correct/--no-gc-correct` flag (default: **True**)
  - Uses corrected counts from Rust before calculating ratios
  - `--verbose` flag for detailed logging
- **WPS GC Correction**: `--gc-correct/--no-gc-correct` flag (default: **True**)
  - `--reference, -r`: Reference FASTA for computing region GC content
  - FASTA-based GC computation using rust-htslib::faidx
  - Graceful fallback if reference not provided

### Changed
- **FSC Rust Backend**: Added `count_fragments_gc_corrected()` function for integrated GC correction
- **WPS Rust Backend**: Updated `calculate_wps()` with `reference_path`, `gc_correct`, `verbose` parameters
- **lowess Dependency**: Updated from 0.3 to 0.4 for improved API

### Documentation
- Updated `docs/features/fsr.md` and `docs/features/wps.md` with GC correction options


## [0.3.0] - 2025-12-16

### Added
- **mFSD Variant Type Support**: Now handles all small variant types:
  - SNV (single nucleotide)
  - MNV (multi-nucleotide)
  - Insertion
  - Deletion
  - Complex (substitution + indel)
- **4-Way Fragment Classification**: Fragments classified as REF, ALT, NonREF, or N (low quality)
- **Comprehensive Statistics**: 6 pairwise KS comparisons (ALT-REF, ALT-NonREF, etc.)
- **Derived Metrics**: VAF_Proxy, Error_Rate, N_Rate, Size_Ratio, Quality_Score
- **Quality Flags**: ALT_Confidence (HIGH/LOW/NONE), KS_Valid
- **Distributions Output**: `--output-distributions` flag generates per-variant size distributions
- **Verbose Logging**: `--verbose` flag for debug-level logging with variant type breakdown
- **MRD Support**: Proper handling of low fragment counts (1-2) common in MRD settings

### Changed
- **BREAKING**: mFSD output format changed from 11 columns to 39 columns
- **Nextflow MFSD**: Distributions and verbose logging enabled by default
- **Fragment Counting**: Now counts fragments (R1 only) instead of reads to avoid double-counting

### Fixed
- **CIGAR Handling**: Improved sequence extraction for INDELs and complex variants

## [0.2.3] - 2025-12-16

### Changed

- **Nextflow Pipeline**: Added `FILTER_MAF` module for multi-sample MAF filtering:
  - New `maf` and `single_sample_maf` columns in samplesheet
  - Filters MAF by `Tumor_Sample_Barcode` matching sample ID (regex: `.*sample_id.*`)
  - `single_sample_maf=true` bypasses filtering for per-sample MAFs
  - Skips MFSD when filtered MAF has zero variants (with warning)
  - Memory-efficient streaming for large MAF files (100s of MBs)
  - Outputs filtered MAF + mFSD results
- **Nextflow Modules**: Updated all container versions to `0.2.3`

## [0.2.2] - 2025-12-15

### Changed
- **Project Structure**: Migrated to recommended maturin "src layout" for better Python/Rust separation:
  - Python code moved from `krewlyzer/` to `src/krewlyzer/`
  - Rust code moved from `krewlyzer-core/` to `rust/`
  - Single `rust/Cargo.toml` (removed duplicate root `Cargo.toml`)
  - Rust module now imports as `krewlyzer._core` (private)
- **Dockerfile**: Rewritten as multi-stage build for smaller image size (~200MB vs ~1GB); amd64 only
- **OCI Labels**: Added container metadata for GitHub Container Registry

### Fixed
- **Distribution Compatibility**: Updated release workflow to build `manylinux_2_28` wheels (compatible with RHEL 8+, AlmaLinux 8, etc.)
- **Source Builds**: Included `clang` and `llvm-devel` in build environment for `bindgen`/`hts-sys`
- **Docker Build**: Added `patchelf` to maturin installation for Linux wheel building
- **Test Imports**: Updated test files to use new `krewlyzer._core` import path
- **CI OpenSSL**: Use system OpenSSL (`OPENSSL_NO_VENDOR=1`) instead of building from source
- **CI Python Versions**: Build wheels for Python 3.10, 3.11, 3.12 in release workflow
- **Documentation**: Updated logo paths for new `src/` layout

## [0.2.1] - 2025-12-15

### Fixed
- **Rust Compilation**: Resolved cross-platform build issues with `coitrees` metadata types (`usize` vs `&usize`) by using explicit `.to_owned()` conversion. 
- **CI Build**: Added `gfortran`, `clang`, and `libclang-dev` to CI workflows to fix `scikit-misc` and `rust-htslib` compilation failures.
- **Permission Errors**: CI scripts now robustly handle `sudo` permissions when installing system dependencies.

## [0.2.0] - 2025-12-12

### Added
- **Unified Engine**: New high-performance Rust core (`krewlyzer-core`) that processes Extract, Motif, FSC, FSD, WPS, and OCF in a single parallelized pass.
- **Fragment Extraction**: dedicated `extract` command (via Rust) to convert BAM to BED with configurable filters.
- **Extract Documentation**: New `docs/features/extract.md` detailing extraction logic and JSON metadata.
- **Calculation Details**: Comprehensive formulas and interpretation guides added to all feature documentation.
- **Root Cargo.toml**: Added to support standard `maturin` builds for the hybrid Python-Rust package.

### Changed
- **Performance**: Significant speedup (3-4x) for end-to-end analysis due to multi-threaded processing and single-pass I/O.
- **Build System**: Migrated to `maturin` backend for robust Rust extension compilation.
- **CLI (`run-all`)**: Now defaults to the Unified Engine.
- **CLI Filters**: Added `--mapq`, `--minlen`, `--maxlen`, `--skip-duplicates`, `--require-proper-pair` flags to `run-all`, `extract`, and `motif`.
- **Motif Outputs**: Renamed output files to use `.tsv` extension consistently (e.g., `{sample}.EndMotif.tsv`).
- **Data Handling**: `motif` now uses the unified engine, eliminating the need for `bedtools` binary entirely.
- **Documentation**: Updated `README.md`, `usage.md`, and `pipeline.md` to reflect the new workflow.
    - Corrected `pipeline.md` samplesheet format documentation to match Nextflow schema.
    - Updated `usage.md` and feature docs to correctly specify output directory arguments.

### Fixed
- **Test Suite**: Cleaned up `tests/` directory, removing obsolete scripts and fixing integration tests (`test_science.py`, `test_run_all_unified_wrapper.py`).
- **Validation**: Fixed BAM header issues in tests.

### Removed
- **Legacy Python Backends**: Removed pure Python implementations of `extract`, `motif`, `fsc`, `fsd`, ensuring all paths use the unified Rust core.
- **Validation Artifacts**: Deleted temporary validation scripts and data.

## [0.1.7] - 2025-11-26

### Fixed
- **PyPI Metadata**: Added `readme` and `license` fields to `pyproject.toml` to ensure the long description is correctly displayed on PyPI.

## [0.1.6] - 2025-11-26

### Fixed
- **Docker Build**: Removed `libatlas-base-dev` dependency from `Dockerfile` as it is not available in the `python:3.10-slim` (Debian Trixie) base image.

## [0.1.5] - 2025-11-26

### Fixed
- **Docker Publishing**: Switched to `GITHUB_TOKEN` for GHCR authentication to fix permission issues.
- **PyPI Publishing**: Added `skip-existing: true` to handle existing versions gracefully.
- **CI/CD**: Added build checks for Python package and Docker image to PR workflows.

## [0.1.4] - 2025-11-26

### Fixed
- **Test Dependencies**: Removed unused `pybedtools` imports from `fsr.py`, `fsd.py`, `uxm.py`, and `fsc.py` which were causing `ImportError` in CI environments where `pybedtools` is not installed.

## [0.1.3] - 2025-11-26

### Changed
- **Dependency Reduction**: Removed `pybedtools` dependency.
- **Refactor**: `motif.py` now uses `pandas` for blacklist filtering and sorting, removing the need for `bedtools` binary.
- **CI/CD**: Added `pytest` and `pytest-mock` to `test` optional dependencies in `pyproject.toml`.

## [0.1.2] - 2025-11-26

### Added
- **Mutant Fragment Size Distribution (`mfsd`)**: New module to compare fragment size distributions of mutant vs. wild-type reads using VCF/MAF input.
- **Enhanced Fragment Size Ratios (`fsr`)**: Added "Ultra-Short" (<100bp) ratio bin.
- **Documentation**: Comprehensive MkDocs website (`docs/`) with material theme.
- **Pipeline**: `run-all` command now supports `--variant-input` for `mfsd` analysis.
- **Nextflow**: Pipeline updated to support optional variant input in samplesheet.

### Changed
- Updated `README.md` to point to the new documentation site.
- Added `mkdocs` and `mkdocs-material` as optional dependencies.
