# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Fixed

- **`describe-output -o page.html` wrote unrendered Markdown.** The file was a
  valid HTML document — doctype, head, styles — whose body was the Markdown
  source escaped into a `<pre>`, so it displayed `##` headings and `|` tables
  as literal text. Nothing said why, and the CLI reported `wrote <path>`.

  `render_html_page` converts the Markdown only when the `markdown` package is
  importable, and `markdown` was declared in **no dependency group at all** —
  not core, not `report`, not `dev`. No install of krewlyzer could satisfy it,
  so every user took the fallback; it worked only where something unrelated
  (mkdocs) had pulled the package in.

  `markdown>=3.4` now sits in the `report` extra, alongside `plotly` — that
  group is human-readable output, which is what both commands produce. The
  fallback also announces itself now, in the page and on the console, matching
  `krewlyzer report`, whose charts each state when plotly is absent rather than
  quietly rendering nothing. Asking for `.html` and getting Markdown is
  confusing precisely because the file *is* HTML; nothing looks wrong except
  the contents.

- **`KREWLYZER_DATA_DIR` is now honoured by every bundled-data lookup, not just
  one.** The README prescribes it for pip installs — the wheel excludes
  `data/` to stay under PyPI's 100 MB limit, so those users keep a separate
  clone and point at it. `AssetManager` honoured it; `get_bundled_gene_bed` and
  `build-gc-reference` did not, and each resolved against the installed
  package instead.

  For a pip user, that meant panel assets loaded correctly while
  `get_bundled_gene_bed` returned `None` — and `None` is not an error here, so
  gene-level FSC and region-MDS were simply absent from the output with only a
  debug-level log to say so. `build-gc-reference` degraded more loudly,
  warning that regions "may include problematic areas", but degraded all the
  same.

  The resolution rule now lives once, in `assets.bundled_data_dir()`. Where the
  variable is set it is the whole answer — resolve there or return nothing,
  never fall back to the bundled copy, matching what `AssetManager` already
  did. A fallback would demote the override to a suggestion and mask a
  misconfigured data clone behind whatever the package happened to ship.

  Nothing caught this because the tests that would have were skipped in CI for
  precisely the same reason: `conftest.DATA_AVAILABLE` asked the dataless
  installed package too. It now resolves the way the runtime does, which puts
  16 previously-skipped tests into every CI run — and they found this on the
  first one.

## [0.9.1] - 2026-08-15

### Added

- **`validate-output --expect`** — reconcile the samples you meant to run
  against the ones that exist. Everything else in the gate validates the
  samples it *finds*, and two failure modes leave nothing to find:

  | | before |
  |---|---|
  | No output directory at all | invisible |
  | Directory exists, no Parquet written | **invisible** — `discover_samples` skips it |
  | Parquet present, no completion marker | ERROR ✓ |

  A job killed between `mkdir` and its first write was indistinguishable from
  one never submitted. At a handful of samples that is obvious; across a
  16,000-sample cohort it is not, and the consumer reads a short cohort without
  knowing it was short.

  Takes a Nextflow samplesheet (using its `sample` column) or a plain list. The
  two absences are reported separately, because an empty directory means the
  job ran and its logs exist while a missing one usually means it was never
  submitted. An unexpected sample is a WARN, not an ERROR — it corrupts
  nothing, but it usually means an identifier that does not round-trip, which
  makes every other count unreliable. A CSV without a `sample` column is an
  error rather than a guess at the first column: reconciling against the wrong
  field reports every sample missing, which reads as a catastrophic run failure
  instead of a malformed input.

  The summary table shows counts; the `--json-report` names every affected
  sample, which is the split that keeps it readable at cohort scale.

- **A check that a release actually published.** Every other gate verifies the
  repository; nothing asked the registries whether the release landed. 0.9.0
  passed every workflow while the docs site served 0.8.3 and no GitHub Release
  existed — both invisible, because green CI only means the configured steps
  ran, not that the steps you need exist.

  `scripts/check_release_artifacts.py` compares the newest tag against PyPI,
  GHCR, the docs `stable` alias, and the Release object. It runs weekly, and
  after the `Release` workflow finishes via `workflow_run` — not a tag trigger
  with a `sleep`, which would idle a runner for 40 minutes and cry wolf
  whenever publishing ran long. The weekly schedule is the load-bearing one: a
  post-release check alone would not have caught the docs, because that deploy
  had been deleted three weeks *before* 0.9.0.

  Three exit codes, because "could not reach the registry" is not "in sync":
  0 matches, 1 stale, 2 unreachable. A network failure reported as success is
  the silent pass the whole script exists to prevent.

- **`check_phi_guard.sh` assertion 6** — every tracked hook exists, matches the
  staged version, and `pre-push` still delegates to git-lfs. Compared against
  the index rather than `HEAD`, so a deliberate edit passes while an
  out-of-band overwrite fails.

- **A stale-extension warning at pytest session start.** Python imports whatever
  `maturin develop` built last, so after a `rust/src/` change the suite tests the
  previous binary — passing, and proving nothing. In `conftest.py` rather than a
  per-tool hook so it fires for a human, for CI, and for any agent. A warning
  and never a failure: the comparison is mtime-based and `git checkout` rewrites
  mtimes, so a branch switch can trip it.

- **`release.yml` now creates the GitHub Release**, using the CHANGELOG section
  for the tag as the body. Until now the Releases page was hand-maintained,
  which is why 0.9.0's tag went green with nothing on it. A missing section
  fails the step; a section over 60,000 characters is cut at a heading boundary
  and linked (0.9.0's was 103,376 against GitHub's 125,000 limit).

- **`MIN_PON_VERSION` and `KREWLYZER_ALLOW_OLD_PON` pinned in the claims
  registry**, now that three documents quote them.

### Changed

- **FILTER_MAF no longer runs once per sample when the MAF is shared.** It
  scanned the entire cohort MAF, line by line in Python, to extract one
  sample's rows — so a 16,552-sample run performed 16,552 full scans of the
  same file, each as its own SLURM job whose scheduling overhead exceeded the
  work it did. `maxForks queue_size.intdiv(2)` then reserved **250 of 500
  slots** for it, upstream of RUNALL, so the opening wave of a run was nothing
  but MAF filtering and no real work could start. This was the observed queue
  clog on the 0.8.3 cohort.

  The new `SPLIT_MAF` does one pass per distinct MAF: **1 job instead of
  16,552**, and the whole queue freed for RUNALL.

  Batching keys on **(maf, single_sample)**, not the MAF alone — grouping a
  pass-through sample with a filtered one would apply a single mode to both and
  silently corrupt one. A sheet with per-sample MAFs yields one task each and
  behaves exactly as before; samples with no MAF keep FILTER_MAF and its
  header-only placeholder.

  Every requested sample gets a file, including those with zero matching
  variants, and `SPLIT_MAF` exits non-zero if it wrote a short set. A sample
  with no file would be dropped by the caller's re-pairing and RUNALL would
  never run for it — silent cohort loss reintroduced one stage earlier than
  `--expect` can see it. The metas ride through the process so re-pairing is a
  `map`, preserving the join-free streaming this workflow is built around.
- **Reverted the RUNALL retry ladder added moments earlier in this cycle.** It
  was built on the claim that the pipeline had no retry policy. That claim was
  wrong: `queue`, `errorStrategy` and `maxRetries` come from the institutional
  config `-profile iris` loads over HTTPS at parse time, so they appear nowhere
  in this repository and a `grep` finds nothing. iris sets `maxRetries = 3`,
  retries any failure, falls back to `ignore`, and selects the partition with a
  closure honouring `--isolated` and `--partition`.

  Because `withName` outranks a generic `process` block, the ladder silently
  replaced all of that — cutting retries from 3 to 1, ignoring `--partition`
  and `--isolated`, and pinning attempt 1 to the smaller of two partitions. The
  `task.attempt` ramps it claimed to revive had been live all along.

  The pre-existing behaviour is better and needs no queue names: a growing time
  request makes SLURM pick the partition, since attempt 2's 4h no longer fits a
  3h queue.

  `tests/unit/test_runall_resources.py` now pins the *absence* of `queue`,
  `errorStrategy` and `maxRetries` from that block, so the next reader who
  greps, finds nothing and "fixes" it lands on the explanation instead.
- **`check_phi_guard.sh` runs in 2.3s instead of ~100s.** `git archive` applies
  the smudge filter, so assertion 1 was unpacking and scanning 660 MB of
  parquet, gzipped BED and FASTA (75s), and assertion 5 grepped the same
  binaries for `/Users/` (13s). `GIT_LFS_SKIP_SMUDGE` exports the 132-byte
  pointers instead, and assertion 5 skips LFS paths. This gives up scanning
  inside LFS payloads: accepted, because they are machine-built compressed
  reference data where a plaintext regex could never have matched.

- **Release guide Phase 5 merges `main` into `develop`**, not the release branch.
  The old form left each branch holding a merge commit the other lacked, so
  `main` sat permanently one ahead; merging `main` makes it an ancestor and
  `git merge-base --is-ancestor main develop` becomes a valid post-release check.

- **The gates table in `AGENTS.md` has a "Read first" column.** `.agents/rules/`
  was described as read "on demand" with nothing creating demand, so those files
  went unopened for a whole release.

### Fixed

- **A 16,552-sample run produced zero samples.** Three faults compounded:

  1. **SPLIT_MAF interpolated the sample ids into its script.** Nextflow strips
     a script block's common indentation *after* interpolating, so 16,552
     zero-indent id lines left no common prefix, nothing was stripped, and
     python saw `    import os` on line 2. `.command.sh` was 16,648 lines: 96
     of code, 16,552 of ids. Exit 1. The ids are now a **staged file** — a
     filename is one token and cannot defeat the dedent at any cohort size.
  2. **Retries could not submit.** The scheduler runs `EnforcePartLimits=ALL`,
     where a partition list enforces the *intersection* of limits and offers
     the union of nodes. `2.h * task.attempt` asked 4h on attempt 2 against a
     list containing a 3h partition and was refused — *"Requested time limit is
     invalid"* — not rerouted. New `--long_partition` sends retries to a
     partition whose limit covers the grown request.
  3. **`ignore` on a fan-in step lost everything.** SPLIT_MAF is one task for
     the whole cohort, so ignoring it left RUNALL with no input and the run
     completed having done nothing. It is now `terminate`: a fan-in failure
     stops the run in ninety seconds instead of after the queue drains.

  The unit tests could not have caught (1): the harness dedents *before*
  substituting, Nextflow does the reverse, so it tested well-formed Python while
  production shipped broken. The stub run could not either — its ids are joined
  with spaces, one line, nothing to defeat. `test_module_interpolations.py` now
  refuses any newline-joined value in a `.nf`.

- **Nextflow 26.04 could not parse `nextflow.config` at all.** The nf-core
  `try { includeConfig } catch` idiom is rejected outright — *"Try-catch blocks
  cannot be mixed with config statements"* — so the pipeline refused to start.
  25.10.3 accepted the same file, which is how it went unnoticed: the only
  signal was a version nobody had run. Replaced with the ternary nf-core moved
  to, which is why their own configs already parse under 26.

  The behaviour change is the better half. The old `catch` swallowed **any**
  failure, so a transient network blip left the run going with no institutional
  config — no `process_single` label, no queue closure honouring
  `--isolated`/`--partition`, no `maxRetries` — and three processes silently
  falling back to 8 CPU / 32 GB / 24h, a request a short partition refuses at
  submission. A whole run's resource policy changed, reported as one line of
  stderr. Skipping now requires the operator to set `NXF_OFFLINE`.

  `tests/unit/test_nextflow_config_parses_on_26.py` pins the construct. It is
  not a parser and cannot claim the file is 26-clean — the parser stops at the
  first error, so more may sit behind any given one. Only a real
  `NXF_VER=26.04.6 ... -stub-run` can establish that.

- **SPLIT_MAF's stub block still referenced the pre-rename input.** Its input
  became `metas` so the metas could ride through and the caller could re-pair
  with a `map` instead of a `join`; the script block and tag were updated and
  the `stub:` block was not. Groovy resolves `${...}` at task runtime, so the
  file parsed, every unit test passed, and `nextflow -stub-run` failed
  immediately with `No such variable: sample_ids`.

  `tests/unit/test_module_interpolations.py` now asserts, for every module,
  that the stub block references no name absent from its inputs and script —
  which is exactly the shape of a rename applied to one block and not the
  other. Written as a subset comparison between the two blocks rather than a
  declaration check: parsing Nextflow's input syntax precisely (`path(x)`, bare
  `path x`, multi-element tuples) produced false positives on seven modules,
  and a test that cries wolf is worse than no test.

- **The Nextflow pipeline reported itself as 0.8.3.** `manifest.version` sat at
  the previous release through all of 0.9.0. Nothing functional depended on it,
  which is why it drifted — but Nextflow writes it into every execution report,
  trace file and Tower entry, so a 16,000-sample cohort would have been labelled
  0.8.3 in its own provenance while running 0.9.0 containers.

  The Phase 2 version bump is a `sed` over container tags in the modules and
  never matched this line. The release guide compounded it by listing
  `nextflow.config` under "container tag" — there is no container pin in that
  file at all. Both corrected, and `tests/unit/test_nextflow_version.py` now
  pins `manifest.version` and every module container tag to `__version__`.

  The container half matters more than the manifest: a stale pin there does not
  mislabel a run, it executes the previous release's defects under this
  version's name.

- **The documentation site stopped publishing on 2026-07-31 and nobody noticed
  for a whole release.** `ba13fd4` (#32) replaced the docs workflow with the
  aggregate check alone, dropping the `mike` deploys and the strict build with
  it. Nothing has published since: `stable` still served **0.8.3** — the
  documentation for the release whose defects 0.9.0 exists to fix — and `dev`
  froze at the same date. Tagging 0.9.0 did not help, because the old workflow
  only deployed on pushes to `main`/`develop`, never on tags.

  Restored, with three changes. **A tag publishes `stable`**, not a push to
  `main`: the tag is the authoritative "released" signal, and hanging it on a
  `main` push is what let PyPI and the container ship while the site sat still.
  **No `paths` filter on push**, because GitHub requires the ref filter and the
  path filter to both match, so `tags: ['*']` beside `paths:` would drop any tag
  whose commits missed `docs/`. And **`workflow_dispatch` takes a version**, so
  0.9.0 can be backfilled without re-tagging.

  `mkdocs build --strict` returns as a PR check. Its absence is why three links
  to non-existent anchors survived into the release.

- **Git LFS had silently disabled the PHI guard.** `git lfs install` writes its
  hooks into `core.hooksPath` — `.githooks/` here — and only declines when it
  recognises the hook already present. Checking out a commit predating
  `.githooks/` removes the hooks, and the next LFS command installs a three-line
  `pre-push` stub in their place. That happened during the 0.9.0 release:
  `pre-commit` was gone and `pre-push` was the stub, so nothing scanned staged
  content for patient identifiers, and nothing reported it.

  `pre-push` now reads the ref list once and delegates to `git lfs pre-push`
  after its own checks pass, so LFS never needs to install anything and the two
  cannot clobber each other. The three LFS hooks are tracked for the same
  reason. A side effect: `git push` uploads LFS objects by itself again, rather
  than needing a separate `git lfs push --all`.

- **The version guard hardcoded `0.9.0`** in one of its two messages instead of
  formatting `MIN_PON_VERSION`, so raising the floor would have left that branch
  naming the old one.

- **Three documentation links pointed at anchors that never existed**, including
  one whose target is an `!!! note` admonition, which generates no anchor at all.

### Documentation

- **The PON version guard is in the README.** 0.9.0 refuses any PON built below
  the floor — the first thing an upgrading user with a self-built model meets —
  and neither the refusal nor its escape hatch was documented there.

- **`build-pon --from-outputs` is in the README and `pon-models.md`**, not just
  the developer guide. It turns a ~15h build into minutes and is how all four
  bundled models were rebuilt for 0.9.0. Both examples previously shown for it
  used a `--pon-variant` flag that `build-pon` does not have.

- **`BreakPointMotif`** added to the README feature table.

## [0.9.0] - 2026-08-11

### Added

- **The four GRCh37 PONs re-aggregated with the σ fix and the breakpoint
  blocks**, and stamped 0.9.0. Rebuilt via `--from-outputs` from the same
  `run-all` directories, so all four cohort digests are unchanged —
  `1e7158a60f3e3ca7`, `a0fed96df4c948ae`, `5c8c873ec1def225`,
  `655a5962b075cea2` — and 22 blocks become 24.

  The rebuild is provably a NaN-ing of residue and nothing else. Against the
  previous build: identical anchor sets, **bit-identical means**, **bit-identical
  σ at all 24.7 M surviving positions**, 1,177,647 positions newly NaN whose old
  σ topped out at 7.9 × 10⁻¹³, and **zero** positions newly finite.

  Residue σ falls to zero in three models and to 2 positions of 15.8 M in
  `xs2.all_unique` (σ 2.2 × 10⁻⁷ against a mean of 9.7 × 10⁻⁸ — the same
  phenomenon two decades above the 10⁻⁹ floor, left rather than chased with a
  reactive threshold change).

  `validate-pon` passes all four with **zero findings** — the first clean pass.
  On a real XS1 plasma sample, `BreakPointMotif` median |z| falls from 5.85 to
  **2.36**, max from 159 to 8.5, and motifs beyond |z| = 10 from 70 to **none**,
  putting it alongside `EndMotif` at 1.82 rather than 3× worse.

- **The product now records what produced it.** `{sample}.metadata.parquet` —
  the consumer's completion marker — gains five columns:
  `krewlyzer_version`, `pon_applied`, `pon_model`, `pon_cohort_digest` and
  `pon_krewlyzer_version`.

  The version *was* recorded, in `{sample}.features.json`, which invariant #2
  says downstream never reads. So the one fact needed to tell two runs apart
  lived precisely where the product ignores it. That mattered the moment this
  release changed what several columns mean: given a results directory there
  was no way to answer *was this produced before or after the FSD bin fix?*
  except by re-running it.

  `pon_model` is the **basename** only. The full path would put the operator's
  home directory into a shipped table.

  Stamped by `run_features`, so a single-feature CLI run records the same
  provenance as `run-all` (invariant #6), and stamped from `pon_for_zscore`
  rather than `pon` — a model loaded but withheld by `--skip-pon` scored
  nothing, and recording it as applied would be a lie about the columns present.

- **`validate-output` requires the PON-derived columns when, and only when, a
  PON was applied.** The contract declared none of them, so the seven WPS
  columns and FSD's `pon_stability` could vanish entirely and the gate would
  pass — which is exactly what happened twice in this release. They cannot be
  required unconditionally, because `--skip-pon`, no PON, and a PON the version
  guard refused are all legitimately unscored; `pon_applied` above is what
  tells the cases apart.

  A directory with no `pon_applied` at all was written by a build older than
  0.9.0, and is reported rather than assumed either way.

- **`no_collided_columns` in `validate-output`** — a `.1`-suffixed or repeated
  column name means a frame collision was written to disk, and no reader can
  tell which copy is current. Generic and universal: it runs on every table, so
  the next step that appends instead of replacing is caught without anyone
  thinking to look for it. It would have found the FSD duplication above
  unprompted.

- **`load_pon_model` refuses a PON older than 0.9.0.** Recording
  `krewlyzer_version` was only half a guard: `stamp-pon` writes it and
  `validate-pon` flags it missing, but nothing stopped a run scoring against a
  model built before the meaning changed — a fabricated `wps_background`, six
  floored σ, a region-MDS fitted over 65–400 bp while samples are measured over
  65–1000 bp. Every z-score against such a model is wrong in a way no schema
  check can see, because the columns are present and finite.

  Refusing, not warning. A warning in a log nobody reads is not a guard, and a
  plausible wrong number is the failure this release exists to remove. Set
  `KREWLYZER_ALLOW_OLD_PON=1` to override deliberately — without a documented
  way out, someone edits the parquet instead and the version stops meaning
  anything.

  Writing it exposed two holes it would otherwise not have covered:

  - `motif` and `region-mds` called `PonModel.load` directly, so those
    subcommands would have scored against a model `run-all` refuses —
    invariant #6. Both now use the guarded loader, and a test forbids any
    other caller.
  - The Rust scorers take a parquet *path* and read it themselves, and that
    path was set regardless of whether the Python load succeeded. A refused
    PON still reached FSD, WPS, OCF and region entropy: the Python half
    stopped and the Rust half carried on.

  The floor is a tuple, `(0, 9, 0)`, not a string — it is a compatibility
  floor, not the package version, and stays put at krewlyzer 1.5.0.

- **The 0.9.0 GRCh37 PONs.** All four rebuilt from `run-all` output via
  `--from-outputs`, and the first to pass `validate-pon` with **zero findings**:
  no fabricated baseline, no non-positive σ, every required block present, and
  a cohort digest identifying what each was built from.

  `test_the_gate_rejects_the_currently_shipped_models` has flipped to
  `test_the_shipped_models_pass_their_own_gate`. Its original docstring said
  that flip would be the acceptance record for the rebuild; this is it. Two
  further assertions join it: every shipped model records its provenance, and
  no two share a `wps_background` — the byte-identical baseline across four
  models is the defect that started this release.

- **`plausible_z_scores` in `validate-output`** — a z-score above 100 is a
  near-zero divisor, not biology. Every fabricated-σ defect in 0.9.0 produced
  *plausible* numbers and survived every schema check; this catches the
  opposite failure, where σ is so small the z explodes, which no schema notices
  either because the column is present, typed and finite.

  Runs on every table rather than a listed few, so a new output cannot opt out
  by omission, and its columns are read even where the contract does not
  declare them. Covers vector z columns element-wise — `wps_nuc_z` carries 200
  values per row, and a per-position σ is exactly where a near-zero divisor
  hides. For scale: real WPS z-scores against these models reach 6.7, with
  nothing above 10 across 6.9M positions.

- **`scripts/check_pon_env.py`** — behavioural probes that answer "does this
  install actually have the PON fixes?". `--version` reports the previous
  release until the bump, and a CLI flag only proves the *Python* side is
  current; neither can see a stale compiled extension, and `git pull` does not
  rebuild one. Each probe drives the real backend on input whose correct answer
  is known, and fails with the remedy. Run it before any PON build.

- **`validate-pon` checks a packing list.** It read the blocks present in the
  file and skipped the rest, so a block that vanished entirely was invisible to
  every check — indistinguishable from one never expected. Demonstrated on the
  real xs2 model: deleting `region_mds`, or `fsc_gene_baseline`, or almost
  every baseline at once, produced an identical finding list each time.

  That is the hole `region_mds` fell through when a Rust reader was handed a
  format it could not parse. New findings: `PON.BLOCK_MISSING` for a core block
  and `PON.BLOCK_MISSING_BAM_ONLY` for one that needs BAM input.

  `fsd_baseline` and `gc_bias` also join the σ checks. Their σ lives in per-bin
  columns rather than one `*_std`, which is why they were skipped — but both
  are applied, and a degenerate `gc_bias` reaches every FSC log-ratio.

- **`input_kind` in PON provenance** — `"bam"`, `"bed"`, `"mixed"` or
  `"outputs"`. `mds_baseline` and `region_mds` need a BAM, so their absence is
  legitimate for a fragment-BED cohort and a defect for a BAM one; without this
  the gate could only warn either way. Models built before it record nothing
  and still warn.

- **`scripts/build_pon_array.sh`** — one SLURM task per sample, feeding
  `build-pon --from-outputs`. Extraction across a cohort is embarrassingly
  parallel and the cluster already has a scheduler for it, so doing it four at
  a time inside one process leaves most of the machine idle:

  | route | 47 samples |
  |---|---:|
  | in-process, 4 parallel (`build_pon.sh`) | ~15.5 h |
  | in-process, `--parallel-samples 0` | ~6.2 h |
  | this array, 12 concurrent | ~5.9 h |
  | this array, unthrottled | ~1.5 h |

  Nextflow would fan out the same way, but every module pins
  `krewlyzer:0.8.3` and the SLURM profile enables singularity — a run today
  would silently execute 0.8.3 and reproduce the defects 0.9.0 removes. That
  needs a 0.9.0 container, which needs the release, which needs this PON. An
  sbatch array needs neither. `--from-outputs` will accept Nextflow's output
  unchanged once a container exists; it is the same layout.

  A finished sample is skipped on resubmission, keyed on the same last-written
  marker `--from-outputs` checks, so a partial directory is never mistaken for
  a complete one. `SAMPLE_LIST` is required rather than derived, as in
  `build_pon.sh`, and a task index past the end of the list is refused —
  a short array would otherwise build a quietly smaller cohort.

- **`krewlyzer build-pon --from-outputs DIR`** — build a PON by aggregating
  per-sample `run-all` output directories, reading only files. No BAM is
  opened.

  The in-process builder re-runs feature extraction itself, four samples at a
  time on one node: 55–97 minutes per sample, roughly 15 hours for a
  47-sample cohort. `run-all` already computes every one of those features, and
  its output directory is a strict superset of what the builder consumes. Every
  defect found in the 0.9.0 models so far has been in *aggregation*, not
  extraction — so each one cost a full rebuild to re-check. Re-aggregating now
  costs minutes.

  Both routes stay and converge on the same collectors, so neither can drift
  into producing a model the other cannot read. `SAMPLE_LIST` and
  `--from-outputs` are mutually exclusive: accepting both would leave it
  ambiguous which cohort the digest describes.

  A directory is aggregated only if it holds the *last* file its pipeline
  writes (`MDS.gene`) plus the inputs no baseline can do without. An incomplete
  directory is named and stops the build unless `--allow-failures` is given —
  silently aggregating a half-written one makes the cohort smaller than its own
  metadata claims.

  One difference is unavoidable: `EndMotif` stores frequencies rounded to six
  decimals, while the in-process path keeps them unrounded in memory. At a
  typical k-mer frequency of 0.0039 against a cross-donor σ of 2.7e-4, that is
  0.37% of σ — the two routes agree on σ to within 0.08% and on z to within
  about 0.002. Since Parquet is the product and the product is what rounds,
  that is the precision that actually exists.

  **A 0.8.3 output directory cannot seed the `wps_background` block.** It
  predates `nrl_at_band_limit`, and without that flag there is no way to tell a
  repeat length from the edge of the window it was searched in. The refusal
  says so. Every other block builds from 0.8.3 output.

- **`krewlyzer stamp-pon`** — records the release a built PON ships with.

  A PON is built from `develop`, where the version still reads the previous
  release, so the model records that however new the code is. Bumping before a
  four-hour build would fix it at the cost of putting a release number on
  unreleased code; this is the other order — build, then stamp when cutting the
  release. Added to the release guide as Phase 2.7, since the `sed`-based
  version update cannot reach inside a Parquet file.

  Afterwards `krewlyzer_version` means *the release this model is published
  with*, not the code that produced it — which is the definition a
  compatibility guard needs. `build_date` is untouched.

  **It refuses to stamp a model that fails `validate-pon`.** Without that it
  would be the shortest path to laundering: run it on one of the models
  carrying the fabricated `167.0 / 5.0` baseline and it would claim exactly the
  compatibility the guard exists to deny. `PON.NO_VERSION` alone does not
  block, since that is the condition being fixed.

- **`krewlyzer validate-pon`** — gates the reference, not the results.
  `validate-output` checks what a run produced; nothing checked the model
  those results are measured against. Every PON defect fixed in this release
  sat in four shipped files, visible in the file and invisible to every check.

  The load-bearing assertion is the same invariant the output gate enforces one
  level up: **a baseline that cannot vary with its cohort was not fitted to
  one.** A single check that σ differs between groups would have caught the
  fabricated `wps_background` in March.

  Also checks: σ positive and finite, `krewlyzer_version` and `cohort_digest`
  recorded, and every entry backed by ≥ 3 samples. Exit codes match
  `validate-output` — `0` satisfied, `1` violation, `2` structural.

  **All four currently-shipped PONs fail it**, which is the acceptance
  criterion: a gate that passes the models we know are wrong is not a gate. A
  test asserts that failure, and flipping it is the record that the rebuild
  worked.

  Wired into the Nextflow pipeline as `KREWLYZER_VALIDATE_PON`, running once
  per run rather than once per sample. Wired at the same time it was added —
  an unwired module is the defect #34 fixed, and it reads as coverage while
  providing none.

- **PON provenance: `krewlyzer_version`, `cohort_digest`, `cohort_label`.**
  The four models in this repository record `n_samples` and nothing else, so
  none can be reproduced, audited, or compared against a rebuild — and the
  build script that produced them has a stale comment claiming 47 where the
  metadata says 21, so even the integer is uncorroborated.

  The cohort is a **salted, non-reversible digest** of the sample IDs, stable
  across paths and ordering. Two builds from the same cohort match; the digest
  reveals nothing about who is in it. A PON ships in this repository, in the
  Docker image and on PyPI, so it is the last place an identifier may appear
  (invariant #4). `--cohort-label` adds a free-text name for humans.

  Writing is inert: `PonModel.load` reads through `meta.get(key, default)`, so
  unknown keys are ignored and existing models still load unchanged.

- **Per-motif z-scores (`EndMotif.frequency_z`, `BreakPointMotif.frequency_z`).**
  `mds_baseline` carries 625 k-mer means and σ — every 4-mer over ACGTN — and
  nothing read them; the tables shipped `Motif, Frequency` alone, so a shift in
  one motif was invisible unless it moved the whole-sample MDS enough to notice.

  The join is on the motif string, never position: 625 baseline keys against
  256 ACGT motifs in the output. And the sample is put on the baseline's scale
  first — sample frequencies sum to 1.0 across the 256, the baseline's
  expectations for those same 256 sum to **0.972**, with the missing 2.79% in
  the N-containing k-mers the output never reports. Comparing them directly
  biased every z upward: measured on a real sample, median **+0.37** naive
  against **−0.21** corrected.

- **WPS z-scoring — the largest PON baseline, previously read by nothing.**
  `wps_baseline` is ~128k anchors of 200-element mean and σ vectors, roughly
  90% of every PON file, and its only consumer was a log line appending
  `"WPS"` to a list of available components.

  `WPS.parquet` now carries `wps_nuc_z` / `wps_tf_z` (per-position z vectors)
  plus three derived shape quantities with their own z-scores, from a new
  `wps_shape_baseline` block.

  **Measurement drove every choice here, and ruled out the obvious ones:**

  | rejected | why, measured |
  |---|---|
  | mean of z over positions | lag-1 autocorrelation **0.986** — a fragment spans ~167 bp and touches many positions, so such a mean has nothing like σ/√200 precision |
  | max of \|z\| over positions | expected max under pure noise is **2.97**, so \|z\|>2 would flag nearly every anchor |
  | centre-minus-flank amplitude | TSS dips at the centre (−6.8 vs −3.4 flanks), CTCF does the reverse — any fixed window is backwards for one |
  | raw peak-to-trough range | correlates **+0.512** with `local_depth`; `log1p` drops it to −0.036 |
  | raw shape correlation | bounded at 1.0; mean 0.844 σ 0.099 means **302/400 anchors cannot reach +2** — Fisher `arctanh` removes the ceiling |

  **Displacement is measured but deliberately not z-scored.**
  `wps_phase_shift_bp` ships because it is cheap and genuinely non-redundant
  (corr −0.24 and −0.28 with the two scored statistics), but it gets no
  baseline: per-sample mean lag varies by **0.26 bp** against a within-sample
  spread of **8.43**, so there is no whole-sample phasing signal, and per
  anchor the intraclass correlation is **0.479** — about half of any lag is
  noise, optimistically, since that estimate used a baseline containing the
  samples being scored. It is also integer-valued, so on a small cohort its σ
  bottoms out at `std([0,0,0,0,0,1]) = 0.408` and a 1 bp shift scores z = 2.4.

  `wps_phase_at_search_limit` marks anchors where the ±30 search ended on its
  own edge (1.8% measured) — `nrl_at_band_limit` one level down.

- **Per-exon MDS baseline (`region_mds_exon`) and `MDS.exon.mds_z`.**
  `MDS.exon` is the finest localisation krewlyzer produces — 1,006 rows on
  xs2, 1,725 on xs1 — and was the only feature table with no baseline at all,
  so it shipped a raw score with nothing to compare it against.

  Keyed on `(gene, name)`: `name` alone is not unique across genes, and
  coordinates would break whenever the panel BED is regenerated. Aggregated in
  Python like the FSC gene and region baselines rather than adding a second
  Rust entry point for a groupby over ~1.7k rows.

  Measured before building it, which refuted the expected shape: exon data is
  **not** sparse. Every exon appears in every sample of its assay (7/7 xs1,
  19/19 xs2) with a measurable spread, and under 0.25% carry fewer than 10
  fragments. So no sparsity handling was needed — only the same NaN-not-floor
  rule as everywhere else. Verified on the real cohort: 1,006/1,006 xs2 exons
  fitted, median σ 0.0079, none unmeasurable.

  A PON built before 0.9.0 has no such block; `mds_z` is then `NaN` and a line
  says so.

- **`tests/real_data/` — a local cohort gate.** Point
  `KREWLYZER_TEST_CORPUS` at a scored cohort and it runs; leave it unset and it
  skips, so `pytest tests/` stays green and CI never sees it. No cohort data is
  committed, and no patient identifier may appear in a test name, parameter id,
  assertion message or artifact — `sample_label()` gives a non-reversible
  handle where output needs one.

  It exists because every defect the PON audit found needed a cohort to see: a
  single fixture cannot reveal a constant, cannot match a real PON arm, and
  cannot show per-anchor sample support. The unit suite proves a scorer *can*
  run; only a cohort proves it *did*.

  The four PON blocks that are built and consumed by nothing are marked
  `xfail(strict=True)`, so Phase B wiring them turns the test red and forces
  the marker off.

- **`read_exact_table`** in `core/output_utils.py` — reads the given path with
  no sibling resolution, the counterpart to the deliberately parquet-first
  `read_table`.

- **`tests/invariants/test_read_after_write.py`** — pins both readers'
  behaviour and an AST check that no function reads with the resolving reader
  after writing, with a reviewed allowlist that must state why each entry is
  safe. Also fails on a stale allowlist entry, so it cannot rot.

- **`tests/invariants/test_pon_format_parity.py`** — asserts every PON-derived
  column is present, populated and equal across `tsv` / `parquet` / `both`.
  One test that would have caught FSD, FSC gene, FSC region and the entropy
  processor together. `both` is covered explicitly: it leaves a TSV *and* a
  parquet on disk, which is what let the parquet-first reader pick the wrong
  one.

- **FSD was never PON-normalised under `--output-format parquet`.** Under
  `tsv` the bins came out as log2 ratios; under `parquet`, raw counts — same
  column names, no warning either way. **0.9.0 makes parquet the Nextflow
  default**, which would have turned this from affecting nobody into affecting
  every pipeline run.

  Two independent bugs, each sufficient on its own:

  1. The caller guarded on `outputs.fsd.exists()`, where `outputs.fsd` names
     the `.tsv`. The Rust writer honours `--output-format`, so under parquet no
     `.tsv` was ever written and the whole post-processing block — PON
     normalisation included — was skipped. Same class as `b77909e`, missed
     site; now uses `resolve_table_path`.
  2. `_write_fsd_output` read back through `read_table`, which is
     **parquet-first**: given `x.FSD.tsv` it prefers `x.FSD.parquet` — the raw
     table the single-pass writer emitted *before* normalisation. So the
     log-ratios were computed, logged as `41 arms normalized`, and then
     overwritten with the raw counts.

  `FSD PON: no arms matched` is now a warning rather than a debug line: zero
  arms means the output is raw counts wearing the same column names as
  log-ratios, and nothing downstream can tell the difference.

- **The WPS background PON baseline was never fitted.** `pon/build.py` looked
  for columns named `nrl`/`nucleosome_repeat_length` and
  `periodicity`/`period_score`; `WPS_background` writes `nrl_bp` and
  `periodicity_score`. Neither ever matched, so the fallback fired on every
  build and all four shipped PONs carry an identical hardcoded
  `167.0 / 5.0 / 0.0 / 1.0` across all 28 groups — from cohorts of 21 and 47
  samples, logging `28 groups` as though it had fitted them. Combined with the
  0.8.x `nrl_bp ≡ 150.0` degeneracy, **every sample ever produced scored
  `nrl_z = -3.4`**: a constant, moderately extreme z presented as a
  measurement. A missing source column is now fatal.

- **Six σ floors replaced with NaN.** `max(std, 0.001)`, `max(nrl_std, 0.1)`,
  `.max(0.01)` and three siblings do not make a z-score conservative — they
  make it arbitrarily large, because the divisor is a number nothing measured.
  That is `nrl_at_band_limit` one level down: a boundary value reported as a
  result. NaN propagates to an absent z, which is what "we could not measure
  this" should look like. A single-sample WPS anchor previously got σ = 0.1 at
  every position.

- **The WPS baseline had no minimum-sample requirement.** FSC gene and FSC
  region have required ≥3 since they were written; WPS — 141k anchors, the
  largest block by 100× — required nothing, so an anchor seen in one sample
  still produced a baseline. Measured on the shipped models that was 1.6% of
  anchors for all_unique/xs1 and 28.8% for duplex/xs1. Now ≥3, and the count
  skipped is logged.

- **`np.std` → `ddof=1` in the PON baselines.** These are a sample of healthy
  donors, not the population. The population form understates the spread and
  inflates every z built from it — by 2.5% at n=21.

- **`build-pon --keep-sample-outputs DIR`.** Per-sample feature outputs were
  extracted to a temp directory and deleted on completion, so every rebuild
  re-ran extraction over every BAM from scratch (~4 h for 47 samples) and
  leave-one-out calibration would have cost *n* of those. Kept outputs also
  survive a failed build, so the rerun can skip what finished.

- **`_log_baseline_quality`** reports how many entries in a block carry a
  measured spread, warns on those that do not, and warns when every entry
  shares one σ — the signature of a fabricated baseline.

- **`docs/cli/index.md` covers every command**, and two tests keep it that way:
  one cross-checks the registered commands against the reference and the
  README, the other resolves every README documentation link against the pages
  mkdocs actually publishes.

- **`krewlyzer report SAMPLE_DIR -o report.html`** — a single-sample report
  with a cross-axis verdict, 16 charts and the tables behind them. Internal
  use: it contains one patient's measurements, so it is generated on demand
  rather than committed or published. `describe-output` remains the shareable,
  structural view.

  **Verdict.** Four independent axes — fragment size, nuclease cutting, tissue
  shedding, chromatin accessibility — each showing its value, source column and
  direction, summarised as *"N of M assessable axes agree"*. Never a composite
  score: one would hide exactly the disagreement worth seeing. An axis with no
  PON z-score is *not assessable*, never counted as disagreement, because that
  would make a thinner run read as a healthier one.

  **Organised by table, not by chart.** Each section carries one output's
  chart, its **why / how / what**, and its columns together. A chart three
  screens from its data is a chart nobody connects to the data.

  **PON state is first-class.** Without one, every z-score and log-ratio is
  absent — most of the interpretable surface — so the report leads with a
  banner and marks each affected column.

  Sixteen charts across the five-act structure from the output EDA notebooks,
  including three taken from the existing notebook's specialised panels: the
  per-region FSD heatmap, which shows whether the size shift is uniform across
  chromosome arms or concentrated on a few — something the summed density curve
  cannot; the WPS foreground profile around TSS and CTCF anchors; and the GC
  correction curve, since everything GC-corrected downstream is multiplied by
  it.

  Plotly, self-contained (~5.6 MB, no network). The theme follows the docs site,
  with an Auto/Light/Dark toggle that the figures respect. `plotly` is an
  optional `[report]` extra; without it the command still runs and each chart
  states that it is missing.

  Three rules the charts hold to, each learned from a real defect here: a
  constant column is reported as constant rather than drawn as a flat line; no
  threshold is drawn as a cut-off, only labelled literature anchors; and a
  value that means "nothing was observed" is never placed on the measurement
  scale — an mFSD variant with no ALT fragment has `ALT_MeanSize = 0`, and
  plotting `0 − REF` would put a fabricated point at the most tumour-like end
  of the very axis being read. Those variants appear at the origin with their
  own marker, because their absence is itself a finding.

- **`krewlyzer describe-output SAMPLE_DIR`** — says what a sample's output
  files contain, which is the question people ask before "is it correct?".

  Row and column counts, dtypes, value ranges, distinct counts and null counts
  per column, plus total size and which tables are read downstream. Markdown by
  default; `-o report.html` renders a self-contained page.

  Everything is either measured from the file or read from `contract.py`.
  Nothing about the biology is restated — that stays in
  `docs/reference/output-files.md` and is linked — so a table that gains a
  column or stops being consumed changes the report automatically and the two
  cannot disagree.

  Each table carries a one-line **what it measures** and **which way cancer
  moves it**, from a new `validate/meaning.py`. Direction is the point: it
  differs per axis, and getting one backwards is the commonest misreading of
  this output — `MDS` was documented the opposite of its own threshold table
  for a year precisely because no test and no report ever stated one. A test
  requires a direction for every table except `metadata`, which is provenance.

  No thresholds are recorded. Every numeric band examined turned out to be a
  display default or refuted outright: the documented ATAC/TFBS entropy range
  flags a perfectly healthy N(167,30) distribution as abnormal. Directions are
  robust, magnitudes are cohort-specific, and a number in the registry would
  acquire an authority it has not earned — a test enforces their absence.

  **Identifier columns are redacted.** Sample directories here are named for
  the patient and several tables carry the sample id as a *column value*, so
  the first version of this report leaked a real identifier through `Sample`
  and `sample_id` examples even though it was generated from deliberately
  renamed files — the filenames were clean and the contents were not. Knowing a
  column holds an identifier is the useful fact; which one is not, and this
  report is meant to be hosted.

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

### Changed

- **WPS PON scoring moved from Python to Rust** (`rust/src/wps.rs::apply_pon_zscore`).
  It was the last PON z-score still computed in Python, and it was the largest:
  89k anchors, a ±30 lag search each, about 5.4M correlations.
  `.agents/rules/architecture.md` puts "PON z-score", "loops over >1000 rows"
  and "row-level computation" on the Rust side; this was all three.

  Measured on a real 76,595-anchor output against a shipped PON: **9.5 s
  against ~17 min**, roughly 107×.

  Written in Python first on purpose, and that is not an apology — three of its
  decisions changed under measurement and would not have survived being written
  in the faster language first: the `log1p` amplitude (raw range correlates
  +0.512 with depth), the Fisher transform (bounded *r* left 302 of 400 anchors
  unable to reach +2), and the deliberate absence of a phase-shift baseline
  (intraclass correlation 0.479). The port is **bug-for-bug**, including a tie
  in the lag search resolving to the most negative lag — measured incidence 0
  in 317 anchors, and correcting it would have destroyed the oracle's ability
  to tell a porting slip from an intended change.

  The Python reference is frozen as that oracle in
  `tests/unit/test_rust_python_equivalence.py`. Nothing imports it outside
  tests; the shipping module is now a ~40-line call. Tolerance was fixed at
  `1e-6` relative *before* the first comparison; measured worst case **2.2e-13**
  on synthetic anchors and **9.0e-14** against real data.

  `apply_wps_pon` now takes the PON **path** rather than a loaded model, and
  `baseline_attr` is renamed `baseline_table` — the argument names a table in
  the parquet, not an attribute on an object.

- **Two silent-read traps found while porting**, both of the "present,
  plausible, returns nothing" class that this release is mostly about:

  Every shipped PON stores `table` and `region_id` as Arrow `large_string`,
  never `string`. A reader that downcasts to `StringArray` gets `None` and
  yields an *empty* baseline — a legitimate state that degrades silently
  instead of raising. `pon_model.rs::PonModel::load` does exactly that and is
  blind to every PON this project has ever shipped. It survived because it has
  **no callers**: each reader loads what it needs, and `pub` items are exempt
  from dead-code warnings. Recorded in `architecture.md` rather than deleted,
  pending an explicit decision on scope.

  Vector columns are `list<double>` from the Python builder and `list<float>`
  from the Rust one. The new reader accepts both and logs the type it could not
  read, rather than reporting no anchors.

- **Four `eprintln!("DEBUG: …")` calls in `wps.rs` now go through the logger.**
  They wrote to stderr on every WPS run, below no level and past every filter.

- **The two PON blocks that were built, stored, and read by nothing are now
  applied.** Found by tracing all 21 baseline blocks from build → parquet →
  read → apply; 19 were reaching scoring code and two were not.

  `wps_baseline_panel` (326 anchors in the xs2 model) — `apply_wps_pon` ran
  only on the genome-wide WPS, so `{sample}.WPS.panel.parquet` shipped raw.
  The anchors closest to the targeted regions were the one WPS output with no
  comparison to a healthy cohort. They now carry the same derived and z-scored
  columns as the genome-wide file. The shape statistics borrow the genome-wide
  `wps_shape_baseline`: a few hundred panel anchors is too few to fit a second
  one, and they overlap the genome-wide set.

  `gc_bias_ontarget` (25 bins) — `PonModel.get_mean` and `get_variance` read
  `self.gc_bias` only, so on-target FSC and FSR normalised against the
  genome-wide GC curves. Capture enrichment shifts the GC profile, which is the
  entire reason panel mode fits a second block. `{sample}.FSC.ontarget`
  `*_log2` and `*_reliability` values change as a result. A PON without the
  block falls back to the genome-wide curves rather than dropping the column.

- **A scored column is NaN when its baseline is absent, not 0.0.** Six
  fabrications in the apply path, after `zscore_or_nan` removed nine from the
  baseline classes:

  | where | when | was | now |
  |---|---|---|---|
  | `fsc_processor` | no `gc_bias` | `{ch}_log2 = 0.0` | NaN |
  | `fsc_processor` | no variance | `{ch}_reliability = 1.0` | NaN |
  | `fsr_processor` | no long fragments | `ratio = short_count` | NaN |
  | `fsr_processor` | ratio ≤ 0 | `short_long_log2 = 0.0` | NaN |
  | `region_entropy.rs` | σ unmeasurable | `z = 0.0` | NaN |
  | `region_entropy.rs` | σ column absent | σ ← 1.0, mean ← 0.0 | NaN |

  Zero is never the cautious choice. A log2 ratio or z-score of zero says
  "this sample sits exactly at the healthy baseline" — the most confident
  statement either column can make, asserted precisely when nothing could be
  compared. Measured: with `gc_bias` absent, three windows whose raw values
  differed fourfold all read `0.0`.

  The Rust pair defeated the Python fix until both were done. `entropy_std`
  defaulted to `1.0` and `entropy_mean` to `0.0` — a standard normal — so a
  builder storing NaN for a single-donor label produced the raw difference
  rather than no score. And `position(...).unwrap_or(0)` read *column zero*
  when a column was not found, which in a PON parquet is `table`.

  `region_entropy_processor` was the fifth site of the `sample_std_or_nan`
  defect and the only one outside `pon/build.py`.

- **`scripts/build_pon_unfiltered.sh` → `scripts/build_pon.sh`**, taking assay
  and variant as arguments. The old script was hardcoded to xs2 and carried a
  header claiming "47 samples" inherited from the xs1 copy it was made from —
  while the xs2 model it produced records 21. Four models built from four
  edited copies of one script is how that happens.

  It also passes `--keep-sample-outputs` and `--cohort-label`, and **ends by
  running `validate-pon` on what it built**. A build that produces a model the
  gate rejects has not succeeded, whatever `build-pon`'s exit code said.

  `docs/guides/building-pon.md` gains the rebuild runbook: what to watch in the
  log, why the sample outputs are kept, and the LFS ordering that has to be
  right.

- **`MIN_SAMPLES_PER_KEY`** names the ≥3 floor once in `pon/build.py`. FSC gene
  and FSC region have required it since they were written and WPS acquired it
  in `4cd634b`; the four now agree by construction rather than by three
  separate literals happening to match.

- **OCF PON logging distinguishes "absent from the baseline" from "matched but
  unscoreable".** Both correctly write NaN, but collapsing them reported
  `0/1 regions matched` for a region that matched perfectly well and simply has
  no variance in the cohort. That is the difference between "rebuild the PON"
  and "expected", and only the log can say which.

- **`tests/integration/test_pon.py`'s PON fixture used a schema that has never
  existed** — `version` / `genome` / `sample_count` where real PONs write
  `schema_version` / `assay` / `n_samples`. Since `PonModel.load` reads through
  `meta.get(key, default)`, it loaded as a completely empty model and the tests
  asserted only `is not None`. `test_pon_model_loading` had been marked
  `skip("PON parquet schema requires production format")` — true of the
  fixture, not the loader. Fixture corrected against the bundled PON, test
  unskipped, and both now assert real field values.

- **PON scoring was skipped under `--output-format parquet` in three more
  places.** Same shape as the FSD defect below: callers name a hardcoded
  `.tsv`, the Rust writer honours `--output-format`, and a bare `.exists()` is
  then False.

  - `apply_fsc_gene_pon` returned `None` with "file not found" and added no
    `depth_zscore`
  - `apply_fsc_region_pon` and the e1-only filter, identically
  - `process_region_entropy` skipped the whole step

  All four now resolve with `resolve_table_path`. Measured before and after on
  the bundled PON: `depth_zscore` absent → present with 5/5 populated.

- **TFBS and ATAC emitted `z_score = 0.0` when there was no baseline.** Zero is
  not a neutral placeholder — it is the most confident claim the column can
  make ("this sample sits exactly at the healthy baseline"), asserted on the
  strength of having no baseline at all, and indistinguishable from a measured
  zero. Now `NaN`, with a WARNING naming how many labels are affected. Applies
  both when no `--pon-model` was given and when the PON lacks the table.

### Removed

- **271 lines of dead Rust.** `_core.wps.apply_pon_zscore` (173 lines) plus
  `load_wps_baseline_from_parquet`, `phase_shift`, `PHASE_MAX_LAG` and the
  parquet imports that existed only to serve them. The crate now builds with
  zero warnings.

  `_core.wps.apply_pon_zscore` in detail: It emitted one scalar z
  per anchor from `wps_long_std`/`wps_short_std` — v1.0 baseline field names,
  while every shipped PON is v2.0 vector format — and pushed `0.0` where σ was
  not positive. It was also unreachable: the call was gated on
  `pon._source_path`, an attribute nothing in the codebase sets. Dead,
  schema-obsolete, fabricating, and implementing the statistic measurement
  refuted.

- **`PonModel.save`** — a second serializer that wrote only the metadata block,
  producing a PON with no baselines at all, while `build-pon` used
  `_save_pon_model`. Nothing in production called it. Removed rather than
  completed: a second writer is a second thing to keep in step with every new
  block.

### Fixed

- **`pon_stability` was an unnormalised inverse variance and collapsed to zero.**
  `1 / (mean(σ²) + 0.01)` only ever looked sane because it was fed the wrong σ —
  `get_stats` returned the first size bin's spread for every size. Once each bin
  got its own, the true magnitudes showed: FSD σ is in fragment counts,
  111–4,573 on xs1 and 947–35,961 on xs2, so `mean(σ²)` reached 1.6×10⁶ and
  8.3×10⁷ and the reciprocal wrote `0.000000` at six decimals for **every arm**.
  The `0.01` floor assumed σ ≈ O(1); it never was. Caught by the
  anti-degeneracy check on a real XS2 sample — invariant #1 working.

  Now `1 / (1 + mean(CV²))` with CV = σ/expected, bounded in **(0, 1]**.
  Dividing by `expected` is justified, not convenient: measured on both models,
  `d log σ / d log expected` is **0.83** and **0.95** — close to 1, so the
  spread is between-donor variation rather than Poisson counting noise (which
  would give 0.5), and the ratio is genuinely scale-free.

  On real output: 41 distinct values per assay at six decimals, against 5 and 1
  before. xs1 spans 0.714–0.758, xs2 0.865–0.926 — and that difference is
  interpretable, the 21-donor xs2 cohort agreeing more tightly than the
  47-donor xs1 one. The between-arm spread is only ~6% but it is real rather
  than estimation noise: split-half over odd/even bins gives r = 0.834 and
  0.997, a full-length reliability of **0.91** and **0.999**.

  **Old and new values are not comparable** — the scale changed from an
  unbounded inverse variance to a bounded agreement score.

- **`pon_stability` was declared `vary=BOTH`, which it cannot satisfy.** It is
  computed from the PON's σ alone, with no sample value entering it, so every
  sample scored against one model gets the same 41 numbers. Now `vary=WITHIN`,
  which still requires it to differ down the rows — the check that caught the
  collapse above.

- **The gate could not see four of the six tables Parquet mode dropped.**
  `.OCF.parquet` and all three `.sync` tables were absent from the output
  contract; only the two on/off-target summaries were declared, and they are
  what reported the loss. All four are declared now, with `ocf_z` marked
  `requires_pon` — the `.sync` tables carry the per-position profiles the
  summaries are reduced from, which is the large half of OCF.

- **`FSC.regions.is_e1` / `is_alt_e1` were declared `string` and written
  `int64`.** `gene_bed.py` parses `fields[8] == "1"` into a Python bool, which
  lands as int64 0/1 — so the contract was wrong, not the writer. The
  synthetic cohort wrote strings too, so fixture and contract agreed with each
  other and neither with the code, which is how it survived. These columns are
  new in 0.9.0, so nothing older disagrees.

- **`fsc_gene_ratios_sum_to_one` flagged regions that have no fragments.** All
  six channels are zero, so the sum is zero and the row reads as a maximal
  partition failure — "worst deviation 1.000000" — when nothing is wrong with
  it. 3 of 1725 regions on a real XS1 plasma sample, 13 of 1725 on a shallow
  one, making the count a function of sequencing depth and burying any real
  break among them. Zero-fragment regions are now excluded; a region with even
  one counted fragment must still partition it.

- **`BreakPointMotif` was scored against the `EndMotif` baseline** — and so
  were the two on-target motif tables. `apply_motif_pon` used
  `pon.mds_baseline` for all four, an **end-motif, genome-wide** block.

  An end motif is the 4-mer at the fragment's 5′ terminus, what the nuclease
  left. A breakpoint motif spans the cut site and includes reference bases
  *not present in the fragment*. The two frequency vectors correlate only
  **0.696** on the same sample, so `frequency_z` there measured the offset
  between two definitions rather than any departure from healthy:

  | | `EndMotif` median / >10 | `BreakPointMotif` median / >10 |
  |---|---|---|
  | XS1 | 1.82 / 0 of 256 | 5.85 / **70** of 256 |
  | XS2 | 4.47 / 39 of 256 | 11.25 / **136** of 256 |

  A correctly fitted baseline gives a median near 0.67 — the half-normal
  median. Demonstrated on synthetic cohorts with the real correlation
  structure: median |z| 6.75 with 86 of 256 beyond 10 against the end-motif
  baseline, **0.75 with none beyond 10** against its own.

  The PON gains `breakpoint_motif_baseline` and
  `breakpoint_motif_baseline_ontarget`, built by both routes. `mds_mean` and
  `mds_std` are NaN in those blocks: MDS is defined on end motifs, and a number
  there would be a different statistic under the same name.

  `mds_baseline_ontarget` needed no new work — it already existed and was
  already used for the whole-sample MDS z, just not for the per-motif
  frequencies. A PON without the breakpoint blocks (every one built before
  0.9.0) now yields no `frequency_z` on those tables rather than a wrong one.

- **A σ of 10⁻¹⁷ was used as a divisor.** Every PON built to date carries
  positions whose baseline spread is floating-point residue: 4.6 % of positive
  σ in xs1.all_unique, 12.0 % in xs2.all_unique, 47.2 % and 55.4 % in the
  duplex pair. On a real XS2 plasma sample that produced **728,007**
  `wps_nuc_z` above 100 and **354,260** above 10⁶, peaking at 6.1 × 10¹⁸.

  WPS is (fragments spanning a position) − (fragments ending near it), and each
  fragment carries a fractional GC-correction weight. Where no donor had
  coverage, or the two terms cancel, the true value is zero — but it computes
  as ~10⁻¹⁷ rather than `0.0`, so the exact `min == max` test added earlier in
  this release never fired. `wps_tf` corroborates it: 0 % residue and 185,079
  honest NaN, because TF-sized fragments are rare enough that those positions
  sum to exact zero.

  `SIGMA_FLOOR = 1e-9` now applies in the builder (`element_wise_std`,
  `mean_and_sd`) *and* in the scorers, so the four models already in the wild
  are safe without a rebuild. Not a taste threshold: across the 24.7 M positive
  σ of the shipped xs1.all_unique baseline, 1,177,647 sit below 10⁻¹², **none**
  sit between 10⁻¹² and 10⁻⁶, and the rest are measurements — 10⁻⁹ is the
  log-space midpoint of six empty decades. Mirrored in
  `pon/model.py::SIGMA_FLOOR` and asserted equal in `validate/claims.py`.

  Measured on a real XS1 plasma sample against the **unchanged** shipped PON:
  max |z| falls from 2.8 × 10¹⁷ to 3.7 × 10⁴ and values above 10⁶ from 18,398
  to **zero**.

  `scripts/check_pon_env.py` gains a matching probe. The existing one asserts
  that *identical* donor vectors give NaN and passed throughout, because the
  real case is *near*-identical; the new one was confirmed to fail against the
  pre-fix binary before it passed against the fixed one.

- **`--output-format parquet` silently discarded every OCF table.** All six —
  `OCF`, `OCF.sync`, `OCF.{on,off}target`, `OCF.{on,off}target.sync`. Measured
  on a real XS1 plasma BAM: 44 logical tables under `both`, **38** under
  `parquet`, exactly the OCF six missing.

  The OCF temp directory is an internal intermediate — Python moves the files
  to their final names and converts them, because `_core.ocf.apply_pon_zscore`
  reads TSV line by line. `pipeline.rs` nonetheless dispatched that
  intermediate write on the user's `--output-format`, so Parquet mode wrote
  `all.ocf.parquet` while the mover looked for `all.ocf.tsv`. Nothing matched,
  nothing moved, the temp directory was deleted. No warning, because an absent
  file is indistinguishable from "OCF was not requested".

  The FSC counts three lines above already carried the right rule — *"an
  internal intermediate; keep as TSV"*. OCF follows it now, and the mover
  reports an empty intermediate as an error. The `.sync` tables are the large
  ones (483 KB TSV vs 56 KB Parquet), so this removed orientation-aware
  fragmentation from precisely the format the release ships in (invariant #2).

- **On-target FSD was never PON-scored under Parquet.** The same defect one
  branch below its own fix: `.exists()` on a hardcoded `.FSD.ontarget.tsv`,
  fifteen lines under the comment explaining why the genome-wide branch had
  stopped doing that. Confirmed on real XS1 and XS2 plasma samples — 69 columns
  and zero `_logR` against the genome-wide table's 137. Resolved with
  `resolve_table_path`; an absent table is now warned about.

- **On/off-target OCF left its intermediate `.tsv` beside the Parquet.**
  `apply_ocf_python_pon` wrote the final format but never called
  `cleanup_intermediate_tsv`, which the genome-wide path does — so a Parquet
  run shipped two files for one table with no way to tell which was current.

- **`fsd_only_size_bins` fired on every PON-scored FSD table.** The check's
  reserved column set was written for the raw table and never updated when PON
  scoring began writing into the same file, so it reported all 67 `{bin}_logR`
  columns plus `pon_stability` as stray — 68 spurious errors on every scored
  sample. A gate that fires on correct output is a gate nobody reads. It went
  unnoticed because the synthetic cohort carried no PON columns until now, so
  the check had never once seen the shape it runs against in production.

- **The WPS column reference documented five columns that do not exist, missed
  eight that do, and typed the profiles as scalars.** Found by diffing the
  reference against a real scored table instead of reading it.

  `wps_nuc_smooth`, `wps_tf_smooth`, `wps_nuc_mean`, `wps_tf_mean` and
  `wps_tf_z` were never written. `capture_mask`, `local_depth` and the six
  derived PON columns were. And `wps_nuc`, `wps_tf`, `prot_frac_nuc`,
  `prot_frac_tf` and `wps_nuc_z` are 200-element **lists**, documented as
  `float` — so the three worked examples indexed columns that do not exist and
  raised `KeyError`. All three are rewritten and were run against a real
  76,595-anchor table before being committed.

  A sweep now checks every `df["col"]` reference in every Python block of
  `output-files.md` against real schemas from a full run: 0 of 40 table
  sections reference a column that does not exist.

- **FSD's histogram range was documented as the filter range.** `meaning.py`
  said "5 bp bins over 65–1000 bp"; the histogram is 67 bins over `[65, 400)`.
  Fragments of 400 bp and above are excluded from the bins **and** from
  `total`, so a long-fragment fraction computed against `total` had a
  denominator that silently excluded them. Corrected, with the range pinned in
  `validate/claims.py` so the two cannot drift again (invariant #5).

- **The FSD PON columns were documented nowhere.** `{bin}_logR` and
  `pon_stability` — the entire product of PON scoring for FSD — appear in
  neither `output-files.md` nor `meaning.py`. Both are now described, including
  what `pon_stability` means and why it is NaN rather than 0.990 when no σ
  could be measured.

- **WPS z-scores were written to a file nothing reads, then destroyed their own
  input when that was fixed.** Two defects stacked, both introduced by the Rust
  port earlier in this release and both caught by reviewing the *product*
  rather than the code.

  `run-all` passes `{sample}.WPS` as the output base — a stem carrying a
  compound extension. `Path("s.WPS").with_suffix(".parquet")` reads `.WPS` as
  the suffix and replaces it, so the scored 18-column table landed on
  `{sample}.parquet` while `{sample}.WPS.parquet` stayed the raw 11-column
  profile. Downstream reads Parquet only (invariant #2) — the same failure the
  `--output-format` fix removed earlier in this release, in a second disguise,
  under a comment that asserted the correct behaviour while the code did the
  opposite.

  Fixing the path exposed the second: the normal case *is* in place, and the
  scorer streams its input in batches while `File::create` truncated that same
  file out from under the reader. It now writes a sibling temp and renames
  atomically, so a crash leaves the original intact rather than a half-written
  product, and a return of 0 still writes nothing.

  Every earlier test used a simple stem (`tmp_path / "o"`), which has no
  compound extension and so could not see either. The new tests use the shape
  the caller actually passes.

- **The four shipped PONs are stamped `0.9.0`.** They were built from
  `develop`, where the version still read `0.8.3`, so every model recorded the
  previous release however new the code was. `krewlyzer_version` now means the
  release a model is *published* with — which is what the compatibility guard
  needs — while `build_date` still records when it was built. Only the metadata
  row changed: all four cohort digests, block counts and row counts are
  unchanged, and `validate-pon` passes on the stamped files, which are the
  files that ship.

- **`validate-pon` reported `0 sample(s)` after checking four models.** The
  summary line is shared with `validate-output`, which counts samples; the PON
  gate fills no samples, so it printed a line saying nothing was inspected
  directly above one saying the contract was satisfied. It now reports what it
  actually checked: `4 model(s)`.

- **Every FSD log-ratio was scored against the wrong size bin.** The largest
  output defect found in this release, and it reached every sample ever
  normalised against a PON.

  `size_bin` is an integer in every sense that matters, and it is stored as a
  **double**: a PON is one long-format table, so every row belonging to another
  block carries a null there and the union column comes out `float64`. The
  reader called `row.get_int()`, which *errors* on a Double rather than
  coercing, and the error fell through to `unwrap_or(0)`.

  So every `size_bin` became 0. An arm's `size_bins` was 67 zeros,
  `size >= *size_bins.last()` held for every size, and `get_expected` returned
  the **last row's** expectation for all 67 bins.

  Measured on a shipped PON and a real sample: **41/41 arms** matched the
  last-bin baseline exactly and the bin-matched baseline not at all. After the
  fix, 41/41 match bin-matched and 0/41 match the old behaviour. The correction
  is large — median log2 shift −1.05, maximum |Δ| 5.10.

  Nothing could have caught it downstream. The log-ratios still varied across
  bins, because the sample numerator varies, so they were present, finite,
  non-degenerate and wrong — invariant #1's failure mode one level subtler than
  the constant it was written for.

  **The PONs are not affected and do not need rebuilding**; the reader was
  broken, not the model. Samples normalised against a PON do need re-running.

- **`pon_stability` used one size bin's sigma for the whole arm.**
  `FsdBaseline::get_stats` returned `std.first()` regardless of the size asked
  for. Measured on the shipped xs1 all-unique PON: sigma varies **41.6×**
  (median) across an arm's 67 bins, up to 56.4×, so `pon_stability` was wrong
  by a median of **4709%** on all 41 arms. Sigma is now interpolated at the
  requested size, exactly as `expected` already was.

- **The same column-lookup defect closed in `ocf.rs`, the third and last copy.**
  `position(...).unwrap_or(0)` made an absent column read column zero, and
  `ocf_mean`/`ocf_std` defaulted to `0.0`/`1.0` — together making z the raw
  difference from zero. OCF's values parse correctly today (verified against a
  shipped PON: an observation 2σ above the mean scores exactly 2.0), so this
  closes a landmine rather than fixing a wrong number. It would have fired the
  moment a column was renamed. The `table` lookups in `fsd.rs` and
  `region_entropy.rs` are by name now too.

  Two `unwrap_or(0)` lookups remain, both on **TSV headers** rather than PON
  columns, where column zero genuinely is the label in most files. Left alone
  deliberately.

- **Three more fabricated defaults on the FSD baseline path**, the same family
  removed from nine baseline classes as `zscore_or_nan` this release: an
  unreadable sigma defaulted to `1.0` (which with the 0.01 floor reads as a
  `pon_stability` of 0.990, indistinguishable from a measurement), an
  unreadable expectation to `0.0`, and `position(...).unwrap_or(0)` made a
  column that could not be found read **column zero** — `table` — the same
  defect already removed from `region_entropy.rs` and missed here. Rows with an
  unreadable `size_bin` or `arm` are now skipped and counted in a warning
  rather than silently collapsing onto bin 0 or an empty key.


- **WPS z-scores were written to a file nobody reads.** `apply_wps_pon`
  honoured `--output-format`, whose default is `tsv`, while the Rust step
  writes `.WPS.parquet` unconditionally. A default run therefore produced two
  files: a `.WPS.tsv` carrying all seven PON columns and a `.WPS.parquet`
  carrying none. Downstream reads Parquet only (invariant #2), so WPS scoring
  was present on disk and absent from the product.

  Measured on a real 89,034-anchor run: the TSV had 18 columns at 928 MB, the
  parquet 11 at 144 MB, and the parquet was the older file. WPS is Parquet-only
  by contract precisely because of that size ratio — 200-element vectors are
  ~6× larger as text.

  The writer no longer takes an `output_format` at all, so the knob cannot be
  set wrong; when a run requests another format, the WPS exception is logged
  rather than silently applied. Introduced in `8265ad3`, and missed because
  every unit test passed `output_format="parquet"` explicitly.

- **FSD log-ratios accumulated on every re-run.** `apply_pon_logratio` found its
  size bins by taking any header whose text before the first `-` parsed as an
  integer — so `65-69_logR`, a column it had written itself, came back as bin 65
  on the next pass. Each run then appended a fresh set of log-ratios *and*
  log-ratios of the previous ones.

  Measured on a real 67-bin sample: **69 columns raw, 137 after one pass, 273
  after two**, every `_logR` name repeated. `_write_fsd_output` reads that back
  with pandas, which renames the collisions to `_logR.1` / `_logR.2` and writes
  those to disk. The correct answer was in the file three times with no way to
  say which was current — and `read_csv` returns the first, which is the oldest.

  Nothing marked the work as done, so a rerun into the same directory, a
  Nextflow retry or a resumed job all triggered it. The normaliser now drops the
  columns it previously produced and rebuilds them, so running it twice gives
  byte-identical output; a size bin is recognised as exactly `{int}-{int}`; and
  it logs when it is replacing rather than appending.

  **The PONs are unaffected** — verified, not assumed. Their cohorts were built
  with `--skip-pon`, so no `_logR` column ever existed for the builder to
  mistake for a bin, and all four models carry exactly 41 arms × 67 bins with no
  duplicates.

- **Two more fabrications in the same function.** `log_ratio = 0.0` when the arm
  is absent from the baseline, and `pon_stability = 1.0` when no σ was measured.
  Zero says "exactly at the healthy baseline" and 1.0 is what an average
  variance of 0.99 would give — both specific readings asserted where nothing
  was measured. Now NaN, as in `zscore_or_nan`, the FSC log2, the FSR ratio and
  the entropy z.

- **The vector σ had the same residue defect, in a different function.** The
  previous fix went to `mean_and_sd`, which serves the *scalar* shape
  statistics. `wps_baseline`'s 200-element `wps_nuc_std` comes from
  `element_wise_std`, a separate implementation carrying its own `sd > 0.0`
  guard — so the fix left every vector σ untouched, and a re-aggregated cohort
  came back byte-identical.

  This is the one that matters. Measured in the shipped models: residue σ
  reaches 4.5% of usable positions in `xs1.all_unique`, 11.8% in
  `xs2.all_unique`, 47.2% in `xs1.duplex` and 55.4% in `xs2.duplex`, touching
  38–74% of anchors. A typical real σ there is 0.4–1.2, so a one-unit deviation
  scores z ≈ 10¹¹.

- **The Rust `mean_and_sd` reported floating-point residue as a spread.**
  `if sd > 0.0 { sd } else { NAN }` is the same hole `sample_std_or_nan` closed
  on the Python side, and it misses the case it exists for: summation error
  grows with n, so identical donors give residue rather than a clean zero.
  21 copies of `0.1` — the xs2 cohort size — yield `7.6e-9`; 8 copies of
  `-0.3333333` yield `3.2e-8`. Both pass `> 0.0` and become divisors. The
  existing test used `[3.0, 3.0, 3.0]`, which cancels exactly, so it passed
  throughout.

  Measured in the freshly rebuilt models, where these reach `compute_z_vector`
  through the same `std > 0` guard:

  | model | positions affected | anchors touched |
  |---|---:|---:|
  | `xs1.all_unique` | 4.5% | 38% |
  | `xs2.all_unique` | 12.2% | 57% |
  | `xs2.duplex` | 55.4% | 74% |

  A typical real σ there is 0.4–1.2, so a one-unit deviation against a residue
  σ scores **z ≈ 10¹¹**. Those are positions where every donor looked the same,
  which is exactly where a sample doing something unusual should have been
  reported as absent rather than astronomically significant.

  Identity is now tested on the values (`min == max`), which is exact and
  invents no tolerance — values one ULP apart keep their spread, pinned by a
  test.

  Found by inspecting `wps_baseline`, the block that is 99% of the file and the
  one `validate-pon` cannot check, because its σ lives inside 200-element
  vectors rather than a column.

- **`build-pon` ran region-MDS at a different fragment cap than `run-all`.**
  It hardcoded `max_len=400`; `run-all` uses 1000, and `run-all` is how every
  sample is measured at scoring time. MDS is normalised entropy, so the wider
  window sees more fragments and reads higher — a baseline fitted over
  65–400 bp and a sample measured over 65–1000 bp are not the same quantity.

  Measured on the xs2 duplex cohort, the same 21 donors aggregated both ways:
  all 146 genes shifted up by 0.0043 ± 0.0016. Against the baseline's own σ of
  0.0043 that is a **median z bias of +1.15** — 86 genes past one σ and 11 past
  two — in every healthy sample, from the cap alone.

  Both call sites now use `sample_params.maxlen`, which defaults to the same
  1000. Invariant #6, caught by diffing a rebuilt PON against its predecessor
  rather than by any test, because nothing crashed.

  The four models rebuilt via `--from-outputs` are unaffected: they aggregate
  `run-all` output, so they already carry the 1000 bp measurement.

- **A kept `build-pon` cache could not be re-aggregated.** `--keep-sample-outputs`
  exists so a rebuild does not re-read every BAM — 55–97 minutes per sample,
  ~15 hours for 47 — and every defect found in the 0.9.0 models has been in
  *aggregation*, so re-checking one should cost minutes.

  It could not, because `process_sample` never called `write_motif_outputs`:
  `run-all` did it at its own call site and `build-pon` did not. A kept cache
  therefore had no `EndMotif` or `MDS` table, even though the same k-mer counts
  reached the model through memory. That is the one input `--from-outputs`
  cannot reconstruct, and `krewlyzer motif` needs a BAM, so it could not be
  backfilled either. Measured against the real 0.9.0 cluster cache: all 47 xs1
  duplex directories were refused, each with `missing .EndMotif`.

  `process_sample` now takes `write_motif_files`, which `build-pon` passes.
  It defaults to `False` so `run-all` keeps writing them itself and nothing is
  written twice. The call has to live there because `write_motif_outputs` takes
  an `ExtractionResult`, which is live inside `process_sample` and never
  returned.

  The resume story is deliberately *"re-aggregate the cache"* rather than
  *"skip finished samples in the extraction loop"*. The in-process path
  collects from live `SampleOutputs` objects, not from the directory, so
  skipping extraction would skip collection too and drop the sample from every
  baseline without saying so. On success the build now prints the exact
  `--from-outputs` command instead of an unactionable "reuse them".

- **`--from-outputs` required extraction assets it never reads.** No BAM, no
  reference, no bin file — but `build-pon` refused when the bundled bin file
  was missing, before reaching the branch that would not have opened it. The
  bin file is LFS-backed, so every `--from-outputs` build failed in CI while
  passing locally.

- **Two Rust baseline readers parse plain TSV and were being handed parquet.**
  `compute_fsd_baseline` and `compute_region_mds_baseline` read with
  `BufReader::lines()` — no gzip, no parquet. `File::open` succeeds on both
  anyway, the header parse then finds no usable columns, and every sample is
  skipped: three samples in, zero out.

  Not an edge case. A `run-all` directory writes `.parquet` and `.tsv.gz` and
  *no* plain `.tsv`, so it is the normal case for the layout `--from-outputs`
  reads. Found against the real 0.8.3 corpus, where `region_mds` came back
  empty from the same files `region_mds_exon` read fine.

  The two failed differently and the quieter one was worse. FSD raised, but
  blamed the backend — "no data returned from Rust", the same misleading
  wording already fixed once in `_compute_wps_baseline`. `region_mds` logged
  one warning and returned `None`, so the block was simply *absent* from the
  model; `validate-pon` iterates the blocks present in the file and skips empty
  ones, so a block that vanished entirely cannot be checked. Build, gate and
  downstream all reported success.

  Inputs are now normalised to plain TSV where the constraint lives, in the two
  callers that know about it. Unreadable inputs are reported rather than
  dropped. Both empty-result paths raise and name an input file.

- **`region_mds` fell back to a standard normal.** `data.get("mds_mean", 0.0),
  data.get("mds_std", 1.0)` makes z equal the raw MDS — about 0.95, an
  ordinary-looking number for a statistic that was never computed. Third
  instance after `_compute_mds_baseline` and `get_periodicity_stats`.

- **A zero spread was reported as a measurement in three of the four PON
  builders.** `pandas.std(ddof=1)` returns `0.0` for identical values, not NaN
  — NaN only when `count < 2`. Only `_compute_fsc_gene_baseline` converted
  that; the other three did not, while carrying comments claiming they did
  ("NaN where pandas could not measure a spread", "an unmeasurable spread
  yields NaN"). A comment asserting behaviour the code lacks is worse than no
  comment: it stops the next reader checking. The Rust twin `mean_and_sd` was
  correct, which is why the WPS blocks passed `validate-pon` and the Python
  ones did not.

  One `sample_std_or_nan` helper now serves all six sites. Identity is tested
  on the values rather than the result: `std([0.95, 0.95, 0.95])` is
  `1.36e-16`, not `0.0` — cancellation in the variance sum leaves residue, so
  a `sd <= 0` guard misses the exact case it exists for, and a z divided by
  `1.36e-16` is `1e16`, worse than the zero it replaced.

- **The PON averaged the edge of the NRL search band in with real
  measurements.** `nrl_bp = 250` is the top of the FFT window, and
  `nrl_at_band_limit` says so on every affected row — the builder ignored the
  flag. Across the xs2 duplex cohort all 174 band-limited rows carry exactly
  250.0, one unique value with zero variance, against 194.9 ± 24.0 for the
  other 414. Twelve of 28 groups are partly limited, four entirely, so those
  groups reported the edge of the search as the healthy expectation.

  The NRL is now fitted from the rows that measured one, and the group keeps
  its row either way — absent, not dropped, because a group xs1 can measure
  and xs2 cannot is information when the two models are compared. New
  `n_at_band_limit` and `n_nrl_fitted` columns record why a baseline is
  missing. Below `MIN_SAMPLES_PER_KEY` measured rows the mean goes too, not
  just the spread: a cohort baseline averaged over one donor is the same
  fabrication as a hardcoded one.

  Periodicity is deliberately still fitted from *all* rows. Those same 174
  band-limited rows hold 174 distinct periodicity values spanning 0.37–0.86 —
  only the peak's position hit the edge, its strength was measured. `nrl_bp`,
  `periodicity_score` and now `nrl_at_band_limit` are all required.

- **`WpsBackgroundBaseline.get_periodicity_stats` invented a standard normal.**
  `row.get("periodicity_mean", 0), row.get("periodicity_std", 1)` made a
  missing column produce a z equal to the raw periodicity — about 0.47, an
  unremarkable number nobody would have questioned. The same fabricated
  baseline the block itself was rebuilt to remove, one level down in the
  reader. Now indexed, so a missing column raises.

  Verified against all three models already built on the cluster by rebuilding
  only the affected blocks from cached per-sample outputs — no BAM re-read.
  `xs1.duplex` went from three `PON.NONPOSITIVE_SIGMA` findings to none;
  `xs2.duplex` and `xs2.all_unique` from one each to none. All three now exit
  0 under `validate-pon`.

- **A `region-mds` failure during a PON build was swallowed by `except: pass`.**
  `region_mds` and `region_mds_exon` are both built by globbing the files that
  step writes, so a failure means two baselines silently missing from the
  model — and the build reported success. Found while inspecting a live
  rebuild: 4 of 4 completed samples had no `MDS.gene.tsv`, and `grep -i mds`
  over 1,803 log lines returned nothing, because there was nothing to find.

  Now warns with the exception type and message, and separately says when the
  gate did not open at all — no gene BED resolved, or the input is a fragment
  BED rather than a BAM. At aggregation, a BAM cohort that produced no MDS
  files at all is an error rather than an empty block.

- **`scripts/build_pon.sh` refuses to run on a pre-0.9.0 krewlyzer.** A 0.8.x
  `build-pon` rebuilds every defect this release fixes, exits 0, and produces
  logs that look clean — the failure is invisible until someone opens the
  model. Caught in practice: the first 0.9.0 rebuild attempt ran 18 minutes on
  `v0.8.3` before the banner was noticed.

  Also corrects the sample-list default from `{assay}_all_unique_pon.txt` to
  `{assay}_allUniq_pon.txt`, which is what the lists are actually called — the
  script had invented a convention rather than using the existing one.

- **`nrl_z` and `periodicity_z` never reached a single output file** — and the
  machinery was fully plumbed. `compute_nrl_zscore` defaulted to
  `group_id="all"`, and no PON has ever held a group by that name: they are
  `Global_All`, `Chr1_All` … `Family_AluY`. Every lookup missed, the function
  returned `None`, and `None` was indistinguishable from "this PON has no
  baseline". The default is now a named `GENOME_WIDE_GROUP` constant, defined
  once, and a group that matches nothing logs what it did find.

  Two defects had to stack for the column to be absent: the fabricated
  `wps_background` baseline (`4cd634b`) *and* this. With the old PON the value
  is `−3.4` for every sample, which is what the fabricated `167.0 / 5.0` gives
  for an `nrl_bp` pinned at 150.

- **27 of 28 NRL baselines were built and never used.** Only `Global_All` was
  scored, though the output carries an `nrl_bp` per chromosome and per Alu
  family and the baseline has one for each. All 28 groups are now scored;
  per-chromosome NRL drift is the reason those baselines exist.

- **An unmeasurable spread produced `z = 0.0`, in twelve places.** Ten baseline
  classes in `pon/model.py` independently ended with
  `if std > 0: return (x - mean) / std` followed by `return 0.0` — nine copies
  of the same three lines, and the same defect nine times. Zero is not a
  cautious answer: it is the most confident claim the column can make ("sits
  exactly at the healthy baseline"), asserted precisely when the baseline
  measured no spread, and indistinguishable from a genuine zero.

  Replaced by one `zscore_or_nan` helper, so it can only be wrong once. Also:
  `FsdBaseline.get_expected`/`get_std` returned `0.0` for an unknown arm (a
  caller dividing by that gets infinity, not an absence), and
  `compute_shape_score` returned `0.0` for an undefined correlation, which is
  a real claim about shape agreement.

- **`WpsBaseline.compute_z_vector` substituted `1.0` for an unusable σ.** Since
  `4cd634b` the builder writes NaN where it could not measure spread;
  `np.where(std > 0, std, 1.0)` is False for NaN, so every one of those became
  a *finite* z — undoing the builder's honesty at the read side, which is the
  harder place to notice. This is the function WPS z-scoring will use.

- **Reading back a file just written could return an older run's data.**
  `read_table` is parquet-first: given `x.tsv` it prefers `x.parquet`. That is
  right when asking "what was produced for this output" and wrong immediately
  after writing a specific file. The Rust backends write plain TSV and Python
  reads it back to honour `--output-format`, so a `.parquet` sibling left by an
  earlier run into the same directory won — which is what a Nextflow retry or
  `-resume` produces.

  Reproduced end to end for region entropy: with a stale parquet present, the
  output carried the **previous run's** `z_score` and `entropy` rather than
  this run's. Fixed at every read-after-write site — region entropy, mFSD, UXM,
  region-MDS gene and exon, on-target MDS — via a new
  `read_exact_table` helper.

  Also reverts an over-correction from earlier in this branch: the post-Rust
  existence check in `process_region_entropy` must be exact, not resolving, or
  a stale sibling satisfies it and the degradation branch never fires.

- **A whole-genome run lost an entire verdict axis.** A panel run splits OCF
  into on- and off-target; a WGS run writes neither and emits a plain
  `.OCF.parquet`. `verdict.py`, `plots.py` and the report's run-facts strip all
  read only the panel pair, so every WGS sample reported OCF absent with the
  file sitting in the directory. Losing an axis makes the verdict read
  *stronger*, not weaker — it shrinks the denominator — which is the dangerous
  direction. The preference order now lives once, in `contract.py`, and all
  three read it.

- **Six optional outputs rendered with no interpretation.** `mFSD`, `UXM`,
  `OCF`, `FSC.regions.e1only`, `fsc_counts` and `correction_factors` got a
  section with a measured shape and nothing else — no lede, no direction, no
  why/how/what — because `test_meaning_registry.py` keyed on `CONTRACT` alone
  and those six are `NOT_CONSUMED`. mFSD in particular had a chart with no
  explanation beside it. The registry now covers everything either command can
  render.

- **Numeric identifier columns were not redacted.** `describe-output` prefers a
  min…max range over the example value whenever one exists, so an integer
  patient key or accession was printed in full by the branch meant to hide it —
  a range over one distinct value *is* the value.

- **`describe-output` and `report` disagreed about the same fact.** The report
  says *gated* / *inventory*; `describe-output` still said *read downstream*,
  a claim about another repository that nothing here keeps in sync. Aligned on
  the report's wording, with the distinction spelled out in the output.

- **README linked to six pages that do not exist.** `/introduction/`,
  `/installation/`, `/usage/`, `/features/extract/`, `/pipeline/` and
  `/citation/` all 404 — they point at a docs layout reorganised out of
  existence. Absolute URLs, so neither mkdocs nor the internal link checker
  ever resolved them.

- **`tests/unit/test_cli_documented.py` checked 8 of 19 commands.** Its regex
  matched only `app.command(name="x")` and silently skipped the eleven
  registered as bare `app.command()(fn)`. It is now cross-checked against typer
  itself, which immediately surfaced a fifth undocumented command, `validate`.

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

  **`output_format` now defaults to `parquet` (was `tsv`).** kreview reads
  Parquet only, and a tsv-only run additionally skips per-sample validation
  because the contract describes the Parquet surface — so the old default
  produced a cohort that was both invisible downstream and unchecked, silently,
  since every reader there swallows exceptions and yields an empty feature
  dict. The pipeline warns at start whenever the selection still excludes
  Parquet.

  **Breaking for anyone relying on the pipeline emitting TSV by default**; pass
  `--output_format tsv` or `both` to restore it.

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
