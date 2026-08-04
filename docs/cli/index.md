# CLI Reference

Krewlyzer provides a command-line interface for all feature extraction tools.

## Main Command

```bash
krewlyzer run-all -i sample.bam --reference hg19.fa --output results/
```

See [run-all](run-all.md) for the unified command that runs all features.

## Individual Commands

| Command | Description |
|---------|-------------|
| `extract` | Extract fragments from BAM |
| `fsc` | Fragment Size Coverage |
| `fsd` | Fragment Size Distribution |
| `fsr` | Fragment Size Ratio |
| `wps` | Windowed Protection Score |
| `motif` | End Motif extraction |
| `ocf` | Orientation-aware Fragmentation |
| `region-entropy` | TFBS/ATAC entropy |
| `region-mds` | Per-gene MDS |
| `mfsd` | Mutant Fragment Size Distribution |
| `uxm` | Fragment-level Methylation |
| `build-pon` | Build Panel of Normals |
| `build-gc-reference` | Build GC reference assets |

## Inspection and Validation

These inspect inputs or a finished output directory rather than producing
features.

!!! warning "`validate` and `validate-output` are different commands"
    `validate` checks the **assets going in** — is this BED well-formed, does
    it have the columns the Rust core expects. `validate-output` checks the
    **results coming out** against the downstream contract. One hyphen apart
    and neither substitutes for the other.

| Command | Purpose |
|---------|---------|
| [`validate`](#validate) | Check input **assets** before a run — BEDs, anchors, GC factors |
| [`describe-output`](#describe-output) | What is in each output file — shape, columns, ranges |
| [`report`](#report) | Single-sample HTML report: verdict, charts, interpretation |
| [`validate-output`](#validate-output) | Check one sample against the output contract |
| [`validate-cohort`](#validate-cohort) | Cross-sample degeneracy checks over fingerprints |
| [`validate-pon`](#validate-pon) | Check a PON **before** anything is scored against it |
| [`stamp-pon`](#stamp-pon) | Record the release a built PON ships with |

---

### `validate` {#validate}

Checks the reference assets a run depends on, before the run rather than after
it. A malformed BED does not usually crash — it produces a table with fewer
rows than it should, which looks exactly like a quiet sample.

```bash
krewlyzer validate --genome hg19                    # every bundled asset
krewlyzer validate --genome hg19 --assay xs2        # ...plus assay-specific
krewlyzer validate --gene-bed my_genes.bed          # one custom file
```

| Option | Checks |
|--------|--------|
| `--genome` / `-G` | All bundled assets for `hg19` or `hg38` |
| `--assay` / `-A` | Assay-specific assets too (requires `--genome`) |
| `--gene-bed` | Gene BED: chrom, start, end, gene, \[name\] |
| `--targets-bed` | Panel target regions |
| `--arms-bed` | Chromosome arms: chrom, start, end, arm |
| `--wps-anchors` | WPS anchors, BED6 |
| `--ocr-file` | OCF open-chromatin regions |
| `--gc-factors` | GC correction factors TSV |
| `--bin-file` | Genomic bins, BED3 |

Use `-G` before reporting a bundled-asset problem: a partial `git lfs pull`
leaves pointer files in place of the real assets, and this is what says so.

---

### `describe-output` {#describe-output}

Answers the question people ask before "is this correct?" — what are these
files and what is in them.

```bash
krewlyzer describe-output RESULTS/{sample_id}/            # Markdown to stdout
krewlyzer describe-output RESULTS/{sample_id}/ -o out.md  # or to a file
krewlyzer describe-output RESULTS/{sample_id}/ -o out.html
```

Per table: rows, columns, size, and whether it is gated by the contract. Per
column: dtype, numeric range or example value, distinct count, null count.

Everything is measured from the file or read from the output contract, so a
table that gains a column changes the description automatically.

!!! note "Identifier columns are redacted"
    Sample directories are named for the patient and several tables carry the
    sample id as a *column value*. Knowing a column holds an identifier is the
    useful fact; which identifier is not. Everything else is the sample's real
    data — treat the output accordingly.

---

### `report` {#report}

A single-sample HTML report: cross-axis verdict, 16 charts, and each table's
data with its interpretation beside it.

```bash
pip install 'krewlyzer[report]'          # plotly, needed for charts
krewlyzer report RESULTS/{sample_id}/ -o report.html
krewlyzer report RESULTS/{sample_id}/ -o report.html --z-threshold 2.5
```

| Option | Default | Description |
|--------|---------|-------------|
| `--output` / `-o` | `report.html` | Where to write |
| `--sample-id` | directory name | Override the sample id |
| `--z-threshold` | `2.0` | \|z\| at which a verdict axis is flagged |

**Verdict.** Four independent axes — fragment size, nuclease cutting, tissue
shedding, chromatin accessibility — reported as *"N of M assessable axes
agree"*, never a single composite score. Direction differs per axis, and an
axis with no PON z-score is *not assessable* rather than counted as
disagreement.

**Organised by table.** Each section carries one output's chart, its why / how
/ what, and its columns together.

!!! warning "Internal use"
    The report contains one patient's measurements. It is generated on demand
    and is not published, committed, or shipped with the documentation. Use
    `describe-output` for anything structural that needs to leave the machine.

!!! info "Without a panel of normals"
    Every z-score and PON log-ratio is absent, which is most of the
    interpretable surface — three of the four verdict axes cannot be assessed.
    The report leads with a banner saying so and marks each affected column.
    Re-run the *feature extraction* with `run-all --pon-model ...` and
    report the new output; `report` has no PON option of its own, because
    it reads what is already on disk and never recomputes a feature.

`--z-threshold` is a **convention, not a clinical cut-off**, and the report
says so. No chart draws a threshold as a cut-off; reference lines appear only
where the value is a labelled literature anchor.

---

### `validate-output` {#validate-output}

Checks a finished directory against the contract its consumers rely on.

```bash
krewlyzer validate-output RESULTS/                       # a cohort directory
krewlyzer validate-output RESULTS/ --json-report out.json
krewlyzer validate-output RESULTS/{sample}/ --fingerprint-out fp.json
```

Three layers: every consumed table present and shaped correctly; domain
invariants; and **anti-degeneracy** — a metric identical across every sample is
an error, not a pass.

| Exit | Meaning |
|------|---------|
| `0` | Contract satisfied |
| `1` | Contract violation |
| `2` | Structural failure — missing directory, unreadable Parquet |

A workflow should retry on `2` and escalate on `1`.

---

### `validate-cohort` {#validate-cohort}

The gather half. Degeneracy is inherently cross-sample: every sample can pass
alone while a metric is constant across all of them.

```bash
krewlyzer validate-cohort FINGERPRINTS/ --min-samples 3
```

Reads the ~20 KB fingerprints written by `validate-output --fingerprint-out`,
never the outputs themselves, so a large cohort reduces from megabytes rather
than terabytes. Below `--min-samples` the check reports SKIP, never PASS.

## See Also

- [Nextflow Pipeline](../nextflow/index.md) - Batch processing
- [Features](../features/index.md) - Feature documentation

---

### `validate-pon` {#validate-pon}

`validate-output` gates the results of a run. This gates the **reference those
results are measured against**.

```bash
krewlyzer validate-pon model.pon.parquet
krewlyzer validate-pon src/krewlyzer/data/pon/GRCh37/*/*.parquet --json-report pon.json
```

Same exit codes as `validate-output`: `0` satisfied, `1` violation, `2`
structural.

| Check | Why |
|---|---|
| **σ is not one value repeated** | A baseline that cannot vary with its cohort was not fitted to one |
| σ is positive and finite | A z divided by zero is infinite, not conservative |
| `krewlyzer_version` recorded | 0.9.0 changes what every feature *means* |
| `cohort_digest` recorded | Otherwise the model cannot be reproduced or compared |
| entries backed by ≥ 3 samples | Below that a baseline entry is an anecdote |

!!! warning "Every PON shipped before 0.9.0 fails this"
    They carry a `wps_background` block hardcoded to `167.0 / 5.0 / 0.0 / 1.0`
    — identical across all 28 groups and all four models, from cohorts of 21
    and 47 samples — and no provenance at all. That is not a false positive;
    it is the reason the command exists. Rebuild with `build-pon`.

**Provenance contains no identifiers.** The cohort is recorded as a salted,
non-reversible digest of the sample IDs, plus an optional free-text
`--cohort-label`. Two builds from the same cohort produce the same digest; the
digest reveals nothing about who is in it. A PON ships in this repository, in
the Docker image and on PyPI, so it is the last place a patient identifier may
appear.

---

### `stamp-pon` {#stamp-pon}

```bash
krewlyzer stamp-pon model.pon.parquet --version 0.9.0
krewlyzer stamp-pon src/krewlyzer/data/pon/GRCh37/*/*.parquet --version 0.9.0 --dry-run
```

A PON is built from `develop`, where `pyproject.toml` still reads the previous
release — so the model records that version however new the code is. Rather
than bump before a four-hour build, set it here when cutting the release.

**This changes what the field means.** Afterwards `krewlyzer_version` is *the
release this model is published with*, not the code that produced it. That is
the definition a compatibility guard needs. `build_date` is untouched, so the
two together still say when it was actually built.

!!! warning "A stamp is an assertion, so it has to be earned"
    `stamp-pon` runs `validate-pon` first and **refuses on failure**. Without
    that it would be the shortest path to laundering: run it on one of the
    models carrying the fabricated `167.0 / 5.0` baseline and it would claim
    exactly the compatibility the guard exists to deny.

    `PON.NO_VERSION` alone does not block — that is the condition being fixed.
    `--force` exists for re-stamping a model that already passed, and says so.

Only the metadata row changes; every baseline is copied through unchanged and
the cohort digest is unaffected. **Re-run `validate-pon` afterwards** — the
file that ships should be the file that was checked.
