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

These read a finished output directory rather than producing features.

| Command | Purpose |
|---------|---------|
| [`describe-output`](#describe-output) | What is in each output file — shape, columns, ranges |
| [`report`](#report) | Single-sample HTML report: verdict, charts, interpretation |
| [`validate-output`](#validate-output) | Check one sample against the output contract |
| [`validate-cohort`](#validate-cohort) | Cross-sample degeneracy checks over fingerprints |

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
    Re-run with `--pon-model` for anything comparative.

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
