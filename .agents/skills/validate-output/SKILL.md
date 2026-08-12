---
name: validate-output
description: Check a finished krewlyzer output directory against the downstream contract
---

# Validating output

```bash
krewlyzer validate-output RESULTS_DIR --json-report report.json
```

`RESULTS_DIR` is a cohort laid out as `{dir}/{sample_id}/`, or a single sample
directory. Exit codes: **0** satisfied, **1** contract violation, **2**
structural — the inputs are not comparable (missing directory, unreadable
Parquet). Retry on 2; escalate on 1. A partial network read produces 2, which
is why the two are distinct.

At cohort scale, scatter and gather instead:

```bash
krewlyzer validate-output SAMPLE_DIR --fingerprint-out fp/SAMPLE.json   # per sample
krewlyzer validate-cohort fp/ --json-report cohort.json                 # once
```

The fingerprint is ~20 KB — a hash and two counts per column. A 14k-sample
cohort reduces from ~270 MB instead of re-reading ~21 TB. `run-all` writes one
automatically on any Parquet run, which is what keeps the cohort step
affordable.

## Reading the result

Findings come in four categories. Take them differently:

- **completion** — `.metadata.parquet` absent. Highest severity despite looking
  trivial: downstream drops the sample from the cohort with no warning.
- **missing / schema** — a table or column is absent or the wrong type. Usually
  an interrupted run; check the sample completed before assuming a code defect.
- **domain** — an invariant broke (frequencies not summing to 1, bare
  chromosome names, channels not partitioning the total). Usually a real
  defect.
- **degeneracy** — a metric is identical across every sample. It carries no
  information and cannot be a function of the input.

Degeneracy is reported **SKIP, never PASS** below two samples. One sample
cannot demonstrate variation, and reporting "pass" there is how a gate ends up
certifying a constant.

## The anti-pattern

**Do not add a column to an exception list to make the gate green.** If a
column is legitimately constant, declare it `Vary.NEVER` in `contract.py` with
a written `constant_reason` — the dataclass refuses the declaration without
one. That reason is what lets the next person judge whether it is still true.

Four `WPS_background` columns were constant across an entire production cohort
while passing every schema check that existed. Silencing the finding would have
preserved exactly that state.
