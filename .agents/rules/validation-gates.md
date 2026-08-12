---
description: What each gate catches, what it cannot, and the order to run them in
alwaysApply: false
---

# Validation gates

`AGENTS.md` lists which gate applies to which change. This file is about the
part that keeps going wrong: the **order** they must run in, and what each one
is blind to.

## Order matters more than the list

```bash
maturin develop --release --manifest-path rust/Cargo.toml   # 1. rebuild
pytest tests/ -q                                            # 2. then test
cargo test --manifest-path rust/Cargo.toml                  # 3. reads source
```

**Step 1 is not optional after any `rust/src/` change.** `pytest` imports
whatever `.so` was built last, so skipping the rebuild means the Python suite is
testing the *previous* binary and passing for the wrong reason. `cargo test`
reads the source directly and is unaffected — which is exactly why the two can
disagree, and why **Python passing is the weaker signal**.

This bites hardest when switching branches. A suite that passed on branch A and
fails on branch B — or vice versa — is usually a stale extension, not a real
difference. Rebuild before believing either result.

## What each gate cannot catch

| Gate | Blind to |
|---|---|
| `pytest tests/ -q` | Anything in Rust you did not rebuild. Also every test that passes vacuously — see below |
| `cargo test` | Any path only reachable through the Python API, and anything about the *bundled assets*, which it never loads |
| `black` / `ruff` / `mypy` | Only `src/krewlyzer/` and `tests/` are checked in CI. `scripts/` is **not** linted |
| `check_output_format.py` | Whether a CLI *exposes* an option — it only verifies internal call sites forward it |
| `test_claims_registry.py` | A number that appears in a doc but has no registry entry. The registry only checks what it knows about |
| `validate-output` | A single sample. Degeneracy is cross-sample by nature; use `validate-cohort` |
| `check_phi_guard.sh` | Identifiers inside `.gitleaks.toml` itself — gitleaks skips its own config |

## A test that passes proves nothing until it has failed

Every significant defect in this repo was *present, plausible, and passing*.
Adding a test that passes on fixed code tells you nothing about whether it would
have caught the bug.

**Revert the fix, watch the test fail, restore the fix.** It takes one rebuild
and it is the only reliable check. Real examples from this repo:

- A conservation suite reported *"exact match: True"* on **zero fragments**,
  because the BAM's proper-pair flag filtered everything out. Caught only by
  printing counts alongside the verdict.
- A first pass at the FSC band tests had **10 of 12 passing against the broken
  code** — a fragment in the wrong channel still sums correctly, so partition
  and ratio checks cannot see a routing bug.
- Three E1 tests built `RegionInfo` by hand and never called `parse_gene_bed`,
  so the entire format-detection path was untested and a column-count change
  silently disabled strand handling.

If you cannot make a new test fail, you have not established what it tests.

## Strict xfail as a forcing function

Pin a known defect with `@pytest.mark.xfail(strict=True)`. When the fix lands
the test **XPASSes and fails the build**, forcing the marker out. A non-strict
xfail rots silently instead.

## Assets are not covered by the Python gates

`pyproject.toml` excludes `data/` from the wheel, so the installed package has
no assets and anything resolving paths through `krewlyzer.__file__` finds
nothing in CI. Tests that validate the *committed* asset files must read from
the checkout (`Path(__file__).parents[2] / "src" / ...`); `conftest.requires_data`
is for tests that need the *runtime* to load them, and skipping is wrong for
the former.

After regenerating anything under `src/krewlyzer/data/`, `git lfs push --all`
before pushing the branch, or the push is rejected for unknown LFS objects.

## Before blessing a cohort

```bash
krewlyzer validate-output RESULTS_DIR     # one sample, or a directory of them
krewlyzer validate-cohort FINGERPRINTS    # cross-sample degeneracy
```

Exit `0` satisfied, `1` contract violation, `2` structural. A workflow should
retry on 2 and escalate on 1.

**Do not add a column to an exception list to make the gate green.** Declaring a
column legitimately constant requires a written `constant_reason`, and that cost
is the point.

## The local cohort gate

```bash
KREWLYZER_TEST_CORPUS=/path/to/scored/cohort python -m pytest tests/real_data -q
```

**Local only.** Unset the variable and it skips, so `pytest tests/` stays green
and CI never sees it. No cohort data is committed, and nothing here may put a
patient identifier in a test name, a parameter id, an assertion message or a
written artifact — sample directories are named for the patient (invariant #4).
`sample_label()` gives a non-reversible handle when output needs one.

**What it catches that nothing else does.** Every defect the PON audit found
needed a cohort to see:

| defect | what a single fixture showed |
|---|---|
| `wps_background` fabricated | nothing — one sample cannot reveal a constant |
| FSD unnormalised under parquet | nothing — the fixture PON matches no arm |
| single-sample WPS anchors | nothing — needs per-anchor `n` across samples |
| `wps_baseline` consumed by nothing | nothing — the code reads fine |

The unit suite proves a scorer *can* run; only a cohort proves it *did*, and
only a cohort can tell a varying metric from a constant one.

**Blind spot.** The corpus is tumour samples, not the healthy cohort the PON
was fitted on, so a shifted z-score mean is expected biology. The calibration
checks here are deliberately loose — they reject a scale error (a floored
sigma, a mismatched baseline), not a shifted one.
