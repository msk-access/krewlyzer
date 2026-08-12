# krewlyzer — standing context (canonical, cross-tool)

> Canonical agent rules for this repo. Read by Antigravity, Cursor, Codex, Copilot,
> and (via `CLAUDE.md` import) Claude Code. Every line is read on every turn — keep it
> tight. Depth lives in `.agents/rules/*.md` and is read on demand, not here.

**What this is:** `krewlyzer` — a cfDNA fragmentomics toolkit (MSKCC / MSK-ACCESS).
Extracts fragment-level features from BAMs: size distributions (FSD/FSR/FSC), end
motifs and motif diversity (MDS), nucleosome protection (WPS), orientation-aware
fragmentation (OCF), region entropy (TFBS/ATAC), fragment methylation (UXM), and
mutant fragment size distributions (mFSD). Ships as a pip package, a Docker image,
and a Nextflow pipeline.

**Architecture:** Rust core (`rust/src/`, pyo3 `extension-module`) behind a Python
CLI (`src/krewlyzer/`). Both the individual tools and `run-all` funnel through
`core/unified_processor.run_features()`, which drives a **single pass** over the
fragment BED for all enabled features — that single pass is why `run-all` does not
simply call the tools in turn. Full module map: `.agents/rules/architecture.md`.

## Load-bearing invariants (don't break these)

<!-- protected:start reason="each of these encodes a defect that shipped; do not weaken without review" -->

1. **A metric that cannot vary with the data is worse than a missing one.**
   It is present, plausible, and survives every schema check, so it gets
   modelled. `nrl_bp` was 150.0 for every sample ever produced;
   `periodicity_score` 0.3333; `adjusted_score` 0.0. For every derived metric,
   assert that **two different inputs produce two different outputs**
   (`tests/invariants/test_antidegeneracy.py`). That single assertion would
   have caught the largest defect in this codebase.

2. **Parquet is the product.** The downstream consumer reads Parquet only and
   never `features.json`. `{sample}.metadata.parquet` is its completion marker
   — a sample without it is dropped from the cohort **silently**. Every reader
   there swallows exceptions and yields an empty feature dict, so a malformed
   table becomes missing columns, not an error. Nothing downstream can tell us
   we shipped something wrong: assert it here or not at all
   (`krewlyzer validate-output`).

3. **A boundary value is not a measurement.** When a search hits its own
   window edge, say so in a separate column rather than reporting the edge as
   a result (`nrl_at_band_limit`). Same failure as (1), one level down.

4. **Patient identifiers must never enter this repository.** Sample directories
   are named for the patient. Use the `P-0000000` and `C-000000-L000`
   placeholders in docs and tests, and never a real identifier inside
   `.gitleaks.toml` itself — gitleaks skips its own config, so only
   `scripts/check_phi_guard.sh` would catch it. Enable the hooks once per
   clone: `git config core.hooksPath .githooks`.

5. **Constants live in one place and are asserted.** Documentation drifted from
   code for a year — a tolerance of 20 documented as 50, a search band of
   140-200 documented as 150-250. Any numeric constant a doc quotes belongs in
   `src/krewlyzer/validate/claims.py`, which fails when the two disagree.

6. **Single-tool and `run-all` output must match.** Every feature CLI accepts
   `--output-format` and `--compress`; eight of eleven once did not, so they
   could only write TSV while `run-all` wrote Parquet.

<!-- protected:end -->

## Gates

The **Read** column is not optional reading — it is the trigger. `.agents/rules/`
is described as read "on demand", and nothing else in this repository creates
demand, so a rule nobody opens is a rule that does not exist.

| When | Run | Read first |
|---|---|---|
| Any change | `pytest tests/ -q` · `black` · `ruff` · `mypy` | — |
| Rust change | `cargo test --manifest-path rust/Cargo.toml --lib` · `cargo clippy` | `rules/architecture.md` — the boundary table, and the three traps in reading a PON parquet |
| Output-writing change | `python scripts/check_output_format.py` | `rules/output-contract.md` — what downstream relies on |
| Docs or constants | `pytest tests/unit/test_claims_registry.py` | — |
| Before blessing a cohort | `krewlyzer validate-output RESULTS_DIR` | `rules/validation-gates.md` — what each gate is blind to |
| Before blessing a cohort or a release | `KREWLYZER_TEST_CORPUS=<dir> pytest tests/real_data` (local only) | `docs/development/release-guide.md` |
| PHI rules touched | `bash scripts/check_phi_guard.sh` | invariant #4 below |
| Build or QC tooling | — | `rules/development.md` |

Order matters: rebuild the extension before `pytest`, or you are testing the
previous binary. What each gate is blind to, and why a passing test proves
nothing until it has failed, is in `.agents/rules/validation-gates.md`.

`cargo test` runs on macOS in CI: `krewlyzer_core` is a pyo3 `extension-module`
cdylib and does not link libpython, so the test binary cannot resolve Python
symbols at link time on Linux.

## Conventions

- **Shared context** beyond these invariants lives in `.agents/rules/`
  (normative, read on demand) and `.agents/memory/` (conventions, modelling
  philosophy, and decision records). Anything machine- or operator-specific
  belongs in your own memory directory, not the repo.

- Git Flow (`docs/development/release-guide.md`): branch off `develop`, one
  `fix/*` or `feature/*` per concern, conventional commits, docs and CHANGELOG
  in the **same** commit as the change they describe.
- Never commit data artifacts. Reference data belongs in `src/krewlyzer/data/`
  (excluded from the wheel to stay under PyPI's 100 MB limit; shipped via
  Docker or an LFS clone), fixtures in `tests/data/fixtures/`.
- `git lfs pull` before running the full suite, or LFS-dependent tests fail
  with "Corrupt footer" or "Not a gzipped file".
- A cold `cargo build` takes ~10 minutes. It is not hung.
