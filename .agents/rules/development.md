---
description: Krewlyzer development guide — project structure, build commands, QC checklist, code quality standards
alwaysApply: true
---

# Krewlyzer Development Guide

## Quick Start

```bash
# Install in editable mode with Rust build
maturin develop --release

# Run tests
pytest tests/

# Run CLI
krewlyzer run-all -i sample.bam -r hg19.fa -o output/
```

## Project Structure

```
krewlyzer/
├── src/krewlyzer/       # Python package
│   ├── wrapper.py       # run-all CLI
│   ├── data/            # Bundled assets (PON, targets, bins)
│   ├── core/            # Processors (sample, motif, fsc, fsd, wps, ocf)
│   ├── pon/             # PON model and build logic
│   └── utils/           # Shared utilities
├── rust/                # Rust extension (top-level)
│   ├── src/
│   │   ├── pipeline.rs  # Single-pass FSC/FSD/WPS/OCF
│   │   ├── mfsd.rs      # Mutant fragment size distribution
│   │   ├── wps.rs       # Windowed protection score
│   │   ├── ocf.rs       # Orientation-aware fragmentation
│   │   └── ...          # Other modules
│   └── Cargo.toml
├── nextflow/            # Nextflow pipeline
│   ├── main.nf
│   ├── nextflow.config
│   ├── modules/local/   # Per-tool NF modules
│   └── subworkflows/    # Subworkflows
├── tests/               # Test suite
└── docs/                # MkDocs documentation
```

## Building (Unified Python + Rust)

### `cargo build` is not a gate

A bare `cargo build --manifest-path rust/Cargo.toml` **fails on macOS**, and a
failure there says nothing about your change:

```
error: linking with `cc` failed: exit status: 1
... -o .../libkrewlyzer_core.dylib ... -dynamiclib -nodefaultlibs
```

`krewlyzer_core` is a pyo3 `extension-module` cdylib, so it deliberately does
not link libpython. Linking the dylib on its own therefore leaves the `Py*`
symbols unresolved, and the macOS linker refuses that; maturin supplies the
extra link arguments that make it work. Verified on macOS on an untouched
checkout, so it is pre-existing rather than anything you did. Not checked on
Linux, where the linker is laxer about undefined symbols in a shared library
and it may well succeed -- which is its own reason not to treat it as a gate:
its result depends on the host, so it is not a signal either way.

The gates that do mean something are
`cargo test --manifest-path rust/Cargo.toml --lib`, `cargo clippy`, and
`maturin develop --release` from the project root.

Note that `cargo clippy` must be run with `--tests`. Without it only the lib
target is linted and test code is skipped entirely -- which in this repo means
the invariant assertions are the one place the linter never looks.

### Toolchain prerequisite: libclang

The Rust build needs **clang and libclang headers** present, not just a Rust
toolchain. `rust-htslib` enables `hts-sys/bindgen` unconditionally, so every
build regenerates htslib's FFI bindings from the C headers using the host
libclang (`apt install clang libclang-dev`, or Xcode's on macOS).

Two consequences worth knowing before you lose an hour to either:

- **A libclang newer than the bindgen in `rust/Cargo.lock` does not error --
  it silently produces opaque structs**, so the build fails inside
  *rust-htslib's own source*, not in ours:

  ```
  error[E0609]: no field `n_targets` on type `sam_hdr_t`
    = note: available field is: `_address`
  ```

  The fix is `cargo update -p hts-sys`, never a source change. See the comment
  at `rust-htslib` in `rust/Cargo.toml`.
- **CI only exercises this on Linux**, and `test.yml` skips documentation-only
  changes (`docs/**`, `**.md`, `mkdocs.yml`). A docs-only stretch of merges can
  leave the Linux build untested for months while a runner-image bump breaks
  it; that is exactly how the clang 21 break reached a release workflow
  unnoticed. Use `workflow_dispatch` on **Run Unit Tests** after a long
  docs-only run.

### Which CI runs on which change

Two workflows, split by cost rather than by topic:

| workflow | path filter | runs |
|---|---|---|
| `quick-checks.yml` | **none** | PHI / secret scan; the tests that read markdown |
| `test.yml` | skips `docs/**`, `**.md`, `mkdocs.yml` | maturin build, full pytest, wheel, Docker, `cargo test`, lint |
| `docs.yml` | docs paths on PRs | aggregate check, `mkdocs build --strict` |

The filter on `test.yml` is only safe because `quick-checks.yml` runs everything
a markdown-only change can break. Two consequences:

- **Never path-filter the PHI scan.** It lived in `test.yml` until it was
  noticed that #78, which changed only `CLAUDE.md`, merged with no scan at all.
  Prose is where an identifier is most likely to be pasted by hand.
- **A test that reads markdown must run in `quick-checks.yml`.** That is
  asserted, not remembered: `tests/unit/test_quick_checks_cover_markdown_readers.py`
  walks the AST of every test for paths built to markdown and fails if one is
  not listed in the workflow. Mentioning a `.md` file in a docstring or a
  message does not count as reading it.

`*` in a GitHub path filter does not cross `/`, so `'*.md'` matches root-level
markdown only. Use `'**.md'` for markdown anywhere.

A PR's path filter is evaluated against the **whole PR diff** (a three-dot diff
against the merge base), not the latest push. Resolving a conflict in a PR that
changes code reruns the full suite, and should: the merge brought in whatever
landed on the base since, so the tree is one nobody has tested.

```bash
# Development build (editable install with Rust)
maturin develop --release

# Build wheel
maturin build --release --out dist

# Build Docker image
maturin build --release --out dist
docker build -t msk-access/krewlyzer:dev .
```

> [!CAUTION]
> **`maturin develop --release` MUST be run from the project root**, not from the `rust/` subdirectory.
> Running from `rust/` builds a `krewlyzer_core` wheel that does NOT update the `.so` at `src/krewlyzer/_core.cpython-*.so`.
> Running from the project root builds the `krewlyzer` wheel which correctly installs the `.so`.
>
> Verify the installed build timestamp:
> ```python
> python -c "import krewlyzer._core as c; import os, datetime; print(datetime.datetime.fromtimestamp(os.path.getmtime(c.__file__)))"
> ```

> [!WARNING]
> **Rust `Path::with_extension()` MUST NOT be used on compound extensions** like `.correction_factors.tsv` or `.ontarget.tsv`.
> Rust treats only the text after the **last** dot as the extension, so `with_extension("")` on `sample.correction_factors.tsv` strips `.tsv` to get `sample.correction_factors`, then `with_extension("tsv")` replaces `.correction_factors` → `sample.tsv` (wrong!).
> Use **string-based suffix stripping** instead: `output_path.strip_suffix(".tsv")`.

## Lint & Type Checking — Run Before Every Commit

```bash
bash scripts/lint.sh          # what CI runs, at the versions CI runs it
bash scripts/lint.sh --fix    # same, but let black and ruff rewrite
```

> **Agent rule:** run `scripts/lint.sh`, not the bare tools. Fix everything
> before pushing.

**Why the script and not `black`/`mypy` from your shell.** CI installs
`.[dev]`, so it runs the `==` pins in `pyproject.toml`. Your shell runs
whatever your environment happens to have, and the two drift without saying so.
Measured on one machine in August 2026:

| tool | pinned | local |
|---|---|---|
| `black` | 26.1.0 | 26.5.1 |
| `ruff` | 0.15.4 | 0.15.4 ✓ |
| `mypy` | 1.19.1 | **2.3.1** |

That mypy gap is not theoretical: 1.19.1 reports
`error: Library stubs not installed for "markdown" [import-untyped]` and 2.3.1
does not, so #71 passed locally and failed CI on it. The script reads the pins
from `pyproject.toml` and runs each tool through `uvx` at that exact version,
so there is one source of truth. Without `uvx` it still runs, but prints the
version it used, because a silently different version is the whole problem.

`tests/unit/test_lint_script_matches_ci.py` asserts the script keeps matching
the workflow — same tools, same flags — so it cannot quietly become a laxer
check that people trust.

### The individual commands

Useful for running one thing in isolation. **These are the commands, not the
versions** — prefer the script above when the answer needs to match CI.

```bash
# 1. Auto-fix import order and style issues
ruff check --fix src/krewlyzer/ tests/

# 2. Format code (auto-reformat)
black src/krewlyzer/ tests/

# 3. Check formatting is clean (what CI runs — no auto-fix)
black --check src/krewlyzer/ tests/

# 4. Type checking — full 0-error baseline (stub covers _core)
mypy src/krewlyzer/ --ignore-missing-imports --no-error-summary

# 5. Rust lints — auto-fix all fixable warnings
cargo clippy --fix --allow-dirty --manifest-path rust/Cargo.toml
# Then verify remaining warnings match allowed categories only:
cargo clippy --manifest-path rust/Cargo.toml -- \
  -D warnings \
  -A clippy::too_many_arguments \
  -A clippy::type_complexity
```

**What each step enforces:**

| Tool | Standard |
|---|---|
| `ruff` | PEP 8 + unused imports + undefined names |
| `black` | Consistent formatting (line-length=88) |
| `mypy` | Full type annotation compliance (0 errors — `_core.pyi` stub covers Rust extension) |
| `cargo clippy` | Rust idioms; `-A too_many_arguments` deferred (large refactor) |

**Baseline (as of 2026-03-02):** `pytest 244/4 ✅ · ruff 0 · black clean · mypy 0 errors ✅`

> [!IMPORTANT]
> **`_core.pyi` must be kept in sync with the Rust extension.**
> When adding or changing a `#[pyfunction]` in `rust/src/*.rs`, update
> `src/krewlyzer/_core.pyi` with the matching Python signature. Run
> `mypy src/krewlyzer/` to confirm 0 errors before committing.

---

## Code Quality Standards

- All modules have `__all__` exports
- All public functions have docstrings
- Use `logging` module (not print)
- Type hints throughout
- Use Rich for CLI progress/status (via `console.*`)

---

## QC Review Checklist (Per Feature)

**Before implementing any feature, review affected files:**

### 1. Logging
- [ ] All user-facing messages use proper logging levels (`debug`/`info`/`warning`/`error`)
- [ ] No `print()` statements for status (use `logger` or `console`)
- [ ] Critical operations have timing/performance logging
- [ ] Errors include context for debugging

### 2. Commenting
- [ ] All modules have docstrings describing purpose
- [ ] All public functions have docstrings with Args/Returns
- [ ] Complex logic has inline comments explaining "why"
- [ ] No commented-out code (delete or use version control)

### 3. Monitoring
- [ ] Pipeline operations have timing metrics
- [ ] Stats/counters for processed items (fragments, regions)
- [ ] Progress indicators for long operations
- [ ] Error counts and success rates tracked

### 4. No Code Duplication
- [ ] Common patterns extracted to helper functions
- [ ] Shared logic in `utils/` or `core/` modules
- [ ] Similar code across files consolidated

### 5. No Unused Code
- [ ] All imports are used
- [ ] No dead functions/methods
- [ ] No unused variables
- [ ] `__all__` exports match public API

### 6. Output Format Enforcement (Cross-Cutting Rule)

> **Every function that writes output files MUST accept and forward `output_format` and `compress` parameters.**

This is enforced by a static analysis gate. **Run it before committing any change that touches output writers:**

```bash
python scripts/check_output_format.py
# Exit 0 = all call sites correct
# Exit 1 = lists exact file:line violations
```

**Functions covered by the gate:**

| Function | File |
|----------|------|
| `write_table()` | `core/output_utils.py` |
| `process_region_entropy()` | `core/region_entropy_processor.py` |
| `compute_and_write_gc_factors()` | Rust via `_core.gc` |
| `write_extraction_outputs()` | `core/sample_processor.py` |
| `_core.run_unified_pipeline()` | Rust via `_core` |

**If you add a new output-writing function**, update `scripts/check_output_format.py` to add it to `RULES` before implementing callers. This puts the gate in place before any call sites exist.

**Known intentional out-of-scope exception:**
- `pon/build.py` → `process_sample()`: PON building always writes TSV (no user-facing format choice). Excluded in the script via `RULES` exceptions.

---

### 7. Cross-Cutting Parameter Discipline

When a parameter (like `output_format`/`compress`) must propagate to every call site of a function across many files:

1. **Write the enforcement script first**, before implementing call sites
2. **Run the script after every change**, not just at the end
3. **Add the script to the QC checklist** (this section)
4. **Never rely on manual read-through alone** for cross-cutting concerns — human review misses things at scale

This rule was introduced after discovering that manual verification passes consistently missed 1–3 call sites per pass across the codebase despite careful review.


## Feature Implementation Workflow

```
1. Create feature branch from develop
2. Review all affected files against QC checklist
3. Document findings in implementation plan
4. Fix QC issues FIRST before adding new features
5. Implement changes
6. Verify all tests pass
7. Re-check QC standards
8. Commit with detailed message
9. Create PR for review
```

---

## Testing

```bash
# Unit tests
pytest tests/unit/ -v

# Integration tests (requires test data)
pytest tests/integration/ -v

# End-to-end tests
pytest tests/e2e/ -v

# Run specific test file
pytest tests/unit/test_fsc.py -v
```

---

## Nextflow Development

```bash
# Run with test profile
nextflow run main.nf -profile test,docker

# Tower/Lilac runs
nextflow run main.nf -profile iris,lilac --samplesheet samples.csv
```
