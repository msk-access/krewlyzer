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

> **Agent rule:** Always run the full lint suite before committing. CI enforces exactly these
> same commands as a blocking `lint` job. Fix all failures before pushing.

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
