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

---

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
