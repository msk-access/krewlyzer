# Changelog

All notable changes to this project will be documented in this file.

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
