# Changelog

All notable changes to this project will be documented in this file.

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
