# Changelog

All notable changes to this project will be documented in this file.

## [0.5.0] - 2026-02-02

### Added

#### PON Framework
- **`build-pon` command**: Generate PON models from cohort samples with FSD/WPS/OCF/MDS baselines
- **Bundled PON assets**: Pre-computed xs1/xs2 PON models for `all_unique` and `duplex` variants
- **`--pon-variant` flag**: Select between `all_unique` (max coverage) or `duplex` (highest accuracy) PON models
- **Duplex tag warning**: mFSD warns when `--duplex` is used but no cD tags found in BAM

#### Panel Mode (MSK-ACCESS)
- **Assay-aware asset resolution**: Auto-load targets, anchors, PON via `-A/--assay xs1|xs2`
- **On/Off-target splitting**: All tools produce separate `.ontarget.tsv` outputs
- **Gene-level FSC**: New `{sample}.FSC.gene.tsv` and `{sample}.FSC.regions.e1only.tsv`
- **Bait edge masking**: WPS `--bait-padding` to remove capture artifacts
- **`--skip-target-regions`**: Force WGS mode for panel assays

#### Duplex Sequencing
- **`--duplex` flag (mFSD)**: Enable family size (cD tag) weighting for duplex BAMs
- **LLR scoring**: Log-likelihood ratio classification for cross-species support
- **GC-weighted mFSD**: 5 new GC-corrected mean fragment size columns

#### Region-Based Analysis
- **Region Entropy**: TFBS/ATAC dual-output architecture with per-cluster entropy
- **Region MDS**: Per-gene Motif Diversity Score with E1-only filtering
- **Rust backend**: New `region_entropy.rs` for high-performance calculation

#### Feature Enhancements
- **Jagged Index**: 1-mer End Motif analysis with C-end fraction
- **WPS v2.0**: Hierarchical stacking, extended anchors, panel-specific normalization
- **Output formats**: `--format` flag (tsv/parquet/json) for all tools
- **`--generate-json`**: Unified JSON export with all features

#### Infrastructure
- **Git LFS**: Large files (.gz, .parquet, .bed) tracked via git-lfs
- **BGZF compression**: All BED outputs use block-gzip format
- **GC references**: Pre-computed Parquet format for GRCh37/GRCh38
- **Unified processor**: Single-pass Rust pipeline for Extract→Motif→FSC→FSD→WPS→OCF

### Changed

#### Nextflow Pipeline
- **nf-core compliance**: Refactored to shared INPUT_CHECK subworkflow
- **KREWLYZER_RUNALL**: Unified module with 30+ output channels
- **`pon_variant` parameter**: New pipeline parameter (default: `all_unique`)

#### CLI
- **9 tools with `--pon-variant`**: run-all, fsc, fsd, fsr, wps, ocf, motif, region-entropy, region-mds
- **Centralized asset resolution**: `resolve_pon_model()` and `resolve_target_regions()`
- **Filter flags**: `--mapq`, `--minlen`, `--maxlen`, `--skip-duplicates`, `--require-proper-pair`
- **Parallel processing**: `--threads` option across all tools

#### Data & Assets
- **Terminology**: Renamed "blacklist" to "exclude_regions" for inclusive language
- **Directory structure**: `data/{type}/{genome}/{variant}/` for organized assets
- **Memory optimization**: Parallel sample processing with spawn context in build-pon

### Fixed
- mFSD: Filter discordant reads with extreme TLEN values
- mFSD: Fix verbose mode hanging by moving debug logging outside parallel loop
- Proper pair detection for legacy BAMs without proper pair flags
- GC correction for gene-level FSC in panel mode
- CIGAR handling improvements for INDELs and complex variants
- BAM extraction for v1 ACCESS and bgzip output format

### Documentation
- **11 docs updated**: PON Variant Selection across pon.md, pipeline.md, usage.md, feature docs
- **Release Guide**: New `docs/advanced/release-guide.md` for version release process
- **Math rendering**: LaTeX formulas in all feature documentation
- **Glossary**: New terms and definitions

### Tests
- **28 new test files**: Unit, integration, and e2e coverage
- **4 PON variant tests**: test_asset_resolution.py
- **Real data tests**: test_real_data.py for end-to-end validation

## [0.3.2] - 2025-12-18

### Fixed
- **CI Build**: Removed `gfortran` and `scikit-misc` - GC correction now fully in Rust
- **FSC GC Correction**: Added `--gc-correct/--no-gc-correct` and `--verbose` flags to FSC CLI
- **Python LOESS Removed**: Removed `gc_correct()` from `helpers.py` and `postprocess.py`

### Changed
- **Dockerfile**: Removed `gfortran` dependency (no longer needed)
- **GitHub Actions**: Removed `gfortran` from CI workflows
- **Nextflow Modules**: Updated all container versions to `0.3.2`

## [0.3.1] - 2025-12-18

### Added
- **Rust LOESS GC Correction**: New `rust/src/gc_correction.rs` module using the `lowess` crate
  - Per-fragment-type correction (short, intermediate, long)
  - Configurable LOESS parameters (fraction, iterations, delta)
- **FSR GC Correction**: `--gc-correct/--no-gc-correct` flag (default: **True**)
  - Uses corrected counts from Rust before calculating ratios
  - `--verbose` flag for detailed logging
- **WPS GC Correction**: `--gc-correct/--no-gc-correct` flag (default: **True**)
  - `--reference, -r`: Reference FASTA for computing region GC content
  - FASTA-based GC computation using rust-htslib::faidx
  - Graceful fallback if reference not provided

### Changed
- **FSC Rust Backend**: Added `count_fragments_gc_corrected()` function for integrated GC correction
- **WPS Rust Backend**: Updated `calculate_wps()` with `reference_path`, `gc_correct`, `verbose` parameters
- **lowess Dependency**: Updated from 0.3 to 0.4 for improved API

### Documentation
- Updated `docs/features/fsr.md` and `docs/features/wps.md` with GC correction options


## [0.3.0] - 2025-12-16

### Added
- **mFSD Variant Type Support**: Now handles all small variant types:
  - SNV (single nucleotide)
  - MNV (multi-nucleotide)
  - Insertion
  - Deletion
  - Complex (substitution + indel)
- **4-Way Fragment Classification**: Fragments classified as REF, ALT, NonREF, or N (low quality)
- **Comprehensive Statistics**: 6 pairwise KS comparisons (ALT-REF, ALT-NonREF, etc.)
- **Derived Metrics**: VAF_Proxy, Error_Rate, N_Rate, Size_Ratio, Quality_Score
- **Quality Flags**: ALT_Confidence (HIGH/LOW/NONE), KS_Valid
- **Distributions Output**: `--output-distributions` flag generates per-variant size distributions
- **Verbose Logging**: `--verbose` flag for debug-level logging with variant type breakdown
- **MRD Support**: Proper handling of low fragment counts (1-2) common in MRD settings

### Changed
- **BREAKING**: mFSD output format changed from 11 columns to 39 columns
- **Nextflow MFSD**: Distributions and verbose logging enabled by default
- **Fragment Counting**: Now counts fragments (R1 only) instead of reads to avoid double-counting

### Fixed
- **CIGAR Handling**: Improved sequence extraction for INDELs and complex variants

## [0.2.3] - 2025-12-16

### Changed

- **Nextflow Pipeline**: Added `FILTER_MAF` module for multi-sample MAF filtering:
  - New `maf` and `single_sample_maf` columns in samplesheet
  - Filters MAF by `Tumor_Sample_Barcode` matching sample ID (regex: `.*sample_id.*`)
  - `single_sample_maf=true` bypasses filtering for per-sample MAFs
  - Skips MFSD when filtered MAF has zero variants (with warning)
  - Memory-efficient streaming for large MAF files (100s of MBs)
  - Outputs filtered MAF + mFSD results
- **Nextflow Modules**: Updated all container versions to `0.2.3`

## [0.2.2] - 2025-12-15

### Changed
- **Project Structure**: Migrated to recommended maturin "src layout" for better Python/Rust separation:
  - Python code moved from `krewlyzer/` to `src/krewlyzer/`
  - Rust code moved from `krewlyzer-core/` to `rust/`
  - Single `rust/Cargo.toml` (removed duplicate root `Cargo.toml`)
  - Rust module now imports as `krewlyzer._core` (private)
- **Dockerfile**: Rewritten as multi-stage build for smaller image size (~200MB vs ~1GB); amd64 only
- **OCI Labels**: Added container metadata for GitHub Container Registry

### Fixed
- **Distribution Compatibility**: Updated release workflow to build `manylinux_2_28` wheels (compatible with RHEL 8+, AlmaLinux 8, etc.)
- **Source Builds**: Included `clang` and `llvm-devel` in build environment for `bindgen`/`hts-sys`
- **Docker Build**: Added `patchelf` to maturin installation for Linux wheel building
- **Test Imports**: Updated test files to use new `krewlyzer._core` import path
- **CI OpenSSL**: Use system OpenSSL (`OPENSSL_NO_VENDOR=1`) instead of building from source
- **CI Python Versions**: Build wheels for Python 3.10, 3.11, 3.12 in release workflow
- **Documentation**: Updated logo paths for new `src/` layout

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
