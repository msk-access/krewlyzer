# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Fixed
- **mFSD Silent Error Swallowing**: Replaced `Err(_) => continue` in BAM record
  iterator with a logged error breaker (max 1000 consecutive errors). Previously,
  corrupt BAM regions could cause infinite silent loops.
- **FILTER_MAF Comment Lines**: Stripped `#` comment lines from filtered MAF output
  in both multi-sample and single-sample modes. Downstream Rust parser already
  skipped them, but output files were unnecessarily large and confusing.
- **mFSD 0-Variant Guard**: Added early exit at both Python (`wrapper.py`) and Rust
  layers when MAF has 0 data lines. Produces header-only TSV instead of attempting
  BAM access. Prevents unnecessary resource allocation for samples without variants.
- **mFSD GC Correction Fallback**: When reference FASTA is unavailable or GC lookup
  fails for a region, GC correction is now skipped (weight=1.0) instead of silently
  using a hardcoded 50% GC to look up correction factors.

### Added
- **mFSD BAM I/O Diagnostics**: Per-variant timing, BAM open/fetch latency logging,
  record counts, and slow-variant warnings (>30s). Enables production debugging of
  the 17% job failure rate on IRIS HPC.
- **mFSD Header Constant**: Extracted 46-column TSV header into `MFSD_HEADER` module
  constant, eliminating duplication between normal output and 0-variant early-exit.

### Tests
- Added `test_mfsd_zero_variants` — verifies 0-variant input produces header-only TSV.
- Added `test_mfsd_maf_with_comment_lines` — verifies MAFs with `#` comment headers
  are parsed correctly.

## [0.8.2] - 2026-03-26

### Fixed
- **Region MDS Output Collision**: Fixed `Path.with_suffix()` compound-extension bug that
  caused `MDS.exon.tsv` and `MDS.gene.tsv` to both resolve to `MDS.tsv`, silently
  destroying per-exon and per-gene data. Affects standalone `region-mds`, `motif`, and
  `run-all` commands.
- **Motif Tracking Paths**: Fixed 5 additional `with_suffix()` bugs in motif output
  tracking (EndMotif, BreakPointMotif, MDS, MDS.ontarget, EndMotif1mer) that caused
  path collisions and silent PON z-score normalization failures.
- **Gene Format Conversion**: Region MDS gene output now correctly undergoes format
  conversion (Parquet/gzip) when no PON is provided.
- **MDS Z-score Logging**: Z-score append failures upgraded from `debug` to `warning`
  level to surface issues in production runs.

### Added
- **Compound Extension Tests**: 57 new parametrized tests (`test_compound_extension.py`)
  covering all 13 compound base names across TSV, Parquet, both formats, and roundtrip.
- **Developer Guide**: "Known Gotchas" section documenting the `Path.with_suffix()` anti-
  pattern with safe alternative and complete table of affected names.

### Changed
- **Nextflow Outputs**: Added missing Parquet emit declarations for TFBS sync (2), ATAC
  sync (2), ATAC/TFBS ontarget sync (2), and fsc_counts (4 — TSV+Parquet, off/on-target).
- **SLURM Script**: Tuned head process memory (16G→32G), queue size (100→200), added
  `--output_format both --compress_tsv true --verbose true` for 14K+ sample production runs.

## [0.8.1] - 2026-03-24

### Fixed
- **WPS Panel GC Correction**: Panel WPS now uses on-target GC correction factors instead
  of off-target. Panel anchors overlap capture regions, making on-target correction more
  accurate. Falls back to off-target factors when on-target unavailable.
- **Rust WPS Background**: Fixed `coitrees` metadata type mismatch in WPS background consumer
  (`wps.rs:1684`) — uses `.clone()` for cross-platform compatibility (returns `&usize` on
  macOS but `usize` on CI Linux).

### Changed
- **Dead Code Cleanup**: Removed 11 unused Python z-score functions from `pon_integration.py`
  (all replaced by Rust equivalents). Module reduced from 448 to 98 lines, 14 to 3 functions.
  Remaining: `load_pon_model`, `compute_nrl_zscore`, `compute_periodicity_zscore`.
- **GC Factor Resolution**: `gc_str`/`gc_ontarget_str` path resolution hoisted to a shared
  section in `unified_processor.py`, eliminating duplication between panel WPS and TFBS/ATAC.

### Added
- **build-pon Logging**: OCF on-target/off-target baseline status now logged during PON build.
- **sbatch Script**: `scripts/build_pon_unfiltered.sh` for building PON from high-coverage
  unfiltered BAMs on SLURM clusters.

### Data
- **PON Models Updated**: Rebuilt xs1/xs2 PON models for both `all_unique` and `duplex`
  variants with krewlyzer 0.8.x.

## [0.8.0] - 2026-03-17

### Added
- **On-target PON Z-scores**: Panel mode now computes on-target/off-target PON baselines
  for MDS, OCF, and FSD features, providing clinically-scoped z-scores in panel assays.
- **FSR On-target Output**: FSR now emits a separate `.FSR.ontarget.tsv` file in panel mode.
- **FSR Real Genomic Coordinates**: Region labels now reflect true genomic window coordinates
  instead of internal indices.

### Fixed
- **Output Format / Compress**: All output files now correctly respect `--output-format` and
  `--compress` flags (gzip compression, path handling, GC correction factors loading).
- **Rust GC Correction**: Added missing `PathBuf` import; fixed path handling for correction
  factor files.
- **Rust Clippy**: Replaced manual string strip with `strip_suffix` (`manual_strip` lint).

### Documentation
- Corrected stale FSR column names in `concepts.md` and `json-output.md`.
- Updated FSR on-target output docs and window/panel-mode descriptions.
- Documented `_core.pyi` stub maintenance requirements.

### CI
- **GitHub Actions Node.js 24**: Bumped `actions/checkout` v4→v5, `actions/setup-python` v5→v6,
  `actions/cache` v4→v5 to address Node.js 20 deprecation (enforced June 2, 2026).

## [0.7.0] - 2026-03-02

### Added
- **Configurable Output Formats**: `--output-format tsv|parquet|both` (default: `tsv`) controls
  all tabular feature outputs. `--compress` gzip-compresses TSV outputs (`.tsv.gz`).
  WPS outputs remain always-Parquet regardless of setting.

### Fixed
- **build-pon Intermediate Files**: Explicitly force `output_format="tsv"` and `compress=False`
  in all `process_sample()` calls within `build-pon` (both parallel and sequential paths).
  Prevents silent failure if default output format changes — intermediate files are internal
  scratch consumed by `pd.read_csv(sep="\t")`.
- **Feature Serializer**: Include `mds_z` in JSON output for the `from_outputs()` code path.
- **OCF Base File**: `OCF.tsv` = all reads (on + off combined), not off-target.
  `OCF.offtarget.tsv` is the true panel off-target score. Corrected in docs and code comments.
- **Rust wps.rs**: Remove erroneous `*` dereference on `node.metadata` (E0614 — `usize` is Copy).
- **Rust gc_correction.rs**: Prefix unused `valid_regions_path` param with `_` to silence
  compiler warning; parameter retained for API symmetry.
- **MkDocs Snippets**: Fix `--8<-- "CHANGELOG.md"` / `--8<-- "CONTRIBUTING.md"` broken includes
  by changing `pymdownx.snippets.base_path` from scalar `docs` to list `['.', 'docs']` so
  repo-root files resolve without path traversal blocked by CI deploy sandbox.

### Documentation
- **Post-0.6.0 Docs Audit** (12 files, 7 issue categories):
  - Fixed broken `--8<-- "../CHANGELOG.md"` / `--8<-- "../CONTRIBUTING.md"` snippet includes
  - Corrected outdated "no global `--output-format` flag" note in `cli/run-all.md`
  - Updated `metadata.json` → `metadata.tsv` across 7 files (8 references total)
  - Added WPS always-Parquet exception note to `reference/output-files.md`
  - Added build-pon intermediate TSV format note to `guides/building-pon.md`
  - Updated test count (248 → 244 + 4 skipped) and added CI lint steps to `developer-guide.md`
  - Added `--output_format` and `--compress` parameters to `nextflow/parameters.md`
- **Output Format Options section**: New section in `reference/output-files.md` documenting
  `--output-format`, `--compress`, WPS always-Parquet exception, and `--generate-json`
- **OCF Variant Clarification**: Added 3-variant table and note block in `reference/output-files.md`
  explaining `OCF.tsv` (all reads) vs `OCF.ontarget.tsv` vs `OCF.offtarget.tsv`
- **docs/index.md**: Replaced `:latest` Docker tag with explicit `:0.7.0` per release policy

### CI
- **Lint Job**: Added parallel `lint` job (`ruff · black · mypy · cargo clippy -- -D warnings`)
  running alongside tests on all push/PR events

## [0.6.0] - 2026-02-28

### Added
- **mFSD Base Quality Filtering**: `--min-baseq` / `-Q` (default 20) gates variant evidence by base quality
- **mFSD GC Correction**: Rust-native LOESS GC bias correction for variant fragment size distributions
- **mFSD Duplex Weighting**: Proper consensus fragment handling via `--duplex`
- **Region MDS `--sample-name`**: Consistent output naming without post-hoc rename
- **Feature Serializer**: Auto-load `fsc_counts`, `region_mds`, `uxm` in `from_outputs()`
- **IRIS Batch Submission**: `scripts/run_krewlyzer_iris.sh` for SLURM/IRIS cluster runs with `--generate_json`
- **nf-core Institutional Configs**: `custom_config_base` param and IRIS profile support
- **Versioned Documentation**: Implemented `mike` for dev/stable doc versions
- **Nextflow mfsd Module**: Full standalone params (`--reference`, `--correction-factors`, `--mapq`, `--minlen`, `--maxlen`, `--min-baseq`, `--duplex`, `--no-skip-duplicates`)
- **Nextflow runall**: `fsc_counts.tsv` output declaration, `--min-baseq` wired
- **mFSD Integration Tests**: 161 lines of new test coverage

### Fixed
- **mFSD MAF Parsing**: Header-based column lookup (fixes column-index mismatch with different MAF flavors)
- **Nextflow FILTER_MAF**: Complete overhaul — eliminated join operator blocking, replaced regex with substring match, fixed SyntaxError in versions.yml, dynamic maxForks for SLURM
- **Nextflow Workflow Streaming**: Fixed RUNALL blocking from `remainder:true`, `failOnMismatch`, channel round-robin; used `multiMap` for proper routing
- **Nextflow RUNALL Outputs**: Added 14 output declarations, fixed BreakPointMotif casing, explicit publishDir
- **Region MDS Nextflow**: `--sample-name` replaces `mv` workaround
- **Nextflow Config**: Executor queueSize placement, `-qs` CLI flag, global publishDir removal
- **WPS CLI Tests**: Fixed `--input` flag (was positional arg) — recovered 2 skipped tests
- **Pandas FutureWarning**: Fixed `pd.concat` with all-NA columns in PON test fixture

### Changed
- **Code Quality**: Black formatted 71 files, ruff fixed 129 lint errors, cargo clippy applied
- **Ruff Config**: Added `[tool.ruff]` to `pyproject.toml` with documented E402/F821 ignores
- **Agent Config**: `.agent/` → `.agent/rules/` with `alwaysApply` frontmatter

### Documentation
- **45-item Audit**: Corrections across 25 doc files including `.csv`→`.tsv` (7 files), `.WPS.tsv.gz`→`.parquet` (3 files), phantom `--output-format` removed, Docker versions→`X.Y.Z`, parameters.md 12→28, outputs.md 14→41, JSON schema corrected, developer guide Rust table 10→19, architecture pipeline signature updated
- **PDF Embedding**: Fixed rendering with mkdocs-pdf plugin

## [0.5.3] - 2026-02-06

### Documentation
- **Admonition Formatting**: Converted blockquotes to proper MkDocs admonitions across 10+ docs files for consistent styling
- **Table Rendering**: Fixed tables not rendering by adding required blank lines before table headers
- **LaTeX Formulas**: Fixed underscore escaping in math blocks (mFSD, FSR, WPS formulas)
- **Abbreviations**: Added glossary with auto-append for consistent acronym tooltips (cfDNA, WPS, etc.)
- **Docs CI**: Added PR validation trigger with `mkdocs build --strict` to catch issues before deployment

## [0.5.2] - 2026-02-05

### Added
- **Dual BAM Support (mFSD)**: New `--mfsd-bam` parameter in `run-all` to use a dedicated duplex BAM for mFSD while other features use the main all_unique BAM
  - Auto-duplex detection from filename patterns (`-duplex`, `_duplex`, `.duplex`)
  - Startup banner shows mFSD BAM path when specified
  - Nextflow: Added `mfsd_bam` column to samplesheet schema

### Fixed
- **Correction Factors Loading**: Fixed delimiter detection for `.correction_factors.tsv` files (CSV format with TSV extension)

### Documentation
- **Release Guide**: Added version format note (use `0.5.2` not `v0.5.2` in code)
- **mFSD docs**: Updated with dual BAM workflow examples
- **Samplesheet docs**: Added `mfsd_bam` column documentation
- **Footer**: Added site footer with author attribution and Antigravity acknowledgment
- **Theme**: Changed color scheme to blood-red (#ef5552) to match "krew" (blood) branding

## [0.5.1] - 2026-02-04

### Fixed
- **Docker Image**: Data folder now bundled in container (was missing in 0.5.0)
- **CI Tests**: Tests requiring bundled data now skip in PyPI installs via `@requires_data` marker

### Changed
- **Docker Tags**: Use versioned tags only (e.g., `:0.5.1`); no `:latest` tag published

### Documentation
- **Installation**: Added Singularity/Apptainer section for HPC clusters
- **Installation**: Clarified Docker uses release tags, not `:latest`
- **Nextflow Examples**: Added IRIS HPC cluster profile example (`-profile iris`)
- **Testing Guide**: Added Data Availability section explaining PyPI vs source
- **Structure**: Moved Changelog to Development section; removed duplicate files
- **Mermaid Diagrams**: Upgraded to `mkdocs-mermaid2` plugin for theme-aware rendering
- **Removed**: Full-site PDF export (unreliable rendering)

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
