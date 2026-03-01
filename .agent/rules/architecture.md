---
description: Krewlyzer architecture reference — Rust/Python boundary, module table, PON patterns, parallelization
alwaysApply: true
---

# Krewlyzer Architecture Reference

## Rust/Python Boundary Rules

| Rust ✓ | Python ✓ |
|--------|----------|
| File I/O (BAM, BED, BED.gz, TSV, Parquet) | CLI / Typer commands |
| Row-level computation | Asset path resolution |
| Loops over >1000 rows | PON model building (metadata) |
| GC correction (LOESS) | Workflow orchestration |
| PON z-score / log-ratio | High-level coordination |
| Interval tree lookups | |
| Fragment counting/aggregation | |

## Rust Modules (19 total)

| Module | Purpose |
|--------|---------|
| `lib.rs` | PyO3 module exports and thread config |
| `pipeline.rs` | Single-pass FSC/FSD/WPS/OCF |
| `engine.rs` | Core engine utilities |
| `bed.rs` | BGZF/gzip BED reader |
| `extract_motif.rs` | BAM parsing, fragment extraction |
| `motif_utils.rs` | Shared 4-mer encoding, MDS calculation, GC content |
| `fsc.rs` | Fragment size coverage + gene aggregation |
| `fsd.rs` | Size distribution + PON log-ratio |
| `wps.rs` | Windowed protection score + PON z-score |
| `ocf.rs` | Orientation-aware fragmentation + PON z-score |
| `mfsd.rs` | Mutant fragment size distribution |
| `region_mds.rs` | Per-region MDS at exon/gene level |
| `region_entropy.rs` | TFBS/ATAC entropy + PON z-score |
| `uxm.rs` | Fragment-level methylation |
| `gc_correction.rs` | LOESS GC bias correction |
| `gc_reference.rs` | Pre-computed GC reference generation |
| `pon_model.rs` | PON model loading and hybrid correction |
| `pon_builder.rs` | PON model construction |
| `filters.rs` | Fragment filtering logic |

## PON Z-Score Pattern

All PON z-score functions read baselines directly from Parquet:

```rust
// Rust signature pattern
pub fn apply_pon_zscore(
    input_path: PathBuf,
    output_path: PathBuf,
    pon_parquet: Option<PathBuf>,
    baseline_table: &str,  // "tfbs_baseline", "atac_baseline", etc.
) -> PyResult<bool>
```

```python
# Python call pattern
_core.module.apply_pon_zscore(input, output, pon_parquet, "baseline_table")
```

## Chromosome Normalization

Always strip `chr` prefix before `chrom_map.get_id()`:

```rust
let chrom_norm = chrom_str.trim_start_matches("chr");
let chrom_id = chrom_map.get_id(chrom_norm);
```

## Gzip Support

`bed::get_reader()` uses flate2 - works with both standard gzip and BGZF.

---

## Parallel Processing Architecture

### Tier-1: Sample-Level Parallelization (`build-pon -P`)

For PON building, samples are processed in parallel using Python `ProcessPoolExecutor`:

```python
# build.py
with ProcessPoolExecutor(max_workers=parallel_samples) as executor:
    futures = {executor.submit(process_sample, ...): sample for sample in samples}
```

### Tier-2: Region-Level Parallelization (Rust)

Within each sample, regions are processed in parallel using Rayon:

| Module | Pattern |
|--------|---------|
| `wps.rs` | `.par_iter()` over anchor regions |
| `region_entropy.rs` | `.par_iter()` over TFBS/ATAC labels |
| `fsc.rs` | `.par_iter()` for gene aggregation |
| `region_mds.rs` | `.par_iter()` over gene regions |

```rust
// Example: WPS region parallelization
regions.par_iter().for_each(|region| {
    // Process each region independently
});
```

### Thread Configuration

```python
_core.configure_threads(num_threads)  # Sets Rayon thread pool size
```

---

## PON Baseline Requirements

| Baseline | MIN_SAMPLES | Notes |
|----------|:-----------:|-------|
| **fsc_gene_baseline** | 3 | Per-gene normalized depth |
| **fsc_region_baseline** | 3 | Per-exon normalized depth |
| **region_mds_baseline** | 3 | Per-gene MDS |
| Other baselines | 1+ | No minimum |

---

## FSC Gene Aggregation (`fsc.rs`)

| Function | Purpose |
|----------|---------|
| `aggregate_by_gene()` | Counts fragments per gene/region with GC correction |
| Output | `FSC.gene.tsv`, `FSC.regions.tsv` |
| Uses | On-target GC factors for panel mode |
