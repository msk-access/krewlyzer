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
| `pon_model.rs` | `FsdBaseline` structs for `fsd.rs`. Its `PonModel::load` has **no callers** and downcasts `table`/`region_id` to `StringArray`, which every shipped PON stores as `large_string` -- so it would read zero rows if anything did call it. Do not reach for it; each reader loads what it needs. |
| `pon_builder.rs` | PON model construction |
| `filters.rs` | Fragment filtering logic |

## PON Z-Score Pattern

Every PON z-score reads its baseline directly from the Parquet, taking the
*path* rather than a loaded model -- the Python side never materialises a PON
just to hand it over. The exact signatures differ per module; read the one you
are calling. `wps.rs::apply_pon_zscore` is the current shape:

```rust
pub fn apply_pon_zscore(
    wps_path: PathBuf,
    pon_parquet_path: PathBuf,
    output_path: PathBuf,
    baseline_table: &str,  // "wps_baseline", "wps_baseline_panel"
    column: &str,          // "wps_nuc", "wps_tf"
) -> PyResult<usize>       // anchors scored; 0 means nothing was written
```

**Returning 0 must mean writing nothing.** The Python caller degrades to the
raw table on a zero, so a scorer that half-writes on a missing baseline
destroys the product it could not improve (invariant #2).

### Two traps when reading a PON parquet

**String width.** The builder writes `large_string`, not `string`; a bare
`downcast_ref::<StringArray>()` returns `None` on every shipped PON and yields
an empty baseline -- which is a legitimate state, so it degrades silently
rather than raising. Handle both widths, or use the row API
(`SerializedFileReader::get_row_iter`), which is logical-typed and sees both.
The row API reads *every* column of every row, so it is only viable for the
small blocks; `wps.rs` projects by root index instead, because the wide PON
table is ~120 columns across ~130k rows.

**Vector element type.** The Python builder writes `list<double>`; the Rust
builder's own fixtures write `list<float>`. A reader that accepts one silently
finds no anchors in the other. `wps.rs` accepts both and logs the type it
could not read rather than returning nothing.

### Which language

Scoring belongs in Rust, but the *design* of a statistic may not start there.
The WPS scoring shipped in Python through 0.8.x on purpose: three of its
decisions -- the `log1p` amplitude, the Fisher transform, and the absent phase
baseline -- were settled by measurement on a real cohort and changed several
times. It moved to `wps.rs` in 0.9.0 once settled, measured at 9.5 s against
~17 min for the Python, on the same 76,595-anchor output. The Python survives
as a frozen equivalence oracle in
`tests/unit/test_rust_python_equivalence.py`; the port is bug-for-bug, so a
divergence there is a finding rather than a sync.

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

---

## Gene BED assets (`scripts/build_gene_bed.py`)

Build-time, not runtime. The GTF is parsed once when the asset is regenerated;
nothing in `src/` or `rust/` reads a GTF, and no GENCODE dependency exists at
run time.

**Canonical transcript policy**, decided here and nowhere else:

1. MANE Select
2. Ensembl canonical
3. Longest CDS

The build **fails** if the GTF carries no MANE tags. Ensembl's GRCh37 line is
frozen at release 87 (2016) and predates MANE, so it would silently fall
through to longest-CDS for every gene. GRCh37 assets must come from GENCODE's
`lift37`.

**Anything derived from strand or exon order belongs in the asset, not the
runtime.** `is_e1` is precomputed for exactly this reason: deriving "first
exon" from coordinates at run time gets every minus-strand gene wrong, and did
— in two separate implementations that disagreed with each other. The runtime
should read the column, not re-derive the answer.

Panel gene symbols are matched through an explicit alias table
(`LEGACY_GENE_ALIASES`); three MSK-ACCESS symbols were renamed by HGNC and
match nothing in a current GTF. An unmatched gene is an error, not a warning,
because the failure mode is silent loss of annotation.
