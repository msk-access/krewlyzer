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

## Key Rust Modules

| Module | Purpose |
|--------|---------|
| `extract_motif.rs` | BAM parsing, fragment extraction |
| `pipeline.rs` | Single-pass FSC/FSD/WPS/OCF |
| `gc_correction.rs` | LOESS GC bias correction |
| `fsc.rs` | Fragment size coverage + gene aggregation |
| `fsd.rs` | Size distribution + PON log-ratio |
| `wps.rs` | Windowed protection score + PON z-score |
| `ocf.rs` | Orientation-aware fragmentation + PON z-score |
| `region_entropy.rs` | TFBS/ATAC entropy + PON z-score |

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
