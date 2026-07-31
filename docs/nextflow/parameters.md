# Nextflow Parameters

All parameters for the Krewlyzer Nextflow pipeline. See `nextflow.config` for defaults.

## Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--samplesheet` | CSV with sample information ([format](samplesheet.md)) |
| `--ref` | Reference genome FASTA (indexed) |

## General Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `./results` | Output directory |
| `--asset_dir` | | Base directory for PON/targets (enables assay resolution) |
| `--targets` | | Global target BED (fallback if not in samplesheet) |
| `--genome` | `hg19` | Genome build (`hg19` or `hg38`) |
| `--threads` | `8` | Threads per process |
| `--verbose` | `false` | Enable verbose logging |

## Mode Selection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--use_runall` | `true` | `true` = unified `run-all` (default), `false` = tool-level subworkflow |

## PON Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--pon_model` | | Global PON model path (fallback) |
| `--pon_variant` | `all_unique` | PON variant: `all_unique` or `duplex` |
| `--skip_pon` | `false` | Skip PON z-score normalization |

## Panel & WPS Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bait_padding` | `50` | Bait edge padding for WPS (bp) |
| `--wps_anchors` | | Custom WPS anchors BED (auto-loaded from assay if not set) |
| `--wps_background` | | Custom WPS background Alu BED (auto-loaded if not set) |

## Filter & Algorithm Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mapq` | `20` | Minimum mapping quality |
| `--minlen` | `65` | Minimum fragment length |
| `--maxlen` | `1000` | Maximum fragment length |
| `--min_baseq` | `20` | Minimum base quality for mFSD variant calling |
| `--duplex` | `true` | Enable duplex weighting for mFSD (graceful fallback) |
| `--skip_duplicates` | `true` | Skip duplicate reads |
| `--require_proper_pair` | `false` | Require proper pairs (disable for duplex BAMs) |
| `--exclude_regions` | | BED file of regions to exclude |

## Feature Toggles

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--no_tfbs` | `false` | Disable TFBS region entropy |
| `--no_atac` | `false` | Disable ATAC region entropy |
| `--skip_target_regions` | `false` | Force WGS mode (ignore bundled targets) |
| `--disable_e1_aggregation` | `false` | Skip E1-only FSC aggregation |
| `--region_mds_e1_only` | `false` | Run region-MDS on E1 (first exon) only |

## Output Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--generate_json` | `true` | Generate unified `features.json` for ML pipelines |
| `--output_format` | `parquet` | Feature output format: `tsv`, `parquet`, or `both`. **Default changed from `tsv` in 0.9.0.** WPS outputs are always Parquet regardless of this setting. |
| `--compress_tsv` | `false` | Gzip-compress all TSV outputs (`.tsv.gz`). Applies only when `output_format` is `tsv` or `both`. Maps to the `--compress` flag in the Python CLI. |
| `--strict_validation` | `false` | Fail a sample when its output violates the contract. The report and fingerprint are written either way; this only controls whether a violation stops the run. |
| `--validate_min_samples` | `3` | Minimum samples before cross-sample degeneracy checks are meaningful. Below this the cohort step reports SKIP, never PASS. |
| `--gc_correct` | `true` | Apply GC bias correction during extraction. |
| `--queue_size` | `100` | Maximum concurrent executor jobs; also derives `FILTER_MAF` maxForks. |

!!! warning "Selecting `tsv` produces a cohort the downstream consumer cannot read"
    kreview reads **Parquet only**. A tsv-only run additionally **skips output
    validation**, because the contract describes the Parquet surface.

    Both failures are silent: every reader downstream swallows exceptions and
    yields an empty feature dict, so a cohort produced with the default simply
    does not appear — no error, no warning, no missing-file complaint.

    This is why the default changed to `parquet` in 0.9.0. Set `tsv` or
    `both` only when you need the text tables; the pipeline logs a warning at
    start whenever the selection excludes Parquet.

## See Also

- [Samplesheet Format](samplesheet.md)
- [CLI Parameters](../cli/run-all.md)
- [Pipeline Outputs](outputs.md)
