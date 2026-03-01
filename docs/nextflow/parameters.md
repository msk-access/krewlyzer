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

## See Also

- [Samplesheet Format](samplesheet.md)
- [CLI Parameters](../cli/run-all.md)
- [Pipeline Outputs](outputs.md)
