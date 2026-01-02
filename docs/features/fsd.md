# Fragment Size Distribution (FSD)

**Command**: `krewlyzer fsd`

## Purpose
Computes high-resolution (5bp bins) fragment length distributions per chromosome arm. Produces ML-ready features with log-ratio normalization and on/off-target split for panel data.

## Key Features

| Feature | Description |
|---------|-------------|
| **67 bins** | 5bp resolution from 65-400bp |
| **GC-weighted** | Corrects for sequencing bias |
| **On/off-target split** | Separate outputs for panel data |
| **Log-ratio normalization** | log2(sample / PoN_expected) |

## Usage

```bash
# WGS
krewlyzer fsd sample.bed.gz -o output_dir/ --genome hg19

# Panel (with target split via run-all)
krewlyzer run-all sample.bam -g ref.fa -o out/ \
    --target-regions panel_targets.bed
```

## Options

| Option | Short | Description |
|--------|-------|-------------|
| `--output` | `-o` | Output directory (required) |
| `--arms-file` | `-a` | Chromosome arms BED file |
| `--pon-model` | `-P` | PON model for log-ratio computation |
| `--genome` | `-G` | Genome build: hg19/hg38 |
| `--gc-correct` | | Apply GC bias correction (default: True) |

## Output Files

### `{sample}.FSD.tsv` (off-target / default)

| Column | Type | Description |
|--------|------|-------------|
| `region` | str | Chromosome arm (chr:start-end) |
| `65-69`, `70-74`, ... | float | Raw GC-weighted counts (67 bins) |
| `total` | float | Sum of all bins |
| `65-69_logR`, ... | float | log2(sample / PoN_expected) *(with PoN)* |
| `pon_stability` | float | 1 / (variance + k) *(with PoN)* |

### `{sample}.FSD.ontarget.tsv` (panel mode only)

Same schema, for fragments overlapping target regions (capture-biased).

## Biological Context

| Signal | Healthy | Cancer |
|--------|---------|--------|
| Modal peak | ~166bp | Left-shifted (~145bp) |
| 10bp periodicity | Clear nucleosome signal | May be altered |
| Arm variation | Minimal | Increased (correlates with CNAs) |

## Normalization Order

1. **GC-weighting** (Rust): Raw counts × correction factor
2. **Log-ratio** (Python): log2(sample / PoN_expected) when PoN provided

> **Note**: For panel data (MSK-ACCESS), use `--target-regions` in `run-all` to separate capture-biased on-target reads from unbiased off-target reads.
