# Fragment Size Coverage (FSC)

**Command**: `krewlyzer fsc`

## Purpose
Computes GC-corrected coverage of cfDNA fragments in 5 biologically-meaningful size channels per genomic bin. Designed for ML feature extraction in cancer detection.

## 5-Channel Fragment Size Categories

| Channel | Size Range | Biological Meaning |
|---------|------------|-------------------|
| **ultra_short** | 65-100bp | Di-nucleosomal debris, early apoptosis |
| **core_short** | 101-149bp | Sub-nucleosomal, specific chromatin states |
| **mono_nucl** | 150-220bp | Mono-nucleosomal (classic cfDNA peak) |
| **di_nucl** | 221-260bp | Di-nucleosomal, transitional |
| **long** | 261-400bp | Multi-nucleosomal, necrosis-associated |

> **Note**: These 5 channels are **non-overlapping** to avoid double-counting in ML models.

## Usage
```bash
krewlyzer fsc sample.bed.gz -o output_dir/ --sample-name SAMPLE [options]
```

## Key Options

| Option | Short | Description |
|--------|-------|-------------|
| `--output` | `-o` | Output directory (required) |
| `--sample-name` | `-s` | Override sample name |
| `--bin-input` | `-b` | Custom bin file (default: 100kb bins) |
| `--pon-model` | `-P` | PON model for log2 ratio normalization |
| `--genome` | `-G` | Genome build: hg19/hg38 |
| `--gc-correct` | | Apply GC bias correction (default: True) |

## Output Format

Output: `{sample}.FSC.tsv`

### Base Columns (always present)

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome |
| `start` | int | Window start (0-based) |
| `end` | int | Window end |
| `ultra_short` | float | GC-weighted count (65-100bp) |
| `core_short` | float | GC-weighted count (101-149bp) |
| `mono_nucl` | float | GC-weighted count (150-220bp) |
| `di_nucl` | float | GC-weighted count (221-260bp) |
| `long` | float | GC-weighted count (261-400bp) |
| `total` | float | GC-weighted total (65-400bp) |

### PoN Columns (when `--pon-model` provided)

| Column | Type | Description |
|--------|------|-------------|
| `*_log2` | float | log2(channel / PoN_mean) |
| `*_reliability` | float | 1 / (PoN_variance + k) |

## Normalization Order

1. **GC-weighting** (Rust): Raw counts × correction factor per (length, GC) bin
2. **PoN log-ratio** (Python): log2(sample / PoN mean) when PoN model provided

## Panel Data Support

For targeted panels (e.g., MSK-ACCESS), use `--target-regions` in `run-all` to compute GC correction from off-target reads only:

```bash
krewlyzer run-all sample.bam -g ref.fa -o out/ \
    --target-regions panel_targets.bed
```
