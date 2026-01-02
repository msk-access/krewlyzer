# Fragment Size Ratio (FSR)

**Command**: `krewlyzer fsr`

## Purpose
Computes short/long fragment ratios for cancer biomarker analysis. Uses PoN-normalization **before** ratio calculation for accurate cross-sample comparison.

## Biological Context

The ratio of short to long fragments is a key indicator of tumor burden in cfDNA:

- **Short Fragments (65-149bp)**: Tumor DNA is typically shorter (~145bp) than healthy DNA (~166bp)
- **Long Fragments (221-400bp)**: Di/multi-nucleosomes, representing stable healthy chromatin

**Key Biomarker**: `short_long_ratio` - Higher ratio = higher probability of tumor DNA

## Usage
```bash
krewlyzer fsr sample.bed.gz -o output_dir/ --sample-name SAMPLE [options]
```

## Key Options

| Option | Short | Description |
|--------|-------|-------------|
| `--output` | `-o` | Output directory (required) |
| `--sample-name` | `-s` | Override sample name |
| `--bin-input` | `-b` | Custom bin file |
| `--pon-model` | `-P` | PON model for count normalization |
| `--genome` | `-G` | Genome build: hg19/hg38 |

## Output Format

Output: `{sample}.FSR.tsv`

| Column | Type | Description |
|--------|------|-------------|
| `region` | str | Genomic region (chr:start-end) |
| `short_count` | int | Raw short fragment count (65-149bp) |
| `long_count` | int | Raw long fragment count (221-400bp) |
| `total_count` | int | Total fragments (65-400bp) |
| `short_norm` | float | short / PoN_short_mean |
| `long_norm` | float | long / PoN_long_mean |
| `short_long_ratio` | float | **short_norm / long_norm** (primary biomarker) |
| `short_long_log2` | float | log2(short_long_ratio) for ML |
| `short_frac` | float | short / total |
| `long_frac` | float | long / total |

## Normalization Order (Critical)

> [!IMPORTANT]
> FSR normalizes counts to PoN **BEFORE** computing ratios:
> 1. **Normalize**: short_norm = short / PoN_short_mean
> 2. **Normalize**: long_norm = long / PoN_long_mean  
> 3. **THEN Ratio**: short_long_ratio = short_norm / long_norm

This ensures accurate cross-sample comparison by removing batch effects before ratio calculation.

## Clinical Interpretation

| Metric | Healthy Plasma | Cancer (ctDNA present) |
|--------|----------------|------------------------|
| Modal fragment size | ~166bp | Left-shifted (~145bp) |
| `short_long_ratio` | Low (baseline) | **Elevated** |
| `short_long_log2` | ~0 | **Positive** |

## Panel Data Support

For targeted panels, use `--target-regions` in `run-all`:

```bash
krewlyzer run-all sample.bam -g ref.fa -o out/ \
    --target-regions panel_targets.bed
```

> **Reference:** See [Citation & Scientific Background](../citation.md#fsr) for DELFI paper details.
