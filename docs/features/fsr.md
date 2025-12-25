# Fragment Size Ratio (FSR)

**Command**: `krewlyzer fsr`

## Purpose
Computes ratios of fragment size classes per genomic window. Unlike FSC, this metric is self-normalizing (ratio) and focuses on the **proportion** of fragment sizes, which is a powerful biomarker for tumor fraction estimation.

## Biological Context
The ratio of short to long fragments is a key indicator of tumor burden in cfDNA ("fragmentomics"). See [DELFI method](../citation.md#fsr) for details.

- **Short Fragments (65-149bp)**: Tumor DNA is typically shorter (~145bp) than healthy DNA (~166bp).
- **Ultra-short Fragments (65-100bp)**: Associated with transcription factor binding sites and open chromatin.
- **Long Fragments (261-399bp)**: Often di-nucleosomes, representing stable, healthy chromatin (Leukocytes).

**Key Biomarkers:**
- **Short/Long Ratio**: The primary metric. Higher ratio = higher probability of tumor DNA.
- **Ultra-short Ratio**: Indicates active gene regulation/transcription factor binding.
- **Nucleosome Footprints**: The 10bp precision of the ranges (150bp vs 167bp) helps separate mono-nucleosomes.

## Usage
```bash
krewlyzer fsr sample.bed.gz -o output_dir/ --sample-name SAMPLE [options]
```

## Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--output` | `-o` | PATH | *required* | Output directory |
| `--sample-name` | `-s` | TEXT | | Override sample name |
| `--bin-input` | `-b` | PATH | | Bin file (default: hg19 100kb bins) |
| `--pon-model` | `-P` | PATH | | PON model for hybrid GC correction |
| `--windows` | `-w` | INT | 100000 | Window size |
| `--continue-n` | `-c` | INT | 50 | Consecutive window count |
| `--genome` | `-G` | TEXT | hg19 | Genome build (hg19/hg38) |
| `--gc-correct` | | FLAG | True | Apply GC bias correction |
| `--verbose` | `-v` | FLAG | | Enable verbose logging |
| `--threads` | `-t` | INT | 0 | Number of threads (0=all) |


## Output Format

Output: `{sample}.FSR.tsv`

| Column | Description | Biological Relevance |
|--------|-------------|----------------------|
| `region` | Genomic region | |
| `ultra_short_count` | Count 65-100bp | TF footprints |
| `short_count` | Count 65-149bp | Tumor enriched |
| `inter_count` | Count 151-259bp | Mono-nucleosome |
| `long_count` | Count 261-399bp | Healthy / Di-nucleosome |
| `total_count` | Count 65-399bp | Total fragments |
| `short_ratio` | Short / Total | Tumor fraction proxy |
| `inter_ratio` | Intermediate / Total | |
| `long_ratio` | Long / Total | Healthy fraction proxy |
| `short_long_ratio` | Short / Long | **Primary Cancer Biomarker** |
| `ultra_short_ratio` | Ultra-short / Total | TF activity indicator |

## Clinical Interpretation

### Healthy vs Cancer

| Metric | Healthy Plasma | Cancer (ctDNA present) |
|--------|----------------|------------------------|
| Modal fragment size | ~166bp | Left-shifted (~145bp) |
| `short_long_ratio` | Low (baseline) | **Elevated** |
| Genome-wide ratio variability | Minimal | Increased aberrations |

### Interpretation Guide
- **High `short_long_ratio`**: Indicates higher tumor fraction (ctDNA) or open chromatin.
- **High `ultra_short_ratio`**: Indicates regions of high transcription factor activity.
- **High `long_ratio`**: Indicates stable, nucleosomal DNA (usually healthy background).

### Clinical Utility
- **Tumor fraction proxy**: Compare `short_long_ratio` to healthy baseline
- **Treatment monitoring**: Track ratio changes over time
- **Performance**: DELFI achieves 57-99% sensitivity at 98% specificity

> **Reference:** See [Citation & Scientific Background](../citation.md#fsr) for detailed paper summary.
