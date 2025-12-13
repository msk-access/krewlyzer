# Fragment Size Ratio (FSR)

**Command**: `krewlyzer fsr`

## Purpose
Computes ratios of fragment size classes per genomic window. Unlike FSC, this metric is self-normalizing (ratio) and focuses on the **proportion** of fragment sizes, which is a powerful biomarker for tumor fraction estimation.

## Biological Context
The ratio of short to long fragments is a key indicator of tumor burden in cfDNA ("fragmentomics").

- **Short Fragments (65-149bp)**: Tumor DNA is typically shorter (~145bp) than healthy DNA (~166bp).
- **Ultra-short Fragments (65-100bp)**: Associated with transcription factor binding sites and open chromatin.
- **Long Fragments (261-399bp)**: Often di-nucleosomes, representing stable, healthy chromatin (Leukocytes).

**Key Biomarkers:**
- **Short/Long Ratio**: The primary metric. Higher ratio = higher probability of tumor DNA.
- **Ultra-short Ratio**: Indicates active gene regulation/transcription factor binding.
- **Nucleosome Footprints**: The 10bp precidion of the ranges (150bp vs 167bp) helps separate mono-nucleosomes.

## Usage
```bash
krewlyzer fsr sample.bed.gz -o output_dir/ --sample-name SAMPLE [options]
```

## Options
- `--bin-input, -b`: Bin file (default: `data/ChormosomeBins/hg19_window_100kb.bed`)
- `--windows, -w`: Window size (default: 100000)
- `--continue-n, -c`: Consecutive window number (default: 50)
- `--threads, -p`: Number of threads.

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

## Interpretation Guide

### Ratios
- **High `short_long_ratio`**: Indicates higher tumor fraction (ctDNA) or open chromatin.
- **High `ultra_short_ratio`**: Indicates regions of high transcription factor activity.
- **High `long_ratio`**: Indicates stable, nucleosomal DNA (usually healthy background).

### Counts
Provided for reference and custom normalization. Note that counts depend on sequencing depth, whereas ratios are generally depth-independent.
