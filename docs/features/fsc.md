# Fragment Size Coverage (FSC)

**Command**: `krewlyzer fsc`

## Purpose
Computes z-scored coverage of cfDNA fragments in different size ranges per genomic bin, with GC correction. This helps identify copy number variations (CNVs) and coverage anomalies specific to certain fragment sizes.

## Biological Context
cfDNA fragment size profiles are informative for cancer detection and tissue-of-origin. FSC measures the normalized coverage depth of:

- **Short (65-149bp)**: Enriched for tumor-derived cfDNA (ctDNA) in cancer patients.
- **Intermediate (151-259bp)**: Represents mono-nucleosomal fragments.
- **Long (261-399bp)**: Represents di-nucleosomal fragments, often from healthy cells.
- **Total (65-399bp)**: Overall coverage.

Differences in coverage patterns between size classes can reveal:
- **Copy Number Alterations (CNAs)**: Detected via Total and size-specific coverage.
- **Chromatin Structure**: Open chromatin (more short fragments) vs. closed chromatin (more long fragments).

## Usage
```bash
krewlyzer fsc sample.bed.gz -o output_dir/ --sample-name SAMPLE [options]
```

## Options
- `--bin-input, -b`: Bin file (default: `data/ChormosomeBins/hg19_window_100kb.bed`)
- `--windows, -w`: Window size (default: 100000)
- `--continue-n, -c`: Consecutive window number (default: 50) - aggregates adjacent bins.
- `--threads, -p`: Number of threads.

## Output Format

Output: `{sample}.FSC.tsv`

| Column | Description |
|--------|-------------|
| `region` | Genomic region (chr:start-end) |
| `short-fragment-zscore` | Z-score of short fragment coverage |
| `itermediate-fragment-zscore` | Z-score of intermediate fragment coverage |
| `long-fragment-zscore` | Z-score of long fragment coverage |
| `total-fragment-zscore` | Z-score of total fragment coverage |

## Interpretation Guide

| Metric | High Z-Score (>2) | Low Z-Score (<-2) |
|--------|-------------------|-------------------|
| **Total** | Copy number gain / Accesssible region | Copy number loss / Closed chromatin |
| **Short** | **Tumor-enriched** / Open chromatin | Depleted ctDNA / Closed chromatin |
| **Long** | Healthy/Leukocyte DNA enriched | Fragmentation / Open chromatin |

**Note**: Counts are GC-corrected before Z-score calculation to remove sequencing bias.

## Calculation Details

1.  **Binning**: Fragments are counted in 100kb genomic bins (or custom regions).
2.  **GC Correction**: Raw counts are adjusted for GC content bias using Loess regression (GC vs Count).
3.  **Aggregation**: Adjusted counts are summed over `N=50` consecutive bins (5Mb effective window) to smooth noise.
4.  **Z-Score Normalization**:
    $$ Z = \frac{X - \mu}{\sigma} $$
    Where $X$ is the summed count for the window, and $\mu, \sigma$ are the genome-wide mean and standard deviation.

