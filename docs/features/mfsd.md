# Mutant Fragment Size Distribution (mFSD)

**Command**: `krewlyzer mfsd`

## Purpose
Compares the size distribution of mutant vs. wild-type reads at variant sites.

## Biological Context
Mutant ctDNA fragments are typically shorter than wild-type cfDNA. This module quantifies this difference using high-depth targeted sequencing data, providing a sensitive marker for ctDNA presence. It performs a Kolmogorov-Smirnov (KS) test to statistically compare the distributions.

## Usage
```bash
krewlyzer mfsd sample.bam --input variants.vcf --output output_dir/ [options]
```

## Options
- `--input`, `-i`: VCF or MAF file (required).
- `--format`, `-f`: Input format ('auto', 'vcf', 'maf'). Default: 'auto'.
- `--map-quality`, `-q`: Minimum mapping quality (default: 20).

## Output Format

Output file: `{sample}.mFSD.tsv`

| Column | Description |
|--------|-------------|
| `sample` | Sample ID |
| `mutation_count` | Number of reads supporting the variant (Mutant) |
| `wt_count` | Number of reads supporting reference (Wild-type) |
| `mutation_mean_size` | Mean fragment length of mutant reads |
| `wt_mean_size` | Mean fragment length of wild-type reads |
| `size_diff` | `wt_mean_size - mutation_mean_size` (Positive = Mutant is shorter) |
| `ks_stat` | Kolmogorov-Smirnov test statistic |
| `p_value` | P-value of the KS test |

## Calculation Details

1.  **Variant Locus**: For each variant in the input VCF/MAF.
2.  **Read Collection**: Reads covering the variant site are classified as **Mutant** (carry alt allele) or **Wild-type** (carry ref allele).
3.  **Size Distribution**: The fragment lengths of both populations are collected.
4.  **Statistical Test**: A two-sample Kolmogorov-Smirnov test compares the cumulative distribution functions (CDFs) of lengths.
    - **Hypothesis**: The size distributions differ.
    - **Expected**: Mutant fragments are shorter, shifting the CDF left.

