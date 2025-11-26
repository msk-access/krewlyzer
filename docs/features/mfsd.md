# Mutant Fragment Size Distribution (mFSD)

**Command**: `krewlyzer mfsd`

## Purpose
Compares the size distribution of mutant vs. wild-type reads at variant sites.

## Biological Context
Mutant ctDNA fragments are typically shorter than wild-type cfDNA. This module quantifies this difference using high-depth targeted sequencing data, providing a sensitive marker for ctDNA presence. It performs a Kolmogorov-Smirnov (KS) test to statistically compare the distributions.

## Usage
```bash
krewlyzer mfsd sample.bam --input variants.vcf --output output.tsv [options]
```

## Options
- `--input`, `-i`: VCF or MAF file (required).
- `--format`, `-f`: Input format ('auto', 'vcf', 'maf'). Default: 'auto'.
- `--map-quality`, `-q`: Minimum mapping quality (default: 20).
