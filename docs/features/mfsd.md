# Mutant Fragment Size Distribution (mFSD)

**Command**: `krewlyzer mfsd`

## Purpose
Compares the size distribution of mutant vs. wild-type fragments at variant sites, with support for all small variant types and 4-way fragment classification.

## Biological Context
Mutant ctDNA fragments are typically shorter (~145bp) than wild-type cfDNA (~166bp). This module quantifies this difference using high-depth targeted sequencing data, providing a sensitive marker for ctDNA presence. It performs a Kolmogorov-Smirnov (KS) test to statistically compare the distributions.

## Variant Types Supported

| Type | Example | Description |
|------|---------|-------------|
| **SNV** | A>T | Single nucleotide variant |
| **MNV** | AT>GC | Multi-nucleotide variant |
| **Insertion** | A>ATG | Pure insertion |
| **Deletion** | ATG>A | Pure deletion |
| **Complex** | ATG>CT | Mixed substitution + indel |

## Fragment Classification

Fragments are classified into 4 categories based on the allele they support:

| Category | Description | Use Case |
|----------|-------------|----------|
| **REF** | Supports reference allele | Healthy cfDNA baseline |
| **ALT** | Supports alternate allele | Tumor signal |
| **NonREF** | Non-REF, non-ALT, non-N | Sequencing errors, subclones |
| **N** | Contains N at variant position | Low quality, uninformative |

## Usage
```bash
krewlyzer mfsd sample.bam --input-file variants.vcf --output output_dir/ [options]
```

## Options
- `--input-file`, `-i`: VCF or MAF file (required)
- `--mapq`, `-q`: Minimum mapping quality (default: 20)
- `--output-distributions`, `-d`: Output per-variant size distributions
- `--sample-name`, `-s`: Sample name for output file
- `--threads`, `-t`: Number of threads (0=all cores)

## Output Format

### Main Output: `{sample}.mFSD.tsv` (39 columns)

#### Variant Info (5)
| Column | Description |
|--------|-------------|
| `Chrom` | Chromosome |
| `Pos` | 1-based position |
| `Ref` | Reference allele |
| `Alt` | Alternate allele |
| `VarType` | SNV, MNV, INS, DEL, or COMPLEX |

#### Counts (5)
| Column | Description |
|--------|-------------|
| `REF_Count` | Fragments supporting REF |
| `ALT_Count` | Fragments supporting ALT |
| `NonREF_Count` | Fragments with errors |
| `N_Count` | Fragments with N bases |
| `Total_Count` | Total fragments |

#### Mean Sizes (4)
| Column | Description |
|--------|-------------|
| `REF_MeanSize` | Mean size of REF fragments |
| `ALT_MeanSize` | Mean size of ALT fragments |
| `NonREF_MeanSize` | Mean size of error fragments |
| `N_MeanSize` | Mean size of N fragments |

#### Pairwise Comparisons (18)
For each pair (6 combinations): Delta, KS statistic, KS p-value

- `Delta_ALT_REF`, `KS_ALT_REF`, `KS_Pval_ALT_REF` (Primary)
- `Delta_ALT_NonREF`, `KS_ALT_NonREF`, `KS_Pval_ALT_NonREF`
- `Delta_REF_NonREF`, `KS_REF_NonREF`, `KS_Pval_REF_NonREF`
- `Delta_ALT_N`, `KS_ALT_N`, `KS_Pval_ALT_N`
- `Delta_REF_N`, `KS_REF_N`, `KS_Pval_REF_N`
- `Delta_NonREF_N`, `KS_NonREF_N`, `KS_Pval_NonREF_N`

#### Derived Metrics (5)
| Column | Description |
|--------|-------------|
| `VAF_Proxy` | ALT / (REF + ALT) |
| `Error_Rate` | NonREF / Total |
| `N_Rate` | N / Total |
| `Size_Ratio` | ALT_MeanSize / REF_MeanSize |
| `Quality_Score` | 1 - N_Rate - Error_Rate |

#### Quality Flags (2)
| Column | Description |
|--------|-------------|
| `ALT_Confidence` | HIGH (≥5), LOW (1-4), NONE (0) |
| `KS_Valid` | TRUE if both REF and ALT ≥2 fragments |

### Optional: `{sample}.mFSD.distributions.tsv`

With `--output-distributions`, outputs per-variant size distributions:

```tsv
Chrom  Pos  Ref  Alt  Category  Size  Count
chr1   12345  A   T   REF       145   3
chr1   12345  A   T   REF       166   12
chr1   12345  A   T   ALT       142   2
```

## Clinical Interpretation

### Healthy vs Cancer

| Metric | Healthy | Cancer (ctDNA) |
|--------|---------|----------------|
| `Delta_ALT_REF` | ~0 | Negative (ALT shorter) |
| `Size_Ratio` | ~1.0 | <1.0 |
| `VAF_Proxy` | 0 | >0 (correlates with tumor fraction) |

### MRD Settings
- Low fragment counts (1-2) produce `NA` for statistics
- `ALT_Confidence` indicates reliability
- `KS_Valid` flags when statistical test is meaningful
