# Fragment-level Methylation (UXM)

**Command**: `krewlyzer uxm`

## Purpose
Computes the proportions of Unmethylated (U), Mixed (X), and Methylated (M) fragments per region.

## Biological Context
Fragment-level methylation (UXM, [Loyfer et al., 2022](../citation.md#uxm)) reveals cell-of-origin and cancer-specific methylation patterns in cfDNA.

## Usage
```bash
# Single-end (default)
krewlyzer uxm /path/to/bam_folder --output uxm_out [options]

# Paired-end mode
krewlyzer uxm /path/to/bam_folder --output uxm_out --type PE [options]
```

## Output
Output: `{sample}.UXM.tsv`

## Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--output` | `-o` | PATH | *required* | Output directory |
| `--sample-name` | `-s` | TEXT | | Override sample name |
| `--mark-input` | `-m` | PATH | | Path to genomic marker file |
| `--genome` | `-G` | TEXT | hg19 | Genome build (hg19/hg38) |
| `--verbose` | `-v` | FLAG | | Enable verbose logging |
| `--threads` | `-t` | INT | 0 | Number of threads (0=all) |

## Clinical Interpretation

### Fragment Classification
| Class | Threshold | Interpretation |
|-------|-----------|----------------|
| **U** (Unmethylated) | ≤25% methylated CpGs | Cell-type specific unmethylated regions |
| **X** (Mixed) | 25-75% | Heterogeneous/mosaic methylation |
| **M** (Methylated) | ≥75% methylated CpGs | Stably methylated regions |

### Healthy cfDNA Composition
Based on the Human Methylation Atlas:

| Cell Type | Contribution |
|-----------|--------------|
| Megakaryocytes | ~31% |
| Granulocytes | ~30% |
| Monocytes/Macrophages | ~20% |
| Endothelial | ~6% |
| Hepatocytes | ~3% |
| Lymphocytes | ~3% |

### Cancer Detection
| Pattern | Interpretation |
|---------|----------------|
| Altered tissue proportions | Tumor DNA shifts composition |
| Non-hematopoietic increase | Possible tumor-derived cfDNA |
| Resolution | Can detect ~0.1% tumor fractions |

> **Reference:** See [Citation & Scientific Background](../citation.md#uxm) for detailed paper summary.
