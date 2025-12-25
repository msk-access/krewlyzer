# Fragment Size Distribution (FSD)

**Command**: `krewlyzer fsd`

## Purpose
Computes high-resolution (5bp bins) fragment length distributions per chromosome arm.

## Biological Context
cfDNA fragmentation patterns at chromosome arms can reflect nucleosome positioning, chromatin accessibility, and cancer-specific fragmentation signatures. See [DELFI method](../citation.md#fsr) for details.

## Usage
```bash
krewlyzer fsd sample.bed.gz --arms-file krewlyzer/data/ChormosomeArms/hg19_arms.bed --output output_dir/ [options]
```
## Output
- `{sample}.FSD.tsv`: Frequency of fragment lengths per chromosome arm.

## Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--output` | `-o` | PATH | *required* | Output directory |
| `--sample-name` | `-s` | TEXT | | Override sample name |
| `--arms-file` | `-a` | PATH | | Chromosome arms BED file |
| `--pon-model` | `-P` | PATH | | PON model for z-score computation |
| `--genome` | `-G` | TEXT | hg19 | Genome build (hg19/hg38) |
| `--gc-correct` | | FLAG | True | Apply GC bias correction |
| `--verbose` | `-v` | FLAG | | Enable verbose logging |
| `--threads` | `-t` | INT | 0 | Number of threads (0=all) |

## Clinical Interpretation

### Healthy vs Cancer

| Metric | Healthy Plasma | Cancer (ctDNA present) |
|--------|----------------|------------------------|
| Modal peak | ~166bp | Left-shifted (~145bp) |
| 10bp periodicity | Clear nucleosome signal | May be altered |
| Arm-level variation | Minimal | Increased (correlates with CNAs) |

### What to Look For
- **Left-shift in distribution**: Tumor DNA is typically shorter
- **Arm-specific changes**: May correlate with copy number alterations
- **Loss of nucleosome periodicity**: Indicates altered chromatin structure

> **Reference:** See [Citation & Scientific Background](../citation.md#fsr) for detailed paper summary.
