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
- `--arms-file`, `-a`: Chromosome arms BED (required)
- `--threads`, `-t`: Number of processes

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
