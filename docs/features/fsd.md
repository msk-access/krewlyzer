# Fragment Size Distribution (FSD)

**Command**: `krewlyzer fsd`

## Purpose
Computes high-resolution (5bp bins) fragment length distributions per chromosome arm.

## Biological Context
cfDNA fragmentation patterns at chromosome arms can reflect nucleosome positioning, chromatin accessibility, and cancer-specific fragmentation signatures.

## Usage
```bash
krewlyzer fsd sample.bed.gz --arms-file krewlyzer/data/ChormosomeArms/hg19_arms.bed --output fsd_out.txt [options]
```
## Output
- `{sample}.FSD.txt`: Frequency of fragment lengths per chromosome arm.

## Options
- `--arms-file`, `-a`: Chromosome arms BED (required)
- `--threads`, `-t`: Number of processes
