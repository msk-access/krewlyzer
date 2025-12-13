# Fragment-level Methylation (UXM)

**Command**: `krewlyzer uxm`

## Purpose
Computes the proportions of Unmethylated (U), Mixed (X), and Methylated (M) fragments per region.

## Biological Context
Fragment-level methylation (UXM, Sun et al., Nature 2023) reveals cell-of-origin and cancer-specific methylation patterns in cfDNA.

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
- `--mark-input`, `-m`: Marker BED file (default: `data/MethMark/Atlas.U25.l4.hg19.bed`)
- `--map-quality`, `-q`: Minimum mapping quality (default: 30)
- `--min-cpg`, `-c`: Minimum CpG per fragment (default: 4)
- `--methy-threshold`, `-tM`: Methylation threshold (default: 0.75)
- `--unmethy-threshold`, `-tU`: Unmethylation threshold (default: 0.25)
- `--type`: Fragment type: SE or PE (default: SE)
- `--threads`, `-t`: Number of processes
