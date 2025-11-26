# Motif-based Feature Extraction

**Command**: `krewlyzer motif`

## Purpose
Extracts end motif, breakpoint motif, and Motif Diversity Score (MDS) from sequencing fragments.

## Biological Context
Motif analysis of cfDNA fragment ends can reveal tissue-of-origin, nucleosome positioning, and mutational processes. MDS quantifies motif diversity, which may be altered in cancer.

## Usage
```bash
krewlyzer motif path/to/input.bam -g path/to/reference.fa -o path/to/output_dir \
    --minlen 65 --maxlen 400 -k 3 --verbose
```

## Output
- `EDM/`: End Motif frequencies
- `BPM/`: Breakpoint Motif frequencies
- `MDS/`: Motif Diversity Scores
- `*.bed.gz`: Intermediate file used for other modules (FSC, FSR, FSD, WPS, OCF).
