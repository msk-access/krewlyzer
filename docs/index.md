# Welcome to Krewlyzer

<p align="center">
  <img src="https://raw.githubusercontent.com/msk-access/krewlyzer/main/src/krewlyzer/logo.svg" alt="Krewlyzer logo" width="120"/>
</p>

**Krewlyzer** is a robust, user-friendly command-line toolkit for extracting a wide range of biological features from cell-free DNA (cfDNA) sequencing data. It is designed for cancer genomics, liquid biopsy research, and clinical bioinformatics, providing high-performance, reproducible feature extraction from BAM files.

Krewlyzer draws inspiration from [cfDNAFE](https://github.com/Cuiwanxin1998/cfDNAFE) and implements state-of-the-art methods for fragmentation, motif, and methylation analysis, all in a modern Pythonic interface with rich parallelization and logging.

---

## TL;DR

Krewlyzer extracts **fragmentomics features** from cfDNA sequencing data:

1. **Input**: BAM file (aligned reads from blood sample)
2. **Output**: Tables of numerical features for ML/statistical analysis
3. **Use case**: Cancer detection, treatment monitoring, tissue-of-origin

**One command does it all:**
```bash
krewlyzer run-all -i sample.bam --reference hg19.fa --output results/
```

> **New to cfDNA?** Start with [What is Cell-Free DNA?](getting-started/concepts.md)  
> **Visual learner?** See the [Overview PDF](resources/overview.md)  
> **Need terminology help?** See the [Glossary](reference/glossary.md)

---

## Key Features

*   **Motif Analysis**: End motifs, breakpoint motifs, and diversity scores.
*   **Fragment Size Analysis**: Coverage (FSC), Ratios (FSR), and Distributions (FSD).
*   **Nucleosome Positioning**: Windowed Protection Scores (WPS).
*   **Tissue of Origin**: Orientation-aware Fragmentation (OCF).
*   **Methylation**: Fragment-level methylation patterns (UXM).
*   **Mutant Analysis**: Mutant vs. Wild-type fragment size comparison (mFSD).

## System Requirements

- Linux or macOS (tested on Ubuntu 20.04, macOS 12+)
- Python 3.8+
- ≥16GB RAM recommended for large BAM files
- [Docker](https://www.docker.com/) (optional, for easiest setup)

## Installation

### With Docker (Recommended)
```bash
docker pull ghcr.io/msk-access/krewlyzer:latest
# Example usage:
docker run --rm -v $PWD:/data ghcr.io/msk-access/krewlyzer:latest run-all -i /data/sample.bam --reference /data/hg19.fa --output /data/output_dir
```

### With uv / pip
```bash
uv venv .venv
source .venv/bin/activate
uv pip install krewlyzer
```
Or install from source:
```bash
git clone https://github.com/msk-access/krewlyzer.git
cd krewlyzer
uv pip install .
```
