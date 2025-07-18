<p align="center">
  <img src="krewlyzer/logo.svg" alt="Krewlyzer logo" width="120"/>
</p>

# Krewlyzer - Circulating Tumor DNA Feature Extractor

<p align="center">
  <a href="https://pypi.org/project/krewlyzer/"><img src="https://img.shields.io/pypi/v/krewlyzer.svg?color=blue" alt="PyPI version"></a>
  <a href="https://github.com/YOUR_GITHUB_USER/YOUR_REPO_NAME/actions"><img src="https://github.com/YOUR_GITHUB_USER/YOUR_REPO_NAME/workflows/CI/badge.svg" alt="GitHub Actions"></a>
  <a href="https://github.com/YOUR_GITHUB_USER/YOUR_REPO_NAME/pkgs/container/krewlyzer"><img src="https://img.shields.io/badge/docker-ready-blue.svg" alt="Docker"></a>
</p>

A command-line tool for extracting features from circulating tumor DNA (ctDNA) BAM files aligned to GRCh37.

## Installation

This project uses [uv](https://github.com/astral-sh/uv), a fast Python package manager. If you don't have uv, install it from their documentation.

### Using a uv Virtual Environment (Recommended)

1. **Create and activate a uv virtual environment:**

```bash
uv venv .venv
source .venv/bin/activate
```

2. **Install the package and dependencies locally:**

```bash
uv pip install .
```

If krewlyzer is published to PyPI, you can use:

```bash
uv pip install krewlyzer
```

## Usage

### Motif-based Feature Extraction

Extract end motif, breakpoint motif, and Motif Diversity Score (MDS) from a BAM file:

```bash
krewlyzer motif path/to/input.bam -g path/to/reference.fa -o path/to/output_dir \
    --minlen 65 --maxlen 400 -k 3 --verbose
```

- Output files are written to EDM, BPM, and MDS subfolders in the output directory.
- Uses rich logging and progress bars for user-friendly feedback.

### Quality Control

```bash
krewlyzer quality-control path/to/bam_file.bam path/to/qc_metrics.txt
```

### Basic Feature Extraction (legacy)

```bash
krewlyzer extract-features path/to/bam_file.bam path/to/output_dir --reference-genome path/to/reference.fa
```

### Process Specific Region

```bash
krewlyzer process-region chr1:1000000-2000000 path/to/bam_file.bam path/to/region_features.txt
```

### Show Version

```bash
krewlyzer version
```

## Input Data

The tool expects:
- BAM files aligned to GRCh37 reference genome
- Reference genome FASTA file
- Genomic regions in format chr:start-end

## Output

The tool generates:
- Motif feature files: End Motif (EDM), Breakpoint Motif (BPM), and Motif Diversity Score (MDS) per sample
- Quality control metrics
- Region-specific feature files

Motif outputs are written to subfolders `EDM/`, `BPM/`, and `MDS/` under your output directory.

## Requirements

- Python 3.8+
- uv
- pysam
- pandas
- biopython
- rich (for logging)
- pytest (for testing)
- typer

## Testing

Run CLI unit tests using pytest:

```bash
pytest tests/
```

This will test CLI argument parsing and help output. For full integration tests, provide small test BAM and reference files in `tests/data/`.


## Reference

This project was inspired by and references the [cfDNAFE](https://github.com/Cuiwanxin1998/cfDNAFE) repository. Please visit their repo for more information and tools related to cfDNA feature extraction.

## License

MIT License
