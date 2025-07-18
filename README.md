# cfDNAFE - Circulating Tumor DNA Feature Extractor

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

If cfdnafe is published to PyPI, you can use:

```bash
uv pip install cfdnafe
```

## Usage

### Motif-based Feature Extraction

Extract end motif, breakpoint motif, and Motif Diversity Score (MDS) from a BAM file:

```bash
cfdnafe motif path/to/input.bam -g path/to/reference.fa -o path/to/output_dir \
    --minlen 65 --maxlen 400 -k 3 --verbose
```

- Output files are written to EDM, BPM, and MDS subfolders in the output directory.
- Uses rich logging and progress bars for user-friendly feedback.

### Quality Control

```bash
cfdnafe quality-control path/to/bam_file.bam path/to/qc_metrics.txt
```

### Basic Feature Extraction (legacy)

```bash
cfdnafe extract-features path/to/bam_file.bam path/to/output_dir --reference-genome path/to/reference.fa
```

### Process Specific Region

```bash
cfdnafe process-region chr1:1000000-2000000 path/to/bam_file.bam path/to/region_features.txt
```

### Show Version

```bash
cfdnafe version
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

## License

MIT License
