<p align="center">
  <img src="krewlyzer/logo.svg" alt="Krewlyzer logo" width="120"/>
</p>

# Krewlyzer - Circulating Tumor DNA Feature Extractor

<p align="center">
  <a href="https://pypi.org/project/krewlyzer/"><img src="https://img.shields.io/pypi/v/krewlyzer.svg?color=blue" alt="PyPI version"></a>
  <a href="https://github.com/msk-acess/krewlyzer/actions"><img src="https://github.com//msk-acess/krewlyzer/workflows/CI/badge.svg" alt="GitHub Actions"></a>
  <a href="https://github.com//msk-acess/krewlyzer/pkgs/container/krewlyzer"><img src="https://img.shields.io/badge/docker-ready-blue.svg" alt="Docker"></a>
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

### Fragment Size Coverage (FSC) Feature Calculation

Calculate fragment size coverage (FSC) features from motif output:

```bash
krewlyzer fsc path/to/motif_output_dir -o path/to/fsc_output_dir
```

- `path/to/motif_output_dir` should be the directory containing `.bed.gz` files output by the `motif` command.
- FSC features will be written to the specified output directory, with one result file per input `.bed.gz` file.
- By default, the command uses the provided 100kb bin file (`data/ChormosomeBins/hg19_window_100kb.bed`); you can override this with `--bin-input`.
- Additional options: `--windows` (window size), `--continue-n` (number of consecutive bins), `--threads` (parallelism).

Each FSC output file contains region-based z-scores for short, intermediate, long, and total fragment sizes, GC-corrected.


### Show Version

```bash
krewlyzer version
```

## Input Data

- For motif-based feature extraction: BAM files aligned to GRCh37 reference genome and a reference genome FASTA file.
- For FSC feature extraction: Directory containing `.bed.gz` files produced by the `motif` command (motif output directory).

## Output

The tool generates:
- Motif feature files: End Motif (EDM), Breakpoint Motif (BPM), and Motif Diversity Score (MDS) per sample (written to subfolders `EDM/`, `BPM/`, and `MDS/` under your motif output directory)
- Fragment Size Coverage (FSC) feature files: one per input `.bed.gz` file, written to your specified FSC output directory
- Fragment Size Ratio (FSR) feature files: one per input `.bed.gz` file, written to your specified FSR output directory
- Fragment Size Distribution (FSD) feature files: one per input `.bed.gz` file, written to your specified FSD output directory

---

## FSC, FSR, and FSD CLI Usage

### FSC: Fragment Size Coverage
Calculates z-score normalized fragment size coverage for each region/bin.

```bash
krewlyzer fsc <motif_output_dir> --output <fsc_output_dir>
```
- Input: Directory of `.bed.gz` files from `motif` command
- Output: `.FSC.txt` per sample
- Options: `--bin-input`, `--windows`, `--continue-n`

### FSR: Fragment Size Ratio
Calculates the ratio of short, intermediate, and long fragments per region/bin.

```bash
krewlyzer fsr <motif_output_dir> --output <fsr_output_dir>
```
- Input: Directory of `.bed.gz` files from `motif` command
- Output: `.FSR.txt` per sample
- Options: `--bin-input`, `--windows`, `--continue-n`

### FSD: Fragment Size Distribution
Calculates the distribution of fragment sizes (in 5bp bins from 65-399bp) per region.

```bash
krewlyzer fsd <motif_output_dir> --arms-file <arms.bed> --output <fsd_output_dir>
```
- Input: Directory of `.bed.gz` files from `motif` command
- Output: `.FSD.txt` per sample
- Options: `--arms-file` (required)

See `krewlyzer <command> --help` for all options and details.

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

This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0). See the [LICENSE](./LICENSE) file for full terms.
