# Krewlyzer - Circulating Tumor DNA Feature Extractor

<p align="center">
  <img src="krewlyzer/logo.svg" alt="Krewlyzer logo" width="120"/>
</p>

<p align="center">
  <a href="https://pypi.org/project/krewlyzer/"><img src="https://img.shields.io/pypi/v/krewlyzer.svg?color=blue" alt="PyPI version"></a>
  <a href="https://github.com/msk-access/krewlyzer/actions"><img src="https://github.com/msk-access/krewlyzer/workflows/CI/badge.svg" alt="GitHub Actions"></a>
  <a href="https://github.com/msk-access/krewlyzer/pkgs/container/krewlyzer"><img src="https://img.shields.io/badge/docker-ready-blue.svg" alt="Docker"></a>
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

Calculate fragment size coverage (FSC) features for all `.bed.gz` files in a folder. The input folder should be the output directory produced by `motif`, containing the `.bed.gz` files. Output files are written to the output directory, one per `.bed.gz` file.

**Parallel Processing:**
- FSC now supports process-based parallel execution (multiprocessing) for per-file processing. Use the `--threads`/`-t` option to set the number of parallel processes (cores) to use. This can significantly improve performance on multi-core systems.

**Usage:**

```bash
krewlyzer fsc BEDGZ_PATH --output OUTDIR [options]
```

**Options:**

| Option | Description | Default |
| --- | --- | --- |
| `--bin-input`, `-b` | Path to bin file | `data/ChormosomeBins/hg19_window_100kb.bed` |
| `--windows`, `-w` | Window size | 100000 |
| `--continue-n`, `-c` | Consecutive window number | 50 |
| `--threads`, `-t` | Number of processes | 1 |

### Fragment Size Ratio (FSR) Feature Calculation

Calculate fragment size ratio (FSR) features for all `.bed.gz` files in a folder. The input folder should be the output directory produced by `motif`, containing the `.bed.gz` files. Output files are written to the output directory, one per `.bed.gz` file.

**Parallel Processing:**
- FSR now supports process-based parallel execution (multiprocessing) for per-file processing. Use the `--threads`/`-t` option to set the number of parallel processes (cores) to use. This can significantly improve performance on multi-core systems.

**Usage:**

```bash
krewlyzer fsr BEDGZ_PATH --output OUTDIR [options]
```

**Options:**

| Option | Description | Default |
| --- | --- | --- |
| `--bin-input`, `-b` | Path to bin file | `data/ChormosomeBins/hg19_window_100kb.bed` |
| `--windows`, `-w` | Window size | 100000 |
| `--continue-n`, `-c` | Consecutive window number | 50 |
| `--threads`, `-t` | Number of processes | 1 |

### Fragment Size Distribution (FSD) Feature Calculation

Calculate fragment size distribution (FSD) features for all `.bed.gz` files in a folder. The input folder should be the output directory produced by `motif`, containing the `.bed.gz` files. Output files are written to the output directory, one per `.bed.gz` file.

**Parallel Processing:**
- FSD now supports process-based parallel execution (multiprocessing) for per-file processing. Use the `--threads`/`-t` option to set the number of parallel processes (cores) to use. This can significantly improve performance on multi-core systems.

**Usage:**

```bash
krewlyzer fsd BEDGZ_PATH --arms-file ARMS_FILE --output OUTDIR [options]
```

**Options:**

| Option | Description | Default |
| --- | --- | --- |
| `--arms-file`, `-a` | Path to arms/region file (BED format) |  |
| `--threads`, `-t` | Number of processes | 1 |

### Windowed Protection Score (WPS) Feature Calculation

Calculate the Windowed Protection Score (WPS) for each region in a transcript/region file.

```bash
krewlyzer wps <motif_output_dir> --output <wps_output_dir>
```

- Input: Directory of `.bed.gz` files from `motif` command
- Output: `.WPS.tsv.gz` per region/sample
- Options:
  - `--tsv-input` (optional, default: `data/TranscriptAnno/transcriptAnno-hg19-1kb.tsv`)
  - `--wpstype` (`L` for long [default], `S` for short)
  - `--empty` (keep files of empty regions)
  - `--threads` (number of threads)

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

### Output Structure Example

```
output_dir/
├── EDM/
│   ├── sample1.edm
│   ├── sample2.edm
│   └── ...
├── BPM/
│   ├── sample1.bpm
│   ├── sample2.bpm
│   └── ...
├── MDS/
│   ├── sample1.mds
│   ├── sample2.mds
│   └── ...
├── FSC/
│   ├── sample1.fsc
│   ├── sample2.fsc
│   └── ...
├── FSR/
│   ├── sample1.fsr
│   ├── sample2.fsr
│   └── ...
├── FSD/
│   ├── sample1.fsd
│   ├── sample2.fsd
│   └── ...
└── WPS/
    ├── region1.wps.tsv.gz
    ├── region2.wps.tsv.gz
    └── ...

```bash
krewlyzer fsr <motif_output_dir> --output <fsr_output_dir> [options]
```

**Options:**
- `--bin-input`, `-b`: Path to bin file (default: `data/ChormosomeBins/hg19_window_100kb.bed`)
- `--windows`, `-w`: Window size (default: 100000)
- `--continue-n`, `-c`: Consecutive window number (default: 50)
- `--threads`, `-t`: Number of processes (default: 1; increase for faster parallel processing)
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

### WPS: Windowed Protection Score
Calculates the Windowed Protection Score (WPS) for each region in a transcript/region file.

```bash
krewlyzer wps <motif_output_dir> --output <wps_output_dir>
```
- Input: Directory of `.bed.gz` files from `motif` command
- Output: `.WPS.tsv.gz` per region/sample
- Options:
  - `--tsv-input` (optional, default: `data/TranscriptAnno/transcriptAnno-hg19-1kb.tsv`)
  - `--wpstype` (`L` for long [default], `S` for short)
  - `--empty` (keep files of empty regions)
  - `--threads` (number of threads)

If `--tsv-input` is not specified, the default transcript annotation file (`data/TranscriptAnno/transcriptAnno-hg19-1kb.tsv`) will be used.

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
