# Krewlyzer: Comprehensive cfDNA Feature Extraction Toolkit

<p align="center">
  <img src="src/krewlyzer/logo.svg" alt="Krewlyzer logo" width="120"/>
</p>

<p align="center">
  <a href="https://pypi.org/project/krewlyzer/"><img src="https://img.shields.io/pypi/v/krewlyzer.svg?color=blue" alt="PyPI version"></a>
  <a href="https://github.com/msk-access/krewlyzer/actions"><img src="https://github.com/msk-access/krewlyzer/workflows/Release/badge.svg" alt="GitHub Actions"></a>
  <a href="https://github.com/msk-access/krewlyzer/pkgs/container/krewlyzer"><img src="https://img.shields.io/badge/docker-ready-blue.svg" alt="Docker"></a>
</p>

**Krewlyzer** is a robust, user-friendly command-line toolkit for extracting a wide range of biological features from cell-free DNA (cfDNA) sequencing data. It is designed for cancer genomics, liquid biopsy research, and clinical bioinformatics, providing high-performance, reproducible feature extraction from BAM files. Krewlyzer draws inspiration from [cfDNAFE](https://github.com/Cuiwanxin1998/cfDNAFE) and implements state-of-the-art methods for fragmentation, motif, and methylation analysis, all in a modern Pythonic interface with rich parallelization and logging.

> [!TIP]
> **Full Documentation**: For detailed usage, feature descriptions, and pipeline integration, visit our [Documentation Site](https://msk-access.github.io/krewlyzer/).

---

## Table of Contents
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Reference Data](#reference-data)
- [Command Summary](#command-summary)
- [Typical Workflow](#typical-workflow)
- [Feature Details & Usage](#feature-details--usage)
  - [Fragment Extraction](#fragment-extraction)
  - [Motif-based Feature Extraction](#motif-based-feature-extraction)
  - [Fragment Size Coverage (FSC)](#fragment-size-coverage-fsc)
  - [Fragment Size Ratio (FSR)](#fragment-size-ratio-fsr)
  - [Fragment Size Distribution (FSD)](#fragment-size-distribution-fsd)
  - [Windowed Protection Score (WPS)](#windowed-protection-score-wps)
  - [Orientation-aware Fragmentation (OCF)](#orientation-aware-fragmentation-ocf)
  - [Fragment-level Methylation (UXM)](#fragment-level-methylation-uxm)
  - [Mutant Fragment Size Distribution (mFSD)](#mutant-fragment-size-distribution-mfsd)
  - [Run All Features](#run-all-features)
- [Output Structure Examples](#output-structure-examples)
- [Troubleshooting](#troubleshooting)
- [Citation & Acknowledgements](#citation--acknowledgements)

---

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
docker run --rm -v $PWD:/data ghcr.io/msk-access/krewlyzer:latest motif /data/sample.bam -g /data/hg19.fa -o /data/motif_out
```

### With uv Virtual Environment
```bash
uv venv .venv
source .venv/bin/activate
uv pip install .
```
Or install from PyPI:
```bash
uv pip install krewlyzer
```

## Reference Data
- **Reference Genome (FASTA):**
  - Download GRCh37/hg19 from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/)
  - BAMs must be sorted, indexed, and aligned to the same build
- **Bin/Region/Marker Files:**
  - Provided in `src/krewlyzer/data/` (see options for each feature)

## Command Summary

| Command   | Description                                  |
|-----------|----------------------------------------------|
| extract   | Extract fragments from BAM to BED            |
| motif     | Motif-based feature extraction               |
| fsc       | Fragment size coverage                       |
| fsr       | Fragment size ratio                          |
| fsd       | Fragment size distribution                   |
| wps       | Windowed protection score                    |
| ocf       | Orientation-aware fragmentation              |
| uxm       | Fragment-level methylation (SE/PE)           |
| mfsd      | Mutant fragment size distribution (4-way)    |
| run-all   | Run all features for a BAM                   |

## Typical Workflow
The recommended way to run krewlyzer is using the **Unified Pipeline** via `run-all`, which processes the BAM file in a single pass for maximum efficiency.

```bash
# Optimized Unified Pipeline
krewlyzer run-all sample.bam --reference hg19.fa --output output_dir \
    --variants variants.maf --bin-input targets.bed --threads 4
```

Alternatively, you can run tools individually. Note that most tools require a fragment BED file (`.bed.gz`) produced by the `extract` command.

```bash
# 1. Extract fragments (BAM -> BED.gz)
krewlyzer extract sample.bam -g hg19.fa -o output_dir

# 2. Run feature tools using the BED file
krewlyzer fsc output_dir/sample.bed.gz --output fsc_out.txt
# ... (wps, fsd, ocf, etc.)

# 3. Motif analysis (Independent of BED, uses BAM directly)
krewlyzer motif sample.bam -g hg19.fa -o output_dir 
```

---

## Feature Details & Usage

### Fragment Extraction
**Purpose:** Extracts cfDNA fragments from BAM to BED format with configurable filters.

**Usage:**
```bash
krewlyzer extract sample.bam -g reference.fa -o output_dir/ [options]
```
- Output: `{sample}.bed.gz` (tabix indexed)
- Options:
  - `--mapq`, `-q`: Minimum mapping quality (default: 20)
  - `--minlen`: Minimum fragment length (default: 65)
  - `--maxlen`: Maximum fragment length (default: 400)
  - `--exclude-regions`, `-x`: BED file of regions to exclude

### Motif-based Feature Extraction
**Purpose:** Extracts end motif, breakpoint motif, and Motif Diversity Score (MDS) from sequencing fragments.

**Biological context:** Motif analysis of cfDNA fragment ends can reveal tissue-of-origin, nucleosome positioning, and mutational processes. MDS quantifies motif diversity, which may be altered in cancer.

**Usage:**
```bash
krewlyzer motif path/to/input.bam -g path/to/reference.fa -o path/to/output_dir \
    --minlen 65 --maxlen 400 -k 4 --verbose
```
- Output: `{sample}.EndMotif.tsv`, `{sample}.BreakPointMotif.tsv`, `{sample}.MDS.tsv`
- Rich logging and progress bars for user-friendly feedback.

### Fragment Size Coverage (FSC)
**Purpose:** Computes z-scored coverage of cfDNA fragments in different size ranges, per genomic bin (default: 100kb), with GC correction.

**Biological context:** cfDNA fragment size profiles are informative for cancer detection and tissue-of-origin. FSC quantifies the coverage of short (65-150bp), intermediate (151-260bp), long (261-400bp), and total (65-400bp) fragments, normalized to genome-wide means.

**Usage:**
```bash
krewlyzer fsc motif_out --output fsc_out [options]
```
- Input: `.bed.gz` files from `motif` command
- Output: One `.FSC` file per sample
- Options:
  - `--bin-input`, `-b`: Bin file (default: `data/ChormosomeBins/hg19_window_100kb.bed`)
  - `--windows`, `-w`: Window size (default: 100000)
  - `--continue-n`, `-c`: Super-bin size (default: 50)
  - `--threads`, `-t`: Number of processes

### Fragment Size Ratio (FSR)
**Purpose:** Calculates the ratio of ultra-short/short/intermediate/long fragments per bin.

**Biological context:** The DELFI method (Mouliere et al., 2018) showed that cfDNA fragment size ratios are highly informative for cancer detection. Krewlyzer uses ultra-short (65-100bp), short (65-150bp), intermediate (151-260bp), and long (261-400bp) bins. The ultra-short bin is a highly specific marker for ctDNA.

**Usage:**
```bash
krewlyzer fsr motif_out --output fsr_out [options]
```
- Input: `.bed.gz` files from `motif` command
- Output: One `.FSR` file per sample
- Options: Same as FSC

### Fragment Size Distribution (FSD)
**Purpose:** Computes high-resolution (5bp bins) fragment length distributions per chromosome arm.

**Biological context:** cfDNA fragmentation patterns at chromosome arms can reflect nucleosome positioning, chromatin accessibility, and cancer-specific fragmentation signatures.

**Usage:**
```bash
krewlyzer fsd motif_out --arms-file krewlyzer/data/ChormosomeArms/hg19_arms.bed --output fsd_out [options]
```
- Input: `.bed.gz` files from `motif` command
- Output: One `.FSD` file per sample
- Options:
  - `--arms-file`, `-a`: Chromosome arms BED (required)
  - `--threads`, `-t`: Number of processes

### Windowed Protection Score (WPS)
**Purpose:** Computes nucleosome protection scores (WPS) for each region in a transcript/region file.

**Biological context:** The WPS (Snyder et al., 2016) quantifies nucleosome occupancy and chromatin accessibility by comparing fragments spanning a window to those ending within it. High WPS indicates nucleosome protection; low WPS, open chromatin.

**Usage:**
```bash
krewlyzer wps motif_out --output wps_out [options]
```
- Input: `.bed.gz` files from `motif` command
- Output: `.WPS.tsv.gz` per region/sample
- Options:
  - `--tsv-input`: Transcript region file (default: `data/TranscriptAnno/transcriptAnno-hg19-1kb.tsv`)
  - `--wpstype`: WPS type (`L` for long [default], `S` for short)
  - `--threads`, `-t`: Number of processes

### Orientation-aware Fragmentation (OCF)
**Purpose:** Computes orientation-aware cfDNA fragmentation (OCF) values in tissue-specific open chromatin regions.

**Biological context:** OCF (Sun et al., Genome Res 2019) measures the phasing of upstream (U) and downstream (D) fragment ends in open chromatin, informing tissue-of-origin of cfDNA.

**Usage:**
```bash
krewlyzer ocf motif_out --output ocf_out [options]
```
- Input: `.bed.gz` files from `motif` command
- Output: `.sync.end` files per tissue and summary `all.ocf.csv` per sample
- Options:
  - `--ocr-input`, `-r`: Open chromatin region BED (default: `data/OpenChromatinRegion/7specificTissue.all.OC.bed`)
  - `--threads`, `-t`: Number of processes

### Fragment-level Methylation (UXM)
**Purpose:** Computes the proportions of Unmethylated (U), Mixed (X), and Methylated (M) fragments per region, supporting both single-end (SE) and paired-end (PE) BAMs.

**Biological context:** Fragment-level methylation (UXM, Sun et al., Nature 2023) reveals cell-of-origin and cancer-specific methylation patterns in cfDNA. Krewlyzer supports both SE and PE mode, pairing reads as in cfDNAFE.

**Usage:**
```bash
# Single-end (default)
krewlyzer uxm /path/to/bam_folder --output uxm_out [options]
# Paired-end mode
krewlyzer uxm /path/to/bam_folder --output uxm_out --type PE [options]
```
- Input: Folder of sorted, indexed BAMs
- Output: `.UXM.tsv` file per sample
- Options:
  - `--mark-input`, `-m`: Marker BED file (default: `data/MethMark/Atlas.U25.l4.hg19.bed`)
  - `--map-quality`, `-q`: Minimum mapping quality (default: 30)
  - `--min-cpg`, `-c`: Minimum CpG per fragment (default: 4)
  - `--methy-threshold`, `-tM`: Methylation threshold (default: 0.75)
  - `--unmethy-threshold`, `-tU`: Unmethylation threshold (default: 0.25)
  - `--type`: Fragment type: SE or PE (default: SE)
  - `--threads`, `-t`: Number of processes

### Mutant Fragment Size Distribution (mFSD)
**Purpose:** Compares the size distribution of mutant vs. wild-type fragments at variant sites with comprehensive 4-way classification.

**Biological context:** Mutant ctDNA fragments are typically shorter (~145bp) than wild-type cfDNA (~166bp). This module quantifies this difference for all small variant types, providing sensitive markers for ctDNA detection and MRD monitoring.

**Variant Types Supported:** SNV, MNV, Insertion, Deletion, Complex

**Fragment Classification:** REF (reference), ALT (alternate), NonREF (errors), N (low quality)

**Usage:**
```bash
krewlyzer mfsd sample.bam --input-file variants.vcf --output output_dir/ [options]
```
- Input: BAM file and VCF/MAF file containing variants.
- Output: 
  - `{sample}.mFSD.tsv` - 39-column summary (counts, mean sizes, 6 pairwise KS tests, derived metrics)
  - `{sample}.mFSD.distributions.tsv` - Per-variant size distributions (with `-d`)
- Options:
  - `--input-file`, `-i`: VCF or MAF file (required)
  - `--mapq`, `-q`: Minimum mapping quality (default: 20)
  - `--output-distributions`, `-d`: Generate per-variant size distributions
  - `--verbose`, `-v`: Enable debug logging
  - `--threads`, `-t`: Number of threads (0=all cores)


### Run All Features
Runs all feature extraction commands (extract, motif, fsc, fsr, fsd, wps, ocf) for a single BAM file in one unified pass. Optional: uxm (requires bisulfite BAM), mfsd (requires variants file).

**Usage:**
```bash
krewlyzer run-all sample.bam --reference hg19.fa --output output_dir/ \
    [--variants variants.maf] [--bisulfite-bam bs.bam] [--threads 4] [--debug]
```
- Options:
  - `--reference`, `-g`: Reference FASTA (required)
  - `--variants`, `-v`: VCF/MAF for mFSD (optional)
  - `--bisulfite-bam`: Bisulfite BAM for UXM (optional)
  - `--mapq`, `-q`: Minimum mapping quality (default: 20)
  - `--threads`, `-t`: Number of threads (0=all cores)
  - `--debug`: Enable debug logging

---

## Output Structure Examples

After `krewlyzer run-all`:
```
output_dir/
├── sample.bed.gz               # Fragment file (Tabix indexed)
├── sample.bed.gz.tbi
├── sample.EndMotif.tsv         # End motif frequencies
├── sample.BreakPointMotif.tsv
├── sample.MDS.tsv              # Motif Diversity Score
├── sample.FSC.tsv              # Fragment Size Coverage
├── sample.FSR.tsv              # Fragment Size Ratio
├── sample.FSD.tsv              # Fragment Size Distribution
├── sample.WPS.tsv.gz           # Windowed Protection Score
├── sample.OCF.tsv              # Orientation-aware Fragmentation summary
├── sample.OCF.sync.tsv         # OCF details
├── sample.mFSD.tsv             # Mutant Fragment Size Distribution (39 columns)
└── sample.mFSD.distributions.tsv  # Per-variant size distributions (optional)
```

---


## Troubleshooting
- **FileNotFoundError:** Ensure all input files/paths exist and are readable. Use absolute paths if possible.
- **PermissionError:** Check output directory permissions.
- **Missing dependencies:** Use Docker or follow [Installation](#installation) for all requirements.
- **Reference mismatch:** BAM and reference FASTA must be from the same genome build.
- **Memory errors:** Use ≥16GB RAM for large BAMs or process in batches.

---

## Citation & Acknowledgements
If you use Krewlyzer in your work, please cite this repository and cfDNAFE. Krewlyzer implements or adapts methods from the following primary literature:

- **DELFI (FSR):** Mouliere F, Chandrananda D, Piskorz AM, et al. Enhanced detection of circulating tumor DNA by fragment size analysis. Sci Transl Med. 2018;10(466):eaat4921. [https://doi.org/10.1126/scitranslmed.aat4921](https://doi.org/10.1126/scitranslmed.aat4921)
- **WPS:** Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA Comprises an In Vivo Nucleosome Footprint that Informs Its Tissues-Of-Origin. Cell. 2016;164(1-2):57-68. [https://doi.org/10.1016/j.cell.2015.11.050](https://doi.org/10.1016/j.cell.2015.11.050)
- **OCF:** Sun K, Jiang P, Chan KC, et al. Orientation-aware plasma cell-free DNA fragmentation analysis in open chromatin regions informs tissue of origin. Genome Res. 2019;29(3):418-427. [https://doi.org/10.1101/gr.242719.118](https://doi.org/10.1101/gr.242719.118)
- **UXM:** Sun K, et al. Fragment-level methylation measures cell-of-origin and cancer-specific signals in cell-free DNA. Nature. 2023;616(7956):563-571. [https://doi.org/10.1038/s41586-022-05580-6](https://doi.org/10.1038/s41586-022-05580-6)

- **cfDNAFE:**
```
@misc{cfDNAFE,
  author = {Wanxin Cui et al.},
  title = {cfDNAFE: A toolkit for comprehensive cell-free DNA fragmentation feature extraction},
  year = {2022},
  howpublished = {\url{https://github.com/Cuiwanxin1998/cfDNAFE}}
}
```
- Developed by the MSK-ACCESS team at Memorial Sloan Kettering Cancer Center.

## References
1. Mouliere F, Chandrananda D, Piskorz AM, et al. Enhanced detection of circulating tumor DNA by fragment size analysis. Sci Transl Med. 2018;10(466):eaat4921. [https://doi.org/10.1126/scitranslmed.aat4921](https://doi.org/10.1126/scitranslmed.aat4921)
2. Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA Comprises an In Vivo Nucleosome Footprint that Informs Its Tissues-Of-Origin. Cell. 2016;164(1-2):57-68. [https://doi.org/10.1016/j.cell.2015.11.050](https://doi.org/10.1016/j.cell.2015.11.050)
3. Sun K, Jiang P, Chan KC, et al. Orientation-aware plasma cell-free DNA fragmentation analysis in open chromatin regions informs tissue of origin. Genome Res. 2019;29(3):418-427. [https://doi.org/10.1101/gr.242719.118](https://doi.org/10.1101/gr.242719.118)
4. Sun K, et al. Fragment-level methylation measures cell-of-origin and cancer-specific signals in cell-free DNA. Nature. 2023;616(7956):563-571. [https://doi.org/10.1038/s41586-022-05580-6](https://doi.org/10.1038/s41586-022-05580-6)


## License

This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0). See the [LICENSE](./LICENSE) file for full terms.
