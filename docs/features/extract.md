# Fragment Extraction

**Command**: `krewlyzer extract`

## Purpose
The `extract` module serves as the entry point for most analysis workflows. It processes a BAM file to extract valid cell-free DNA (cfDNA) fragments and saves them in a standardized, compressed BED format. It also generates a JSON metadata file containing processing statistics and configuration.

## Biological Context
Raw sequencing data (BAM) contains reads that must be paired and filtered to reconstruct physical DNA fragments. This step standardizes the data, removing PCR duplicates and low-quality mappings, ensuring downstream analysis focuses on high-confidence unique molecules.

## Usage
```bash
krewlyzer extract sample.bam -g hg19.fa -o output_dir/ [options]
```

## Options
- `--genome`, `-g`: Reference genome FASTA (required).
- `--output`, `-o`: Output directory.
- `--mapq`, `-q`: Minimum mapping quality (default: 20).
- `--minlen`: Minimum fragment length (default: 65).
- `--maxlen`: Maximum fragment length (default: 400).
- `--exclude-regions`, `-x`: BED file of regions to blacklist (e.g., centromeres).

## Output Files

### 1. Fragment File (`{sample}.bed.gz`)
A block-gzipped, tabix-indexed BED file containing the coordinates of extracted fragments.
- **Format**: BED3+3 (chrom, start, end, name, score, strand) - *Note: currently simplified to BED3 or similar depending on implementation.*
- **Coordinates**: 0-based, half-open (standard BED).

### 2. Metadata File (`{sample}.metadata.json`)
A JSON file capturing run provenance and statistics.

**Example Structure:**
```json
{
  "sample_id": "CasePlasma",
  "total_fragments": 8123456,
  "filters": {
    "mapq": 20,
    "min_length": 65,
    "max_length": 400
  },
  "timestamp": "2023-10-27T10:30:00.123456"
}
```

**Fields:**
- `sample_id`: Identifier derived from input filename or user argument.
- `total_fragments`: Number of valid fragments extracted after filtering.
- `filters`: Configuration parameters used for the run.
- `timestamp`: Execution time.
