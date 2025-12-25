# Fragment Extraction

**Command**: `krewlyzer extract`

## Purpose
The `extract` module serves as the entry point for most analysis workflows. It processes a BAM file to extract valid cell-free DNA (cfDNA) fragments and saves them in a standardized, compressed BED format with GC content. It also generates metadata and optional GC correction factors.

## Biological Context
Raw sequencing data (BAM) contains reads that must be paired and filtered to reconstruct physical DNA fragments. This step standardizes the data, removing PCR duplicates and low-quality mappings, ensuring downstream analysis focuses on high-confidence unique molecules.

## Usage
```bash
krewlyzer extract sample.bam -r hg19.fa -o output_dir/ [options]
```

## Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--reference` | `-r` | PATH | *required* | Reference genome FASTA (indexed) |
| `--output` | `-o` | PATH | *required* | Output directory |
| `--genome` | `-G` | TEXT | hg19 | Genome build for GC assets |
| `--gc-correct` | | FLAG | True | Compute GC correction factors |
| `--exclude-regions` | `-x` | PATH | | BED file of regions to exclude |
| `--mapq` | `-q` | INT | 20 | Minimum mapping quality |
| `--minlen` | | INT | 65 | Minimum fragment length |
| `--maxlen` | | INT | 400 | Maximum fragment length |
| `--skip-duplicates` | | FLAG | True | Skip duplicate reads |
| `--require-proper-pair` | | FLAG | True | Require proper read pairs |
| `--chromosomes` | | TEXT | | Comma-separated chromosomes to process |
| `--sample-name` | `-s` | TEXT | | Override sample name |
| `--verbose` | `-v` | FLAG | | Enable verbose logging |
| `--threads` | `-t` | INT | 0 | Number of threads (0=all) |

## Output Files

### 1. Fragment File (`{sample}.bed.gz`)
A block-gzipped, tabix-indexed BED file containing fragment coordinates with GC content.

- **Format**: BED4 (chrom, start, end, gc_content)
- **Coordinates**: 0-based, half-open (standard BED)

### 2. Tabix Index (`{sample}.bed.gz.tbi`)
Index for fast random access to fragment regions.

### 3. Metadata File (`{sample}.metadata.json`)
JSON file with run statistics and configuration.

```json
{
  "sample_id": "CasePlasma",
  "total_fragments": 8123456,
  "filters": {
    "mapq": 20,
    "min_length": 65,
    "max_length": 400
  },
  "timestamp": "2024-12-25T10:30:00.123456"
}
```

### 4. GC Correction Factors (`{sample}.correction_factors.csv`)
Per-GC-bin correction factors for downstream normalization (generated when `--gc-correct` is enabled).

| Column | Description |
|--------|-------------|
| gc_bin | GC content bin (0.00-1.00) |
| short_factor | Correction for short fragments |
| intermediate_factor | Correction for intermediate fragments |
| long_factor | Correction for long fragments |

