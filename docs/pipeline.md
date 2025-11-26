# Pipeline Integration

## Run All Features

The `run-all` command allows you to execute all feature extraction modules for a single BAM file in one go.

### Usage
```bash
krewlyzer run-all sample.bam --reference hg19.fa --output all_features_out \
    [--variant-input variants.vcf] \
    [--threads N] \
    [--type SE|PE]
```

### Arguments
- `sample.bam`: Input BAM file (sorted, indexed).
- `--reference`, `-g`: Reference genome FASTA.
- `--output`, `-o`: Output directory.
- `--variant-input`, `-v`: (Optional) VCF or MAF file for `mfsd` analysis.
- `--threads`, `-t`: Number of threads (default: 1).
- `--type`: Fragment type for UXM (SE or PE).

## Nextflow Pipeline

Krewlyzer includes a Nextflow pipeline (`krewlyzer.nf`) for processing multiple samples in parallel on HPC clusters or local machines.

### Usage
```bash
nextflow run krewlyzer.nf --samplesheet samplesheet.csv --ref /path/to/reference.fa
```

### Samplesheet Format
A CSV file with the following columns:
```csv
sample_id,bam,vcf
sample1,/path/to/sample1.bam,/path/to/sample1.vcf
sample2,/path/to/sample2.bam,
```
- `sample_id`: Unique identifier.
- `bam`: Path to BAM file.
- `vcf` (or `maf`): Optional path to variant file.

### Profiles
- `-profile lsf`: For LSF clusters.
- `-profile slurm`: For SLURM clusters.
- `-profile docker`: Run using Docker container.
