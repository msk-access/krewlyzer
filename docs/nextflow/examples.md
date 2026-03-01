# Nextflow Examples

Common workflows for batch processing with Krewlyzer.

## Basic Run

```bash
nextflow run msk-access/krewlyzer \
    --samplesheet samples.csv \
    --ref /path/to/hg19.fa \
    --outdir results/
```

## MSK-ACCESS Panel

```bash
nextflow run msk-access/krewlyzer \
    --samplesheet samples.csv \
    --ref /path/to/hg19.fa \
    --asset_dir /path/to/krewlyzer/data/ \
    --outdir results/
```

With samplesheet:
```csv
sample,bam,mfsd_bam,meth_bam,vcf,bed,maf,single_sample_maf,assay,pon,targets
ACCESS_001,/data/sample1.bam,,,,,,false,XS2,,
ACCESS_002,/data/sample2.bam,,,,,,false,XS2,,
```

## With Variant Analysis

```bash
nextflow run msk-access/krewlyzer \
    --samplesheet samples.csv \
    --ref /path/to/hg19.fa \
    --outdir results/
```

With samplesheet:
```csv
sample,bam,mfsd_bam,meth_bam,vcf,bed,maf,single_sample_maf,assay,pon,targets
SAMPLE_001,/data/sample1.bam,,,/data/variants.vcf,,,,XS2,,
SAMPLE_002,/data/sample2.bam,,,,,,/data/cohort.maf,false,XS2,,
```

## Using Docker

```bash
nextflow run msk-access/krewlyzer \
    -profile docker \
    --samplesheet samples.csv \
    --ref /path/to/hg19.fa \
    --outdir results/
```

## SLURM Cluster

```bash
nextflow run msk-access/krewlyzer \
    -profile slurm \
    --samplesheet samples.csv \
    --ref /path/to/hg19.fa \
    --outdir results/
```

## IRIS HPC (MSK)

For running on MSKCC's IRIS cluster, use the [nf-core/configs](https://github.com/nf-core/configs) institutional profile:

```bash
nextflow run msk-access/krewlyzer \
    -profile iris \
    --samplesheet samples.csv \
    --ref /data1/ref/hg19/hg19.fa \
    --outdir /scratch/$GROUP/results/
```

!!! tip
    The `iris` profile automatically configures:
    - SLURM executor with proper queue settings
    - Singularity with pre-cached images at `/data1/core006/resources/singularity_image_library`
    - Scratch paths and work directories

For preemptable (faster) queue:

```bash
nextflow run msk-access/krewlyzer \
    -profile iris \
    --preemptable true \
    --samplesheet samples.csv \
    --ref /data1/ref/hg19/hg19.fa \
    --outdir /scratch/$GROUP/results/
```

## Resume Failed Run

```bash
nextflow run msk-access/krewlyzer \
    -resume \
    --samplesheet samples.csv \
    --ref /path/to/hg19.fa \
    --outdir results/
```
