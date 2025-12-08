# Usage Guide

## Command Summary

| Command   | Description                                  |
|-----------|----------------------------------------------|
| `motif`   | Motif-based feature extraction               |
| `fsc`     | Fragment size coverage                       |
| `fsr`     | Fragment size ratio                          |
| `fsd`     | Fragment size distribution                   |
| `wps`     | Windowed protection score                    |
| `ocf`     | Orientation-aware fragmentation              |
| `uxm`     | Fragment-level methylation (SE/PE)           |
| `mfsd`    | Mutant fragment size distribution            |
| `run-all` | Run all features for a BAM                   |

## Reference Data
- **Reference Genome (FASTA):**
  - Download GRCh37/hg19 from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/)
  - BAMs must be sorted, indexed, and aligned to the same build
- **Bin/Region/Marker Files:**
  - Provided in `krewlyzer/data/` (see options for each feature)

## Typical Workflow

```bash
# 1. Motif extraction (produces .bed.gz files)
krewlyzer motif sample.bam -g hg19.fa -o motif_out

# 2. Extract additional features from motif output:
krewlyzer fsc motif_out --output fsc_out
krewlyzer fsr motif_out --output fsr_out
krewlyzer fsd motif_out --arms-file krewlyzer/data/ChormosomeArms/hg19_arms.bed --output fsd_out
krewlyzer wps motif_out --output wps_out
krewlyzer ocf motif_out --output ocf_out
krewlyzer uxm /path/to/bam_folder --output uxm_out
krewlyzer mfsd sample.bam --input variants.vcf --output mfsd_out/sample.mfsd.tsv

# 3. Run all features in one call:
krewlyzer run-all sample.bam --reference hg19.fa --output all_features_out --variant-input variants.vcf
```

## Targeted Panel Usage (ACCESS, etc.)

For targeted sequencing panels (e.g., MSK-ACCESS), FSC/FSR require a custom regions BED file instead of the default genome-wide 100kb bins:

```bash
# Using run-all with custom target regions
krewlyzer run-all sample.bam --reference hg19.fa --output out/ \
  --bin-input /path/to/MSK-ACCESS-v2_canonicaltargets.bed

# Or run FSC/FSR individually with target regions
krewlyzer fsc motif_out -b targets.bed -w 1 -c 1 --output fsc_out
krewlyzer fsr motif_out -b targets.bed -w 1 -c 1 --output fsr_out
```

> **Note:** Without `--bin-input`, FSC/FSR will produce zeros for targeted panels since data only covers specific gene regions, not genome-wide bins.

