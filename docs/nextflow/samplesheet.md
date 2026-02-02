# Samplesheet Format

The Nextflow pipeline accepts a CSV samplesheet with the following columns:

## Columns

```csv
sample,bam,meth_bam,vcf,bed,maf,single_sample_maf,assay,pon,targets
```

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `sample` | TEXT | ✓ | Sample identifier |
| `bam` | PATH | ✓ | WGS/Panel BAM file |
| `meth_bam` | PATH | | Bisulfite BAM for UXM |
| `vcf` | PATH | | VCF for mFSD |
| `bed` | PATH | | Pre-extracted .bed.gz |
| `maf` | PATH | | MAF for mFSD |
| `single_sample_maf` | BOOL | | Skip MAF filtering if `true` |
| `assay` | TEXT | | Assay code: `XS1`, `XS2`, `WGS` |
| `pon` | PATH | | Sample-specific PON (overrides assay) |
| `targets` | PATH | | Sample-specific targets (overrides assay) |

## Input Logic

| Input Combination | Workflow |
|-------------------|----------|
| `bam` only | Full run-all (extract → features) |
| `bam` + `vcf`/`maf` | run-all + mFSD |
| `meth_bam` only | UXM methylation |
| `bed` only | FSC, FSR, FSD, WPS, OCF (no extract) |

## Assay Resolution

When `assay` is set and `--asset_dir` is provided:

| Assay | PON | Targets | Genes |
|-------|-----|---------|:-----:|
| `XS1` | `xs1.pon.parquet` | `xs1.targets.bed` | 128 |
| `XS2` | `xs2.pon.parquet` | `xs2.targets.bed` | 146 |
| `WGS` | `wgs.pon.parquet` | None | - |

## Example

```csv
sample,bam,meth_bam,vcf,bed,maf,single_sample_maf,assay,pon,targets
# MSK-ACCESS samples (auto-resolve PON/targets)
ACCESS_001,/data/sample1.bam,,,,,false,XS1,,
ACCESS_002,/data/sample2.bam,,,,/data/cohort.maf,false,XS2,,

# WGS (no targets)
WGS_001,/data/wgs.bam,,/data/wgs.vcf,,,,WGS,,

# Custom PON/targets
CUSTOM,/data/custom.bam,,,,,,,/data/custom.pon.parquet,/data/custom.bed
```
