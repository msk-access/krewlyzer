# Samplesheet Format

The Nextflow pipeline accepts a CSV samplesheet with the following columns:

## Columns

```csv
sample,bam,mfsd_bam,meth_bam,vcf,bed,maf,single_sample_maf,assay,pon,targets
```

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `sample` | TEXT | ✓ | Sample identifier |
| `bam` | PATH | ✓ | WGS/Panel BAM file (all_unique for panels) |
| `mfsd_bam` | PATH | | Duplex BAM for mFSD (falls back to `bam` if empty) |
| `meth_bam` | PATH | | Bisulfite BAM for UXM |
| `vcf` | PATH | | VCF for mFSD, used when `maf` is empty |
| `bed` | PATH | | Pre-extracted .bed.gz |
| `maf` | PATH | | MAF for mFSD (takes precedence over `vcf`) |
| `single_sample_maf` | BOOL | | Skip MAF filtering if `true` |
| `assay` | TEXT | | Assay code: `XS1`, `XS2`, `WGS` |
| `pon` | PATH | | Sample-specific PON (overrides assay) |
| `targets` | PATH | | Sample-specific targets (overrides assay) |

## Input Logic

| Input Combination | Workflow |
|-------------------|----------|
| `bam` only | Full run-all (extract → features) |
| `bam` + `vcf`/`maf` | run-all + mFSD |
| `bam` + `mfsd_bam` + `maf` | run-all with dual BAM support |
| `meth_bam` only | UXM methylation |
| `bed` only | FSC, FSR, FSD, WPS, OCF (no extract) |

### Variants for mFSD

mFSD reads one variant file per sample, in both `--use_runall true` and
tool-level mode:

- **`maf`**: filtered to the sample's `Tumor_Sample_Barcode` rows, unless
  `single_sample_maf` is `true`.
- **`vcf`**: used as given, every record. The pipeline does not filter a VCF
  by sample, so give each row a VCF of that sample's calls.
- **Both set**: the MAF is used and the VCF is ignored, with a warning in the
  log naming the sample.
- **Neither set**: run-all receives a header-only placeholder MAF and records
  `mfsd_status=no_variants` in `{sample}.metadata`; tool-level mode skips the
  mFSD process.

The file type is detected from the extension: `.vcf` is read as VCF, anything
else as MAF. Use an uncompressed `.vcf` in run-all mode for now: run-all counts
the variant lines with a plain-text read before mFSD starts, and a `.vcf.gz`
fails that count.

## Assay Resolution

When `assay` is set and `--asset_dir` is provided:

| Assay | PON | Targets | Genes |
|-------|-----|---------|:-----:|
| `XS1` | `xs1.pon.parquet` | `xs1.targets.bed` | 128 |
| `XS2` | `xs2.pon.parquet` | `xs2.targets.bed` | 146 |
| `WGS` | `wgs.pon.parquet` | None | - |

## Example

```csv
sample,bam,mfsd_bam,meth_bam,vcf,bed,maf,single_sample_maf,assay,pon,targets
# MSK-ACCESS V1 samples (auto-resolve PON/targets)
ACCESS_001,/data/sample1.bam,,,,,,,XS1,,
ACCESS_002,/data/sample2.bam,,,,,,/data/cohort.maf,false,XS2,,

# MSK-ACCESS V2 with duplex BAM for mFSD
ACCESS_V2,/data/sample.all_unique.bam,/data/sample.duplex.bam,,,,,maf.tsv,true,XS2,,

# WGS (no targets)
WGS_001,/data/wgs.bam,,,/data/wgs.vcf,,,,WGS,,

# Custom PON/targets
CUSTOM,/data/custom.bam,,,,,,,,,/data/custom.pon.parquet,/data/custom.bed
```

