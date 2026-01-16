# Target Regions Directory

This directory contains panel target regions for different assays.

## Expected Files

| Assay Code | Filename | Description |
|------------|----------|-------------|
| `XS1` | `XS1_targets.bed` | MSK-ACCESS V1 panel targets |
| `XS2` | `XS2_targets.bed` | MSK-ACCESS V2 panel targets |

## File Format

Standard BED format (tab-separated, no header):

```
chr1    100000    100500    gene1_exon1
chr1    105000    105300    gene1_exon2
chr2    200000    200800    gene2_exon1
```

## How Targets Are Used

When `--target-regions` is provided:

1. **GC Correction** - Model trained on **off-target** reads only (unbiased)
2. **Output Splitting** - Generates:
   - `*.tsv` → Off-target (recommended for ML)
   - `*.ontarget.tsv` → On-target (capture-biased)
3. **FSC/FSR Aggregation** - Disabled for panel data

## Nextflow Auto-Resolution

The pipeline resolves targets based on the `assay` column in samplesheet:

```csv
sample,bam,...,assay
sample1,sample1.bam,...,XS1
sample2,sample2.bam,...,XS2
```

Resolution: `{asset_dir}/targets/{assay}_targets.bed`

## PON Files

PON models are in the adjacent `pon/` directory:
- `msk-access-v1.pon.parquet` (alias: XS1)
- `msk-access-v2.pon.parquet` (alias: XS2)
