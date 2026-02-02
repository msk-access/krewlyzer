# Nextflow Outputs

Output files produced by the Nextflow pipeline.

## Output Directory Structure

```
results/
├── extract/
│   ├── {sample}.bed.gz
│   └── {sample}.bed.gz.tbi
├── motif/
│   ├── {sample}.EndMotif.tsv
│   └── {sample}.MDS.tsv
├── fsc/
│   ├── {sample}.FSC.tsv
│   └── {sample}.FSC.gene.tsv
├── fsd/
│   └── {sample}.FSD.tsv
├── fsr/
│   └── {sample}.FSR.tsv
├── wps/
│   ├── {sample}.WPS.parquet
│   └── {sample}.WPS_background.parquet
├── ocf/
│   └── {sample}.OCF.tsv
├── region_entropy/
│   ├── {sample}.TFBS.tsv
│   └── {sample}.ATAC.tsv
├── region_mds/
│   ├── {sample}.MDS.exon.tsv
│   └── {sample}.MDS.gene.tsv
├── mfsd/
│   └── {sample}.mFSD.tsv
└── uxm/
    └── {sample}.UXM.tsv
```

## Panel Mode Outputs

When `--targets` is provided, outputs are split:

| File | Content |
|------|---------|
| `{sample}.FSC.tsv` | Off-target features |
| `{sample}.FSC.ontarget.tsv` | On-target features |

## Available Modules

| Module | Description |
|--------|-------------|
| `KREWLYZER_EXTRACT` | Fragment extraction |
| `KREWLYZER_FSC` | Fragment Size Coverage |
| `KREWLYZER_FSR` | Fragment Size Ratio |
| `KREWLYZER_FSD` | Fragment Size Distribution |
| `KREWLYZER_WPS` | Windowed Protection Score |
| `KREWLYZER_OCF` | Orientation cfDNA Fragmentation |
| `KREWLYZER_MOTIF` | End Motif & MDS |
| `KREWLYZER_REGION_ENTROPY` | TFBS/ATAC entropy |
| `KREWLYZER_REGION_MDS` | Per-gene MDS |
| `KREWLYZER_UXM` | Methylation |
| `KREWLYZER_MFSD` | Mutant Fragment Size |
| `KREWLYZER_RUNALL` | Full pipeline |
| `KREWLYZER_BUILD_PON` | Build PON |
