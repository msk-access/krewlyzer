# Nextflow Outputs

Output files produced by the Krewlyzer Nextflow pipeline.

## WGS Output Directory

```
results/
├── {sample}.bed.gz                         # Extracted fragments
├── {sample}.bed.gz.tbi                     # Tabix index
├── {sample}.metadata.json                  # Run metadata and QC metrics
├── {sample}.correction_factors.tsv         # GC correction factors
├── {sample}.EndMotif.tsv                   # 4-mer end motif frequencies
├── {sample}.EndMotif1mer.tsv               # 1-mer end motif + Jagged Index
├── {sample}.BreakPointMotif.tsv            # Breakpoint motif frequencies
├── {sample}.MDS.tsv                        # Motif Diversity Score
├── {sample}.FSC.tsv                        # Fragment size coverage (off-target)
├── {sample}.FSR.tsv                        # Fragment size ratio
├── {sample}.FSD.tsv                        # Fragment size distribution (per arm)
├── {sample}.WPS.parquet                    # WPS nucleosome profiles (foreground)
├── {sample}.WPS_background.parquet         # WPS Alu stacking (background)
├── {sample}.OCF.tsv                        # OCF tissue-of-origin scores
├── {sample}.OCF.sync.tsv                   # OCF sync scores (detail)
├── {sample}.TFBS.tsv                       # TFBS entropy (808 TFs)
├── {sample}.ATAC.tsv                       # ATAC entropy (23 cancer types)
└── {sample}.features.json                  # Unified ML features (--generate_json)
```

## Panel Mode Outputs (with `--assay` or `--targets`)

When assay or targets are provided, additional on-target and panel-specific files are generated:

```
results/
├── ... (all WGS outputs above) ...
│
│── # GC Correction
├── {sample}.correction_factors.ontarget.tsv    # On-target GC factors
│
│── # Motif (on-target split)
├── {sample}.EndMotif.ontarget.tsv
├── {sample}.BreakPointMotif.ontarget.tsv
├── {sample}.MDS.ontarget.tsv
│
│── # FSC (gene-centric + regions)
├── {sample}.FSC.ontarget.tsv                   # On-target FSC
├── {sample}.FSC.gene.tsv                       # Gene-level FSC (e.g., 146 genes for xs2)
├── {sample}.FSC.regions.tsv                    # Per-exon/target FSC
├── {sample}.FSC.regions.e1only.tsv             # E1-only FSC (first exon per gene)
│
│── # FSD (on-target split)
├── {sample}.FSD.ontarget.tsv
│
│── # WPS (panel-specific anchors)
├── {sample}.WPS.panel.parquet                  # Panel gene WPS profiles
│
│── # OCF (on/off-target + panel-filtered OCRs)
├── {sample}.OCF.ontarget.tsv
├── {sample}.OCF.ontarget.sync.tsv
├── {sample}.OCF.offtarget.tsv
├── {sample}.OCF.offtarget.sync.tsv
│
│── # TFBS/ATAC (on-target + sync)
├── {sample}.TFBS.ontarget.tsv
├── {sample}.TFBS.sync.tsv
├── {sample}.TFBS.ontarget.sync.tsv
├── {sample}.ATAC.ontarget.tsv
├── {sample}.ATAC.sync.tsv
├── {sample}.ATAC.ontarget.sync.tsv
│
│── # Region MDS (per-gene/exon)
├── {sample}.MDS.exon.tsv                       # Per-exon MDS scores
└── {sample}.MDS.gene.tsv                       # Gene-level aggregated MDS
```

## mFSD Variant Outputs (with VCF/MAF in samplesheet)

```
results/
├── {sample}.mFSD.tsv                       # Per-variant mFSD summary
└── {sample}.mFSD.distributions.tsv         # Per-variant size distributions (optional)
```

## UXM Methylation Outputs (with `meth_bam` in samplesheet)

```
results/
└── {sample}.UXM.tsv                        # Fragment-level methylation
```

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
