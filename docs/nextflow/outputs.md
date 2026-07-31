# Nextflow Outputs

Output files produced by the Krewlyzer Nextflow pipeline.

## WGS Output Directory

```
results/
├── {sample}.bed.gz                         # Extracted fragments
├── {sample}.bed.gz.tbi                     # Tabix index
├── {sample}.metadata.tsv                   # Run metadata and QC metrics (tabular)
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

## Validation Artifacts (Parquet runs)

Written per sample by `run-all`, and gathered once per cohort:

| File | Scope | Contents |
|------|-------|----------|
| `{sample}.validation.json` | sample | Contract findings for that sample |
| `{sample}.fingerprint.json` | sample | ~20 KB summary — a hash and two counts per column |
| `cohort.validation.json` | cohort | Cross-sample degeneracy findings |

The split exists because **degeneracy is inherently cross-sample**: every sample
can pass on its own while a metric is constant across all of them. A sample
directory is ~1.5 GB, so the gather step compares fingerprints rather than
re-reading tables.

Skipped entirely for tsv-only runs — see `--output_format` in
[Parameters](parameters.md#output-parameters).

!!! note "Tool-level mode"
    `--use_runall false` produces no fingerprints, so cohort validation is
    skipped rather than run on a partial set, which would report degeneracy
    findings that are artefacts of the missing samples.

---

## Output Changes in 0.9.0

Six output families changed value semantics. Values are **not comparable**
across this boundary, and PON baselines built on uncollapsed input should be
rebuilt.

| Output | Change |
|--------|--------|
| Every positional family — `WPS`, `OCF`, TFBS/ATAC, `MDS.exon/gene`, `FSC.gene/regions` | Fragment coordinates corrected when R1 is the rightmost mate (~48% of reads on uncollapsed input) |
| `FSC.gene.*`, `FSC.regions.*` | Size bands aligned to the genome-bin bands; new `ultra_long` and `ultra_long_ratio` |
| `FSC.regions.*` | New `strand`, `is_e1`, `is_alt_e1` columns |
| `FSC.regions.e1only.*` | Selects on the E1 flags; genes with no annotated first exon are omitted rather than represented by an internal exon |
| `MDS.gene.mds_e1`, `mds_e1_z` | Strand-aware on panel data for the first time; `NaN` where the gene has no E1, instead of a fabricated `0.0` |
| `WPS_background.*` | NRL family is data-dependent; new `nrl_at_band_limit` marks a right-censored estimate |

---

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
