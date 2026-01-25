# Krewlyzer Reference Data

This folder contains bundled reference genome and annotation files required for cfDNA feature extraction. Assets are auto-loaded based on the `--genome` flag (hg19/hg38).

> **See Also**: [Installation Guide](https://msk-access.github.io/krewlyzer/installation/) | [Glossary](https://msk-access.github.io/krewlyzer/glossary/)

---

## Directory Structure

```
data/
├── ChromosomeArms/      # Chromosome arm definitions (p/q arms)
├── ChromosomeBins/      # 100kb genomic windows for FSC/FSR
├── OpenChromatinRegion/ # OCF tissue-specific regions
├── WpsAnchors/          # WPS anchor regions (TSS + CTCF)
│   └── GRCh37/
│       ├── hg19.wps_anchors.bed.gz     # Genome-wide
│       ├── xs1.wps_anchors.bed.gz      # MSK-ACCESS v1 (1,611 anchors)
│       └── xs2.wps_anchors.bed.gz      # MSK-ACCESS v2 (1,820 anchors)
├── WpsBackground/       # Alu elements for background stacking
├── MethMark/            # Methylation markers for UXM
├── gc/                  # GC reference files for correction
├── exclude-regions/     # Blacklist regions (centromeres, etc.)
├── genes/               # Gene-grouped BED files for panels
│   └── GRCh37/
│       ├── xs1.genes.bed.gz            # MSK-ACCESS v1 (128 genes)
│       └── xs2.genes.bed.gz            # MSK-ACCESS v2 (146 genes)
├── targets/             # Target region BEDs for panels
│   └── GRCh37/
│       ├── xs1.targets.bed             # MSK-ACCESS v1 targets
│       └── xs2.targets.bed             # MSK-ACCESS v2 targets
├── pon/                 # Panel of Normals models
├── TFBS/                # Transcription factor binding sites for region entropy
└── ATAC/                # ATAC-seq peaks for region entropy
```

---

## Asset Descriptions

### ChromosomeArms
Chromosome arm boundary BED files for FSD arm-level analysis.
- `hg19.arms.bed` - GRCh37/hg19 arm coordinates
- `hg38.arms.bed` - GRCh38/hg38 arm coordinates

### ChromosomeBins
100kb genomic windows for FSC/FSR bin-based counting.
- `hg19_window_100kb.bed`
- `hg38_window_100kb.bed`

### WpsAnchors
Combined TSS + CTCF anchor regions for WPS foreground profiling.
- `hg19.wps_anchors.bed.gz`
- `hg38.wps_anchors.bed.gz`

### WpsBackground
Alu element coordinates for WPS background hierarchical stacking.
- `hg19.alu_consensus.bed.gz`
- `hg38.alu_consensus.bed.gz`

### OpenChromatinRegion
Tissue-specific open chromatin regions for OCF analysis.
- `7specificTissue.all.OC.bed` - 7 tissue types (hg19 only)

> **Note**: OCF regions currently hg19 only. Use `--ocr-input` for custom hg38 regions.

### MethMark
Methylation marker BED files for UXM cell-type deconvolution.
- `Atlas.U25.l4.hg19.bed` / `..hg38.bed` - 25 marker regions
- `Atlas.U250.l4.hg19.bed` / `..hg38.bed` - 250 marker regions

### gc
GC content reference files for LOESS bias correction.
- `ref_genome_GC_hg19.parquet` / `ref_genome_GC_hg38.parquet`

### exclude-regions
Blacklist regions excluded from analysis (centromeres, telomeres, low-mappability).
- `hg19-blacklist.v2.bed.gz`
- `hg38-blacklist.v2.bed.gz`

### pon
Panel of Normals (PON) models for z-score normalization.
- `GRCh37/xs1.pon.parquet` - MSK-ACCESS v1 (xs1)
- `GRCh37/xs2.pon.parquet` - MSK-ACCESS v2 (xs2)

### targets
Target region BED files for panel sequencing.
- `xs1.targets.bed` - MSK-ACCESS v1 capture regions
- `xs2.targets.bed` - MSK-ACCESS v2 capture regions

### genes
Gene-grouped target regions for gene-centric analysis.
- `xs1.genes.bed.gz` - MSK-ACCESS v1 (128 genes)
- `xs2.genes.bed.gz` - MSK-ACCESS v2 (146 genes)

---

## Auto-Resolution

When using `--genome hg19` or `--genome hg38`, Krewlyzer automatically resolves:

```bash
# These are equivalent:
krewlyzer fsc -i sample.bed.gz -o out/ --genome hg19
krewlyzer fsc -i sample.bed.gz -o out/ \
    --bin-input data/ChromosomeBins/hg19_window_100kb.bed
```

---

## Custom Data Files

Override bundled assets with your own files:

```bash
# Custom bins
krewlyzer fsc -i sample.bed.gz -o out/ --bin-input /path/to/custom.bed

# Custom PON
krewlyzer fsr -i sample.bed.gz -o out/ -P /path/to/custom.pon.parquet

# Custom targets (panel mode)
krewlyzer run-all -i sample.bam -r ref.fa -o out/ \
    --target-regions /path/to/panel_targets.bed
```

---

## Data Sources

- **Blacklists**: [ENCODE Blacklist v2](https://github.com/Boyle-Lab/Blacklist)
- **Chromosome arms**: [UCSC cytoBand](https://hgdownload.soe.ucsc.edu/goldenPath/)
- **Alu elements**: [RepeatMasker](http://www.repeatmasker.org/)
- **Methylation markers**: [Loyfer et al., Nature 2022](https://doi.org/10.1038/s41586-022-04480-x)

---

*Reference: Asset structure inspired by [cfDNAFE](https://github.com/Cuiwanxin1998/cfDNAFE).*
