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
├── TranscriptAnno/      # TSS annotations for WPS
├── WpsAnchors/          # WPS anchor regions (TSS + CTCF)
├── WpsBackground/       # Alu elements for background stacking
├── MethMark/            # Methylation markers for UXM
├── gc/                  # GC reference files for correction
├── exclude-regions/     # Blacklist regions (centromeres, etc.)
├── pon/                 # Panel of Normals models
├── targets/             # Target region BEDs (MSK-ACCESS)
└── CNVdependency/       # Mappability and GC wig files
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
- `hg19_wps_anchors.bed.gz`
- `hg38_wps_anchors.bed.gz`

### WpsBackground
Alu element coordinates for WPS background hierarchical stacking.
- `hg19_alu_background.bed.gz`
- `hg38_alu_background.bed.gz`

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
- `hg19.gc_reference.parquet` / `hg38.gc_reference.parquet`

### exclude-regions
Blacklist regions excluded from analysis (centromeres, telomeres, low-mappability).
- `hg19-blacklist.v2.bed`
- `hg38-blacklist.v2.bed`

### pon
Panel of Normals (PON) models for z-score normalization.
- `msk-access-v1.pon.parquet` - MSK-ACCESS v1 (XS1)
- `msk-access-v2.pon.parquet` - MSK-ACCESS v2 (XS2)
- `wgs.pon.parquet` - Whole genome sequencing

### targets
Target region BED files for panel sequencing.
- `MSK-ACCESS-v2_targets.bed` - MSK-ACCESS v2 capture regions

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
