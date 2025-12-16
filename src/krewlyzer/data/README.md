# Krewlyzer Reference Data

This folder contains reference genome and annotation files required for motif and feature extraction in Krewlyzer. The structure and file recommendations are inspired by the [cfDNAFE](https://github.com/Cuiwanxin1998/cfDNAFE) project.

## Table of Contents
- [BlackList](#blacklist)
- [CNVdependency](#cnvdependency)
- [ChormosomeArms](#chormosomearms)
- [ChormosomeBins](#chormosomebins)
- [MethMark](#methmark)
- [OpenChromatinRegion](#openchromatinregion)
- [TranscriptAnno](#transcriptanno)

## Subfolders and Contents

### BlackList
Contains blacklist BED files for hg19 and hg38 genomes:
- `hg19-blacklist.v2.bed`
- `hg38-blacklist.v2.bed`

### CNVdependency
Contains GC content and mappability wig files, centromere locations, and sequence info for CNV analysis:
- GC content: `gc_hg19_*.wig`, `gc_hg38_*.wig`
- Mappability: `map_hg19_*.wig`, `map_hg38_*.wig`
- Centromere: `GRCh37.p13_centromere_UCSC-gapTable.txt`, `GRCh38.GCA_000001405.2_centromere_acen.txt`
- Sequence info: `seqinfo_hg19_ucsc.rds`, `seqinfo_hg38_ucsc.rds`

### ChormosomeArms
Chromosome arm annotation BED files:
- `hg19.arms.bed`
- `hg38.arms.bed`

### ChormosomeBins
Genome bin BED files (e.g., 100kb bins):
- `hg19_window_100kb.bed`
- `hg38_window_100kb.bed`

### MethMark
Contains methylation marker BED files:
- `Atlas.U25.l4.hg19.bed`, `Atlas.U25.l4.hg19.full.bed`, `Atlas.U25.l4.hg38.bed`, `Atlas.U25.l4.hg38.full.bed`
- `Atlas.U250.l4.hg19.bed`, `Atlas.U250.l4.hg19.full.bed`, `Atlas.U250.l4.hg38.full.bed`
- `Markers.U250.hg19.bed`, `Markers.U250.hg38.bed`

### OpenChromatinRegion
Contains open chromatin region BED files:
- `7specificTissue.all.OC.bed`

### TranscriptAnno
Transcript annotation tables for hg19 and hg38:
- `transcriptAnno-hg19-10kb.tsv`, `transcriptAnno-hg19-1kb.tsv`
- `transcriptAnno-hg38-10kb.tsv`, `transcriptAnno-hg38-1kb.tsv`

## Downloading Data
You can download reference files from:
- [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/)
- [ENCODE Project](https://www.encodeproject.org/)
- Or use your own custom files as needed.

**Note:** Large files are not included in this repository. Place them here as needed for your analyses.

---

*Reference: This folder structure and file recommendations are inspired by [cfDNAFE](https://github.com/Cuiwanxin1998/cfDNAFE).*
