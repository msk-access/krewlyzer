# Orientation-aware Fragmentation (OCF)

**Command**: `krewlyzer ocf`

## Purpose
Computes orientation-aware cfDNA fragmentation (OCF) values in tissue-specific open chromatin regions.

## Biological Context
OCF (Sun et al., Genome Res 2019) measures the phasing of upstream (U) and downstream (D) fragment ends in open chromatin, informing tissue-of-origin of cfDNA.

## Usage
```bash
krewlyzer ocf motif_out --output ocf_out [options]
```

## Options
- `--ocr-input`, `-r`: Open chromatin region BED (default: `data/OpenChromatinRegion/7specificTissue.all.OC.bed`)
- `--threads`, `-t`: Number of processes
