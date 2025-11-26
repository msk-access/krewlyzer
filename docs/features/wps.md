# Windowed Protection Score (WPS)

**Command**: `krewlyzer wps`

## Purpose
Computes nucleosome protection scores (WPS) for each region in a transcript/region file.

## Biological Context
The WPS (Snyder et al., 2016) quantifies nucleosome occupancy and chromatin accessibility by comparing fragments spanning a window to those ending within it.
- **High WPS**: Nucleosome protection.
- **Low WPS**: Open chromatin / nucleosome depletion.

## Usage
```bash
krewlyzer wps motif_out --output wps_out [options]
```

## Options
- `--tsv-input`: Transcript region file (default: `data/TranscriptAnno/transcriptAnno-hg19-1kb.tsv`)
- `--wpstype`: WPS type (`L` for long [default], `S` for short)
- `--threads`, `-t`: Number of processes
