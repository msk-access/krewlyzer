# Fragment Size Coverage (FSC)

**Command**: `krewlyzer fsc`

## Purpose
Computes z-scored coverage of cfDNA fragments in different size ranges, per genomic bin (default: 100kb), with GC correction.

## Biological Context
cfDNA fragment size profiles are informative for cancer detection and tissue-of-origin. FSC quantifies the coverage of:
- **Short**: 65-150bp
- **Intermediate**: 151-260bp
- **Long**: 261-400bp
- **Total**: 65-400bp

Values are normalized to genome-wide means.

## Usage
```bash
krewlyzer fsc motif_out --output fsc_out [options]
```

## Options
- `--bin-input`, `-b`: Bin file (default: `data/ChormosomeBins/hg19_window_100kb.bed`)
- `--windows`, `-w`: Window size (default: 100000)
- `--continue-n`, `-c`: Super-bin size (default: 50)
- `--threads`, `-t`: Number of processes
