# Fragment Size Ratio (FSR)

**Command**: `krewlyzer fsr`

## Purpose
Calculates the ratio of fragment counts in specific size bins relative to the total count per genomic window.

## Biological Context
The DELFI method (Mouliere et al., 2018) showed that cfDNA fragment size ratios are highly informative for cancer detection. Krewlyzer calculates ratios for:
- **Ultra-Short**: 65-100bp (Highly specific for ctDNA)
- **Short**: 65-150bp
- **Intermediate**: 151-260bp
- **Long**: 261-400bp

## Usage
```bash
krewlyzer fsr motif_out --output fsr_out [options]
```

## Options
- `--bin-input`, `-b`: Bin file (default: `data/ChormosomeBins/hg19_window_100kb.bed`)
- `--windows`, `-w`: Window size (default: 100000)
- `--continue-n`, `-c`: Super-bin size (default: 50)
- `--threads`, `-t`: Number of processes

## Targeted Panel Usage
For targeted sequencing (ACCESS, etc.), use custom regions with `windows=1`:
```bash
krewlyzer fsr motif_out -b targets.bed -w 1 -c 1 --output fsr_out
```

