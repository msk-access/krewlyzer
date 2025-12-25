# Panel of Normals (PON) Models

This directory contains pre-built PON models for supported assays.

## Purpose

PON models provide baseline fragmentomics patterns from healthy cfDNA samples.
These are used for:
- **GC correction**: Assay-specific GC bias curves
- **FSD normalization**: Expected size distributions per chromosome arm
- **WPS baseline**: Expected nucleosome positioning per transcript region

## Available Models

| Model | Assay | Samples | Status |
|-------|-------|---------|--------|
| `msk-access-v1.pon.parquet` | MSK-ACCESS v1 | TBD | 🚧 Pending |
| `msk-access-v2.pon.parquet` | MSK-ACCESS v2 | TBD | 🚧 Pending |

## Usage

```bash
krewlyzer run-all sample.bam --pon-model src/krewlyzer/data/pon/msk-access-v2.pon.parquet -o out/
```

## Building Custom PON

```bash
krewlyzer build-pon \
    --sample-list healthy_plasma.txt \
    --assay my-assay \
    --reference hg19.fa \
    --output my-assay.pon.parquet
```

See `krewlyzer build-pon --help` for details.

## Model Format

PON models use Parquet format containing:
- `gc_bias`: GC correction curves for short/intermediate/long fragments
- `fsd_baseline`: Size distribution per chromosome arm
- `wps_baseline`: Mean/std WPS per transcript region
