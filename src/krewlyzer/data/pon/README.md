# Panel of Normals (PON) Models

This directory contains pre-built PON models for MSK-ACCESS cfDNA assays.

## Purpose

PON models provide baseline fragmentomics patterns from healthy cfDNA samples.
These are used for:
- **GC correction**: Assay-specific GC bias curves
- **FSD normalization**: Expected size distributions per chromosome arm
- **WPS baseline**: Expected nucleosome positioning per transcript region
- **OCF baseline**: Open chromatin fragmentation patterns
- **MDS baseline**: Motif diversity score baselines
- **TFBS/ATAC baseline**: Regulatory region entropy baselines

## Available Models

| Model | Assay | Read Type | Size | Tables |
|-------|-------|-----------|------|--------|
| `xs1.duplex.pon.parquet` | MSK-ACCESS v1 | Duplex only | 0.35 MB | gc_bias, ocf, mds, atac |
| `xs1.all_unique.pon.parquet` | MSK-ACCESS v1 | Duplex + Simplex + Singletons | 0.39 MB | gc_bias, ocf, mds, atac, tfbs |
| `xs2.duplex.pon.parquet` | MSK-ACCESS v2 | Duplex only | 0.38 MB | gc_bias, ocf, mds, atac, tfbs |
| `xs2.all_uniq.pon.parquet` | MSK-ACCESS v2 | Duplex + Simplex + Singletons | 0.39 MB | gc_bias, ocf, mds, atac, tfbs |

### Read Type Guidance

| Read Type | Description | Use Case |
|-----------|-------------|----------|
| **Duplex** | High-accuracy reads with both strands sequenced | Low tumor fraction detection, high specificity |
| **All Unique** | Duplex + Simplex + Singletons | Maximum sensitivity, higher noise tolerance |

## Usage

```bash
# Duplex BAM with duplex PON
krewlyzer run-all sample.duplex.bam -r hg19.fa -o out/ \
    --pon-model src/krewlyzer/data/pon/GRCh37/xs2.duplex.pon.parquet

# All-unique BAM with all-unique PON
krewlyzer run-all sample.all_unique.bam -r hg19.fa -o out/ \
    --pon-model src/krewlyzer/data/pon/GRCh37/xs2.all_uniq.pon.parquet
```

## Building Custom PON

```bash
krewlyzer build-pon \
    --sample-list healthy_plasma.txt \
    --reference hg19.fa \
    --output custom.pon.parquet
```

See `krewlyzer build-pon --help` for details.

## Model Format

PON models use Parquet format containing these baseline tables:

| Table | Description |
|-------|-------------|
| `metadata` | PON build info (samples, date, version) |
| `gc_bias` | GC correction curves per fragment length bin |
| `fsd_baseline` | Size distribution per chromosome arm |
| `wps_baseline` | Mean/std WPS per transcript region |
| `ocf_baseline` | OCF mean/std per tissue type |
| `mds_baseline` | Motif diversity mean/std per gene |
| `tfbs_baseline` | TFBS entropy mean/std per TF |
| `atac_baseline` | ATAC entropy mean/std per cancer type |
