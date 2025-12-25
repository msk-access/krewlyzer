# Panel of Normals (PON)

PON models provide baseline fragmentomics patterns from healthy cfDNA samples for normalization and z-score computation.

## Purpose

PON models enable:

1. **GC Correction** - Assay-specific GC bias curves
2. **FSD Normalization** - Expected size distributions per chromosome arm
3. **WPS Baseline** - Mean/std nucleosome positioning per transcript

---

## Available Models

| Model | Assay | Status |
|-------|-------|--------|
| `msk-access-v1.pon.parquet` | MSK-ACCESS v1 | 🚧 Pending |
| `msk-access-v2.pon.parquet` | MSK-ACCESS v2 | 🚧 Pending |

Models are stored in `src/krewlyzer/data/pon/`.

---

## Using a PON Model

### With run-all

```bash
krewlyzer run-all sample.bam \
    --reference hg19.fa \
    --pon-model msk-access-v2.pon.parquet \
    --output results/
```

### With Individual Tools

```bash
krewlyzer fsc sample.bed.gz \
    --pon-model msk-access-v2.pon.parquet \
    -o output/

krewlyzer wps sample.bed.gz \
    --pon-model msk-access-v2.pon.parquet \
    -o output/
```

---

## Building a Custom PON

### Step 1: Prepare Sample List

Create a text file with BAM paths (one per line):

```
/path/to/healthy_plasma_001.bam
/path/to/healthy_plasma_002.bam
...
```

### Step 2: Build PON

```bash
krewlyzer build-pon \
    --sample-list healthy_samples.txt \
    --assay my-assay \
    --reference hg19.fa \
    --output my-assay.pon.parquet
```

### Options

| Option | Description |
|--------|-------------|
| `--sample-list` | Text file with BAM paths |
| `--assay` | Assay identifier (e.g., "MSK-ACCESS-v2") |
| `--reference` | Reference FASTA |
| `--output` | Output Parquet file |
| `--genome` | Genome build (hg19/hg38) |
| `--threads` | Number of threads |

---

## PON Model Format

PON models use Parquet format with multiple row groups:

### Metadata Row Group

| Field | Type | Description |
|-------|------|-------------|
| `schema_version` | string | Model format version |
| `assay` | string | Assay identifier |
| `build_date` | string | ISO timestamp |
| `n_samples` | int | Number of samples |
| `reference` | string | Genome build |

### GC Bias Row Group

| Column | Description |
|--------|-------------|
| `gc_bin` | GC content bin (0.0-1.0) |
| `short_expected` | Expected coverage for short fragments |
| `short_std` | Standard deviation |
| `intermediate_expected` | Expected for intermediate |
| `intermediate_std` | Standard deviation |
| `long_expected` | Expected for long fragments |
| `long_std` | Standard deviation |

### FSD Baseline Row Group

| Column | Description |
|--------|-------------|
| `arm` | Chromosome arm (e.g., "chr1p") |
| `size_bin` | Fragment size (5bp bins) |
| `proportion_mean` | Mean proportion |
| `proportion_std` | Standard deviation |

### WPS Baseline Row Group

| Column | Description |
|--------|-------------|
| `region_id` | Transcript/region ID |
| `wps_long_mean` | Mean long WPS |
| `wps_long_std` | Standard deviation |
| `wps_short_mean` | Mean short WPS |
| `wps_short_std` | Standard deviation |

---

## Hybrid GC Correction

When a PON model is provided, Krewlyzer applies **hybrid correction**:

```
Algorithm:
1. PON correction: corrected = observed / pon_expected[gc]
2. Residual LOESS: residual = loess(gc, corrected)  
3. Final: final = corrected / residual
```

This corrects for:
- **Assay-specific bias** (from PON)
- **Sample-specific residual** (within-sample LOESS)

---

## Assay Compatibility

PON models are assay-specific. If your sample's assay differs from the PON:

!!! warning "Assay Mismatch"
    Krewlyzer logs a warning but continues processing. Results may be less accurate.

For best results, use a PON built from the same assay and protocol.
