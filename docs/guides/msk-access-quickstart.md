# MSK-ACCESS Quickstart

This guide walks through analyzing MSK-ACCESS panel data with Krewlyzer.

## Prerequisites

- Krewlyzer installed (`pip install krewlyzer` or Docker)
- MSK-ACCESS BAM file (sorted, indexed)
- Reference genome (hg19.fa)

---

## Basic Panel Analysis

```bash
krewlyzer run-all -i ACCESS_sample.bam \
    --reference hg19.fa \
    --output results/ \
    --assay xs2 \
    --target-regions MSK-ACCESS-v2_targets.bed
```

### What This Does

| Flag | Effect |
|------|--------|
| `--assay xs2` | Enables MSK-ACCESS v2 optimizations |
| `--target-regions` | Splits outputs into on/off-target |

---

## Output Files

### Standard Outputs

| File | Description |
|------|-------------|
| `sample.bed.gz` | Extracted fragments |
| `sample.FSC.tsv` | Fragment size coverage |
| `sample.FSC.ontarget.tsv` | FSC for on-target regions |
| `sample.FSR.tsv` | Fragment size ratios |
| `sample.FSD.tsv` | Size distribution by arm |
| `sample.WPS.parquet` | Nucleosome protection (genome-wide) |
| `sample.OCF.tsv` | Tissue of origin footprint |
| `sample.EndMotif.tsv` | End motif frequencies |
| `sample.MDS.tsv` | Motif Diversity Score |

### Panel-Specific Outputs (with `--assay`)

| File | Description |
|------|-------------|
| `sample.FSC.gene.tsv` | **Gene-level FSC** (146 genes for xs2) |
| `sample.WPS.panel.parquet` | **Panel-focused WPS** (~2k anchors) |

---

## With PON Normalization

For z-score computation against healthy baselines:

```bash
krewlyzer run-all -i ACCESS_sample.bam \
    --reference hg19.fa \
    --output results/ \
    --assay xs2 \
    --target-regions MSK-ACCESS-v2_targets.bed \
    --pon-model xs2.pon.parquet
```

This adds z-score columns to all outputs (e.g., `z_core_short`, `wps_nuc_z`, `mds_z`).

---

## JSON Output for ML

Generate a unified JSON file for machine learning pipelines:

```bash
krewlyzer run-all -i ACCESS_sample.bam \
    --reference hg19.fa \
    --output results/ \
    --assay xs2 \
    --target-regions MSK-ACCESS-v2_targets.bed \
    --pon-model xs2.pon.parquet \
    --generate-json
```

This creates `sample.features.json` with all features. See [JSON Output](../features/output/json-output.md) for schema details.

---

## Assay Options

| Assay Code | Panel | Genes | WPS Anchors |
|------------|-------|:-----:|:-----------:|
| `xs1` | MSK-ACCESS v1 | 128 | 1,611 |
| `xs2` | MSK-ACCESS v2 | 146 | 1,820 |

---

## Standalone Tool Examples

### FSC with Gene Aggregation

```bash
krewlyzer fsc -i sample.bed.gz -o results/ \
    --assay xs2 \
    --target-regions targets.bed
```

Outputs: `sample.FSC.tsv` + `sample.FSC.gene.tsv`

### WPS with Dual Output

```bash
krewlyzer wps -i sample.bed.gz -o results/ \
    --assay xs2 \
    --target-regions targets.bed
```

Outputs: `sample.WPS.parquet` (genome-wide) + `sample.WPS.panel.parquet` (panel genes)

---

## Building Your Own PON

Create a PON from your healthy plasma samples:

```bash
# Create sample list file
ls /path/to/healthy/*-duplex.bam > healthy_samples.txt

# Build PON (HPC with custom temp directory)
krewlyzer build-pon healthy_samples.txt \
    --assay xs2 \
    --reference hg19.fa \
    --target-regions MSK-ACCESS-v2_targets.bed \
    --temp-dir /scratch/$USER/pon_tmp \
    --output xs2.pon.parquet \
    --threads 16 \
    --verbose
```

> **Tip**: Use `--temp-dir` to specify a directory with more disk space than `/tmp`. Each BAM extraction creates temporary BED.gz files (~100MB each).

---

## Nextflow Batch Processing

For multiple samples, use the Nextflow pipeline:

```bash
nextflow run nextflow/main.nf \
    --samplesheet samples.csv \
    --ref /path/to/hg19.fa \
    --outdir results/
```

### Samplesheet Example

```csv
sample,bam,meth_bam,vcf,bed,maf,single_sample_maf,assay,pon,targets
ACCESS_001,/data/sample1.bam,,,,,false,XS2,,
ACCESS_002,/data/sample2.bam,,,,,false,XS2,,
ACCESS_003,/data/sample3.bam,,,,,false,XS2,,
```

When `assay=XS2` is set, the pipeline automatically resolves:
- PON: `{asset_dir}/pon/GRCh37/xs2.pon.parquet`
- Targets: `{asset_dir}/targets/GRCh37/xs2.targets.bed`

---

## Key Differences from WGS

| Aspect | WGS | MSK-ACCESS |
|--------|-----|------------|
| Coverage | Uniform | Targeted |
| GC Model | All fragments | **Off-target only** |
| FSC Output | Windows | + **Gene-level** |
| WPS Output | Genome-wide | + **Panel-specific** |
| Fragments used | All | Split on/off-target |

---

## Troubleshooting

### "0% of reads pass filters"

Duplex/consensus BAMs need:
```bash
krewlyzer run-all ... --no-require-proper-pair
```

### Missing Gene BED

Ensure bundled data is present:
```bash
ls $(python -c "from krewlyzer.assets import AssetManager; a=AssetManager('hg19'); print(a.base_path)")/genes/GRCh37/
# Should show: xs1.genes.bed.gz, xs2.genes.bed.gz
```

---

## Next Steps

- [Panel Mode Details](panel-mode.md) - How on/off-target splitting works
- [PON Building](building-pon.md) - Create your own PON
- [JSON Schema](../features/output/json-output.md) - ML integration
- [Troubleshooting](../resources/troubleshooting.md) - Common issues
