# Panel of Normals (PON)

The Panel of Normals (PON) is a unified model built from healthy plasma samples that enables:

1. **GC bias correction** - Per-fragment correction for GC content bias
2. **Z-score normalization** - Detect deviations from healthy baseline for all features
3. **Panel mode support** - Dual on/off-target baselines for capture panels

## Quick Start

```bash
# Build PON from healthy samples
krewlyzer build-pon samples.txt --assay msk-access-v2 -r hg19.fa -o pon.parquet

# Use PON for sample processing
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ -P pon.parquet
```

## Auto-PON Loading

When you specify an assay with `-A`, krewlyzer automatically loads the bundled PON:

```bash
# Auto-loads bundled PON for xs2 assay
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ -A xs2 -G hg19
```

This is equivalent to explicitly passing `-P` with the bundled PON path.

## Skipping Z-Score Normalization (`--skip-pon`)

For **ML training workflows** where PON samples are used as true negatives, use `--skip-pon` to output raw features without z-score normalization:

```bash
# Process PON samples as ML negatives (no z-scores)
krewlyzer run-all -i pon_sample.bam -r hg19.fa -o out/ -A xs2 --skip-pon
```

!!! warning
    `-P` and `--skip-pon` are **mutually exclusive**. If you specify an explicit PON model, you want z-scores applied. Use `--skip-pon` only with `-A` (assay) for the ML negatives workflow.

The `--skip-pon` flag:
- Works with `-A/--assay` (auto-loads bundled PON but skips z-scores)
- Available on all tools: `run-all`, `fsc`, `fsd`, `fsr`, `wps`, `ocf`, `region-entropy`, `motif`
- Logs which tools are skipping normalization

## PON Variant Selection (`--pon-variant`)

For duplex sequencing workflows (fgbio/Marianas), use `--pon-variant duplex` to select PONs built from duplex consensus reads:

```bash
# Default: all_unique PON (maximum coverage for standard cfDNA)
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ -A xs2

# Duplex PON (highest accuracy for duplex sequencing data)
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ -A xs2 --pon-variant duplex
```

| Variant | Description | Best For |
|---------|-------------|----------|
| `all_unique` | Built from all unique reads | Standard cfDNA (default) |
| `duplex` | Built from duplex consensus reads | Duplex sequencing workflows |

!!! tip
    The `--pon-variant` flag is independent of the `--duplex` flag for mFSD. Use `--duplex` for mFSD weighting (enables cD tag usage), and `--pon-variant` for PON selection across all tools.

The `--pon-variant` flag:
- Defaults to `all_unique` (maximum coverage PON)
- Available on all PON-using tools: `run-all`, `fsc`, `fsd`, `fsr`, `wps`, `ocf`, `motif`, `region-entropy`, `region-mds`
- File structure: `pon/{genome}/{variant}/{assay}.{variant}.pon.parquet`

## PON Components

A PON is a single long-format Parquet with a `table` column; each block below is
one value of it. All 24 are present in the shipped GRCh37 models.

| Block | Description | Used by |
|---|---|---|
| `metadata` | Assay, cohort size, build date, and the provenance fields below | every reader |
| `gc_bias` | Expected coverage by GC content per fragment type | FSC, FSR, WPS |
| `fsd_baseline` | Size distribution per chromosome arm | FSD |
| `wps_baseline` | 200-position mean and σ profiles per anchor | WPS |
| `wps_baseline_panel` | As `wps_baseline`, over the assay's panel anchors | WPS panel |
| `wps_shape_baseline` | Per-anchor baselines for the derived shape statistics | WPS |
| `wps_background` | NRL and periodicity per Alu group | WPS background |
| `ocf_baseline` | Open-chromatin scores per tissue | OCF |
| `mds_baseline` | End-motif k-mer frequencies and motif diversity | EndMotif, MDS |
| `breakpoint_motif_baseline` | **Breakpoint** k-mer frequencies | BreakPointMotif |
| `tfbs_baseline` | Per-TF entropy mean/σ | Region entropy |
| `atac_baseline` | Per-cancer-type entropy mean/σ | Region entropy |
| `region_mds` | Per-gene MDS mean/σ | Region MDS |
| `region_mds_exon` | Per-exon MDS mean/σ | Region MDS |
| `fsc_gene_baseline` | Per-gene normalised depth mean/σ | FSC gene |
| `fsc_region_baseline` | Per-exon normalised depth mean/σ | FSC region |

!!! note "Breakpoint motifs are not end motifs"
    An end motif is the 4-mer at the fragment's 5′ terminus — what the nuclease
    left. A breakpoint motif spans the cut site and includes reference bases
    *not present in the fragment*. The two frequency vectors correlate only
    **0.696** on the same sample, so they need separate baselines. Before 0.9.0
    `BreakPointMotif` was scored against `mds_baseline` and its `frequency_z`
    measured the offset between two definitions rather than any departure from
    healthy.

### What a PON records about itself

`metadata` carries provenance, so a model can be identified after the fact:

| Field | Meaning |
|---|---|
| `krewlyzer_version` | The release the model **ships with**, set by `stamp-pon` — not the version of the code that happened to build it |
| `cohort_digest` | A salted, non-reversible fingerprint of the donor set. Two builds from the same cohort match; it reveals nothing about who is in it |
| `cohort_label` | Optional free-text name |
| `input_kind` | `bam`, `bed`, `mixed` or `outputs` — what the cohort was made of, so the gate can tell a block that was never asked for from one that failed |
| `n_samples`, `build_date`, `assay`, `panel_mode` | As named |

### Version compatibility

**A PON older than 0.9.0 is refused.** That release changed what several
features *mean* — a fabricated `wps_background`, floored σ, a region-MDS fitted
over a different fragment range — so an older model is not merely stale, it is
measuring something else, and every z-score against it is wrong in a way no
schema check can see.

The floor is a compatibility floor, not the package version: it stays at 0.9.0
in later releases, because what changed happened once.

```bash
# Deliberately score against an older model anyway
KREWLYZER_ALLOW_OLD_PON=1 krewlyzer run-all ...
```

The escape hatch exists on purpose. Without a documented way out, someone who
genuinely wants an old model for comparison edits the Parquet instead, and then
the version stops meaning anything.

### σ that is not a measurement

Where a baseline could not measure a spread — every donor identical, or no donor
covering that position — the model stores **NaN**, and the scorers emit no
z-score rather than dividing by a placeholder.

This is enforced at both ends. Values below `1e-9` are treated as floating-point
residue rather than a spread: WPS is a difference of GC-weighted counts, so a
true zero computes as ~1e-17 instead of `0.0`, and dividing by that produced
z-scores up to 6.1 × 10¹⁸ before 0.9.0.

## Panel Mode

For capture panels (like MSK-ACCESS), use `--target-regions` when building the PON:

```bash
krewlyzer build-pon samples.txt -a msk-access-v2 -r hg19.fa -T targets.bed -o pon.parquet
```

This enables:

- **GC model trained on off-target only** — Avoids capture bias
- **Separate on/off-target baselines** — FSD, MDS, OCF each get on-target variants
- **Panel mode detection** — Sample processing auto-detects matching PON mode

Panel mode generates these additional baselines:

| On-Target Baseline | Description |
|---|---|
| `gc_bias_ontarget` | GC correction from on-target reads |
| `fsd_baseline_ontarget` | Size distribution from on-target fragments |
| `mds_baseline_ontarget` | End-motif diversity from on-target fragments |
| `breakpoint_motif_baseline_ontarget` | Breakpoint motifs from on-target fragments |
| `wps_baseline_panel` | WPS profiles over the assay's panel anchors |
| `ocf_baseline_ontarget` | OCF scores from on-target reads |
| `ocf_baseline_offtarget` | OCF scores from off-target reads |
| `tfbs_baseline_ontarget` | TFBS entropy from on-target reads |
| `atac_baseline_ontarget` | ATAC entropy from on-target reads |

## Building a PON

See [build-pon CLI](../guides/building-pon.md) for detailed options.

**Requirements:**
- Minimum 10 healthy samples recommended
- Same assay/panel as samples to be processed
- Same reference genome

### Rebuilding without re-reading the BAMs (`--from-outputs`)

A from-scratch build re-extracts every fragment and takes roughly 15 hours for
a 47-donor cohort. Most rebuilds do not need that: **aggregation** is what
changed, not the per-sample features. `--from-outputs` re-aggregates a
directory of existing `run-all` outputs in minutes per model.

```bash
krewlyzer build-pon --from-outputs /path/to/runall_dirs \
    --assay xs1 --genome hg19 -o xs1.all_unique.pon.parquet
krewlyzer validate-pon xs1.all_unique.pon.parquet
```

This is how all four bundled models were rebuilt for 0.9.0 when the σ floor and
the BreakPointMotif baseline changed — both aggregation-stage fixes.

Use it when the **builder** changed. Use a full build when anything upstream of
it did: the fragment filters, the GC model, or the feature code itself. The
per-sample directories must still exist and must contain every table the model
needs, so confirm before relying on it — a cohort cleaned up after its last
build forces the full path.

> [!NOTE]
> The cohort digest is computed from the sample list, so re-aggregating the same
> cohort must leave it **unchanged**. A digest that moves means the inputs moved,
> not the aggregation — investigate before shipping the model.

## Using PON in Processing

When `--pon-model` is provided to `run-all`:

1. PON is loaded once and passed to all processors
2. Each feature computes z-scores against healthy baseline
3. Output includes both raw values and PON-normalized columns

## Output Columns

With PON, additional columns are added to outputs:

| Feature | PON Column(s) | Description |
|---------|--------------|-------------|
| FSC | `*_log2` | Log2 ratio vs PON expected |
| FSC Gene | `depth_zscore` | Gene-level depth z-score |
| FSC Region | `depth_zscore` | Exon-level depth z-score |
| FSD | `ratio_log2` | Size distribution log ratio |
| WPS | `wps_zscore` | Z-score vs region baseline |
| OCF | `ocf_zscore` | Z-score vs OCF baseline |
| Motif | `mds_z` | Z-score for MDS |
| TFBS | `entropy_z` | Z-score per TF |
| ATAC | `entropy_z` | Z-score per cancer type |
| Region MDS | `mds_z`, `mds_e1_z` | Gene-level and E1 z-scores |

## API Reference

```python
from krewlyzer.pon.model import PonModel

# Load existing PON
pon = PonModel.load("path/to/pon.parquet")

# Access components
gc_expected = pon.get_mean("short")  # Expected at median GC
variance = pon.get_variance("short")  # For reliability scoring

# Check panel mode
if pon.panel_mode:
    print(f"Built with: {pon.target_regions_file}")
```

---

## PON Baselines in Detail

### GC Bias Model (`gc_bias`)

Stores expected fragment coverage per GC content (0-100%) for each fragment type:

| Fragment Type | Size Range | Purpose |
|---------------|------------|---------|
| `short` | 65-149bp | Short fragment correction |
| `intermediate` | 150-259bp | Mono-nucleosomal |
| `long` | 260-400bp | Di-nucleosomal |
| `wps_long` | 120-180bp | WPS nucleosomal |
| `wps_short` | 35-80bp | WPS TF footprint |

### FSD Baseline (`fsd_baseline` / `fsd_baseline_ontarget`)

Size distribution per chromosome arm (46 arms):
- `expected`: Mean proportion at each size bin
- `std`: Standard deviation across PON samples

In panel mode, `fsd_baseline_ontarget` provides a separate baseline for on-target fragment distributions. This is used automatically when processing `{s}.FSD.ontarget.tsv`.

### WPS Baseline (`wps_baseline`)

Per-region nucleosome positioning metrics.

**Schema v1.0 (Scalar):**
- `wps_long_mean/std`: Single nucleosomal WPS value per region
- `wps_short_mean/std`: Single TF footprint value per region

**Schema v2.0 (Vector):**
- `wps_nuc_mean/std`: 200-element vector (nucleosomal footprint)
- `wps_tf_mean/std`: 200-element vector (TF footprint)

!!! tip
    v2.0 enables position-specific z-scores and **Shape Correlation Score** for cancer detection.

#### Shape Score Interpretation

| Score | Interpretation |
|-------|---------------|
| 0.9-1.0 | Healthy nucleosome positioning |
| 0.5-0.9 | Mild chromatin disorganization |
| <0.5 | Significant disruption (cancer signal) |

See [WPS Features](../features/core/wps.md) for output column details.

### OCF Baseline (`ocf_baseline` / `ocf_baseline_ontarget` / `ocf_baseline_offtarget`)

Per-region open chromatin footprint:
- `ocf_mean/std`: OCF score baseline
- `sync_mean/std`: Synchronization score baseline

In panel mode, separate `ocf_baseline_ontarget` and `ocf_baseline_offtarget` baselines are built, enabling z-score normalization for `{s}.OCF.ontarget.tsv` and `{s}.OCF.offtarget.tsv` independently.

### MDS Baseline (`mds_baseline` / `mds_baseline_ontarget`)

Motif diversity expectations:
- `kmer_expected`: 256 4-mer frequencies from healthy samples
- `kmer_std`: Variability per k-mer
- `mds_mean/std`: Expected Motif Diversity Score

In panel mode, `mds_baseline_ontarget` provides a separate baseline computed from on-target fragments only. The on-target MDS z-score is written to `{s}.MDS.ontarget.tsv`.

### TFBS Baseline (`tfbs_baseline`)

Per-TF size entropy:
- `label_stats`: Mean/std entropy per TF (808 transcription factors)
- Enables z-score per TF for detailed regulatory analysis

### ATAC Baseline (`atac_baseline`)

Per-cancer-type size entropy:
- `label_stats`: Mean/std entropy per cancer type (23 types)
- Enables tissue-of-origin scoring

### Region MDS Baseline (`region_mds`)

Per-gene MDS expectations:
- `gene_baseline`: Dict of gene → {mds_mean, mds_std, mds_e1_mean, mds_e1_std}
- Enables gene-level anomaly detection
- E1 (first exon) tracked separately for promoter-proximal sensitivity

### FSC Gene Baseline (`fsc_gene_baseline`)

Per-gene normalized depth baseline (panel mode only):
- `data`: Dict of gene → (mean_depth, std_depth, n_samples)
- Requires minimum 3 samples for reliable statistics
- Clinical use: z-score >> 0 = amplification, z-score << 0 = deletion

### FSC Region Baseline (`fsc_region_baseline`)

Per-exon/probe normalized depth baseline (panel mode only):
- `data`: Dict of region_id → (mean_depth, std_depth, n_samples)
- Region IDs formatted as "chrom:start-end"
- Covers all exons (no filtering by variance)
- Enables detection of focal copy number changes affecting single exons

---

## Interpreting Z-Scores

Z-scores measure how many standard deviations a sample differs from the healthy PON baseline:

$$
z = \frac{x_{\text{sample}} - \mu_{\text{PON}}}{\sigma_{\text{PON}}}
$$

### Clinical Interpretation

| Z-Score Range | Interpretation | Action |
|---------------|----------------|--------|
| **-2 to +2** | Normal range | Within healthy variation |
| **|z| = 2-3** | Mild deviation | Monitor, may be noise |
| **|z| > 3** | Significant | Investigate for ctDNA |
| **|z| > 5** | Extreme | High tumor burden likely |

### Per-Feature Z-Score Meaning

| Feature | Z-Score Column | Positive Z Means | Negative Z Means |
|---------|----------------|------------------|------------------|
| **FSC** | `z_core_short` | More short fragments | Fewer short fragments |
| **FSD** | - | Shifted size distribution | - |
| **WPS** | `wps_nuc_z` | Stronger nucleosome signal | Disrupted nucleosomes |
| **OCF** | `ocf_z` | More open chromatin | Less accessible |
| **MDS** | `mds_z` | More diverse motifs | Less diverse |
| **TFBS** | `entropy_z` | Higher entropy (diverse sizes) | Lower entropy (restricted) |
| **ATAC** | `entropy_z` | Higher entropy | Lower entropy |
| **Region MDS** | `mds_z`, `mds_e1_z` | More diverse at gene | Restricted motifs (aberrant) |

### ML Feature Usage

```python
# Extract z-score features for classification
features = {
    "fsc_short_z": sample_fsc["z_core_short"].mean(),
    "wps_nuc_z": sample_wps["wps_nuc_z"].mean(),
    "mds_z": sample_motif["mds_z"],
}

# Higher |z| = more likely to be tumor
combined_signal = sum(abs(z) for z in features.values())
```

!!! tip
    **Combine z-scores across features** - Single extreme values may be noise, but consistent deviations across FSC, WPS, and MDS are highly indicative of ctDNA.

