# Output Files Reference

Complete reference for every file Krewlyzer produces — what it contains, what each column means, when to use it, and how to apply it in ML models.

---

## Quick Reference

| File | Feature | Resolution | ML Signal |
|------|---------|-----------|-----------|
| [`{s}.FSD.tsv`](#fsd-fragment-size-distribution) | FSD | Per chr arm | Arm-level fragmentation shift |
| [`{s}.FSR.tsv`](#fsr-fragment-size-ratio) | FSR | 5 Mb window | PON-normalized short/long ratio |
| [`{s}.FSC.tsv`](#fsc-fragment-size-coverage) | FSC | 5 Mb window | Multi-channel size coverage |
| [`{s}.FSC.gene.tsv`](#fsc-gene-level) | FSC | Per gene | Gene fragmentation composition |
| [`{s}.FSC.regions.tsv`](#fsc-region-level) | FSC | Per exon/target | Exon fragmentation composition |
| [`{s}.FSC.regions.e1only.tsv`](#fsc-e1-only) | FSC-E1 | Per gene (E1) | Promoter-proximal signal |
| [`{s}.fsc_counts.tsv`](#fsc-counts-pre-correction) | FSC-raw | Per GC bin | GC correction diagnostics |
| [`{s}.WPS.parquet`](#wps-windowed-protection-score) | WPS | Per anchor | Nucleosome positioning |
| [`{s}.WPS.panel.parquet`](#wps-panel) | WPS | Panel anchors | Gene-level nucleosome |
| [`{s}.WPS_background.parquet`](#wps-background) | WPS-bg | Alu stacks | Global chromatin state |
| [`{s}.EndMotif.tsv`](#endmotif) | Motif | Global (1 row) | 256 4-mer end frequencies |
| [`{s}.EndMotif1mer.tsv`](#endmotif-1-mer) | Motif | Global (4 rows) | Base composition at ends |
| [`{s}.BreakPointMotif.tsv`](#breakpointmotif) | Motif | Global (1 row) | 256 4-mer break frequencies |
| [`{s}.MDS.tsv`](#mds-motif-diversity-score) | MDS | Global (1 row) | Motif diversity scalar |
| [`{s}.MDS.exon.tsv`](#mds-exon-level) | Region-MDS | Per exon | Per-exon motif diversity |
| [`{s}.MDS.gene.tsv`](#mds-gene-level) | Region-MDS | Per gene | Per-gene MDS + E1 |
| [`{s}.OCF.tsv`](#ocf-orientation-aware-fragmentation) | OCF | Per tissue | Tissue-of-origin score |
| [`{s}.OCF.sync.tsv`](#ocf-sync) | OCF | Positional | Strand-phased profiles |
| [`{s}.TFBS.tsv`](#tfbs-transcription-factor-binding-site-entropy) | TFBS | Per TF | TF footprint entropy |
| [`{s}.TFBS.sync.tsv`](#tfbs-sync) | TFBS | Per TF × size | Size-resolved TF profiles |
| [`{s}.ATAC.tsv`](#atac-chromatin-accessibility-entropy) | ATAC | Per tissue | ATAC region entropy |
| [`{s}.ATAC.sync.tsv`](#atac-sync) | ATAC | Per tissue × size | Size-resolved ATAC profiles |
| [`{s}.mFSD.tsv`](#mfsd-mutant-fragment-size-distribution) | mFSD | Per variant | Variant fragment size metrics |
| [`{s}.mFSD.distributions.tsv`](#mfsd-distributions) | mFSD | Per variant × size | Raw size histograms |
| [`{s}.UXM.tsv`](#uxm-methylation) | UXM | Per region | U/X/M methylation fractions |
| [`{s}.correction_factors.tsv`](#gc-correction-factors) | GC | Per (len, GC) bin | GC bias weights |
| [`{s}.metadata.json`](#metadata-json) | Meta | Global | Run parameters + QC |
| [`{s}.features.json`](#features-json) | All | All | Unified ML feature export |

> `{s}` = sample name. On-target variants of most files follow `{s}.*.ontarget.tsv`.

---

## Core Fragmentomics

### FSD (Fragment Size Distribution)

**File:** `{sample}.FSD.tsv` / `{sample}.FSD.ontarget.tsv`

FSD measures how many fragments of each size (65–400 bp, 5 bp bins) come from each chromosomal arm. Each row is one arm.

#### Columns

| Column | Type | Description |
|--------|------|-------------|
| `region` | str | Arm label, e.g. `chr1p`, `chr17q` |
| `65-69`, `70-74`, … `395-399` | float | GC-corrected fragment count in that 5 bp size bin |
| `total` | float | Total GC-corrected fragment count for this arm |

#### Purpose & Use Cases

- **ARM-LEVEL fragmentation fingerprint**: Each arm's histogram reflects the chromatin state of that chromosomal region
- **Aneuploidy / CNV detection**: Arms with deletions or amplifications show altered absolute counts in `total`
- **Tumor-specific size shift**: Cancer plasma shows a systematic shift toward shorter fragments genome-wide

#### ML Use Case

```python
import pandas as pd
import numpy as np

df = pd.read_csv("sample.FSD.tsv", sep="\t", index_col="region")
size_cols = [c for c in df.columns if c != "total"]

# Feature vector: proportion of each size class per arm
props = df[size_cols].div(df["total"] + 1e-9, axis=0)

# Reduce to 3-class: short (65-149), mono (150-259), long (260-400)
props["short_frac"] = props[[c for c in size_cols if int(c.split("-")[0]) < 150]].sum(axis=1)
props["mono_frac"]  = props[[c for c in size_cols if 150 <= int(c.split("-")[0]) < 260]].sum(axis=1)
props["long_frac"]  = props[[c for c in size_cols if int(c.split("-")[0]) >= 260]].sum(axis=1)
```

**Best for**: Arm-level `short_frac` as 46-feature input (2 arms × 23 chromosomes) per sample for pan-cancer models.

!!! note "On-target variant"
    `FSD.ontarget.tsv` uses only reads overlapping target capture regions. Useful for panel-specific copy number analysis, but capture-biased — prefer off-target for fragmentation ML features.

---

### FSR (Fragment Size Ratio)

**File:** `{sample}.FSR.tsv` / `{sample}.FSR.ontarget.tsv`

FSR computes the **PON-normalized short-to-long ratio** per 5 Mb window. This is the primary genome-wide tumor fraction biomarker.

#### Columns

| Column | Type | Description |
|--------|------|-------------|
| `region` | str | 5 Mb window, e.g. `chr1:0-5000000` |
| `short_count` | float | GC-corrected count of short frags (65–149 bp) |
| `long_count` | float | GC-corrected count of long frags (261–400 bp) |
| `total_count` | float | Total GC-corrected count |
| `short_norm` | float | `short_count / PON_short_mean` — PON-normalized short |
| `long_norm` | float | `long_count / PON_long_mean` — PON-normalized long |
| `short_long_ratio` | float | `short_norm / long_norm` — primary biomarker |
| `short_long_log2` | float | `log2(short_long_ratio)` — ML-ready signed metric |
| `short_frac` | float | `short_count / total_count` — raw proportion |
| `long_frac` | float | `long_count / total_count` — raw proportion |

#### Purpose & Use Cases

- **Tumor fraction estimation**: `short_long_ratio` increases proportionally with ctDNA fraction
- **Genome-wide cancer screen**: Profile of 500+ windows per sample captures focal and arm-level alterations
- **PON comparison**: PON normalization (`short_norm`, `long_norm`) removes batch effects before ratio — critical for cross-batch comparison

!!! important "Why not just use `short_frac` from FSC?"
    `short_frac` is a raw proportion — it conflates true biology with library prep and GC bias. `short_long_ratio` divides PON-normalized values, canceling technical noise. Use FSR for any cross-sample comparison.

#### ML Use Case

```python
df = pd.read_csv("sample.FSR.tsv", sep="\t")

# Primary feature vector: ~500 windows × 1 scalar
feature_vec = df["short_long_log2"].values  # signed, mean ~0 in healthy

# Genome-wide statistics as compact features
features = {
    "fsr_mean":   df["short_long_log2"].mean(),
    "fsr_std":    df["short_long_log2"].std(),
    "fsr_q90":    df["short_long_log2"].quantile(0.9),
    "fsr_n_elevated": (df["short_long_log2"] > 0.3).sum(),
}
```

**Typical range**: healthy ~0.0 ± 0.15; high ctDNA > +0.4 (more short frags than PON)

---

### FSC (Fragment Size Coverage)

**File:** `{sample}.FSC.tsv` / `{sample}.FSC.ontarget.tsv`

FSC counts fragments in 6 non-overlapping size channels across 5 Mb windows — the foundational multi-channel coverage feature.

#### Columns

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome |
| `start` | int | Window start (0-based) |
| `end` | int | Window end |
| `ultra_short` | float | GC-corrected count, 65–100 bp |
| `core_short` | float | GC-corrected count, 101–149 bp |
| `mono_nucl` | float | GC-corrected count, 150–220 bp |
| `di_nucl` | float | GC-corrected count, 221–260 bp |
| `long` | float | GC-corrected count, 261–400 bp |
| `ultra_long` | float | GC-corrected count, 401–1000 bp |
| `total` | float | GC-corrected total, 65–1000 bp |
| `mean_gc` | float | Mean GC fraction of fragments in window |
| `*_log2` | float | log₂(channel / PON_mean), with `--pon-model` |
| `*_reliability` | float | 1/(PON_variance + k) — weight for PON columns |

#### Purpose & Use Cases

- **Multi-channel fragmentation profile**: Each channel represents fragments sharing a biological origin (nucleosomal, sub-nucleosomal, apoptotic)
- **CNV proxy**: `total` counts per arm reflect coverage depth — useful for detecting large-scale copy number events
- **PON log2 ratios**: `core_short_log2`, `mono_nucl_log2` etc. are analogous to CNV log-ratio tracks

#### ML Use Case

```python
df = pd.read_csv("sample.FSC.tsv", sep="\t")

# 6-channel feature matrix: windows × channels
channels = ["ultra_short", "core_short", "mono_nucl", "di_nucl", "long", "ultra_long"]
X = df[channels].values  # shape: (n_windows, 6)

# Normalize to proportions (remove depth variation)
X_prop = X / (df["total"].values[:, None] + 1e-9)

# With PON: use log2 ratios directly (already depth-normalized)
log2_cols = [c + "_log2" for c in channels if c + "_log2" in df.columns]
X_pon = df[log2_cols].values  # shape: (n_windows, 6) — centered near 0 in healthy
```

**Best for**: Input to CNV callers, genome-wide fragmentation classifiers, and tumor fraction regression.

---

### FSC Gene-Level

**File:** `{sample}.FSC.gene.tsv`  
**Requires:** `--assay xs2` (or other assay code) / `run-all`

Aggregates FSC across all exons for each panel gene. Rows = genes.

#### Columns

| Column | Type | Description |
|--------|------|-------------|
| `gene` | str | HGNC symbol (e.g. `ATM`, `TP53`) |
| `n_regions` | int | Number of exons/targets captured |
| `total_bp` | int | Total base pairs covered |
| `ultra_short` … `long` | float | GC-corrected count per size class |
| `total` | float | Total GC-corrected count |
| `ultra_short_ratio` … `long_ratio` | float | `channel / total` — size composition |
| `normalized_depth` | float | RPKM-like: `(total × 10⁹) / (total_bp × total_frags)` |

#### Purpose & Use Cases

- **Gene-level copy number**: `normalized_depth` enables comparing coverage across genes in the same sample
- **Gene fragmentation composition**: `*_ratio` columns show whether a gene's reads are enriched for short (tumor) or long (normal) fragments
- **Panel-level feature matrix**: 146 genes × 6 channels = 876 features per sample

#### ML Use Case

```python
df = pd.read_csv("sample.FSC.gene.tsv", sep="\t").set_index("gene")

# Gene-level short enrichment
df["tumor_signal"] = df["ultra_short_ratio"] + df["core_short_ratio"]

# Full channel composition feature matrix: 146 genes × 5 ratios
ratio_cols = ["ultra_short_ratio", "core_short_ratio", "mono_nucl_ratio", "di_nucl_ratio", "long_ratio"]
X = df[ratio_cols].values  # shape: (146, 5) per sample

# Normalized depth for CNV
cnv_proxy = df["normalized_depth"]  # pivot across samples → CNV log-ratio
```

**Best for**: Gene-specific models, tissue-of-origin (which genes show altered fragmentation), copy number inference.

---

### FSC Region-Level

**File:** `{sample}.FSC.regions.tsv`  
**Requires:** `--assay` or `run-all`

Per-exon/target fragment size coverage. Most granular FSC output.

#### Columns

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome |
| `start` / `end` | int | Exon/target coordinates |
| `gene` | str | Gene symbol |
| `region_name` | str | Unique exon/target identifier |
| `region_bp` | int | Region size in bp |
| `ultra_short` … `long` | float | GC-corrected counts |
| `total` | float | Total count |
| `ultra_short_ratio` … `long_ratio` | float | `channel / total` |
| `normalized_depth` | float | RPKM-like depth |

#### Purpose & Use Cases

- **Exon-level resolution**: Useful when only specific exons (e.g., hotspot exons in `TP53`) show altered fragmentation
- **Fine-grained CNV**: Detect focal amplifications or deletions at single-exon scale
- **Input for PON building**: Use to compute exon-level expected depth profiles

#### ML Use Case

```python
df = pd.read_csv("sample.FSC.regions.tsv", sep="\t")

# Pivot: regions × channels, one row per region
pivot = df.pivot_table(index="region_name", values=["core_short_ratio", "mono_nucl_ratio"])

# Filter to well-covered regions
well_covered = df[df["total"] > 50]["region_name"]
df_filtered = df[df["region_name"].isin(well_covered)]
```

---

### FSC E1-Only

**File:** `{sample}.FSC.regions.e1only.tsv`  
**Requires:** `--assay` or `run-all` (disable with `--disable-e1-aggregation`)

First exon (E1) per gene only. Same columns as `FSC.regions.tsv`.

#### Purpose & Use Cases

- **Promoter-proximal fragmentation**: E1 = first exon = nucleosome-depleted region (NDR) near TSS. NDRs have the most cancer-specific fragmentation patterns (Helzer et al. 2025)
- **Highest cancer signal**: E1 consistently outperforms whole-gene FSC in early cancer detection tasks
- **Compact feature set**: 146 genes × 1 exon = compact, interpretable vector

#### ML Use Case

```python
df = pd.read_csv("sample.FSC.regions.e1only.tsv", sep="\t").set_index("gene")

# Primary cancer signal: promoter short enrichment
df["promoter_short"] = df["ultra_short_ratio"] + df["core_short_ratio"]

# 146-gene feature vector — best single FSC feature for early detection
X = df["promoter_short"].values
```

!!! tip "Use E1 over gene-level for early detection models"
    `e1only` routinely achieves lower AUC for early-stage cancer vs `FSC.gene.tsv` because E1 captures NDR-specific fragmentation that is washed out by whole-gene averaging.

---

### FSC Counts (Pre-Correction)

**File:** `{sample}.fsc_counts.tsv` / `{sample}.correction_factors.ontarget.tsv`

Raw bin-level fragment counts before GC correction — used internally for GC model training.

#### Purpose & Use Cases

- **GC bias diagnostics**: Compare observed vs expected counts per (length, GC) bin
- **Batch QC**: Libraries with systematic GC bias show large correction factors in specific bins
- **Not for ML features**: GC-corrected values in `FSC.tsv` are the right input; `fsc_counts` is pre-correction

---

## Nucleosome Positioning

### WPS (Windowed Protection Score)

**File:** `{sample}.WPS.parquet`

Per-anchor nucleosome protection profiles. Each row is one genomic anchor (gene TSS or CTCF site). WPS captures how protected (nucleosome-covered) a region is to fragments of two sizes.

#### Columns (Parquet)

| Column | Type | Description |
|--------|------|-------------|
| `region_id` | str | Anchor identifier (e.g. `ENSG00000142611_TSS`) |
| `chrom` | str | Chromosome |
| `center` | int | Anchor midpoint |
| `strand` | str | `+` / `-` |
| `region_type` | str | `TSS`, `CTCF`, etc. |
| `wps_nuc` | float | Raw nucleosomal WPS (120–180 bp fragments) |
| `wps_tf` | float | Raw TF footprint WPS (35–80 bp fragments) |
| `wps_nuc_smooth` | float | Savitzky-Golay smoothed nucleosomal WPS |
| `wps_tf_smooth` | float | Savitzky-Golay smoothed TF WPS |
| `wps_nuc_mean` | float | Mean WPS across anchor window |
| `wps_tf_mean` | float | Mean TF WPS across anchor window |
| `prot_frac_nuc` | float | Fraction of window with WPS > 0 (nucleosome-covered) |
| `prot_frac_tf` | float | Fraction of window with TF WPS > 0 |
| `wps_nuc_z` | float | Z-score vs PON baseline (with `--pon-model`) |
| `wps_tf_z` | float | TF WPS Z-score vs PON |

#### Purpose & Use Cases

- **Nucleosome positioning**: `wps_nuc_smooth` detects nucleosome phasing at TSS/CTCF anchors
- **TF accessibility**: `wps_tf` / `prot_frac_tf` reflects accessible chromatin at TF binding sites
- **Cancer signal**: Tumor DNA shows flattened/disrupted WPS profiles at TSS of active genes
- **NRL (Nucleosome Repeat Length)**: Computed in `WPS_background.parquet` from Alu stacking

#### ML Use Case

```python
import pandas as pd

df = pd.read_parquet("sample.WPS.parquet")

# Feature vector per sample: mean WPS across all anchors
features = {
    "wps_nuc_global_mean": df["wps_nuc_mean"].mean(),
    "wps_nuc_global_std":  df["wps_nuc_mean"].std(),
    "prot_frac_nuc_mean":  df["prot_frac_nuc"].mean(),
    "prot_frac_tf_mean":   df["prot_frac_tf"].mean(),
}

# Per-anchor feature matrix for gene-level models
X = df[["wps_nuc_mean", "wps_tf_mean", "prot_frac_nuc", "prot_frac_tf"]].values
# shape: (~15,000 anchors, 4) — one row per TSS/CTCF
```

**Best for**: Deep learning input (WPS profiles as 1D signals), nucleosome periodicity score, TSS accessibility classifier.

---

### WPS Panel

**File:** `{sample}.WPS.panel.parquet`  
**Requires:** `--assay`

Same columns as `WPS.parquet`, filtered to panel gene anchors (~1,820 for xs2). Rows are panel-gene TSS and CTCF sites.

#### ML Use Case

```python
df = pd.read_parquet("sample.WPS.panel.parquet")

# Compact panel feature: 1820 anchors × 4 = 7280 features
X = df[["wps_nuc_mean", "wps_tf_mean", "prot_frac_nuc", "prot_frac_tf"]].values

# Gene-indexed lookup
df_gene = df.set_index("region_id")
tp53_wps = df_gene.loc["TP53_TSS", "wps_nuc_mean"]
```

---

### WPS Background

**File:** `{sample}.WPS_background.parquet`

Hierarchical Alu element stacking analysis capturing **global chromatin state and nucleosome repeat length (NRL)**.

#### Columns

| Column | Type | Description |
|--------|------|-------------|
| `group_id` | str | `Global_All`, `Family_AluY/S/J/Other`, `Chr{N}_All` |
| `stacked_wps_nuc` | float[] | 30-position binned stacked WPS (nucleosomal) |
| `stacked_wps_tf` | float[] | 30-position binned stacked WPS (TF) |
| `alu_count` | int | Number of Alu elements in this group |
| `mean_wps_nuc` | float | Mean WPS amplitude |
| `nrl_bp` | float | Estimated Nucleosome Repeat Length in bp (~190 in healthy) |
| `nrl_deviation_bp` | float | Deviation from expected 190 bp NRL |
| `periodicity_score` | float | Signal-to-noise ratio of periodicity (0–1) |
| `adjusted_score` | float | Periodicity score penalized by NRL deviation |
| `fragment_ratio` | float | Ratio of short/long fragments at Alu sites |

#### Purpose & Use Cases

- **Global chromatin compaction**: `nrl_bp` shortens in cancer (chromatin opens globally)
- **NRL as tumor biomarker**: Healthy plasma NRL ~190 bp; cancer < 185 bp
- **Background correction**: Used to compute "global_pon" baseline for WPS z-scores

#### ML Use Case

```python
df = pd.read_parquet("sample.WPS_background.parquet")
global_row = df[df["group_id"] == "Global_All"].iloc[0]

features = {
    "nrl_bp": global_row["nrl_bp"],
    "periodicity_score": global_row["periodicity_score"],
    "adjusted_score": global_row["adjusted_score"],
    "fragment_ratio_bg": global_row["fragment_ratio"],
}
```

**Best for**: Global chromatin state features, cancer screening, NRL as continuous tumor fraction predictor.

---

## Motif & Tissue-of-Origin

### EndMotif

**File:** `{sample}.EndMotif.tsv` / `{sample}.EndMotif.ontarget.tsv`

4-mer frequencies at fragment 5′ ends. One row per sample, 256 columns (one per AAAA→TTTT 4-mer).

#### Columns

| Column | Description |
|--------|-------------|
| `AAAA` … `TTTT` | Frequency of that 4-mer at fragment ends (sums to 1.0) |

#### Purpose & Use Cases

- **Tissue-of-origin**: Different tissues have distinct end-motif preferences based on DNASE1L3 activity
- **Cancer detection**: DNASE1L3 is suppressed in cancer, producing a flattened, less-specific motif profile
- **MDS input**: Raw material for Motif Diversity Score calculation

#### ML Use Case

```python
df = pd.read_csv("sample.EndMotif.tsv", sep="\t")
# Single row: 256 4-mer frequencies
X = df.iloc[0].values  # 256-dimensional feature vector

# Reduce by GC content group (64 → 5 groups)
gc_groups = {"AT-rich": [k for k in df.columns if k.count("A") + k.count("T") >= 3], ...}
```

---

### EndMotif 1-mer

**File:** `{sample}.EndMotif1mer.tsv`

Single-base (A/C/G/T) composition at fragment ends.

#### Columns: `base`, `fraction`

#### Purpose & Use Cases

- **GC bias QC**: Should be roughly balanced; extreme GC skew indicates library quality issues
- **DNASE1L3 proxy**: Healthy cfDNA has distinct strand-asymmetric base preferences

---

### BreakPointMotif

**File:** `{sample}.BreakPointMotif.tsv` / `{sample}.BreakPointMotif.ontarget.tsv`

4-mer frequencies at internal fragment **breakpoints** (rather than ends). Same 256-column format as EndMotif.

#### Purpose vs EndMotif

| | EndMotif | BreakPointMotif |
|---|---------|-----------------|
| **Measures** | DNASE1L3 cutting preference | Mechanical fragmentation patterns |
| **Cancer signal** | DNASE1L3 suppression | Chromatin compaction / MNase-like cleavage |
| **Correlation** | r ~ 0.6 with BPM | Complementary signal |

**Best for**: Combining both in multimodal ML models — they capture distinct biological processes.

---

### MDS (Motif Diversity Score)

**File:** `{sample}.MDS.tsv` / `{sample}.MDS.ontarget.tsv`

Single-number summary of end-motif randomness (Shannon entropy of 256 4-mers).

#### Columns: `MDS` (float)

**Range**: Healthy plasma ~0.80–0.85; cancer (DNASE1L3 suppressed) < 0.75

#### ML Use Case

```python
mds = float(pd.read_csv("sample.MDS.tsv", sep="\t")["MDS"].iloc[0])
# Single scalar — highly interpretable cancer feature
```

---

### MDS Exon-Level

**File:** `{sample}.MDS.exon.tsv`  
**Requires:** `region-mds` command or `run-all`

MDS calculated per exon/target from BAM reads overlapping that region.

#### Columns

| Column | Description |
|--------|-------------|
| `gene` | Gene symbol |
| `name` | Exon identifier (`gene:exonN` for WGS; target name for panel) |
| `chrom` | Chromosome |
| `start` / `end` | Exon coordinates |
| `strand` | Strand |
| `n_fragments` | Fragments overlapping this exon |
| `mds` | Motif Diversity Score for this exon |

#### ML Use Case

```python
df = pd.read_csv("sample.MDS.exon.tsv", sep="\t")

# Filter low-coverage exons
df_high = df[df["n_fragments"] >= 20]

# Per-exon MDS as feature matrix (rows = exons)
X = df_high[["mds"]].values  # or pivot into gene × exon matrix
```

---

### MDS Gene-Level

**File:** `{sample}.MDS.gene.tsv`  
**Requires:** `region-mds` command or `run-all`

Gene-level aggregation of per-exon MDS, plus E1 (first exon) MDS as the promoter-proximal signal.

#### Columns

| Column | Description |
|--------|-------------|
| `gene` | Gene symbol |
| `n_exons` | Number of exons with data |
| `n_fragments` | Total fragments across all exons |
| `mds_mean` | Mean MDS across all exons |
| `mds_e1` | MDS of E1 (first exon) only |
| `mds_std` | Standard deviation of per-exon MDS |
| `mds_z` | Z-score vs PON (with `--pon-model`) |
| `mds_e1_z` | E1 MDS z-score vs PON |

#### ML Use Case

```python
df = pd.read_csv("sample.MDS.gene.tsv", sep="\t").set_index("gene")

# E1 MDS is highest-signal feature (promoter-proximal NDR)
X_e1 = df["mds_e1"].values           # 146 features for panel
X_z  = df["mds_e1_z"].fillna(0).values  # PON-normalized — zero-centered in healthy
```

**Best for**: Gene-level cancer classifiers, promoter aberration detection, combining with FSC-E1.

---

### OCF (Orientation-aware cfDNA Fragmentation)

**File:** `{sample}.OCF.tsv` / `{sample}.OCF.ontarget.tsv` / `{sample}.OCF.offtarget.tsv`

Tissue-of-origin scores based on strand asymmetry of fragment ends at tissue-specific open chromatin regions.

#### Columns

| Column | Description |
|--------|-------------|
| `tissue` | Tissue type (Liver, Lung, Colon, Placenta, etc.) |
| `OCF` | Raw score = `(U - D) / (U + D)` — upstream/downstream asymmetry |
| `ocf_z` | OCF z-score vs PON (with `--pon-model`) |

#### Purpose & Use Cases

- **Tissue-of-origin**: Which tissue contributed most cfDNA — useful for cancer site-of-origin
- **Multi-tissue mixture deconvolution**: Multiple elevated `ocf_z` values suggest multi-tissue contribution
- **ctDNA fraction proxy**: Overall OCF magnitude correlates with cfDNA purity

#### ML Use Case

```python
df = pd.read_csv("sample.OCF.tsv", sep="\t").set_index("tissue")

# OCF score vector across tissues (e.g. 10 tissues = 10 features)
X_ocf = df["OCF"].values
X_z   = df["ocf_z"].fillna(0).values  # zero-centered in healthy

# Tissue with max signal
top_tissue = df["ocf_z"].idxmax()
```

**Best for**: Primary tumor site-of-origin classification, multi-class tissue deconvolution.

---

### OCF Sync

**File:** `{sample}.OCF.sync.tsv`

Positional strand-specific protection profiles — detailed positional data underlying the OCF summary score.

#### Columns: `label`, `count`, `mean_size`, `entropy`

**Use**: Raw data for visualizing strand phasing; input to advanced nucleosome positioning models.

---

### TFBS (Transcription Factor Binding Site Entropy)

**File:** `{sample}.TFBS.tsv` / `{sample}.TFBS.ontarget.tsv`

Fragment size entropy at TFBS regions for 808 transcription factors. Reflects chromatin accessibility and TF binding.

#### Columns

| Column | Description |
|--------|-------------|
| `label` | TF name (e.g. `CTCF`, `SP1`, `E2F1`) |
| `count` | Fragment count |
| `mean_size` | Mean fragment size at this TF's sites |
| `entropy` | Shannon entropy of fragment size distribution |

#### Purpose & Use Cases

- **TF accessibility**: Low entropy = dominated by one size class (nucleosomal); high entropy = mixed, accessible
- **TF-specific cancer signal**: Some TFs (E2F family, SP1) show altered fragmentation in cancer
- **808-feature vector**: One entropy value per TF = large, rich feature set

#### ML Use Case

```python
df = pd.read_csv("sample.TFBS.tsv", sep="\t").set_index("label")

# 808 TF entropy values — high-dimensional feature vector
X_tfbs = df["entropy"].values

# Mean size as complementary feature
X_size = df["mean_size"].values

# Combine
X = np.stack([X_tfbs, X_size], axis=1)  # shape: (808, 2)
```

---

### TFBS Sync

**File:** `{sample}.TFBS.sync.tsv`

Per-TF × per-size distribution — the raw size histogram for each TF.

#### Columns: `label`, `size`, `count`, `proportion`

**Use**: Size-resolved TF footprinting; input for detecting nucleosome-footprint transitions.

---

### ATAC (Chromatin Accessibility Entropy)

**File:** `{sample}.ATAC.tsv` / `{sample}.ATAC.ontarget.tsv`

Fragment size entropy at ATAC-seq accessible regions for 23 cancer-relevant tissue types.

#### Columns

Same as TFBS: `label` (tissue type), `count`, `mean_size`, `entropy`

#### Purpose & Use Cases

- **Tissue-specific accessibility**: Which tissue's ATAC peaks show altered fragmentation
- **Cancer type inference**: Different cancer types show tissue-specific ATAC entropy patterns
- **Complementary to OCF**: OCF uses strand asymmetry; ATAC uses size entropy — different signal axis

#### ML Use Case

```python
df = pd.read_csv("sample.ATAC.tsv", sep="\t").set_index("label")

# 23 tissue entropy values — compact tissue-of-origin feature
X_atac = df["entropy"].values

# Combine OCF + ATAC tissue vectors for multimodal tissue classifier
X_combined = np.concatenate([ocf_scores, atac_entropy])
```

---

### ATAC Sync

**File:** `{sample}.ATAC.sync.tsv`

Per-tissue × per-size fragment distributions.

#### Columns: `label`, `size`, `count`, `proportion`

---

## Variant-Level

### mFSD (Mutant Fragment Size Distribution)

**File:** `{sample}.mFSD.tsv`  
**Requires:** `--maf` / `run-all` with MAF input

Per-variant fragment size analysis. Compares ALT-bearing fragments vs REF, NonREF, and N (uncertain) allele classes.

#### Column Groups (46 total)

| Group | Columns | Description |
|-------|---------|-------------|
| **Variant** (5) | `Chrom`, `Pos`, `Ref`, `Alt`, `VarType` | Variant coordinates and type |
| **Raw Counts** (5) | `REF_Count`, `ALT_Count`, `NonREF_Count`, `N_Count`, `Total_Count` | Fragments per allele class |
| **GC-Weighted** (5) | `REF_Weighted`, `ALT_Weighted`, `NonREF_Weighted`, `N_Weighted`, `VAF_GC_Corrected` | GC-corrected counts and VAF |
| **Log-Likelihood** (2) | `ALT_LLR`, `REF_LLR` | For low-N variants (ALT_Count < 5) |
| **Mean Sizes** (4) | `REF_MeanSize`, `ALT_MeanSize`, `NonREF_MeanSize`, `N_MeanSize` | Mean fragment size per class |
| **KS Tests** (18) | `Delta_*/KS_*/KS_Pval_*` × 6 pairings | KS distance + p-value for each allele pair |
| **Derived** (5) | `VAF_Proxy`, `Error_Rate`, `N_Rate`, `Size_Ratio`, `Quality_Score` | Summary biomarkers |
| **Flags** (2) | `ALT_Confidence`, `KS_Valid` | Quality indicators |

#### Purpose & Use Cases

- **MRD detection**: `VAF_GC_Corrected` + `KS_Pval_ALT_REF` identify true somatic mutations vs noise
- **Fragment size as orthogonal evidence**: `Delta_ALT_REF` negative = ALT fragments shorter than REF = tumor-derived ctDNA
- **Duplex support**: `ALT_LLR` provides statistically valid evidence even at 1–2 fragment counts
- **Multi-variant aggregation**: Combine evidence across many variants to estimate ctDNA fraction

#### ML Use Case

```python
df = pd.read_csv("sample.mFSD.tsv", sep="\t")

# Filter to high-quality variants
df_hq = df[(df["ALT_Confidence"] == "HIGH") & (df["KS_Valid"] == True)]

# Per-variant feature vector
features = df_hq[["VAF_GC_Corrected", "Delta_ALT_REF", "KS_ALT_REF",
                   "Size_Ratio", "Quality_Score", "ALT_LLR"]].values

# Sample-level aggregation: weighted average across variants
sample_vaf = (df_hq["VAF_GC_Corrected"] * df_hq["Quality_Score"]).sum() / df_hq["Quality_Score"].sum()
```

**Best for**: MRD detection models, ctDNA fraction estimation from targeted panels, variant-level cancer classifiers.

---

### mFSD Distributions

**File:** `{sample}.mFSD.distributions.tsv`  
**Requires:** `--output-distributions` flag

Per-variant raw size histograms for manual inspection.

#### Columns: `Chrom`, `Pos`, `Ref`, `Alt`, `Category`, `Size`, `Count`

**Use**: Visualizing the ALT vs REF size shift for individual variants; QC and debugging.

---

## Methylation

### UXM (Methylation)

**File:** `{sample}.UXM.tsv`  
**Requires:** Methylation-enabled BAM input

CpG methylation classification per genomic region. Each fragment is classified as Unmethylated (U), Partially-methylated (X), or fully Methylated (M).

#### Columns

| Column | Description |
|--------|-------------|
| `region` | Genomic region identifier |
| `U` | Fraction of unmethylated fragments |
| `X` | Fraction of partially methylated fragments |
| `M` | Fraction of fully methylated fragments |

#### Purpose & Use Cases

- **Tumor DNA methylation**: Cancer shows global hypomethylation (↑U) and focal hypermethylation (↑M at TSGs)
- **Tissue-of-origin**: Tissue-specific methylation patterns enable ctDNA source attribution
- **Multi-modal fusion**: Combine with FSR and MDS for maximum sensitivity in early detection

#### ML Use Case

```python
df = pd.read_csv("sample.UXM.tsv", sep="\t").set_index("region")

# 3-class methylation state per region
X = df[["U", "X", "M"]].values  # shape: (n_regions, 3); rows sum to 1.0

# Sample-level statistics
features = {
    "global_U": df["U"].mean(),   # high = hypomethylation (cancer signal)
    "global_M": df["M"].mean(),   # varies by tissue
    "U_std": df["U"].std(),       # heterogeneity across regions
}
```

---

## Diagnostic / QC Outputs

### GC Correction Factors

**File:** `{sample}.correction_factors.tsv` / `{sample}.correction_factors.ontarget.tsv`

GC-bias correction weights per (fragment length bin, GC content) pair.

#### Columns

| Column | Description |
|--------|-------------|
| `length_bin_min` | Fragment length bin lower bound (bp) |
| `length_bin_max` | Fragment length bin upper bound (bp) |
| `gc_percent` | GC content percentage (0–100) |
| `observed` | Observed fragment count |
| `expected` | Expected count from GC model |
| `correction_factor` | `expected / observed` — multiply raw counts by this |

#### Purpose & Use Cases

- **Library QC**: Values far from 1.0 indicate GC bias in sequencing/capture
- **Batch effect debugging**: Compare correction factors across runs to detect systematic issues
- **Input to GC model PON**: Used to build PON correction factor baselines

!!! note "Not for ML features"
    GC correction is already applied to FSC, FSR, FSD, and WPS outputs. Use those corrected values. `correction_factors.tsv` is diagnostic data.

---

### Metadata JSON

**File:** `{sample}.metadata.json`

Run parameters, QC metrics, and processing provenance.

```json
{
  "sample_id": "sample_001",
  "krewlyzer_version": "0.6.0",
  "genome": "hg19",
  "assay": "xs2",
  "total_fragments": 8234567,
  "on_target_rate": 0.44,
  "mean_fragment_size": 168.3,
  "duplication_rate": 0.12,
  "processing_time_s": 342.1
}
```

**Use**: Filter samples by QC thresholds before ML training (e.g. `total_fragments > 5M`, `on_target_rate > 0.3`).

---

### Features JSON

**File:** `{sample}.features.json`  
**Requires:** `--generate-json`

Unified export of all features above in a single JSON file for ML pipelines. See [JSON Output Reference](../features/output/json-output.md) for the complete schema.

```python
import json
with open("sample.features.json") as f:
    features = json.load(f)

# Access any output without reading individual TSVs
fsr_windows = features["fsr"]["off_target"]
mds = features["motif"]["mds"]
region_mds = features["region_mds"]["gene"]
```

---

## Feature Selection Guide for ML Models

| Model type | Recommended features | File(s) |
|-----------|---------------------|---------|
| Pan-cancer screen (WGS) | FSR `short_long_log2`, MDS, WPS `prot_frac_nuc_mean`, FSD `short_frac` | FSR, MDS, WPS_background, FSD |
| Pan-cancer screen (panel) | FSC-E1 `promoter_short`, MDS-E1, OCF `ocf_z`, WPS panel | FSC.e1only, MDS.gene, OCF, WPS.panel |
| Tumor fraction regression | FSR `short_long_log2` genome-wide, WPS `nrl_bp` | FSR, WPS_background |
| Cancer type / site-of-origin | OCF tissue vector, ATAC tissue vector, TFBS vector | OCF, ATAC, TFBS |
| Gene-level amplification | FSC gene `normalized_depth`, FSD arm `total` | FSC.gene, FSD |
| MRD (residual disease) | mFSD `VAF_GC_Corrected`, `Quality_Score`, `Delta_ALT_REF` | mFSD |
| Promoter accessibility | WPS E1 `wps_nuc_mean`, MDS E1, FSC-E1 ratios | WPS.panel, MDS.gene, FSC.e1only |
| Methylation-augmented | UXM `U`/`M` + FSR + MDS | UXM, FSR, MDS |

---

## See Also

- [JSON Output Reference](../features/output/json-output.md) — unified JSON schema for all features
- [FSR vs FSC ratios](../features/core/fsc.md#fsc-ratios-vs-fsr--not-redundant) — why they differ
- [PON Models](../reference/pon-models.md) — normalization baselines
- [GC Correction](../guides/gc-correction.md) — correction methodology
