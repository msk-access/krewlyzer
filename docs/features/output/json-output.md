# JSON Feature Output

The `--generate-json` flag produces a single JSON file containing all features for ML integration.

```bash
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ \
    --assay xs2 \
    --generate-json
```

This generates `{sample}.features.json` alongside the standard TSV/Parquet outputs.

---

## Top-Level Structure

```json
{
  "schema_version": "1.0",
  "sample_id": "sample_001",
  "krewlyzer_version": "0.6.0",
  "timestamp": "2026-02-28T19:00:00",
  "metadata": {
    "genome": "hg19",
    "assay": "xs2",
    "panel_mode": true,
    "on_target_rate": 0.45
  },
  "features": {
    "fsc": { ... },
    "fsc_gene": [ ... ],
    "fsc_region": [ ... ],
    "fsc_region_e1": [ ... ],
    "fsc_counts": [ ... ],
    "fsr": { ... },
    "fsd": { ... },
    "wps": { ... },
    "wps_panel": { ... },
    "wps_background": { ... },
    "motif": { ... },
    "ocf": { ... },
    "tfbs": { ... },
    "atac": { ... },
    "gc_factors": { ... },
    "mfsd": { ... },
    "region_mds": { ... },
    "uxm": { ... }
  },
  "qc": { ... }
}
```

---

## Feature Schemas

### FSC (Fragment Size Coverage)

Window-based fragment size counts. Split into `off_target` (genome-wide) and `on_target` (panel capture regions).

```json
"fsc": {
  "off_target": [
    {
      "region": "chr1:0-500000",
      "ultra_short": 123,
      "core_short": 4567,
      "mono_nucl": 8901,
      "di_nucl": 2345,
      "long": 678,
      "total": 16614,
      "log_ratio_core_short": -0.15,
      "zscore_core_short": -1.23
    }
  ],
  "on_target": [ ... ]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `region` | string | Genomic window coordinates |
| `ultra_short` | int | Fragments 65–99 bp |
| `core_short` | int | Fragments 100–149 bp |
| `mono_nucl` | int | Fragments 150–259 bp |
| `di_nucl` | int | Fragments 260–399 bp |
| `long` | int | Fragments 400+ bp |
| `log_ratio_*` | float | Log₂(observed/expected) vs PON |
| `zscore_*` | float | Z-score vs PON (if PON provided) |

---

### FSC Gene (Panel Mode Only)

Gene-level fragment size aggregation. Only present with `--assay`.

```json
"fsc_gene": [
  {
    "gene": "ATM",
    "n_regions": 62,
    "total_bp": 8432,
    "ultra_short": 1234,
    "core_short": 5678,
    "mono_nucl": 9012,
    "core_short_ratio": 0.282,
    "normalized_depth": 1245.67,
    "z_core_short": -0.45
  }
]
```

---

### FSC Region (Panel Mode Only)

Per-exon/target fragment size data. More granular than gene-level.

```json
"fsc_region": [
  {
    "chrom": "1",
    "start": 11168235,
    "end": 11168345,
    "gene": "MTOR",
    "region_name": "MTOR_target_02",
    "region_bp": 110,
    "ultra_short": 8.0,
    "core_short": 229.0,
    "mono_nucl": 804.0,
    "di_nucl": 88.0,
    "total": 1129.0,
    "normalized_depth": 1272.71
  }
]
```

---

### FSC Region E1 (Panel Mode Only)

First exon (E1) per gene. E1 = promoter-proximal region with stronger cancer signal (Helzer et al. 2025).

```json
"fsc_region_e1": [
  {
    "chrom": "14",
    "start": 105238685,
    "end": 105238805,
    "gene": "AKT1",
    "region_name": "exon_AKT1_15a_1",
    "region_bp": 120,
    "ultra_short": 19.0,
    "mono_nucl": 635.0,
    "total": 1082.0,
    "normalized_depth": 3039.77
  }
]
```

!!! tip
    Use `fsc_region_e1` for **early cancer detection** models where promoter fragmentation changes are a primary signal.

---

### FSC Counts (Raw GC-Binned Counts)

Raw per-GC-bin, per-size-class fragment counts. Source: `{sample}.fsc_counts.tsv`.

```json
"fsc_counts": [
  {
    "gc_bin": 0.40,
    "len_bin": 150,
    "count": 1234,
    "expected": 1012.5,
    "correction_factor": 1.218
  }
]
```

!!! note
    `fsc_counts` contains pre-correction data. The GC bias correction is already applied to `fsc`, `fsr`, and `fsd` values.

---

### FSR (Fragment Size Ratios)

Biomarker ratios for tumor detection. Split into `off_target` and `on_target`.

```json
"fsr": {
  "off_target": [
    {
      "region": "chr1:0-500000",
      "ultra_short_ratio": 0.0074,
      "core_short_ratio": 0.275,
      "mono_nucl_ratio": 0.536,
      "di_nucl_ratio": 0.141,
      "long_ratio": 0.041,
      "core_short_long_ratio": 6.73
    }
  ],
  "on_target": [ ... ]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `*_ratio` | float | Fraction of total fragments in that size class |
| `core_short_long_ratio` | float | Primary cancer biomarker (higher = more tumor content) |

---

### FSD (Fragment Size Distribution)

Per-arm size distribution profiles. Split into `off_target` and `on_target`.

```json
"fsd": {
  "off_target": {
    "arms": ["1p", "1q", "2p", "..."],
    "size_bins": ["65-69", "70-74", "...", "395-399"],
    "counts": [[123, 456, ...], [...]],
    "total": [12345, 13456, ...]
  },
  "on_target": { ... }
}
```

| Field | Type | Description |
|-------|------|-------------|
| `arms` | string[] | Chromosomal arm labels |
| `size_bins` | string[] | Fragment size range labels (bp) |
| `counts` | int[][] | Count matrix: arms × size bins |
| `total` | int[] | Total fragment count per arm |

---

### WPS (Windowed Protection Score)

Nucleosome positioning profiles. Stored as **columnar arrays** (one value per region), not row records.

```json
"wps": {
  "regions": ["ENSG00000142611_TSS", "ENSG00000157764_TSS", "..."],
  "chrom": ["1", "7", "..."],
  "center": [11166302, 140453136, "..."],
  "wps_nuc": [24.5, 18.2, "..."],
  "wps_tf": [-3.2, -1.8, "..."],
  "wps_nuc_smooth": [23.1, 17.9, "..."],
  "wps_tf_smooth": [-3.0, -1.7, "..."],
  "wps_nuc_mean": [22.8, 17.5, "..."],
  "wps_tf_mean": [-2.9, -1.6, "..."],
  "wps_nuc_z": [1.2, 0.8, "..."],
  "wps_tf_z": [-0.4, -0.2, "..."],
  "prot_frac_nuc": [0.62, 0.55, "..."],
  "prot_frac_tf": [0.41, 0.38, "..."]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `regions` | string[] | Region IDs (one entry per anchor) |
| `wps_nuc` | float[] | Nucleosomal WPS (120–180 bp fragments) |
| `wps_tf` | float[] | TF footprint WPS (35–80 bp fragments) |
| `wps_*_smooth` | float[] | Savitzky-Golay smoothed profiles |
| `wps_*_z` | float[] | Z-scores vs PON baseline |
| `prot_frac_*` | float[] | Protected fraction (values > 0) |

!!! tip "Reconstructing a DataFrame"
    ```python
    import pandas as pd
    wps = features["wps"]
    df = pd.DataFrame({"region_id": wps["regions"], "wps_nuc": wps["wps_nuc"], ...})
    ```

---

### WPS Panel (Panel Mode Only)

Same schema as standard TSV/Parquet WPS output, but stored as row records and filtered to panel gene anchors.

```json
"wps_panel": {
  "n_anchors": 1820,
  "data": [
    {
      "region_id": "ATM_TSS",
      "chrom": "11",
      "center": 108093559,
      "strand": "+",
      "wps_nuc": 22.1,
      "wps_tf": -2.8,
      "prot_frac_nuc": 0.59,
      "prot_frac_tf": 0.38,
      "wps_nuc_z": 0.9,
      "wps_tf_z": -0.3
    }
  ]
}
```

---

### WPS Background

Alu element hierarchical stacking profiles (global, family, per-chromosome groups).

```json
"wps_background": {
  "n_elements": 27,
  "data": [
    {
      "group_id": "Global_All",
      "stacked_wps_nuc": [0.12, 0.08, "..."],
      "stacked_wps_tf": [-0.03, -0.01, "..."],
      "alu_count": 142567,
      "mean_wps_nuc": 0.11,
      "mean_wps_tf": -0.02,
      "nrl_bp": 192.4,
      "nrl_deviation_bp": 2.4,
      "periodicity_score": 0.78,
      "adjusted_score": 0.69,
      "fragment_ratio": 0.31
    }
  ]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `group_id` | string | `Global_All`, `Family_AluY/S/J/Other`, `Chr{N}_All` |
| `stacked_wps_nuc` | float[] | 30-bin binned WPS nucleosomal profile |
| `nrl_bp` | float | Nucleosome Repeat Length in bp (expected ~190 bp) |
| `periodicity_score` | float | SNR-based quality 0–1 |
| `adjusted_score` | float | Periodicity score penalized by NRL deviation |

---

### Motif

End motif (EDM) and breakpoint motif (BPM) 4-mer frequencies, plus MDS scores for off-target and on-target.

```json
"motif": {
  "edm": { "AAAA": 0.0042, "AAAC": 0.0039, "...": "..." },
  "bpm": { "AAAA": 0.0038, "AAAC": 0.0041, "...": "..." },
  "edm_1mer": { "A": 0.197, "C": 0.345, "G": 0.321, "T": 0.138 },
  "mds": 0.8234,
  "mds_z": -1.23,
  "edm_on_target": { "AAAA": 0.0051, "...": "..." },
  "bpm_on_target": { "AAAA": 0.0043, "...": "..." },
  "mds_on_target": 0.7918,
  "mds_z_on_target": -0.91
}
```

| Field | Type | Description |
|-------|------|-------------|
| `edm` | dict | 256 4-mer frequencies at fragment **ends** (off-target) |
| `bpm` | dict | 256 4-mer frequencies at **breakpoints** (off-target) |
| `edm_1mer` | dict | Single-base composition at fragment ends (A/C/G/T) |
| `mds` | float | Motif Diversity Score (off-target) |
| `mds_z` | float | MDS z-score vs PON (only with `--pon-model`) |
| `edm_on_target` | dict | On-target EDM frequencies (panel mode only) |
| `bpm_on_target` | dict | On-target BPM frequencies (panel mode only) |
| `mds_on_target` | float | On-target MDS (panel mode only) |
| `mds_z_on_target` | float | On-target MDS z-score vs PON (panel + `--pon-model`) |

---

### OCF (Orientation-aware cfDNA Fragmentation)

Open chromatin footprint scores by tissue type. Includes off-target, on-target, and positional sync profiles.

```json
"ocf": {
  "off_target": [ { "tissue": "Liver", "score": 0.42, "n_fragments": 12345 } ],
  "on_target":  [ { "tissue": "Liver", "score": 0.51, "n_fragments": 3421 } ],
  "offtarget":  [ { "tissue": "Liver", "score": 0.39, "..." } ],
  "sync":           [ { "pos": -150, "strand_ratio": 0.61, "..." } ],
  "sync_offtarget": [ { ... } ],
  "sync_ontarget":  [ { ... } ]
}
```

| Sub-key | When present | Description |
|---------|-------------|-------------|
| `off_target` | Always | Genome-wide OCF scores |
| `on_target` | Panel mode | On-target capture region OCF |
| `offtarget` | Panel mode | Panel-specific off-target scores |
| `sync` | Always | Positional strand-specific profiles |
| `sync_offtarget` | Panel mode | Sync profiles for off-target |
| `sync_ontarget` | Panel mode | Sync profiles for on-target |

---

### TFBS (Transcription Factor Binding Site Entropy)

Fragment size entropy at TFBS regions. Includes sync (per-TF × per-size) profiles.

```json
"tfbs": {
  "off_target": [
    { "region": "CTCF_chr1_12345", "entropy": 3.45, "n_fragments": 234, "mean_size": 167.5 }
  ],
  "on_target": [ ... ],
  "sync": [ { "tf": "CTCF", "size_bin": 150, "count": 234, "fraction": 0.052 } ],
  "sync_ontarget": [ ... ]
}
```

---

### ATAC (Chromatin Accessibility Regions)

Fragment size entropy at ATAC-seq accessible regions. Includes sync (per-tissue × per-size) profiles.

```json
"atac": {
  "off_target": [
    { "region": "peak_chr1_23456", "entropy": 3.21, "n_fragments": 189, "mean_size": 145.2 }
  ],
  "on_target": [ ... ],
  "sync": [ { "tissue": "BRCA", "size_bin": 150, "count": 189, "fraction": 0.041 } ],
  "sync_ontarget": [ ... ]
}
```

---

### GC Factors (Diagnostic)

GC bias correction factors used internally during processing.

!!! note
    **Not recommended for ML features.** The GC correction is already applied to FSC/FSR/FSD. Use those corrected values instead. GC factors are useful for QC, batch effect detection, and panel development.

```json
"gc_factors": {
  "off_target": [
    { "len_bin": 100, "gc_pct": 45, "correction_factor": 1.12 }
  ],
  "on_target": [ ... ]
}
```

---

### mFSD (Mutant Fragment Size Distribution)

Per-variant fragment size distribution metrics. Only present when a MAF file is provided.

```json
"mfsd": {
  "enabled": true,
  "n_variants": 47,
  "variants": [
    {
      "Chrom": "17", "Pos": 7577548, "Ref": "C", "Alt": "T", "VarType": "SNP",
      "REF_Count": 234, "ALT_Count": 12, "NonREF_Count": 8, "N_Count": 45, "Total_Count": 299,
      "REF_Weighted": 221.4, "ALT_Weighted": 11.3, "NonREF_Weighted": 7.6, "N_Weighted": 42.8, "VAF_GC_Corrected": 0.048,
      "ALT_LLR": 3.41, "REF_LLR": -3.41,
      "REF_MeanSize": 168.3, "ALT_MeanSize": 142.7, "NonREF_MeanSize": 155.1, "N_MeanSize": 171.2,
      "Delta_ALT_REF": -25.6, "KS_ALT_REF": 0.31, "KS_Pval_ALT_REF": 0.003,
      "Delta_ALT_NonREF": -12.4, "KS_ALT_NonREF": 0.18, "KS_Pval_ALT_NonREF": 0.041,
      "Delta_REF_NonREF": 13.2, "KS_REF_NonREF": 0.15, "KS_Pval_REF_NonREF": 0.089,
      "Delta_ALT_N": -28.5, "KS_ALT_N": 0.34, "KS_Pval_ALT_N": 0.001,
      "Delta_REF_N": -3.1, "KS_REF_N": 0.08, "KS_Pval_REF_N": 0.621,
      "Delta_NonREF_N": -16.1, "KS_NonREF_N": 0.19, "KS_Pval_NonREF_N": 0.033,
      "VAF_Proxy": 0.049, "Error_Rate": 0.027, "N_Rate": 0.150, "Size_Ratio": 0.848, "Quality_Score": 0.731,
      "ALT_Confidence": "HIGH", "KS_Valid": true
    }
  ]
}
```

When no MAF is provided: `"mfsd": { "enabled": false }`

| Field Group | Columns | Description |
|------------|---------|-------------|
| Variant info | `Chrom`, `Pos`, `Ref`, `Alt`, `VarType` | Variant coordinates and type |
| Counts | `REF_Count`, `ALT_Count`, `NonREF_Count`, `N_Count`, `Total_Count` | Raw fragment counts per allele class |
| GC-Weighted | `REF_Weighted`, `ALT_Weighted`, `NonREF_Weighted`, `N_Weighted`, `VAF_GC_Corrected` | GC-bias corrected counts/VAF |
| Log-Likelihood | `ALT_LLR`, `REF_LLR` | Log-likelihood ratios (duplex/low-N) |
| Mean Sizes | `REF_MeanSize`, `ALT_MeanSize`, `NonREF_MeanSize`, `N_MeanSize` | Mean fragment size per allele class |
| KS Tests | `Delta_*/KS_*/KS_Pval_*` | Kolmogorov–Smirnov distance + p-value (6 pairings) |
| Derived | `VAF_Proxy`, `Error_Rate`, `N_Rate`, `Size_Ratio`, `Quality_Score` | Summary biomarker values |
| Flags | `ALT_Confidence`, `KS_Valid` | Quality flags |

---

### Region MDS (Per-Exon/Gene Motif Diversity Score)

Per-exon and gene-level MDS from `region-mds` command. Present when `--region-mds` is run.

```json
"region_mds": {
  "n_exons": 1820,
  "mds_exon_mean": 0.512,
  "mds_exon_std": 0.041,
  "exon": [
    {
      "gene": "ATM",
      "name": "ATM:exon1",
      "chrom": "11",
      "start": 108093558,
      "end": 108093795,
      "strand": "+",
      "n_fragments": 234,
      "mds": 0.531
    }
  ],
  "n_genes": 146,
  "mds_e1_mean": 0.504,
  "gene": [
    {
      "gene": "ATM",
      "n_exons": 23,
      "n_fragments": 5210,
      "mds_mean": 0.519,
      "mds_e1": 0.531,
      "mds_std": 0.038
    }
  ]
}
```

| Field | Level | Description |
|-------|-------|-------------|
| `exon[]` | Exon | Per-exon records from `{sample}.MDS.exon.tsv` |
| `exon[].mds` | Exon | MDS for this exon/target region |
| `exon[].name` | Exon | Exon identifier (`gene:exonN` for WGS, target name for panel) |
| `mds_exon_mean` | Summary | Mean MDS across all exons |
| `mds_exon_std` | Summary | Std dev of per-exon MDS |
| `gene[]` | Gene | Per-gene records from `{sample}.MDS.gene.tsv` |
| `gene[].mds_mean` | Gene | Mean MDS across all exons of this gene |
| `gene[].mds_e1` | Gene | MDS of the first exon (E1) — promoter proxy |
| `mds_e1_mean` | Summary | Mean E1 MDS across all genes |

---

### UXM (Methylation)

CpG methylation features. Only present when methylation BAM input is provided.

```json
"uxm": {
  "enabled": true,
  "data": [
    { "region": "chr1:0-500000", "U_fraction": 0.82, "X_fraction": 0.05, "M_fraction": 0.13 }
  ]
}
```

When no methylation input: `"uxm": { "enabled": false }`

---

## ML Integration Example

```python
import json
import pandas as pd

with open("sample.features.json") as f:
    features = json.load(f)

# FSC off-target (genome-wide)
df_fsc = pd.DataFrame(features["fsc"]["off_target"])

# WPS — reconstruct DataFrame from columnar arrays
wps = features["wps"]
df_wps = pd.DataFrame({"region_id": wps["regions"], "wps_nuc": wps["wps_nuc"], "wps_tf": wps["wps_tf"]})
print(f"WPS anchors: {len(df_wps)}")

# Region MDS — per-gene E1 MDS
df_mds_gene = pd.DataFrame(features["region_mds"]["gene"])
print(f"Mean E1 MDS: {features['region_mds']['mds_e1_mean']:.3f}")

# mFSD — variant-level fragment size shift
if features.get("mfsd", {}).get("enabled"):
    df_mfsd = pd.DataFrame(features["mfsd"]["variants"])
    print(f"mFSD variants: {features['mfsd']['n_variants']}")

# Motif diversity
mds_z = features["motif"].get("mds_z", 0)
print(f"MDS: {features['motif']['mds']:.4f}")
```

---

## JSON Generation with Panel Mode

For MSK-ACCESS panels:

```bash
krewlyzer run-all -i sample.bam -r hg19.fa -o out/ \
    --assay xs2 \
    --target-regions targets.bed \
    --pon-model xs2.pon.parquet \
    --generate-json
```

Panel mode adds these additional features to the JSON:
- `fsc_gene`: 146 genes (gene-level FSC)
- `fsc_region`: per-exon FSC
- `fsc_region_e1`: first-exon FSC
- `wps_panel`: 1,820 anchors
- `wps_background`: Alu stacking profiles
- `ocf.on_target`, `tfbs.on_target`, `atac.on_target`: panel-specific entropy
- PON z-scores across all features

---

## Output File Structure

Krewlyzer generates TSV/Parquet files alongside the optional unified JSON:

### Core Fragmentomics

```
out/
├── sample.FSD.tsv                   # Fragment size distribution (arm-level)
├── sample.FSD.ontarget.tsv          # Panel mode: on-target FSD
├── sample.FSR.tsv                   # Fragment size ratio (short/long)
├── sample.FSR.ontarget.tsv          # Panel mode: on-target FSR
├── sample.FSC.tsv                   # Fragment size coverage (bin-level)
├── sample.FSC.ontarget.tsv          # Panel mode: on-target FSC
├── sample.FSC.gene.tsv              # Gene-level FSC (with --assay)
├── sample.FSC.regions.tsv           # Exon-level FSC (aggregate_by='region')
├── sample.FSC.regions.e1only.tsv    # E1-only FSC (first exon per gene)
├── sample.fsc_counts.tsv            # Raw GC-binned fragment counts (pre-correction)
└── sample.correction_factors.tsv    # GC correction factors
```

### WPS (Windowed Protection Score)

```
out/
├── sample.WPS.parquet               # Per-region WPS profiles (foreground)
├── sample.WPS.panel.parquet         # Panel-specific anchors (with --assay)
└── sample.WPS_background.parquet    # Alu stacking profiles (background)
```

### Motif & Tissue-of-Origin

```
out/
├── sample.EndMotif.tsv              # 4-mer end motif frequencies (off-target)
├── sample.EndMotif.ontarget.tsv     # 4-mer end motif frequencies (on-target)
├── sample.EndMotif1mer.tsv          # Single-base end compositions (A/C/G/T)
├── sample.BreakPointMotif.tsv       # 4-mer breakpoint motif frequencies
├── sample.BreakPointMotif.ontarget.tsv
├── sample.MDS.tsv                   # Motif Diversity Score (off-target)
├── sample.MDS.ontarget.tsv          # Motif Diversity Score (on-target)
├── sample.MDS.exon.tsv              # Per-exon MDS (region-mds command)
├── sample.MDS.gene.tsv              # Per-gene aggregated MDS (region-mds command)
├── sample.OCF.tsv                   # Orientation-aware fragmentation
├── sample.OCF.ontarget.tsv          # Panel mode: on-target OCF
├── sample.OCF.offtarget.tsv         # Panel mode: off-target OCF
└── sample.OCF.sync.tsv              # Positional strand-specific profiles
```

### Region Entropy (TFBS/ATAC)

```
out/
├── sample.TFBS.tsv                  # TF binding site entropy (808 factors)
├── sample.TFBS.ontarget.tsv         # Panel mode: on-target TFBS
├── sample.TFBS.sync.tsv             # Per-TF × per-size profiles
├── sample.TFBS.ontarget.sync.tsv
├── sample.ATAC.tsv                  # ATAC-seq peak entropy (23 cancer types)
├── sample.ATAC.ontarget.tsv         # Panel mode: on-target ATAC
├── sample.ATAC.sync.tsv             # Per-tissue × per-size profiles
└── sample.ATAC.ontarget.sync.tsv
```

### Variant (mFSD)

```
out/
├── sample.mFSD.tsv                  # Per-variant fragment size metrics (46 columns)
└── sample.mFSD.distributions.tsv    # Per-variant size histograms (optional)
```

### Methylation (UXM)

```
out/
└── sample.UXM.tsv                   # CpG methylation U/X/M fractions
```

### Unified Output

```
out/
├── sample.metadata.json             # Run metadata and QC metrics
└── sample.features.json             # All features (with --generate-json)
```

!!! note
    The `--generate-json` flag produces the unified JSON **in addition to** the standard TSV/Parquet outputs.

---

## See Also

- [Pipeline Integration](../../nextflow/index.md) - `run-all` command and `--generate-json` flag
- [Input File Formats](../../reference/input-formats.md) - Custom input file specifications
- [Troubleshooting](../../resources/troubleshooting.md) - Common format issues
