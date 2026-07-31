<!--
GENERATED FILE -- DO NOT EDIT.

Produced by scripts/build_example_report.py from fabricated tables. Edit the
generator or `krewlyzer describe-output`, then run:

    python scripts/build_example_report.py
-->

# Example output report

What [`krewlyzer describe-output`](../cli/index.md) produces for one sample.

```bash
krewlyzer describe-output RESULTS/{sample_id}/ -o report.html
```

!!! info "This example uses fabricated data"
    Every number below is invented. The tables have the right shapes, columns
    and dtypes, but the values carry no biology — this page exists to show the
    *form* of the report, not to be read as a result.

    A report from a real sample is not committed here. Identifier columns are
    redacted, but the measurements themselves are still that patient's data,
    and this site is public.

---

`RESULTS/P-0000000-T01-XS1`

**27 of 27 contracted tables present**, 118.5 KB total.

## Which way is the tumour signal?

Direction differs per axis, and getting one backwards is the commonest misreading of this output — `MDS` was documented the wrong way round for a year. A single elevated metric is a hypothesis; several agreeing across independent axes is evidence.

| Table | Measures | Cancer direction |
|---|---|---|
| `FSD` | Fragment-size histogram, 5 bp bins over 65–1000 bp, per chromosome arm. | mass shifts from the ~166 bp mononucleosomal mode toward ~145 bp |
| `FSD.ontarget` | As FSD, restricted to captured regions. | as FSD, but compare on-target to on-target only |
| `FSR` | Short:long ratio per window, with each side already divided by its PON mean — so the ratio is healthy-relative and cancels batch effects. | higher `short_long_log2` |
| `FSR.ontarget` | As FSR, over captured regions. | higher `short_long_log2` |
| `FSC` | GC-corrected coverage split into six non-overlapping size channels per genomic bin. | higher short-fragment fraction in affected bins |
| `FSC.ontarget` | As FSC, over captured regions. | higher short-fragment fraction |
| `FSC.gene` | The six size channels aggregated per gene, plus an RPKM-like `normalized_depth` usable as a copy-number proxy. | higher short-fragment fraction at affected genes |
| `FSC.regions` | The same channels per exon or capture tile — the finest FSC resolution. | higher short-fragment fraction at affected regions |
| `EndMotif` | Frequency of each of the 256 4-mers at fragment 5′ ends. | a few motifs dominate — the distribution narrows |
| `EndMotif.ontarget` | As EndMotif, over captured regions. | distribution narrows |
| `BreakPointMotif` | 4-mer frequencies at the breakpoint rather than the fragment end — the sequence context spanning the cut site. | distribution narrows |
| `BreakPointMotif.ontarget` | As BreakPointMotif, over captured regions. | distribution narrows |
| `EndMotif1mer` | Single-base composition at fragment ends. | higher C-end fraction |
| `MDS` | Motif Diversity Score: Shannon entropy of the 256 4-mer end motifs, normalised by log2(256) so it lands in [0, 1]. | **lower** MDS — stereotyped cutting |
| `MDS.ontarget` | As MDS, over captured regions. | **lower** MDS |
| `MDS.gene` | MDS per gene, plus `mds_e1` for the first exon specifically — the promoter's nucleosome-depleted region, where accessibility differences are largest and early-cancer signal is strongest (Helzer 2025). | **lower** `mds_e1` at a driver locus |
| `MDS.exon` | MDS per exon or capture tile — localises aberrant cutting. | **lower** MDS |
| `OCF.ontarget` | Orientation-aware fragmentation at tissue-specific open chromatin. | the shedding tissue rises; a **fall** in T-cell OCF is equally informative, since it reflects dilution of normal haematopoietic background by tumour DNA |
| `OCF.offtarget` | As OCF, outside captured regions — the unbiased view for panel data. | as OCF |
| `TFBS` | Shannon entropy of the fragment-size distribution at the binding sites of 808 transcription factors. | shifts either way depending on the factor — interpret per TF, and prefer z-scores |
| `TFBS.ontarget` | As TFBS, over captured regions. | per-factor; z-scores only |
| `ATAC` | The same size-entropy computation over each of 23 TCGA cancer types' accessible peaks — a panel-compatible accessibility readout, orthogonal to OCF's orientation signal. | shifts at the peaks of the matching cancer type; z-scores only |
| `ATAC.ontarget` | As ATAC, over captured regions. | z-scores only |
| `WPS` | Windowed protection score around TSS and CTCF anchors: how strongly each position is protected from cutting, which traces nucleosome placement. | protection at active promoters weakens as positioning is disrupted |
| `WPS.panel` | As WPS, over the panel's anchors. | as WPS |
| `WPS_background` | WPS stacked over Alu elements, giving a genome-wide chromatin quality readout and the nucleosome repeat length from its periodicity. | NRL drifts from the ~190 bp healthy repeat length |

*No thresholds are given. Every numeric band examined turned out to be a display default or refuted outright — the documented ATAC/TFBS entropy range flags a perfectly healthy distribution as abnormal. Directions are robust; magnitudes are cohort-specific.*

## Contents

| Table | Rows | Cols | Size | Read downstream |
|---|---:|---:|---:|:---:|
| [`FSD`](#fsd) | 4 | 8 | 5.1 KB | yes |
| [`FSD.ontarget`](#fsdontarget) | 4 | 8 | 5.1 KB | yes |
| [`FSR`](#fsr) | 4 | 6 | 4.2 KB | yes |
| [`FSR.ontarget`](#fsrontarget) | 4 | 6 | 4.2 KB | yes |
| [`FSC`](#fsc) | 4 | 16 | 10.0 KB | yes |
| [`FSC.ontarget`](#fscontarget) | 4 | 16 | 10.0 KB | yes |
| [`FSC.gene`](#fscgene) | 4 | 9 | 6.1 KB | yes |
| [`FSC.regions`](#fscregions) | 4 | 12 | 7.5 KB | yes |
| [`EndMotif`](#endmotif) | 256 | 2 | 5.2 KB | yes |
| [`EndMotif.ontarget`](#endmotifontarget) | 256 | 2 | 5.2 KB | yes |
| [`BreakPointMotif`](#breakpointmotif) | 256 | 2 | 5.2 KB | yes |
| [`BreakPointMotif.ontarget`](#breakpointmotifontarget) | 256 | 2 | 5.2 KB | yes |
| [`EndMotif1mer`](#endmotif1mer) | 4 | 3 | 2.3 KB | yes |
| [`MDS`](#mds) | 1 | 1 | 1.2 KB | yes |
| [`MDS.ontarget`](#mdsontarget) | 1 | 1 | 1.2 KB | yes |
| [`MDS.gene`](#mdsgene) | 4 | 6 | 4.0 KB | yes |
| [`MDS.exon`](#mdsexon) | 4 | 2 | 1.7 KB | yes |
| [`OCF.ontarget`](#ocfontarget) | 4 | 3 | 2.3 KB | yes |
| [`OCF.offtarget`](#ocfofftarget) | 4 | 3 | 2.3 KB | yes |
| [`TFBS`](#tfbs) | 4 | 5 | 3.4 KB | yes |
| [`TFBS.ontarget`](#tfbsontarget) | 4 | 5 | 3.4 KB | yes |
| [`ATAC`](#atac) | 4 | 5 | 3.4 KB | yes |
| [`ATAC.ontarget`](#atacontarget) | 4 | 5 | 3.4 KB | yes |
| [`WPS`](#wps) | 4 | 5 | 5.0 KB | yes |
| [`WPS.panel`](#wpspanel) | 4 | 6 | 5.6 KB | yes |
| [`WPS_background`](#wps_background) | 4 | 6 | 4.3 KB | yes |
| [`metadata`](#metadata) | 1 | 3 | 2.2 KB | yes |

## Tables

### FSD

4 rows × 8 columns · 5.1 KB · read downstream

Fragment-size histogram, 5 bp bins over 65–1000 bp, per chromosome arm. Summed across arms this is the sample's genome-wide size density — the most direct picture of the whole thesis.

**Cancer direction:** mass shifts from the ~166 bp mononucleosomal mode toward ~145 bp

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `region` | object | `chr1:0-100000` | 4 | 0 | constant by design — identifier column; its values are the join key, not a measurement |
| `65-69` | float64 | 2.03248 … 6.03842 | 4 | 0 |  |
| `70-74` | float64 | 1.91036 … 7.53682 | 4 | 0 |  |
| `75-79` | float64 | 2.48373 … 8.63722 | 4 | 0 |  |
| `80-84` | float64 | 1.89857 … 9.83996 | 4 | 0 |  |
| `85-89` | float64 | 2.03913 … 8.98808 | 4 | 0 |  |
| `90-94` | float64 | 1.53932 … 9.74764 | 4 | 0 |  |
| `total` | float64 | 28.3173 … 35.7955 | 4 | 0 | must differ between rows and between samples |

### FSD.ontarget

4 rows × 8 columns · 5.1 KB · read downstream

As FSD, restricted to captured regions. Panel capture biases fragment size, so on-target and off-target are not interchangeable.

**Cancer direction:** as FSD, but compare on-target to on-target only

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `region` | object | `chr1:0-100000` | 4 | 0 | constant by design — identifier column; its values are the join key, not a measurement |
| `65-69` | float64 | 2.03248 … 6.03842 | 4 | 0 |  |
| `70-74` | float64 | 1.91036 … 7.53682 | 4 | 0 |  |
| `75-79` | float64 | 2.48373 … 8.63722 | 4 | 0 |  |
| `80-84` | float64 | 1.89857 … 9.83996 | 4 | 0 |  |
| `85-89` | float64 | 2.03913 … 8.98808 | 4 | 0 |  |
| `90-94` | float64 | 1.53932 … 9.74764 | 4 | 0 |  |
| `total` | float64 | 28.3173 … 35.7955 | 4 | 0 | must differ between rows and between samples |

### FSR

4 rows × 6 columns · 4.2 KB · read downstream

Short:long ratio per window, with each side already divided by its PON mean — so the ratio is healthy-relative and cancels batch effects. `short_long_log2` is the ML-ready signed form, ~0 in healthy plasma.

**Cancer direction:** higher `short_long_log2`

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `region` | object | `chr1:0-100000` | 4 | 0 | constant by design — identifier column; its values are the join key, not a measurement |
| `total_count` | float64 | 68.0089 … 88.5765 | 4 | 0 | must differ between rows and between samples |
| `short_long_ratio` | float64 | 0.539347 … 6.32175 | 4 | 0 | must differ between rows and between samples |
| `short_long_log2` | float64 | -0.890714 … 2.66032 | 4 | 0 | must differ between rows and between samples |
| `short_frac` | float64 | 0.298855 … 0.755407 | 4 | 0 | must differ between rows and between samples |
| `long_frac` | float64 | 0.119493 … 0.554105 | 4 | 0 | must differ between rows and between samples |

### FSR.ontarget

4 rows × 6 columns · 4.2 KB · read downstream

As FSR, over captured regions.

**Cancer direction:** higher `short_long_log2`

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `region` | object | `chr1:0-100000` | 4 | 0 | constant by design — identifier column; its values are the join key, not a measurement |
| `total_count` | float64 | 68.0089 … 88.5765 | 4 | 0 | must differ between rows and between samples |
| `short_long_ratio` | float64 | 0.539347 … 6.32175 | 4 | 0 | must differ between rows and between samples |
| `short_long_log2` | float64 | -0.890714 … 2.66032 | 4 | 0 | must differ between rows and between samples |
| `short_frac` | float64 | 0.298855 … 0.755407 | 4 | 0 | must differ between rows and between samples |
| `long_frac` | float64 | 0.119493 … 0.554105 | 4 | 0 | must differ between rows and between samples |

### FSC

4 rows × 16 columns · 10.0 KB · read downstream

GC-corrected coverage split into six non-overlapping size channels per genomic bin. Ties fragment size to *location*, so focal signal survives a genome average that would hide it.

**Cancer direction:** higher short-fragment fraction in affected bins

> ⚠️ Not PON-normalised — the healthy expectation is not zero. Read the spread across bins, not the offset.

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `chrom` | object | `chr1` | 4 | 0 | constant by design — identifier column; its values are the join key, not a measurement |
| `start` | int64 | 0 … 3000 | 4 | 0 |  |
| `end` | int64 | 1000 … 4000 | 4 | 0 |  |
| `ultra_short` | float64 | 20.32 … 60.38 | 4 | 0 |  |
| `core_short` | float64 | 19.1 … 75.37 | 4 | 0 |  |
| `mono_nucl` | float64 | 24.84 … 86.37 | 4 | 0 |  |
| `di_nucl` | float64 | 18.99 … 98.4 | 4 | 0 |  |
| `long` | float64 | 20.39 … 89.88 | 4 | 0 |  |
| `ultra_long` | float64 | 15.39 … 97.48 | 4 | 0 |  |
| `total` | float64 | 283.17 … 357.95 | 4 | 0 | must differ between rows and between samples |
| `ultra_short_log2` | float64 | 4.41414 … 5.9397 | 4 | 0 |  |
| `core_short_log2` | float64 | 4.32912 … 6.25493 | 4 | 0 |  |
| `mono_nucl_log2` | float64 | 4.69153 … 6.44907 | 4 | 0 |  |
| `di_nucl_log2` | float64 | 4.32121 … 6.63517 | 4 | 0 |  |
| `long_log2` | float64 | 4.41886 … 6.50589 | 4 | 0 |  |
| `ultra_long_log2` | float64 | 4.03474 … 6.62176 | 4 | 0 |  |

### FSC.ontarget

4 rows × 16 columns · 10.0 KB · read downstream

As FSC, over captured regions.

**Cancer direction:** higher short-fragment fraction

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `chrom` | object | `chr1` | 4 | 0 | constant by design — identifier column; its values are the join key, not a measurement |
| `start` | int64 | 0 … 3000 | 4 | 0 |  |
| `end` | int64 | 1000 … 4000 | 4 | 0 |  |
| `ultra_short` | float64 | 20.32 … 60.38 | 4 | 0 |  |
| `core_short` | float64 | 19.1 … 75.37 | 4 | 0 |  |
| `mono_nucl` | float64 | 24.84 … 86.37 | 4 | 0 |  |
| `di_nucl` | float64 | 18.99 … 98.4 | 4 | 0 |  |
| `long` | float64 | 20.39 … 89.88 | 4 | 0 |  |
| `ultra_long` | float64 | 15.39 … 97.48 | 4 | 0 |  |
| `total` | float64 | 283.17 … 357.95 | 4 | 0 | must differ between rows and between samples |
| `ultra_short_log2` | float64 | 4.41414 … 5.9397 | 4 | 0 |  |
| `core_short_log2` | float64 | 4.32912 … 6.25493 | 4 | 0 |  |
| `mono_nucl_log2` | float64 | 4.69153 … 6.44907 | 4 | 0 |  |
| `di_nucl_log2` | float64 | 4.32121 … 6.63517 | 4 | 0 |  |
| `long_log2` | float64 | 4.41886 … 6.50589 | 4 | 0 |  |
| `ultra_long_log2` | float64 | 4.03474 … 6.62176 | 4 | 0 |  |

### FSC.gene

4 rows × 9 columns · 6.1 KB · read downstream

The six size channels aggregated per gene, plus an RPKM-like `normalized_depth` usable as a copy-number proxy.

**Cancer direction:** higher short-fragment fraction at affected genes

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `gene` | object | `GENE0` | 4 | 0 | constant by design — identifier column; its values are the join key, not a measurement |
| `total` | float64 | 101.624 … 301.921 | 4 | 0 |  |
| `ultra_short_ratio` | float64 | 0.05355 … 0.229732 | 4 | 0 | must differ between rows and between samples |
| `core_short_ratio` | float64 | 0.106942 … 0.278661 | 4 | 0 | must differ between rows and between samples |
| `mono_nucl_ratio` | float64 | 0.030423 … 0.336053 | 4 | 0 | must differ between rows and between samples |
| `di_nucl_ratio` | float64 | 0.022248 … 0.278156 | 4 | 0 | must differ between rows and between samples |
| `long_ratio` | float64 | 0.056897 … 0.288693 | 4 | 0 | must differ between rows and between samples |
| `ultra_long_ratio` | float64 | 0.131837 … 0.321962 | 4 | 0 | must differ between rows and between samples |
| `normalized_depth` | float64 | 0.573109 … 2.26104 | 4 | 0 | must differ between rows and between samples |

### FSC.regions

4 rows × 12 columns · 7.5 KB · read downstream

The same channels per exon or capture tile — the finest FSC resolution. Carries `strand`, `is_e1` and `is_alt_e1` so promoter-proximal regions can be selected without re-deriving them.

**Cancer direction:** higher short-fragment fraction at affected regions

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `gene` | object | `GENE0` | 4 | 0 | constant by design — identifier column; its values are the join key, not a measurement |
| `total` | float64 | 101.624 … 301.921 | 4 | 0 | must differ between rows and between samples |
| `strand` | object | `+` | 2 | 0 | constant by design — copied from the gene BED asset rather than measured from the sample, so it is identical across every sample by construction |
| `is_e1` | object | `1` | 2 | 0 | constant by design — copied from the gene BED asset rather than measured from the sample, so it is identical across every sample by construction |
| `is_alt_e1` | object | `1` | 2 | 0 | constant by design — copied from the gene BED asset rather than measured from the sample, so it is identical across every sample by construction |
| `ultra_short_ratio` | float64 | 0.05355 … 0.229732 | 4 | 0 | must differ between rows and between samples |
| `core_short_ratio` | float64 | 0.106942 … 0.278661 | 4 | 0 | must differ between rows and between samples |
| `mono_nucl_ratio` | float64 | 0.030423 … 0.336053 | 4 | 0 | must differ between rows and between samples |
| `di_nucl_ratio` | float64 | 0.022248 … 0.278156 | 4 | 0 | must differ between rows and between samples |
| `long_ratio` | float64 | 0.056897 … 0.288693 | 4 | 0 | must differ between rows and between samples |
| `ultra_long_ratio` | float64 | 0.131837 … 0.321962 | 4 | 0 | must differ between rows and between samples |
| `normalized_depth` | float64 | 0.573109 … 2.26104 | 4 | 0 |  |

### EndMotif

256 rows × 2 columns · 5.2 KB · read downstream

Frequency of each of the 256 4-mers at fragment 5′ ends. Nuclease cutting preference is sequence-specific, and tumour cfDNA is cut differently.

**Cancer direction:** a few motifs dominate — the distribution narrows

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `Motif` | object | `AAAA` | 256 | 0 | constant by design — the 256 4-mers are a fixed alphabet |
| `Frequency` | float64 | 0.000351 … 0.006977 | 249 | 0 | must differ between rows and between samples |

### EndMotif.ontarget

256 rows × 2 columns · 5.2 KB · read downstream

As EndMotif, over captured regions.

**Cancer direction:** distribution narrows

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `Motif` | object | `AAAA` | 256 | 0 | constant by design — the 256 4-mers are a fixed alphabet |
| `Frequency` | float64 | 0.000351 … 0.006977 | 249 | 0 | must differ between rows and between samples |

### BreakPointMotif

256 rows × 2 columns · 5.2 KB · read downstream

4-mer frequencies at the breakpoint rather than the fragment end — the sequence context spanning the cut site.

**Cancer direction:** distribution narrows

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `Motif` | object | `AAAA` | 256 | 0 | constant by design — the 256 4-mers are a fixed alphabet |
| `Frequency` | float64 | 0.000351 … 0.006977 | 249 | 0 | must differ between rows and between samples |

### BreakPointMotif.ontarget

256 rows × 2 columns · 5.2 KB · read downstream

As BreakPointMotif, over captured regions.

**Cancer direction:** distribution narrows

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `Motif` | object | `AAAA` | 256 | 0 | constant by design — the 256 4-mers are a fixed alphabet |
| `Frequency` | float64 | 0.000351 … 0.006977 | 249 | 0 | must differ between rows and between samples |

### EndMotif1mer

4 rows × 3 columns · 2.3 KB · read downstream

Single-base composition at fragment ends. DNASE1L3 leaves 5′ overhangs; filling them in during library prep shifts the terminal base. A cheap complement to MDS — MDS measures diversity, this measures the chemistry of the cut end.

**Cancer direction:** higher C-end fraction

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `base` | object | `A` | 4 | 0 | constant by design — the four DNA bases |
| `count` | int64 | 10 … 40 | 4 | 0 |  |
| `fraction` | float64 | 0.057142 … 0.590047 | 4 | 0 | must differ between rows and between samples |

### MDS

1 rows × 1 columns · 1.2 KB · read downstream

Motif Diversity Score: Shannon entropy of the 256 4-mer end motifs, normalised by log2(256) so it lands in [0, 1].

**Cancer direction:** **lower** MDS — stereotyped cutting

> ⚠️ The direction is the opposite of intuition and was documented backwards for a year. Lower is the abnormal end.

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `MDS` | float64 | 5.21386 … 5.21386 | 1 | 0 | must differ between samples |

### MDS.ontarget

1 rows × 1 columns · 1.2 KB · read downstream

As MDS, over captured regions.

**Cancer direction:** **lower** MDS

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `MDS` | float64 | 5.21386 … 5.21386 | 1 | 0 | must differ between samples |

### MDS.gene

4 rows × 6 columns · 4.0 KB · read downstream

MDS per gene, plus `mds_e1` for the first exon specifically — the promoter's nucleosome-depleted region, where accessibility differences are largest and early-cancer signal is strongest (Helzer 2025).

**Cancer direction:** **lower** `mds_e1` at a driver locus

> ⚠️ `mds_e1` is NaN where the gene has no captured first exon, which on a targeted panel is most genes. NaN is not zero.

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `gene` | object | `gene_0` | 4 | 0 | constant by design — identifier column; its values are the join key, not a measurement |
| `mds_mean` | float64 | 2.03248 … 6.03842 | 4 | 0 | must differ between rows and between samples |
| `mds_e1` | float64 | 1.91036 … 7.53682 | 4 | 0 | must differ between rows and between samples |
| `mds_std` | float64 | 2.48373 … 8.63722 | 4 | 0 | must differ between rows and between samples |
| `mds_z` | float64 | 1.89857 … 9.83996 | 4 | 0 | must differ between rows and between samples |
| `mds_e1_z` | float64 | 2.03913 … 8.98808 | 4 | 0 | must differ between rows and between samples |

### MDS.exon

4 rows × 2 columns · 1.7 KB · read downstream

MDS per exon or capture tile — localises aberrant cutting.

**Cancer direction:** **lower** MDS

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `gene` | object | `gene_0` | 4 | 0 | constant by design — identifier column; its values are the join key, not a measurement |
| `mds` | float64 | 2.03248 … 6.03842 | 4 | 0 | must differ between rows and between samples |

### OCF.ontarget

4 rows × 3 columns · 2.3 KB · read downstream

Orientation-aware fragmentation at tissue-specific open chromatin. Fragments ending upstream versus downstream of a region are counted separately; a tissue actively shedding DNA produces a phased excess.

**Cancer direction:** the shedding tissue rises; a **fall** in T-cell OCF is equally informative, since it reflects dilution of normal haematopoietic background by tumour DNA

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `tissue` | object | `tissue_0` | 4 | 0 | constant by design — fixed tissue/label vocabulary from the bundled atlas |
| `OCF` | float64 | 2.03248 … 6.03842 | 4 | 0 | must differ between rows and between samples |
| `ocf_z` | float64 | 1.91036 … 7.53682 | 4 | 0 | must differ between rows and between samples |

### OCF.offtarget

4 rows × 3 columns · 2.3 KB · read downstream

As OCF, outside captured regions — the unbiased view for panel data.

**Cancer direction:** as OCF

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `tissue` | object | `tissue_0` | 4 | 0 | constant by design — fixed tissue/label vocabulary from the bundled atlas |
| `OCF` | float64 | 2.03248 … 6.03842 | 4 | 0 | must differ between rows and between samples |
| `ocf_z` | float64 | 1.91036 … 7.53682 | 4 | 0 | must differ between rows and between samples |

### TFBS

4 rows × 5 columns · 3.4 KB · read downstream

Shannon entropy of the fragment-size distribution at the binding sites of 808 transcription factors. A proxy for occupancy and local accessibility, and the highest-dimensional single feature emitted.

**Cancer direction:** shifts either way depending on the factor — interpret per TF, and prefer z-scores

> ⚠️ Absolute entropy bands are refuted: a healthy-like N(167,30) sample already exceeds the documented 'normal' range. Use z-scores only.

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `label` | object | `label_0` | 4 | 0 | constant by design — fixed tissue/label vocabulary from the bundled atlas |
| `count` | float64 | 2.03248 … 6.03842 | 4 | 0 | must differ between samples |
| `mean_size` | float64 | 1.91036 … 7.53682 | 4 | 0 | must differ between samples |
| `entropy` | float64 | 2.48373 … 8.63722 | 4 | 0 | must differ between samples |
| `z_score` | float64 | 1.89857 … 9.83996 | 4 | 0 | must differ between samples |

### TFBS.ontarget

4 rows × 5 columns · 3.4 KB · read downstream

As TFBS, over captured regions.

**Cancer direction:** per-factor; z-scores only

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `label` | object | `label_0` | 4 | 0 | constant by design — fixed tissue/label vocabulary from the bundled atlas |
| `count` | float64 | 2.03248 … 6.03842 | 4 | 0 | must differ between samples |
| `mean_size` | float64 | 1.91036 … 7.53682 | 4 | 0 | must differ between samples |
| `entropy` | float64 | 2.48373 … 8.63722 | 4 | 0 | must differ between samples |
| `z_score` | float64 | 1.89857 … 9.83996 | 4 | 0 | must differ between samples |

### ATAC

4 rows × 5 columns · 3.4 KB · read downstream

The same size-entropy computation over each of 23 TCGA cancer types' accessible peaks — a panel-compatible accessibility readout, orthogonal to OCF's orientation signal.

**Cancer direction:** shifts at the peaks of the matching cancer type; z-scores only

> ⚠️ Absolute entropy bands are refuted — see TFBS.

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `label` | object | `label_0` | 4 | 0 | constant by design — fixed tissue/label vocabulary from the bundled atlas |
| `count` | float64 | 2.03248 … 6.03842 | 4 | 0 | must differ between samples |
| `mean_size` | float64 | 1.91036 … 7.53682 | 4 | 0 | must differ between samples |
| `entropy` | float64 | 2.48373 … 8.63722 | 4 | 0 | must differ between samples |
| `z_score` | float64 | 1.89857 … 9.83996 | 4 | 0 | must differ between samples |

### ATAC.ontarget

4 rows × 5 columns · 3.4 KB · read downstream

As ATAC, over captured regions.

**Cancer direction:** z-scores only

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `label` | object | `label_0` | 4 | 0 | constant by design — fixed tissue/label vocabulary from the bundled atlas |
| `count` | float64 | 2.03248 … 6.03842 | 4 | 0 | must differ between samples |
| `mean_size` | float64 | 1.91036 … 7.53682 | 4 | 0 | must differ between samples |
| `entropy` | float64 | 2.48373 … 8.63722 | 4 | 0 | must differ between samples |
| `z_score` | float64 | 1.89857 … 9.83996 | 4 | 0 | must differ between samples |

### WPS

4 rows × 5 columns · 5.0 KB · read downstream

Windowed protection score around TSS and CTCF anchors: how strongly each position is protected from cutting, which traces nucleosome placement. Per-position vectors, not scalars.

**Cancer direction:** protection at active promoters weakens as positioning is disrupted

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `region_type` | object | `TSS` | 2 | 0 | constant by design — fixed anchor vocabulary (TSS, CTCF) |
| `wps_nuc` | object | `[8 values] 0.5214 …` | — | 0 | must differ between rows and between samples |
| `wps_tf` | object | `[8 values] 0.2607 …` | — | 0 | must differ between rows and between samples |
| `prot_frac_nuc` | object | `[8 values] 0.1564 …` | — | 0 | must differ between rows and between samples |
| `prot_frac_tf` | object | `[8 values] 0.1043 …` | — | 0 | must differ between rows and between samples |

### WPS.panel

4 rows × 6 columns · 5.6 KB · read downstream

As WPS, over the panel's anchors.

**Cancer direction:** as WPS

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `region_type` | object | `TSS` | 2 | 0 | constant by design — fixed anchor vocabulary (TSS, CTCF) |
| `wps_nuc` | object | `[8 values] 0.5214 …` | — | 0 | must differ between rows and between samples |
| `wps_tf` | object | `[8 values] 0.2607 …` | — | 0 | must differ between rows and between samples |
| `prot_frac_nuc` | object | `[8 values] 0.1564 …` | — | 0 | must differ between rows and between samples |
| `prot_frac_tf` | object | `[8 values] 0.1043 …` | — | 0 | must differ between rows and between samples |
| `local_depth` | float64 | 5.93257 … 81.1995 | 4 | 0 | must differ between rows and between samples |

### WPS_background

4 rows × 6 columns · 4.3 KB · read downstream

WPS stacked over Alu elements, giving a genome-wide chromatin quality readout and the nucleosome repeat length from its periodicity.

**Cancer direction:** NRL drifts from the ~190 bp healthy repeat length

> ⚠️ Check `nrl_at_band_limit`: where it is set, `nrl_bp` is the edge of the search band rather than a measurement — no peak was found.

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `group_id` | object | `AluY_0` | 4 | 0 | constant by design — identifier column; its values are the join key, not a measurement |
| `nrl_bp` | float64 | 184.065 … 192.077 | 4 | 0 | must differ between samples |
| `nrl_deviation_bp` | float64 | 1.91036 … 7.53682 | 4 | 0 | must differ between samples |
| `periodicity_score` | float64 | 0.248373 … 0.863722 | 4 | 0 | must differ between samples |
| `adjusted_score` | float64 | 0.189857 … 0.983996 | 4 | 0 | must differ between samples |
| `fragment_ratio` | float64 | 0.00203913 … 0.00898808 | 4 | 0 | must differ between rows and between samples |

### metadata

1 rows × 3 columns · 2.2 KB · read downstream

Run provenance and fragment totals. Also the **completion marker** — a sample without this file is dropped from the cohort silently.

| Column | Type | Range / example | Distinct | Null | Contract |
|---|---|---|---:|---:|---|
| `sample_id` | object | `<identifier — redacted>` | 1 | 0 | must differ between samples |
| `total_fragments` | int64 | 1e+06 … 1e+06 | 1 | 0 | must differ between samples |
| `genome` | object | `hg19` | 1 | 0 |  |

---

Shape, ranges and examples are measured from this sample. The one-line interpretation and direction for each table come from the output contract, so they cannot drift from it; the longer narrative — how to use each table, worked examples — is in the [output reference](https://msk-access.github.io/krewlyzer/reference/output-files/).

Values in identifier columns are redacted. Everything else is this sample's real data.
