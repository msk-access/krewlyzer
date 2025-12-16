# Citation & Scientific Background

If you use Krewlyzer in your work, please cite this repository and the relevant methods papers below.

## Primary Literature

Krewlyzer implements or adapts methods from the following foundational papers in cfDNA fragmentomics:

---

### OCF — Orientation-aware Fragmentation {#ocf}

> **Sun K, Jiang P, Chan KC, et al.** Orientation-aware plasma cell-free DNA fragmentation analysis in open chromatin regions informs tissue of origin. *Genome Res.* 2019;29(3):418-427. [DOI](https://doi.org/10.1101/gr.242719.118)

**Key Concept:** OCF measures differentially phased fragment ends (Upstream/Downstream) at tissue-specific open chromatin regions to infer tissue-of-origin.

**Mechanism:**
- In open chromatin → nucleosomes are evicted → longer linker DNA exposed
- During apoptosis → endonuclease cuts exposed linker DNA  
- Creates characteristic pattern: **U ends peak ~60bp right, D ends peak ~60bp left** of OCR center

**Healthy Baseline:**
- **T-cells:** Highest OCF (dominant cfDNA source)
- **Liver:** Second highest
- **Other tissues:** Near zero

**Cancer Pattern:**
| Cancer Type | OCF Change |
|-------------|------------|
| HCC (liver) | ↑ Liver OCF, correlates with tumor fraction (R=0.36) |
| Colorectal | ↑ Intestine OCF (R=0.89), ↓ T-cell OCF |
| Lung | ↑ Lung OCF, ↓ T-cell OCF |

---

### WPS — Windowed Protection Score {#wps}

> **Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J.** Cell-free DNA Comprises an In Vivo Nucleosome Footprint that Informs Its Tissues-Of-Origin. *Cell.* 2016;164(1-2):57-68. [DOI](https://doi.org/10.1016/j.cell.2015.11.050)

**Key Concept:** WPS quantifies nucleosome occupancy by comparing fragments that span a protection window vs those ending within it.

**Formula:**
$$WPS(k) = N_{spanning}(k) - N_{ends}(k)$$

**Interpretation:**
| WPS Value | Meaning |
|-----------|---------|
| **Positive** | Nucleosome present (DNA protected) |
| **~Zero** | Transitional region |
| **Negative** | Open chromatin (nucleosome-free) |

**Healthy vs Cancer:**
- Nucleosome patterns are cell-type specific → infer tissue-of-origin
- Cancer: Aberrant nucleosome positioning at oncogene/TSG promoters
- Loss of 10bp periodicity at dysregulated genes

---

### FSC/FSR — Fragment Size Coverage & Ratio (DELFI) {#fsr}

> **Cristiano S, et al.** Genome-wide cell-free DNA fragmentation in patients with cancer. *Nature.* 2019;570(7761):385-389. [DOI](https://doi.org/10.1038/s41586-019-1272-6)

> **Mouliere F, Chandrananda D, et al.** Enhanced detection of circulating tumor DNA by fragment size analysis. *Sci Transl Med.* 2018;10(466):eaat4921. [DOI](https://doi.org/10.1126/scitranslmed.aat4921)

**Key Concept:** DELFI (DNA Evaluation of Fragments for earLy Interception) analyzes short/long fragment ratios genome-wide for cancer detection.

**Fragment Classes:**
| Class | Size Range | Origin |
|-------|------------|--------|
| Short | 100-150bp | Enriched in tumor cfDNA |
| Long | 151-220bp | Healthy/mono-nucleosomal |

**Healthy vs Cancer:**
| Metric | Healthy | Cancer |
|--------|---------|--------|
| Modal peak | ~166bp | Left-shifted (~145bp) |
| Short/Long ratio | Low (baseline) | **Elevated** |
| Genome-wide variability | Minimal | Increased aberrations |

**Performance:** 57-99% sensitivity across 7 cancer types at 98% specificity (AUC=0.94)

---

### UXM — Fragment-level Methylation {#uxm}

> **Loyfer N, et al.** A DNA methylation atlas of normal human cell types. *Nature.* 2022;613(7943):355-364. [DOI](https://doi.org/10.1038/s41586-022-05580-6)

**Key Concept:** Classify each cfDNA fragment as Unmethylated (U), Mixed (X), or Methylated (M) to deconvolve cell-type contributions.

**Classification Thresholds:**
- **U:** ≤25% methylated CpGs
- **M:** ≥75% methylated CpGs  
- **X:** Between 25-75%

**Healthy cfDNA Composition:**
| Cell Type | Contribution |
|-----------|--------------|
| Megakaryocytes | ~31% |
| Granulocytes | ~30% |
| Monocytes/Macrophages | ~20% |
| Endothelial | ~6% |
| Hepatocytes | ~3% |

**Resolution:** Achieves ~0.1% detection (10x better than array-based methods)

---

### Motif / Jagged Ends {#motif}

> **Zhou Q, et al.** Detection and characterization of jagged ends of double-stranded DNA in plasma. *Genome Res.* 2020;30(8):1144-1153. [DOI](https://doi.org/10.1101/gr.261396.120)

**Key Concept:** cfDNA fragments have single-stranded "jagged" ends that vary by tissue origin and health status.

**Key Findings:**
- **87.8%** of cfDNA molecules have jagged ends
- Jaggedness relates to nuclease activity (DNASE1/DNASE1L3)
- End motif diversity reflects fragmentation patterns

**Healthy vs Cancer:**
| Metric | Healthy | Cancer (ctDNA) |
|--------|---------|----------------|
| Jaggedness | Lower | **Higher** |
| Fetal vs Maternal | Fetal has higher jaggedness | — |
| Tumor vs Wild-type | — | Tumor-derived has higher jaggedness |

**MDS (Motif Diversity Score):**
- **High (~1.0):** Random/diverse fragmentation (healthy-like)
- **Low:** Stereotyped fragmentation (possible tumor signal)

---

## Acknowledgements

Krewlyzer was developed by the **MSK-ACCESS** team at Memorial Sloan Kettering Cancer Center.

The fragmentomics methods implemented here build upon foundational work from laboratories worldwide including Dennis Lo (CUHK), Jay Shendure (UW), Victor Velculescu (JHU), and others.

