# Motif-based Feature Extraction

**Command**: `krewlyzer motif`

!!! info "Plain English"
    Motif analysis looks at the 4-letter DNA sequences at fragment ends.
    Different enzymes cut DNA at different sequences—tumors have more diverse cutting patterns.

    **Key metric**: MDS (Motif Diversity Score) - **higher MDS = more abnormal cutting = potential tumor signal**

---

## Purpose
Extracts end motif, breakpoint motif, and Motif Diversity Score (MDS) from sequencing fragments.

---

## Processing Flowchart

```mermaid
flowchart LR
    BAM[BAM File] --> RUST[Rust Backend]
    REF[Reference FASTA] --> RUST
    RUST --> EM["EndMotif.tsv"]
    RUST --> BM["BreakPointMotif.tsv"]
    RUST --> MDS["MDS.tsv"]
    
    subgraph "With --target-regions"
        RUST --> EM_ON["EndMotif.ontarget.tsv"]
        RUST --> BM_ON["BreakPointMotif.ontarget.tsv"]
    end
```

---

## Biological Context
Motif analysis of cfDNA fragment ends reveals tissue-of-origin, nucleosome positioning, and nuclease activity. MDS quantifies motif diversity, which may be altered in cancer. See [Zhou et al., 2020](../../resources/citation.md#motif) for details.

---

## Usage
```bash
krewlyzer motif -i /path/to/input.bam -r /path/to/reference.fa -o /path/to/output_dir \
    -k 4 --threads 4
```

## Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | PATH | *required* | Input BAM file |
| `--reference` | `-r` | PATH | *required* | Reference genome FASTA (indexed) |
| `--output` | `-o` | PATH | *required* | Output directory |
| `--genome` | `-G` | TEXT | hg19 | Genome build (hg19/hg38) |
| `--target-regions` | `-T` | PATH | | Target BED (for on/off-target motifs) |
| `--skip-target-regions` | | FLAG | | Force WGS mode (ignore bundled targets) |
| `--assay` | `-A` | TEXT | | Assay code (xs1/xs2) for bundled assets |
| `--pon-model` | `-P` | PATH | | PON model for MDS z-score computation |
| `--pon-variant` | | TEXT | all_unique | PON variant: `all_unique` or `duplex` |
| `--skip-pon` | | FLAG | | Skip PON z-score normalization |
| `--kmer` | `-k` | INT | 4 | K-mer size for motif extraction |
| `--chromosomes` | | TEXT | | Comma-separated chromosomes to process |
| `--sample-name` | `-s` | TEXT | | Override sample name |
| `--require-proper-pair` | | FLAG | True | Require proper pairs (disable for duplex) |
| `--verbose` | `-v` | FLAG | | Enable verbose logging |
| `--threads` | `-t` | INT | 0 | Number of threads (0=all) |

---

## Output Files

| File | Description |
|------|-------------|
| `{sample}.EndMotif.tsv` | K-mer frequencies at fragment 5' ends |
| `{sample}.BreakPointMotif.tsv` | K-mer frequencies flanking breakpoints |
| `{sample}.MDS.tsv` | Motif Diversity Score |
| `{sample}.EndMotif1mer.tsv` | 1-mer frequencies with C-end fraction (Jagged Index) |

---

## Jagged Index (1-mer End Motifs)

The Jagged Index captures the fraction of fragment ends terminating with Cytosine (C), which is elevated in tumor-derived cfDNA with "jagged" single-stranded overhangs.

### Output: EndMotif1mer.tsv

```tsv
base    count     fraction
A       661428    0.188214
C       1261449   0.358953    # ← C-end fraction
G       849433    0.241712
T       741931    0.211121
# c_fraction    0.358953
# entropy       1.952998
# c_bias        0.108953
```

### Metrics

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| `c_fraction` | C_count / total_count | Fraction of C-ending fragments |
| `entropy` | Shannon entropy (0-2 bits) | Randomness of 1-mer distribution |
| `c_bias` | c_fraction - 0.25 | Deviation from expected 25% |

### C-End Fraction Interpretation

| c_fraction | c_bias | Interpretation |
|------------|--------|----------------|
| 0.25-0.30 | 0 to +0.05 | Normal / healthy-like |
| 0.30-0.35 | +0.05 to +0.10 | Mildly elevated (possible tumor) |
| >0.35 | >0.10 | **Elevated** (tumor signal) |

### Biological Basis

cfDNA fragmentation by DNASE1L3 produces fragments with single-stranded 5' overhangs ("jagged ends"). Tumor-derived cfDNA shows:
- **~87.8% jagged ends** (vs lower in healthy)
- **Higher C-end fraction** due to preferential C-terminal cutting

!!! tip
    The C-end fraction complements MDS for detecting tumor-derived cfDNA.
    Use both metrics together for improved sensitivity.

---

## Formulas

### Motif Diversity Score (MDS)

MDS quantifies the randomness of 4-mer end motifs using normalized Shannon entropy:

$$
\text{MDS} = \frac{-\sum_{i} p_i \times \log_2(p_i)}{\log_2(4^k)}
$$

**Variables:**
- $p_i$ = frequency of the i-th motif
- $k$ = k-mer length (default: 4)
- Result range: $[0, 1]$

**Interpretation:**

| MDS Value | Meaning |
|-----------|---------|
| ~1.0 | Random/diverse (healthy-like) |
| < 0.8 | Stereotyped (possible tumor signal) |

---

## PON Normalization

When `--pon-model` is provided, MDS output includes z-score normalization:

```bash
krewlyzer motif -i sample.bam -r hg19.fa -o output/ \
    --pon-model healthy_cohort.pon.parquet
```

### Output Columns with PON

| Column | Formula | Description |
|--------|---------|-------------|
| `MDS` | Shannon entropy | Raw score (0-1) |
| `mds_z` | `(MDS - PON_mean) / PON_std` | Z-score vs healthy baseline |

### Z-Score Interpretation

| mds_z | Meaning |
|-------|---------|
| -2 to +2 | Within normal range |
| < -2 | **Abnormally low diversity** (possible tumor) |
| > +2 | Rare (check data quality) |

!!! tip
    Lower MDS (negative z-score) indicates stereotyped cutting patterns often associated with tumor-derived cfDNA.

---

## Panel Mode (--target-regions)

When `--target-regions` is provided, motif analysis produces **separate outputs** for on-target and off-target fragments:

```bash
krewlyzer motif -i sample.bam -r hg19.fa -o output/ \
    --target-regions MSK-ACCESS_targets.bed
```

### Outputs in Panel Mode

| File | Contents | Use Case |
|------|----------|----------|
| `{sample}.EndMotif.tsv` | **Off-target** fragments | Unbiased global motif signal |
| `{sample}.EndMotif.ontarget.tsv` | **On-target** fragments | Local capture region analysis |
| `{sample}.MDS.tsv` | Off-target MDS | Primary biomarker |

!!! important
    **Off-target = unbiased** – preferred for fragmentomics biomarkers.  
    **On-target = capture-biased** – use cautiously; reflects library prep artifacts.

### Why Split?

```mermaid
flowchart TB
    subgraph "On-Target (Capture Bias)"
        CAP["Hybridization probes"] --> BIAS["GC & motif bias"]
        BIAS --> LOCAL["Local signal only"]
    end
    
    subgraph "Off-Target (Unbiased)"
        WGS["Random cfDNA"] --> PURE["True biological signal"]
        PURE --> GLOBAL["Global fragmentomics"]
    end
```

For panel data (MSK-ACCESS), on-target fragments have capture bias from hybridization probes. Off-target reads represent unbiased cfDNA and should be used for motif biomarkers.

---

## Clinical Interpretation

| Metric | Healthy | Cancer (ctDNA) |
|--------|---------|----------------|
| MDS | Higher (random) | Lower (stereotyped) |
| Jagged ends | Lower | **Higher** (87.8% jagged) |
| Specific motifs | Baseline | Cancer-associated enriched |

### Biological Basis
- cfDNA fragmentation driven by nucleases (DNASE1, DNASE1L3)
- ~87.8% of cfDNA molecules have jagged (single-stranded) ends
- Tumor-derived fragments show higher jaggedness than wild-type

---

## See Also

- [Citation & Scientific Background](../../resources/citation.md#motif) - Zhou et al. paper
- [Troubleshooting](../../resources/troubleshooting.md) - Common issues
