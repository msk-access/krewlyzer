# Motif-based Feature Extraction

**Command**: `krewlyzer motif`

## Purpose
Extracts end motif, breakpoint motif, and Motif Diversity Score (MDS) from sequencing fragments.

## Biological Context
Motif analysis of cfDNA fragment ends can reveal tissue-of-origin, nucleosome positioning, and nuclease activity. MDS quantifies motif diversity, which may be altered in cancer. See [Zhou et al., 2020](../citation.md#motif) for details on jagged ends and nuclease biology.

## Usage
```bash
krewlyzer motif /path/to/input.bam -g /path/to/reference.fa -o /path/to/output_dir \
    -k 4 --threads 4
```

## Output
- `{Sample}.EndMotif.tsv`: Frequency of K-mer sequences at fragment 5' ends.
- `{Sample}.BreakPointMotif.tsv`: Frequency of K-mer sequences flanking the breakpoint (context + fragment start).
- `{Sample}.MDS.tsv`: Motif Diversity Score (Shannon entropy of end motifs).

## Method
Motif frequencies are calculated by analyzing the K-mer sequences at the 5'-ends of fragments.
- **End Motif**: Raw sequence at the end of the fragment.
- **Breakpoint Motif**: Sequence context from the reference genome flanking the fragment start point.
- **MDS**: Normalized Shannon entropy of the End Motif distribution.

## Calculation Details

### Motif Diversity Score (MDS)
MDS quantifies the randomness of the 4-mer end motifs.

$$ MDS = \frac{ - \sum_{i} p_i \log_2(p_i) }{ \log_2(4^k) } $$

Where:
- $p_i$ is the frequency of the $i$-th motif.
- $k$ is the k-mer length (default 4).
- The denominator normalizes the score to $[0, 1]$.

## Clinical Interpretation

### MDS Values
| MDS | Interpretation |
|-----|----------------|
| **High (~1.0)** | Random/diverse fragmentation (healthy-like) |
| **Low** | Stereotyped fragmentation (possible tumor signal) |

### Healthy vs Cancer
| Metric | Healthy | Cancer (ctDNA) |
|--------|---------|----------------|
| MDS | Higher (random) | Lower (stereotyped) |
| Jagged ends | Lower | **Higher** |
| Specific motifs | Baseline | Enriched cancer-associated motifs |

### Biological Basis
- cfDNA fragmentation is driven by nucleases (DNASE1, DNASE1L3)
- ~87.8% of cfDNA molecules have jagged (single-stranded) ends
- Tumor-derived fragments show higher jaggedness than wild-type

> **Reference:** See [Citation & Scientific Background](../citation.md#motif) for detailed paper summary.
