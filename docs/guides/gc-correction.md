# GC Bias Correction

GC content bias is a major source of variability in cfDNA sequencing. Krewlyzer implements LOESS-based correction in Rust for high performance.

## Why GC Correction?

DNA fragments with extreme GC content (very AT-rich or GC-rich) are under-represented due to:

1. **PCR amplification bias** - GC-rich regions amplify less efficiently
2. **Sequencing bias** - Certain GC ranges have lower quality scores
3. **Capture bias** - Hybridization efficiency varies with GC

Without correction, fragment counts are confounded by GC content, obscuring biological signal.

---

## How It Works

### LOESS Regression

Krewlyzer uses **Locally Estimated Scatterplot Smoothing (LOESS)** to model the relationship between GC content and fragment count:

```
For each fragment type (short, intermediate, long):
1. Bin fragments by GC content (0.00-1.00 in 0.01 steps)
2. Fit LOESS curve: count = f(gc)
3. Compute correction factor: factor[gc] = median(count) / loess_fit[gc]
4. Apply: corrected_count = raw_count × factor[gc]
```

### Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| LOESS fraction | 0.3 | Fraction of data for local fitting |
| LOESS iterations | 3 | Robustness iterations |
| Delta | 0.01 | Smoothing delta |

---

## Fragment Length Bins

GC correction is applied per fragment length bin (17 bins, 20bp width):

| Bin | Range | Fragment Type |
|-----|-------|---------------|
| 0 | 60-79bp | Ultra-short |
| 1 | 80-99bp | Ultra-short |
| 2-4 | 100-159bp | Short |
| 5-7 | 160-219bp | Intermediate |
| 8-16 | 220-400bp | Long |

---

## Usage

### Automatic (Default)

GC correction is **enabled by default** for most tools:

```bash
krewlyzer extract -i sample.bam -r hg19.fa -o output/
# Generates: sample.correction_factors.tsv
```

### Disable GC Correction

```bash
krewlyzer fsc -i sample.bed.gz --no-gc-correct -o output/
```

### Using Pre-computed Factors

```bash
# mFSD can use factors from extract. It computes each fragment's GC from the
# reference, so the reference is required with them.
krewlyzer mfsd -i sample.bam -V variants.vcf -g reference.fa \
    --correction-factors output/sample.correction_factors.tsv \
    -o output/
```

---

## Per-Tool GC Correction

| Tool | GC Option | Source | Notes |
|------|-----------|--------|-------|
| **extract** | `--gc-correct` | Computes factors | Generates `.correction_factors.tsv` |
| **FSC** | `--gc-correct` | From extract | Via `run_unified_pipeline` |
| **FSR** | `--gc-correct` | From extract | Via `run_unified_pipeline` |
| **FSD** | `--gc-correct` | From extract | Via `run_unified_pipeline` |
| **WPS** | `--gc-correct` | From extract | Via `run_unified_pipeline` |
| **OCF** | `--gc-correct` | From extract | Via `run_unified_pipeline` |
| **Region Entropy** | `--gc-correct`, `--gc-factors` | From extract | TFBS/ATAC fragment weighting |
| **Region MDS** | | N/A | Uses raw motif counts (no GC correction) |
| **mFSD** | `--correction-factors` | Manual input | Uses pre-computed CSV |
| **motif** | N/A | N/A | No GC correction |
| **UXM** | N/A | N/A | No GC correction |

### Panel Data (--target-regions)

When `--target-regions` is provided to `extract`:
- GC model is built from **off-target** fragments only
- Avoids capture bias contamination
- Generates both `.correction_factors.*` (off-target) and `.correction_factors.ontarget.*` (on-target), in the run's `--output-format`

**Which table a feature uses.** Every consumer finds its table through one resolver
(`core.gc_assets.resolve_gc_factors`): beside the BED first, then in the output directory;
Parquet, then `.tsv.gz`, then `.tsv`. TSV keeps full precision, so all three carry identical
factors and GC correction does not depend on the output format. Each sample's metadata records
the kind of table each population used (`gc_factors_offtarget`, `gc_factors_ontarget`).

!!! warning "0.9.x cohorts"
    Up to 0.9.2, panel WPS, every TFBS and ATAC table, and gene FSC looked only for a plain
    `.tsv` factor table. Runs with `--output-format parquet`, or with `--compress` (the Nextflow
    default), write none, so those features ran **without GC correction**, while PON donors
    (built with TSV forced) were corrected. Gene FSC was uncorrected for samples and donors
    alike. Re-run those cohorts; `build-pon --from-outputs` refuses sample directories from
    before the fix.

---

## Building GC Reference Assets

The `build-gc-reference` command pre-computes reference GC data:

```bash
# Standard (WGS) mode
krewlyzer build-gc-reference hg19.fa -o data/gc/ -e hg19-blacklist.bed

# Panel mode (generates both WGS and on-target assets)
krewlyzer build-gc-reference hg19.fa -o data/gc/ -T msk_targets.bed
```

### Output Files

| Mode | File | Description |
|------|------|-------------|
| Standard | `valid_regions_{genome}.bed.gz` | 100kb bins for GC estimation |
| Standard | `ref_genome_GC_{genome}.parquet` | Expected fragment counts |
| Panel | `valid_regions_{genome}.ontarget.bed.gz` | Bins overlapping targets |
| Panel | `ref_genome_GC_{genome}.ontarget.parquet` | On-target expected counts |

### Options

| Option | Description |
|--------|-------------|
| `-o, --output` | Output directory (required) |
| `-e, --exclude-regions` | Blacklist BED to exclude |
| `-T, --target-regions` | Target BED for panel mode |
| `--bin-size` | Bin size in bp (default: 100kb) |
| `-n, --genome-name` | Genome name (default: from FASTA) |

## Correction Factors File

The `extract` command generates `{sample}.correction_factors.tsv` (or `.tsv.gz` with `--compress`):

```tsv
length_bin_min	length_bin_max	gc_percent	observed	expected	correction_factor
75	80	36	53	1281884	1.0000
75	80	39	56	1259020	1.0000
...
```

| Column | Description |
|--------|-------------|
| `length_bin_min` | Fragment length bin lower bound (bp) |
| `length_bin_max` | Fragment length bin upper bound (bp) |
| `gc_percent` | GC content percentage (0–100) |
| `observed` | Raw fragment count |
| `expected` | Expected count from LOESS model |
| `correction_factor` | `expected / observed` — multiply raw counts by this |

---

## PON-based Hybrid Correction

When a PON model is provided, correction uses a **hybrid approach**:

```mermaid
flowchart LR
    obs[Observed Counts] --> pon[PON Correction]
    pon --> |"observed / pon_expected"| res[Residual LOESS]
    res --> |"pon_corrected / residual"| final[Final Counts]
```

**Algorithm:**
1. **PON correction**: Divide by assay-specific expected coverage
2. **Residual LOESS**: Fit sample-specific residual bias
3. **Final**: Divide PON-corrected counts by residual

This removes both assay-wide and sample-specific GC effects.

---

## Implementation Details

### Rust Module: `gc_correction.rs`

Key structures:

```rust
/// Fragment length bins (17 bins, 60-400bp)
pub struct LengthBin(u8);

/// GC correction factors lookup
pub struct CorrectionFactors {
    factors: HashMap<(LengthBin, u8), f64>,
    stats: HashMap<(LengthBin, u8), CorrectionBinStats>,
}
```

### Python Interface

```python
from krewlyzer import _core

# Compute and save correction factors
_core.gc.compute_and_write_gc_factors(
    bed_path="sample.bed.gz",
    gc_reference_path="gc_reference.parquet",
    output_path="correction_factors.tsv",
    output_format="both",   # "tsv", "parquet", or "both"
    compress=True,           # gzip-compress TSV output
)
```

---

## The redesigned GC model (library only)

!!! warning "Not yet used by any feature (spec Step 2)"
    Everything in this section is library code in `rust/src/gc_model.rs` and
    `rust/src/gc_reference.rs`. Every GC-weighted feature still reads the correction factors
    described above, until spec Step 2 moves the GC tally and the lookups onto the shared fragment
    record. Nothing here changes an output.

**One GC definition.** A fragment's cell is its exact length and the integer GC count of its
insert with `k` bp trimmed from each end (`gc_key`; `k` defaults to 10 until Step 3 chooses it
from {0, 2, 5, 10}). A, C, G and T count in either case. Any other base in the trimmed window
(N, IUPAC) means no key, and so does a fragment shorter than `2k + 20` bp. The observed tally,
the expected counts and the lookup all use this one rule. Lengths run from 20 to 600 bp.

**Expected counts** (`ExpectedTable`) are exhaustive: every start position at every exact length
in a `Territory`, counted from prefix sums. A window counts only if it lies wholly inside one
interval. The territory's contigs resolve through `contigs::ContigNames`, and its SHA-256
includes each contig's length, so the same coordinates on another assembly have another digest.

**The support** comes from the reference, not from donors. At each length it is the fewest GC
cells holding 99.5% of the territory's expected windows.

**The fit** (`gc_model::fit`):

1. At each length, `r = p_obs(GC | length) / p_exp(GC | length)` over the support.
2. `r` is smoothed over (length, GC fraction) by a Gaussian-kernel **local-linear** fit,
   weighted by the expected count at the sample's depth. The bandwidths are configuration
   (defaults 5 bp and 0.02 GC fraction, truncated at 4 SD) until Step 3 chooses them by held-out
   residual. The remaining curvature bias, about `(κh)²/2`, is estimated and reported.
3. Raw weight = `1 / smoothed r`, capped to [1/10, 10].
4. A cell outside the support takes its length's nearest in-support value. Its fragments never
   enter the smoother.
5. Each length's weights are divided by their fragment-weighted mean over the whole tally, so
   every length keeps its total weight: `Σ obs·w = Σ obs`.

Below 100,000 fragments, or when a length has no support cell, the fit returns "cannot fit". It
never returns weights of 1. `weight(length, gc)` returns no value, rather than 1.0, for an
unresolved GC or a cell outside the model; what to do then is Step 3's policy. `FitQc` reports
the shares of weight in clamped and capped cells, the lengths left unobserved or unfitted, and
the fragments outside the length range, which are not fitted (Step 2 decides what they get).

**Provenance.** A model records its genome, territory and support digests, trim, GC definition
version, fit configuration and krewlyzer version. Its compatibility hash covers all of these
except the krewlyzer version. The Parquet loader refuses a table with no provenance or with
another GC definition version: it does not reinterpret it.

---

## When to Disable

Disable GC correction (`--no-gc-correct`) when:

- Comparing raw signal across samples with known GC differences
- Validating uncorrected patterns
- Working with very small regions where LOESS may be unstable

For most applications, **keep correction enabled** (default).
