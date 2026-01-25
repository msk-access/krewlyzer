# Troubleshooting

## Common Issues

### File Not Found Error
**Error**: `FileNotFoundError: [Errno 2] No such file or directory`

**Solution**: Ensure all input files (BAM, FASTA, BED) exist and paths are correct. Use absolute paths.

---

### Permission Error
**Error**: `PermissionError: [Errno 13] Permission denied`

**Solution**: Check write permissions for the output directory.

---

### Missing Dependencies
**Error**: `ModuleNotFoundError: No module named '...'`

**Solution**: Ensure Krewlyzer is installed correctly:
```bash
uv pip install krewlyzer
```
Or use the Docker image for a complete environment.

---

### Reference Mismatch
**Issue**: Results look wrong or empty.

**Solution**: Ensure BAM files and reference FASTA are from the **same genome build** (both hg19 or both hg38). Krewlyzer defaults to hg19 for bundled data files.

---

### Memory Errors
**Issue**: Process crashes on large BAM files.

**Solutions**:
1. Increase available RAM (≥16GB recommended)
2. Reduce thread count: `--threads 4`
3. Process chromosomes separately:
```bash
krewlyzer extract -i sample.bam -r hg19.fa -o output/ --chromosomes chr1,chr2
```

---

## Duplex/Consensus BAM Issues {#duplex-bam}

### "FILTER COMPATIBILITY WARNING: 0% of reads pass"

**Cause**: Duplex or consensus BAMs don't have proper pair flags set correctly. The default `--require-proper-pair` filter removes all reads.

**What's happening**:
- Standard BAMs: R1+R2 reads are flagged as "proper pairs"
- Duplex BAMs: Consensus reads are often single-ended or unpaired
- Result: `--require-proper-pair` (default) filters out everything

**Solution**: Disable proper pair requirement:
```bash
krewlyzer extract -i duplex.bam -r hg19.fa -o output/ \
    --no-require-proper-pair

krewlyzer run-all -i duplex.bam -r hg19.fa -o output/ \
    --no-require-proper-pair
```

### Auto-detection in run-all

The `run-all` command automatically detects this situation:

```
⚠️ FILTER COMPATIBILITY WARNING
   Only 0.00% of sampled reads would pass current filters.
   Only 0.0% of reads are marked as proper pairs.

   Suggested command:
   krewlyzer run-all -i sample.bam ... --no-require-proper-pair
```

If you see this warning, re-run with `--no-require-proper-pair`.

### Which BAM types need --no-require-proper-pair?

| BAM Type | Proper Pairs? | Need Flag? |
|----------|---------------|------------|
| Standard WGS | ✅ Yes | No |
| Standard Panel | ✅ Yes | No |
| Duplex/UMI | ❌ No | **Yes** |
| Consensus | ❌ No | **Yes** |
| Single-end | ❌ No | **Yes** |

---

## GC Correction Issues {#gc-correction}

### "GC correction assets not found"
**Cause**: Missing GC reference files for the specified genome.

**Solutions**:
1. Verify genome build matches bundled assets:
   ```bash
   krewlyzer extract -i sample.bam -r hg19.fa -o output/ -G hg19
   ```
2. Use `--no-gc-correct` to skip GC correction:
   ```bash
   krewlyzer fsc -i sample.bed.gz -o output/ --no-gc-correct
   ```

### Correction factors look wrong
**Cause**: Insufficient coverage or extreme GC bias.

**Diagnosis**: Check `correction_factors.csv`:
- Factors should be ~0.5-2.0 for most bins
- Extreme factors (>10 or <0.1) indicate problems

**Solutions**:
1. Increase coverage (>1M fragments recommended)
2. Check BAM for quality issues (duplicates, low MAPQ)

---

## hg38 / GRCh38 Issues {#hg38}

### "OCF regions not available for hg38"
**Cause**: Bundled OCF regions only exist for hg19/GRCh37.

**Solutions**:
1. Provide a custom OCR file:
   ```bash
   krewlyzer ocf -i sample.bed.gz -o output/ -G hg38 \
       -r custom_ocr_regions.bed.gz
   ```
2. In `run-all`, OCF is automatically skipped for hg38 with a warning

### Missing assets for hg38
**Cause**: Some bundled assets only exist for hg19.

**Affected**:
- OCF regions (hg19 only)
- Some methylation markers

**Solution**: Use `-G hg19` if your data supports both, or provide custom files.

---

## Asset Validation {#asset-validation}

### Validating custom files before running

If you're using custom override files (e.g., `--arms-file`, `--wps-anchors`), validate their format first:

```bash
# Validate specific files
krewlyzer validate --gene-bed my_genes.bed
krewlyzer validate --arms-bed my_arms.bed --wps-anchors my_anchors.bed
```

Validation checks column counts, data types, and format requirements. If a file fails validation, you'll see a detailed error with:
- Expected format and columns
- Line number of the first error
- Example of correct format
- Link to documentation

See [Input File Formats](advanced/input-formats.md) for complete format specifications.

### Validating bundled assets

After updating krewlyzer or modifying data files, validate bundled assets:

```bash
# Validate all bundled assets for a genome
krewlyzer validate --genome hg19

# Also validate assay-specific assets
krewlyzer validate --genome hg19 --assay xs2
```

### Validation during analysis

For production pipelines, use `--validate-assets` to verify bundled files before each run:

```bash
# Recommended: validate once when setting up a new environment
krewlyzer run-all -i sample.bam -r hg19.fa -o output/ \
    --validate-assets

# Or with build-pon
krewlyzer build-pon samples.txt -a xs2 -r hg19.fa -o pon.parquet \
    --validate-assets
```

> [!TIP]
> Use `--validate-assets` after krewlyzer updates or when debugging unexpected errors.
> For routine processing, validation overhead (~100ms) can be skipped.

### Common format errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Expected at least N columns" | Missing columns | Ensure file is tab-delimited with all required columns |
| "Invalid arm format" | Wrong arm name | Use patterns like `1p`, `1q`, `22p` (not `Arm1`) |
| "Expected numeric" | Non-numeric coordinates | Check start/end columns are integers |

---

## PON Model Issues {#pon}

### "Assay mismatch warning"
**Cause**: PON model built with different assay than sample.

**Impact**: Results may be less accurate but processing continues.

**Solution**: Use a PON model built from the same assay:
```bash
krewlyzer fsc -i sample.bed.gz -o output/ \
    -P matching_assay.pon.parquet
```

### PON normalization looks wrong
**Diagnosis**: Check log-ratio columns:
- `*_logR` should be centered around 0 for healthy samples
- Extreme values (>5 or <-5) indicate mismatch

**Solutions**:
1. Verify PON matches assay and genome build
2. Run without PON to get raw values:
   ```bash
   krewlyzer fsc -i sample.bed.gz -o output/
   ```

---

## Panel Data Issues {#panel}

### FSC/FSR output is all zeros
**Cause**: Default 100kb bins don't overlap panel targets.

**Solution**: Provide custom bins matching your panel:
```bash
krewlyzer fsc -i sample.bed.gz -o output/ \
    --bin-input panel_targets.bed
```

### On-target vs off-target confusion
**Solution**: Use `--target-regions` to generate separate outputs:
```bash
krewlyzer run-all -i sample.bam -r ref.fa -o output/ \
    --target-regions panel_targets.bed
```
- `.tsv` files = off-target (unbiased)
- `.ontarget.tsv` files = on-target (capture-biased)

---

## Getting Help

If your issue isn't listed:

1. Check [GitHub Issues](https://github.com/msk-access/krewlyzer/issues)
2. Run with `--verbose` for detailed logging
3. Open a new issue with:
   - Command used
   - Error message
   - Krewlyzer version (`krewlyzer --version`)
