#!/bin/bash
#SBATCH --job-name=mfsd_diag
#SBATCH --partition=cmobic_cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=mfsd_diag_%j.log
#SBATCH --error=mfsd_diag_%j.err

# ============================================================================
# mFSD Diagnostic Test — IRIS SLURM Script
# ============================================================================
# Purpose: Run mFSD standalone on the known-hanging sample P-0000302-T03-XS1
#          with full diagnostic logging to identify the root cause of 16h hangs.
#
# Usage:
#   1. Push feature branch:  git push origin feature/mfsd-diagnostics
#   2. On IRIS:
#        cd /usersoftware/shahr2/github/krewlyzer
#        git fetch origin && git checkout feature/mfsd-diagnostics
#        pip install -e .
#   3. Submit:  sbatch scripts/test_mfsd_iris.sh
#   4. Monitor: tail -f mfsd_diag_*.log
# ============================================================================

set -euo pipefail

# --- Configuration ---
SAMPLE_ID="P-0000302-T03-XS1"
BAM="/data1/share001/share/access_12_245/X/T/XT375722-T-duplex.bam"
MAF="/data1/core006/access/production/resources/cbioportal/current/msk_solid_heme/data_mutations_extended.txt"
REF="/data1/core006/access/production/resources/reference/versions/hg19/Homo_sapiens_assembly19.fasta"
OUTDIR="$PWD/mfsd_diag_output"
KREWLYZER="/usersoftware/shahr2/github/krewlyzer"

# --- Environment ---
eval "$(micromamba shell hook --shell bash)"
micromamba activate krewlyzer  # or whichever env has krewlyzer installed

echo "=== mFSD Diagnostic Test ==="
echo "Date:      $(date)"
echo "Host:      $(hostname)"
echo "Sample:    ${SAMPLE_ID}"
echo "BAM:       ${BAM}"
echo "MAF:       ${MAF}"
echo "REF:       ${REF}"
echo "Krewlyzer: $(which krewlyzer) ($(krewlyzer --version 2>&1 || true))"
echo "Git:       $(cd ${KREWLYZER} && git log --oneline -1)"
echo ""

# --- Validate inputs ---
if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found: ${BAM}"
    exit 1
fi
if [[ ! -f "${BAM%.bam}.bai" ]] && [[ ! -f "${BAM}.bai" ]]; then
    echo "ERROR: BAM index not found for ${BAM}"
    exit 1
fi

# --- Step 1: Filter MAF (same as FILTER_MAF Mode 1) ---
mkdir -p "${OUTDIR}"
FILTERED_MAF="${OUTDIR}/${SAMPLE_ID}.filtered.maf"

echo ">>> Step 1: Filtering MAF for ${SAMPLE_ID}..."
python3 -c "
sample_id = '${SAMPLE_ID}'
sample_id_lower = sample_id.lower()
variant_count = 0
with open('${MAF}') as infile, open('${FILTERED_MAF}', 'w') as outfile:
    header_written = False
    tsb_col = None
    for line in infile:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if not header_written:
            outfile.write(line)
            header_written = True
            for i, col in enumerate(fields):
                if col == 'Tumor_Sample_Barcode':
                    tsb_col = i
                    break
            continue
        if tsb_col is not None and len(fields) > tsb_col:
            if sample_id_lower == fields[tsb_col].lower():
                outfile.write(line)
                variant_count += 1
print(f'Filtered: {variant_count} variants for {sample_id}')
"

echo ">>> Filtered MAF: ${FILTERED_MAF}"
echo ">>> Variant count: $(tail -n +2 ${FILTERED_MAF} | wc -l | tr -d ' ')"
echo ""

# --- Step 2: Run mFSD with full diagnostics ---
echo ">>> Step 2: Running mFSD with RUST_LOG=info ..."
echo ">>> Start time: $(date)"
echo ""

RUST_LOG=info krewlyzer mfsd \
    -i "${BAM}" \
    -V "${FILTERED_MAF}" \
    -o "${OUTDIR}" \
    -s "${SAMPLE_ID}" \
    -g "${REF}" \
    --duplex \
    --verbose \
    --threads 4

echo ""
echo ">>> End time: $(date)"
echo ">>> Output:"
ls -lh "${OUTDIR}/"

# --- Step 3: Verify output ---
MFSD_TSV="${OUTDIR}/${SAMPLE_ID}.mFSD.tsv"
if [[ -f "${MFSD_TSV}" ]]; then
    echo ">>> mFSD output columns: $(head -1 ${MFSD_TSV} | tr '\t' '\n' | wc -l | tr -d ' ')"
    echo ">>> mFSD data rows: $(tail -n +2 ${MFSD_TSV} | wc -l | tr -d ' ')"
    echo ">>> SUCCESS"
else
    echo ">>> FAILURE: ${MFSD_TSV} not produced"
fi
