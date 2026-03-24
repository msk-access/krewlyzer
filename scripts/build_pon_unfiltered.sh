#!/bin/bash
#SBATCH --job-name=krewlyzer_pon_xs2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --partition=cmobic_cpu
#SBATCH --output=pon_xs2_unfiltered_%j.log
#SBATCH --error=pon_xs2_unfiltered_%j.err

# ==============================================================================
# Build PON from unfiltered BAMs (~3000-6000x, 5-10GB each, 47 samples)
#
# Resource math:
#   - 200GB RAM / ~25GB per sample = 4 parallel samples (with headroom)
#   - 50 CPUs  / 4 parallel        = ~12 threads per sample
#   - 47 samples / 4 parallel × ~20 min each ≈ 4 hours
# ==============================================================================

set -euo pipefail

# Activate environment with krewlyzer
eval "$(micromamba shell hook --shell bash)"
micromamba activate pygbcms

echo "=== Krewlyzer PON Build (Unfiltered) ==="
echo "Job ID:    ${SLURM_JOB_ID}"
echo "Node:      $(hostname)"
echo "CPUs:      ${SLURM_CPUS_PER_TASK}"
echo "Memory:    ${SLURM_MEM_PER_NODE}MB"
echo "Start:     $(date)"
echo "========================================="

# --- Configuration ---
SAMPLE_LIST="/data1/shahr2/shahr2/test/krewlyzer/xs2_allUniq_pon.txt"
ASSAY="xs2"
REFERENCE="/data1/core006/access/production/resources/reference/versions/hg19/Homo_sapiens_assembly19.fasta"
OUTPUT_DIR="/data1/shahr2/shahr2/test/krewlyzer/pon/xs2"
OUTPUT="${OUTPUT_DIR}/xs2.all_unique.pon.parquet"

if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir -p "${OUTPUT_DIR}"
    echo "Created output directory: ${OUTPUT_DIR}"
fi

# Log sample count
N_SAMPLES=$(grep -cv '^\s*$' "${SAMPLE_LIST}" || true)
echo "Samples:   ${N_SAMPLES}"
echo "Output:    ${OUTPUT}"
echo ""

# --- Run build-pon ---
# 4 parallel samples × ~12 threads each = 48 threads used
# Each sample peaks at ~20-25GB for 3000-6000x unfiltered BAMs
krewlyzer build-pon \
    "${SAMPLE_LIST}" \
    --assay "${ASSAY}" \
    -r "${REFERENCE}" \
    -o "${OUTPUT}" \
    --threads 50 \
    --parallel-samples 4 \
    --memory-per-sample 25 \
    --sample-timeout 0 \
    -v

echo ""
echo "=== PON Build Complete ==="
echo "Output: ${OUTPUT}"
echo "End:    $(date)"