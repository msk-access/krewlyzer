#!/bin/bash
#SBATCH --job-name=krewlyzer_pon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --partition=cmobic_cpu
#SBATCH --output=pon_%x_%j.log
#SBATCH --error=pon_%x_%j.err

# ==============================================================================
# Build one PON from unfiltered BAMs.
#
#   sbatch scripts/build_pon.sh xs2 all_unique
#   sbatch scripts/build_pon.sh xs1 duplex
#
# Replaces build_pon_unfiltered.sh, which was hardcoded to xs2 and carried a
# header claiming "47 samples" inherited from the xs1 copy it was made from --
# while the xs2 model it produced records 21. Four models built from four
# edited copies of one script is how that happens; this takes the two things
# that differ as arguments.
#
# Resource maths (unchanged, they were right):
#   200GB RAM / ~25GB per sample = 4 parallel samples, with headroom
#   50 CPUs   / 4 parallel       = ~12 threads each
#   ~20 min per sample, so 47 samples in roughly 4 hours
# ==============================================================================

set -euo pipefail

ASSAY="${1:?usage: build_pon.sh <assay: xs1|xs2> <variant: all_unique|duplex>}"
VARIANT="${2:?usage: build_pon.sh <assay: xs1|xs2> <variant: all_unique|duplex>}"

case "$ASSAY" in xs1|xs2) ;; *) echo "unknown assay: $ASSAY" >&2; exit 2 ;; esac
case "$VARIANT" in all_unique|duplex) ;; *) echo "unknown variant: $VARIANT" >&2; exit 2 ;; esac

# --- Configuration -----------------------------------------------------------
BASE="/data1/shahr2/shahr2/test/krewlyzer"
SAMPLE_LIST="${SAMPLE_LIST:-${BASE}/${ASSAY}_${VARIANT}_pon.txt}"
REFERENCE="${REFERENCE:-/data1/core006/access/production/resources/reference/versions/hg19/Homo_sapiens_assembly19.fasta}"
OUTPUT_DIR="${OUTPUT_DIR:-${BASE}/pon/${ASSAY}}"
OUTPUT="${OUTPUT_DIR}/${ASSAY}.${VARIANT}.pon.parquet"

# Per-sample feature outputs, kept.
#
# Without this they are extracted to a temp directory and deleted on success,
# so every rebuild re-runs extraction over every BAM -- hours -- and
# leave-one-out calibration costs n of those. They also survive a failed build,
# so a rerun skips what finished.
KEEP_DIR="${KEEP_DIR:-${OUTPUT_DIR}/${VARIANT}.sample_outputs}"

# Recorded in the PON alongside a salted digest of the sample IDs. Free text;
# no identifier is stored either way.
COHORT_LABEL="${COHORT_LABEL:-${ASSAY}-${VARIANT}-healthy-donors}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate "${KREWLYZER_ENV:-pygbcms}"

mkdir -p "${OUTPUT_DIR}"

if [ ! -f "${SAMPLE_LIST}" ]; then
    echo "Sample list not found: ${SAMPLE_LIST}" >&2
    echo "Set SAMPLE_LIST=... to override." >&2
    exit 2
fi

N_SAMPLES=$(grep -cv '^\s*$' "${SAMPLE_LIST}" || true)

echo "=== Krewlyzer PON build ==="
echo "Assay:     ${ASSAY}"
echo "Variant:   ${VARIANT}"
echo "Samples:   ${N_SAMPLES}   (${SAMPLE_LIST})"
echo "Cohort:    ${COHORT_LABEL}"
echo "Output:    ${OUTPUT}"
echo "Keeping:   ${KEEP_DIR}"
echo "Version:   $(krewlyzer --version)"
echo "Start:     $(date)"
echo "==========================="

krewlyzer build-pon \
    "${SAMPLE_LIST}" \
    --assay "${ASSAY}" \
    -r "${REFERENCE}" \
    -o "${OUTPUT}" \
    --cohort-label "${COHORT_LABEL}" \
    --keep-sample-outputs "${KEEP_DIR}" \
    --threads 50 \
    --parallel-samples 4 \
    --memory-per-sample 25 \
    --sample-timeout 0 \
    -v

echo ""
echo "=== Checking the model before anyone scores against it ==="
# Non-negotiable: this is the gate that fails on a fabricated baseline, and
# every PON shipped before 0.9.0 fails it. A build that produces a model this
# rejects has not succeeded, whatever build-pon's exit code said.
krewlyzer validate-pon "${OUTPUT}" --json-report "${OUTPUT%.parquet}.validation.json"

echo ""
echo "=== Complete ==="
echo "Output:    ${OUTPUT}"
echo "End:       $(date)"
