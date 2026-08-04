#!/bin/bash
#SBATCH --job-name=krewlyzer_runall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --partition=cmobic_cpu
#SBATCH --output=runall_%x_%A_%a.log
#SBATCH --error=runall_%x_%A_%a.err

# ==============================================================================
# Extract features for one PON cohort, one SLURM task per sample.
#
#   SAMPLE_LIST=/path/to/list.txt sbatch --array=1-47%12 \
#       scripts/build_pon_array.sh xs1 duplex
#
# Then aggregate them -- minutes, on a login node, no BAM re-read:
#
#   krewlyzer build-pon --from-outputs <OUTPUT_DIR>/<variant>.sample_outputs \
#       --assay xs1 -r <REFERENCE> -o <OUTPUT_DIR>/xs1.duplex.pon.parquet
#   krewlyzer validate-pon <OUTPUT_DIR>/xs1.duplex.pon.parquet
#
# Why this exists
# ---------------
# `build_pon.sh` does extraction and aggregation in one process, four samples
# at a time on a single node: 55-97 min per sample, ~15 h for 47. The work is
# embarrassingly parallel across samples and the cluster already has the
# scheduler for it.
#
#   in-process, 4 parallel (build_pon.sh)   ~15.5 h
#   in-process, --parallel-samples 0        ~6.2 h
#   this array, 12 concurrent               ~5.9 h
#   this array, unthrottled                 ~1.5 h
#
# Nextflow would do the same fan-out, but every module pins
# `krewlyzer:0.8.3` and the SLURM profile enables singularity -- so a run today
# would silently execute 0.8.3 and reproduce the defects 0.9.0 removes. That
# needs a 0.9.0 container, which needs the release, which needs this PON.
# An sbatch array needs no container and no new infrastructure.
#
# `--from-outputs` will accept Nextflow's output unchanged once it exists; it
# is the same per-sample layout.
#
# The array index picks the line
# ------------------------------
# `--array=1-N` maps task k to line k of SAMPLE_LIST, so N must match the line
# count. The guard below refuses an index past the end rather than processing
# an empty path, which is how a short array silently builds a smaller cohort.
# ==============================================================================

set -euo pipefail

ASSAY="${1:?usage: build_pon_array.sh <assay: xs1|xs2> <variant: all_unique|duplex>}"
VARIANT="${2:?usage: build_pon_array.sh <assay: xs1|xs2> <variant: all_unique|duplex>}"

case "$ASSAY" in xs1|xs2) ;; *) echo "unknown assay: $ASSAY" >&2; exit 2 ;; esac
case "$VARIANT" in all_unique|duplex) ;; *) echo "unknown variant: $VARIANT" >&2; exit 2 ;; esac

# Required, never derived -- the same rule as build_pon.sh. Composing this name
# from the variant was guessed twice and wrong both times; the naming belongs
# to whoever assembled the cohort.
if [ -z "${SAMPLE_LIST:-}" ]; then
    echo "usage: SAMPLE_LIST=<file> sbatch --array=1-<N> build_pon_array.sh <assay> <variant>" >&2
    echo "" >&2
    echo "  SAMPLE_LIST must name the BAM list for this assay and variant." >&2
    echo "  <N> must equal its line count." >&2
    exit 2
fi

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "FATAL: no SLURM_ARRAY_TASK_ID. Submit with --array=1-<N>, e.g." >&2
    echo "       sbatch --array=1-\$(grep -cv '^[[:space:]]*\$' \"\$SAMPLE_LIST\") ..." >&2
    exit 2
fi

REFERENCE="${REFERENCE:?set REFERENCE to the genome FASTA}"
OUTPUT_DIR="${OUTPUT_DIR:-$(pwd)/pon/${ASSAY}}"
SAMPLE_OUTPUTS="${SAMPLE_OUTPUTS:-${OUTPUT_DIR}/${VARIANT}.sample_outputs}"

if [ ! -f "${SAMPLE_LIST}" ]; then
    echo "Sample list not found: ${SAMPLE_LIST}" >&2
    exit 2
fi

# Blank and commented lines are dropped first, so the index means what it says.
N_SAMPLES=$(grep -cv '^[[:space:]]*\(#.*\)\?$' "${SAMPLE_LIST}")
if [ "${SLURM_ARRAY_TASK_ID}" -gt "${N_SAMPLES}" ]; then
    echo "FATAL: task ${SLURM_ARRAY_TASK_ID} is past the end of the list" >&2
    echo "       (${N_SAMPLES} samples). Submit --array=1-${N_SAMPLES}." >&2
    exit 2
fi

BAM=$(grep -v '^[[:space:]]*\(#.*\)\?$' "${SAMPLE_LIST}" | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ ! -f "${BAM}" ]; then
    echo "FATAL: sample ${SLURM_ARRAY_TASK_ID} does not exist: ${BAM}" >&2
    exit 2
fi

SAMPLE=$(basename "${BAM}")
SAMPLE="${SAMPLE%.bam}"
SAMPLE="${SAMPLE%.cram}"
SAMPLE_DIR="${SAMPLE_OUTPUTS}/${SAMPLE}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate "${KREWLYZER_ENV:-krewlyzer}"

# Capability, not version. `krewlyzer --version` reports the previous release
# on develop and will until the bump, so a version comparison rejects exactly
# the build we want. See build_pon.sh for the 18-minute run that preceded this.
if ! krewlyzer validate-pon --help >/dev/null 2>&1; then
    echo "FATAL: this krewlyzer has no 'validate-pon', so it predates the PON" >&2
    echo "       rebuild work. Install from a develop checkout:" >&2
    echo "       maturin develop --release" >&2
    exit 3
fi

# Skip what is already finished.
#
# The marker is the *last* output a sample writes, matching what
# `--from-outputs` checks. Anything earlier would call a half-written
# directory complete and let the aggregation read it.
if compgen -G "${SAMPLE_DIR}/${SAMPLE}.MDS.gene.*" >/dev/null 2>&1 \
   && [ -z "${FORCE_REPROCESS:-}" ]; then
    echo "${SAMPLE}: already complete, skipping (set FORCE_REPROCESS=1 to redo)"
    exit 0
fi

mkdir -p "${SAMPLE_DIR}"

echo "=== krewlyzer run-all: task ${SLURM_ARRAY_TASK_ID}/${N_SAMPLES} ==="
echo "Sample:    ${SAMPLE}"
echo "Assay:     ${ASSAY} / ${VARIANT}"
echo "Output:    ${SAMPLE_DIR}"
echo "Version:   $(krewlyzer --version 2>&1 | head -1)"
echo "Start:     $(date)"
echo "==============================================="

# --skip-pon: this cohort is what a PON will be built *from*, so scoring it
# against an existing one would be both circular and pointless.
krewlyzer run-all \
    -i "${BAM}" \
    -r "${REFERENCE}" \
    -o "${SAMPLE_DIR}" \
    --assay "${ASSAY}" \
    --skip-pon \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

# Fail the task if the marker is missing, rather than leaving a directory the
# aggregation will refuse later with no clue which step went wrong.
if ! compgen -G "${SAMPLE_DIR}/${SAMPLE}.MDS.gene.*" >/dev/null 2>&1; then
    echo "FATAL: ${SAMPLE} finished but wrote no MDS.gene table, so the run is" >&2
    echo "       incomplete and --from-outputs will refuse it. Check the" >&2
    echo "       region-mds warnings above." >&2
    exit 4
fi

echo ""
echo "=== Complete: ${SAMPLE} ==="
echo "End:       $(date)"
