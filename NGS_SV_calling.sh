#!/bin/bash
set -Eeuo pipefail

# ==============================================================================
# Short-read WGS structural-variant calling
#
# Callers:
#   1. Manta
#   2. smoove
#   3. DELLY
#   4. dysgu
#
# Input:
#   - reference FASTA
#   - coordinate-sorted, duplicate-marked and indexed BAM
#
# Usage:
#   bash Short_read_SV_calling.sh reference.fa sample.sorted.dedup.bam SAMPLE_ID 8
#
# Conda environments used by this script:
#   sv_py27   : legacy environment for Manta + smoove
#   sv_py310  : modern environment for DELLY + dysgu
#
# Example environment setup:
#
#   mamba create -n sv_py27 -c conda-forge -c bioconda \
#       python=2.7 manta=1.6.0 smoove=0.2.8 samtools bcftools htslib -y
#
#   mamba create -n sv_py310 -c conda-forge -c bioconda \
#       python=3.10 delly dysgu samtools bcftools -y
#
# Notes:
#   - Manta 1.6.0 is kept in the Python 2.7 environment.
#   - smoove is also run in the same legacy environment here to match the
#     tested workflow and its legacy dependency stack.
#   - Current dysgu releases require a modern Python 3 environment; Python 3.10
#     is recommended here rather than an unspecified "Python 3.x".
# ==============================================================================

if [ "$#" -lt 3 ]; then
    echo "Usage: bash $0 <reference.fa> <sample.sorted.dedup.bam> <sample_id> [threads]"
    exit 1
fi

REF=$(readlink -f "$1")
BAM=$(readlink -f "$2")
SAMPLE="$3"
THREADS="${4:-8}"

PY27_ENV="${PY27_ENV:-sv_py27}"
PY3_ENV="${PY3_ENV:-sv_py310}"

OUTDIR=$(readlink -f "${SAMPLE}.shortread_SV")
MANTA_DIR="${OUTDIR}/Manta"
SMOOVE_DIR="${OUTDIR}/Smoove"
DELLY_DIR="${OUTDIR}/Delly"
DYSGU_DIR="${OUTDIR}/Dysgu"

mkdir -p "${OUTDIR}" "${MANTA_DIR}" "${SMOOVE_DIR}" "${DELLY_DIR}" "${DYSGU_DIR}"

# ------------------------------------------------------------------------------
# Load conda
# ------------------------------------------------------------------------------

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda was not found in PATH."
    echo "Please initialize conda before running this script."
    exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"

# ------------------------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------------------------

[ -f "${REF}" ] || { echo "ERROR: reference not found: ${REF}"; exit 1; }
[ -f "${BAM}" ] || { echo "ERROR: BAM not found: ${BAM}"; exit 1; }

# Reference FASTA index
if [ ! -f "${REF}.fai" ]; then
    conda activate "${PY3_ENV}"
    samtools faidx "${REF}"
fi

# BAM index
if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM%.bam}.bai" ]; then
    conda activate "${PY3_ENV}"
    samtools index -@ "${THREADS}" "${BAM}"
fi

echo "============================================================"
echo " Short-read SV calling"
echo " Sample   : ${SAMPLE}"
echo " Reference: ${REF}"
echo " BAM      : ${BAM}"
echo " Threads  : ${THREADS}"
echo " Output   : ${OUTDIR}"
echo "============================================================"


# ==============================================================================
# 1. Manta
# ==============================================================================

echo
echo "===== [1/4] Manta ====="

conda activate "${PY27_ENV}"

MANTA_RUN="${MANTA_DIR}/run"

if [ ! -s "${MANTA_DIR}/Manta.candidateSV.vcf.gz" ]; then

    rm -rf "${MANTA_RUN}"

    configManta.py \
        --bam "${BAM}" \
        --referenceFasta "${REF}" \
        --runDir "${MANTA_RUN}"

    python "${MANTA_RUN}/runWorkflow.py" \
        -m local \
        -j "${THREADS}"

    cp "${MANTA_RUN}/results/variants/candidateSV.vcf.gz" \
       "${MANTA_DIR}/Manta.candidateSV.vcf.gz"

    if [ -f "${MANTA_RUN}/results/variants/candidateSV.vcf.gz.tbi" ]; then
        cp "${MANTA_RUN}/results/variants/candidateSV.vcf.gz.tbi" \
           "${MANTA_DIR}/Manta.candidateSV.vcf.gz.tbi"
    fi

else
    echo "Manta output exists; skipping."
fi


# ==============================================================================
# 2. smoove
# ==============================================================================

echo
echo "===== [2/4] smoove ====="

conda activate "${PY27_ENV}"

if [ ! -s "${SMOOVE_DIR}/${SAMPLE}-smoove.genotyped.vcf.gz" ]; then

    # smoove and several of its dependencies use temporary files.
    # Keep TMPDIR inside the analysis directory instead of system /tmp.
    SMOOVE_TMP="${SMOOVE_DIR}/tmp"
    mkdir -p "${SMOOVE_TMP}"

    export TMPDIR="${SMOOVE_TMP}"
    export TMP="${SMOOVE_TMP}"
    export TEMP="${SMOOVE_TMP}"

    smoove call \
        --outdir "${SMOOVE_DIR}" \
        --name "${SAMPLE}" \
        --fasta "${REF}" \
        -p "${THREADS}" \
        --genotype \
        "${BAM}"

else
    echo "smoove output exists; skipping."
fi


# ==============================================================================
# 3. DELLY
# ==============================================================================

echo
echo "===== [3/4] DELLY ====="

conda activate "${PY3_ENV}"

cd "${DELLY_DIR}"

if [ ! -s "Delly.vcf" ]; then

    # Keep the type-specific DELLY workflow used in the original analysis.
    # SV types: deletion, duplication, inversion, translocation and insertion.
    TYPES=("DEL" "DUP" "INV" "TRA" "INS")

    rm -f delly.*.bcf delly.*.vcf Delly.vcf \
          Delly.header.txt Delly.body.txt Delly.body.sorted.txt

    FIRST_VCF=""

    for TYPE in "${TYPES[@]}"; do

        echo "DELLY calling: ${TYPE}"

        delly call \
            -t "${TYPE}" \
            -o "delly.${TYPE}.bcf" \
            -g "${REF}" \
            "${BAM}"

        bcftools view \
            "delly.${TYPE}.bcf" \
            > "delly.${TYPE}.vcf"

        if [ -z "${FIRST_VCF}" ] && [ -s "delly.${TYPE}.vcf" ]; then
            FIRST_VCF="delly.${TYPE}.vcf"
        fi
    done

    [ -n "${FIRST_VCF}" ] || {
        echo "ERROR: DELLY produced no VCF output."
        exit 1
    }

    grep '^#' "${FIRST_VCF}" > Delly.header.txt

    : > Delly.body.txt
    for TYPE in "${TYPES[@]}"; do
        if [ -s "delly.${TYPE}.vcf" ]; then
            grep -v '^#' "delly.${TYPE}.vcf" >> Delly.body.txt
        fi
    done

    LC_ALL=C sort -k1,1V -k2,2n \
        Delly.body.txt \
        > Delly.body.sorted.txt

    cat Delly.header.txt Delly.body.sorted.txt > Delly.vcf

else
    echo "DELLY output exists; skipping."
fi

cd "${OUTDIR}"


# ==============================================================================
# 4. dysgu
# ==============================================================================

echo
echo "===== [4/4] dysgu ====="

conda activate "${PY3_ENV}"

DYSGU_TMP="${DYSGU_DIR}/tmp"
mkdir -p "${DYSGU_TMP}"

if [ ! -s "${DYSGU_DIR}/Dysgu.vcf" ]; then

    dysgu run \
        -x \
        -p "${THREADS}" \
        "${REF}" \
        "${DYSGU_TMP}" \
        "${BAM}" \
        > "${DYSGU_DIR}/Dysgu.vcf"

else
    echo "dysgu output exists; skipping."
fi


# ==============================================================================
# Summary
# ==============================================================================

echo
echo "============================================================"
echo " Finished: ${SAMPLE}"
echo "============================================================"
echo "Manta : ${MANTA_DIR}/Manta.candidateSV.vcf.gz"
echo "smoove: ${SMOOVE_DIR}/${SAMPLE}-smoove.genotyped.vcf.gz"
echo "DELLY : ${DELLY_DIR}/Delly.vcf"
echo "dysgu : ${DYSGU_DIR}/Dysgu.vcf"
echo "============================================================"
