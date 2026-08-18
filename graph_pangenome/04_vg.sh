#!/bin/bash
set -Eeuo pipefail

# ==============================================================================
# VG: whole-genome graph construction from reference + VCF
#
# Usage:
#   bash 04_vg_from_vcf_direct.sh reference.fa variants.vcf.gz graph.vg \
#       --threads 32 --memory 256G [--gfa] [--no-handle-sv]
#
# This script intentionally DOES NOT split by chromosome.
#
# Note:
#   --memory is a documented target for the current node. `vg construct` has a
#   thread option but no whole-run RAM reservation option.
# ==============================================================================

if [ "$#" -lt 3 ]; then
    echo "Usage: bash $0 <reference.fa> <variants.vcf.gz> <output.vg> [--threads N] [--memory 256G] [--gfa] [--no-handle-sv]"
    exit 1
fi

REF=$(readlink -f "$1")
VCF=$(readlink -f "$2")
OUT="$3"
shift 3

THREADS=32
MEMORY="256G"
MAKE_GFA=0
HANDLE_SV=1

while [ "$#" -gt 0 ]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --memory) MEMORY="$2"; shift 2 ;;
        --gfa) MAKE_GFA=1; shift ;;
        --no-handle-sv) HANDLE_SV=0; shift ;;
        -h|--help)
            echo "Usage: bash $0 <reference.fa> <variants.vcf.gz> <output.vg> [--threads N] [--memory 256G] [--gfa] [--no-handle-sv]"
            exit 0 ;;
        *) echo "ERROR: unknown option: $1"; exit 1 ;;
    esac
done

for CMD in vg samtools tabix; do
    command -v "$CMD" >/dev/null 2>&1 || {
        echo "ERROR: $CMD not found"
        exit 1
    }
done

[ -f "$REF" ] || { echo "ERROR: reference not found: $REF"; exit 1; }
[ -f "$VCF" ] || { echo "ERROR: VCF not found: $VCF"; exit 1; }

case "$VCF" in
    *.vcf.gz|*.vcf.bgz) ;;
    *)
        echo "ERROR: input VCF must be bgzip-compressed (.vcf.gz or .vcf.bgz)"
        exit 1
        ;;
esac

if [ ! -s "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

if [ ! -s "${VCF}.tbi" ] && [ ! -s "${VCF}.csi" ]; then
    tabix -p vcf "$VCF"
fi

ARGS=(
    vg construct
    -r "$REF"
    -v "$VCF"
    -t "$THREADS"
)

if [ "$HANDLE_SV" -eq 1 ]; then
    ARGS+=( -S )
fi

echo "============================================================"
echo "VG whole-genome construction"
echo "Reference  : $REF"
echo "VCF        : $VCF"
echo "Threads    : $THREADS"
echo "Memory target for current node: $MEMORY"
echo "Handle SV  : $HANDLE_SV"
echo "Output     : $OUT"
echo "============================================================"

"${ARGS[@]}" > "$OUT"

if [ ! -s "$OUT" ]; then
    echo "ERROR: VG output is empty: $OUT"
    exit 1
fi

echo "Finished VG graph: $OUT"

if [ "$MAKE_GFA" -eq 1 ]; then
    GFA="${OUT%.*}.gfa"
    vg convert -f -t "$THREADS" "$OUT" > "$GFA"
    echo "GFA output: $GFA"
fi
