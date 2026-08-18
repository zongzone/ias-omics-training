#!/bin/bash
set -Eeuo pipefail

# ==============================================================================
# PGGB: partition first, then run each community sequentially
#
# Usage:
#   bash 03_pggb.sh all.fa.gz outdir \
#       --threads 32 --poa-threads 16 --memory 256G \
#       --identity 98 --segment 10k --min-match 47 \
#       [--haps 61] [--vcf-ref REF]
#
# This script intentionally does NOT run whole-genome PGGB in one job.
# It uses the official partition-before-pggb workflow and then executes the
# generated community commands one by one.
#
# Note:
#   --memory is a documented target for the current node. PGGB does not expose a
#   single whole-pipeline RAM allocation flag.
# ==============================================================================

if [ "$#" -lt 2 ]; then
    echo "Usage: bash $0 <all_assemblies.fa[.gz]> <outdir> [options]"
    exit 1
fi

INPUT=$(readlink -f "$1")
OUTDIR="$2"
shift 2

THREADS=32
POA_THREADS=16
MEMORY="256G"
IDENTITY=98
SEGMENT="10k"
MIN_MATCH=47
HAPS=""
VCF_REF=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --poa-threads) POA_THREADS="$2"; shift 2 ;;
        --memory) MEMORY="$2"; shift 2 ;;
        --identity) IDENTITY="$2"; shift 2 ;;
        --segment) SEGMENT="$2"; shift 2 ;;
        --min-match) MIN_MATCH="$2"; shift 2 ;;
        --haps) HAPS="$2"; shift 2 ;;
        --vcf-ref) VCF_REF="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: bash $0 <all_assemblies.fa[.gz]> <outdir> --threads N --poa-threads N --memory 256G --identity 98 --segment 10k --min-match 47 [--haps N] [--vcf-ref REF]"
            exit 0 ;;
        *) echo "ERROR: unknown option: $1"; exit 1 ;;
    esac
done

for CMD in partition-before-pggb pggb samtools; do
    command -v "$CMD" >/dev/null 2>&1 || {
        echo "ERROR: $CMD not found"
        exit 1
    }
done

[ -f "$INPUT" ] || { echo "ERROR: input FASTA not found: $INPUT"; exit 1; }

mkdir -p "$OUTDIR"
OUTDIR=$(readlink -f "$OUTDIR")

if [ ! -s "${INPUT}.fai" ]; then
    samtools faidx "$INPUT"
fi

PARTITION_STDOUT="${OUTDIR}/partition.stdout.log"
COMMANDS="${OUTDIR}/pggb_commands.sh"

ARGS=(
    partition-before-pggb
    -i "$INPUT"
    -o "$OUTDIR"
    -t "$THREADS"
    --poa-threads "$POA_THREADS"
    -p "$IDENTITY"
    -s "$SEGMENT"
    -k "$MIN_MATCH"
    --skip-viz
    --resume
)

if [ -n "$HAPS" ]; then
    ARGS+=( -n "$HAPS" )
fi

if [ -n "$VCF_REF" ]; then
    ARGS+=( -V "${VCF_REF}:1000" )
fi

echo "============================================================"
echo "PGGB partition-first direct run"
echo "Input       : $INPUT"
echo "Threads     : $THREADS"
echo "POA threads : $POA_THREADS"
echo "Memory target for current node: $MEMORY"
echo "Identity    : $IDENTITY"
echo "Segment     : $SEGMENT"
echo "Min match   : $MIN_MATCH"
echo "Output      : $OUTDIR"
echo "============================================================"

# The official partition script prints the generated per-community PGGB commands
# to stdout. Capture them while retaining the complete stdout log.
"${ARGS[@]}" | tee "$PARTITION_STDOUT"

# Extract only the generated multi-line `pggb -i ...` command blocks.
awk '
    /^pggb -i / {capture=1}
    capture {print}
    capture && /--poa-threads/ {
        print ""
        capture=0
    }
' "$PARTITION_STDOUT" > "$COMMANDS"

if ! grep -q '^pggb -i ' "$COMMANDS"; then
    echo "ERROR: no generated PGGB community commands were found."
    echo "Check: $PARTITION_STDOUT"
    exit 1
fi

chmod +x "$COMMANDS"

N=$(grep -c '^pggb -i ' "$COMMANDS")

echo
echo "Generated PGGB community commands: $N"
echo "Commands file: $COMMANDS"
echo
echo "Running communities sequentially to control peak memory..."
echo

bash "$COMMANDS"

echo
echo "Finished all PGGB communities."
echo "Community outputs are under: $OUTDIR"
