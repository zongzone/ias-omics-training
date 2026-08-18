#!/bin/bash
set -Eeuo pipefail

# ==============================================================================
# Minigraph-Cactus: large-genome workflow
#
# Usage:
#   bash 02_cactus_pangenome.sh \
#       seqfile.txt reference_sample outdir outname \
#       --threads 64 --memory 850G --workdir /path/to/scratch [--giraffe]
#
# This is one cactus-pangenome run. The script does NOT manually split
# genomes by chromosome. Cactus retains its default internal graphmap-split step,
# which reduces downstream memory use.
# ==============================================================================

if [ "$#" -lt 4 ]; then
    echo "Usage: bash $0 <seqfile.txt> <reference_sample> <outdir> <outname> --threads N --memory 850G [--workdir DIR] [--giraffe]"
    exit 1
fi

SEQFILE=$(readlink -f "$1")
REFERENCE="$2"
OUTDIR="$3"
OUTNAME="$4"
shift 4

THREADS=64
MEMORY="850G"
WORKDIR=""
GIRAFFE=0

while [ "$#" -gt 0 ]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --memory) MEMORY="$2"; shift 2 ;;
        --workdir) WORKDIR="$2"; shift 2 ;;
        --giraffe) GIRAFFE=1; shift ;;
        -h|--help)
            echo "Usage: bash $0 <seqfile.txt> <reference_sample> <outdir> <outname> --threads N --memory 850G [--workdir DIR] [--giraffe]"
            exit 0 ;;
        *) echo "ERROR: unknown option: $1"; exit 1 ;;
    esac
done

command -v cactus-pangenome >/dev/null 2>&1 || {
    echo "ERROR: cactus-pangenome not found"
    exit 1
}

[ -f "$SEQFILE" ] || { echo "ERROR: seqfile not found: $SEQFILE"; exit 1; }

if ! awk 'NF && $1 !~ /^#/ {print $1}' "$SEQFILE" | grep -Fxq "$REFERENCE"; then
    echo "ERROR: reference sample '$REFERENCE' is not present in $SEQFILE"
    exit 1
fi

mkdir -p "$OUTDIR"
OUTDIR=$(readlink -f "$OUTDIR")

JOBSTORE="${OUTDIR}.jobstore"
LOGFILE="${OUTDIR}/${OUTNAME}.cactus.log"

MG_CORES="${CACTUS_MG_CORES:-$THREADS}"
MAP_CORES="${CACTUS_MAP_CORES:-8}"
CONS_CORES="${CACTUS_CONS_CORES:-$THREADS}"

if [ "$THREADS" -gt 1 ]; then
    INDEX_DEFAULT=$((THREADS - 1))
else
    INDEX_DEFAULT=1
fi
INDEX_CORES="${CACTUS_INDEX_CORES:-$INDEX_DEFAULT}"

# Do not allow stage core settings to exceed the whole-run ceiling.
[ "$MG_CORES" -le "$THREADS" ] || MG_CORES="$THREADS"
[ "$MAP_CORES" -le "$THREADS" ] || MAP_CORES="$THREADS"
[ "$CONS_CORES" -le "$THREADS" ] || CONS_CORES="$THREADS"
[ "$INDEX_CORES" -le "$THREADS" ] || INDEX_CORES="$THREADS"

ARGS=(
    cactus-pangenome
    "$JOBSTORE"
    "$SEQFILE"
    --outDir "$OUTDIR"
    --outName "$OUTNAME"
    --reference "$REFERENCE"

    --gfa
    --gbz
    --vcf

    --maxCores "$THREADS"
    --maxMemory "$MEMORY"
    --doubleMem true
    --retryCount 5

    --mgCores "$MG_CORES"
    --mapCores "$MAP_CORES"
    --consCores "$CONS_CORES"
    --indexCores "$INDEX_CORES"

    --logFile "$LOGFILE"
)

# Optional stage-memory overrides. Leave unset unless a stage is known to be
# underestimated on the current dataset.
if [ -n "${CACTUS_MG_MEMORY:-}" ]; then
    ARGS+=(--mgMemory "$CACTUS_MG_MEMORY")
fi

if [ -n "${CACTUS_CONS_MEMORY:-}" ]; then
    ARGS+=(--consMemory "$CACTUS_CONS_MEMORY")
fi

if [ -n "${CACTUS_INDEX_MEMORY:-}" ]; then
    ARGS+=(--indexMemory "$CACTUS_INDEX_MEMORY")
fi

if [ -n "$WORKDIR" ]; then
    mkdir -p "$WORKDIR"
    WORKDIR=$(readlink -f "$WORKDIR")
    ARGS+=(--workDir "$WORKDIR")
fi

if [ "$GIRAFFE" -eq 1 ]; then
    ARGS+=(--giraffe)
fi

echo "============================================================"
echo "Minigraph-Cactus direct run"
echo "Reference  : $REFERENCE"
echo "Threads    : $THREADS"
echo "Max memory : $MEMORY"
echo "mg/map     : $MG_CORES / $MAP_CORES cores"
echo "cons/index : $CONS_CORES / $INDEX_CORES cores"
echo "Output     : $OUTDIR"
echo "JobStore   : $JOBSTORE"
echo "============================================================"

printf 'Command:'
printf ' %q' "${ARGS[@]}"
printf '\n'

exec "${ARGS[@]}"
