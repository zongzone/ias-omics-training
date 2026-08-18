#!/bin/bash
set -Eeuo pipefail

# ==============================================================================
# Minigraph: whole-genome direct construction
#
# Usage:
#   bash 01_minigraph_direct.sh reference.fa assemblies.list output.gfa \
#       --threads 24 --memory 300G
#
# Note:
#   --memory is a documented target for the current node. Minigraph has no
#   whole-run memory-allocation option; Bash itself cannot reserve cluster RAM.
# ==============================================================================

if [ "$#" -lt 3 ]; then
    echo "Usage: bash $0 <reference.fa> <assemblies.list> <output.gfa> [--threads N] [--memory 300G]"
    exit 1
fi

REF=$(readlink -f "$1")
LIST=$(readlink -f "$2")
OUT="$3"
shift 3

THREADS=24
MEMORY="300G"

while [ "$#" -gt 0 ]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --memory)  MEMORY="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: bash $0 <reference.fa> <assemblies.list> <output.gfa> [--threads N] [--memory 300G]"
            exit 0 ;;
        *) echo "ERROR: unknown option: $1"; exit 1 ;;
    esac
done

command -v minigraph >/dev/null 2>&1 || { echo "ERROR: minigraph not found"; exit 1; }
[ -f "$REF" ] || { echo "ERROR: reference not found: $REF"; exit 1; }
[ -f "$LIST" ] || { echo "ERROR: assembly list not found: $LIST"; exit 1; }

ASSEMBLIES=()
while IFS= read -r FA; do
    [ -z "$FA" ] && continue
    [[ "$FA" =~ ^[[:space:]]*# ]] && continue
    ABS=$(readlink -f "$FA")
    [ -f "$ABS" ] || { echo "ERROR: assembly not found: $FA"; exit 1; }
    ASSEMBLIES+=("$ABS")
done < "$LIST"

[ "${#ASSEMBLIES[@]}" -gt 0 ] || { echo "ERROR: no assemblies found in $LIST"; exit 1; }

echo "============================================================"
echo "Minigraph whole-genome construction"
echo "Reference : $REF"
echo "Assemblies: ${#ASSEMBLIES[@]}"
echo "Threads   : $THREADS"
echo "Memory target for current node: $MEMORY"
echo "Output    : $OUT"
echo "============================================================"

minigraph \
    -cxggs \
    -t "$THREADS" \
    "$REF" \
    "${ASSEMBLIES[@]}" \
    > "$OUT"

echo "Finished: $OUT"
