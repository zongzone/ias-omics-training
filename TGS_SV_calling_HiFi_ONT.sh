#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================================
# Third-generation sequencing (TGS) structural variant calling
# PacBio HiFi and Oxford Nanopore (ONT)
#
# Usage:
#   bash TGS_SV_calling_HiFi_ONT.sh hifi reference.fa reads.fastq.gz SAMPLE [threads] [outdir]
#   bash TGS_SV_calling_HiFi_ONT.sh ont  reference.fa reads.fastq.gz SAMPLE [threads] [outdir]
#
# Main software:
#   Common : samtools, cuteSV, SVIM, Sniffles2
#   HiFi   : pbmm2, pbsv
#   ONT    : minimap2
#   Optional ONT caller: NanoSV (requires a NanoSV config file)
#
# Default minimum SV size in this public template: 50 bp
# Override, for example: MIN_SV=30 bash TGS_SV_calling_HiFi_ONT.sh ...
# ============================================================================

usage() {
    cat <<'USAGE'
Usage:
  bash TGS_SV_calling_HiFi_ONT.sh hifi REF.fa READS.fastq.gz SAMPLE [THREADS] [OUTDIR]
  bash TGS_SV_calling_HiFi_ONT.sh ont  REF.fa READS.fastq.gz SAMPLE [THREADS] [OUTDIR]

Examples:
  bash TGS_SV_calling_HiFi_ONT.sh hifi reference.fa sample.hifi.fastq.gz cattle01 16
  bash TGS_SV_calling_HiFi_ONT.sh ont  reference.fa sample.ont.fastq.gz  cattle02 16

Optional NanoSV for ONT:
  export NANOSV_CONFIG=/path/to/config_nanosv.ini
  bash TGS_SV_calling_HiFi_ONT.sh ont reference.fa sample.ont.fastq.gz cattle02 16
USAGE
}

if [[ $# -lt 4 ]]; then
    usage
    exit 1
fi

MODE="${1,,}"
REF="$2"
READS="$3"
SAMPLE="$4"
THREADS="${5:-16}"
OUTDIR="${6:-${SAMPLE}_${MODE}_SV}"
MIN_SV="${MIN_SV:-50}"
MAX_SV="${MAX_SV:-100000}"
MIN_MAPQ="${MIN_MAPQ:-20}"

[[ "$MODE" == "hifi" || "$MODE" == "ont" ]] || {
    echo "ERROR: mode must be hifi or ont" >&2
    usage
    exit 1
}

[[ -f "$REF" ]] || { echo "ERROR: reference not found: $REF" >&2; exit 1; }
[[ -f "$READS" ]] || { echo "ERROR: reads not found: $READS" >&2; exit 1; }

mkdir -p "$OUTDIR"
REF=$(readlink -f "$REF")
READS=$(readlink -f "$READS")
OUTDIR=$(readlink -f "$OUTDIR")

# Reference FASTA index
samtools faidx "$REF"

# ============================================================================
# PacBio HiFi
#   Alignment: pbmm2 (PacBio-recommended aligner)
#   Callers  : pbsv + cuteSV + SVIM + Sniffles2
# ============================================================================
run_hifi() {
    local BAM="$OUTDIR/${SAMPLE}.sort.bam"

    echo "===== HiFi: pbmm2 alignment ====="
    pbmm2 align \
        "$REF" \
        "$READS" \
        "$BAM" \
        --sort \
        --preset CCS \
        --rg "@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}" \
        --num-threads "$THREADS"

    # pbmm2 normally creates a BAM index for sorted output; create one if absent.
    [[ -f "${BAM}.bai" ]] || samtools index -@ "$THREADS" "$BAM"

    echo "===== HiFi: pbsv ====="
    pbsv discover \
        "$BAM" \
        "$OUTDIR/${SAMPLE}.pbsv.svsig.gz"

    pbsv call \
        --ccs \
        "$REF" \
        "$OUTDIR/${SAMPLE}.pbsv.svsig.gz" \
        "$OUTDIR/${SAMPLE}.pbsv.vcf"

    echo "===== HiFi: cuteSV ====="
    mkdir -p "$OUTDIR/cutesv_tmp"
    cuteSV \
        "$BAM" \
        "$REF" \
        "$OUTDIR/${SAMPLE}.cutesv.vcf" \
        "$OUTDIR/cutesv_tmp" \
        --max_cluster_bias_INS 1000 \
        --diff_ratio_merging_INS 0.9 \
        --max_cluster_bias_DEL 1000 \
        --diff_ratio_merging_DEL 0.5 \
        --genotype \
        --sample "$SAMPLE" \
        --threads "$THREADS" \
        --min_mapq "$MIN_MAPQ" \
        --min_size "$MIN_SV" \
        --max_size "$MAX_SV"

    echo "===== HiFi: SVIM ====="
    rm -rf "$OUTDIR/svim_output"
    svim alignment \
        "$OUTDIR/svim_output" \
        "$BAM" \
        "$REF" \
        --min_mapq "$MIN_MAPQ" \
        --min_sv_size "$MIN_SV" \
        --max_sv_size "$MAX_SV" \
        --minimum_depth 3 \
        --sample "$SAMPLE"

    echo "===== HiFi: Sniffles2 ====="
    sniffles \
        --input "$BAM" \
        --vcf "$OUTDIR/${SAMPLE}.sniffles.vcf.gz" \
        --snf "$OUTDIR/${SAMPLE}.sniffles.snf" \
        --threads "$THREADS" \
        --reference "$REF"

    echo "===== HiFi calling finished ====="
}

# ============================================================================
# Oxford Nanopore (ONT)
#   Alignment: minimap2 -x map-ont
#   Callers  : cuteSV + SVIM + Sniffles2
#   Optional : NanoSV if NANOSV_CONFIG is provided
# ============================================================================
run_ont() {
    local BAM="$OUTDIR/${SAMPLE}.sort.bam"

    echo "===== ONT: minimap2 alignment ====="
    minimap2 \
        -ax map-ont \
        --MD \
        -R "@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}" \
        -t "$THREADS" \
        "$REF" \
        "$READS" \
        | samtools sort \
            -@ "$THREADS" \
            -o "$BAM"

    samtools index -@ "$THREADS" "$BAM"

    echo "===== ONT: cuteSV ====="
    mkdir -p "$OUTDIR/cutesv_tmp"
    cuteSV \
        "$BAM" \
        "$REF" \
        "$OUTDIR/${SAMPLE}.cutesv.vcf" \
        "$OUTDIR/cutesv_tmp" \
        --max_cluster_bias_INS 100 \
        --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 \
        --diff_ratio_merging_DEL 0.3 \
        --genotype \
        --sample "$SAMPLE" \
        --threads "$THREADS" \
        --min_mapq "$MIN_MAPQ" \
        --min_size "$MIN_SV" \
        --max_size "$MAX_SV"

    echo "===== ONT: SVIM ====="
    rm -rf "$OUTDIR/svim_output"
    svim alignment \
        "$OUTDIR/svim_output" \
        "$BAM" \
        "$REF" \
        --min_mapq "$MIN_MAPQ" \
        --min_sv_size "$MIN_SV" \
        --max_sv_size "$MAX_SV" \
        --minimum_depth 3 \
        --sample "$SAMPLE"

    echo "===== ONT: Sniffles2 ====="
    sniffles \
        --input "$BAM" \
        --vcf "$OUTDIR/${SAMPLE}.sniffles.vcf.gz" \
        --snf "$OUTDIR/${SAMPLE}.sniffles.snf" \
        --threads "$THREADS" \
        --reference "$REF"

    # NanoSV is kept as an optional caller because it requires a separate config file.
    if [[ -n "${NANOSV_CONFIG:-}" ]]; then
        if ! command -v NanoSV >/dev/null 2>&1; then
            echo "WARNING: NANOSV_CONFIG was provided, but NanoSV is not in PATH. Skipping NanoSV." >&2
        elif [[ ! -f "$NANOSV_CONFIG" ]]; then
            echo "WARNING: NanoSV config not found: $NANOSV_CONFIG. Skipping NanoSV." >&2
        else
            echo "===== ONT: NanoSV (optional) ====="
            NanoSV \
                -t "$THREADS" \
                -c "$NANOSV_CONFIG" \
                "$BAM" \
                -o "$OUTDIR/${SAMPLE}.NanoSV.vcf" \
                > "$OUTDIR/${SAMPLE}.NanoSV.log" 2>&1
        fi
    fi

    echo "===== ONT calling finished ====="
}

case "$MODE" in
    hifi) run_hifi ;;
    ont)  run_ont ;;
esac

echo
echo "Output directory: $OUTDIR"
echo "Sample: $SAMPLE"
echo "Mode: $MODE"
echo "Minimum SV size: ${MIN_SV} bp"
echo "Minimum mapping quality: ${MIN_MAPQ}"
