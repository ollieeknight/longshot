#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs"
WORK_DIR="${BASE_DIR}/msk"
BARCODE_DIR="/data/cephfs-1/work/groups/romagnani/users/knighto_c/data/adaptive_nk/objects/longread/barcodes"
OUT_DIR="${WORK_DIR}/downstream/01_population_bams"
DONORS=(HC01 HC02 HC03)
NPROC=$(nproc)

command -v sinto samtools gawk >/dev/null 2>&1 || { echo "Missing sinto, samtools, or gawk"; exit 1; }
mkdir -p "${OUT_DIR}"

for donor in "${DONORS[@]}"; do
    bam="${WORK_DIR}/workspace/${donor}_NK/mapped.bam"
    barcodes="${BARCODE_DIR}/${donor}.tsv"
    class="${WORK_DIR}/workspace/${donor}_NK/collapsed_classification.txt"
    out="${OUT_DIR}/${donor}"

    [[ -f "$bam" && -f "$barcodes" && -f "$class" ]] || { echo "Missing BAM, barcodes, or classification for $donor"; exit 1; }
    
    mkdir -p "$out"
    [[ -f "${bam}.bai" ]] || samtools index -@ "$NPROC" "$bam"
    
    skip_file=$(mktemp)
    tail -n +2 "$class" | gawk -F'\t' '$5 ~ /^(RPL|RPS|MRPL|MRPS|MALAT1)/ {id=$1; sub(/^molecule\//, "", id); print id}' | sort -u > "$skip_file"
    
    if [[ -s "$skip_file" ]]; then
        filt_bam=$(mktemp --suffix=".bam")
        samtools view -h "$bam" | \
            gawk -F'\t' "NR==FNR {skip[\$0]; next} /^@/ {print; next} {match(\$0, /zm:i:([^\t]+)/, a); if (!(a[1] in skip)) print}" "$skip_file" - | \
            samtools view -b -o "$filt_bam" - || { rm -f "$skip_file" "$filt_bam"; exit 1; }
        samtools index -@ "$NPROC" "$filt_bam"
        bam="$filt_bam"
    fi
    rm -f "$skip_file"
    
    sinto filterbarcodes -b "$bam" -c "$barcodes" --barcodetag CB \
        --outdir "$out" --nproc "$NPROC"
    
    [[ -n "$(ls "$out"/*.bam 2>/dev/null)" ]] || { echo "No BAMs for $donor"; exit 1; }
    for popbam in "$out"/*.bam; do samtools index -@ "$NPROC" "$popbam"; done
    
    [[ "$bam" != "${WORK_DIR}/workspace/${donor}_NK/mapped.bam" ]] && rm -f "$bam" "$bam.bai"
done