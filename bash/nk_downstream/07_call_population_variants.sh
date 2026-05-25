#!/usr/bin/bash

WORK_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
OUT_DIR="${WORK_DIR}/downstream/07_population_variants"
LOG_DIR="${OUT_DIR}/logs"
DONORS=(HC01 HC02 HC03)

mkdir -p "$LOG_DIR"

for donor in "${DONORS[@]}"; do
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=clair3_${donor}
#SBATCH --output=${LOG_DIR}/%x_%A_%a.log
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=199G

set -euo pipefail
shopt -s nullglob

REF="/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/pacbio/reference/human_GRCh38_no_alt_analysis_set.fasta"
SPLIT_DIR="${WORK_DIR}/downstream/01_population_bams"
OUT_DIR="${WORK_DIR}/downstream/07_population_variants"
CLAIR3_SIF="/data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/clair3_latest.sif"
CLAIR3_MODEL="/opt/models/hifi"
donor_dir="\${SPLIT_DIR}/${donor}"
# Keep chrM/chrY out of the default call set to avoid Clair3 NaN PL edge cases.
CTGS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"

for bam in "\${donor_dir}"/*.bam; do
    [[ -f "\${bam}.bai" ]] || samtools index -@ "\$(nproc)" "\$bam"

    source /data/cephfs-1/home/users/knighto_c/work/bin/miniforge3/etc/profile.d/conda.sh
    conda activate sinto

    n_reads=\$(samtools view -c "\$bam")

    pop=\$(basename "\$bam" .bam)
    vcf_dir="\${OUT_DIR}/${donor}/clair3/\${pop}"
    mkdir -p "\$vcf_dir"

    run_clair3_cmd=(
        apptainer exec -B /data "\$CLAIR3_SIF" /opt/bin/run_clair3.sh
        --bam_fn="\$(realpath "\$bam")"
        --ref_fn="\$(realpath "\$REF")"
        --threads="\$(nproc)"
        --platform=hifi
        --model_path="\$CLAIR3_MODEL"
        --output="\$(realpath "\$vcf_dir")"
        --gvcf
        --ctg_name="\$CTGS"
        --no_phasing_for_fa
        --pileup_only
        --sample_name="${donor}_\${pop}"
    )

    "\${run_clair3_cmd[@]}" || {
        echo "[WARN] Clair3 failed for ${donor}_\${pop}; retrying in clean output dir"
        rm -rf "\${vcf_dir:?}"/*
        "\${run_clair3_cmd[@]}"
    }

    if [[ -f "\${vcf_dir}/merge_output.gvcf.gz" ]]; then
        mv "\${vcf_dir}/merge_output.gvcf.gz" "\${vcf_dir}/${donor}_\${pop}.g.vcf.gz"
        [[ -f "\${vcf_dir}/merge_output.gvcf.gz.tbi" ]] && \
            mv "\${vcf_dir}/merge_output.gvcf.gz.tbi" "\${vcf_dir}/${donor}_\${pop}.g.vcf.gz.tbi" || \
            bcftools index -t "\${vcf_dir}/${donor}_\${pop}.g.vcf.gz"
    fi
done
EOF
done