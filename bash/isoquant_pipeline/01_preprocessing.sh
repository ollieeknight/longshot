#!/bin/bash
set -euo pipefail

BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
REFERENCE_DIR="/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/pacbio/reference"

PACBIO_DIR="${BASE_DIR}/workspace/pacbio"
BAM_DIR="${BASE_DIR}/workspace/pacbio/bams"
LOG_DIR="${BASE_DIR}/workspace/pacbio/logs"

mkdir -p "${BAM_DIR}"
mkdir -p "${LOG_DIR}"

MAPPING_FILE="${BASE_DIR}/workspace/pacbio/sample_mapping.csv"
echo "sample_id,suffix" > "${MAPPING_FILE}"

IDX=1
for SAMPLE in HC01_NK HC02_NK HC03_NK; do
    IDX_STR=$(printf "%02d" $IDX)
    echo "${SAMPLE},_${IDX_STR}" >> "${MAPPING_FILE}"
    case $SAMPLE in
        HC01_NK) SAMPLE_DIR="1_A01"; MOVIE="m84094_260226_200655_s1"; MUX_BC="bcM0001" ;;
        HC02_NK) SAMPLE_DIR="1_B01"; MOVIE="m84094_260226_220949_s2"; MUX_BC="bcM0002" ;;
        HC03_NK) SAMPLE_DIR="1_C01"; MOVIE="m84094_260227_001254_s3"; MUX_BC="bcM0003" ;;
    esac

    INPUT_BAM="${BASE_DIR}/raw/${SAMPLE_DIR}/hifi_reads/${MOVIE}.hifi_reads.${MUX_BC}.bam"
    WORKSPACE_DIR="${PACBIO_DIR}/${SAMPLE}"
    mkdir -p "${WORKSPACE_DIR}"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${SAMPLE}_preprocess
#SBATCH --output=${LOG_DIR}/${SAMPLE}.log
#SBATCH --error=${LOG_DIR}/${SAMPLE}.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=48:00:00

set -euo pipefail

source /data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/miniforge3/etc/profile.d/conda.sh
conda activate pacbio

# 1. Segment, Lima, Tag, Refine, Correct, Dedup
skera split -j \${SLURM_CPUS_PER_TASK} ${INPUT_BAM} ${REFERENCE_DIR}/../adapters/mas16_primers.fasta ${WORKSPACE_DIR}/segmented.bam
lima --isoseq -j \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/segmented.bam ${REFERENCE_DIR}/../primers/10x_3kit_primers.fasta ${WORKSPACE_DIR}/fl.bam
isoseq tag --design T-12U-16B -j \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/fl.5p--3p.bam ${WORKSPACE_DIR}/flt.bam
isoseq refine --require-polya -j \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/flt.bam ${REFERENCE_DIR}/../primers/10x_3kit_primers.fasta ${WORKSPACE_DIR}/fltnc.bam
isoseq correct -j \${SLURM_CPUS_PER_TASK} --barcodes ${REFERENCE_DIR}/../barcodes/3M-3pgex-may-2023.REVCOMP.txt.gz ${WORKSPACE_DIR}/fltnc.bam ${WORKSPACE_DIR}/fltnc.corrected.bam
samtools sort -@ \${SLURM_CPUS_PER_TASK} -t CB ${WORKSPACE_DIR}/fltnc.corrected.bam -o ${WORKSPACE_DIR}/fltnc.corrected.sorted.bam
isoseq groupdedup -j \${SLURM_CPUS_PER_TASK} --keep-non-real-cells ${WORKSPACE_DIR}/fltnc.corrected.sorted.bam ${WORKSPACE_DIR}/dedup.bam

# 2. Align to genome
pbmm2 align --preset ISOSEQ --sort -j \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/dedup.bam ${REFERENCE_DIR}/human_GRCh38_no_alt_analysis_set.fasta ${WORKSPACE_DIR}/mapped.bam

# 3. Inject Donor Suffix into CB tag to prevent merging collisions
SUFFIX="_${IDX_STR}"
MAPPED_MOD="${BAM_DIR}/${SAMPLE}_mapped_mod.bam"

samtools view -h "${WORKSPACE_DIR}/mapped.bam" | awk -v OFS="\t" -v suf="\${SUFFIX}" '{
    if (\$0 ~ /^@/) { print }
    else {
        for(i=12; i<=NF; i++) {
            if (\$i ~ /^CB:Z:/) { \$i = \$i suf }
        }
        print \$0
    }
}' | samtools view -@ \${SLURM_CPUS_PER_TASK} -b - > "\${MAPPED_MOD}"

samtools index -@ \${SLURM_CPUS_PER_TASK} "\${MAPPED_MOD}"

EOF

    ((IDX++))
done