#!/bin/bash
BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs"
REFERENCE_DIR="/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/pacbio/reference"

PROJECT_ID="pacbio"

WORKING_DIR="${BASE_DIR}/${PROJECT_ID}"

for SAMPLE in PBMC1 PBMC2; do
    WORKSPACE_DIR="${WORKING_DIR}/workspace/${SAMPLE}"
    OUTS_DIR="${WORKING_DIR}/outs/${SAMPLE}"

    mkdir -p "${WORKSPACE_DIR}"
    mkdir -p "${OUTS_DIR}"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${SAMPLE}_pigeon
#SBATCH --output=${WORKING_DIR}/workspace/logs/${SAMPLE}_pigeon.log
#SBATCH --error=${WORKING_DIR}/workspace/logs/${SAMPLE}_pigeon.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=48:00:00

set -euo pipefail

source /data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/miniforge3/etc/profile.d/conda.sh
conda activate pacbio

# Change to workspace directory for pigeon
cd ${WORKSPACE_DIR}

# Map to reference genome
pbmm2 align --preset ISOSEQ --sort -j \${SLURM_CPUS_PER_TASK} dedup.bam ${REFERENCE_DIR}/human_GRCh38_no_alt_analysis_set.fasta mapped.bam

# Collapse into unique isoforms
isoseq collapse --num-threads \${SLURM_CPUS_PER_TASK} mapped.bam collapsed.gff

# Prepare collapsed GFF
pigeon prepare collapsed.gff

# Classify isoforms
pigeon classify \
    --fl collapsed.abundance.txt \
    --cage-peak ${REFERENCE_DIR}/refTSS_v3.3_human_coordinate.hg38.sorted.bed \
    --poly-a ${REFERENCE_DIR}/polyA.list.txt \
    -j \${SLURM_CPUS_PER_TASK} \
    collapsed.sorted.gff \
    ${REFERENCE_DIR}/gencode.v39.annotation.sorted.gtf \
    ${REFERENCE_DIR}/human_GRCh38_no_alt_analysis_set.fasta \
    collapsed

# Filter isoforms
pigeon filter -j \${SLURM_CPUS_PER_TASK} collapsed_classification.txt --isoforms collapsed.sorted.gff

# Saturation report
pigeon report -j \${SLURM_CPUS_PER_TASK} collapsed_classification.filtered_lite_classification.txt saturation.txt

# Seurat-compatible count matrices
pigeon make-seurat \
    -j \${SLURM_CPUS_PER_TASK} \
    --dedup dedup.fasta \
    --group collapsed.group.txt \
    --keep-ribo-mito-genes \
    -d ${OUTS_DIR} \
    collapsed_classification.filtered_lite_classification.txt

EOF

done