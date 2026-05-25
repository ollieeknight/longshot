#!/bin/bash
BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs"
REFERENCE_DIR="/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/pacbio/reference"

PROJECT_ID="pacbio"

WORKING_DIR="${BASE_DIR}/${PROJECT_ID}"

# PBMC1 merge
mkdir -p "${WORKING_DIR}/workspace/PBMC1"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=PBMC1_merge
#SBATCH --output=${WORKING_DIR}/workspace/logs/PBMC1_merge.log
#SBATCH --error=${WORKING_DIR}/workspace/logs/PBMC1_merge.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=48:00:00

set -euo pipefail

source /data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/miniforge3/etc/profile.d/conda.sh
conda activate pacbio

# Merge PBMC1_s2 and PBMC1_s3
samtools merge -@ \${SLURM_CPUS_PER_TASK} ${WORKING_DIR}/workspace/PBMC1/fltnc.bam \
    ${WORKING_DIR}/workspace/PBMC1_s2/fltnc.bam \
    ${WORKING_DIR}/workspace/PBMC1_s3/fltnc.bam

# Barcode correction
isoseq correct --num-threads \${SLURM_CPUS_PER_TASK} \
    --barcodes ${REFERENCE_DIR}/../barcodes/3M-3pgex-may-2023.REVCOMP.txt.gz \
    ${WORKING_DIR}/workspace/PBMC1/fltnc.bam \
    ${WORKING_DIR}/workspace/PBMC1/fltnc.corrected.bam

# Sort by cell barcode
samtools sort -@ \${SLURM_CPUS_PER_TASK} -t CB \
    ${WORKING_DIR}/workspace/PBMC1/fltnc.corrected.bam \
    -o ${WORKING_DIR}/workspace/PBMC1/fltnc.corrected.sorted.bam

# Deduplication
isoseq groupdedup --num-threads \${SLURM_CPUS_PER_TASK} --log-level INFO --keep-non-real-cells \
    ${WORKING_DIR}/workspace/PBMC1/fltnc.corrected.sorted.bam \
    ${WORKING_DIR}/workspace/PBMC1/dedup.bam

EOF

# PBMC2 merge
mkdir -p "${WORKING_DIR}/workspace/PBMC2"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=PBMC2_merge
#SBATCH --output=${WORKING_DIR}/workspace/logs/PBMC2_merge.log
#SBATCH --error=${WORKING_DIR}/workspace/logs/PBMC2_merge.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=48:00:00

set -euo pipefail

source /data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/miniforge3/etc/profile.d/conda.sh
conda activate pacbio

# Merge PBMC2_s4 and PBMC2_s1
samtools merge -@ \${SLURM_CPUS_PER_TASK} ${WORKING_DIR}/workspace/PBMC2/fltnc.bam \
    ${WORKING_DIR}/workspace/PBMC2_s4/fltnc.bam \
    ${WORKING_DIR}/workspace/PBMC2_s1/fltnc.bam

# Barcode correction
isoseq correct --num-threads \${SLURM_CPUS_PER_TASK} \
    --barcodes ${REFERENCE_DIR}/../barcodes/3M-3pgex-may-2023.REVCOMP.txt.gz \
    ${WORKING_DIR}/workspace/PBMC2/fltnc.bam \
    ${WORKING_DIR}/workspace/PBMC2/fltnc.corrected.bam

# Sort by cell barcode
samtools sort -@ \${SLURM_CPUS_PER_TASK} -t CB \
    ${WORKING_DIR}/workspace/PBMC2/fltnc.corrected.bam \
    -o ${WORKING_DIR}/workspace/PBMC2/fltnc.corrected.sorted.bam

# Deduplication
isoseq groupdedup --num-threads \${SLURM_CPUS_PER_TASK} --log-level INFO --keep-non-real-cells \
    ${WORKING_DIR}/workspace/PBMC2/fltnc.corrected.sorted.bam \
    ${WORKING_DIR}/workspace/PBMC2/dedup.bam

EOF