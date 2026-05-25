#!/bin/bash
set -euo pipefail

BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
WORK_DIR="${BASE_DIR}/workspace/isoquant"
BAM_DIR="${BASE_DIR}/workspace/pacbio/bams"

mkdir -p "${WORK_DIR}/logs"

sbatch << "EOF"
#!/bin/bash
#SBATCH --job-name=2_isoquant_joint
#SBATCH --output=/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk/workspace/isoquant/logs/isoquant_joint.%j.log
#SBATCH --error=/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk/workspace/isoquant/logs/isoquant_joint.%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=600G
#SBATCH --time=48:00:00

set -euo pipefail

BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
REFERENCE_DIR="/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/pacbio/reference"
WORK_DIR="${BASE_DIR}/workspace/isoquant"
BAM_DIR="${BASE_DIR}/workspace/pacbio/bams"

source /data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/miniforge3/etc/profile.d/conda.sh
conda activate pacbio

cd ${WORK_DIR}

INPUT_BAMS="${BAM_DIR}/HC01_NK_mapped_mod.bam ${BAM_DIR}/HC02_NK_mapped_mod.bam ${BAM_DIR}/HC03_NK_mapped_mod.bam"

chrs=$(echo chr{1..22} chrX chrY chrM)
# Run IsoQuant in Joint Single-Cell Mode (using BAM tags)
isoquant \
    --reference ${REFERENCE_DIR}/human_GRCh38_no_alt_analysis_set.fasta \
    --genedb ${REFERENCE_DIR}/gencode.v39.annotation.sorted.gtf \
    --complete_genedb \
    --process_only_chr ${chrs} \
    --bam ${INPUT_BAMS} \
    --data_type pacbio \
    --barcoded_bam \
    --barcode_tag CB \
    --umi_tag XM \
    --read_group barcode \
    --threads ${SLURM_CPUS_PER_TASK} \
    --prefix joint \
    --output ${WORK_DIR}/isoquant_out
EOF