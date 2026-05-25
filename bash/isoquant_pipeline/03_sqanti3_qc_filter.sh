#!/bin/bash
set -euo pipefail

BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
WORK_DIR="${BASE_DIR}/workspace/sqanti3"
ISOQUANT_DIR="${BASE_DIR}/workspace/isoquant"

mkdir -p "${WORK_DIR}/logs"

sbatch << "EOF"
#!/bin/bash
#SBATCH --job-name=sqanti3_joint
#SBATCH --output=/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk/workspace/sqanti3/logs/sqanti3_joint.%j.log
#SBATCH --error=/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk/workspace/sqanti3/logs/sqanti3_joint.%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=600G
#SBATCH --time=24:00:00

set -euo pipefail

BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
REFERENCE_DIR="/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/pacbio/reference"
WORK_DIR="${BASE_DIR}/workspace/sqanti3"
ISOQUANT_DIR="${BASE_DIR}/workspace/isoquant"
SQANTI_SIF="/data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/sqanti3_latest.sif"

mkdir -p "${WORK_DIR}/sqanti_qc"
mkdir -p "${WORK_DIR}/sqanti_filter"

cd "${WORK_DIR}"

ISOQUANT_GTF="${ISOQUANT_DIR}/isoquant_out/joint/joint.transcript_models.gtf"

apptainer run -B /data "${SQANTI_SIF}" sqanti3_qc.py \
    --isoforms "${ISOQUANT_GTF}" \
    --refGTF "${REFERENCE_DIR}/gencode.v39.annotation.sorted.gtf" \
    --refFasta "${REFERENCE_DIR}/human_GRCh38_no_alt_analysis_set.fasta" \
    --CAGE_peak "${REFERENCE_DIR}/refTSS_v3.3_human_coordinate.hg38.sorted.bed" \
    --polyA_motif_list "${REFERENCE_DIR}/polyA.list.txt" \
    --coverage "${REFERENCE_DIR}/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv" \
    -t ${SLURM_CPUS_PER_TASK} \
    -d "${WORK_DIR}/sqanti_qc" \
    -o joint \
    --report pdf

apptainer run -B /data "${SQANTI_SIF}" sqanti3_filter.py rules \
    --sqanti_class "${WORK_DIR}/sqanti_qc/joint_classification.txt" \
    --filter_isoforms "${WORK_DIR}/sqanti_qc/joint_corrected.fasta" \
    --filter_gtf "${WORK_DIR}/sqanti_qc/joint_corrected.gtf" \
    -d "${WORK_DIR}/sqanti_filter" \
    -o joint_filtered \
    --skip_report
EOF