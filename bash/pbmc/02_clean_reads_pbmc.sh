#!/bin/bash
BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs"
REFERENCE_DIR="/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/pacbio"

PROJECT_ID="pacbio"

WORKING_DIR="${BASE_DIR}/${PROJECT_ID}"

mkdir -p "${WORKING_DIR}/workspace/logs"

for SAMPLE in PBMC1_s2 PBMC1_s3 PBMC2_s4 PBMC2_s1; do
    case $SAMPLE in
        PBMC1_s2)
            SAMPLE_DIR="1_PBMC1"
            MOVIE="m84039_250411_025511_s2"
            MUX_BC="bcM0002"
            ;;
        PBMC1_s3)
            SAMPLE_DIR="1_PBMC1"
            MOVIE="m84039_250411_045809_s3"
            MUX_BC="bcM0002"
            ;;
        PBMC2_s4)
            SAMPLE_DIR="1_PBMC2"
            MOVIE="m84039_250411_070102_s4"
            MUX_BC="bcM0003"
            ;;
        PBMC2_s1)
            SAMPLE_DIR="1_PBMC2"
            MOVIE="m84039_250412_044518_s1"
            MUX_BC="bcM0003"
            ;;
    esac

    INPUT_BAM="${WORKING_DIR}/raw/${SAMPLE_DIR}/${MOVIE}.hifi_reads.${MUX_BC}.bam"
    WORKSPACE_DIR="${WORKING_DIR}/workspace/${SAMPLE}"

    mkdir -p "${WORKSPACE_DIR}"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${SAMPLE}_clean
#SBATCH --output=${WORKING_DIR}/workspace/logs/${SAMPLE}_clean.log
#SBATCH --error=${WORKING_DIR}/workspace/logs/${SAMPLE}_clean.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=48:00:00

set -euo pipefail

source /data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/miniforge3/etc/profile.d/conda.sh
conda activate pacbio

# segment reads
skera split --num-threads \${SLURM_CPUS_PER_TASK} ${INPUT_BAM} ${REFERENCE_DIR}/adapters/mas16_primers.fasta ${WORKSPACE_DIR}/segmented.bam

# primer removal
lima  --isoseq --num-threads \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/segmented.bam ${REFERENCE_DIR}/primers/10x_3kit_primers.fasta ${WORKSPACE_DIR}/fl.bam

# tag UMI and cell barcode
isoseq tag --design T-12U-16B --num-threads \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/fl.5p--3p.bam ${WORKSPACE_DIR}/flt.bam

# trim polyA and remove concatemers
isoseq refine --require-polya --num-threads \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/flt.bam ${REFERENCE_DIR}/primers/10x_3kit_primers.fasta ${WORKSPACE_DIR}/fltnc.bam
EOF

done