#!/bin/bash
BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs"
PROJECT_ID="pacbio"
WORKING_DIR="${BASE_DIR}/${PROJECT_ID}"
RAW_DIR="${WORKING_DIR}/raw"

mkdir -p "${WORKING_DIR}/raw/logs"

for SAMPLE in PBMC1_s2 PBMC1_s3 PBMC2_s4 PBMC2_s1; do
    case $SAMPLE in
        PBMC1_s2)
            SAMPLE_DIR="1_PBMC1"
            MOVIE="m84039_250411_025511_s2"
            MUX_BC="bcM0002"
            BASE_URL="https://downloads.pacbcloud.com/public/dataset/MAS-Seq/DATA-Revio-Kinnex-PBMC-20kcells-10xGEMX3p-rep1/0-CCS"
            ;;
        PBMC1_s3)
            SAMPLE_DIR="1_PBMC1"
            MOVIE="m84039_250411_045809_s3"
            MUX_BC="bcM0002"
            BASE_URL="https://downloads.pacbcloud.com/public/dataset/MAS-Seq/DATA-Revio-Kinnex-PBMC-20kcells-10xGEMX3p-rep1/0-CCS"
            ;;
        PBMC2_s4)
            SAMPLE_DIR="1_PBMC2"
            MOVIE="m84039_250411_070102_s4"
            MUX_BC="bcM0003"
            BASE_URL="https://downloads.pacbcloud.com/public/dataset/MAS-Seq/DATA-Revio-Kinnex-PBMC-20kcells-10xGEMX3p-rep2/0-CCS"
            ;;
        PBMC2_s1)
            SAMPLE_DIR="1_PBMC2"
            MOVIE="m84039_250412_044518_s1"
            MUX_BC="bcM0003"
            BASE_URL="https://downloads.pacbcloud.com/public/dataset/MAS-Seq/DATA-Revio-Kinnex-PBMC-20kcells-10xGEMX3p-rep2/0-CCS"
            ;;
    esac

    BAM_FILE="${MOVIE}.hifi_reads.${MUX_BC}.bam"
    PBI_FILE="${BAM_FILE}.pbi"
    OUTPUT_DIR="${RAW_DIR}/${SAMPLE_DIR}"

    mkdir -p "${OUTPUT_DIR}"

    # Download BAM files
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=dl_${SAMPLE}_bam
#SBATCH --output=${WORKING_DIR}/raw/logs/download_${BAM_FILE}.log
#SBATCH --error=${WORKING_DIR}/raw/logs/download_${BAM_FILE}.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=24:00:00

set -euo pipefail
wget -c --progress=dot:giga -O ${OUTPUT_DIR}/${BAM_FILE} ${BASE_URL}/${BAM_FILE}
EOF

    # Download PBI index files
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=dl_${SAMPLE}_pbi
#SBATCH --output=${WORKING_DIR}/logs/download_${PBI_FILE}.log
#SBATCH --error=${WORKING_DIR}/logs/download_${PBI_FILE}.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=24:00:00

set -euo pipefail
wget -c --progress=dot:giga -O ${OUTPUT_DIR}/${PBI_FILE} ${BASE_URL}/${PBI_FILE}
EOF

done