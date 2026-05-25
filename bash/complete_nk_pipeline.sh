#!/bin/bash
BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs"
REFERENCE_DIR="/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/pacbio/reference"

PROJECT_ID="msk"

WORKING_DIR="${BASE_DIR}/${PROJECT_ID}"

mkdir -p "${WORKING_DIR}/workspace/logs"

for SAMPLE in HC01_NK HC02_NK HC03_NK; do
    case $SAMPLE in
        HC01_NK)
            SAMPLE_DIR="1_A01"
            MOVIE="m84094_260226_200655_s1"
            MUX_BC="bcM0001"
            ;;
        HC02_NK)
            SAMPLE_DIR="1_B01"
            MOVIE="m84094_260226_220949_s2"
            MUX_BC="bcM0002"
            ;;
        HC03_NK)
            SAMPLE_DIR="1_C01"
            MOVIE="m84094_260227_001254_s3"
            MUX_BC="bcM0003"
            ;;
    esac

    INPUT_BAM="${WORKING_DIR}/raw/${SAMPLE_DIR}/hifi_reads/${MOVIE}.hifi_reads.${MUX_BC}.bam"
    WORKSPACE_DIR="${WORKING_DIR}/workspace/${SAMPLE}"
    OUTS_DIR="${WORKING_DIR}/outs/${SAMPLE}"

    mkdir -p "${WORKSPACE_DIR}"
    mkdir -p "${OUTS_DIR}"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${SAMPLE}_isoseq
#SBATCH --output=${WORKING_DIR}/workspace/logs/${SAMPLE}_isoseq.log
#SBATCH --error=${WORKING_DIR}/workspace/logs/${SAMPLE}_isoseq.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=48:00:00

set -euo pipefail

source /data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/miniforge3/etc/profile.d/conda.sh
conda activate pacbio

cd ${WORKSPACE_DIR}

# segment reads
skera split -j \${SLURM_CPUS_PER_TASK} ${INPUT_BAM} ${REFERENCE_DIR}/../adapters/mas16_primers.fasta ${WORKSPACE_DIR}/segmented.bam

# primer removal
lima  --isoseq -j \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/segmented.bam ${REFERENCE_DIR}/../primers/10x_3kit_primers.fasta ${WORKSPACE_DIR}/fl.bam

# tag UMI and cell barcode
isoseq tag --design T-12U-16B -j \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/fl.5p--3p.bam ${WORKSPACE_DIR}/flt.bam

# trim polyA and remove concatemers
isoseq refine --require-polya -j \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/flt.bam ${REFERENCE_DIR}/../primers/10x_3kit_primers.fasta ${WORKSPACE_DIR}/fltnc.bam

# barcode correction
isoseq correct -j \${SLURM_CPUS_PER_TASK} --barcodes ${REFERENCE_DIR}/../barcodes/3M-3pgex-may-2023.REVCOMP.txt.gz ${WORKSPACE_DIR}/fltnc.bam ${WORKSPACE_DIR}/fltnc.corrected.bam

# sort by cell barcode
samtools sort -@ \${SLURM_CPUS_PER_TASK} -t CB ${WORKSPACE_DIR}/fltnc.corrected.bam -o ${WORKSPACE_DIR}/fltnc.corrected.sorted.bam

# deduplication
isoseq groupdedup -j \${SLURM_CPUS_PER_TASK} --keep-non-real-cells ${WORKSPACE_DIR}/fltnc.corrected.sorted.bam ${WORKSPACE_DIR}/dedup.bam

# Map to reference genome
pbmm2 align --preset ISOSEQ --sort -j \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/dedup.bam ${REFERENCE_DIR}/human_GRCh38_no_alt_analysis_set.fasta ${WORKSPACE_DIR}/mapped.bam

# Collapse into unique isoforms
isoseq collapse -j \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/mapped.bam ${WORKSPACE_DIR}/collapsed.gff

# Prepare collapsed GFF
pigeon prepare ${WORKSPACE_DIR}/collapsed.gff

# Classify isoforms
pigeon classify \
    --fl ${WORKSPACE_DIR}/collapsed.abundance.txt \
    --cage-peak ${REFERENCE_DIR}/refTSS_v3.3_human_coordinate.hg38.sorted.bed \
    --poly-a ${REFERENCE_DIR}/polyA.list.txt \
    -j \${SLURM_CPUS_PER_TASK} \
    ${WORKSPACE_DIR}/collapsed.sorted.gff \
    ${REFERENCE_DIR}/gencode.v39.annotation.sorted.gtf \
    ${REFERENCE_DIR}/human_GRCh38_no_alt_analysis_set.fasta \
    ${WORKSPACE_DIR}/collapsed

# Filter isoforms
pigeon filter -j \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/collapsed_classification.txt --isoforms ${WORKSPACE_DIR}/collapsed.sorted.gff

# Saturation report
pigeon report -j \${SLURM_CPUS_PER_TASK} ${WORKSPACE_DIR}/collapsed_classification.filtered_lite_classification.txt ${WORKSPACE_DIR}/saturation.txt

# Change to workspace directory and run pigeon make-seurat with relative paths
pigeon make-seurat \
    -j \${SLURM_CPUS_PER_TASK} \
    --dedup dedup.fasta \
    --group collapsed.group.txt \
    --keep-ribo-mito-genes \
    -d ${OUTS_DIR} \
    collapsed_classification.filtered_lite_classification.txt

# Remove intermediate files
rm -f ${WORKSPACE_DIR}/segmented.bam \
      ${WORKSPACE_DIR}/segmented.bam.pbi \
      ${WORKSPACE_DIR}/segmented.non_passing.bam \
      ${WORKSPACE_DIR}/segmented.non_passing.bam.pbi \
      ${WORKSPACE_DIR}/fl.5p--3p.bam \
      ${WORKSPACE_DIR}/fl.5p--3p.bam.pbi \
      ${WORKSPACE_DIR}/flt.bam \
      ${WORKSPACE_DIR}/flt.bam.pbi \
      ${WORKSPACE_DIR}/fltnc.bam \
      ${WORKSPACE_DIR}/fltnc.bam.pbi \
      ${WORKSPACE_DIR}/fltnc.corrected.bam \
      ${WORKSPACE_DIR}/fltnc.corrected.bam.pbi \
      ${WORKSPACE_DIR}/fltnc.corrected.sorted.bam \
      ${WORKSPACE_DIR}/fltnc.corrected_intermediate.bam.pbi
rm -f ${WORKSPACE_DIR}/*.consensusreadset.xml
rm -f ${WORKSPACE_DIR}/segmented.found_adapters.csv.gz \
      ${WORKSPACE_DIR}/segmented.read_lengths.csv \
      ${WORKSPACE_DIR}/segmented.ligations.csv
EOF

done