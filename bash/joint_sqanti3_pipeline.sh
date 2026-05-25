#!/bin/bash
#SBATCH --job-name=joint_nk_sqanti3
#SBATCH --output=/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk/workspace/logs/joint_sqanti3.log
#SBATCH --error=/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk/workspace/logs/joint_sqanti3.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=600G
#SBATCH --time=48:00:00

set -euo pipefail

# --- Setup Paths ---
BASE_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
WORKSPACE_DIR="${BASE_DIR}/workspace"
JOINT_DIR="${WORKSPACE_DIR}/joint_analysis"
REFERENCE_DIR="/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/pacbio/reference"
SQANTI_SIF="/data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/sqanti3_latest.sif"

mkdir -p "${JOINT_DIR}"
cd "${JOINT_DIR}"

source /data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/miniforge3/etc/profile.d/conda.sh
conda activate pacbio

if false; then # SKIPPING STEPS 1 AND 2

echo "======================================================"
echo " Step 1: Append Suffix to Cell Barcodes & Merge BAMs"
echo "======================================================"
MAPPED_BAM_LIST=""
DEDUP_BAM_LIST=""

for SAMPLE in HC01_NK HC02_NK HC03_NK; do
    echo "Processing $SAMPLE..."
    MAPPED_BAM="${WORKSPACE_DIR}/${SAMPLE}/mapped.bam"
    DEDUP_BAM="${WORKSPACE_DIR}/${SAMPLE}/dedup.bam"
    
    # We will suffix the CB column with _HC01, etc.
    SUFFIX="_${SAMPLE}"
    
    MAPPED_MOD="${JOINT_DIR}/${SAMPLE}_mapped_mod.bam"
    DEDUP_MOD="${JOINT_DIR}/${SAMPLE}_dedup_mod.bam"
    DEDUP_FASTA="${JOINT_DIR}/${SAMPLE}_dedup_mod.fasta"

    # Use awk to inject the sample suffix into the CB (Cell Barcode) tags for mapped BAM
    if [ ! -f "${MAPPED_MOD}" ]; then
        samtools view -h "${MAPPED_BAM}" | awk -v OFS="\t" -v suf="${SUFFIX}" '{
            if ($0 ~ /^@/) { print }
            else {
                for(i=12; i<=NF; i++) {
                    if ($i ~ /^CB:Z:/) { $i = $i suf }
                }
                print $0
            }
        }' | samtools view -@ ${SLURM_CPUS_PER_TASK} -b - > "${MAPPED_MOD}"
    fi

    # Do the exact same string manipulation for the dedup BAM (needed by pigeon later for counting UMIs)
    if [ ! -f "${DEDUP_MOD}" ]; then
        samtools view -h "${DEDUP_BAM}" | awk -v OFS="\t" -v suf="${SUFFIX}" '{
            if ($0 ~ /^@/) { print }
            else {
                for(i=12; i<=NF; i++) {
                    if ($i ~ /^CB:Z:/) { $i = $i suf }
                }
                print $0
            }
        }' | samtools view -@ ${SLURM_CPUS_PER_TASK} -b - > "${DEDUP_MOD}"
        
        # We also need the dedup fasta for Pigeon
        samtools fasta -@ ${SLURM_CPUS_PER_TASK} "${DEDUP_MOD}" > "${DEDUP_FASTA}"
    fi

    MAPPED_BAM_LIST="${MAPPED_BAM_LIST} ${MAPPED_MOD}"
    DEDUP_BAM_LIST="${DEDUP_BAM_LIST} ${DEDUP_FASTA}"
done

echo "Merging mapped BAM files..."
samtools merge -@ ${SLURM_CPUS_PER_TASK} -f "${JOINT_DIR}/merged_mapped.bam" ${MAPPED_BAM_LIST}

echo "Merging dedup FASTA files..."
cat ${DEDUP_BAM_LIST} > "${JOINT_DIR}/merged_dedup.fasta"

echo "======================================================"
echo " Step 2: Joint IsoSeq Collapse"
echo "======================================================"
isoseq collapse -j ${SLURM_CPUS_PER_TASK} \
    "${JOINT_DIR}/merged_mapped.bam" \
    "${JOINT_DIR}/joint_collapsed.gff"

echo "======================================================"
echo " Step 3: Comprehensive QC with SQANTI3"
echo "======================================================"
# SQANTI3's built-in CSV parser crashes on the massive "read names" column generated 
# by the single-cell isoseq collapse. We extract just the Isoform ID and Count columns 
# to bypass the 131KB Python field size limit.
cut -f1,2 "${JOINT_DIR}/joint_collapsed.abundance.txt" > "${JOINT_DIR}/clean_abundance.txt"

# Note: we use 'apptainer run' instead of 'exec' so the container's entrypoint
# (which activates conda) is used. If it requires an env name, we pass `-n sqanti3`. 
apptainer run -B /data "${SQANTI_SIF}" sqanti3_qc.py \
    --isoforms "${JOINT_DIR}/joint_collapsed.gff" \
    --refGTF "${REFERENCE_DIR}/gencode.v39.annotation.sorted.gtf" \
    --refFasta "${REFERENCE_DIR}/human_GRCh38_no_alt_analysis_set.fasta" \
    --CAGE_peak "${REFERENCE_DIR}/refTSS_v3.3_human_coordinate.hg38.sorted.bed" \
    --polyA_motif_list "${REFERENCE_DIR}/polyA.list.txt" \
    --coverage "${REFERENCE_DIR}/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv" \
    -t ${SLURM_CPUS_PER_TASK} \
    -d "${JOINT_DIR}/sqanti_qc" \
    -o joint_nk \
    --report pdf \
    --fl_count "${JOINT_DIR}/clean_abundance.txt"

echo "======================================================"
echo " Step 4: Strict Filtering (SQANTI3 Rules Filter)"
echo "======================================================"
apptainer run -B /data "${SQANTI_SIF}" sqanti3_filter.py rules \
    --sqanti_class "${JOINT_DIR}/sqanti_qc/joint_nk_classification.txt" \
    --filter_isoforms "${JOINT_DIR}/sqanti_qc/joint_nk_corrected.fasta" \
    --filter_gtf "${JOINT_DIR}/sqanti_qc/joint_nk_corrected.gtf" \
    -d "${JOINT_DIR}/sqanti_filter" \
    -o joint_nk_filtered \
    --skip_report

fi # END SKIP

echo "======================================================"
echo " Step 5: Generate Single-Cell Count Matrices (Pigeon)"
echo "======================================================"
# Pigeon requires 'filtered_collapsed.group.txt' and a classification format input.
# We create a mapping file that Pigeon natively understands based on SQANTI3's passing IDs.

cat << EOF > "${JOINT_DIR}/make_pigeon_compat.py"
import pandas as pd
import sys

# SQANTI3 produces an inclusion-list.txt containing exactly the IDs that passed the filter.
filter_res = "${JOINT_DIR}/sqanti_filter/joint_nk_filtered_inclusion-list.txt"
group_in   = "${JOINT_DIR}/joint_collapsed.group.txt"
sq_class   = "${JOINT_DIR}/sqanti_qc/joint_nk_classification.txt"

print(f"Reading filtered classify file: {filter_res}")
# 1. Get IDs that passed the filter
try:
    with open(filter_res) as f:
        pass_ids = set([line.strip() for line in f if line.strip() != "isoform" and line.strip() != ""])
    print(f"Successfully read inclusion list. Found {len(pass_ids)} passed isoforms.")
except Exception as e:
    print("Warning: Could not read inclusion list:", e)
    pass_ids = set()

print("First 5 pass_ids:", list(pass_ids)[:5])

print(f"Reading group file: {group_in}")
# 2. Filter collapsed group file
# The grouping file contains the mapping: PB.XX.Y \t cell_CB1,cell_CB2...
df_group = pd.read_csv(group_in, sep='\t', header=None, names=["isoform", "reads"])
print(f"Found {len(df_group)} rows in group file.")
print("First 5 isoforms in group file:", df_group['isoform'].head().tolist())

df_group_filtered = df_group[df_group['isoform'].isin(pass_ids)]
print(f"After filtering, {len(df_group_filtered)} rows remain.")

if len(df_group_filtered) == 0:
    print("ERROR: No rows passed filtering. Check isoform ID formatting between files!")

df_group_filtered.to_csv("${JOINT_DIR}/joint_filtered_collapsed.group.txt", sep='\t', index=False, header=False)

# 3. Create a Pigeon-compatible minimal classification file
print("Creating pigeon compatible class file...")

# Pigeon's 'make-seurat' completely crashes due to strings like 'NA' in specific columns that it expects to be purely integer.
# SQANTI3 introduces 'NA' strings, so we must clean these up.
df_cls = pd.read_csv(sq_class, sep='\t', low_memory=False)
df_cls = df_cls[df_cls['isoform'].isin(pass_ids)]

# Pigeon expects structural categories to either be integer columns, or completely NA free.
# Let's replace 'NA' strings with empty strings or 0 for columns that are traditionally numeric.
for col in df_cls.columns:
    if df_cls[col].dtype == object:
        # Check if the column consists mostly of digits mixed with 'NA's
        non_na = df_cls[col][df_cls[col] != 'NA'].dropna()
        if len(non_na) > 0 and non_na.astype(str).str.isdigit().all():
            df_cls[col] = df_cls[col].replace('NA', '0').fillna('0').astype(int)

df_cls.to_csv("${JOINT_DIR}/pigeon_compatible_class.txt", sep='\t', index=False)
print("make_pigeon_compat.py finished successfully.")
EOF

python "${JOINT_DIR}/make_pigeon_compat.py"

mkdir -p "${JOINT_DIR}/seurat_out"
pigeon make-seurat -j ${SLURM_CPUS_PER_TASK} \
    --group "${JOINT_DIR}/joint_filtered_collapsed.group.txt" \
    --dedup "${JOINT_DIR}/merged_dedup.fasta" \
    --keep-ribo-mito-genes \
    -d "${JOINT_DIR}/seurat_out" \
    "${JOINT_DIR}/pigeon_compatible_class.txt"

echo "Pipeline successfully completed. Filtered Seurat objects are in ${JOINT_DIR}/seurat_out"
