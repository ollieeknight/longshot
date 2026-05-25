#!/usr/bin/bash

WORK_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
OUT_DIR="${WORK_DIR}/downstream/07_population_variants"
LOG_DIR="${OUT_DIR}/logs"
DONORS=(HC01 HC02 HC03)

mkdir -p "$LOG_DIR"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=joint_call
#SBATCH --output=${LOG_DIR}/%x_%A.log
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail
shopt -s nullglob

WORK_DIR="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
OUT_DIR="\${WORK_DIR}/downstream/07_population_variants"
JOINT_DIR="\${OUT_DIR}/joint"

command -v bcftools >/dev/null 2>&1 || { echo "Missing bcftools"; exit 1; }
mkdir -p "\$JOINT_DIR"

echo "Creating per-donor master VCFs..."

# Find all donors dynamically
donors=(\$(ls -d "\${OUT_DIR}"/*/ | xargs -n1 basename))

for donor in "\${donors[@]}"; do
    echo "Processing \${donor}..."
    donor_dir="\${JOINT_DIR}/\${donor}"
    mkdir -p "\$donor_dir"
    
    # Collect all population VCFs for this donor
    all_vcfs=\$(mktemp)
    for pop_dir in "\${OUT_DIR}/\${donor}/clair3"/*/; do
        pop=\$(basename "\$pop_dir")
        vcf="\${OUT_DIR}/\${donor}/clair3/\${pop}/pileup.vcf.gz"
        [[ -f "\$vcf" ]] && echo "\$vcf" >> "\$all_vcfs"
    done
    
    [[ -s "\$all_vcfs" ]] || { echo "  No VCFs found for \${donor}"; rm "\$all_vcfs"; continue; }
    n_pops=\$(wc -l < "\$all_vcfs")
    
    echo "  Merging \${n_pops} populations into master VCF..."
    bcftools merge -l "\$all_vcfs" --threads 4 -O z -o "\${donor_dir}/\${donor}_master.vcf.gz"
    bcftools index -t "\${donor_dir}/\${donor}_master.vcf.gz"
    
    rm "\$all_vcfs"
done

echo "Done. Master VCFs per donor: \${JOINT_DIR}/*/*.vcf.gz"
EOF

