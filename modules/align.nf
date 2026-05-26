process SAMTOOLS_MERGE_FLTNC {
    tag "${meta.sample_id}"
    label 'process_medium'
    container "${params.container_samtools}"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.sample_id}_fltnc_merged.bam"), emit: merged_bam
    path "versions.yml",                                         emit: versions

    script:
    def bam_list = bams instanceof List ? bams : [bams]
    def merge_cmd = bam_list.size() == 1
        ? "cp ${bam_list[0]} ${meta.sample_id}_fltnc_merged.bam"
        : "samtools merge -f -@ ${task.cpus} ${meta.sample_id}_fltnc_merged.bam ${bam_list.join(' ')}"
    """
    ${merge_cmd}
    samtools index -@ ${task.cpus} ${meta.sample_id}_fltnc_merged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}


process PREPARE_WHITELIST {
    tag "${meta.sample_id}"
    label 'process_low'
    container "${params.container_multiqc}"

    input:
    tuple val(meta), path(sr_barcodes)

    output:
    tuple val(meta), path("${meta.sample_id}_whitelist.txt"), emit: whitelist
    path "versions.yml",                                       emit: versions

    script:
    """
    python3 -c "
import gzip
import lzma
import sys

def rev_comp(seq):
    tb = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(tb)[::-1]

input_file = '${sr_barcodes}'
output_file = '${meta.sample_id}_whitelist.txt'

open_fun = lzma.open if input_file.endswith('.xz') else (gzip.open if input_file.endswith('.gz') else open)
with open_fun(input_file, 'rt') as f_in:
    barcodes = [line.strip().split('-')[0] for line in f_in if line.strip()]

with open(output_file, 'w') as f_out:
    for bc in barcodes:
        f_out.write(rev_comp(bc) + '\\n')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | grep -oP '(?<=Python )\\S+')
    END_VERSIONS
    """
}


process ISOSEQ_CORRECT {
    tag "${meta.sample_id}"
    label 'process_high'
    container "${params.container_isoseq}"

    input:
    tuple val(meta), path(fltnc_bam), path(whitelist)

    output:
    tuple val(meta), path("${meta.sample_id}_corrected.bam"), emit: corrected_bam
    path "versions.yml",                                      emit: versions

    script:
    def wl = whitelist.name.endsWith('.xz') ? 'whitelist_isoseq.txt.gz' : whitelist
    def decomp = whitelist.name.endsWith('.xz') ? "xz -dkc ${whitelist} | gzip -c > whitelist_isoseq.txt.gz" : ''
    """
    ${decomp}
    isoseq correct \\
        -j ${task.cpus} \\
        --barcodes ${wl} \\
        ${fltnc_bam} \\
        ${meta.sample_id}_corrected.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoseq: \$(isoseq --version 2>&1 | grep -oP '(?<=isoseq )\\S+')
    END_VERSIONS
    """
}


process SAMTOOLS_SORT_CB {
    tag "${meta.sample_id}"
    label 'process_medium'
    container "${params.container_samtools}"

    input:
    tuple val(meta), path(corrected_bam)

    output:
    tuple val(meta), path("${meta.sample_id}_sorted_cb.bam"), emit: sorted_bam
    path "versions.yml",                                      emit: versions

    script:
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -t CB \\
        ${corrected_bam} \\
        -o ${meta.sample_id}_sorted_cb.bam

    samtools index -@ ${task.cpus} ${meta.sample_id}_sorted_cb.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}


process ISOSEQ_GROUPDEDUP {
    tag "${meta.sample_id}"
    label 'process_high'
    container "${params.container_isoseq}"

    input:
    tuple val(meta), path(sorted_bam)

    output:
    tuple val(meta), path("${meta.sample_id}_dedup.bam"), emit: dedup_bam
    path "versions.yml",                                  emit: versions

    script:
    """
    isoseq groupdedup \\
        -j ${task.cpus} \\
        --keep-non-real-cells \\
        ${sorted_bam} \\
        ${meta.sample_id}_dedup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoseq: \$(isoseq --version 2>&1 | grep -oP '(?<=isoseq )\\S+')
    END_VERSIONS
    """
}


process PBMM2_ALIGN {
    tag "${meta.sample_id}"
    label 'process_high'
    container "${params.container_pbmm2}"

    input:
    tuple val(meta), path(dedup_bam)

    output:
    tuple val(meta), path("${meta.sample_id}_aligned.bam"), path("${meta.sample_id}_aligned.bam.bai"), emit: aligned_bam
    path "versions.yml",                                                                                emit: versions

    script:
    """
    pbmm2 align \\
        --preset ISOSEQ \\
        --sort \\
        -j ${task.cpus} \\
        ${dedup_bam} \\
        ${params.ref_fasta} \\
        ${meta.sample_id}_aligned.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmm2: \$(pbmm2 --version 2>&1 | grep -oP '(?<=pbmm2 )\\S+')
    END_VERSIONS
    """
}


process CB_SUFFIX_INJECT {
    tag "${meta.sample_id}"
    label 'process_medium'
    container "${params.container_samtools}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample_id}_suffixed.bam"), path("${meta.sample_id}_suffixed.bam.bai"), emit: bam
    path "versions.yml",                                                                                  emit: versions

    script:
    """
    samtools view -h ${bam} \\
    | awk -v OFS="\\t" -v suf="${meta.suffix}" '{
        if (\$0 ~ /^@/) { print }
        else {
            for (i = 12; i <= NF; i++) {
                if (\$i ~ /^CB:Z:/) { \$i = \$i suf }
            }
            print \$0
        }
    }' \\
    | samtools view -@ ${task.cpus} -b -o ${meta.sample_id}_suffixed.bam

    samtools index -@ ${task.cpus} ${meta.sample_id}_suffixed.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}


process GENERATE_CRAM {
    tag "${meta.sample_id}"
    label 'process_low'
    container "${params.container_samtools}"
    publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/alignment" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample_id}.cram"), path("${meta.sample_id}.cram.crai"), emit: cram
    path "versions.yml",                                                                    emit: versions

    script:
    """
    samtools view \\
        -@ ${task.cpus} \\
        -C \\
        -T ${params.ref_fasta} \\
        -o ${meta.sample_id}.cram \\
        ${bam}
    samtools index ${meta.sample_id}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
