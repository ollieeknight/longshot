process ISOQUANT_DISCOVERY {
    tag "${experiment}"
    label 'process_ultra'
    container "${params.container_isoquant}"
    publishDir { "${params.outdir}/${experiment}/joint/transcript_model" }, mode: 'copy'

    input:
    tuple val(experiment), val(metas), path(bams), path(bais)

    output:
    tuple val(experiment), path("isoquant_out/${experiment}/${experiment}.transcript_models.gtf"), emit: transcript_gtf
    path "versions.yml",                                                                           emit: versions

    script:
    def chrs = ((1..22).collect { "chr${it}" } + ['chrX', 'chrY', 'chrM']).join(' ')
    def bam_args = bams instanceof List ? bams.join(' ') : bams.toString()
    """
    isoquant \\
        --reference ${params.ref_fasta} \\
        --genedb ${params.ref_gtf} \\
        --complete_genedb \\
        --process_only_chr ${chrs} \\
        --bam ${bam_args} \\
        --data_type pacbio \\
        --barcoded_bam \\
        --barcode_tag CB \\
        --umi_tag XM \\
        --read_group barcode \\
        --threads ${task.cpus} \\
        --prefix ${experiment} \\
        --output isoquant_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(isoquant --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1)
    END_VERSIONS
    """
}


process SQANTI3_QC {
    tag "${experiment}"
    label 'process_high'
    container "${params.container_sqanti3}"
    publishDir { "${params.outdir}/${experiment}/joint/sqanti3" }, mode: 'copy',
               saveAs: { fn -> fn.startsWith("sqanti_qc/") ? fn.replace("sqanti_qc/", "") : fn }

    input:
    tuple val(experiment), path(isoforms_gtf)

    output:
    tuple val(experiment), path("sqanti_qc/${experiment}_classification.txt"),
                           path("sqanti_qc/${experiment}_corrected.fasta"),
                           path("sqanti_qc/${experiment}_corrected.gtf"),   emit: sqanti_results
    path "sqanti_qc/${experiment}_junctions.txt",                           emit: junctions
    path "versions.yml",                                                    emit: versions

    script:
    def intropolis_arg = params.intropolis ? "--coverage ${params.intropolis}" : ""
    def polya_peak_arg = params.polya_peak ? "--polyA_peak ${params.polya_peak}" : ""
    """
    sqanti3_qc.py \\
        --isoforms ${isoforms_gtf} \\
        --refGTF ${params.ref_gtf} \\
        --refFasta ${params.ref_fasta} \\
        --CAGE_peak ${params.cage_peaks} \\
        --polyA_motif_list ${params.polya_list} \\
        ${polya_peak_arg} \\
        ${intropolis_arg} \\
        -t ${task.cpus} \\
        -d sqanti_qc \\
        -o ${experiment} \\
        --report pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_qc.py --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1)
    END_VERSIONS
    """
}


process SQANTI3_FILTER {
    tag "${experiment}"
    label 'process_medium'
    container "${params.container_sqanti3}"
    publishDir { "${params.outdir}/${experiment}/joint/sqanti3" }, mode: 'copy',
               saveAs: { fn -> fn.startsWith("sqanti_filter/") ? fn.replace("sqanti_filter/", "") : fn }

    input:
    tuple val(experiment), path(classification), path(corrected_fasta), path(corrected_gtf)

    output:
    tuple val(experiment), path("sqanti_filter/${experiment}_filtered_classification.txt"), emit: filtered_class
    tuple val(experiment), path("sqanti_filter/${experiment}_filtered.gtf"),               emit: filtered_gtf
    tuple val(experiment), path("sqanti_filter/${experiment}_filtered.fasta"),             emit: filtered_fasta
    path "versions.yml",                                                                   emit: versions

    script:
    """
    sqanti3_filter.py rules \\
        --sqanti_class ${classification} \\
        --filter_isoforms ${corrected_fasta} \\
        --filter_gtf ${corrected_gtf} \\
        -d sqanti_filter \\
        -o ${experiment}_filtered \\
        --skip_report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_filter.py --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1)
    END_VERSIONS
    """
}


process SQANTI3_RESCUE {
    tag "${experiment}"
    label 'process_medium'
    container "${params.container_sqanti3}"
    publishDir { "${params.outdir}/${experiment}/joint/sqanti3" }, mode: 'copy',
               saveAs: { fn -> fn.startsWith("sqanti_rescue/") ? fn.replace("sqanti_rescue/", "") : fn }

    input:
    tuple val(experiment), path(filter_classification), path(corrected_fasta), path(filtered_gtf)

    output:
    tuple val(experiment), path("sqanti_rescue/${experiment}_rescued.gtf"),   emit: rescued_gtf
    tuple val(experiment), path("sqanti_rescue/${experiment}_rescued.fasta"), emit: rescued_fasta
    path "sqanti_rescue/${experiment}_rescue_inclusion_list.tsv",             emit: inclusion_list
    path "versions.yml",                                                      emit: versions

    script:
    """
    sqanti3_rescue.py rules \\
        --filter_class ${filter_classification} \\
        --corrected_isoforms_fasta ${corrected_fasta} \\
        --filtered_isoforms_gtf ${filtered_gtf} \\
        -rg ${params.ref_gtf} \\
        -rf ${params.ref_fasta} \\
        --mode automatic \\
        -d sqanti_rescue \\
        -o ${experiment}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_rescue.py --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1)
    END_VERSIONS
    """
}


process ISOQUANT_QUANTIFY {
    tag "${meta.sample_id}"
    label 'process_high'
    container "${params.container_isoquant}"
    publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/counts" }, mode: 'copy',
               saveAs: { fn -> fn.startsWith("isoquant_out/") ? fn.replaceFirst("isoquant_out/[^/]+/", "") : fn }

    input:
    tuple val(meta), path(bam), path(bai), path(filtered_gtf)

    output:
    tuple val(meta), path("isoquant_out/${meta.library_id}/"), emit: counts_dir
    path "versions.yml",                                       emit: versions

    script:
    def chrs = ((1..22).collect { "chr${it}" } + ['chrX', 'chrY', 'chrM']).join(' ')
    """
    isoquant \\
        --reference ${params.ref_fasta} \\
        --genedb ${filtered_gtf} \\
        --complete_genedb \\
        --process_only_chr ${chrs} \\
        --bam ${bam} \\
        --data_type pacbio \\
        --barcoded_bam \\
        --barcode_tag CB \\
        --umi_tag XM \\
        --read_group barcode \\
        --threads ${task.cpus} \\
        --prefix ${meta.library_id} \\
        --output isoquant_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(isoquant --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1)
    END_VERSIONS
    """
}
