process COLLECT_INSTRUMENT_STATS {
    tag "${meta.id}"
    label 'process_low'
    container "${params.container_samtools}"

    input:
    tuple val(meta), path(stats_dir)

    output:
    path "instrument_stats/*", optional: true, emit: stats
    path "versions.yml",                        emit: versions

    script:
    """
    mkdir -p instrument_stats
    if [ -d "${stats_dir}" ]; then
        for f in "${stats_dir}"/*.hifi_reads.lima_summary.txt \
                 "${stats_dir}"/*.hifi_reads.lima_counts.txt \
                 "${stats_dir}"/*.ccs_report.json; do
            [ -f "\$f" ] && cp "\$f" instrument_stats/ || true
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}


process SAMTOOLS_FLAGSTAT {
    tag "${meta.sample_id} (${stage})"
    label 'process_low'
    container "${params.container_samtools}"

    input:
    tuple val(meta), val(stage), path(bam)

    output:
    path "${meta.sample_id}_${stage}.flagstat", emit: flagstat
    path "versions.yml",                         emit: versions

    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${meta.sample_id}_${stage}.flagstat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}


process MOSDEPTH {
    tag "${meta.sample_id}"
    label 'process_low'
    container "${params.container_mosdepth}"
    publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/qc/mosdepth" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample_id}.*"),    emit: results
    path "${meta.sample_id}.mosdepth.summary.txt",   emit: summary
    path "versions.yml",                              emit: versions

    script:
    """
    mosdepth \\
        --threads ${task.cpus} \\
        --fast-mode \\
        --no-abbrev \\
        ${meta.sample_id} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | grep -oP '(?<=mosdepth )\\S+')
    END_VERSIONS
    """
}


process NANOSTAT {
    tag "${meta.sample_id}"
    label 'process_low'
    container "${params.container_nanostat}"
    publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/qc/nanostat" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "${meta.sample_id}_NanoStats.txt", emit: stats
    path "versions.yml",                    emit: versions

    script:
    """
    NanoStat \\
        --bam ${bam} \\
        --name ${meta.sample_id} \\
        --outdir . \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanostat: \$(NanoStat --version 2>&1 | sed 's/NanoStat //')
    END_VERSIONS
    """
}


process MULTIQC {
    label 'process_low'
    container "${params.container_multiqc}"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path reports, stageAs: 'reports/?/*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data

    script:
    def config_arg = params.multiqc_config ? "--config ${params.multiqc_config}" : ""
    """
    multiqc \\
        ${config_arg} \\
        --force \\
        reports/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | grep -oP '(?<=version )\\S+')
    END_VERSIONS
    """
}
