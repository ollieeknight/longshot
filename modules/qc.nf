process COLLECT_INSTRUMENT_STATS {
    tag "${meta.id}"
    container "${params.container_samtools}"
    publishDir { "${params.outdir}/${meta.id}/qc/instrument" }, mode: 'copy'

    input:
    tuple val(meta), path(stats_dir)

    output:
    path "instrument_stats/*", optional: true, emit: stats

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
    """
}


process SAMTOOLS_FLAGSTAT {
    tag "${meta.sample_id} ${stage}"
    container "${params.container_samtools}"
    publishDir { "${params.outdir}/${meta.experiment ?: meta.id}/qc/flagstat" }, mode: 'copy'

    input:
    tuple val(meta), val(stage), path(bam)

    output:
    path "${meta.sample_id}_${stage}.flagstat", emit: flagstat

    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${meta.sample_id}_${stage}.flagstat
    """
}


process CRAMINO {
    tag "${meta.sample_id ?: meta.id}"
    container "${params.container_cramino}"
    publishDir { "${params.outdir}/${meta.experiment ?: meta.id}/qc/cramino" }, mode: 'copy'

    input:
    tuple val(meta), val(stage), val(extra_flags), path(bam)

    output:
    path "${meta.sample_id ?: "${meta.id}_${meta.run_id}"}_${stage}.cramino.txt", emit: stats

    script:
    // sample_id is set post-merge (unique across runs); smrt_meta (pre-demux) has no
    // sample_id, so fall back to id+run_id to keep filenames unique per SMRT cell.
    def out_id = meta.sample_id ?: "${meta.id}_${meta.run_id}"
    """
    cramino \\
        --threads ${task.cpus} \\
        ${extra_flags} \\
        ${bam} \\
        > ${out_id}_${stage}.cramino.txt
    """
}


process MULTIQC {
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
    """
}
