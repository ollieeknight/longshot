process ISOQUANT_DISCOVERY {
    tag "${experiment}"
    container "${params.container_isoquant}"
    publishDir { "${params.outdir}/${experiment}/joint/transcript_model" }, mode: 'copy'

    input:
    tuple val(experiment), val(metas), path(bams), path(bais)

    output:
    tuple val(experiment), path("isoquant_out/${experiment}/${experiment}.transcript_models.gtf"), emit: transcript_gtf

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
        --fl_data \\
        --check_canonical \\
        --count_exons \\
        --barcoded_bam \\
        --barcode_tag CB \\
        --umi_tag XM \\
        --read_group barcode \\
        --threads ${task.cpus} \\
        --prefix ${experiment} \\
        --output isoquant_out
    """
}


process SQANTI3_QC {
    tag "${experiment}"
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
    """
}


process SQANTI3_FILTER {
    tag "${experiment}"
    container "${params.container_sqanti3}"
    publishDir { "${params.outdir}/${experiment}/joint/sqanti3" }, mode: 'copy',
               saveAs: { fn -> fn.startsWith("sqanti_filter/") ? fn.replace("sqanti_filter/", "") : fn }

    input:
    tuple val(experiment), path(classification), path(corrected_fasta), path(corrected_gtf)

    output:
    tuple val(experiment), path("sqanti_filter/${experiment}_filtered_classification.txt"), emit: filtered_class
    tuple val(experiment), path("sqanti_filter/${experiment}_filtered.gtf"),               emit: filtered_gtf
    tuple val(experiment), path("sqanti_filter/${experiment}_filtered.fasta"),             emit: filtered_fasta

    script:
    """
    sqanti3_filter.py rules \\
        --sqanti_class ${classification} \\
        --filter_isoforms ${corrected_fasta} \\
        --filter_gtf ${corrected_gtf} \\
        -d sqanti_filter \\
        -o ${experiment}_filtered \\
        --skip_report
    """
}


process SQANTI3_RESCUE {
    tag "${experiment}"
    container "${params.container_sqanti3}"
    publishDir { "${params.outdir}/${experiment}/joint/sqanti3" }, mode: 'copy',
               saveAs: { fn -> fn.startsWith("sqanti_rescue/") ? fn.replace("sqanti_rescue/", "") : fn }

    input:
    tuple val(experiment), path(filter_classification), path(corrected_fasta), path(filtered_gtf)

    output:
    tuple val(experiment), path("sqanti_rescue/${experiment}_rescued.gtf"),   emit: rescued_gtf
    tuple val(experiment), path("sqanti_rescue/${experiment}_rescued.fasta"), emit: rescued_fasta
    path "sqanti_rescue/${experiment}_rescue_inclusion_list.tsv",             emit: inclusion_list

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
    """
}


// ── Plan 03: IsoQuant chromosome sharding ────────────────────────────────────

process ISOQUANT_DISCOVERY_SHARD {
    tag "${experiment} shard${shard.id} [${shard.chrs.join(',')}]"
    // ponytail: downgraded from process_ultra — each shard is ~1/17th the data
    container "${params.container_isoquant}"

    input:
    tuple val(experiment), val(shard), path(bams), path(bais)

    output:
    tuple val(experiment), val(shard),
          path("isoquant_out/${experiment}_shard${shard.id}/${experiment}_shard${shard.id}.transcript_models.gtf"),
          emit: shard_gtf

    script:
    def chr_list = shard.chrs.join(' ')
    def bam_args = bams instanceof List ? bams.join(' ') : bams.toString()
    def prefix   = "${experiment}_shard${shard.id}"
    """
    isoquant \\
        --reference ${params.ref_fasta} \\
        --genedb ${params.ref_gtf} \\
        --complete_genedb \\
        --process_only_chr ${chr_list} \\
        --bam ${bam_args} \\
        --data_type pacbio \\
        --fl_data \\
        --check_canonical \\
        --count_exons \\
        --barcoded_bam \\
        --barcode_tag CB \\
        --umi_tag XM \\
        --read_group barcode \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        --output isoquant_out
    """
}


process MERGE_SHARD_GTFS {
    tag "${experiment}"
    container "${params.container_multiqc}"
    publishDir { "${params.outdir}/${experiment}/joint/transcript_model" }, mode: 'copy'

    input:
    tuple val(experiment), val(shards), path(gtfs)

    output:
    tuple val(experiment), path("${experiment}.transcript_models.gtf"), emit: transcript_gtf

    script:
    """
    python3 -c "
import re, sys

gtf_files   = sorted('${gtfs}'.split())
out_file    = '${experiment}.transcript_models.gtf'
header_seen = set()
entry_seen  = set()

with open(out_file, 'w') as out:
    for gtf in gtf_files:
        with open(gtf) as f:
            for line in f:
                if line.startswith('#'):
                    if line not in header_seen:
                        header_seen.add(line)
                        out.write(line)
                    continue

                chrom = line.split('\\t')[0]

                # Prefix novel IDs with chromosome name to prevent collision across shards
                line = re.sub(r'(gene_id \\\"novel_gene_)',       r'\\g<1>' + chrom + r'_', line)
                line = re.sub(r'(transcript_id \\\"novel_transcript_)', r'\\g<1>' + chrom + r'_', line)
                line = re.sub(r'(gene_name \\\"novel_gene_)',      r'\\g<1>' + chrom + r'_', line)

                key = line.strip()
                if key not in entry_seen:
                    entry_seen.add(key)
                    out.write(line)
"
    """
}


// ── Plan 03 Part B: SQANTI3 chunking ─────────────────────────────────────────

process SQANTI3_SPLIT_GTF {
    tag "${experiment} (${params.sqanti_chunks} chunks)"
    container "${params.container_multiqc}"

    input:
    tuple val(experiment), path(gtf)

    output:
    tuple val(experiment), path("chunk_*.gtf"), emit: chunks

    script:
    """
    python3 -c "
import math

chunks = ${params.sqanti_chunks}
header_lines = []
transcript_blocks = []
current_block = []

with open('${gtf}') as f:
    for line in f:
        if line.startswith('#'):
            header_lines.append(line)
            continue
        if '\\t' not in line:
            continue
        fields = line.split('\\t')
        feature = fields[2]
        if feature == 'transcript':
            if current_block:
                transcript_blocks.append(current_block)
            current_block = [line]
        elif current_block:
            current_block.append(line)
if current_block:
    transcript_blocks.append(current_block)

per_chunk = math.ceil(len(transcript_blocks) / chunks)
for i in range(chunks):
    chunk_blocks = transcript_blocks[i*per_chunk:(i+1)*per_chunk]
    if not chunk_blocks:
        continue
    with open(f'chunk_{i+1:02d}.gtf', 'w') as out:
        out.writelines(header_lines)
        for block in chunk_blocks:
            out.writelines(block)
"
    """
}


process SQANTI3_QC_CHUNK {
    tag "${experiment} chunk${chunk_id}"
    // ponytail: downgraded from process_high — 1/8th the transcripts per chunk
    container "${params.container_sqanti3}"

    input:
    tuple val(experiment), val(chunk_id), path(chunk_gtf)

    output:
    tuple val(experiment), val(chunk_id),
          path("sqanti_qc/${experiment}_chunk${chunk_id}_classification.txt"),
          path("sqanti_qc/${experiment}_chunk${chunk_id}_corrected.fasta"),
          path("sqanti_qc/${experiment}_chunk${chunk_id}_corrected.gtf"), emit: chunk_results

    script:
    def intropolis_arg = params.intropolis ? "--coverage ${params.intropolis}" : ""
    def polya_peak_arg = params.polya_peak ? "--polyA_peak ${params.polya_peak}" : ""
    """
    sqanti3_qc.py \\
        --isoforms ${chunk_gtf} \\
        --refGTF ${params.ref_gtf} \\
        --refFasta ${params.ref_fasta} \\
        --CAGE_peak ${params.cage_peaks} \\
        --polyA_motif_list ${params.polya_list} \\
        ${polya_peak_arg} \\
        ${intropolis_arg} \\
        -t ${task.cpus} \\
        -d sqanti_qc \\
        -o ${experiment}_chunk${chunk_id} \\
        --report skip
    """
}


process SQANTI3_MERGE_CHUNKS {
    tag "${experiment}"
    container "${params.container_multiqc}"
    publishDir { "${params.outdir}/${experiment}/joint/sqanti3" }, mode: 'copy',
               saveAs: { fn -> fn.startsWith("sqanti_qc/") ? fn.replace("sqanti_qc/", "") : fn }

    input:
    tuple val(experiment), val(chunk_ids), path(classifications), path(fastas), path(gtfs)

    output:
    tuple val(experiment), path("sqanti_qc/${experiment}_classification.txt"),
                           path("sqanti_qc/${experiment}_corrected.fasta"),
                           path("sqanti_qc/${experiment}_corrected.gtf"), emit: sqanti_results

    script:
    """
    mkdir -p sqanti_qc

    # Merge classification TSVs: header from first chunk, data rows from all
    cls_arr=( ${classifications} )
    head -1 "\${cls_arr[0]}" > sqanti_qc/${experiment}_classification.txt
    for f in "\${cls_arr[@]}"; do
        tail -n+2 "\$f" >> sqanti_qc/${experiment}_classification.txt
    done

    # Concatenate corrected FASTAs
    cat ${fastas} > sqanti_qc/${experiment}_corrected.fasta

    # Concatenate corrected GTFs
    cat ${gtfs} > sqanti_qc/${experiment}_corrected.gtf
    """
}


process ISOQUANT_QUANTIFY {
    tag "${meta.sample_id}"
    container "${params.container_isoquant}"
    publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/counts" }, mode: 'copy',
               saveAs: { fn -> fn.startsWith("isoquant_out/") ? fn.replaceFirst("isoquant_out/[^/]+/", "") : fn }

    input:
    tuple val(meta), path(bam), path(bai), path(filtered_gtf)

    output:
    tuple val(meta), path("isoquant_out/${meta.library_id}/"), emit: counts_dir

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
    """
}
