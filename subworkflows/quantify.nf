include { ISOQUANT_QUANTIFY  } from '../modules/quantify'

workflow QUANTIFY {
    take:
    ch_aligned_bam  // [meta, bam, bai] — per library
    ch_filtered_gtf // [experiment, gtf]

    main:
    // ── Step 1: Per-library quantification using filtered GTF ─────────────────
    // Combine each library's original aligned BAM with its experiment's filtered GTF
    ch_aligned_bam
        .map { meta, bam, bai -> [ meta.experiment, meta, bam, bai ] }
        .combine(ch_filtered_gtf, by: 0)
        .map { exp, meta, bam, bai, gtf -> [ meta, bam, bai, gtf ] }
        | ISOQUANT_QUANTIFY

    emit:
    counts_dir = ISOQUANT_QUANTIFY.out.counts_dir
    versions   = ISOQUANT_QUANTIFY.out.versions
}
