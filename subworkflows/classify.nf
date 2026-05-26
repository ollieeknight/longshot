include { CB_SUFFIX_INJECT   } from '../modules/align'
include { ISOQUANT_DISCOVERY } from '../modules/quantify'
include { SQANTI3_QC         } from '../modules/quantify'
include { SQANTI3_FILTER     } from '../modules/quantify'

workflow CLASSIFY {
    take:
    ch_aligned_bam // [meta, bam, bai] — per library

    main:
    // ── Step 1: Assign CB suffixes deterministically (sorted by library_id) ────
    ch_aligned_bam
        .map { meta, bam, bai -> [ meta.experiment, meta, bam, bai ] }
        .groupTuple(by: 0)
        .flatMap { experiment, metas, bams, bais ->
            def items = [metas, bams, bais].transpose()
            def sorted = items.sort { a, b -> a[0].library_id <=> b[0].library_id }
            sorted.withIndex().collect { item, idx ->
                def (meta, bam, bai) = item
                [ meta + [suffix: String.format("_%02d", idx + 1)], bam, bai ]
            }
        }
        .set { ch_with_suffix }

    // ── Step 2: Inject CB suffix into BAM tags ────────────────────────────────
    CB_SUFFIX_INJECT(ch_with_suffix)

    // ── Step 3: Group suffixed BAMs by experiment for joint discovery ─────────
    CB_SUFFIX_INJECT.out.bam
        .map { meta, bam, bai -> [ meta.experiment, meta, bam, bai ] }
        .groupTuple(by: 0)
        .set { ch_experiment_bams }

    // ── Step 4: Joint transcript model discovery ──────────────────────────────
    ISOQUANT_DISCOVERY(ch_experiment_bams)

    // ── Step 5: SQANTI3 QC + Filter ───────────────────────────────────────────
    SQANTI3_QC(ISOQUANT_DISCOVERY.out.transcript_gtf)
    SQANTI3_FILTER(SQANTI3_QC.out.sqanti_results)

    emit:
    filtered_gtf    = SQANTI3_FILTER.out.filtered_gtf
    filtered_fasta  = SQANTI3_FILTER.out.filtered_fasta
    filtered_class  = SQANTI3_FILTER.out.filtered_class
    sqanti_class    = SQANTI3_QC.out.sqanti_results.map { exp, cls, fa, gtf -> cls }
    versions        = Channel.empty()
                      .mix(CB_SUFFIX_INJECT.out.versions)
                      .mix(ISOQUANT_DISCOVERY.out.versions)
                      .mix(SQANTI3_QC.out.versions)
                      .mix(SQANTI3_FILTER.out.versions)
}
