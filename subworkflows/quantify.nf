include { CB_SUFFIX_INJECT   } from '../modules/align'
include { ISOQUANT_DISCOVERY } from '../modules/quantify'
include { SQANTI3_QC         } from '../modules/quantify'
include { SQANTI3_FILTER     } from '../modules/quantify'
include { ISOQUANT_QUANTIFY  } from '../modules/quantify'
include { EXPORT_LIBRARY_MTX; GENERATE_SHARED_CATALOG; CALCULATE_SATURATION } from '../modules/exporter'


workflow QUANTIFY {
    take:
    ch_aligned_bam  // [meta, bam, bai] — per library; meta has experiment + library_id

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

    // ── Step 6: Per-library quantification using filtered GTF ─────────────────
    // Combine each library's original aligned BAM with its experiment's filtered GTF
    ch_aligned_bam
        .map { meta, bam, bai -> [ meta.experiment, meta, bam, bai ] }
        .combine(SQANTI3_FILTER.out.filtered_gtf, by: 0)
        .map { exp, meta, bam, bai, gtf -> [ meta, bam, bai, gtf ] }
        | ISOQUANT_QUANTIFY

    // ── Step 7: Export 10x-style QC Matrices enriched with SQANTI3 Metadata ──
    ISOQUANT_QUANTIFY.out.counts_dir
        .map { meta, counts -> [ meta.experiment, meta, counts ] }
        .combine(SQANTI3_FILTER.out.filtered_class, by: 0)
        .map { exp, meta, counts, class_file -> [ meta, counts, class_file ] }
        | EXPORT_LIBRARY_MTX

    // ── Step 8: Build Unified Cohort Shared Catalog ──────────────────────────
    SQANTI3_FILTER.out.filtered_gtf
        .join(SQANTI3_FILTER.out.filtered_class)
        | GENERATE_SHARED_CATALOG

    // ── Step 9: Calculate Sequencing Saturation Curves ───────────────────────
    ISOQUANT_QUANTIFY.out.counts_dir
        | CALCULATE_SATURATION

    emit:
    counts_dir      = ISOQUANT_QUANTIFY.out.counts_dir
    filtered_gtf    = SQANTI3_FILTER.out.filtered_gtf
    filtered_fasta  = SQANTI3_FILTER.out.filtered_fasta
    filtered_class  = SQANTI3_FILTER.out.filtered_class
    sqanti_class    = SQANTI3_QC.out.sqanti_results.map { exp, cls, fa, gtf -> cls }
    qc_export       = EXPORT_LIBRARY_MTX.out.mtx_export
    shared_catalog  = GENERATE_SHARED_CATALOG.out.catalog
    shared_map      = GENERATE_SHARED_CATALOG.out.map
    saturation      = CALCULATE_SATURATION.out.report
    versions        = Channel.empty()
                      .mix(CB_SUFFIX_INJECT.out.versions)
                      .mix(ISOQUANT_DISCOVERY.out.versions)
                      .mix(SQANTI3_QC.out.versions)
                      .mix(SQANTI3_FILTER.out.versions)
                      .mix(ISOQUANT_QUANTIFY.out.versions)
                      .mix(EXPORT_LIBRARY_MTX.out.versions)
                      .mix(GENERATE_SHARED_CATALOG.out.versions)
                      .mix(CALCULATE_SATURATION.out.versions)
}
