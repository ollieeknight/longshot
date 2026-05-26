include { EXPORT_LIBRARY_MTX; GENERATE_SHARED_CATALOG; CALCULATE_SATURATION } from '../modules/exporter'

workflow EXPORT {
    take:
    ch_counts_dir     // [meta, counts_dir] — per library
    ch_filtered_class // [experiment, classification] — per experiment
    ch_filtered_gtf   // [experiment, gtf] — per experiment

    main:
    // ── Step 1: Export 10x-style QC Matrices enriched with SQANTI3 Metadata ──
    ch_counts_dir
        .map { meta, counts -> [ meta.experiment, meta, counts ] }
        .combine(ch_filtered_class, by: 0)
        .map { exp, meta, counts, class_file -> [ meta, counts, class_file ] }
        | EXPORT_LIBRARY_MTX

    // ── Step 2: Build Unified Cohort Shared Catalog ──────────────────────────
    ch_filtered_gtf
        .join(ch_filtered_class)
        | GENERATE_SHARED_CATALOG

    // ── Step 3: Calculate Sequencing Saturation Curves ───────────────────────
    ch_counts_dir
        | CALCULATE_SATURATION

    emit:
    qc_export      = EXPORT_LIBRARY_MTX.out.mtx_export
    shared_catalog = GENERATE_SHARED_CATALOG.out.catalog
    shared_map     = GENERATE_SHARED_CATALOG.out.map
    saturation     = CALCULATE_SATURATION.out.report
    versions       = Channel.empty()
                     .mix(EXPORT_LIBRARY_MTX.out.versions)
                     .mix(GENERATE_SHARED_CATALOG.out.versions)
                     .mix(CALCULATE_SATURATION.out.versions)
}
