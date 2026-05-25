include { SKERA_SPLIT; DETECT_SAMPLE_INDICES; CONSTRUCT_MULTIPLEX_PRIMERS; LIMA_ISOSEQ; LIMA_MULTIPLEX; MERGE_INDEX_BAMS; ISOSEQ_TAG; ISOSEQ_REFINE } from '../modules/preprocess'
include { SAMTOOLS_FLAGSTAT } from '../modules/qc'


workflow PREPROCESS {
    take:
    ch_input  // [meta, bam]  — one entry per samplesheet row (per lane)

    main:
    // ── Group inputs by unique raw BAM path ──────────────────────────────────
    ch_input
        .map { meta, bam -> [ bam, meta ] }
        .groupTuple(by: 0)
        .map { bam, metas ->
            def first_meta = metas[0]
            def smrt_meta = [
                id:                  first_meta.run_id,
                run_id:              first_meta.run_id,
                experiment:          first_meta.experiment
            ]
            [ smrt_meta, bam, metas ]
        }
        .set { ch_grouped_smrtcells }

    // ── Run Skera split once per raw BAM flowcell ─────────────────────────────
    SKERA_SPLIT(ch_grouped_smrtcells.map { smrt_meta, bam, metas -> [ smrt_meta, bam ] })

    // ── Run pre-flight sample index detection QC ─────────────────────────────
    ch_grouped_smrtcells
        .map { smrt_meta, bam, metas ->
            def idx_str = metas.collect { it.tenx_index }.findAll { it != null }.join(',')
            [ smrt_meta, bam, idx_str ]
        }
        | DETECT_SAMPLE_INDICES

    // ── Join split BAM with grouped metadata to branch on demultiplexing ──────
    SKERA_SPLIT.out.segmented_bam
        .join(ch_grouped_smrtcells.map { smrt_meta, bam, metas -> [ smrt_meta, metas ] })
        .branch { smrt_meta, segmented_bam, metas ->
            def unique_indices = metas.collect { it.tenx_index }.findAll { it != null }.unique()
            multiplexed: unique_indices.size() > 1
            standard: true
        }
        .set { ch_branched_segmented }

    // ── Standard single-index or non-multiplexed branch ──────────────────────
    ch_branched_segmented.standard
        .map { smrt_meta, segmented_bam, metas ->
            def meta = metas[0]
            // Auto-detect 5' vs 3' kit from index name prefix (e.g. SI-GA maps to 5' kit)
            def primers = (meta.tenx_index && meta.tenx_index.startsWith("SI-GA")) ? file(params.tenx_5kit_primers) : file(params.tenx_3kit_primers)
            [ meta, segmented_bam, primers ]
        }
        | LIMA_ISOSEQ

    // ── Multiplexed index demultiplexing branch ──────────────────────────────
    // 1. Construct index-specific primers FASTA
    ch_branched_segmented.multiplexed
        .map { smrt_meta, segmented_bam, metas ->
            def mapping = metas.collect { "${it.library_id}:${it.tenx_index}" }.join(',')
            [ smrt_meta, mapping ]
        }
        | CONSTRUCT_MULTIPLEX_PRIMERS

    // 2. Run LIMA multiplexed demultiplexing
    ch_branched_segmented.multiplexed
        .map { smrt_meta, segmented_bam, metas -> [ smrt_meta, segmented_bam ] }
        .join(CONSTRUCT_MULTIPLEX_PRIMERS.out.primers)
        | LIMA_MULTIPLEX

    // 3. Flatten split BAM files and merge them back per individual library ID
    LIMA_MULTIPLEX.out.split_bams
        .join(ch_branched_segmented.multiplexed.map { smrt_meta, segmented_bam, metas -> [ smrt_meta, metas ] })
        .flatMap { smrt_meta, bams, metas ->
            metas.collect { meta ->
                def lib_bams = bams.findAll { it.name.contains(meta.library_id) }
                [ meta, lib_bams ]
            }
        }
        | MERGE_INDEX_BAMS

    // ── Mix standard and demultiplexed/merged BAMs together ──────────────────
    LIMA_ISOSEQ.out.fl_bam
        .mix(MERGE_INDEX_BAMS.out.fl_bam)
        .set { ch_fl_bam }

    // ── downstream tagging & refining per library ID ─────────────────────────
    ISOSEQ_TAG(ch_fl_bam)
    ISOSEQ_REFINE(ISOSEQ_TAG.out.flt_bam)

    // QC: flagstat on FLTNC (post-refine read count)
    ISOSEQ_REFINE.out.fltnc_bam
        .map { meta, bam -> [ meta, 'fltnc', bam ] }
        | SAMTOOLS_FLAGSTAT

    emit:
    fltnc_bam    = ISOSEQ_REFINE.out.fltnc_bam
    lima_reports = LIMA_ISOSEQ.out.lima_reports.mix(LIMA_MULTIPLEX.out.lima_reports)
    flagstat     = SAMTOOLS_FLAGSTAT.out.flagstat
    versions     = Channel.empty()
                    .mix(SKERA_SPLIT.out.versions)
                    .mix(LIMA_ISOSEQ.out.versions)
                    .mix(LIMA_MULTIPLEX.out.versions)
                    .mix(ISOSEQ_TAG.out.versions)
                    .mix(ISOSEQ_REFINE.out.versions)
}
