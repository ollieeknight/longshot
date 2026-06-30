include { SAMTOOLS_MERGE_FLTNC } from '../modules/align'
include { PREPARE_WHITELIST      } from '../modules/align'
include { ISOSEQ_CORRECT       } from '../modules/align'
include { SAMTOOLS_SORT_CB     } from '../modules/align'
include { ISOSEQ_GROUPDEDUP    } from '../modules/align'
include { PBMM2_ALIGN          } from '../modules/align'
include { GENERATE_CRAM        } from '../modules/align'
include { SAMTOOLS_FLAGSTAT; CRAMINO as CRAMINO_ALIGNED } from '../modules/qc'


workflow BARCODE_ALIGN {
    take:
    ch_fltnc_bam  // [meta, bam] — per lane, meta has sample_id for grouping

    main:
    // Merge FLTNC BAMs across SMRT cells for the same library
    ch_fltnc_bam
        .map { meta, bam -> [ meta.sample_id, meta, bam ] }
        .groupTuple(by: 0)
        .map { sample_id, metas, bams ->
            // Use the first meta; all metas for same sample_id are identical except run_id
            [ metas[0], bams ]
        }
        | SAMTOOLS_MERGE_FLTNC

    // Separate libraries with short-read barcodes from those without
    SAMTOOLS_MERGE_FLTNC.out.merged_bam
        .branch { meta, bam ->
            custom: meta.shortread_barcodes != null
            standard: true
        }
        .set { ch_barcode_source_branch }

    // Custom: Prepare the custom reverse-complemented whitelist
    ch_barcode_source_branch.custom
        .map { meta, bam -> [ meta, file(meta.shortread_barcodes) ] }
        | PREPARE_WHITELIST

    // Combine standard BAMs with global static whitelist
    ch_barcode_source_branch.standard
        .map { meta, bam -> [ meta, bam, file(params.tenx_whitelist) ] }
        .set { ch_standard_correct_input }

    // Combine custom BAMs with their custom prepared whitelist
    ch_barcode_source_branch.custom
        .join(PREPARE_WHITELIST.out.whitelist)
        .map { meta, bam, whitelist -> [ meta, bam, whitelist ] }
        .set { ch_custom_correct_input }

    // Mix standard and custom inputs for barcode correction
    ch_standard_correct_input
        .mix(ch_custom_correct_input)
        .set { ch_correct_input }

    ISOSEQ_CORRECT(ch_correct_input)
    SAMTOOLS_SORT_CB(ISOSEQ_CORRECT.out.corrected_bam)
    ISOSEQ_GROUPDEDUP(SAMTOOLS_SORT_CB.out.sorted_bam)

    PBMM2_ALIGN(ISOSEQ_GROUPDEDUP.out.dedup_bam)

    // QC: flagstat post-dedup & post-alignment
    ISOSEQ_GROUPDEDUP.out.dedup_bam
        .map { meta, bam -> [ meta, 'dedup', bam ] }
        .mix(
            PBMM2_ALIGN.out.aligned_bam
                .map { meta, bam, bai -> [ meta, 'aligned', bam ] }
        )
        | SAMTOOLS_FLAGSTAT

    // cramino: merged FLTNC, dedup, and aligned stages
    SAMTOOLS_MERGE_FLTNC.out.merged_bam
        .map { meta, bam -> [ meta, 'merged', '--ubam', bam ] }
        .mix( ISOSEQ_GROUPDEDUP.out.dedup_bam
                .map { meta, bam -> [ meta, 'dedup', '--ubam', bam ] } )
        .mix( PBMM2_ALIGN.out.aligned_bam
                .map { meta, bam, bai -> [ meta, 'aligned', '--spliced --karyotype', bam ] } )
        | CRAMINO_ALIGNED

    // Archive as CRAM
    GENERATE_CRAM(PBMM2_ALIGN.out.aligned_bam)

    emit:
    aligned_bam  = PBMM2_ALIGN.out.aligned_bam
    cram         = GENERATE_CRAM.out.cram
    flagstat     = SAMTOOLS_FLAGSTAT.out.flagstat
    cramino      = CRAMINO_ALIGNED.out.stats
}
