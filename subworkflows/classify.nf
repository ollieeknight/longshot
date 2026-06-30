include { CB_SUFFIX_INJECT           } from '../modules/align'
include { SPLIT_BAM_BY_SHARD         } from '../modules/align'
include { ISOQUANT_DISCOVERY_SHARD   } from '../modules/quantify'
include { MERGE_SHARD_GTFS           } from '../modules/quantify'
include { SQANTI3_SPLIT_GTF          } from '../modules/quantify'
include { SQANTI3_QC_CHUNK           } from '../modules/quantify'
include { SQANTI3_MERGE_CHUNKS       } from '../modules/quantify'
include { SQANTI3_FILTER             } from '../modules/quantify'
include { SQANTI3_RESCUE             } from '../modules/quantify'

def isoquant_shards() {
    return [
        [id: '01', chrs: ['chr1']],
        [id: '02', chrs: ['chr19']],
        [id: '03', chrs: ['chr11']],
        [id: '04', chrs: ['chr2']],
        [id: '05', chrs: ['chr17']],
        [id: '06', chrs: ['chr5',  'chr13']],
        [id: '07', chrs: ['chr6',  'chr21']],
        [id: '08', chrs: ['chr10', 'chr22']],
        [id: '09', chrs: ['chr4',  'chr20']],
        [id: '10', chrs: ['chr14', 'chr8']],
        [id: '11', chrs: ['chr3']],
        [id: '12', chrs: ['chr16', 'chr18']],
        [id: '13', chrs: ['chr12', 'chrM']],
        [id: '14', chrs: ['chrX',  'chrY']],
        [id: '15', chrs: ['chr7']],
        [id: '16', chrs: ['chr9']],
        [id: '17', chrs: ['chr15']],
    ]
}

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

    // ── Step 3: Cross-product every library BAM with every chromosome shard ───
    Channel.from(isoquant_shards()).set { ch_shards }

    CB_SUFFIX_INJECT.out.bam
        .combine(ch_shards)
        // now: [meta, bam, bai, shard]
        | SPLIT_BAM_BY_SHARD

    // ── Step 4: Group shard BAMs by (experiment, shard) for joint discovery ──
    SPLIT_BAM_BY_SHARD.out.shard_bam
        .map { meta, shard, bam, bai ->
            [ "${meta.experiment}__shard${shard.id}", meta.experiment, shard, bam, bai ]
        }
        .groupTuple(by: 0)
        .map { key, experiments, shards, bams, bais ->
            [ experiments[0], shards[0], bams, bais ]
        }
        | ISOQUANT_DISCOVERY_SHARD

    // ── Step 5: Gather all shard GTFs per experiment and merge ────────────────
    ISOQUANT_DISCOVERY_SHARD.out.shard_gtf
        .groupTuple(by: 0)
        // now: [experiment, [shard,...], [gtf,...]]
        | MERGE_SHARD_GTFS

    // ── Step 6: Split merged GTF into SQANTI3 chunks ──────────────────────────
    MERGE_SHARD_GTFS.out.transcript_gtf
        | SQANTI3_SPLIT_GTF

    // Flatten chunk list into individual [experiment, chunk_id, gtf] tuples
    SQANTI3_SPLIT_GTF.out.chunks
        .flatMap { experiment, chunk_gtfs ->
            def chunks = chunk_gtfs instanceof List ? chunk_gtfs : [chunk_gtfs]
            chunks.collect { gtf ->
                def chunk_id = gtf.name.replaceAll(/chunk_(\d+)\.gtf/, '$1')
                [ experiment, chunk_id, gtf ]
            }
        }
        | SQANTI3_QC_CHUNK

    // ── Step 7: Merge SQANTI3 chunk results per experiment ───────────────────
    SQANTI3_QC_CHUNK.out.chunk_results
        .groupTuple(by: 0)
        // now: [experiment, [chunk_ids], [cls_files], [fasta_files], [gtf_files]]
        | SQANTI3_MERGE_CHUNKS

    // ── Step 8: SQANTI3 Filter + Rescue (unchanged) ──────────────────────────
    SQANTI3_FILTER(SQANTI3_MERGE_CHUNKS.out.sqanti_results)

    SQANTI3_FILTER.out.filtered_class
        .join(SQANTI3_MERGE_CHUNKS.out.sqanti_results.map { exp, cls, fa, gtf -> [exp, fa] })
        .join(SQANTI3_FILTER.out.filtered_gtf)
        .map { exp, filter_cls, corr_fa, filt_gtf -> [exp, filter_cls, corr_fa, filt_gtf] }
        | SQANTI3_RESCUE

    emit:
    rescued_gtf     = SQANTI3_RESCUE.out.rescued_gtf
    rescued_fasta   = SQANTI3_RESCUE.out.rescued_fasta
    filtered_gtf    = SQANTI3_FILTER.out.filtered_gtf
    filtered_fasta  = SQANTI3_FILTER.out.filtered_fasta
    filtered_class  = SQANTI3_FILTER.out.filtered_class
    sqanti_class    = SQANTI3_MERGE_CHUNKS.out.sqanti_results.map { exp, cls, fa, gtf -> cls }
    versions        = Channel.empty()
                      .mix(CB_SUFFIX_INJECT.out.versions)
                      .mix(SPLIT_BAM_BY_SHARD.out.versions)
                      .mix(ISOQUANT_DISCOVERY_SHARD.out.versions)
                      .mix(MERGE_SHARD_GTFS.out.versions)
                      .mix(SQANTI3_SPLIT_GTF.out.versions)
                      .mix(SQANTI3_QC_CHUNK.out.versions)
                      .mix(SQANTI3_MERGE_CHUNKS.out.versions)
                      .mix(SQANTI3_FILTER.out.versions)
                      .mix(SQANTI3_RESCUE.out.versions)
}
