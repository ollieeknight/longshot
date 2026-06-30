// GRCh38 chromosome-shard partition used to parallelize joint IsoQuant
// discovery in CLASSIFY. Shards are grouped to roughly balance read volume
// per shard (large chromosomes alone, small ones paired), not by chromosome
// number — chr1 dominates a shard on its own, chr5+chr13 together approximate
// the same load. Re-balance manually if porting to a non-human genome.
def chromosomeShardMap() {
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
