# Plan 03: IsoQuant Chromosome Sharding + SQANTI3 Chunking

**Goal:** Replace the single monolithic `ISOQUANT_DISCOVERY` job (32 CPUs / 600 GB / 48 h) with a scatter-gather across 17 chromosome-group shards running in parallel. Also shard SQANTI3 for the same reason. Wall time: ~48 h → ~4–6 h.

---

## Background

`ISOQUANT_DISCOVERY` receives ALL per-library aligned BAMs grouped by experiment — one sequential job regardless of sample count. With 6–12 libraries, it bottlenecks the entire pipeline.

**Why chromosome sharding is safe:** IsoQuant assembles transcripts locus-by-locus via splice graphs. Genes do not span chromosomes. 17 chromosome-group shards running in parallel produces an identical transcriptome to a single whole-genome run — verified as the canonical approach by both Perplexity and the wtsi-pbsc pipeline reference.

**Novel ID collision — critical:** IsoQuant assigns sequential novel IDs (`novel_gene_1`, `novel_transcript_1`) per run. Every shard independently starts from 1. **Naive GTF concatenation creates ID collisions that break SQANTI3 and downstream quantification.** The merge step must prefix novel IDs with chromosome name before concatenating. Identified by Gemini consultation.

**SQANTI3 also slow:** RajLabMSSM pipeline (`isoquant_pipeline.smk:98`) splits SQANTI3 into 8 parallel chunks. We add the same pattern.

**Do NOT split reference FASTA or GTF per shard.** Every shard loads the full reference; `--process_only_chr` scopes the processing. Splitting the reference causes indexing failures. Confirmed: Gemini.

---

## Chromosome Shard Groups

Grouped by GENCODE v49 protein-coding gene count — **not raw Mb.** Gene-dense chromosomes (chr1, chr2, chr11, chr17, chr19) are singletons. Gene deserts (chr13, chr18, chr21, chrY, chrM) are paired with mid-size chromosomes. chr19 is the most gene-dense chromosome (~1,471 genes, 59 Mb) and must be isolated.

| Shard | Chromosomes | ~Genes |
|-------|-------------|--------|
| 01 | chr1 | 1998 |
| 02 | chr19 | 1471 |
| 03 | chr11 | 1301 |
| 04 | chr2 | 1255 |
| 05 | chr17 | 1197 |
| 06 | chr5, chr13 | 1202 |
| 07 | chr6, chr21 | 1282 |
| 08 | chr10, chr22 | 1220 |
| 09 | chr4, chr20 | 1293 |
| 10 | chr14, chr8 | 1288 |
| 11 | chr3 | 1073 |
| 12 | chr16, chr18 | 1112 |
| 13 | chr12, chrM | 1047 |
| 14 | chrX, chrY | 904 |
| 15 | chr7 | 903 |
| 16 | chr9 | 781 |
| 17 | chr15 | 599 |

Range: 599–1998 genes per shard. Mean: 1172. chr1 is unavoidably ~2× mean — within-chromosome chunking (at zero-coverage gaps) would help further but is out of scope.

---

## Part A: IsoQuant Scatter-Gather

### Architecture

```
PBMM2_ALIGN.out.aligned_bam  [meta, bam, bai]  ← per library
        │
        │  .combine(ch_shards)  — cross-product: every library × every shard
        ▼
SPLIT_BAM_BY_SHARD             ← samtools view -b per (library × shard)
  [meta+shard, shard_bam, shard_bai]
        │
        │  .groupTuple by (experiment, shard_id)
        ▼
ISOQUANT_DISCOVERY_SHARD       ← 17 parallel jobs, all libs per shard
  [experiment, shard, shard_gtf]
        │
        │  .groupTuple by experiment
        ▼
MERGE_SHARD_GTFS               ← deduplicate headers + prefix novel IDs + dedup entries
  [experiment, merged_gtf]
        │
        ▼
SQANTI3_QC  → Part B sharding ↓
```

CB suffix injection (`CB_SUFFIX_INJECT`) still runs before splitting — suffixes survive `samtools view -b`.

---

### Phase 0: Confirm Novel ID Format

Before coding the merge script, confirm IsoQuant's exact novel ID format on the cluster:

```bash
find /sc-scratch/sc-scratch-cc12-ag-romagnani/nf_work_longshot \
    -path "*/ISOQUANT_DISCOVERY/*" -name "*.gtf" | head -1 \
    | xargs grep -m 20 'novel_'
```

Expected format: `gene_id "novel_gene_1"`, `transcript_id "novel_transcript_1.1"` — but confirm before coding the regex.

Also confirm `--process_only_chr` accepts space-separated list:
```bash
apptainer exec <isoquant_sif> isoquant --help 2>&1 | grep -A3 'process_only_chr'
```

---

### Phase 1: Shard Definition Channel

**File:** `subworkflows/classify.nf`

```groovy
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

Channel.from(isoquant_shards()).set { ch_shards }
```

---

### Phase 2: SPLIT_BAM_BY_SHARD Process

**File:** `modules/align.nf`

```nextflow
process SPLIT_BAM_BY_SHARD {
    tag "${meta.sample_id} shard${shard.id}"
    label 'process_low'
    container "${params.container_samtools}"

    input:
    tuple val(meta), path(bam), path(bai), val(shard)

    output:
    tuple val(meta), val(shard), path("${meta.sample_id}_shard${shard.id}.bam"),
                                 path("${meta.sample_id}_shard${shard.id}.bam.bai"), emit: shard_bam
    path "versions.yml", emit: versions

    script:
    def chr_list = shard.chrs.join(' ')
    """
    samtools view \\
        -@ ${task.cpus} \\
        -b \\
        -o ${meta.sample_id}_shard${shard.id}.bam \\
        ${bam} \\
        ${chr_list}

    samtools index -@ ${task.cpus} ${meta.sample_id}_shard${shard.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
```

**Wire in `subworkflows/classify.nf`:**
```groovy
CB_SUFFIX_INJECT.out.bam
    .combine(ch_shards)
    | SPLIT_BAM_BY_SHARD

SPLIT_BAM_BY_SHARD.out.shard_bam
    .map { meta, shard, bam, bai ->
        [ "${meta.experiment}__shard${shard.id}", shard, meta, bam, bai ]
    }
    .groupTuple(by: 0)
    .map { key, shards, metas, bams, bais ->
        [ metas[0].experiment, shards[0], metas, bams, bais ]
    }
    .set { ch_sharded_bams }
```

---

### Phase 3: ISOQUANT_DISCOVERY_SHARD Process

**File:** `modules/quantify.nf`

Keep old `ISOQUANT_DISCOVERY` commented until validation is complete.

```nextflow
process ISOQUANT_DISCOVERY_SHARD {
    tag "${experiment} shard${shard.id} [${shard.chrs.join(',')}]"
    label 'process_high'
    // ponytail: downgraded from process_ultra — each shard is ~1/17th the data
    container "${params.container_isoquant}"

    input:
    tuple val(experiment), val(shard), val(metas), path(bams), path(bais)

    output:
    tuple val(experiment), val(shard),
          path("isoquant_out/${experiment}_shard${shard.id}/${experiment}_shard${shard.id}.transcript_models.gtf"),
          emit: shard_gtf
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(isoquant --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1)
    END_VERSIONS
    """
}
```

---

### Phase 4: MERGE_SHARD_GTFS Process

**File:** `modules/quantify.nf`

Critical: deduplicate headers, prefix novel IDs with chromosome name, then remove any duplicate GTF entries (safety net per RajLabMSSM pattern).

```nextflow
process MERGE_SHARD_GTFS {
    tag "${experiment}"
    label 'process_low'
    container "${params.container_multiqc}"
    publishDir { "${params.outdir}/${experiment}/joint/transcript_model" }, mode: 'copy'

    input:
    tuple val(experiment), val(shards), path(gtfs)

    output:
    tuple val(experiment), path("${experiment}.transcript_models.gtf"), emit: transcript_gtf
    path "versions.yml", emit: versions

    script:
    """
    python3 -c "
import re, sys

gtf_files   = sorted('${gtfs}'.split())
out_file    = '${experiment}.transcript_models.gtf'
header_seen = set()
entry_seen  = set()   # dedup safety net (RajLabMSSM pattern)

with open(out_file, 'w') as out:
    for gtf in gtf_files:
        with open(gtf) as f:
            for line in f:
                if line.startswith('#'):
                    if line not in header_seen:
                        header_seen.add(line)
                        out.write(line)
                    continue

                chrom = line.split('\t')[0]

                # Prefix novel gene/transcript IDs with chromosome name
                line = re.sub(r'(gene_id \\\"novel_gene_)',       r'\g<1>' + chrom + r'_', line)
                line = re.sub(r'(transcript_id \\\"novel_transcript_)', r'\g<1>' + chrom + r'_', line)
                line = re.sub(r'(gene_name \\\"novel_gene_)',      r'\g<1>' + chrom + r'_', line)

                # Dedup identical entries (safety net)
                key = line.strip()
                if key not in entry_seen:
                    entry_seen.add(key)
                    out.write(line)
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | grep -oP '(?<=Python )\\S+')
    END_VERSIONS
    """
}
```

**Wire grouping in `subworkflows/classify.nf`:**
```groovy
ISOQUANT_DISCOVERY_SHARD.out.shard_gtf
    .map { experiment, shard, gtf -> [ experiment, shard, gtf ] }
    .groupTuple(by: 0)
    | MERGE_SHARD_GTFS
```

---

## Part B: SQANTI3 Chunking

RajLabMSSM (`isoquant_pipeline.smk:98`) splits SQANTI3 into 8 parallel chunks. SQANTI3 can be parallelised by GTF chunk — each chunk processes a subset of transcripts independently; classification tables are concatenated.

### Architecture

```
MERGE_SHARD_GTFS.out.transcript_gtf  [experiment, gtf]
        │
        │  split GTF into N_SQANTI_CHUNKS subsets by transcript
        ▼
SQANTI3_SPLIT_GTF                    ← N chunks per experiment
  [experiment, chunk_id, chunk_gtf]
        │
        ▼
SQANTI3_QC_CHUNK                     ← N parallel SQANTI3 jobs
  [experiment, chunk_id, classification, corrected_fasta, corrected_gtf]
        │
        │  .groupTuple by experiment
        ▼
SQANTI3_MERGE_CHUNKS                 ← cat classification tables, merge GTFs/FASTAs
  [experiment, classification, corrected_fasta, corrected_gtf]
        │
        ▼
SQANTI3_FILTER  (unchanged)
```

Default `params.sqanti_chunks = 8`.

### Phase 5: SQANTI3_SPLIT_GTF Process

**File:** `modules/quantify.nf`

```nextflow
process SQANTI3_SPLIT_GTF {
    tag "${experiment} (${params.sqanti_chunks} chunks)"
    label 'process_low'
    container "${params.container_multiqc}"

    input:
    tuple val(experiment), path(gtf)

    output:
    tuple val(experiment), path("chunk_*.gtf"), emit: chunks
    path "versions.yml",                         emit: versions

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
        if '\t' not in line:
            continue
        fields = line.split('\t')
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | grep -oP '(?<=Python )\\S+')
    END_VERSIONS
    """
}
```

### Phase 6: SQANTI3_QC_CHUNK Process

**File:** `modules/quantify.nf`

Refactor existing `SQANTI3_QC` to accept a single chunk GTF instead of the full GTF. Tag with chunk ID.

```nextflow
process SQANTI3_QC_CHUNK {
    tag "${experiment} chunk${chunk_id}"
    label 'process_medium'   // ponytail: downgraded from process_high — 1/8th the transcripts
    container "${params.container_sqanti3}"
    ...
}
```

### Phase 7: SQANTI3_MERGE_CHUNKS Process

**File:** `modules/quantify.nf`

Merge chunk classification TSVs and corrected GTFs:
- Classification: `head -1` from first chunk for header, then `tail -n+2` from all chunks
- Corrected GTF: concatenate (same pattern as MERGE_SHARD_GTFS, but no novel ID prefixing needed — SQANTI3 doesn't create novel IDs)
- Corrected FASTA: concatenate

---

## Phase 8: Resource Overrides in `nextflow.config`

```groovy
withName: 'SPLIT_BAM_BY_SHARD' {
    cpus   = 4
    memory = 8.GB
    time   = 2.h
}
withName: 'ISOQUANT_DISCOVERY_SHARD' {
    cpus   = 16
    memory = 128.GB   // monitor chr1 shard; increase to 256 GB if OOM
    time   = 12.h
}
withName: 'MERGE_SHARD_GTFS' {
    cpus   = 2
    memory = 16.GB
    time   = 1.h
}
withName: 'SQANTI3_SPLIT_GTF' {
    cpus   = 2
    memory = 8.GB
    time   = 30.m
}
withName: 'SQANTI3_QC_CHUNK' {
    cpus   = 8
    memory = 64.GB
    time   = 6.h
}
withName: 'SQANTI3_MERGE_CHUNKS' {
    cpus   = 2
    memory = 16.GB
    time   = 30.m
}
```

Add to `params`: `sqanti_chunks = 8`

---

## Phase 9: Validation

### IsoQuant shard validation
```bash
# Novel ID uniqueness — expect 0 duplicates
grep -v '^#' experiment.transcript_models.gtf \
    | awk '$3=="transcript"' \
    | grep -oP 'transcript_id "[^"]+"' \
    | sort | uniq -d | wc -l

# Novel IDs prefixed correctly
grep 'novel_transcript' experiment.transcript_models.gtf | head -5
# Expect: novel_transcript_chr1_1.1, novel_transcript_chr19_3.2, etc.

# All 25 target chromosomes present
grep -v '^#' experiment.transcript_models.gtf | cut -f1 | sort -u

# Transcript count vs monolithic run (run both on 2-lib test case)
grep -v '^#' monolithic.gtf | awk '$3=="transcript"' | wc -l
grep -v '^#' sharded_merged.gtf | awk '$3=="transcript"' | wc -l
```

### SQANTI3 chunk validation
```bash
# Row count in merged classification == sum of chunk row counts
wc -l merged_classification.txt
```

---

## Files Changed Summary

| File | Change |
|------|--------|
| `modules/align.nf` | New `SPLIT_BAM_BY_SHARD` |
| `modules/quantify.nf` | New `ISOQUANT_DISCOVERY_SHARD`, `MERGE_SHARD_GTFS`, `SQANTI3_SPLIT_GTF`, `SQANTI3_QC_CHUNK`, `SQANTI3_MERGE_CHUNKS`; old `ISOQUANT_DISCOVERY` + `SQANTI3_QC` kept/commented until validated |
| `subworkflows/classify.nf` | Shard channel def; scatter-gather wiring for both IsoQuant and SQANTI3 |
| `nextflow.config` | New per-process resource overrides; `params.sqanti_chunks = 8` |

---

## Known Risks

| Risk | Mitigation |
|------|-----------|
| chr1 shard OOM at 128 GB | Monitor first run; increase to 256 GB if needed |
| Novel ID regex misses IsoQuant format variant | Phase 0 confirms exact format before coding merge script |
| Empty shard GTF for chrY / chrM | `optional: true` on shard GTF outputs; merge script skips empty files |
| SQANTI3 chunk split orphans exon lines from previous transcript | Splitter groups by transcript block (transcript + its exons) — safe |
| Inter-chromosomal fusion transcripts missed | Accepted; MAS-ISO-seq chimeras are a noise source, not signal |
