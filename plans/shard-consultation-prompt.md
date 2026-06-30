# Consultation Prompt: IsoQuant Sharding Strategy for PacBio scRNA-seq Pipeline

## Context

We have a Nextflow DSL2 pipeline (`longshot`) for PacBio MAS-ISO-seq single-cell long-read RNA-seq. It processes HiFi BAMs through adapter splitting (skera), primer trimming (lima), barcode tagging and deduplication (isoseq), alignment (pbmm2), and joint isoform discovery + quantification (IsoQuant).

The bottleneck is **joint isoform discovery** (`ISOQUANT_DISCOVERY`). This step takes all per-library aligned BAMs for an experiment and runs IsoQuant once across all of them to build a consensus transcript model. With 6+ libraries:

- Single job: 32 CPUs, 600 GB RAM, 48 h walltime
- All library BAMs are passed simultaneously: `isoquant --bam lib1.bam lib2.bam ... libN.bam`
- IsoQuant already uses `--process_only_chr chr1 chr2 ... chr22 chrX chrY chrM` to exclude alt contigs
- Output: one `experiment.transcript_models.gtf` used by SQANTI3 downstream

The joint modeling across libraries is required — IsoQuant uses all samples together to build a consensus isoform catalog. This cannot be parallelised per-library.

## The Problem

With 6–12 libraries, this single job becomes the wall-time and memory bottleneck. We want to shard the discovery step to run in parallel.

## Hard Constraint

**IsoQuant assembles transcripts locus by locus.** All reads from the same gene must be in the same shard. You cannot split reads mid-chromosome. Valid shard boundaries are chromosome boundaries only.

## Proposed Approach

**Scatter-gather by chromosome groups.**

1. **Pre-split** each aligned BAM by chromosome group using `samtools view -b`. Each library produces N shard BAMs.
2. **Run IsoQuant per shard** — each shard job receives all libraries' shard BAMs for that chromosome group, plus `--process_only_chr` scoped to just that group's chromosomes.
3. **Merge GTFs** — concatenate the N per-shard GTF outputs into a single experiment GTF (GTF headers deduplicated, chromosome-scoped genes have no cross-shard references).
4. **Continue unchanged** — SQANTI3 receives the merged GTF as before.

**Why group chromosomes rather than use 25 individual chromosomes?**
- chr1 alone is ~249 Mb (~8% of the target genome); chrY + chrM together are <0.1%
- 25 individual shards = highly unequal job sizes; chr1 still dominates wall time
- Grouping ~25 chromosomes into ~25–30 roughly equal-sized buckets by genomic length balances wall time across shards
- Chromosome sizes in hg38 are fixed and known — groups can be hardcoded, no runtime read counting

**Proposed groups (hg38, target chromosomes only: chr1–22, X, Y, M):**

Approximate sizes (Mb): chr1=249, chr2=242, chr3=198, chr4=190, chr5=182, chr6=171, chr7=159, chr8=145, chr9=138, chr10=133, chr11=135, chr12=133, chr13=114, chr14=107, chr15=102, chr16=90, chr17=83, chr18=80, chr19=59, chr20=64, chr21=47, chr22=51, chrX=156, chrY=57, chrM=0.016

Total target genome: ~3,088 Mb. With 25 shards, target ~123 Mb per shard.

One reasonable grouping (25 shards):
- Shard 1: chr1 (249 Mb) — largest, unavoidably so
- Shard 2: chr2 (242 Mb)
- Shard 3: chr3 (198 Mb)
- Shard 4: chr4 (190 Mb)
- Shard 5: chr5 (182 Mb)
- Shard 6: chr6 (171 Mb)
- Shard 7: chr7 (159 Mb)
- Shard 8: chrX (156 Mb)
- Shard 9: chr8 (145 Mb)
- Shard 10: chr9 + chrY (138 + 57 = 195 Mb)
- Shard 11: chr10 (133 Mb)
- Shard 12: chr11 (135 Mb)
- Shard 13: chr12 (133 Mb)
- Shard 14: chr13 (114 Mb)
- Shard 15: chr14 (107 Mb)
- Shard 16: chr15 (102 Mb)
- Shard 17: chr16 + chr21 (90 + 47 = 137 Mb)
- Shard 18: chr17 + chr22 (83 + 51 = 134 Mb)
- Shard 19: chr18 + chr19 (80 + 59 = 139 Mb)
- Shard 20: chr20 + chrM (64 + 0.016 = 64 Mb)

This gives ~20 shards of 64–249 Mb. The top 9 are single chromosomes; smaller ones are paired.

**Expected improvement:**
- Current: 1 × 48 h × 600 GB × 32 CPUs
- Sharded: ~20 parallel jobs × ~4–6 h × ~50 GB × 16 CPUs
- Wall time: 48 h → ~6 h (dominated by chr1/chr2 shards)

## What We Want Feedback On

1. **Is chromosome-grouped sharding the right strategy, or is there a better approach** (e.g. IsoQuant-native parallelism, a different tool, a different split axis)?

2. **Are the proposed chromosome groups sensible**, or should the grouping be done differently (e.g. by gene density rather than raw base pairs, since read coverage tracks gene density better than genomic length)?

3. **GTF merge correctness**: After sharded IsoQuant runs, we concatenate the per-shard GTFs. IsoQuant assigns transcript IDs like `experiment.transcript_models_ENSGXXX.1`. Is there any risk of ID collision or inconsistency across shards that a simple concatenation wouldn't handle?

4. **BAM pre-splitting overhead**: With 6–12 libraries × 20 shards = 120–240 `samtools view -b` jobs. Each is lightweight but adds pipeline steps. Is there a reason NOT to pre-split and instead just pass full BAMs to each IsoQuant shard with `--process_only_chr` scoped to the group? (Concern: IsoQuant would then read all 6–12 full-genome BAMs 20 times, which is massive redundant I/O.)

5. **Any IsoQuant-specific gotchas** when running in single-chromosome or chromosome-group mode that differ from whole-genome mode (e.g. junction evidence from reads spanning chromosome-boundary-adjacent genes, edge cases in the `--complete_genedb` flag behaviour)?

## Pipeline Stack for Reference

- Nextflow DSL2, SLURM + Apptainer, BIH HPC
- IsoQuant 3.13.0 (`quay.io/biocontainers/isoquant:3.13.0--pyh106432d_0`)
- Reads: PacBio HiFi, MAS-ISO-seq (16-molecule concatenates), 10x 3' or 5' GEX barcoded
- Reference: GRCh38 no-alt, GENCODE v49
- Typical cohort: 6–12 libraries per experiment, each library ~5–20 M FLTNC reads post-dedup
