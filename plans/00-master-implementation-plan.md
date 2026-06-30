# Master Implementation Plan: longshot Pipeline Improvements

**For use at the start of a new implementation session.**
Execute plans in order: Plan 02 first (independent, lower risk), then Plan 03 (depends on IsoQuant module from 02).

---

## Context Summary

`longshot` is a Nextflow DSL2 PacBio MAS-ISO-seq single-cell pipeline running on BIH HPC (SLURM + Apptainer). It processes HiFi BAMs through: skera → lima → isoseq tag/refine/correct/groupdedup → pbmm2 → IsoQuant joint discovery → SQANTI3 → IsoQuant quantification → export.

**Input:** `hifi_reads/*.bcM0001.bam` from Revio — post-CCS, post-sample-index-lima. This is the true upstream limit; no subreads or pre-lima BAMs are available from Revio.

**Key files:**
- `main.nf` — entry point and channel orchestration
- `subworkflows/preprocess.nf` — skera → lima → isoseq tag/refine
- `subworkflows/align.nf` — isoseq correct/dedup → pbmm2
- `subworkflows/classify.nf` — CB suffix injection → IsoQuant → SQANTI3
- `modules/preprocess.nf`, `modules/align.nf`, `modules/quantify.nf`, `modules/qc.nf`
- `nextflow.config` — SLURM resource labels and container paths

**Dry-run command (syntax check, no cluster jobs):**
```bash
nextflow run main.nf -preview -profile standard --samplesheet assets/example_samplesheet.csv
```

---

## Plan 02: QC Gap Fixes

**File:** `plans/02-qc-gap-fixes.md`
**Risk:** Low. All changes are additive (new outputs, new publishDir). No logic changes.
**Estimated effort:** 1 session.

### What changes

| Change | File(s) | Why |
|--------|---------|-----|
| `--dump-clips --dump-removed` on LIMA | `modules/preprocess.nf` | Captures clipped/rejected reads for QC. RajLabMSSM pattern. |
| Capture `isoseq refine` filter summary | `modules/preprocess.nf` | Shows poly-A / chimera attrition — currently invisible |
| Capture `isoseq correct` stats | `modules/preprocess.nf` | Shows barcode match rate — currently invisible |
| Capture skera split stats | `modules/preprocess.nf` | Shows segmentation yield — currently invisible |
| `publishDir` on `SAMTOOLS_FLAGSTAT` | `modules/qc.nf` | Flagstat only in MultiQC now; R can't read it |
| `publishDir` on lima reports | `modules/preprocess.nf` | Same |
| `publishDir` on instrument stats | `modules/qc.nf` | Same |
| `--fl_data --check_canonical --count_exons` on IsoQuant | `modules/quantify.nf` | `--fl_data` is correct for post-dedup full-length reads; `--check_canonical` filters artefactual junctions; `--count_exons` adds exon-level output. Dev_Brain_IsoSeq pattern. |

### Before starting Phase 1–3

Run Phase 0 on the cluster to confirm exact sidecar filenames:
```bash
find /sc-scratch/sc-scratch-cc12-ag-romagnani/nf_work_longshot \
    -path "*/ISOSEQ_REFINE/*" \( -name "*.filter_summary*" -o -name "*.report.csv" \)
find /sc-scratch/sc-scratch-cc12-ag-romagnani/nf_work_longshot \
    -path "*/ISOSEQ_CORRECT/*" \( -name "*.json" -o -name "*.csv" -o -name "*.log" \)
find /sc-scratch/sc-scratch-cc12-ag-romagnani/nf_work_longshot \
    -path "*/SKERA_SPLIT/*" -name "*report*" -o -path "*/SKERA_SPLIT/*" -name "*.log"
```

---

## Plan 03: IsoQuant Chromosome Sharding + SQANTI3 Chunking

**File:** `plans/03-isoquant-chromosome-sharding.md`
**Risk:** Medium. Replaces the most resource-intensive process. Requires validation run before production use.
**Estimated effort:** 2 sessions (implementation + validation).
**Depends on:** Plan 02 IsoQuant flag additions (add `--fl_data --check_canonical --count_exons` to `ISOQUANT_DISCOVERY_SHARD`).

### The problem

`ISOQUANT_DISCOVERY`: 1 job × 32 CPUs × 600 GB × 48 h for all libraries. Becomes the wall-time bottleneck at 6+ libraries.

### The solution

**Part A — IsoQuant scatter-gather:**
Split into 17 parallel chromosome-group shards. Each shard receives all library BAMs but scoped to its chromosomes. Merge 17 GTFs into one. **Expected: 48 h → ~6 h wall time.**

**Part B — SQANTI3 chunking:**
Split transcript GTF into 8 chunks, run SQANTI3 in parallel, merge classification tables. Same scatter-gather pattern.

### Critical: Novel ID collision

IsoQuant assigns sequential novel IDs (`novel_gene_1`, `novel_transcript_1`) per run. Every shard starts from 1. **Naive concatenation breaks SQANTI3 and downstream quantification.** The `MERGE_SHARD_GTFS` process prefixes novel IDs with the chromosome name from GTF column 1 before concatenating.

This is not a problem in the current monolithic run, nor in the reference pipelines (RajLabMSSM, Dev_Brain_IsoSeq) — it is introduced specifically by sharding and solved by the prefix step.

### New processes

| Process | File | Replaces |
|---------|------|---------|
| `SPLIT_BAM_BY_SHARD` | `modules/align.nf` | — (new) |
| `ISOQUANT_DISCOVERY_SHARD` | `modules/quantify.nf` | `ISOQUANT_DISCOVERY` |
| `MERGE_SHARD_GTFS` | `modules/quantify.nf` | — (new) |
| `SQANTI3_SPLIT_GTF` | `modules/quantify.nf` | — (new) |
| `SQANTI3_QC_CHUNK` | `modules/quantify.nf` | `SQANTI3_QC` |
| `SQANTI3_MERGE_CHUNKS` | `modules/quantify.nf` | — (new) |

### Before starting Phase 1

Run Phase 0 on the cluster:
```bash
# Confirm novel ID format
find /sc-scratch/sc-scratch-cc12-ag-romagnani/nf_work_longshot \
    -path "*/ISOQUANT_DISCOVERY/*" -name "*.gtf" | head -1 \
    | xargs grep -m 20 'novel_'

# Confirm --process_only_chr syntax
apptainer exec <isoquant_sif> isoquant --help 2>&1 | grep -A3 'process_only_chr'
```

### Chromosome shard groups (by gene density, not raw Mb)

chr19 is the most gene-dense chromosome (~1,471 genes in 59 Mb) — must be a singleton. Do NOT group by raw Mb; chr18 + chr19 would be catastrophically unbalanced.

```
Shard 01: chr1                  ~1998 genes
Shard 02: chr19                 ~1471 genes
Shard 03: chr11                 ~1301 genes
Shard 04: chr2                  ~1255 genes
Shard 05: chr17                 ~1197 genes
Shard 06: chr5, chr13           ~1202 genes
Shard 07: chr6, chr21           ~1282 genes
Shard 08: chr10, chr22          ~1220 genes
Shard 09: chr4, chr20           ~1293 genes
Shard 10: chr14, chr8           ~1288 genes
Shard 11: chr3                  ~1073 genes
Shard 12: chr16, chr18          ~1112 genes
Shard 13: chr12, chrM           ~1047 genes
Shard 14: chrX, chrY             ~904 genes
Shard 15: chr7                   ~903 genes
Shard 16: chr9                   ~781 genes
Shard 17: chr15                  ~599 genes
```

### Resource targets (starting point — adjust after first run)

```
SPLIT_BAM_BY_SHARD:       4 CPUs /   8 GB /  2 h
ISOQUANT_DISCOVERY_SHARD: 16 CPUs / 128 GB / 12 h  (monitor chr1; may need 256 GB)
MERGE_SHARD_GTFS:          2 CPUs /  16 GB /  1 h
SQANTI3_QC_CHUNK:          8 CPUs /  64 GB /  6 h
SQANTI3_MERGE_CHUNKS:      2 CPUs /  16 GB / 30 m
```

---

## Findings from Reference Pipeline Review

Incorporated into the plans above. Summary for context:

### RajLabMSSM isoseq pipeline
- **`--dump-clips --dump-removed` on lima** → added to Plan 02
- **SQANTI3 split into 8 chunks** → Part B of Plan 03
- Confirmed: they run IsoQuant on all samples jointly (no sharding) — novel IDs are consistent because one job assigns them
- Uses `isoseq3` (v3.2.2, 2019-era) — obsolete for our use case; we use `isoseq` v4.3.0

### Dev_Brain_IsoSeq pipeline
- **`--fl_data --check_canonical --count_exons` on IsoQuant** → added to Plan 02 and Plan 03
- Uses TALON (database approach) not IsoQuant — TALON enforces ID consistency via SQLite DB across all samples; not applicable to our stack but confirms the "merge all samples" requirement
- **Multi-donor criterion for novel transcripts:** novel transcripts kept only if seen in ≥2 donors, OR have external CAGE/PolyA support — consider adding to exporter post-processing (not in scope for Plans 02–03)
- **Fraction_As ≤ 0.75 internal priming filter** — partially covered by `--require-polya` in refine; evaluate separately
- **Strand filter:** remove transcripts with `strand == "*"` before export — add to exporter (not in scope here)

### Novel ID architecture (clarification)
Novel ID collision is **not** a problem in unsharded pipelines. Both reference pipelines run a single joint discovery job → one ID assignment pass → no duplicates. The collision exists only in our sharded design and is solved by the chromosome-prefix step in `MERGE_SHARD_GTFS`.

---

## Execution Order for New Session

```
1. Plan 02, Phase 0   — cluster: find sidecar filenames
2. Plan 02, Phase 1   — add --dump-clips --dump-removed to lima
3. Plan 02, Phase 2   — capture refine/correct/skera stats
4. Plan 02, Phase 3   — publishDir for flagstat/lima/instrument stats
5. Plan 02, Phase 4   — add --fl_data --check_canonical --count_exons to IsoQuant
6. Dry-run verify     — nextflow run -preview
7. Plan 03, Phase 0   — cluster: confirm novel ID format + --process_only_chr syntax
8. Plan 03, Phase 1   — shard channel definition
9. Plan 03, Phase 2   — SPLIT_BAM_BY_SHARD process
10. Plan 03, Phase 3  — ISOQUANT_DISCOVERY_SHARD process
11. Plan 03, Phase 4  — MERGE_SHARD_GTFS process (novel ID prefixing)
12. Plan 03, Phase 5–7 — SQANTI3 chunking
13. Plan 03, Phase 8  — nextflow.config resource overrides
14. Plan 03, Phase 9  — validation run on 2-library test case
```
