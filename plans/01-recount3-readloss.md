# Plan: recount3 junction database + per-step read-loss tracking

## Context

Pipeline: `longshot` — PacBio MAS-seq single-cell long-read RNA-seq.

Two goals:
1. Replace intropolis with recount3 as the junction coverage source for SQANTI3.
2. Add read-count tracking at every processing step, with per-cell resolution where feasible.

### Current state

- `params.intropolis` → passed to `SQANTI3_QC` as `--coverage` flag (`modules/quantify.nf:59`).
- `SAMTOOLS_FLAGSTAT` already called: after `ISOSEQ_REFINE` (preprocess) and after `PBMM2_ALIGN` (align). Both fed into MultiQC via `main.nf:159–160`.
- No flagstat at: SKERA_SPLIT, LIMA, ISOSEQ_CORRECT, ISOSEQ_GROUPDEDUP, ISOQUANT output.
- Per-cell tracking: not implemented anywhere.

### SQANTI3 `--coverage` format

From intropolis: tab-separated, columns `chrom start end dot dot strand count1 count2 ...` (one column per sample). SQANTI3 accepts a single-sample file; only the count column matters.

---

## Phase 0 — Documentation discovery (run before touching code)

Deploy a subagent to:

1. **Recount3 junction format**: Read `https://rna.recount.bio/docs/` or run
   ```bash
   Rscript -e "library(recount3); ?read_counts"
   ```
   to confirm the exact column layout of recount3 junction files (`.RDS` vs flat TSV).

2. **SQANTI3 `--coverage` spec**: Read SQANTI3 docs or source (`sqanti3_qc.py --help`) to confirm:
   - Exact column names/order expected
   - Whether multi-sample columns are summed or max'd
   - Whether 0-based or 1-based coordinates

3. **recount3 download endpoints**: Identify the flat-file URL pattern for junction data
   (`recount3.recount.bio/data/{organism}/{project}/...`). Confirm TSV format available without R.

4. **IsoQuant `read_to_transcript.tsv`**: Read IsoQuant docs to confirm columns — especially
   whether CB tag is preserved, what the read-id column looks like, and whether unmapped
   reads appear with a null transcript ID.

Deliverable: fill in the "Allowed APIs / file formats" table below before Phase 1.

```
| Item                          | Source file/URL | Confirmed value |
|-------------------------------|-----------------|-----------------|
| SQANTI3 --coverage col layout | ...             | ...             |
| recount3 junction TSV cols    | ...             | ...             |
| recount3 download URL pattern | ...             | ...             |
| IsoQuant read_to_tx.tsv cols  | ...             | ...             |
```

---

## Phase 1 — Recount3 junction reference preparation

### Goal

Produce a pre-computed junction file in intropolis format from recount3 data for hg38,
stored as a reference file on the cluster (not computed per-run).

### Task

Create `assets/prepare_recount3_junctions.R` (run once offline, not in pipeline):

```r
library(recount3)

# 1. Pull GTEx or SRA project junction data for hg38
# ponytail: GTEx is the closest tissue-matched equivalent to intropolis GTEx data
proj <- available_projects(organism = "human")
rse  <- create_rse_manual(
  project      = "GTEX",          # adjust to chosen project
  project_home = "data_sources/gtex",
  organism     = "human",
  annotation   = "gencode_v29",
  type         = "jxn"
)

# 2. Extract junction counts and reshape to intropolis-compatible TSV
# SQANTI3 --coverage: chrom\tstart\tend\t.\t.\tstrand\tcount
jxn <- rowRanges(rse)
counts <- rowSums(assay(rse, "counts"))  # aggregate across all samples

out <- data.frame(
  chrom  = as.character(seqnames(jxn)),
  start  = start(jxn) - 1,  # convert to 0-based if SQANTI3 expects it (confirm in Phase 0)
  end    = end(jxn),
  dot1   = ".",
  dot2   = ".",
  strand = as.character(strand(jxn)),
  count  = counts
)
out <- out[out$count >= 10, ]  # match intropolis min_count_10 filter

write.table(out, "recount3.gtex.hg38.min10.junctions.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

**Verification:**
- File exists and has `>100k` lines
- Column count = 7 per line: `awk '{print NF}' | sort -u` returns `7`
- No header row
- Chroms are `chr`-prefixed (required for hg38 GRCh38)

### Pipeline param change

In `nextflow.config`, add:
```groovy
params.recount3_junctions = "/sc-projects/.../ref/hs/pacbio/recount3.gtex.hg38.min10.junctions.tsv"
params.intropolis          = null   // deprecated; keep for backward compat
```

In `modules/quantify.nf:59`, change:
```groovy
// OLD:
def intropolis_arg = params.intropolis ? "--coverage ${params.intropolis}" : ""

// NEW:
def junction_db = params.recount3_junctions ?: params.intropolis
def intropolis_arg = junction_db ? "--coverage ${junction_db}" : ""
```

**Anti-pattern guard:** Do NOT rename the Nextflow variable — keeping `intropolis_arg` avoids
touching the shell command block and keeps git diff minimal.

---

## Phase 2 — Per-step read-loss tracking (pipeline-wide flagstat)

### Goal

Every major BAM-producing step emits a flagstat file. All flagstats feed into MultiQC.
MultiQC will render a "read counts per step" table automatically from flagstat files
that have distinct stage names in their filename.

### Current flagstat call sites

| Step | Process | Flagstat? |
|------|---------|-----------|
| After ISOSEQ_REFINE | `preprocess.nf:101` | ✅ yes (stage = "refine") |
| After PBMM2_ALIGN   | `align.nf:70`       | ✅ yes (stage = "align")  |
| After SKERA_SPLIT   | —                   | ❌ missing |
| After LIMA          | —                   | ❌ missing |
| After ISOSEQ_CORRECT| —                   | ❌ missing |
| After ISOSEQ_GROUPDEDUP | —               | ❌ missing |

### Implementation

`SAMTOOLS_FLAGSTAT` in `modules/qc.nf` already accepts `[meta, bam, stage]`. Add calls at
each missing step.

#### `subworkflows/preprocess.nf`

After `SKERA_SPLIT` (around line 80):
```groovy
ch_skera_bam
    .map { meta, bam -> [meta, bam, "skera"] }
    | SAMTOOLS_FLAGSTAT
```

After `LIMA_ISOSEQ` / `LIMA_MULTIPLEX` merge (after the branch/join, before ISOSEQ_TAG):
```groovy
ch_lima_bam
    .map { meta, bam -> [meta, bam, "lima"] }
    | SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_LIMA
```

Collect all new flagstat outputs into `preprocess.nf` emit:
```groovy
emit:
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
                   .mix(SAMTOOLS_FLAGSTAT_LIMA.out.flagstat)
                   .mix(SAMTOOLS_FLAGSTAT_SKERA.out.flagstat)
```

#### `subworkflows/align.nf`

After `ISOSEQ_CORRECT` (insert before SAMTOOLS_SORT_CB):
```groovy
ISOSEQ_CORRECT.out.bam
    .map { meta, bam -> [meta, bam, "correct"] }
    | SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_CORRECT
```

After `ISOSEQ_GROUPDEDUP` (insert before PBMM2_ALIGN):
```groovy
ISOSEQ_GROUPDEDUP.out.bam
    .map { meta, bam -> [meta, bam, "dedup"] }
    | SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_DEDUP
```

Add to `align.nf` emit:
```groovy
flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
               .mix(SAMTOOLS_FLAGSTAT_CORRECT.out.flagstat)
               .mix(SAMTOOLS_FLAGSTAT_DEDUP.out.flagstat)
```

**Verification:**
- `nextflow run . -stub` completes without channel errors
- MultiQC output contains a "General Statistics" table with columns: skera, lima, refine, correct, dedup, align
- Read counts decrease monotonically down the table (sanity check)

---

## Phase 3 — Per-cell read-loss tracking

### Goal

For each library, produce a TSV: `{library_id}.per_cell_readloss.tsv` with columns:
```
CB  reads_post_correct  reads_post_dedup  reads_post_align  reads_quantified
```

This is the most expensive part. Feasibility by step:

| Step | Method | Cost |
|------|--------|------|
| post-ISOSEQ_CORRECT | `samtools view -F 4 \| python count_by_tag.py CB` | fast |
| post-ISOSEQ_GROUPDEDUP | same | fast |
| post-PBMM2_ALIGN | `samtools view -F 4 \| python count_by_tag.py CB` | fast |
| post-ISOQUANT_QUANTIFY | parse `read_to_transcript.tsv`, group by CB tag in read-id | moderate |

Pre-ISOSEQ_CORRECT steps (SKERA, LIMA, REFINE) do not have reliable CB tags — skip
per-cell tracking there. Report total-read loss only (Phase 2 flagstat covers this).

### New process: `COUNT_READS_PER_CELL`

Add to `modules/qc.nf`:

```groovy
process COUNT_READS_PER_CELL {
    tag "${meta.library_id}:${stage}"
    container params.samtools_container

    input:
    tuple val(meta), path(bam), val(stage)

    output:
    tuple val(meta), path("${meta.library_id}_${stage}.per_cell.tsv"), emit: per_cell

    script:
    """
    samtools view -F 4 ${bam} \\
        | python3 -c "
import sys, collections
counts = collections.Counter()
for line in sys.stdin:
    fields = line.split('\\t')
    tags = {f[:2]: f[5:] for f in fields[11:] if len(f) > 5 and f[2] == ':'}
    cb = tags.get('CB', tags.get('CR', 'NA'))
    counts[cb] += 1
print('CB\\treads')
for cb, n in sorted(counts.items()):
    print(f'{cb}\\t{n}')
" > ${meta.library_id}_${stage}.per_cell.tsv
    """
}
```

Call this in `subworkflows/align.nf` after ISOSEQ_CORRECT, ISOSEQ_GROUPDEDUP, and PBMM2_ALIGN
(same insertion points as Phase 2).

### IsoQuant per-cell counts

Add a new process `COUNT_ISOQUANT_PER_CELL` in `modules/exporter.nf`:

```groovy
process COUNT_ISOQUANT_PER_CELL {
    tag "${meta.library_id}"
    container params.python_container

    input:
    tuple val(meta), path(isoquant_dir)

    output:
    tuple val(meta), path("${meta.library_id}_quantified.per_cell.tsv"), emit: per_cell

    script:
    """
    python3 -c "
import sys, collections, gzip, os, glob

# read_to_transcript.tsv: read_id \\t transcript_id (or '.' if unassigned)
# read_id encodes CB as last field after '#' separator (IsoQuant convention — verify in Phase 0)
f = glob.glob('${isoquant_dir}/**/read_to_transcript.tsv', recursive=True)[0]
counts = collections.Counter()
with open(f) as fh:
    for line in fh:
        read_id, tx = line.rstrip().split('\\t')[:2]
        if tx == '.':
            continue
        # CB is encoded in read_id by IsoQuant as: readname#CB_UMI
        cb = read_id.split('#')[1].split('_')[0] if '#' in read_id else 'NA'
        counts[cb] += 1
print('CB\\treads')
for cb, n in sorted(counts.items()):
    print(f'{cb}\\t{n}')
" > ${meta.library_id}_quantified.per_cell.tsv
    """
}
```

**Important:** The read_id format in IsoQuant output must be confirmed in Phase 0 before
implementing the CB extraction logic above.

### Joining per-cell tables

Add a final process `MERGE_PER_CELL_READLOSS` in `modules/exporter.nf` that joins the four
per-cell TSVs by CB for each library:

```groovy
process MERGE_PER_CELL_READLOSS {
    tag "${meta.library_id}"
    container params.python_container

    input:
    tuple val(meta), path(correct_tsv), path(dedup_tsv), path(align_tsv), path(quant_tsv)

    output:
    tuple val(meta), path("${meta.library_id}.per_cell_readloss.tsv"), emit: readloss

    script:
    """
    python3 -c "
import pandas as pd
c = pd.read_csv('${correct_tsv}', sep='\\t').rename(columns={'reads': 'reads_post_correct'})
d = pd.read_csv('${dedup_tsv}',   sep='\\t').rename(columns={'reads': 'reads_post_dedup'})
a = pd.read_csv('${align_tsv}',   sep='\\t').rename(columns={'reads': 'reads_post_align'})
q = pd.read_csv('${quant_tsv}',   sep='\\t').rename(columns={'reads': 'reads_quantified'})
out = c.merge(d, on='CB', how='outer') \\
        .merge(a, on='CB', how='outer') \\
        .merge(q, on='CB', how='outer') \\
        .fillna(0).astype({'reads_post_correct': int, 'reads_post_dedup': int,
                           'reads_post_align': int, 'reads_quantified': int})
out.to_csv('${meta.library_id}.per_cell_readloss.tsv', sep='\\t', index=False)
"
    """
}
```

Wire inputs in `subworkflows/export.nf` or `main.nf` by joining the four per_cell channels
on `meta.library_id`.

**Verification:**
- Output file exists for each library
- `CB` column contains valid 10x barcodes (16-char strings)
- `reads_post_correct >= reads_post_dedup >= reads_post_align >= reads_quantified` for all rows

---

## Phase 4 — Verification

1. `nextflow run . -stub -profile slurm,agilent_v7` — full DAG renders, no channel errors.
2. On a single small library (≤1M reads), run with `-profile slurm,agilent_v7,umi`:
   - MultiQC report shows 6 flagstat stages in General Statistics table
   - `results/{library_id}.per_cell_readloss.tsv` present and non-empty
   - SQANTI3 classification files identical structure to intropolis run (junction column populated)
3. Grep check — no old hardcoded intropolis path in any `.nf` file:
   ```bash
   grep -r "intropolis.v1" modules/ subworkflows/ main.nf
   # should return nothing
   ```

---

## Files touched

| File | Change |
|------|--------|
| `nextflow.config` | Add `params.recount3_junctions`, nullify `params.intropolis` |
| `modules/quantify.nf:59` | Use `recount3_junctions ?: intropolis` fallback |
| `modules/qc.nf` | Add `COUNT_READS_PER_CELL` process |
| `modules/exporter.nf` | Add `COUNT_ISOQUANT_PER_CELL`, `MERGE_PER_CELL_READLOSS` |
| `subworkflows/preprocess.nf` | Add flagstat after SKERA and LIMA |
| `subworkflows/align.nf` | Add flagstat + per-cell count after CORRECT and DEDUP |
| `subworkflows/export.nf` | Wire per-cell join and merge |
| `assets/prepare_recount3_junctions.R` | One-off reference prep script (not in DAG) |

---

## Known risks / open questions

- **IsoQuant read_id CB encoding**: CB position in read_id is assumed `#CB_UMI` — must confirm
  from IsoQuant docs or test output before Phase 3. If wrong, the per-cell quant count is broken.
- **Coordinate system**: recount3 junctions may be 1-based; SQANTI3 may expect 0-based (or vice versa).
  Confirm in Phase 0 or SQANTI3 will silently fail to match junctions.
- **LIMA output channels**: LIMA has two code paths (`LIMA_ISOSEQ` and `LIMA_MULTIPLEX`). The
  flagstat after LIMA must be added to both branches or the channel mix will be unbalanced.
- **Container for COUNT_READS_PER_CELL**: Uses `params.samtools_container` + inline Python.
  If the samtools container lacks Python 3, either inline awk instead or add a python container param.
