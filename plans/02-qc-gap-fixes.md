# Plan 02: QC Gap Fixes

**Goal:** Fix all QC publishing gaps so every attrition step is visible and all intermediate stats land in `results/` for direct use in R. Add lima flags to capture more read-level QC. Harden IsoQuant flags.

---

## Background

Revio exports **only** post-CCS, post-sample-index-lima files (`hifi_reads/*.bcM0001.bam`). No pre-lima combined BAM exists. Current pipeline entry is the true upstream limit.

**QC blind spots identified in audit:**

| Gap | Impact |
|-----|--------|
| `isoseq refine` filter summary not captured | Can't see poly-A / chimera attrition |
| `isoseq correct` stats not captured | Can't see barcode match rate per library |
| Skera split stats not captured | Can't see segmentation yield per SMRT cell |
| `SAMTOOLS_FLAGSTAT` not published to `results/` | Only in MultiQC; can't parse in R |
| Lima reports not published individually | Only in MultiQC |
| Instrument stats (CCS reports) not published individually | Only in MultiQC |
| Lima missing `--dump-clips --dump-removed` | No visibility into clipped/rejected reads |

**Sources:** RajLabMSSM pipeline (`single_cell_pipeline.smk:77–86`) uses `--dump-clips --dump-removed`; Dev_Brain_IsoSeq uses `--fl_data --check_canonical` on IsoQuant.

---

## Phase 0: Stats File Discovery

Run on the cluster before coding. Confirm exact filenames emitted by each tool.

```bash
# isoseq refine sidecar files
find /sc-scratch/sc-scratch-cc12-ag-romagnani/nf_work_longshot \
    -path "*/ISOSEQ_REFINE/*" \
    \( -name "*.filter_summary*" -o -name "*.report.csv" \) | head -20

# isoseq correct stats
find /sc-scratch/sc-scratch-cc12-ag-romagnani/nf_work_longshot \
    -path "*/ISOSEQ_CORRECT/*" \
    \( -name "*.json" -o -name "*.csv" -o -name "*.log" \) | head -20

# skera split stats
find /sc-scratch/sc-scratch-cc12-ag-romagnani/nf_work_longshot \
    -path "*/SKERA_SPLIT/*" -name "*.log" -o \
    -path "*/SKERA_SPLIT/*" -name "*report*" | head -20
```

Document exact suffixes before Phase 1.

---

## Phase 1: Lima Flag Additions

**File:** `modules/preprocess.nf` — `LIMA_ISOSEQ` and `LIMA_MULTIPLEX`

Add to both lima processes:
```bash
lima --isoseq \
    --dump-clips \
    --dump-removed \
    -j ${task.cpus} \
    ...
```

Add clipped/removed BAM outputs to output block:
```nextflow
path "${meta.id}_fl.lima.clips",          optional: true, emit: lima_clips
path "${meta.id}_fl.lima.removed.bam",    optional: true, emit: lima_removed
```

These files show exactly which reads were trimmed vs rejected — useful for diagnosing low yield.

---

## Phase 2: Capture Skera and Isoseq Stats

**File:** `modules/preprocess.nf`

### 2a. SKERA_SPLIT — add report output

Add to output block (adjust suffix per Phase 0):
```nextflow
path "${meta.id}.skera.log", optional: true, emit: skera_report
```

`publishDir`:
```nextflow
publishDir { "${params.outdir}/${meta.experiment}/qc/skera" }, mode: 'copy'
```

### 2b. ISOSEQ_REFINE — add filter summary

```nextflow
path "${meta.id}_fltnc.filter_summary.json", optional: true, emit: filter_summary
path "${meta.id}_fltnc.report.csv",          optional: true, emit: report
```

`publishDir`:
```nextflow
publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/qc/isoseq_refine" }, mode: 'copy'
```

### 2c. ISOSEQ_CORRECT — add stats

If `isoseq correct` emits a stats file (confirm Phase 0), add output + publishDir to `qc/isoseq_correct/`. If stderr only, capture as `.log`:
```bash
isoseq correct ... 2> ${meta.sample_id}_correct.log
```

Emit from subworkflow. Mix all new stat files into `ch_multiqc_reports` in `main.nf`.

---

## Phase 3: Publish Flagstat, Lima Reports, Instrument Stats

**Files:** `modules/qc.nf`, `modules/preprocess.nf`

### 3a. SAMTOOLS_FLAGSTAT
```nextflow
publishDir { "${params.outdir}/${meta.experiment ?: meta.id}/qc/flagstat" }, mode: 'copy'
```
Three stages (fltnc, dedup, aligned) publish to same dir, distinguished by `${stage}` in filename.

### 3b. LIMA_ISOSEQ and LIMA_MULTIPLEX
```nextflow
publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/qc/lima" }, mode: 'copy'
```

### 3c. COLLECT_INSTRUMENT_STATS
```nextflow
publishDir { "${params.outdir}/${meta.id}/qc/instrument" }, mode: 'copy'
```

---

## Phase 4: IsoQuant Flag Hardening

**File:** `modules/quantify.nf` — `ISOQUANT_DISCOVERY` (and `ISOQUANT_DISCOVERY_SHARD` from Plan 03 when implemented)

Add flags identified from Dev_Brain_IsoSeq (`run_isoQuant.bash:28–41`):

```bash
isoquant \
    --reference ${params.ref_fasta} \
    --genedb ${params.ref_gtf} \
    --complete_genedb \
    --process_only_chr ${chrs} \
    --bam ${bam_args} \
    --data_type pacbio \
    --fl_data \           # NEW: mark reads as full-length (correct for post-dedup reads)
    --check_canonical \   # NEW: validate canonical splice sites
    --count_exons \       # NEW: add exon-level counts to output
    --barcoded_bam \
    --barcode_tag CB \
    --umi_tag XM \
    --read_group barcode \
    --threads ${task.cpus} \
    --prefix ${experiment} \
    --output isoquant_out
```

**Why `--fl_data`:** Tells IsoQuant these are full-length reads (post-`isoseq groupdedup`), not noisy short reads. Changes internal scoring — this is the correct flag for our data type.

**Why `--check_canonical`:** Validates GT-AG / AT-AC splice sites. Filters artefactual junctions from MAS-seq chimeras.

**Why `--count_exons`:** Adds exon-level quantification table — useful for downstream exon-skipping analysis.

---

## Phase 5: Verification

```bash
# Dry run
nextflow run main.nf -preview -profile standard --samplesheet assets/example_samplesheet.csv

# Post-cluster-run checks
find results/ -name "*.flagstat" | wc -l         # expect 3 × libraries
find results/ -name "*.filter_summary.json" | head
find results/ -name "*.lima.*" -path "*/qc/lima/*" | head
find results/ -path "*/qc/instrument/*" | head

# R smoke test
Rscript -e "
  jsonlite::fromJSON('results/.../qc/isoseq_refine/sample_fltnc.filter_summary.json')
  readr::read_delim('results/.../qc/flagstat/sample_fltnc.flagstat', delim=' ', col_names=FALSE)
"
```

---

## Files Changed Summary

| File | Change |
|------|--------|
| `modules/preprocess.nf` | `--dump-clips --dump-removed` on lima; add outputs to `SKERA_SPLIT`, `ISOSEQ_REFINE`, `ISOSEQ_CORRECT`; add `publishDir` to lima processes |
| `subworkflows/preprocess.nf` | Emit new stat channels |
| `modules/qc.nf` | Add `publishDir` to `SAMTOOLS_FLAGSTAT`, `COLLECT_INSTRUMENT_STATS` |
| `modules/quantify.nf` | Add `--fl_data --check_canonical --count_exons` to IsoQuant |
| `main.nf` | Mix new stat channels into `ch_multiqc_reports` |
