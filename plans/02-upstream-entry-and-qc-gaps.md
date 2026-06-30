# Plan 02: QC Gap Fixes

**Goal:** Fix all QC publishing gaps so every attrition step is visible and all intermediate stats land in `results/` for direct use in R.

---

## Background

Revio exports **only** the post-CCS, post-sample-index-lima files (`hifi_reads/*.bcM0001.bam`). There is no pre-lima combined HiFi BAM and no subreads BAM available. The `hifi_reads/*.bcM0001.bam` files are the true upstream limit — current pipeline entry is already correct.

The pipeline has the following QC blind spots identified in audit:

| Gap | Impact |
|-----|--------|
| `isoseq refine` filter summary not captured | Can't see poly-A / chimera attrition |
| `isoseq correct` stats not captured | Can't see barcode match rate per library |
| Skera split stats not captured | Can't see segmentation yield per SMRT cell |
| `SAMTOOLS_FLAGSTAT` not published to `results/` | Only in MultiQC; can't parse in R |
| Lima reports not published individually | Only in MultiQC |
| Instrument stats (CCS reports) not published individually | Only in MultiQC |

---

## Phase 0: Stats File Discovery

Before coding, confirm exact filenames emitted by each tool on the cluster.

Run these on the cluster against any existing work directory or via `--help`:

```bash
# skera split output files
apptainer exec <skera_sif> skera split --help 2>&1 | grep -i report

# isoseq refine output files
apptainer exec <isoseq_sif> isoseq refine --help 2>&1 | grep -i report
# Also check work dir for an existing refine run:
find /sc-scratch/sc-scratch-cc12-ag-romagnani/nf_work_longshot -name "*.filter_summary*" -o -name "*.report.csv" 2>/dev/null | head -20

# isoseq correct output files
apptainer exec <isoseq_sif> isoseq correct --help 2>&1 | grep -i report
find /sc-scratch/sc-scratch-cc12-ag-romagnani/nf_work_longshot -path "*/ISOSEQ_CORRECT/*" -name "*.json" -o -path "*/ISOSEQ_CORRECT/*" -name "*.csv" 2>/dev/null | head -20
```

**Expected findings to document before Phase 1:**
- Exact filename suffixes for refine sidecar files (e.g. `*.filter_summary.json`, `*.report.csv`)
- Whether `isoseq correct` writes a stats file at all, or only logs to stderr
- Exact filename suffix for skera split report (e.g. `*.split_counts.tsv`)

---

## Phase 1: Capture Skera and Isoseq Stats

**Files:** `modules/preprocess.nf`

### 1a. SKERA_SPLIT — add report output

Skera emits a split report. Add to output block (adjust filename suffix per Phase 0 findings):

```nextflow
path "${meta.id}_segmented.skera.log", optional: true, emit: skera_report
```

Add `publishDir`:
```nextflow
publishDir { "${params.outdir}/${meta.experiment}/qc/skera" }, mode: 'copy'
```

Emit `skera_report` from `subworkflows/preprocess.nf` and mix into `ch_multiqc_reports` in `main.nf`.

### 1b. ISOSEQ_REFINE — add filter summary output

```nextflow
path "${meta.id}_fltnc.filter_summary.json", optional: true, emit: filter_summary
path "${meta.id}_fltnc.report.csv",          optional: true, emit: report
```

Add `publishDir`:
```nextflow
publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/qc/isoseq_refine" }, mode: 'copy'
```

Emit from subworkflow, mix into MultiQC channel.

### 1c. ISOSEQ_CORRECT — add stats output (if file exists per Phase 0)

If `isoseq correct` writes a stats file:
```nextflow
path "${meta.sample_id}_corrected.report.json", optional: true, emit: correct_report
```

Add `publishDir`:
```nextflow
publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/qc/isoseq_correct" }, mode: 'copy'
```

If it only logs to stderr, capture stderr in script block:
```bash
isoseq correct ... 2> ${meta.sample_id}_correct.log
```
and publish the `.log` file instead.

**Verification:**
```bash
find results/ -path "*/qc/skera/*" | head
find results/ -path "*/qc/isoseq_refine/*" | head
find results/ -path "*/qc/isoseq_correct/*" | head
```

---

## Phase 2: Publish Flagstat, Lima Reports, Instrument Stats

**Files:** `modules/qc.nf`, `modules/preprocess.nf`

### 2a. SAMTOOLS_FLAGSTAT — add publishDir

`modules/qc.nf`. The process uses `meta.sample_id` — path it under experiment using `meta.experiment` if available, else fall back to sample_id:

```nextflow
publishDir { "${params.outdir}/${meta.experiment ?: meta.id}/qc/flagstat" }, mode: 'copy'
```

Note: flagstat is called with three stages (fltnc, dedup, aligned) — all three will publish to the same dir, distinguished by `${stage}` in the filename. That's correct.

### 2b. LIMA_ISOSEQ and LIMA_MULTIPLEX — add publishDir

`modules/preprocess.nf`:

```nextflow
publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/qc/lima" }, mode: 'copy'
```

### 2c. COLLECT_INSTRUMENT_STATS — add publishDir

`modules/qc.nf`:

```nextflow
publishDir { "${params.outdir}/${meta.id}/qc/instrument" }, mode: 'copy'
```

**Verification:**
```bash
find results/ -name "*.flagstat" | wc -l   # expect 3 × number of libraries
find results/ -name "*.lima.*" -path "*/qc/lima/*" | head
find results/ -path "*/qc/instrument/*" | head
```

---

## Phase 3: Verification

1. Dry run passes:
   ```bash
   nextflow run main.nf -preview -profile standard --samplesheet assets/example_samplesheet.csv
   ```

2. On cluster — run one SMRT cell, then confirm:
   ```bash
   # All flagstat stages present
   find results/ -name "*.flagstat" | sort

   # Refine attrition visible
   find results/ -name "*.filter_summary.json" | xargs -I{} python3 -c "import json,sys; d=json.load(open('{}'));print('{}:',d)"

   # Instrument CCS stats present
   find results/ -name "*.ccs_report.json" -path "*/qc/instrument/*"
   ```

3. R smoke test — confirm files parse cleanly:
   ```r
   library(jsonlite)
   library(readr)
   fromJSON("results/.../qc/isoseq_refine/sample_fltnc.filter_summary.json")
   read_delim("results/.../qc/flagstat/sample_fltnc.flagstat", delim=" ", col_names=FALSE)
   ```

---

## Files Changed Summary

| File | Change |
|------|--------|
| `modules/preprocess.nf` | Add outputs + publishDir to `SKERA_SPLIT`, `ISOSEQ_REFINE`, `ISOSEQ_CORRECT`, `LIMA_ISOSEQ`, `LIMA_MULTIPLEX` |
| `subworkflows/preprocess.nf` | Emit new stat channels from each process |
| `modules/qc.nf` | Add publishDir to `SAMTOOLS_FLAGSTAT`, `COLLECT_INSTRUMENT_STATS` |
| `main.nf` | Mix new stat channels into `ch_multiqc_reports` |
