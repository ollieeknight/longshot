# Plan 04 — Replace NanoStat + MosDepth with cramino across all pipeline stages

## Goal

Remove NanoStat and MosDepth entirely. Add a single reusable `CRAMINO` process that runs
at every BAM/CRAM checkpoint in the pipeline, providing a read-loss funnel from raw HiFi
input through to final aligned CRAM. cramino `--spliced --karyotype` on the aligned BAM
supersedes both tools: NanoStat (read-length/quality distributions) and MosDepth (coverage
depth / per-chromosome stats).

---

## Phase 0 — Facts & API Inventory

### cramino CLI (confirmed from cramino/README.md)

```
cramino [OPTIONS] <INPUT>

--threads <N>          parallel decompression threads [default: 4]
--ubam                 metrics for unaligned reads (skip alignment requirement)
--spliced              metrics for spliced/aligned reads (adds junctions etc.)
--karyotype            per-chromosome read distribution
--format text|json|tsv [default: text]
--hist [<FILE>]        read-length histogram
--reference <REF>      required only for CRAM decompression
-V, --version
```

No native MultiQC module confirmed — outputs are plain text files.
NanoStat IS natively parsed by MultiQC; cramino is NOT (as of MultiQC 1.34).
**Gap:** MultiQC will no longer render long-read QC in the HTML report.
**Mitigation:** cramino files still published to disk; MultiQC collects them as
unrecognised custom content. Acceptable for now — flag for follow-up if needed.

### BAM stages available (confirmed from subworkflows + modules)

| # | Stage label | BAM source | cramino flags | Subworkflow |
|---|---|---|---|---|
| 1 | `raw` | samplesheet input `hifi_bam` | `--ubam` | preprocess |
| 2 | `segmented` | `SKERA_SPLIT.out.segmented_bam` | `--ubam` | preprocess |
| 3 | `fl` | `ch_fl_bam` (LIMA_ISOSEQ/LIMA_MULTIPLEX output mix) | `--ubam` | preprocess |
| 4 | `fltnc` | `ISOSEQ_REFINE.out.fltnc_bam` | `--ubam` | preprocess |
| 5 | `merged` | `SAMTOOLS_MERGE_FLTNC.out.merged_bam` | `--ubam` | align |
| 6 | `dedup` | `ISOSEQ_GROUPDEDUP.out.dedup_bam` | `--ubam` | align |
| 7 | `aligned` | `PBMM2_ALIGN.out.aligned_bam` (tuple: meta, bam, bai) | `--spliced --karyotype` | align |

### Existing NanoStat wiring to remove (confirmed)

- `modules/qc.nf:88-114` — NANOSTAT process definition
- `nextflow.config` param `container_nanostat`
- `subworkflows/align.nf:10` — `include { NANOSTAT }` import
- `subworkflows/align.nf:74` — `NANOSTAT(PBMM2_ALIGN.out.aligned_bam)` call
- `subworkflows/align.nf:84` — `nanostat = NANOSTAT.out.stats` emit
- `subworkflows/align.nf:94` — `.mix(NANOSTAT.out.versions)` in versions
- `main.nf:164` — `.mix( ALIGN.out.nanostat )` in MultiQC channel

### Existing MosDepth wiring to remove (confirmed)

- `modules/qc.nf:58-85` — MOSDEPTH process definition
- `nextflow.config` param `container_mosdepth`
- `subworkflows/align.nf:9` — `include { MOSDEPTH }` import
- `subworkflows/align.nf:73` — `MOSDEPTH(PBMM2_ALIGN.out.aligned_bam)` call
- `subworkflows/align.nf:83` — `mosdepth = MOSDEPTH.out.summary` emit
- `subworkflows/align.nf:93` — `.mix(MOSDEPTH.out.versions)` in versions
- `main.nf:163` — `.mix( ALIGN.out.mosdepth )` in MultiQC channel

### Container needed

Add to `nextflow.config` params block:
```
container_cramino = "quay.io/biocontainers/cramino:0.14.5--hdbdd923_0"
```
Verify latest tag at: https://quay.io/repository/biocontainers/cramino?tab=tags
(Use `conda search -c bioconda cramino` on HPC if quay.io tag is uncertain.)

---

## Phase 1 — Add CRAMINO process; remove NANOSTAT

**File:** `modules/qc.nf`

### 1a. Add CRAMINO process (after MOSDEPTH, before MULTIQC)

```nextflow
process CRAMINO {
    tag "${meta.sample_id} (${stage})"
    label 'process_low'
    container "${params.container_cramino}"
    publishDir { "${params.outdir}/${meta.experiment ?: meta.id}/qc/cramino" }, mode: 'copy'

    input:
    tuple val(meta), val(stage), val(extra_flags), path(bam)

    output:
    path "${meta.sample_id}_${stage}.cramino.txt", emit: stats
    path "versions.yml",                            emit: versions

    script:
    """
    cramino \\
        --threads ${task.cpus} \\
        ${extra_flags} \\
        ${bam} \\
        > ${meta.sample_id}_${stage}.cramino.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cramino: \$(cramino --version 2>&1 | sed 's/cramino //')
    END_VERSIONS
    """
}
```

### 1b. Delete NANOSTAT and MOSDEPTH processes

- Remove lines 88–114 of `modules/qc.nf` (the full NANOSTAT process block).
- Remove lines 58–85 of `modules/qc.nf` (the full MOSDEPTH process block).

### 1c. Add container param; remove nanostat + mosdepth params

In `nextflow.config` params block:
- Add: `container_cramino = "quay.io/biocontainers/cramino:0.14.5--hdbdd923_0"`
- Remove: `container_nanostat = "..."`
- Remove: `container_mosdepth = "..."`

**Verification:**
```bash
grep -n 'NANOSTAT\|nanostat\|MOSDEPTH\|mosdepth' modules/qc.nf nextflow.config
# should return nothing
grep -n 'CRAMINO\|cramino\|container_cramino' modules/qc.nf nextflow.config
# should return the new process and param
```

---

## Phase 2 — Wire cramino in PREPROCESS subworkflow

**File:** `subworkflows/preprocess.nf`

### 2a. Add import

```nextflow
include { SAMTOOLS_FLAGSTAT; CRAMINO } from '../modules/qc'
```

### 2b. Add cramino calls

After the existing `ISOSEQ_REFINE` / `SAMTOOLS_FLAGSTAT` block, add cramino at four
preprocess checkpoints. Pattern: `.map { meta, bam -> [ meta, '<stage>', '<flags>', bam ] } | CRAMINO`

```nextflow
// cramino: raw HiFi BAM
ch_input
    .map { meta, bam -> [ meta, 'raw', '--ubam', bam ] }
    | CRAMINO

// cramino: post-Skera segmented BAM
SKERA_SPLIT.out.segmented_bam
    .map { meta, bam -> [ meta, 'segmented', '--ubam', bam ] }
    | CRAMINO

// cramino: post-Lima FL BAM
ch_fl_bam
    .map { meta, bam -> [ meta, 'fl', '--ubam', bam ] }
    | CRAMINO

// cramino: post-Refine FLTNC BAM
ISOSEQ_REFINE.out.fltnc_bam
    .map { meta, bam -> [ meta, 'fltnc', '--ubam', bam ] }
    | CRAMINO
```

Note: `ch_input` in preprocess.nf is the raw samplesheet channel. Check it is
accessible at the point of insertion (it is — it's the `take:` channel, used
throughout the workflow body).

### 2c. Add cramino emit and versions mix

```nextflow
emit:
// ... existing emits ...
cramino_reports  = CRAMINO.out.stats
versions = Channel.empty()
    // ... existing mixes ...
    .mix(CRAMINO.out.versions)
```

Note: CRAMINO is called 4× here, so `.out.stats` emits 4 items per run.
All four wire into the single `cramino_reports` channel; collect at MultiQC time.

**Verification:**
```bash
grep -n 'CRAMINO\|cramino' subworkflows/preprocess.nf
# should show import, 4 calls, emit, versions mix
```

---

## Phase 3 — Wire cramino in ALIGN subworkflow; remove NANOSTAT

**File:** `subworkflows/align.nf`

### 3a. Update import

```nextflow
# Remove:
include { SAMTOOLS_FLAGSTAT; MOSDEPTH; NANOSTAT } from '../modules/qc'
# Replace with:
include { SAMTOOLS_FLAGSTAT; CRAMINO             } from '../modules/qc'
```

### 3b. Remove MOSDEPTH and NANOSTAT calls

Remove both calls:
```nextflow
MOSDEPTH(PBMM2_ALIGN.out.aligned_bam)
NANOSTAT(PBMM2_ALIGN.out.aligned_bam)
```

Add cramino at three align checkpoints:
```nextflow
// cramino: post-merge FLTNC (unaligned)
SAMTOOLS_MERGE_FLTNC.out.merged_bam
    .map { meta, bam -> [ meta, 'merged', '--ubam', bam ] }
    | CRAMINO

// cramino: post-dedup (unaligned)
ISOSEQ_GROUPDEDUP.out.dedup_bam
    .map { meta, bam -> [ meta, 'dedup', '--ubam', bam ] }
    | CRAMINO

// cramino: post-alignment (spliced, with karyotype)
PBMM2_ALIGN.out.aligned_bam
    .map { meta, bam, bai -> [ meta, 'aligned', '--spliced --karyotype', bam ] }
    | CRAMINO
```

### 3c. Update emit block

```nextflow
emit:
aligned_bam  = PBMM2_ALIGN.out.aligned_bam
cram         = GENERATE_CRAM.out.cram
flagstat     = SAMTOOLS_FLAGSTAT.out.flagstat
cramino      = CRAMINO.out.stats      // replaces mosdepth + nanostat emits
versions     = Channel.empty()
                .mix(SAMTOOLS_MERGE_FLTNC.out.versions)
                .mix(PREPARE_WHITELIST.out.versions)
                .mix(ISOSEQ_CORRECT.out.versions)
                .mix(SAMTOOLS_SORT_CB.out.versions)
                .mix(ISOSEQ_GROUPDEDUP.out.versions)
                .mix(PBMM2_ALIGN.out.versions)
                .mix(SAMTOOLS_FLAGSTAT.out.versions)
                .mix(CRAMINO.out.versions)    // replaces MOSDEPTH + NANOSTAT versions
                .mix(GENERATE_CRAM.out.versions)
```

**Verification:**
```bash
grep -n 'NANOSTAT\|nanostat\|MOSDEPTH\|mosdepth' subworkflows/align.nf
# should return nothing
grep -n 'CRAMINO\|cramino' subworkflows/align.nf
# should show import, 3 calls, emit, versions mix
```

---

## Phase 4 — Update main.nf MultiQC wiring

**File:** `main.nf`

### 4a. Replace mosdepth + nanostat lines in MultiQC channel

```nextflow
# Remove (lines 163-164):
        .mix( ALIGN.out.mosdepth )
        .mix( ALIGN.out.nanostat )
# Replace with:
        .mix( ALIGN.out.cramino.flatten() )
        .mix( PREPROCESS.out.cramino_reports.flatten() )
```

The `.flatten()` is needed because each cramino emit carries multiple files (one per stage
invocation); flatten ensures individual path items enter the MultiQC stageAs pattern correctly.

**Verification:**
```bash
grep -n 'nanostat\|NANOSTAT' main.nf
# should return nothing
grep -n 'cramino\|CRAMINO' main.nf
# should return the two new mix lines
```

---

## Phase 5 — Add post-Skera length filter (SAMTOOLS_LENGTH_FILTER)

**Rationale:** Skera can produce very short artefact segments from incomplete MAS16
splitting. Removing these pre-Lima reduces Lima runtime and false demux assignments.
PacBio HiFi reads are already CCS-filtered to Q20, so quality filtering adds nothing.

**Approach:** `samtools view -e` expression filter — lossless (all BAM tags preserved),
no new container (reuses `container_samtools`), no FASTQ round-trip.

> **Why not chopper?** chopper is FASTQ-only; using it requires `samtools fastq | chopper |
> samtools import` which drops PacBio BAM tags (`rq`, `np`, `zm`) that Lima and downstream
> IsoSeq tools rely on. `samtools view -e` achieves the same length filter without tag loss.

### 5a. Add SAMTOOLS_LENGTH_FILTER process to `modules/preprocess.nf`

```nextflow
process SAMTOOLS_LENGTH_FILTER {
    tag "${meta.id}"
    label 'process_low'
    container "${params.container_samtools}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}_filtered.bam"), emit: filtered_bam
    path "versions.yml",                               emit: versions

    script:
    """
    samtools view \\
        -@ ${task.cpus} \\
        -b \\
        -e "length(seq) >= ${params.min_read_length}" \\
        ${bam} \\
        -o ${meta.id}_filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
```

### 5b. Add param to `nextflow.config`

```groovy
min_read_length = 200   // minimum segment length post-Skera (bp)
```

### 5c. Wire into `subworkflows/preprocess.nf`

Add import:

```nextflow
include { SKERA_SPLIT; ...; SAMTOOLS_LENGTH_FILTER } from '../modules/preprocess'
```

Insert filter step between SKERA_SPLIT and the Lima branching logic:

```nextflow
// Filter Skera artefacts by minimum read length (lossless, tags preserved)
SKERA_SPLIT.out.segmented_bam
    | SAMTOOLS_LENGTH_FILTER

// Replace all downstream uses of SKERA_SPLIT.out.segmented_bam
// with SAMTOOLS_LENGTH_FILTER.out.filtered_bam
```

Concretely, the Lima join currently reads:

```nextflow
SKERA_SPLIT.out.segmented_bam
    .join(ch_grouped_smrtcells.map { smrt_meta, bam, metas -> ... })
```

Change `SKERA_SPLIT.out.segmented_bam` → `SAMTOOLS_LENGTH_FILTER.out.filtered_bam` in
that join and anywhere else `segmented_bam` feeds Lima or downstream processes.

Also update the cramino `segmented` call (Phase 2) to run **before** the filter so the
read-loss is visible:

```nextflow
// cramino: post-Skera (before length filter) — shows raw segment count
SKERA_SPLIT.out.segmented_bam
    .map { meta, bam -> [ meta, 'segmented', '--ubam', bam ] }
    | CRAMINO

// cramino: post-filter — shows count after artefact removal
SAMTOOLS_LENGTH_FILTER.out.filtered_bam
    .map { meta, bam -> [ meta, 'filtered', '--ubam', bam ] }
    | CRAMINO
```

This gives a direct read-loss number at the filter step.

Add to versions mix:

```nextflow
.mix(SAMTOOLS_LENGTH_FILTER.out.versions)
```

**Verification:**
```bash
grep -n 'SAMTOOLS_LENGTH_FILTER\|min_read_length' \
    modules/preprocess.nf subworkflows/preprocess.nf nextflow.config
# should show process def, import+call+versions, param
```

---

## Phase 7 — Per-sample adapter kit customization; rename mas16_primers param

**Goal:** Replace hardcoded `params.mas16_primers` with a flexible `params.adapter_primers`
that defaults to MAS16 but can be overridden globally or per-sample via samplesheet
`adapter_kit` column. (`--emit-all` not available in pbskera 1.4.0 — confirmed.)

### 7a. Update `nextflow.config`

```groovy
// Remove:
mas16_primers    = "${projectDir}/assets/adapters/mas16_primers.fasta"

// Add:
adapter_primers  = "${projectDir}/assets/adapters/mas16_primers.fasta"  // default: MAS16 Kinnex kit
```

### 7b. Add adapter_kit resolution in `main.nf` samplesheet map block

After the existing `chemistry` resolution block, add:

```groovy
// Resolve adapter FASTA: samplesheet adapter_kit column overrides params.adapter_primers
def kit_val = row.containsKey('adapter_kit') && row.adapter_kit ? row.adapter_kit.trim() : ""
def is_null_kit = kit_val == "" || kit_val.equalsIgnoreCase("null") || kit_val.equalsIgnoreCase("na") || kit_val.equalsIgnoreCase("none")
def adapter_file = params.adapter_primers  // default
if (!is_null_kit) {
    if (kit_val == 'mas8')  adapter_file = "${projectDir}/assets/adapters/mas8_primers.fasta"
    else if (kit_val == 'mas12') adapter_file = "${projectDir}/assets/adapters/mas12_primers.fasta"
    else if (kit_val == 'mas16') adapter_file = "${projectDir}/assets/adapters/mas16_primers.fasta"
    else adapter_file = kit_val  // treat as direct file path for custom adapters
}
```

Store as `adapter_file` in the meta map returned from this block.

### 7c. Update `modules/preprocess.nf` SKERA_SPLIT process

Change input to accept adapters as a path (not from params):

```nextflow
input:
tuple val(meta), path(hifi_bam), path(adapters)

script:
"""
skera split \\
    -j ${task.cpus} \\
    ${hifi_bam} \\
    ${adapters} \\
    ${meta.id}_segmented.bam \\
    2> ${meta.id}_skera.log
```

> ℹ️ **`--emit-all` does not exist in pbskera 1.4.0** — confirmed via `skera split --help`.
> Only available flags: `-j/--num-threads`, `--log-level`, `--log-file`.
> Read recovery for adapter-not-found reads is not supported at this version.

### 7d. Update `subworkflows/preprocess.nf` SKERA_SPLIT call

The subworkflow calls SKERA_SPLIT with `ch_grouped_smrtcells`. The adapter file must be
passed per SMRT cell. Since all samples in a SMRT cell share the same adapter kit
(physically the same sequencing run), use the first sample's `adapter_file`:

```nextflow
// Map grouped SMRT cells to include resolved adapter file
ch_grouped_smrtcells
    .map { smrt_meta, bam, metas ->
        def adapter = file(metas[0].adapter_file)
        [ smrt_meta, bam, adapter ]
    }
    | SKERA_SPLIT
```

Replace the existing `SKERA_SPLIT(ch_grouped_smrtcells.map { smrt_meta, bam, metas -> [ smrt_meta, bam ] })` call.

### 7e. Update `CLAUDE.md` samplesheet schema table

Add row:

```markdown
| `adapter_kit` | *Optional* | MAS-seq adapter kit for Skera splitting (`mas8`, `mas12`, `mas16`, or custom FASTA path). Defaults to `mas16` if empty. | `mas8` |
```

### 7f. Remove all remaining `mas16_primers` references

```bash
grep -rn 'mas16_primers' modules/ subworkflows/ main.nf nextflow.config
# should return nothing after rename
```

**Verification:**
```bash
grep -rn 'mas16_primers' modules/ subworkflows/ main.nf nextflow.config CLAUDE.md
# Expected: no output

grep -rn 'adapter_primers\|adapter_file\|adapter_kit' \
    modules/ subworkflows/ main.nf nextflow.config CLAUDE.md
# Expected: nextflow.config (param), main.nf (resolution block + meta),
#           preprocess.nf module (input), preprocess.nf subworkflow (call), CLAUDE.md (schema)
```

---

## Phase 6 — Final verification

```bash
# 1. No traces of nanostat or mosdepth anywhere
grep -rn 'nanostat\|NANOSTAT\|NanoStat\|container_nanostat\|mosdepth\|MOSDEPTH\|container_mosdepth' \
    modules/ subworkflows/ main.nf nextflow.config
# Expected: no output

# 2. cramino wired at all expected points (now 5 preprocess stages: raw, segmented, filtered, fl, fltnc)
grep -rn 'CRAMINO\|cramino' \
    modules/ subworkflows/ main.nf nextflow.config
# Expected: qc.nf (process def), nextflow.config (container param),
#           preprocess.nf (import + 5 calls + emit + versions),
#           align.nf (import + 3 calls + emit + versions),
#           main.nf (2 mix lines)

# 3. Length filter wired
grep -rn 'SAMTOOLS_LENGTH_FILTER\|min_read_length' \
    modules/ subworkflows/ nextflow.config
# Expected: process def, import+call+versions, param

# 4. Dry-run syntax check
nextflow run main.nf -profile standard --samplesheet assets/example_samplesheet.csv -preview 2>&1 | tail -20
```

---

## Known gaps / follow-up

- **MultiQC cramino + mosdepth parsing:** cramino text output won't be rendered as a named
  section in the MultiQC HTML (unlike NanoStat and MosDepth which have native parsers).
  Files are published to disk and collected by MultiQC but will appear as unrecognised
  content. Check if a cramino MultiQC plugin ships with a newer MultiQC version; if so,
  update `container_multiqc`. Coverage depth summary (previously from MosDepth) is no
  longer in MultiQC — cramino karyotype provides per-chromosome read counts instead.
- **CRAM stage:** cramino can also run on the final CRAM (`GENERATE_CRAM.out.cram`) with
  `--spliced --karyotype --reference`. Omitted here because it requires the reference FASTA
  path and the CRAM is an archive copy — the aligned BAM cramino report covers the same
  reads. Add if per-archive metrics are needed.
- **Container tag:** verify `cramino:0.14.5--hdbdd923_0` is available on the HPC registry
  before running. Adjust tag as needed.
