# Plan: ASE and Cluster-Resolved Secondary Analysis Sub-Pipelines

## Overview

**Goals:**

1. **Core ASE pipeline** — per-library, haplotype-resolved isoform quantification and variant-level allele-specific expression, branching from `ALIGN.out.aligned_bam` without touching the main quantification path.
2. **Cluster-resolved secondary analyses** — cell-type-stratified analyses (isoform quantification, ASE, somatic variant calling) driven by a `--from_secondary_analysis` flag that accepts cluster-barcode TSV files from upstream Seurat/Scanpy analysis.

**Architecture at a glance:**

```text
ALIGN.out.aligned_bam
    │
    ├─── [PART 1: Core ASE] ──────────────────────────────────────────────────
    │         │
    │    DEEPVARIANT_CALL → LONGPHASE_PHASE → WHATSHAP_HAPLOTAG
    │                                               │
    │                               SPLIT_BY_HAPLOTYPE (H1 / H2)
    │                                    │               │
    │                          ISOQUANT_H1         ISOQUANT_H2
    │                                   └─────────────────┘
    │                          LORALS_CALC_ASE (on haplotagged BAM)
    │
    └─── [PART 2: Cluster Secondary] ────────────────────────────────────────
              │
         (requires --from_secondary_analysis /path/to/cluster_tsvs/)
              │
    SAMTOOLS_SUBSET_CLUSTER (CB tag per cluster)
              │
       ┌──────┴──────────────┬──────────────────┬────────────────────┐
    Branch A            Branch B            Branch C            Branch D
    IsoQuant        LORALS ASE (per      DeepSomatic         Cluster × HP
    per cluster     cluster, req PART 1) somatic             IsoQuant
    → DIE (R)                            per cluster         (combinatorial)
```

**HPC execution note:** All pipeline development and testing happens on the BIH HPC cluster via SLURM. Never test locally on macOS — processes require containers (Apptainer) and the reference data lives on `/sc-projects/`. Always use `nextflow run main.nf -profile slurm ...` for test runs. Dry-run syntax: `nextflow run main.nf -preview -profile standard --samplesheet assets/example_samplesheet.csv` (head-node only, no jobs submitted).

---

## Phase 0: Documentation Discovery (COMPLETED)

### Tool versions and containers

| Tool | Version | Container | Purpose |
| --- | --- | --- | --- |
| DeepVariant | 1.8.0 | `quay.io/biocontainers/deepvariant:1.8.0--py310h5ef7bb4_1` | Germline variant calling |
| LongPhase | 2.0.1 | `quay.io/biocontainers/longphase:2.0--h9ee0642_0` | Haplotype phasing (HiFi) |
| WhatsHap | 2.8 | `quay.io/biocontainers/whatshap:2.3--py39hc16433a_1` | Haplotagging only (adds HP tag) |
| IsoQuant | 3.13.0 | `quay.io/biocontainers/isoquant:3.13.0--pyh106432d_0` | Isoform quantification (reuse existing) |
| LORALS | latest | custom SIF (pip install from GitHub) | Variant-level ASE |
| DeepSomatic | 1.9.0 | `google/deepsomatic:v1.9.0` | Somatic calling per cluster (tumor-normal paired mode) |
| samtools | 1.23.1 | `quay.io/biocontainers/samtools:1.23.1--ha83d96e_0` | BAM subsetting/splitting |

### Key design decisions

1. **DeepVariant for germline.** DeepVariant PACBIO mode achieves the highest SNP F1 (0.9968 at 30x) and best INDEL precision on PacBio HiFi. Actively maintained by Google; same tool family as DeepSomatic.
2. **PEPPER-Margin-DeepVariant is redundant here.** It adds a haplotagging round-trip that LongPhase already performs downstream. Use plain DeepVariant.
3. **LongPhase for phasing, WhatsHap for haplotagging.** LongPhase is 10× faster with better phase-block N50 on HiFi. WhatsHap `haplotag` adds the HP BAM tag that downstream tools need. Do not use WhatsHap `phase`.
4. **IsoQuant cannot stratify by HP tag natively.** Strategy: haplotag BAM → split by HP tag → run IsoQuant per haplotype BAM. CB tags survive haplotagging and splitting (separate read tag, untouched by HP operations).
5. **LORALS on unsplit haplotagged BAM.** More reads than either H1/H2 split; LORALS reads the HP tag internally for variant-level ASE.
6. **DeepSomatic for somatic calling (tumor-normal paired mode).** Mutect2 is broken for HiFi (`WellformedReadFilter` causes a JVM crash when disabled — `--disable-read-filter` is insufficient). DeepSomatic (Google, Nature Biotech 2025) is the choice: same tool family as DeepVariant, Parabricks GPU support, and — critically — **tumor-normal paired mode** (`--model_type PACBIO`, `--reads_tumor`, `--reads_normal`) which is more sensitive and specific than tumor-only. A `somatic_normal_bam` param (matched germline/PBMC BAM) is required; each cluster BAM is the tumor.
7. **CB tags in subset BAMs.** `samtools view -D CB:whitelist.txt` filters by the CB read tag; all other tags (HP, XM, etc.) are preserved.

### Global anti-patterns to avoid

| Anti-pattern | Correct approach |
| --- | --- |
| `--model_type PACBIO_HIFI` | Use `--model_type PACBIO` |
| `whatshap phase` | Use LongPhase for phasing; WhatsHap for haplotag only |
| Mutect2 on HiFi BAMs | Use DeepSomatic; Mutect2 WellformedReadFilter crashes on HiFi |
| `--disable-read-filter NotDuplicateReadFilter` with Mutect2 | Insufficient — WellformedReadFilter also blocks HiFi reads |
| Split BAM before haplotagging | Haplotag first; split after |
| `isoquant --read_group tag:HP` on unsplit BAM | IsoQuant cannot stratify by (CB, HP) pairs — split BAM first |
| Running `nextflow run` locally on macOS | All runs on HPC via SLURM |

### Reference files (all on HPC)

| Purpose | Path |
| --- | --- |
| Reference FASTA | `/sc-projects/sc-proj-cc12-ag-romagnani/ref/hs/pacbio/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` |
| Reference GTF | `/sc-projects/sc-proj-cc12-ag-romagnani/ref/hs/pacbio/gencode.v49.annotation.sorted.gtf` |
| CAGE peaks | `/sc-projects/sc-proj-cc12-ag-romagnani/ref/hs/pacbio/refTSS_v4.1_human_coordinate.hg38.sorted.bed` |
| gnomAD (somatic germline resource) | `/sc-projects/sc-proj-cc12-ag-romagnani/ref/hs/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz` |
| Panel of normals | `/sc-projects/sc-proj-cc12-ag-romagnani/ref/hs/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz` |

---

## PART 1: Core ASE Pipeline (Per-Library, Haplotype-Resolved)

**Module file:** `modules/ase.nf`  
**Subworkflow:** `subworkflows/ase.nf`  
**Param gate:** `params.run_ase = false`  
**Branch point:** `ALIGN.out.aligned_bam` → `[meta, bam, bai]`

---

### Phase 1: Germline Variant Calling

**New process `DEEPVARIANT_CALL` in `modules/ase.nf`.**

**Input:** `[meta, bam, bai]` from `ALIGN.out.aligned_bam`  
**Output:** `[meta, vcf_gz, tbi]`, `versions.yml`

```bash
run_deepvariant \
    --model_type PACBIO \
    --ref ${params.ref_fasta} \
    --reads ${bam} \
    --output_vcf ${meta.sample_id}_deepvariant.vcf.gz \
    --num_shards ${task.cpus} \
    --intermediate_results_dir tmp_dv
```

**Container:** `container_deepvariant = "quay.io/biocontainers/deepvariant:1.8.0--py310h5ef7bb4_1"`  
**Label:** `process_high` (cpus = 16, memory = 64 GB)  
**publishDir:** `${params.outdir}/${meta.experiment}/${meta.library_id}/ase/variants/`

**Verification:**

- [ ] Output `.vcf.gz` is bgzip-compressed with `.tbi` index
- [ ] `grep -c "PACBIO\b" modules/ase.nf` confirms `--model_type PACBIO` (not `PACBIO_HIFI`)

---

### Phase 2: Haplotype Phasing

#### Phase 2a: LongPhase (phasing)

**New process `LONGPHASE_PHASE` in `modules/ase.nf`.**

**Input:** `[meta, vcf_gz, tbi, bam, bai]` (DeepVariant output joined with original BAM by meta)  
**Output:** `[meta, phased_vcf_gz, tbi]`, `versions.yml`

```bash
longphase phase \
    --snp-file ${vcf_gz} \
    --bam-file ${bam} \
    --reference ${params.ref_fasta} \
    --threads ${task.cpus} \
    --out-prefix ${meta.sample_id}_longphase \
    --pb \
    --indels

bgzip ${meta.sample_id}_longphase.vcf
tabix ${meta.sample_id}_longphase.vcf.gz
```

**Container:** `container_longphase = "quay.io/biocontainers/longphase:2.0--h9ee0642_0"`  
**Label:** `process_medium`

#### Phase 2b: WhatsHap (haplotagging)

**New process `WHATSHAP_HAPLOTAG` in `modules/ase.nf`.**

**Input:** `[meta, phased_vcf_gz, tbi, bam, bai]`  
**Output:** `[meta, haplotagged_bam, haplotagged_bai]`, `versions.yml`

```bash
whatshap haplotag \
    --reference ${params.ref_fasta} \
    --output ${meta.sample_id}_haplotagged.bam \
    --output-threads ${task.cpus} \
    ${phased_vcf_gz} \
    ${bam}

samtools index ${meta.sample_id}_haplotagged.bam
```

**Container:** `container_whatshap = "quay.io/biocontainers/whatshap:2.3--py39hc16433a_1"`  
**Label:** `process_medium`

**Verification:**

- [ ] Phased VCF has `0|1` / `1|0` GT notation (not `0/1`)
- [ ] `samtools view -h haplotagged.bam | grep -c "HP:i"` > 0
- [ ] CB tags still present: `samtools view haplotagged.bam | head | grep -o "CB:Z:[^ ]*"`

---

### Phase 3: Haplotype BAM Splitting

**New process `SPLIT_BY_HAPLOTYPE` in `modules/ase.nf`.**

**Input:** `[meta, haplotagged_bam, haplotagged_bai]`  
**Output:** `[meta, h1_bam, h1_bai, h2_bam, h2_bai]`, `versions.yml`

```bash
samtools view -b -d HP:1 -@ ${task.cpus} ${haplotagged_bam} > ${meta.sample_id}_H1.bam
samtools index ${meta.sample_id}_H1.bam

samtools view -b -d HP:2 -@ ${task.cpus} ${haplotagged_bam} > ${meta.sample_id}_H2.bam
samtools index ${meta.sample_id}_H2.bam
```

`-d HP:1/2` requires samtools ≥ 1.12. Pipeline uses 1.23.1 — confirmed compatible.

**Container:** `container_samtools` (existing)  
**Label:** `process_medium`

**Verification:**

- [ ] `samtools view -c H1.bam` > 0 and `samtools view -c H2.bam` > 0
- [ ] H1 + H2 read count ≤ haplotagged BAM count (unphased reads excluded — expected)
- [ ] CB tags present in both split BAMs

---

### Phase 4: Per-Haplotype IsoQuant Quantification

**New process `ISOQUANT_QUANTIFY_HAPLOTYPE` in `modules/ase.nf`.**

**Input:** `[meta, haplotype, bam, bai, filtered_gtf]` where `haplotype` is `"H1"` or `"H2"`  
**Output:** `[meta, haplotype, counts_dir]`, `versions.yml`

```bash
isoquant \
    --reference ${params.ref_fasta} \
    --genedb ${filtered_gtf} \
    --complete_genedb \
    --bam ${bam} \
    --data_type pacbio \
    --barcoded_bam \
    --barcode_tag CB \
    --umi_tag XM \
    --read_group barcode \
    --threads ${task.cpus} \
    --prefix ${meta.library_id}_${haplotype} \
    --output isoquant_out
```

**Container:** `container_isoquant` (existing — reuse)  
**publishDir:** `${params.outdir}/${meta.experiment}/${meta.library_id}/ase/isoquant_${haplotype}/`

**Channel construction (in subworkflow):**

```groovy
SPLIT_BY_HAPLOTYPE.out
    .flatMap { meta, h1, h1i, h2, h2i ->
        [ [meta, "H1", h1, h1i], [meta, "H2", h2, h2i] ]
    }
    .combine(ch_filtered_gtf, by: 0)
    .map { meta, hap, bam, bai, gtf -> [meta, hap, bam, bai, gtf] }
    | ISOQUANT_QUANTIFY_HAPLOTYPE
```

**Verification:**

- [ ] Two output dirs per library: `*_H1/` and `*_H2/`
- [ ] Both contain `*.transcript_counts.tsv` and `*.gene_counts.tsv`
- [ ] H1 + H2 transcript counts sum to ~80% of full library (acceptable haplotype assignment rate)

---

### Phase 5: Variant-Level ASE (LORALS)

**Three chained processes in `modules/ase.nf`:**

All run on the **unsplit haplotagged BAM** (more reads than H1/H2 split; LORALS reads HP tag internally).

**Input for `LORALS_CALC_ASE`:** `[meta, haplotagged_bam, haplotagged_bai, phased_vcf_gz, tbi]`

```bash
# LORALS_CALC_ASE
lorals calc_ase \
    -f ${phased_vcf_gz} \
    -b ${haplotagged_bam} \
    -o ${meta.sample_id}_ase.tsv

# LORALS_ANNOTATE_ASE
lorals annotate_ase \
    -f ${phased_vcf_gz} \
    -a ${meta.sample_id}_ase.tsv \
    -o ${meta.sample_id}_ase_annotated.tsv

# LORALS_CALC_ASTS
lorals calc_asts \
    -f ${phased_vcf_gz} \
    -b ${haplotagged_bam} \
    -g ${params.ref_gtf} \
    -o ${meta.sample_id}_asts.tsv
```

**Container:** `container_lorals` — build custom SIF from `pip install git+https://github.com/LappalainenLab/lorals.git`  
**publishDir:** `${params.outdir}/${meta.experiment}/${meta.library_id}/ase/lorals/`  
**Label:** `process_medium`

**Verification:**

- [ ] `*_ase.tsv` has columns: `contig`, `position`, `variantID`, `refCount`, `altCount`
- [ ] `*_ase_annotated.tsv` row count ≤ `*_ase.tsv` (filtered)
- [ ] `*_asts.tsv` has transcript-level allele counts with gene/transcript IDs

---

### Phase 6: Core Subworkflow Assembly

**New file: `subworkflows/ase.nf`**

```groovy
include { DEEPVARIANT_CALL            } from '../modules/ase'
include { LONGPHASE_PHASE             } from '../modules/ase'
include { WHATSHAP_HAPLOTAG           } from '../modules/ase'
include { SPLIT_BY_HAPLOTYPE          } from '../modules/ase'
include { ISOQUANT_QUANTIFY_HAPLOTYPE } from '../modules/ase'
include { LORALS_CALC_ASE             } from '../modules/ase'
include { LORALS_ANNOTATE_ASE         } from '../modules/ase'
include { LORALS_CALC_ASTS            } from '../modules/ase'

workflow ASE {
    take:
    ch_aligned_bam   // [meta, bam, bai] from ALIGN
    ch_filtered_gtf  // [experiment, gtf] from CLASSIFY

    main:
    DEEPVARIANT_CALL(ch_aligned_bam)

    ch_aligned_bam
        .join(DEEPVARIANT_CALL.out.vcf)
        | LONGPHASE_PHASE

    ch_aligned_bam
        .join(LONGPHASE_PHASE.out.phased_vcf)
        | WHATSHAP_HAPLOTAG

    SPLIT_BY_HAPLOTYPE(WHATSHAP_HAPLOTAG.out.haplotagged_bam)

    SPLIT_BY_HAPLOTYPE.out.haplotype_bams
        .flatMap { meta, h1, h1i, h2, h2i ->
            [ [meta.experiment, meta, "H1", h1, h1i],
              [meta.experiment, meta, "H2", h2, h2i] ]
        }
        .combine(ch_filtered_gtf, by: 0)
        .map { exp, meta, hap, bam, bai, gtf -> [meta, hap, bam, bai, gtf] }
        | ISOQUANT_QUANTIFY_HAPLOTYPE

    WHATSHAP_HAPLOTAG.out.haplotagged_bam
        .join(LONGPHASE_PHASE.out.phased_vcf)
        | LORALS_CALC_ASE

    LORALS_CALC_ASE.out.ase_tsv
        .join(LONGPHASE_PHASE.out.phased_vcf)
        | LORALS_ANNOTATE_ASE

    WHATSHAP_HAPLOTAG.out.haplotagged_bam
        .join(LONGPHASE_PHASE.out.phased_vcf)
        | LORALS_CALC_ASTS

    emit:
    haplotagged_bam  = WHATSHAP_HAPLOTAG.out.haplotagged_bam
    phased_vcf       = LONGPHASE_PHASE.out.phased_vcf
    haplotype_counts = ISOQUANT_QUANTIFY_HAPLOTYPE.out.counts_dir
    ase_annotated    = LORALS_ANNOTATE_ASE.out.ase_annotated_tsv
    asts             = LORALS_CALC_ASTS.out.asts_tsv
    versions         = Channel.empty()
                        .mix(DEEPVARIANT_CALL.out.versions)
                        .mix(LONGPHASE_PHASE.out.versions)
                        .mix(WHATSHAP_HAPLOTAG.out.versions)
                        .mix(SPLIT_BY_HAPLOTYPE.out.versions)
                        .mix(ISOQUANT_QUANTIFY_HAPLOTYPE.out.versions)
                        .mix(LORALS_CALC_ASE.out.versions)
                        .mix(LORALS_ANNOTATE_ASE.out.versions)
                        .mix(LORALS_CALC_ASTS.out.versions)
}
```

**Edit `main.nf`** — after `QUANTIFY` call:

```groovy
if (params.run_ase) {
    ASE(ALIGN.out.aligned_bam, CLASSIFY.out.filtered_gtf)
}
```

**Add to `nextflow.config` params:**

```groovy
run_ase               = false
container_deepvariant = "quay.io/biocontainers/deepvariant:1.8.0--py310h5ef7bb4_1"
container_longphase   = "quay.io/biocontainers/longphase:2.0--h9ee0642_0"
container_whatshap    = "quay.io/biocontainers/whatshap:2.3--py39hc16433a_1"
container_lorals      = null  // custom SIF — path set at runtime
```

---

## PART 2: Cluster-Resolved Secondary Analyses

**Module file:** `modules/secondary.nf`  
**Subworkflow:** `subworkflows/secondary.nf`  
**Entry flag:** `--from_secondary_analysis /path/to/cluster_tsvs/`  
**Branch point:** aligned BAM (or haplotagged BAM for Branch B) → subset by CB tag per cluster → parallel analysis branches

---

### Phase 7: Input Design and Parameter Specification

#### File format

The `--from_secondary_analysis` flag points to a **directory of per-library TSV files**. Each file is named `{library_id}.tsv` and contains one row per cell with a header:

```tsv
cluster    barcode
NK56bright ACGTACGTACGTACGT-1
NK56bright TTGCTTGCTTGCTTGC-1
NK56dim    GCTAGCTAGCTAGCTA-1
```

- `cluster`: string label from Seurat/Scanpy; safe characters only (`[A-Za-z0-9_-]`); no spaces.
- `barcode`: CB tag value exactly as stored in the BAM (no `CB:Z:` prefix).
- One row per cell. Cells absent from the BAM are silently ignored by samtools.
- Typical export from Seurat: `write.table(data.frame(cluster=seu$seurat_clusters, barcode=colnames(seu)), "libA.tsv", sep="\t", quote=FALSE, row.names=FALSE)`

#### Params (add to `nextflow.config`)

```groovy
// ── Cluster-resolved secondary analyses ─────────────────────────────────────
secondary_analysis_dir         = null   // path to dir of {library_id}.tsv files
run_cluster_isoquant           = false  // Branch A: per-cluster IsoQuant → DIE
run_cluster_ase                = false  // Branch B: per-cluster LORALS (requires run_ase = true)
run_cluster_somatic            = false  // Branch C: DeepSomatic tumor-normal per cluster
run_cluster_haplotype_isoquant = false  // Branch D: cluster × haplotype (combinatorial)

somatic_normal_bam  = null  // matched germline/PBMC BAM for DeepSomatic tumor-normal mode (required for Branch C)
container_deepsomatic = "google/deepsomatic:v1.9.0"
```

#### Channel construction in `main.nf`

```groovy
if (params.secondary_analysis_dir) {
    ch_cluster_barcodes = Channel
        .fromPath("${params.secondary_analysis_dir}/*.tsv")
        .flatMap { f ->
            def lib_id = f.baseName
            f.readLines().drop(1)
             .collect { it.split('\t')[0] }.unique()
             .collect { cluster -> [lib_id, cluster, f] }
        }

    SECONDARY(
        ch_cluster_barcodes,
        ALIGN.out.aligned_bam,
        params.run_ase ? ASE.out.haplotagged_bam : Channel.empty(),
        params.run_ase ? ASE.out.phased_vcf      : Channel.empty(),
        CLASSIFY.out.filtered_gtf
    )
}
```

> `f.readLines()` is a Groovy `java.io.File` method. Safe inside `flatMap` on the head node. Do not call inside a process script.

**Verification:**

- [ ] Channel emits one item per unique (library_id, cluster) pair
- [ ] Library IDs in TSV filenames match `meta.library_id` values in `ch_aligned_bam`
- [ ] Cluster names with unsafe characters → validation error in `preflight_samplesheet` before any jobs submit

---

### Phase 8: Per-Cluster BAM Subsetting

**New process `SAMTOOLS_SUBSET_CLUSTER` in `modules/secondary.nf`.**

This is the common entry point for all cluster-resolved branches. Subsets the BAM to reads whose CB tag appears in the cluster's barcode whitelist.

**Input:** `[meta, cluster, bam, bai, cluster_tsv]`  
**Output:** `[meta, cluster, subset_bam, subset_bai]`, `versions.yml`

```bash
awk -F'\t' -v clust="${cluster}" 'NR>1 && $1==clust {print $2}' \
    ${cluster_tsv} > ${meta.library_id}_${cluster}_barcodes.txt

samtools view -b \
    -D CB:${meta.library_id}_${cluster}_barcodes.txt \
    -@ ${task.cpus} \
    ${bam} \
    > ${meta.library_id}_${cluster}.bam

samtools index ${meta.library_id}_${cluster}.bam
```

**Container:** `container_samtools` (existing)  
**Label:** `process_low`  
**publishDir:** none (intermediate; too large to publish)

**Channel join (in subworkflow):**

```groovy
// BAM source: haplotagged if run_cluster_ase, else aligned
def ch_base_bam = params.run_cluster_ase ? ch_haplotagged_bam : ch_aligned_bam

ch_cluster_barcodes
    .combine(
        ch_base_bam.map { meta, bam, bai -> [meta.library_id, meta, bam, bai] },
        by: 0
    )
    .map { lib_id, cluster, tsv, meta, bam, bai -> [meta, cluster, bam, bai, tsv] }
    | SAMTOOLS_SUBSET_CLUSTER
```

**Verification:**

- [ ] `samtools view -c subset.bam` > 0 for a known non-empty cluster
- [ ] All reads have CB tags matching the whitelist
- [ ] HP and XM tags intact if subset came from haplotagged BAM

---

### Phase 9: Analysis Branches

All branches take `SAMTOOLS_SUBSET_CLUSTER.out` as input and are independently gated by params.

---

#### Branch A: Per-Cluster Isoform Quantification → Differential Isoform Expression

**Param gate:** `params.run_cluster_isoquant = false`

**New process `ISOQUANT_QUANTIFY_CLUSTER` in `modules/secondary.nf`.**

**Input:** `[meta, cluster, bam, bai, filtered_gtf]`  
**Output:** `[meta, cluster, counts_dir]`, `versions.yml`

```bash
isoquant \
    --reference ${params.ref_fasta} \
    --genedb ${filtered_gtf} \
    --complete_genedb \
    --bam ${bam} \
    --data_type pacbio \
    --barcoded_bam \
    --barcode_tag CB \
    --umi_tag XM \
    --read_group barcode \
    --threads ${task.cpus} \
    --prefix ${meta.library_id}_${cluster} \
    --output isoquant_out
```

**Container:** `container_isoquant` (existing)  
**publishDir:** `${params.outdir}/${meta.experiment}/${meta.library_id}/secondary/${cluster}/isoquant/`

**Downstream (post-pipeline, in R):** Merge `*.transcript_counts.tsv` across clusters/libraries → DRIMSeq (transcript proportion testing) or satuRn (scalable isoform testing for single-cell) for differential isoform expression. DEXSeq for exon-level differential usage.

---

#### Branch B: Per-Cluster Allele-Specific Expression (LORALS)

**Param gate:** `params.run_cluster_ase = false`  
**Prerequisite:** `params.run_ase = true` — Branch B subsets the **haplotagged** BAM (not raw aligned BAM), so `WHATSHAP_HAPLOTAG` must have run.

**New process `LORALS_CALC_ASE_CLUSTER` in `modules/secondary.nf`.**

**Input:** `[meta, cluster, subset_bam, subset_bai, phased_vcf_gz, tbi]`  
**Output:** `[meta, cluster, ase_annotated_tsv]`, `versions.yml`

```bash
lorals calc_ase \
    -f ${phased_vcf_gz} \
    -b ${subset_bam} \
    -o ${meta.library_id}_${cluster}_ase.tsv

lorals annotate_ase \
    -f ${phased_vcf_gz} \
    -a ${meta.library_id}_${cluster}_ase.tsv \
    -o ${meta.library_id}_${cluster}_ase_annotated.tsv
```

**Container:** `container_lorals` (custom SIF, same as Phase 5)  
**publishDir:** `${params.outdir}/${meta.experiment}/${meta.library_id}/secondary/${cluster}/lorals/`

**Known limitation:** LORALS is bulk-focused. Recommend ≥200 cells per cluster for adequate allelic read depth. Check cell counts before enabling.

---

#### Branch C: Per-Cluster Somatic Variant Calling (DeepSomatic, tumor-normal)

**Param gate:** `params.run_cluster_somatic = false`  
**Prerequisite:** `params.somatic_normal_bam` must be set — path to a matched germline or PBMC BAM. Each cluster subset BAM is the tumor; the normal BAM is shared across all clusters.

DeepSomatic (Google, Nature Biotech 2025) in **tumor-normal paired mode** (`--model_type PACBIO`). This mode is more sensitive and specific than tumor-only. Do not use Mutect2 — `WellformedReadFilter` crashes on HiFi with no workaround.

**New process `DEEPSOMATIC_SOMATIC_CLUSTER` in `modules/secondary.nf`.**

**Input:** `[meta, cluster, tumor_bam, tumor_bai, normal_bam, normal_bai]`  
**Output:** `[meta, cluster, somatic_vcf_gz, tbi]`, `versions.yml`

```bash
run_deepsomatic \
    --model_type PACBIO \
    --ref ${params.ref_fasta} \
    --reads_tumor ${tumor_bam} \
    --reads_normal ${normal_bam} \
    --output_vcf ${meta.library_id}_${cluster}_somatic.vcf.gz \
    --num_shards ${task.cpus}
```

**Channel construction** — inject the normal BAM into each cluster item:

```groovy
SAMTOOLS_SUBSET_CLUSTER.out.subset_bam
    .map { meta, cluster, bam, bai ->
        [meta, cluster, bam, bai,
         file(params.somatic_normal_bam),
         file("${params.somatic_normal_bam}.bai")]
    }
    | DEEPSOMATIC_SOMATIC_CLUSTER
```

**Container:** `container_deepsomatic = "google/deepsomatic:v1.9.0"` — Docker Hub only (no bioconda). Pull to Apptainer SIF on HPC: `apptainer pull deepsomatic.sif docker://google/deepsomatic:v1.9.0`. GPU-accelerated via NVIDIA Clara Parabricks.  
**publishDir:** `${params.outdir}/${meta.experiment}/${meta.library_id}/secondary/${cluster}/somatic/`  
**Label:** `process_high` (32 GB minimum)

**Validation guard in `main.nf`** — fail early if normal BAM missing:

```groovy
if (params.run_cluster_somatic && !params.somatic_normal_bam) {
    error "--run_cluster_somatic requires --somatic_normal_bam <path>"
}
```

**Use case:** Identify cluster-enriched somatic mutations — clonal NK populations, clonal haematopoiesis, NK cell malignancies. Post-pipeline: compare per-cluster VAF distributions against the shared normal.

---

#### Branch D: Cluster × Haplotype IsoQuant (Most Granular)

**Param gate:** `params.run_cluster_haplotype_isoquant = false`

Produces one BAM per (library, cluster, haplotype) triplet. Combines Phase 8 (CB subsetting) with Phase 3 (HP splitting).

```text
ALIGN → WHATSHAP_HAPLOTAG → SAMTOOLS_SUBSET_CLUSTER (CB) → SPLIT_BY_HAPLOTYPE (HP) → ISOQUANT_QUANTIFY_HAPLOTYPE
```

**Combinatorial warning:** `n_libraries × n_clusters × 2 haplotypes` jobs. Example: 10 libraries × 8 clusters = 160 IsoQuant runs. Always check per-(cluster, haplotype) cell counts before enabling; recommend ≥100 cells per BAM. Add a pre-flight count process that reports cell counts and fails fast below threshold.

---

### Phase 10: Secondary Subworkflow Assembly

**New file: `subworkflows/secondary.nf`**

```groovy
include { SAMTOOLS_SUBSET_CLUSTER    } from '../modules/secondary'
include { ISOQUANT_QUANTIFY_CLUSTER  } from '../modules/secondary'
include { LORALS_CALC_ASE_CLUSTER    } from '../modules/secondary'
include { DEEPSOMATIC_SOMATIC_CLUSTER } from '../modules/secondary'

workflow SECONDARY {
    take:
    ch_cluster_barcodes   // [lib_id, cluster, tsv]
    ch_aligned_bam        // [meta, bam, bai] from ALIGN
    ch_haplotagged_bam    // [meta, bam, bai] from ASE.out (may be empty)
    ch_phased_vcf         // [meta, vcf_gz, tbi] from ASE.out (may be empty)
    ch_filtered_gtf       // [experiment, gtf] from CLASSIFY

    main:
    def ch_base_bam = params.run_cluster_ase ? ch_haplotagged_bam : ch_aligned_bam

    ch_for_subset = ch_cluster_barcodes
        .combine(
            ch_base_bam.map { meta, bam, bai -> [meta.library_id, meta, bam, bai] },
            by: 0
        )
        .map { lib_id, cluster, tsv, meta, bam, bai -> [meta, cluster, bam, bai, tsv] }

    SAMTOOLS_SUBSET_CLUSTER(ch_for_subset)

    if (params.run_cluster_isoquant) {
        SAMTOOLS_SUBSET_CLUSTER.out.subset_bam
            .map { meta, cluster, bam, bai -> [meta.experiment, meta, cluster, bam, bai] }
            .combine(ch_filtered_gtf, by: 0)
            .map { exp, meta, cluster, bam, bai, gtf -> [meta, cluster, bam, bai, gtf] }
            | ISOQUANT_QUANTIFY_CLUSTER
    }

    if (params.run_cluster_ase) {
        SAMTOOLS_SUBSET_CLUSTER.out.subset_bam
            .join(ch_phased_vcf, by: [0])
            | LORALS_CALC_ASE_CLUSTER
    }

    if (params.run_cluster_somatic) {
        SAMTOOLS_SUBSET_CLUSTER.out.subset_bam
            .map { meta, cluster, bam, bai ->
                [meta, cluster, bam, bai,
                 file(params.somatic_normal_bam),
                 file("${params.somatic_normal_bam}.bai")]
            }
            | DEEPSOMATIC_SOMATIC_CLUSTER
    }

    emit:
    cluster_isoquant_counts = params.run_cluster_isoquant
        ? ISOQUANT_QUANTIFY_CLUSTER.out.counts_dir : Channel.empty()
    cluster_ase             = params.run_cluster_ase
        ? LORALS_CALC_ASE_CLUSTER.out.ase_annotated_tsv : Channel.empty()
    cluster_somatic         = params.run_cluster_somatic
        ? DEEPSOMATIC_SOMATIC_CLUSTER.out.somatic_vcf : Channel.empty()
}
```

---

## PART 3: Verification and Testing

**All verification runs on the HPC via SLURM.** No local Mac testing — containers and reference files are only available on the cluster.

---

### Core ASE Verification

**Full run (one library):**

```bash
nextflow run main.nf \
    -profile slurm \
    --samplesheet assets/example_samplesheet.csv \
    --run_ase true
```

**Expected output tree:**

```text
results/<experiment>/<library>/
    ase/
        variants/    ← DeepVariant VCF
        isoquant_H1/ ← H1 transcript/gene counts
        isoquant_H2/ ← H2 transcript/gene counts
        lorals/      ← *_ase_annotated.tsv, *_asts.tsv
```

**Spot checks:**

- [ ] H1 + H2 isoform count totals ≈ 80% of full library counts (20% haplotype assignment loss acceptable)
- [ ] Heterozygous SNP allele counts in `*_ase_annotated.tsv` ≈ 50:50 in diploid regions
- [ ] `nextflow run main.nf -preview` channel graph shows ASE branch is parallel to QUANTIFY (not serial)

---

### Secondary Analysis Verification

**Full run with all branches:**

```bash
nextflow run main.nf \
    -profile slurm \
    --samplesheet assets/example_samplesheet.csv \
    --run_ase true \
    --from_secondary_analysis /path/to/cluster_tsvs/ \
    --run_cluster_isoquant true \
    --run_cluster_ase true \
    --run_cluster_somatic true
```

**Expected output tree:**

```text
results/<experiment>/<library>/
    secondary/<cluster>/
        isoquant/  ← per-cluster transcript counts
        lorals/    ← per-cluster *_ase_annotated.tsv
        somatic/   ← per-cluster *_somatic.vcf.gz
```

**Spot checks:**

- [ ] `samtools view -c results/.../secondary/NK56bright/*.bam` > 0
- [ ] No reads from cluster A in cluster B BAM: `samtools view -D CB:clusterB_barcodes.txt clusterA.bam | wc -l` → 0
- [ ] Cluster BAM cell counts match Seurat metadata (n cells per cluster)
- [ ] Running without `--from_secondary_analysis` → SECONDARY workflow absent from `-preview` graph

---

### Global Anti-Pattern Grep Checks

Run after implementing each module file:

```bash
# modules/ase.nf
grep -n "PACBIO_HIFI"         modules/ase.nf   # must be empty
grep -n "whatshap phase"      modules/ase.nf   # must be empty
grep -n "ignore-linked-read"  modules/ase.nf   # must be empty
grep -n "versions.yml"        modules/ase.nf   # one match per process

# modules/secondary.nf
grep -n "Mutect2"             modules/secondary.nf  # must be empty
grep -n "hifi_revio_ss\b"     modules/secondary.nf  # must be empty (use hifi_revio_ssrs)
grep -n "versions.yml"        modules/secondary.nf  # one match per process
```

---

## Appendix: Open Questions and Future Extensions

| Question | Status | Notes |
| --- | --- | --- |
| Per-experiment BAM merge before DeepVariant? | Open | More reads → better phasing. Add optional `SAMTOOLS_MERGE_EXPERIMENT` before `DEEPVARIANT_CALL` if per-library coverage <10x. |
| SQANTI3 rescue integration? | Open | Rescue should run before ASE — ASE uses CLASSIFY-filtered GTF. |
| Population phasing? | Open | LongPhase or WhatsHap `phase` with 1000G reference panel (`1000G_phase1.snps.high_confidence.hg38.vcf.gz` in GATKBundle) could improve short phase blocks. |
| Cell-type-resolved ASE | Planned | Phase 9 Branch B — LORALS per cluster. Recommend ≥200 cells/cluster. |
| Per-cluster somatic calling | Planned | Phase 9 Branch C — DeepSomatic `--model_type PACBIO` tumor-normal paired mode. Requires `somatic_normal_bam`. Not Mutect2. |
| Per-cluster isoform quantification | Planned | Phase 9 Branch A — IsoQuant per cluster → DRIMSeq/satuRn DIE in R. |
| Cluster × haplotype IsoQuant | Planned | Phase 9 Branch D — most granular. Check cell counts before enabling. |
| Differential isoform expression | Planned | Post-pipeline R: DRIMSeq or satuRn on per-cluster count matrices. |
| Trajectory-resolved ASE | Future | Bin cells by pseudotime, treat bins as pseudo-clusters, feed into `--from_secondary_analysis`. No new pipeline code needed. |
| LORALS depth floor | Open | <200 cells/cluster → underpowered ASE. Consider aggregating libraries before per-cluster LORALS. |
| DeepSomatic bioconda packaging | Monitor | Currently Docker Hub only. Watch for bioconda recipe — would simplify Apptainer SIF build. |
