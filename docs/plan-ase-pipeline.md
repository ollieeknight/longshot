# Plan: Allele-Specific Isoform Expression (ASE) Sub-Pipeline

**Goal:** Extend `longshot` with a haplotype-aware sub-pipeline that produces per-haplotype isoform counts and variant-level ASE statistics per NK cell library, branching off the existing aligned BAMs without modifying the main quantification path.

**Branch point:** `ALIGN.out.aligned_bam` → `[meta, bam, bai]` (per library, CB+XM tagged, coordinate-sorted)

---

## Phase 0: Documentation Discovery (COMPLETED)

### Allowed APIs / Tool Versions

| Tool | Version | Container | Key flag |
|---|---|---|---|
| DeepVariant | 1.8.0 | `google/deepvariant:1.8.0` | `--model_type PACBIO` |
| LongPhase | 2.0.1 | binary or build from GitHub | `--pb --indels` |
| WhatsHap | 2.8 | `quay.io/biocontainers/whatshap:2.3--py39hc16433a_1` | haplotag only |
| LORALS | latest | pip install from GitHub | `calc_ase`, `calc_asts` |
| IsoQuant | 3.13.0 | existing `container_isoquant` | `--read_group barcode` |

### Key Design Decisions (from research)

1. **LongPhase for phasing, WhatsHap for haplotagging.** LongPhase is 10× faster with better phase-block N50 on HiFi. WhatsHap `haplotag` adds the HP BAM tag that downstream tools need.
2. **No tool natively does single-cell + haplotype-aware isoform quantification.** Strategy: haplotag BAM → split by HP tag → run IsoQuant per haplotype BAM (H1, H2) using existing `container_isoquant`.
3. **LORALS for variant-level ASE.** Bulk-focused but accepts phased BAM + phased VCF. Run on the unsplit haplotagged BAM; LORALS reads the HP tag internally.
4. **DeepVariant is per-library** (on PBMM2_ALIGN output BAM). More reads per library → better phasing, so consider optional per-experiment merging in Phase 6.
5. **CB tags survive haplotagging and BAM splitting** — samtools writes HP tag at read level; CB is a separate read tag, untouched.

### Anti-Patterns to Avoid

- Do NOT use `google/deepvariant` with `--model_type PACBIO_HIFI` — the correct flag is `--model_type PACBIO`.
- Do NOT run WhatsHap `phase` — use LongPhase for phasing. WhatsHap is only used for `haplotag`.
- Do NOT split BAM by haplotype before haplotagging — HP tags must be added first.
- Do NOT use `--ignore-linked-read` in WhatsHap haplotag — this flag disables BX barcode handling (irrelevant here but avoid accidentally disabling CB handling).
- Do NOT use `isoquant --read_group tag:HP` on the unsplit BAM for haplotype quantification — IsoQuant cannot stratify output by (CB, HP) pairs natively.

---

## Phase 1: Variant Calling with DeepVariant

### What to implement

New process `DEEPVARIANT_CALL` in `modules/ase.nf`.

**Input:** `[meta, bam, bai]` from `ALIGN.out.aligned_bam`  
**Output:** `[meta, vcf_gz, tbi]`

**CLI pattern (from DeepVariant 1.8 PacBio case study):**
```bash
/opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
    --ref ${params.ref_fasta} \
    --reads ${bam} \
    --output_vcf ${meta.sample_id}_deepvariant.vcf.gz \
    --num_shards ${task.cpus} \
    --intermediate_results_dir tmp_dv
```

**Container:** `google/deepvariant:1.8.0` — add to `nextflow.config` as `container_deepvariant`.

**Label:** `process_high` (CPU-intensive; set cpus = 16, memory = 64 GB).

**publishDir:** `${params.outdir}/${meta.experiment}/${meta.library_id}/ase/variants/`

### Verification checklist
- [ ] Output VCF is bgzip-compressed (`.vcf.gz`) with tabix index (`.tbi`)
- [ ] Grep `run_deepvariant` in process script — confirm `--model_type PACBIO` (not `PACBIO_HIFI`)
- [ ] `nextflow.config` has `container_deepvariant = "google/deepvariant:1.8.0"`

---

## Phase 2: Phasing (LongPhase → WhatsHap haplotag)

### Phase 2a: LONGPHASE_PHASE

New process in `modules/ase.nf`.

**Input:** `[meta, vcf_gz, tbi, bam, bai]` (join DeepVariant output with original BAM by meta)  
**Output:** `[meta, phased_vcf]`

**CLI pattern (from LongPhase 2.0 README, HiFi mode):**
```bash
longphase phase \
    --snp-file ${vcf_gz} \
    --bam-file ${bam} \
    --reference ${params.ref_fasta} \
    --threads ${task.cpus} \
    --out-prefix ${meta.sample_id}_longphase \
    --pb \
    --indels
```

Output file: `${meta.sample_id}_longphase.vcf` (uncompressed; bgzip+index after).

**Container:** `container_longphase` — build or pull LongPhase binary container. Check `quay.io/biocontainers/longphase` for availability; fall back to `FROM ubuntu` + binary download if not on bioconda.

**Label:** `process_medium`

### Phase 2b: WHATSHAP_HAPLOTAG

New process in `modules/ase.nf`.

**Input:** `[meta, phased_vcf_gz, tbi, bam, bai]`  
**Output:** `[meta, haplotagged_bam, haplotagged_bai]`

**CLI pattern (from WhatsHap 2.8 docs):**
```bash
whatshap haplotag \
    --reference ${params.ref_fasta} \
    --output ${meta.sample_id}_haplotagged.bam \
    --output-threads ${task.cpus} \
    ${phased_vcf_gz} \
    ${bam}

samtools index ${meta.sample_id}_haplotagged.bam
```

**Container:** `quay.io/biocontainers/whatshap:2.3--py39hc16433a_1` — add as `container_whatshap`. (Check bioconda for 2.8 image; use latest available.)

**Label:** `process_medium`

### Verification checklist
- [ ] LongPhase output VCF has phased GT fields (`0|1`, `1|0` notation)
- [ ] WhatsHap output BAM has HP tag: `samtools view -h output.bam | head -100 | grep -c "HP:i"`
- [ ] CB tags still present after haplotagging: `samtools view output.bam | head | grep -o "CB:Z:[^ ]*"`

---

## Phase 3: Haplotype BAM Splitting

### What to implement

New process `SPLIT_BY_HAPLOTYPE` in `modules/ase.nf`.

**Input:** `[meta, haplotagged_bam, haplotagged_bai]`  
**Output:** `[meta, h1_bam, h1_bai, h2_bam, h2_bai]`

**Logic:** samtools `-d HP:1` / `-d HP:2` filter on the HP tag. Untagged reads (HP absent, e.g. reads spanning non-heterozygous regions) are excluded from per-haplotype BAMs but are available in the unsplit BAM for LORALS.

**CLI pattern:**
```bash
samtools view -b -d HP:1 -@ ${task.cpus} \
    ${haplotagged_bam} > ${meta.sample_id}_H1.bam
samtools index ${meta.sample_id}_H1.bam

samtools view -b -d HP:2 -@ ${task.cpus} \
    ${haplotagged_bam} > ${meta.sample_id}_H2.bam
samtools index ${meta.sample_id}_H2.bam
```

**Container:** `container_samtools` (already defined)

**Label:** `process_medium`

**Note on `-d` flag:** Available in samtools ≥ 1.12. Current pipeline uses `samtools:1.23.1` — confirmed compatible.

### Verification checklist
- [ ] H1 and H2 BAMs non-empty: `samtools view -c H1.bam` > 0
- [ ] CB tags present in split BAMs
- [ ] No HP tag in output reads (reads are split, not filtered on value presence)
- [ ] Sanity: `samtools view -c H1.bam` + `samtools view -c H2.bam` ≤ `samtools view -c haplotagged.bam` (unphased reads excluded)

---

## Phase 4: Per-Haplotype IsoQuant Quantification

### What to implement

New process `ISOQUANT_QUANTIFY_HAPLOTYPE` in `modules/ase.nf`. Nearly identical to existing `ISOQUANT_QUANTIFY` in `modules/quantify.nf` but takes a haplotype label and runs on H1/H2 BAMs.

**Input:** `[meta, haplotype, bam, bai, filtered_gtf]`  
- `haplotype`: `"H1"` or `"H2"` string  
**Output:** `[meta, haplotype, counts_dir]`

**CLI pattern:** Copy from `modules/quantify.nf:ISOQUANT_QUANTIFY` exactly, replacing output prefix:
```bash
isoquant \
    --reference ${params.ref_fasta} \
    --genedb ${filtered_gtf} \
    --complete_genedb \
    --process_only_chr ${chrs} \
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

**Container:** `container_isoquant` (already defined — reuse)

**publishDir:** `${params.outdir}/${meta.experiment}/${meta.library_id}/ase/isoquant_${haplotype}/`

**Channel construction in subworkflow:** Emit H1 and H2 as a mixed channel with haplotype label:
```groovy
SPLIT_BY_HAPLOTYPE.out
    .flatMap { meta, h1, h1i, h2, h2i ->
        [ [meta, "H1", h1, h1i], [meta, "H2", h2, h2i] ]
    }
    .combine(ch_filtered_gtf, by: 0)  // ch_filtered_gtf keyed by meta.experiment
    .map { meta, hap, bam, bai, gtf -> [meta, hap, bam, bai, gtf] }
    | ISOQUANT_QUANTIFY_HAPLOTYPE
```

### Verification checklist
- [ ] Two output directories per library: `*_H1/` and `*_H2/`
- [ ] Both contain `*.transcript_counts.tsv` and `*.gene_counts.tsv`
- [ ] H1 + H2 transcript counts roughly sum to total (allow ~20% loss from unphased reads)

---

## Phase 5: Variant-Level ASE with LORALS

### What to implement

Three chained processes in `modules/ase.nf`:
1. `LORALS_CALC_ASE` — allelic coverage per variant
2. `LORALS_ANNOTATE_ASE` — QC filter (indel ratio, mapping quality, blacklists)
3. `LORALS_CALC_ASTS` — transcript-level allele counts

All take the **unsplit haplotagged BAM** + phased VCF (more reads than either H1/H2 split).

**Input for LORALS_CALC_ASE:** `[meta, haplotagged_bam, haplotagged_bai, phased_vcf_gz, tbi]`

**CLI patterns (from LORALS GitHub):**
```bash
# Step 1
lorals calc_ase \
    -f ${phased_vcf_gz} \
    -b ${haplotagged_bam} \
    -o ${meta.sample_id}_ase.tsv

# Step 2
lorals annotate_ase \
    -f ${phased_vcf_gz} \
    -a ${meta.sample_id}_ase.tsv \
    -o ${meta.sample_id}_ase_annotated.tsv

# Step 3
lorals calc_asts \
    -f ${phased_vcf_gz} \
    -b ${haplotagged_bam} \
    -g ${params.ref_gtf} \
    -o ${meta.sample_id}_asts.tsv
```

**Container:** `container_lorals` — build from `pip install git+https://github.com/LappalainenLab/lorals.git` in a Python 3 base image. Add to `nextflow.config`.

**publishDir:** `${params.outdir}/${meta.experiment}/${meta.library_id}/ase/lorals/`

**Label:** `process_medium`

### Verification checklist
- [ ] `*_ase.tsv` has columns: `contig`, `position`, `variantID`, `refCount`, `altCount`
- [ ] `*_ase_annotated.tsv` row count ≤ `*_ase.tsv` (filtered)
- [ ] `*_asts.tsv` has transcript-level allele counts with gene/transcript IDs

---

## Phase 6: Subworkflow Assembly and main.nf Integration

### What to implement

**New file:** `subworkflows/ase.nf`

```groovy
include { DEEPVARIANT_CALL         } from '../modules/ase'
include { LONGPHASE_PHASE          } from '../modules/ase'
include { WHATSHAP_HAPLOTAG        } from '../modules/ase'
include { SPLIT_BY_HAPLOTYPE       } from '../modules/ase'
include { ISOQUANT_QUANTIFY_HAPLOTYPE } from '../modules/ase'
include { LORALS_CALC_ASE          } from '../modules/ase'
include { LORALS_ANNOTATE_ASE      } from '../modules/ase'
include { LORALS_CALC_ASTS         } from '../modules/ase'

workflow ASE {
    take:
    ch_aligned_bam   // [meta, bam, bai] from ALIGN
    ch_filtered_gtf  // [experiment, gtf] from CLASSIFY

    main:
    DEEPVARIANT_CALL(ch_aligned_bam)

    // Join VCF back with original BAM for phasing
    ch_aligned_bam
        .join(DEEPVARIANT_CALL.out.vcf)
        | LONGPHASE_PHASE

    // Join phased VCF with original BAM for haplotagging
    ch_aligned_bam
        .join(LONGPHASE_PHASE.out.phased_vcf)
        | WHATSHAP_HAPLOTAG

    SPLIT_BY_HAPLOTYPE(WHATSHAP_HAPLOTAG.out.haplotagged_bam)

    // Haplotype IsoQuant
    SPLIT_BY_HAPLOTYPE.out.haplotype_bams
        .flatMap { meta, h1, h1i, h2, h2i ->
            [ [meta.experiment, meta, "H1", h1, h1i],
              [meta.experiment, meta, "H2", h2, h2i] ]
        }
        .combine(ch_filtered_gtf, by: 0)
        .map { exp, meta, hap, bam, bai, gtf -> [meta, hap, bam, bai, gtf] }
        | ISOQUANT_QUANTIFY_HAPLOTYPE

    // LORALS on haplotagged BAM
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

**Edit `main.nf`:** After line ~150 (the `QUANTIFY` call), add:
```groovy
ASE(ALIGN.out.aligned_bam, CLASSIFY.out.filtered_gtf)
```

**Edit `nextflow.config`:** Add to params block:
```groovy
container_deepvariant = "google/deepvariant:1.8.0"
container_longphase   = "quay.io/biocontainers/longphase:2.0--h9ee0642_0"  // verify on bioconda
container_whatshap    = "quay.io/biocontainers/whatshap:2.3--py39hc16433a_1"
container_lorals      = null  // build custom SIF from pip install
```

**Optional future extension:** Add `params.run_ase = false` flag and guard the ASE workflow call with `if (params.run_ase)` — prevents mandatory DeepVariant cost on every run.

### Verification checklist
- [ ] Dry run succeeds: `nextflow run main.nf -preview -profile standard --samplesheet assets/example_samplesheet.csv`
- [ ] Channel graph shows ASE branch is parallel to QUANTIFY, not serial
- [ ] All processes emit `versions.yml`
- [ ] `publishDir` paths follow `${params.outdir}/${experiment}/${library_id}/ase/` convention

---

## Phase 7: Final Verification

### End-to-end check
1. Run `nextflow run main.nf -profile slurm --run_ase true` on one library
2. Confirm `results/<experiment>/<library>/ase/` tree contains:
   - `variants/` — DeepVariant VCF
   - `isoquant_H1/` and `isoquant_H2/` — per-haplotype counts
   - `lorals/` — `*_ase_annotated.tsv`, `*_asts.tsv`
3. Check H1+H2 isoform count totals sum to ~80% of full library counts (acceptable haplotype assignment rate)
4. Spot-check a heterozygous SNP: allele counts in `*_ase_annotated.tsv` should be roughly 50:50 in diploid regions

### Anti-pattern grep checks
```bash
grep -n "PACBIO_HIFI" modules/ase.nf         # must return nothing
grep -n "whatshap phase" modules/ase.nf       # must return nothing
grep -n "ignore-linked-read" modules/ase.nf   # must return nothing
grep -n "versions.yml" modules/ase.nf         # must return one match per process
```

---

## Open Questions / Future Extensions

| Question | Notes |
|---|---|
| Per-experiment BAM merge before DeepVariant? | More reads → better phasing. Add optional `SAMTOOLS_MERGE_EXPERIMENT` step before `DEEPVARIANT_CALL` if per-library coverage is insufficient (<10x). |
| SQANTI3 rescue integration? | Rescue should run before ASE since ASE uses the CLASSIFY-filtered GTF. See plan in `docs/plan-rescue-pipeline.md` (if created). |
| Population phasing? | WhatsHap `phase` supports population-level phasing with a reference panel (e.g., 1000G). Could improve short phase blocks in low-coverage libraries. |
| Cell-type-resolved ASE | Post-processing: group cells by Seurat cluster label, subset H1/H2 BAMs by CB whitelist per cluster, re-quantify. Not in pipeline scope yet. |
| LORALS single-cell awareness | LORALS bulk output can be filtered post-hoc by CB tag using `samtools view -D CB:barcodes.txt`. Script wrapper could produce per-cluster ASE tables. |
