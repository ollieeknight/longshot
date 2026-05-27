# `longshot` Project Guide

Welcome to the `longshot` (formerly `longreadr`) PacBio single-cell long-read transcriptomics Nextflow pipeline. This document serves as a developer and runtime guide for building, modifying, and running the pipeline.

---

## 🚀 Execution & Command Reference

### 1. Dry Run / Testing
Always perform a preflight dry run to check imports and channel structure before submitting to the cluster:
```bash
nextflow run main.nf -preview -profile standard --samplesheet assets/example_samplesheet.csv
```

### 2. Standard Local Profile (Test Runs)
Runs the pipeline locally (processes are executed sequentially on the local machine):
```bash
nextflow run main.nf \
    -profile standard \
    --samplesheet samplesheet.csv \
    --ref_fasta reference/GRCh38.fasta \
    --ref_gtf reference/gencode.gtf \
    --cage_peaks reference/refTSS.bed \
    --polya_list reference/polyA.list.txt
```

### 3. Production HPC Profile (SLURM + Apptainer)
Submits processes to the SLURM scheduler and executes inside Apptainer containers:
```bash
nextflow run main.nf \
    -profile slurm \
    --samplesheet samplesheet.csv \
    --ref_fasta reference/GRCh38.fasta \
    --ref_gtf reference/gencode.gtf \
    --cage_peaks reference/refTSS.bed \
    --polya_list reference/polyA.list.txt \
    --container_sqanti3 path/to/sqanti3_latest.sif
```

---

## 📂 Samplesheet Column Guide

The pipeline expects a comma-separated samplesheet (`.csv`) with the following schema:

| Column | Status | Purpose | Example |
|---|---|---|---|
| `experiment` | **Required** | Cohort or batch label (used for joint modeling) | `nk_activation` |
| `library_id` | **Required** | Unique library or sample identifier | `libA` |
| `run_id` | **Required** | SMRT cell run identifier | `m84094_260226_200655_s1` |
| `bam` | **Required** | Path to the raw PacBio multiplexed BAM | `/data/run1.bam` |
| `10x_index` | *Optional* | 10x sample index (for conditional demultiplexing) | `SI-GA-A1` |
| `chemistry` | *Optional* | 10x chemistry kit type (`3prime` or `5prime`, auto-detected if empty) | `5prime` |
| `shortread_barcodes` | *Optional* | Matched Illumina barcodes TSV path (for cell correction) | `libA_barcodes.tsv.gz` |

* **Fallback behavior:** If `10x_index`, `chemistry`, or `shortread_barcodes` are empty, blank, or set to `NULL`, `null`, `NA`, or `none`, the pipeline automatically bypasses demultiplexing / custom whitelisting and falls back to static global defaults (with `chemistry` defaulting to `params.chemistry` which is `3prime`).

---

## 📂 Codebase Directory Layout

```
longreadr/
  ├── main.nf                 ← Pipeline entry point, validation, and orchestrator
  ├── nextflow.config         ← Cluster-specific, profile-scoped configurations
  ├── assets/                 ← Whitelists, index CSV kits, primers, and samplesheets
  │   ├── indexes/            ← 10x Single and Dual index kit mapping tables
  │   ├── primers/            ← Standard 10x 3' and 5' cDNA primer FASTAs
  │   └── barcodes/           ← 10x cell barcode whitelists (static)
  ├── subworkflows/           ← Multi-step preprocessing, alignment, and quantification flows
  │   ├── preprocess.nf       ← Dynamic index demultiplexing, tags, and polyA trims
  │   ├── align.nf            ← Empirical whitelists, sorting, dedup, and mapping
  │   └── quantify.nf         ← Joint modeling, SQANTI3 filters, and TSV exports
  └── modules/                ← Single-tool process definitions
        ├── preprocess.nf
        ├── align.nf
        ├── quantify.nf
        ├── qc.nf             ← Flagstat, mosdepth, nanostat, and MultiQC
        └── exporter.nf       ← 10x QC matrix, saturation curve, and shared catalog
```

---

## 💻 Developer & Contribution Guidelines

### 1. Nextflow DSL2 Best Practices
* Keep all processes **isolated** in `modules/` and group multi-step logic in `subworkflows/`.
* Every process **must emit a `versions.yml`** file in standard nf-core format.
* Avoid absolute paths at the top level of `nextflow.config`; always scope cluster-specific paths inside the `slurm` profile.
* Use `mode: 'copy'` for light-weight publishing, and consider `mode: 'link'` or `mode: 'symlink'` when publishing massive CRAMs on the same filesystem.

---

## 🔭 Future Plans

### Allele-Specific Isoform Expression (ASE) Sub-Pipeline

Full implementation plan: `docs/plan-ase-pipeline.md`

**Goal:** Per-haplotype isoform quantification and variant-level ASE per NK cell library, branching from `ALIGN.out.aligned_bam` without touching the main quantification path.

**7-phase implementation stack:**

| Phase | Tool | Purpose |
| --- | --- | --- |
| 1 | DeepVariant 1.8.0 (`--model_type PACBIO`) | Per-library variant calling |
| 2a | LongPhase 2.0 (`--pb --indels`) | Haplotype phasing (faster + better N50 than WhatsHap for HiFi) |
| 2b | WhatsHap 2.8 `haplotag` | Add HP tags to reads (CB tags preserved) |
| 3 | `samtools view -d HP:1/2` | Split BAM into H1 and H2 (CB tags intact) |
| 4 | IsoQuant (reuse existing container) | Per-haplotype isoform quantification |
| 5 | LORALS (`calc_ase` → `annotate_ase` → `calc_asts`) | Variant + transcript-level ASE |
| 6 | `subworkflows/ase.nf` + `modules/ase.nf` | Pipeline integration |

**Key anti-patterns:**

* `--model_type PACBIO` not `PACBIO_HIFI`
* WhatsHap for haplotagging only — LongPhase does phasing
* Split BAM only after haplotagging
* IsoQuant cannot stratify by HP tag natively — must split BAM first

**Open question:** LORALS is bulk-only. Cell-type-resolved ASE (the real goal for NK subtypes) requires a post-processing step: subset haplotagged BAM by CB whitelist per Seurat cluster → per-cluster LORALS run. Not in current scope.

**New containers needed:** `container_deepvariant`, `container_longphase`, `container_whatshap`, `container_lorals`

**Reference files needed:** None beyond existing `ref_fasta` and `ref_gtf`. For full-mode rescue (if added later): run `sqanti3_qc.py` on reference transcriptome once to generate reference classification file.

---

### 2. GString Bash & Python Escaping Rules ⚠️
When writing Python scripts inside Nextflow's double-quoted GString blocks (`"""`):
* Always escape bash variables like `\$VAR` or `\$NBAM` with a backslash to prevent Groovy compilation errors.
* Always escape double backslashes in regex or escapes (e.g. use `\\t` and `\\n` instead of `\t` and `\n`) to ensure they are passed correctly to the shell environment.
* Avoid using python variables named like Groovy GString templates; python f-string `{var}` formatting is safe as long as it is not preceded by `$`.
* When writing nested python code, use standard indentation inside the double-quotes for legibility.
