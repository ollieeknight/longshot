# `longshot` Project Guide

Welcome to the `longshot` (formerly `longreadr`) PacBio single-cell long-read transcriptomics Nextflow pipeline. This document serves as a developer and runtime guide for building, modifying, and running the pipeline.

---

## Execution Environment

**All pipeline runs happen on the BIH HPC cluster (SLURM + Apptainer). Never run locally on macOS** — reference data and containers are only available on the cluster. The `standard` profile exists only for Nextflow syntax/channel validation via `-preview`; it does not run any actual processes.

---

## Execution & Command Reference

### 1. Dry Run / Syntax Check
Validates imports and channel structure on the head node. No jobs submitted, no containers needed:
```bash
nextflow run main.nf -preview -profile standard --samplesheet assets/example_samplesheet.csv
```

### 2. Production HPC Profile (SLURM + Apptainer)
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
| `adapter_kit` | *Optional* | MAS-seq adapter kit for Skera splitting: `mas8`, `mas12`, `mas16`, or a custom FASTA path. Defaults to `params.adapter_primers` (MAS16) if empty. | `mas8` |

* **Fallback behavior:** If `10x_index`, `chemistry`, `shortread_barcodes`, or `adapter_kit` are empty, blank, or set to `NULL`, `null`, `NA`, or `none`, the pipeline automatically bypasses demultiplexing / custom whitelisting and falls back to static global defaults (with `chemistry` defaulting to `params.chemistry` which is `3prime`, and `adapter_kit` defaulting to `params.adapter_primers` which is MAS16).

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
  │   ├── barcode_align.nf    ← Empirical whitelists, barcode correction, dedup, and mapping
  │   ├── classify.nf         ← Sharded IsoQuant discovery, SQANTI3 QC/filter/rescue
  │   ├── quantify.nf         ← Per-library IsoQuant quantification against filtered GTF
  │   └── export.nf           ← 10x-style MTX export, shared catalog, saturation curves
  ├── modules/                ← Single-tool process definitions
  │   ├── preprocess.nf
  │   ├── align.nf
  │   ├── quantify.nf
  │   ├── qc.nf             ← Flagstat, cramino (read-loss funnel), and MultiQC
  │   └── exporter.nf       ← 10x QC matrix, saturation curve, and shared catalog
  └── lib/
        └── genome_shards.nf ← Static chromosome-shard partition map used by classify.nf
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

**Not implemented — see `docs/plan-ase-pipeline.md` for the full plan.** No code from this section exists in `modules/` or `subworkflows/` yet.

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
