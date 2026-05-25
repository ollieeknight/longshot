# `longshot` 🚀

> **A Production-Grade Nextflow DSL2 Pipeline for Single-Cell Long-Read Transcriptomics**

`longshot` is an end-to-end, automated Nextflow pipeline designed for processing PacBio HiFi MAS-seq single-cell transcriptomics data. It takes you seamlessly from raw instrument multiplexed runs to fully annotated, ready-to-import Single-Cell count matrices and cohort-wide structural catalogs.

---

## 🌟 Key Features & Highlights

### 🛡️ 1. Intelligent Pre-Flight Index Detector & Sanity Guard
Before launching heavy, expensive HPC compute jobs, `longshot` runs an ultra-fast **pre-flight index detector** (`DETECT_SAMPLE_INDICES`):
* Downsamples your raw BAM file to the first 20,000 reads (in seconds).
* Scans the cDNA adapter regions to empirically verify the presence of your specified 10x sample indices.
* **Wrong Index Guard (Fail-Fast):** Halts immediately and throws a descriptive error if the specified index (e.g., `SI-GA-A1`) is absent, preventing wasted cluster hours due to metadata typos.
* **Single-Index Auto-Override (Read Recovery):** If multiple indices are specified in the samplesheet but the detector finds that >99.5% of reads belong to a single index, it **automatically bypasses index demultiplexing** and runs standard cDNA trimming to maximize read recovery.

### 🔀 2. Automated & Conditional 10x Index Demultiplexing
If a raw flowcell contains multiple verified indices:
* The pipeline groups the run by BAM path to run `SKERA_SPLIT` **exactly once per SMRT cell**, saving massive processing overhead.
* Nextflow dynamically builds a custom multiplexed primers FASTA containing all 4 constituent barcode sequences per index.
* **Chemistry-Aware Kit Selection:** Automatically detects 5' GEX vs 3' GEX kits and modifies the correct primer (`5p` or `3p`) accordingly.
* Runs `lima` to separate the multiplexed libraries, and automatically merges the constituent index BAMs per library using `samtools merge`.

### 🧪 3. Paired Short-Read Barcode Correction (Empirical Whitelist)
Correcting cell barcodes against a static 3-million possible barcodes list can introduce noise by rescuing ambient reads to empty droplets. If you have paired short-read (Illumina) data:
* Provide the CellRanger `barcodes.tsv.gz` path in your samplesheet.
* The pipeline **auto-strips the `-1` gem-group suffix** and **reverse-complements the sequences** to match PacBio's sequencing chemistry orientation completely automatically.
* It uses this high-confidence empirical list for `isoseq correct`, strictly keeping reads belonging to genuine cells.
* If no short-read barcodes are specified, it gracefully falls back to the global static 10x whitelist.

### 📊 4. CellRanger-Style Library QC Exporter (Features Enriched with SQANTI3)
* Exports counts in standard 10x-style format (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`) per library for both **gene-level** and **transcript-level** matrices.
* **SQANTI3 Metadata Integration:** For transcripts, we enrich the `features.tsv.gz` by mapping Ensembl IDs to **Gene Symbols** and merging each transcript directly with its **SQANTI3 structural category** (e.g. `full-splice_match`, `novel_in_catalog`), allowing you to load counts directly into Seurat using standard `Read10X()` and perform instant downstream transcript-QC in R!
* **Seurat Safety:** Replaces underscores (`_`) with dashes (`-`) in transcript IDs to prevent Seurat from crashing.

### 📂 5. Cohort-Wide Shared Catalog
* Because `longshot` performs **joint isoform discovery** across all libraries, all libraries automatically share a unified structural catalog.
* Generates `shared_isoform_catalog.tsv.gz` (cleaned of NA values that crash Pigeon/Seurat) and `shared_isoform_map.tsv.gz` linking transcript IDs to gene names, gene symbols, and structural categories across all donors.

### 📉 6. Sequencing Saturation Curves
* Downsamples IsoQuant read assignments (`*.read_to_transcript.tsv`) in Python in seconds to compute saturation rates (for genes, isoforms, cells, and UMIs) at 10%, 25%, 50%, 75%, and 100% sequencing depth. Saturation plots are natively integrated into the final **MultiQC report**.

---

## 🛠️ Pipeline Architecture

```
main.nf
  ├── subworkflows/preprocess.nf  →  modules/preprocess.nf  (SKERA → DETECT_INDEX → DYNAMIC_PRIMERS → LIMA → MERGE → TAG → REFINE)
  ├── subworkflows/align.nf       →  modules/align.nf        (MERGE → PREP_WHITELIST → CORRECT → SORT_CB → DEDUP → PBMM2 → CRAM)
  ├── subworkflows/quantify.nf    →  modules/quantify.nf     (CB_SUFFIX → ISOQUANT_DISC → SQANTI3 → ISOQUANT_QUANT → EXPORTER)
  └── modules/qc.nf               (FLAGSTAT, MOSDEPTH, NANOSTAT, MULTIQC)
```

---

## 🚀 Quick Start

### 1. Configure your `samplesheet.csv`
The samplesheet is fully generalized to handle standard runs, multiplexed runs, and custom short-read whitelists:

```csv
experiment,library_id,10x_index,run_id,bam,shortread_barcodes
nk_activation,libA,SI-GA-A1,m84094_260226_200655_s1,/data/run1.bam,/data/illumina/libA/barcodes.tsv.gz
nk_activation,libB,SI-GA-A2,m84094_260226_200655_s1,/data/run1.bam,/data/illumina/libB/barcodes.tsv.gz
nk_activation,libC,SI-GA-A3,m84094_260226_220949_s1,/data/run2.bam,NULL
nk_activation,libD,SI-GA-A4,m84094_260226_220949_s1,/data/run2.bam,
```

### 2. Run the Pipeline on your HPC (SLURM + Apptainer)
Launch the pipeline with your HPC profile:
```bash
nextflow run main.nf \
    -profile slurm \
    --samplesheet samplesheet.csv \
    --ref_fasta reference/GRCh38.fasta \
    --ref_gtf reference/gencode.v39.gtf \
    --cage_peaks reference/refTSS.bed \
    --polya_list reference/polyA.list.txt \
    --container_sqanti3 path/to/sqanti3_latest.sif
```

---

## 📂 Output Directory Structure

The final results are highly organized and structured for immediate downstream analysis in R/Seurat:

```
results/
  ├── NK_cohort/
  │   ├── joint/
  │   │   ├── shared_isoform_catalog.tsv.gz    ← Cleaned SQANTI3 Class
  │   │   └── shared_isoform_map.tsv.gz        ← Unified Gene-Transcript Map
  │   │
  │   ├── libA/
  │   │   ├── qc_export/                       ← Ready for R Seurat import
  │   │   │   ├── gene/
  │   │   │   │   ├── matrix.mtx.gz
  │   │   │   │   ├── barcodes.tsv.gz
  │   │   │   │   └── features.tsv.gz
  │   │   │   └── transcript/
  │   │   │       ├── matrix.mtx.gz
  │   │   │       ├── barcodes.tsv.gz
  │   │   │       └── features.tsv.gz          ← Enriched with SQANTI3 categories
  │   │   │
  │   │   └── qc/
  │   │       └── libA_saturation.tsv          ← Rarefaction statistics
  │   │
  │   └── [other libraries]
  └── multiqc/                                 ← Now includes saturation graphs
```
