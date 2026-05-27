# longshot

A Nextflow DSL2 pipeline for PacBio HiFi MAS-seq single-cell long-read transcriptomics, taking raw multiplexed SMRT cell BAMs and producing per-library count matrices with a cohort-wide isoform catalog.

---

## What it does

**Preprocessing** — detects and validates 10x sample indices against a downsampled read set before committing to full demultiplexing. Runs SKERA splitting once per SMRT cell, builds chemistry-aware primer FASTAs, demultiplexes with lima, and merges constituent BAMs per library.

**Barcode correction** — if paired Illumina barcodes are provided, strips gem-group suffixes, reverse-complements sequences for PacBio orientation, and uses the resulting empirical whitelist for `isoseq correct`. Falls back to the static 10x whitelist otherwise.

**Alignment and deduplication** — sorts by cell barcode, deduplicates with `isoseq groupdedup`, and aligns to GRCh38 with `pbmm2`.

**Transcript discovery and QC** — performs joint isoform discovery across all libraries with IsoQuant, runs SQANTI3 QC and rules-based filtering, then rescues valid isoforms discarded by the filter (automatic mode).

**Quantification** — re-runs IsoQuant per library against the rescued transcript model to produce per-cell counts.

**Export** — writes gene-level and transcript-level matrices in 10x format (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`). Transcript features are annotated with SQANTI3 structural categories for direct use with `Read10X()`. Generates a cohort-wide shared isoform catalog and saturation curves at multiple depth fractions.

---

## Samplesheet

```csv
experiment,library_id,run_id,bam,10x_index,chemistry,shortread_barcodes
nk_cohort,libA,m84094_260226_200655_s1,/data/run1.bam,SI-GA-A1,3prime,/data/illumina/libA/barcodes.tsv.gz
nk_cohort,libB,m84094_260226_200655_s1,/data/run1.bam,SI-GA-A2,3prime,/data/illumina/libB/barcodes.tsv.gz
nk_cohort,libC,m84094_260226_220949_s1,/data/run2.bam,SI-GA-A3,,NULL
```

`10x_index`, `chemistry`, and `shortread_barcodes` are optional. Set to empty or `NULL` to use defaults (`chemistry` falls back to `params.chemistry`, shortread correction is skipped).

---

## Running

Dry run (checks channel wiring without submitting jobs):

```bash
nextflow run main.nf -preview -profile standard --samplesheet assets/example_samplesheet.csv
```

Production (SLURM + Apptainer):

```bash
nextflow run main.nf \
    -profile slurm \
    --samplesheet samplesheet.csv
```

Reference files are hardcoded in the `slurm` profile in `nextflow.config`. Override any at the command line with `--ref_fasta`, `--ref_gtf`, `--cage_peaks`, `--polya_list`, etc.

---

## Output structure

```
results/
  {experiment}/
    joint/
      transcript_model/   — IsoQuant discovery GTF
      sqanti3/            — QC reports, filter results, rescued GTF and FASTA
    {library_id}/
      counts/             — IsoQuant per-library output
      qc_export/
        gene/             — matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz
        transcript/       — same, features annotated with SQANTI3 categories
      qc/
        mosdepth/
        nanostat/
        flagstat/
  multiqc/                — aggregated QC including saturation curves
```

---

## Pipeline structure

```
main.nf
  subworkflows/preprocess.nf   — SKERA, index detection, lima, merge, tag, refine
  subworkflows/align.nf        — whitelist prep, isoseq correct, dedup, pbmm2, CRAM
  subworkflows/classify.nf     — CB suffix injection, IsoQuant discovery, SQANTI3 QC/filter/rescue
  subworkflows/quantify.nf     — per-library IsoQuant quantification
  subworkflows/export.nf       — MTX export, shared catalog, saturation curves
  modules/qc.nf                — flagstat, mosdepth, NanoStat, MultiQC
```