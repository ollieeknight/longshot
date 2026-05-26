# Reference Files

All reference files required to run the `longshot` pipeline. Download once and point `nextflow.config` or your run script at the resulting paths.

---

## 1. Genome FASTA — GRCh38 no-alt analysis set

The no-alt analysis set excludes alternative contigs and decoy sequences, which prevents multi-mapping artefacts in pbmm2 and downstream tools.

**File to download:**
```
GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

**Source:** NCBI — GRCh38.p14 alignment pipelines directory
```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/
```

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

**Pipeline param:** `--ref_fasta`

---

## 2. Gene Annotation GTF — GENCODE v49

GENCODE v49 is the current release for GRCh38. It uses `chr`-style chromosome names compatible with pbmm2, IsoQuant, and SQANTI3.

**Source:** https://www.gencodegenes.org/human/

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
gunzip gencode.v49.annotation.gtf.gz

# Sort by chromosome and position (required by IsoQuant and SQANTI3)
grep "^#" gencode.v49.annotation.gtf > gencode.v49.annotation.sorted.gtf
grep -v "^#" gencode.v49.annotation.gtf \
    | sort -k1,1 -k4,4n -k5,5n \
    >> gencode.v49.annotation.sorted.gtf

# Compress and index (for tabix-based tools)
bgzip -c gencode.v49.annotation.sorted.gtf > gencode.v49.annotation.sorted.gtf.gz
tabix -p gff gencode.v49.annotation.sorted.gtf.gz
```

**Pipeline param:** `--ref_gtf` — point at `gencode.v49.annotation.sorted.gtf`

> **Note:** Do not use the RefSeq GTF from the NCBI FASTA directory — it uses `NC_000001.11`-style chromosome names that are incompatible with the pipeline.

---

## 3. CAGE Peaks — refTSS v4.1 (hg38)

Used by SQANTI3 to annotate transcription start sites.

**Source:** SQANTI3 official reference dataset (Figshare)
```
https://figshare.com/articles/dataset/SQANTI3_reference_annotations_-_human/28537373
```

```bash
wget -O refTSS_v4.1_human_coordinate.hg38.bed.txt.gz \
    https://ndownloader.figshare.com/files/52801133
gunzip refTSS_v4.1_human_coordinate.hg38.bed.txt.gz

# Sort
sort -k1,1 -k2,2n refTSS_v4.1_human_coordinate.hg38.bed.txt \
    > refTSS_v4.1_human_coordinate.hg38.sorted.bed
```

**Pipeline param:** `--cage_peaks`

---

## 4. PolyA Motif List — SQANTI3

A short plain-text file listing polyA signal motifs.

**Source:** SQANTI3 official reference dataset (Figshare)
```
https://figshare.com/articles/dataset/SQANTI3_reference_annotations_-_human/28537373
```

```bash
wget -O mouse_and_human.polyA_motif.txt \
    https://ndownloader.figshare.com/files/52801139
```

**Pipeline param:** `--polya_list`

---

## 5. PolyA Site Atlas — atlas.clusters (optional)

Provides polyA cleavage site coordinates for SQANTI3 `--polyA_peak` annotation. Not currently wired as a pipeline parameter but useful for manual SQANTI3 runs.

**Source:** SQANTI3 official reference dataset (Figshare)

```bash
wget -O atlas.clusters.2.0.GRCh38.96_chr.bed.gz \
    https://ndownloader.figshare.com/files/52801130
gunzip atlas.clusters.2.0.GRCh38.96_chr.bed.gz
```

---

## 6. Junction Coverage File — Intropolis (optional)

Improves SQANTI3 splice junction annotation. This parameter is **optional** (`--intropolis null` by default).

**Source:** SQANTI3 official reference dataset (Figshare)

```bash
wget -O intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz \
    https://ndownloader.figshare.com/files/52801127
gunzip intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz

# Sort (required by SQANTI3)
sort -k1,1 -k2,2n intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified \
    > intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.sorted.tsv
```

**Pipeline param:** `--intropolis` (optional)

---

## 6. SQANTI3 Container (SIF)

SQANTI3 is not available as a standard Apptainer pull from BioContainers in all environments. Build or pull the SIF once:

```bash
apptainer pull sqanti3_6.0.1.sif docker://quay.io/biocontainers/sqanti3:6.0.1--hdfd78af_0
```

**Pipeline param:** `--container_sqanti3 /path/to/sqanti3_6.0.1.sif`

---

## Version Summary

| File | Version | Source |
|---|---|---|
| Genome FASTA | GRCh38 / GCA_000001405.15 | NCBI |
| Gene annotation | GENCODE v49 | gencodegenes.org |
| CAGE peaks | refTSS v4.1 | figshare/28537373 |
| PolyA motifs | mouse_and_human.polyA_motif | figshare/28537373 |
| PolyA site atlas | atlas.clusters 2.0 GRCh38 (optional) | figshare/28537373 |
| Junction coverage | intropolis v1 modified (optional) | figshare/28537373 |
| SQANTI3 container | 6.0.1 | quay.io/biocontainers |
