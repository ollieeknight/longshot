#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "case_count",
    "case_group",
    "case_prop",
    "category",
    "category_lower",
    "delta_prop",
    "detected",
    "detected_case",
    "detected_ref",
    "donor",
    "donors_detected_case",
    "donors_detected_ref",
    "donors_negative",
    "donors_positive",
    "gene_count",
    "gene_id",
    "isoform_count",
    "isoform_prop",
    "max_prop",
    "mean_case_prop",
    "mean_delta_prop",
    "mean_ref_prop",
    "n_donors",
    "population",
    "populations_detected",
    "ref_count",
    "ref_group",
    "ref_prop",
    "sample_id",
    "seqname",
    "shared_isoform_id",
    "strand",
    "total_counts",
    "tx_end",
    "tx_start"
  ))
}

# =============================================================================
# Rank recurrent novel isoforms across donors and populations.
# =============================================================================

BASE_DIR <- "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs"
PROJECT_ID <- "msk"

WORKING_DIR <- file.path(BASE_DIR, PROJECT_ID)
COUNT_PATH <- file.path(WORKING_DIR, "downstream", "03_population_pseudobulks", "isoform_pseudobulk_counts.tsv.gz")
SAMPLE_PATH <- file.path(WORKING_DIR, "downstream", "03_population_pseudobulks", "sample_metadata.tsv")
CATALOG_PATH <- file.path(WORKING_DIR, "downstream", "02_shared_isoform_catalog", "shared_isoform_catalog.tsv.gz")
OUTDIR <- file.path(WORKING_DIR, "downstream", "05_novel_isoforms")

COMPARISONS <- list(
  c("adaptive_nk", "cd56dim_nk"),
  c("cd56bright_nk", "cd56dim_nk"),
  c("cycling_nk", "cd56dim_nk")
)

NOVEL_CATEGORIES <- c(
  "novel_in_catalog",
  "novel_not_in_catalog",
  "incomplete-splice_match",
  "genic",
  "antisense"
)

MIN_COUNT <- 3

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

counts_dt <- fread(COUNT_PATH)
sample_dt <- fread(SAMPLE_PATH)
catalog_dt <- fread(CATALOG_PATH)

sample_cols <- intersect(sample_dt$sample_id, names(counts_dt))

long_dt <- melt(
  counts_dt,
  id.vars = c("shared_isoform_id", "gene_id", "category"),
  measure.vars = sample_cols,
  variable.name = "sample_id",
  value.name = "isoform_count"
)

long_dt <- merge(long_dt, sample_dt, by = "sample_id", all.x = TRUE)
long_dt <- merge(
  long_dt,
  catalog_dt[, .(shared_isoform_id, seqname, strand, tx_start, tx_end, n_donors)],
  by = "shared_isoform_id",
  all.x = TRUE
)

long_dt[, category_lower := tolower(category)]
long_dt <- long_dt[category_lower %in% NOVEL_CATEGORIES]

gene_totals <- long_dt[, .(gene_count = sum(isoform_count)), by = .(sample_id, gene_id)]
long_dt <- merge(long_dt, gene_totals, by = c("sample_id", "gene_id"), all.x = TRUE)
long_dt[, isoform_prop := fifelse(gene_count > 0, isoform_count / gene_count, 0)]
long_dt[, detected := isoform_count >= MIN_COUNT]

overall_summary <- long_dt[, .(
  category = category[[1]],
  seqname = seqname[[1]],
  strand = strand[[1]],
  tx_start = tx_start[[1]],
  tx_end = tx_end[[1]],
  donors_detected = uniqueN(donor[detected]),
  populations_detected = uniqueN(population[detected]),
  total_counts = sum(isoform_count),
  max_prop = max(isoform_prop)
), by = .(gene_id, shared_isoform_id)]

setorder(overall_summary, -donors_detected, -total_counts, -max_prop)
fwrite(overall_summary, file.path(OUTDIR, "recurrent_novel_isoforms.tsv.gz"), sep = "\t")

for (comparison in COMPARISONS) {
  case_group <- comparison[[1]]
  ref_group <- comparison[[2]]
  comparison_name <- paste(case_group, "vs", ref_group, sep = "_")

  keep_samples <- sample_dt[population %in% c(case_group, ref_group)]
  donor_counts <- keep_samples[, .N, by = .(donor, population)]
  balanced_donors <- donor_counts[, .N, by = donor][N == 2, donor]

  if (length(balanced_donors) < 2) {
    next
  }

  donor_level <- long_dt[donor %in% balanced_donors & population %in% c(case_group, ref_group), .(
    case_count = isoform_count[population == case_group][1],
    ref_count = isoform_count[population == ref_group][1],
    case_prop = isoform_prop[population == case_group][1],
    ref_prop = isoform_prop[population == ref_group][1]
  ), by = .(gene_id, shared_isoform_id, category, donor, seqname, strand, tx_start, tx_end)]

  donor_level <- donor_level[!is.na(case_count) & !is.na(ref_count)]
  donor_level[, delta_prop := case_prop - ref_prop]
  donor_level[, detected_case := case_count >= MIN_COUNT]
  donor_level[, detected_ref := ref_count >= MIN_COUNT]

  if (nrow(donor_level) == 0) {
    next
  }

  ranked <- donor_level[, .(
    n_donors = .N,
    donors_detected_case = sum(detected_case),
    donors_detected_ref = sum(detected_ref),
    mean_case_prop = mean(case_prop),
    mean_ref_prop = mean(ref_prop),
    mean_delta_prop = mean(delta_prop),
    donors_positive = sum(delta_prop > 0),
    donors_negative = sum(delta_prop < 0),
    category = category[[1]],
    seqname = seqname[[1]],
    strand = strand[[1]],
    tx_start = tx_start[[1]],
    tx_end = tx_end[[1]]
  ), by = .(gene_id, shared_isoform_id)]

  ranked[, consistent_direction := donors_positive == n_donors | donors_negative == n_donors]
  ranked <- ranked[order(-consistent_direction, -abs(mean_delta_prop), -donors_detected_case, donors_detected_ref)]

  fwrite(ranked, file.path(OUTDIR, paste0(comparison_name, "_novel_isoform_candidates.tsv.gz")), sep = "\t")
}