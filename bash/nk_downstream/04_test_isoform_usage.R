#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "adj_p_value",
    "case_gene_count",
    "case_group",
    "case_isoform_count",
    "case_prop",
    "category",
    "consistent_direction",
    "delta_asin",
    "delta_prop",
    "donor",
    "donors_negative",
    "donors_positive",
    "gene_count",
    "gene_id",
    "isoform_count",
    "isoform_prop",
    "mean_case_prop",
    "mean_delta_asin",
    "mean_delta_prop",
    "mean_ref_prop",
    "n_consistent_isoforms",
    "n_donors",
    "n_isoforms_tested",
    "p_value",
    "population",
    "ref_gene_count",
    "ref_group",
    "ref_isoform_count",
    "ref_prop",
    "sample_id",
    "sd_delta_asin",
    "shared_isoform_id"
  ))
}

# =============================================================================
# Test donor-matched isoform usage shifts within genes.
# =============================================================================

BASE_DIR <- "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs"
PROJECT_ID <- "msk"

WORKING_DIR <- file.path(BASE_DIR, PROJECT_ID)
COUNT_PATH <- file.path(WORKING_DIR, "downstream", "03_population_pseudobulks", "isoform_pseudobulk_counts.tsv.gz")
SAMPLE_PATH <- file.path(WORKING_DIR, "downstream", "03_population_pseudobulks", "sample_metadata.tsv")
OUTDIR <- file.path(WORKING_DIR, "downstream", "04_isoform_usage")

COMPARISONS <- list(
  c("adaptive_nk", "cd56dim_nk"),
  c("cd56bright_nk", "cd56dim_nk"),
  c("cycling_nk", "cd56dim_nk")
)

MIN_GENE_COUNT <- 20
MIN_ISOFORM_COUNT <- 5

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

paired_test <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(NA_real_)
  }
  tryCatch(t.test(x)$p.value, error = function(...) NA_real_)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  min(x)
}

counts_dt <- fread(COUNT_PATH)
sample_dt <- fread(SAMPLE_PATH)

sample_cols <- intersect(sample_dt$sample_id, names(counts_dt))
if (length(sample_cols) == 0) {
  stop("No pseudobulk sample columns were found in the isoform count table")
}

long_dt <- melt(
  counts_dt,
  id.vars = c("shared_isoform_id", "gene_id", "category"),
  measure.vars = sample_cols,
  variable.name = "sample_id",
  value.name = "isoform_count"
)

long_dt <- merge(long_dt, sample_dt, by = "sample_id", all.x = TRUE)
long_dt[, isoform_count := as.numeric(isoform_count)]

gene_totals <- long_dt[, .(gene_count = sum(isoform_count)), by = .(sample_id, gene_id)]
long_dt <- merge(long_dt, gene_totals, by = c("sample_id", "gene_id"), all.x = TRUE)
long_dt[, isoform_prop := fifelse(gene_count > 0, isoform_count / gene_count, 0)]

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

  subset_dt <- long_dt[donor %in% balanced_donors & population %in% c(case_group, ref_group)]

  donor_level <- subset_dt[, .(
    case_isoform_count = isoform_count[population == case_group][1],
    ref_isoform_count = isoform_count[population == ref_group][1],
    case_gene_count = gene_count[population == case_group][1],
    ref_gene_count = gene_count[population == ref_group][1],
    case_prop = isoform_prop[population == case_group][1],
    ref_prop = isoform_prop[population == ref_group][1]
  ), by = .(gene_id, shared_isoform_id, category, donor)]

  donor_level <- donor_level[
    !is.na(case_isoform_count) & !is.na(ref_isoform_count) &
      !is.na(case_gene_count) & !is.na(ref_gene_count)
  ]

  donor_level[, delta_prop := case_prop - ref_prop]
  donor_level[, delta_asin := asin(sqrt((case_isoform_count + 0.5) / (case_gene_count + 1))) -
    asin(sqrt((ref_isoform_count + 0.5) / (ref_gene_count + 1)))]

  donor_level <- donor_level[
    pmax(case_gene_count, ref_gene_count) >= MIN_GENE_COUNT &
      pmax(case_isoform_count, ref_isoform_count) >= MIN_ISOFORM_COUNT
  ]

  if (nrow(donor_level) == 0) {
    next
  }

  feature_results <- donor_level[, .(
    n_donors = .N,
    mean_case_prop = mean(case_prop),
    mean_ref_prop = mean(ref_prop),
    mean_delta_prop = mean(delta_prop),
    mean_delta_asin = mean(delta_asin),
    sd_delta_asin = sd(delta_asin),
    p_value = paired_test(delta_asin),
    donors_positive = sum(delta_prop > 0),
    donors_negative = sum(delta_prop < 0)
  ), by = .(gene_id, shared_isoform_id, category)]

  feature_results[, adj_p_value := p.adjust(p_value, method = "BH")]
  feature_results[, consistent_direction := donors_positive == n_donors | donors_negative == n_donors]
  feature_results <- feature_results[order(adj_p_value, -abs(mean_delta_prop))]

  gene_results <- feature_results[, .(
    n_isoforms_tested = .N,
    min_feature_p_value = safe_min(p_value),
    min_feature_adj_p_value = safe_min(adj_p_value),
    max_abs_delta_prop = max(abs(mean_delta_prop), na.rm = TRUE),
    n_consistent_isoforms = sum(consistent_direction, na.rm = TRUE)
  ), by = gene_id]

  fwrite(feature_results, file.path(OUTDIR, paste0(comparison_name, "_isoform_results.tsv.gz")), sep = "\t")
  fwrite(gene_results, file.path(OUTDIR, paste0(comparison_name, "_gene_summary.tsv.gz")), sep = "\t")
}