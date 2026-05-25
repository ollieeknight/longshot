#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".",
    "adj_p_value",
    "case_count",
    "case_gene_count",
    "case_group",
    "case_prop",
    "consistent_direction",
    "count",
    "delta_asin",
    "delta_prop",
    "donor",
    "donors_negative",
    "donors_positive",
    "feature_id",
    "gene_count",
    "gene_for_testing",
    "gene_id",
    "mean_case_prop",
    "mean_delta_prop",
    "mean_ref_prop",
    "n_donors",
    "p_value",
    "polya_feature_id",
    "polya_pos",
    "population",
    "position",
    "prop",
    "ref_count",
    "ref_gene_count",
    "ref_group",
    "ref_prop",
    "sample_id",
    "shared_gene",
    "shared_isoform_id",
    "tss_feature_id",
    "tss_pos"
  ))
}

# =============================================================================
# Test donor-matched TSS and polyA usage shifts within genes.
# =============================================================================

BASE_DIR <- "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs"
PROJECT_ID <- "msk"

WORKING_DIR <- file.path(BASE_DIR, PROJECT_ID)
COUNT_PATH <- file.path(WORKING_DIR, "downstream", "03_population_pseudobulks", "isoform_pseudobulk_counts.tsv.gz")
SAMPLE_PATH <- file.path(WORKING_DIR, "downstream", "03_population_pseudobulks", "sample_metadata.tsv")
CATALOG_PATH <- file.path(WORKING_DIR, "downstream", "02_shared_isoform_catalog", "shared_isoform_catalog.tsv.gz")
OUTDIR <- file.path(WORKING_DIR, "downstream", "06_tss_polya_usage")

COMPARISONS <- list(
  c("adaptive_nk", "cd56dim_nk"),
  c("cd56bright_nk", "cd56dim_nk"),
  c("cycling_nk", "cd56dim_nk")
)

MIN_SITE_COUNT <- 5
MIN_GENE_COUNT <- 20

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

paired_test <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(NA_real_)
  }
  tryCatch(t.test(x)$p.value, error = function(...) NA_real_)
}

aggregate_site_counts <- function(counts_dt, site_map, sample_cols) {
  long_counts <- melt(
    counts_dt,
    id.vars = c("shared_isoform_id", "gene_id", "category"),
    measure.vars = sample_cols,
    variable.name = "sample_id",
    value.name = "count"
  )

  long_counts <- merge(long_counts, site_map, by = "shared_isoform_id", all.x = FALSE)
  long_counts[, .(count = sum(count)), by = .(gene_id = gene_for_testing, feature_id, sample_id)]
}

run_site_test <- function(site_counts, sample_dt, case_group, ref_group, out_prefix) {
  keep_samples <- sample_dt[population %in% c(case_group, ref_group)]
  donor_counts <- keep_samples[, .N, by = .(donor, population)]
  balanced_donors <- donor_counts[, .N, by = donor][N == 2, donor]

  if (length(balanced_donors) < 2) {
    return(invisible(NULL))
  }

  subset_counts <- site_counts[sample_id %in% keep_samples[donor %in% balanced_donors, sample_id]]
  subset_counts <- merge(subset_counts, sample_dt, by = "sample_id", all.x = TRUE)

  gene_totals <- subset_counts[, .(gene_count = sum(count)), by = .(sample_id, gene_id)]
  subset_counts <- merge(subset_counts, gene_totals, by = c("sample_id", "gene_id"), all.x = TRUE)
  subset_counts[, prop := fifelse(gene_count > 0, count / gene_count, 0)]

  donor_level <- subset_counts[, .(
    case_count = count[population == case_group][1],
    ref_count = count[population == ref_group][1],
    case_prop = prop[population == case_group][1],
    ref_prop = prop[population == ref_group][1],
    case_gene_count = gene_count[population == case_group][1],
    ref_gene_count = gene_count[population == ref_group][1]
  ), by = .(gene_id, feature_id, donor)]

  donor_level <- donor_level[
    !is.na(case_count) & !is.na(ref_count) &
      pmax(case_count, ref_count) >= MIN_SITE_COUNT &
      pmax(case_gene_count, ref_gene_count) >= MIN_GENE_COUNT
  ]

  if (nrow(donor_level) == 0) {
    return(invisible(NULL))
  }

  donor_level[, delta_prop := case_prop - ref_prop]
  donor_level[, delta_asin := asin(sqrt((case_count + 0.5) / (case_gene_count + 1))) -
    asin(sqrt((ref_count + 0.5) / (ref_gene_count + 1)))]

  result_dt <- donor_level[, .(
    n_donors = .N,
    mean_case_prop = mean(case_prop),
    mean_ref_prop = mean(ref_prop),
    mean_delta_prop = mean(delta_prop),
    p_value = paired_test(delta_asin),
    donors_positive = sum(delta_prop > 0),
    donors_negative = sum(delta_prop < 0)
  ), by = .(gene_id, feature_id)]

  result_dt[, adj_p_value := p.adjust(p_value, method = "BH")]
  result_dt[, consistent_direction := donors_positive == n_donors | donors_negative == n_donors]
  result_dt <- result_dt[order(adj_p_value, -abs(mean_delta_prop))]

  fwrite(result_dt, paste0(out_prefix, ".tsv.gz"), sep = "\t")
}

counts_dt <- fread(COUNT_PATH)
sample_dt <- fread(SAMPLE_PATH)
catalog_dt <- fread(CATALOG_PATH)

sample_cols <- intersect(sample_dt$sample_id, names(counts_dt))

catalog_dt[, gene_for_testing := shared_gene]
catalog_dt[is.na(gene_for_testing) | gene_for_testing == "", gene_for_testing := paste0("UNASSIGNED_", shared_isoform_id)]

tss_map <- catalog_dt[, .(
  shared_isoform_id,
  gene_for_testing,
  feature_id = tss_feature_id,
  position = tss_pos
)]

polya_map <- catalog_dt[, .(
  shared_isoform_id,
  gene_for_testing,
  feature_id = polya_feature_id,
  position = polya_pos
)]

tss_counts <- aggregate_site_counts(counts_dt, tss_map, sample_cols)
polya_counts <- aggregate_site_counts(counts_dt, polya_map, sample_cols)

fwrite(unique(tss_map), file.path(OUTDIR, "tss_feature_metadata.tsv.gz"), sep = "\t")
fwrite(unique(polya_map), file.path(OUTDIR, "polya_feature_metadata.tsv.gz"), sep = "\t")

for (comparison in COMPARISONS) {
  case_group <- comparison[[1]]
  ref_group <- comparison[[2]]
  comparison_name <- paste(case_group, "vs", ref_group, sep = "_")

  run_site_test(
    tss_counts,
    sample_dt,
    case_group,
    ref_group,
    file.path(OUTDIR, paste0(comparison_name, "_tss_results"))
  )

  run_site_test(
    polya_counts,
    sample_dt,
    case_group,
    ref_group,
    file.path(OUTDIR, paste0(comparison_name, "_polya_results"))
  )
}