#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

WORK_DIR <- "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
CAT_DIR <- file.path(WORK_DIR, "downstream", "02_shared_isoform_catalog")
OUT_DIR <- file.path(CAT_DIR, "exploration")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

catalog_path <- file.path(CAT_DIR, "shared_isoform_catalog.tsv.gz")
map_path <- file.path(CAT_DIR, "shared_isoform_map.tsv.gz")

if (!file.exists(catalog_path)) stop("Missing file: ", catalog_path)
if (!file.exists(map_path)) stop("Missing file: ", map_path)

message("Loading shared catalog and map...")
catalog <- fread(catalog_path)
iso_map <- fread(map_path)

required_catalog <- c("shared_isoform_id", "shared_gene", "modal_category", "n_donors", "tss_feature_id", "polya_feature_id")
required_map <- c("donor", "pbid", "shared_isoform_id")

miss_catalog <- setdiff(required_catalog, names(catalog))
miss_map <- setdiff(required_map, names(iso_map))
if (length(miss_catalog) > 0) stop("Catalog is missing columns: ", paste(miss_catalog, collapse = ", "))
if (length(miss_map) > 0) stop("Map is missing columns: ", paste(miss_map, collapse = ", "))

message("1) Building donor-overlap summary...")
donor_overlap <- catalog[, .N, by = n_donors][order(-n_donors)]
fwrite(donor_overlap, file.path(OUT_DIR, "01_donor_overlap_summary.tsv.gz"), sep = "\t")

message("2) Ranking recurrent novel isoforms...")
novel_idx <- grepl("novel", catalog$modal_category, ignore.case = TRUE)
novel_candidates <- copy(catalog[novel_idx])
setorder(novel_candidates, -n_donors, shared_gene, shared_isoform_id)
novel_top <- novel_candidates[, .(
  shared_isoform_id,
  shared_gene,
  modal_category,
  n_donors,
  n_exons,
  donors,
  exon_signature,
  intron_signature,
  tss_feature_id,
  polya_feature_id
)]
fwrite(novel_top, file.path(OUT_DIR, "02_top_recurrent_novel_isoforms.tsv.gz"), sep = "\t")

message("3) Measuring per-gene TSS/polyA diversity...")
gene_site_diversity <- catalog[!is.na(shared_gene) & shared_gene != "", .(
  n_isoforms = uniqueN(shared_isoform_id),
  n_tss = uniqueN(tss_feature_id),
  n_polya = uniqueN(polya_feature_id),
  max_donor_support = suppressWarnings(max(n_donors, na.rm = TRUE))
), by = shared_gene]

gene_site_diversity[!is.finite(max_donor_support), max_donor_support := NA_real_]
setorder(gene_site_diversity, -n_tss, -n_polya, -n_isoforms)
fwrite(gene_site_diversity, file.path(OUT_DIR, "03_gene_tss_polya_diversity.tsv.gz"), sep = "\t")

message("4) Creating a stats-ready isoform table...")
stats_ready <- merge(
  iso_map,
  catalog[, .(shared_isoform_id, shared_gene, modal_category, n_donors, tss_feature_id, polya_feature_id)],
  by = "shared_isoform_id",
  all.x = TRUE
)
setorder(stats_ready, donor, shared_gene, shared_isoform_id, pbid)
fwrite(stats_ready, file.path(OUT_DIR, "04_stats_ready_isoform_map.tsv.gz"), sep = "\t")

message("Done.")
message("Wrote files in: ", OUT_DIR)
message(" - 01_donor_overlap_summary.tsv.gz")
message(" - 02_top_recurrent_novel_isoforms.tsv.gz")
message(" - 03_gene_tss_polya_diversity.tsv.gz")
message(" - 04_stats_ready_isoform_map.tsv.gz")
