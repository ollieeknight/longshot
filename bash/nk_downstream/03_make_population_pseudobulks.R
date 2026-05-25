#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

WORK_DIR <- "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
BARCODE_DIR <- "/data/cephfs-1/work/groups/romagnani/users/knighto_c/data/adaptive_nk/objects/longread/barcodes"
MAP_PATH <- file.path(WORK_DIR, "downstream/02_shared_isoform_catalog/shared_isoform_map.tsv.gz")
CATALOG_PATH <- file.path(WORK_DIR, "downstream/02_shared_isoform_catalog/shared_isoform_catalog.tsv.gz")
OUT_DIR <- file.path(WORK_DIR, "downstream/03_population_pseudobulks")
DONORS <- c("HC01_NK", "HC02_NK", "HC03_NK")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

normalise_barcode <- function(x) sub("-1$", "", trimws(x))
normalise_pbid <- function(x) sub(":.*$", "", trimws(as.character(x)))

pick_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) stop("Missing: ", paste(paths, collapse = ", "))
  hit[[1]]
}

resolve_sample_name <- function(donor_id) {
  if (dir.exists(file.path(WORK_DIR, "outs", donor_id, "isoforms_seurat"))) return(donor_id)
  if (dir.exists(file.path(WORK_DIR, "outs", paste0(donor_id, "_NK"), "isoforms_seurat"))) return(paste0(donor_id, "_NK"))
  stop("Can't find outs directory for ", donor_id)
}

resolve_map_donor <- function(donor_id, sample_name, shared_donors) {
  candidates <- unique(c(donor_id, sample_name, sub("_NK$", "", c(donor_id, sample_name))))
  hit <- intersect(candidates, shared_donors)
  if (length(hit) == 0) stop("Can't match donor in shared map: ", donor_id)
  hit[[1]]
}

resolve_barcode_path <- function(donor_id, sample_name) {
  candidates <- unique(c(
    file.path(BARCODE_DIR, paste0(sub("_NK$", "", donor_id), ".tsv")),
    file.path(BARCODE_DIR, paste0(sub("_NK$", "", sample_name), ".tsv")),
    file.path(BARCODE_DIR, paste0(c(donor_id, sample_name), ".tsv"))
  ))
  pick_existing(candidates)
}

aggregate_to_shared <- function(mat, shared_ids) {
  idx <- match(shared_ids, unique(shared_ids))
  map <- sparseMatrix(
    i = idx,
    j = seq_along(shared_ids),
    x = 1,
    dims = c(length(unique(shared_ids)), length(shared_ids))
  )
  out <- map %*% mat
  rownames(out) <- unique(shared_ids)
  out
}

align_to_ids <- function(block, all_ids) {
  row_idx <- match(rownames(block), all_ids)
  if (any(is.na(row_idx))) stop("Unexpected row id while aligning blocks")
  row_map <- sparseMatrix(
    i = row_idx,
    j = seq_along(row_idx),
    x = 1,
    dims = c(length(all_ids), length(row_idx))
  )
  out <- row_map %*% block
  rownames(out) <- all_ids
  colnames(out) <- colnames(block)
  out
}

write_sparse_singlecell_outputs <- function(sc_mat, iso_annot, out_prefix) {
  # Persist sparse single-cell counts without densifying to a giant wide table.
  writeMM(sc_mat, paste0(out_prefix, ".mtx"))

  iso_dt <- data.table(shared_isoform_id = rownames(sc_mat))
  iso_dt <- merge(iso_dt, iso_annot, by = "shared_isoform_id", all.x = TRUE)
  setcolorder(iso_dt, c("shared_isoform_id", "gene_id", "category"))
  fwrite(iso_dt, paste0(out_prefix, "_rows.tsv.gz"), sep = "\t")

  col_dt <- data.table(cell_id = colnames(sc_mat), col_index = seq_len(ncol(sc_mat)))
  fwrite(col_dt, paste0(out_prefix, "_cols.tsv.gz"), sep = "\t")
}

shared_map <- fread(MAP_PATH)
shared_catalog <- fread(CATALOG_PATH)
shared_donors <- unique(shared_map$donor)
iso_annot <- unique(shared_catalog[, .(shared_isoform_id, gene_id = shared_gene, category = modal_category)])

bulk_blocks <- list()

for (donor_id in DONORS) {
  sample_name <- resolve_sample_name(donor_id)
  map_donor <- resolve_map_donor(donor_id, sample_name, shared_donors)
  iso_dir <- file.path(WORK_DIR, "outs", sample_name, "isoforms_seurat")
  
  feature_table <- fread(pick_existing(c(file.path(iso_dir, "genes.tsv"), file.path(iso_dir, "features.tsv"))), 
                         sep = "\t", header = FALSE)
  barcode_table <- fread(file.path(iso_dir, "barcodes.tsv"), sep = "\t", header = FALSE)
  count_matrix <- readMM(file.path(iso_dir, "matrix.mtx"))
  
  feature_ids <- normalise_pbid(feature_table[[1]])
  cell_barcodes <- normalise_barcode(as.character(barcode_table[[1]]))
  
  pop_labels <- fread(resolve_barcode_path(donor_id, sample_name), sep = "\t", header = FALSE)
  if (ncol(pop_labels) < 2) stop("Barcode file needs ≥2 columns: ", resolve_barcode_path(donor_id, sample_name))
  
  setnames(pop_labels, c("barcode", "population"))
  pop_labels[, barcode := normalise_barcode(barcode)]
  pop_labels <- unique(pop_labels[, .(barcode, population)])
  
  cell_annot <- data.table(barcode = cell_barcodes, col_id = seq_along(cell_barcodes))
  cell_annot <- merge(cell_annot, pop_labels, by = "barcode", all.x = TRUE)
  matched <- cell_annot[!is.na(population)]
  
  if (nrow(matched) == 0) stop("No annotated cells for ", donor_id)
  
  donor_map <- unique(shared_map[donor == map_donor, .(pbid, shared_isoform_id)])
  donor_features <- data.table(pbid = feature_ids, row_id = seq_along(feature_ids))
  donor_features <- merge(donor_features, donor_map, by = "pbid", all.x = TRUE)[!is.na(shared_isoform_id)]
  
  if (nrow(donor_features) == 0) stop("No donor isoforms mapped: ", donor_id)
  
  count_matrix <- count_matrix[donor_features$row_id, , drop = FALSE]
  shared_ids <- donor_features$shared_isoform_id
  
  # Build per-donor population pseudobulks and collapse repeated shared IDs.
  pop_levels <- sort(unique(matched$population))
  bulk_raw <- do.call(cbind, lapply(pop_levels, function(pop) {
    keep <- matched[population == pop, col_id]
    Matrix::rowSums(count_matrix[, keep, drop = FALSE])
  }))
  colnames(bulk_raw) <- pop_levels
  bulk_collapsed <- aggregate_to_shared(bulk_raw, shared_ids)
  colnames(bulk_collapsed) <- paste(donor_id, colnames(bulk_collapsed), sep = "__")
  bulk_blocks[[donor_id]] <- bulk_collapsed

  # Build donor-specific single-cell matrix at shared isoform level.
  sc_cols <- paste(donor_id, matched$barcode, sep = "__")
  sc_raw <- count_matrix[, matched$col_id, drop = FALSE]
  colnames(sc_raw) <- sc_cols
  sc_collapsed <- aggregate_to_shared(sc_raw, shared_ids)

  write_sparse_singlecell_outputs(
    sc_collapsed,
    iso_annot,
    file.path(OUT_DIR, paste0(donor_id, "_singlecell_isoform_counts"))
  )
}

all_iso_ids <- sort(unique(unlist(lapply(bulk_blocks, rownames))))
bulk_matrix <- do.call(cbind, lapply(names(bulk_blocks), function(donor_id) {
  align_to_ids(bulk_blocks[[donor_id]], all_iso_ids)
}))
rownames(bulk_matrix) <- all_iso_ids

bulk_iso <- data.table(shared_isoform_id = all_iso_ids)
bulk_iso <- merge(bulk_iso, iso_annot, by = "shared_isoform_id", all.x = TRUE)

pop_cols <- sort(unique(sub("^.*__", "", colnames(bulk_matrix))))
for (col in pop_cols) {
  col_ids <- grep(paste0("__", col, "$"), colnames(bulk_matrix))
  if (length(col_ids) == 0) {
    bulk_iso[[col]] <- 0
  } else {
    bulk_iso[[col]] <- as.numeric(Matrix::rowSums(bulk_matrix[, col_ids, drop = FALSE]))
  }
}

gene_ids <- bulk_iso$gene_id
gene_ids[is.na(gene_ids) | gene_ids == ""] <- paste0("UNASSIGNED_", bulk_iso$shared_isoform_id[is.na(gene_ids) | gene_ids == ""])
unique_gene_ids <- unique(gene_ids)
gene_mapping <- sparseMatrix(
  i = match(gene_ids, unique_gene_ids), j = seq_along(gene_ids), x = 1,
  dims = c(length(unique_gene_ids), length(gene_ids))
)
rownames(gene_mapping) <- unique_gene_ids

gene_bulk <- gene_mapping %*% as.matrix(bulk_iso[, ..pop_cols])
bulk_gene <- data.table(gene_id = rownames(gene_mapping))
for (col in pop_cols) {
  bulk_gene[[col]] <- as.numeric(gene_bulk[, col])
}

fwrite(bulk_iso, file.path(OUT_DIR, "bulk_isoform_counts.tsv.gz"), sep = "\t")
fwrite(bulk_gene, file.path(OUT_DIR, "bulk_gene_counts.tsv.gz"), sep = "\t")