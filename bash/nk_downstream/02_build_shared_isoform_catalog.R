#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

invisible(utils::globalVariables(c(
  ".", "attributes", "donor", "end", "exact_structure_key", "exon_end", "exon_signature",
  "exon_start", "feature", "gene_for_sites", "gene_id", "gene_label", "gff_gene_id",
  "intron_signature", "junction_chain_key", "modal_category", "n_exons", "pbid",
  "pigeon_category", "pigeon_gene", "pigeon_subcategory", "pigeon_transcript",
  "polya_cluster_index", "polya_feature_id", "polya_pos", "seqname", "shared_gene",
  "shared_isoform_id", "shared_junction_id", "shared_transcript", "start", "strand",
  "transcript_id", "tss_cluster_index", "tss_feature_id", "tss_pos", "tx_end", "tx_start"
)))

BASE_DIR <- "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs"
WORK_DIR <- file.path(BASE_DIR, "msk")
BARCODE_DIR <- "/data/cephfs-1/work/groups/romagnani/users/knighto_c/data/adaptive_nk/objects/longread/barcodes"
OUT_DIR <- file.path(WORK_DIR, "downstream", "02_shared_isoform_catalog")
CP_DIR <- file.path(OUT_DIR, "checkpoints")
DONORS <- c("HC01_NK", "HC02_NK", "HC03_NK")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CP_DIR, recursive = TRUE, showWarnings = FALSE)

HAS_QS <- requireNamespace("qs", quietly = TRUE)

pick_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) stop("Missing: ", paste(paths, collapse = ", "))
  hit[[1]]
}

extract_attribute <- function(x, key) {
  patterns <- c(
    paste0("(?:^|;)\\s*", key, "=([^;]+)"),
    paste0("(?:^|;)\\s*", key, " \\\"([^\\\"]+)\\\""),
    paste0("(?:^|;)\\s*", key, " '([^']+)'")
  )
  value <- rep(NA_character_, length(x))
  for (pattern in patterns) {
    keep <- is.na(value)
    if (!any(keep)) break
    match <- regexec(pattern, x[keep], perl = TRUE)
    parsed <- regmatches(x[keep], match)
    value[keep] <- vapply(parsed, function(z) if (length(z) >= 2) z[[2]] else NA_character_, character(1))
  }
  trimws(value)
}

pick_column <- function(dt, choices, label) {
  hit <- intersect(choices, names(dt))
  if (length(hit) == 0) stop("Can't find ", label, " in: ", paste(names(dt), collapse = ", "))
  hit[[1]]
}

cluster_positions <- function(x, tolerance = 25L) {
  x <- as.integer(x)
  order_idx <- order(x, na.last = TRUE)
  out <- integer(length(x))
  current_cluster <- 0L
  last_pos <- NA_integer_
  
  for (idx in order_idx) {
    if (is.na(x[[idx]])) {
      out[[idx]] <- NA_integer_
      next
    }
    if (is.na(last_pos) || abs(x[[idx]] - last_pos) > tolerance) {
      current_cluster <- current_cluster + 1L
    }
    out[[idx]] <- current_cluster
    last_pos <- x[[idx]]
  }
  out
}

safe_median_int <- function(x) {
  x <- as.integer(x[!is.na(x)])
  if (length(x) == 0) NA_integer_ else as.integer(round(stats::median(x)))
}

read_or_build_checkpoint <- function(path, build_fun, label) {
  if (HAS_QS && file.exists(path)) {
    message("Loading: ", label)
    return(qs::qread(path))
  }
  value <- build_fun()
  if (HAS_QS) {
    message("Saving: ", label)
    qs::qsave(value, path, preset = "high")
  }
  value
}

read_donor_gff <- function(donor) {
  gff_path <- file.path(WORK_DIR, "workspace", donor, "collapsed.sorted.gff")
  if (!file.exists(gff_path)) stop("Missing GFF: ", gff_path)
  
  gff <- fread(gff_path, sep = "\t", header = FALSE, quote = "", comment.char = "#", fill = TRUE, showProgress = TRUE)
  setnames(gff, c("seqname", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"))
  
  gff[, transcript_id := extract_attribute(attributes, "transcript_id")]
  gff[is.na(transcript_id), transcript_id := extract_attribute(attributes, "ID")]
  gff[, gene_id := extract_attribute(attributes, "gene_id")]
  gff[is.na(gene_id), gene_id := extract_attribute(attributes, "gene")]
  
  transcripts <- gff[feature == "transcript", .(
    donor, pbid = transcript_id, seqname, strand,
    tx_start = as.integer(start), tx_end = as.integer(end), gff_gene_id = gene_id
  )]
  
  exons <- gff[feature == "exon", .(
    donor, pbid = transcript_id, exon_start = as.integer(start), exon_end = as.integer(end)
  )][!is.na(pbid)]
  setorder(exons, donor, pbid, exon_start, exon_end)
  
  exon_summary <- exons[, .(
    n_exons = .N,
    exon_signature = paste(paste0(exon_start, "-", exon_end), collapse = "|"),
    intron_signature = if (.N > 1) paste(paste0(exon_end[-.N], "-", exon_start[-1]), collapse = "|") else "monoexonic"
  ), by = .(donor, pbid)]
  
  merged <- merge(transcripts, exon_summary, by = c("donor", "pbid"), all.x = TRUE)
  merged[, exact_structure_key := paste(seqname, strand, exon_signature, sep = "|")]
  merged[, junction_chain_key := paste(seqname, strand, intron_signature, sep = "|")]
  merged[, tss_pos := fifelse(strand == "+", tx_start, tx_end)]
  merged[, polya_pos := fifelse(strand == "+", tx_end, tx_start)]
  merged
}

read_donor_classification <- function(donor) {
  class_path <- pick_existing_file(c(
    file.path(WORK_DIR, "workspace", donor, "collapsed_classification.filtered_lite_classification.txt"),
    file.path(WORK_DIR, "workspace", donor, "collapsed_classification.filtered_classification.txt"),
    file.path(WORK_DIR, "workspace", donor, "collapsed_classification.txt")
  ))
  
  dt <- fread(class_path, sep = "\t", header = TRUE, showProgress = FALSE)
  
  isoform_col <- pick_column(dt, c("isoform", "pbid", "id"), "isoform")
  gene_col <- pick_column(dt, c("associated_gene", "gene", "gene_name"), "gene")
  transcript_col <- pick_column(dt, c("associated_transcript", "transcript", "transcript_id"), "transcript")
  category_col <- pick_column(dt, c("structural_category", "category"), "category")
  
  out <- dt[, .(
    donor, pbid = as.character(get(isoform_col)), pigeon_gene = as.character(get(gene_col)),
    pigeon_transcript = as.character(get(transcript_col)), pigeon_category = as.character(get(category_col))
  )]
  
  if ("subcategory" %in% names(dt)) {
    out[, pigeon_subcategory := as.character(dt[["subcategory"]])]
  } else {
    out[, pigeon_subcategory := NA_character_]
  }
  out
}

structure_list <- read_or_build_checkpoint(
  file.path(CP_DIR, "structure_list.qs"),
  function() rbindlist(lapply(DONORS, read_donor_gff), fill = TRUE),
  "structure_list"
)

classification_list <- read_or_build_checkpoint(
  file.path(CP_DIR, "classification_list.qs"),
  function() rbindlist(lapply(DONORS, read_donor_classification), fill = TRUE),
  "classification_list"
)

isoform_map <- read_or_build_checkpoint(
  file.path(CP_DIR, "isoform_map_initial.qs"),
  function() merge(structure_list, classification_list, by = c("donor", "pbid"), all.x = TRUE),
  "isoform_map_initial"
)

isoform_map[, gene_label := pigeon_gene]
isoform_map[is.na(gene_label) | gene_label == "", gene_label := gff_gene_id]

exact_levels <- unique(isoform_map$exact_structure_key)
junction_levels <- unique(isoform_map$junction_chain_key)

isoform_map[, shared_isoform_id := sprintf("SHARED.%06d", match(exact_structure_key, exact_levels))]
isoform_map[, shared_junction_id := sprintf("JUNC.%06d", match(junction_chain_key, junction_levels))]

shared_catalog <- read_or_build_checkpoint(
  file.path(CP_DIR, "shared_catalog_pre_sites.qs"),
  function() {
    isoform_map[, {
      gene_votes <- na.omit(gene_label)
      category_votes <- na.omit(pigeon_category)
      transcript_votes <- na.omit(pigeon_transcript)
      
      .(seqname = seqname[[1]], strand = strand[[1]], tx_start = as.integer(min(tx_start, na.rm = TRUE)),
        tx_end = as.integer(max(tx_end, na.rm = TRUE)), tss_pos = safe_median_int(tss_pos),
        polya_pos = safe_median_int(polya_pos), n_exons = as.integer(n_exons[[1]]),
        exon_signature = exon_signature[[1]], intron_signature = intron_signature[[1]],
        shared_junction_id = shared_junction_id[[1]],
        shared_gene = if (length(gene_votes) > 0) names(sort(table(gene_votes), decreasing = TRUE))[1] else NA_character_,
        shared_transcript = if (length(transcript_votes) > 0) names(sort(table(transcript_votes), decreasing = TRUE))[1] else NA_character_,
        modal_category = if (length(category_votes) > 0) names(sort(table(category_votes), decreasing = TRUE))[1] else NA_character_,
        n_donors = as.integer(uniqueN(donor)), donors = paste(sort(unique(donor)), collapse = ","),
        n_gene_labels = as.integer(uniqueN(na.omit(gene_label))))
    }, by = shared_isoform_id]
  },
  "shared_catalog_pre_sites"
)

shared_catalog[, gene_for_sites := shared_gene]
shared_catalog[is.na(gene_for_sites) | gene_for_sites == "", 
               gene_for_sites := paste0("UNASSIGNED_", seqname, "_", strand)]

shared_catalog[, tss_cluster_index := cluster_positions(tss_pos, tolerance = 25L), 
               by = .(gene_for_sites, seqname, strand)]
shared_catalog[, polya_cluster_index := cluster_positions(polya_pos, tolerance = 25L), 
               by = .(gene_for_sites, seqname, strand)]
shared_catalog[, tss_feature_id := sprintf("%s_TSS_%03d", gene_for_sites, tss_cluster_index)]
shared_catalog[, polya_feature_id := sprintf("%s_POLYA_%03d", gene_for_sites, polya_cluster_index)]

isoform_map <- read_or_build_checkpoint(
  file.path(CP_DIR, "isoform_map_final.qs"),
  function() merge(isoform_map, shared_catalog[, .(shared_isoform_id, shared_gene, shared_transcript, modal_category, tss_feature_id, polya_feature_id)], 
                   by = "shared_isoform_id", all.x = TRUE),
  "isoform_map_final"
)

shared_catalog <- shared_catalog[!grepl("^(RPL|RPS|MRPL|MRPS|MALAT1)", shared_gene, ignore.case = TRUE)]
isoform_map <- isoform_map[shared_isoform_id %in% shared_catalog$shared_isoform_id]

fwrite(shared_catalog, file.path(OUT_DIR, "shared_isoform_catalog.tsv.gz"), sep = "\t")
fwrite(isoform_map, file.path(OUT_DIR, "shared_isoform_map.tsv.gz"), sep = "\t")

