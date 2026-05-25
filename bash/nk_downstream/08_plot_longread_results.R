#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =============================================================================
# Build lab-meeting figures from downstream long-read analysis outputs.
#
# By default this script creates ggplot objects and prints them in sequence.
# If --save-dir is provided, plots are also saved as PDF files.
# =============================================================================

BASE_DIR <- "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs"
PROJECT_ID <- "msk"

WORKING_DIR <- file.path(BASE_DIR, PROJECT_ID)
ISOFORM_DIR <- file.path(WORKING_DIR, "downstream", "04_isoform_usage")
NOVEL_DIR <- file.path(WORKING_DIR, "downstream", "05_novel_isoforms")
SITE_DIR <- file.path(WORKING_DIR, "downstream", "06_tss_polya_usage")

read_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  fread(path)
}

extract_comparisons <- function() {
  files <- list.files(ISOFORM_DIR, pattern = "_isoform_results\\.tsv\\.gz$", full.names = FALSE)
  if (length(files) == 0) {
    return(character())
  }
  sub("_isoform_results\\.tsv\\.gz$", "", files)
}

load_comparison_data <- function(comparison) {
  list(
    isoform = read_if_exists(file.path(ISOFORM_DIR, paste0(comparison, "_isoform_results.tsv.gz"))),
    gene = read_if_exists(file.path(ISOFORM_DIR, paste0(comparison, "_gene_summary.tsv.gz"))),
    novel = read_if_exists(file.path(NOVEL_DIR, paste0(comparison, "_novel_isoform_candidates.tsv.gz"))),
    tss = read_if_exists(file.path(SITE_DIR, paste0(comparison, "_tss_results.tsv.gz"))),
    polya = read_if_exists(file.path(SITE_DIR, paste0(comparison, "_polya_results.tsv.gz")))
  )
}

make_hit_flag <- function(dt, p_col = "adj_p_value", effect_col = "mean_delta_prop", p_cutoff = 0.05, effect_cutoff = 0.05) {
  dt <- copy(dt)
  dt[, hit := fifelse(!is.na(get(p_col)) & get(p_col) <= p_cutoff & abs(get(effect_col)) >= effect_cutoff, "hit", "other")]
  dt
}

plot_isoform_volcano <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) {
    return(NULL)
  }
  dt <- make_hit_flag(dt)
  dt[, neglog10_adj_p := -log10(pmax(adj_p_value, 1e-300))]

  ggplot(dt, aes(x = mean_delta_prop, y = neglog10_adj_p, color = hit)) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = c(hit = "#d73027", other = "#636363")) +
    labs(x = "Mean delta isoform proportion (case - ref)", y = "-log10(FDR)", color = "")
}

plot_top_isoforms <- function(dt, n_top = 20) {
  if (is.null(dt) || nrow(dt) == 0) {
    return(NULL)
  }

  dt <- copy(dt)
  dt[, rank_score := ifelse(is.na(adj_p_value), Inf, adj_p_value)]
  dt <- dt[order(rank_score, -abs(mean_delta_prop))]
  dt <- dt[seq_len(min(n_top, .N))]

  dt[, isoform_label := paste(gene_id, shared_isoform_id, sep = " | ")]
  dt[, isoform_label := factor(isoform_label, levels = rev(isoform_label))]
  dt[, direction := fifelse(mean_delta_prop >= 0, "higher_in_case", "higher_in_ref")]

  ggplot(dt, aes(x = isoform_label, y = mean_delta_prop, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(higher_in_case = "#1f78b4", higher_in_ref = "#b2df8a")) +
    labs(x = "Gene | shared isoform", y = "Mean delta isoform proportion", fill = "")
}

plot_top_genes <- function(dt, n_top = 20) {
  if (is.null(dt) || nrow(dt) == 0) {
    return(NULL)
  }

  dt <- copy(dt)
  dt[, rank_score := ifelse(is.na(min_feature_adj_p_value), Inf, min_feature_adj_p_value)]
  dt <- dt[order(rank_score, -max_abs_delta_prop)]
  dt <- dt[seq_len(min(n_top, .N))]

  dt[, gene_id := factor(gene_id, levels = rev(gene_id))]

  ggplot(dt, aes(x = gene_id, y = max_abs_delta_prop, fill = n_consistent_isoforms)) +
    geom_col() +
    coord_flip() +
    labs(x = "Gene", y = "Max |delta isoform proportion|", fill = "Consistent isoforms")
}

plot_novel_candidates <- function(dt, n_top = 20) {
  if (is.null(dt) || nrow(dt) == 0) {
    return(NULL)
  }

  dt <- copy(dt)
  dt <- dt[order(-consistent_direction, -abs(mean_delta_prop), -donors_detected_case, donors_detected_ref)]
  dt <- dt[seq_len(min(n_top, .N))]

  dt[, candidate := paste(gene_id, shared_isoform_id, sep = " | ")]
  dt[, candidate := factor(candidate, levels = rev(candidate))]
  dt[, direction := fifelse(mean_delta_prop >= 0, "higher_in_case", "higher_in_ref")]

  ggplot(dt, aes(x = candidate, y = mean_delta_prop, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(higher_in_case = "#ff7f00", higher_in_ref = "#6a3d9a")) +
    labs(x = "Novel candidate", y = "Mean delta isoform proportion", fill = "")
}

plot_novel_category_counts <- function() {
  dt <- read_if_exists(file.path(NOVEL_DIR, "recurrent_novel_isoforms.tsv.gz"))
  if (is.null(dt) || nrow(dt) == 0) {
    return(NULL)
  }

  counts <- dt[, .(n_isoforms = .N), by = category]
  setorder(counts, -n_isoforms)
  counts[, category := factor(category, levels = rev(category))]

  ggplot(counts, aes(x = category, y = n_isoforms)) +
    geom_col(fill = "#33a02c") +
    coord_flip() +
    labs(x = "Novel category", y = "Number of recurrent isoforms")
}

plot_site_usage <- function(dt, site_label, n_top = 20) {
  if (is.null(dt) || nrow(dt) == 0) {
    return(NULL)
  }

  dt <- copy(dt)
  dt <- make_hit_flag(dt)
  dt[, rank_score := ifelse(is.na(adj_p_value), Inf, adj_p_value)]
  dt <- dt[order(rank_score, -abs(mean_delta_prop))]
  dt <- dt[seq_len(min(n_top, .N))]

  dt[, feature_label := paste(gene_id, feature_id, sep = " | ")]
  dt[, feature_label := factor(feature_label, levels = rev(feature_label))]
  dt[, direction := fifelse(mean_delta_prop >= 0, "higher_in_case", "higher_in_ref")]

  ggplot(dt, aes(x = feature_label, y = mean_delta_prop, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(higher_in_case = "#e31a1c", higher_in_ref = "#a6cee3")) +
    labs(x = paste(site_label, "feature"), y = paste("Mean delta", site_label, "proportion"), fill = "")
}

plot_site_volcano <- function(dt, site_label) {
  if (is.null(dt) || nrow(dt) == 0) {
    return(NULL)
  }
  dt <- make_hit_flag(dt)
  dt[, neglog10_adj_p := -log10(pmax(adj_p_value, 1e-300))]

  ggplot(dt, aes(x = mean_delta_prop, y = neglog10_adj_p, color = hit)) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = c(hit = "#d73027", other = "#636363")) +
    labs(x = paste("Mean delta", site_label, "proportion (case - ref)"), y = "-log10(FDR)", color = "")
}

save_plot_list <- function(plot_list, save_dir, width = 8, height = 5.5) {
  if (is.null(save_dir) || save_dir == "") {
    return(invisible(NULL))
  }
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  for (nm in names(plot_list)) {
    p <- plot_list[[nm]]
    if (is.null(p)) {
      next
    }
    out <- file.path(save_dir, paste0(nm, ".pdf"))
    ggsave(filename = out, plot = p, width = width, height = height, units = "in", device = "pdf")
  }
}

parse_args <- function(args) {
  out <- list(comparison = NULL, n_top = 20L, save_dir = NULL)
  if (length(args) == 0) {
    return(out)
  }

  for (arg in args) {
    if (grepl("^--comparison=", arg)) {
      out$comparison <- sub("^--comparison=", "", arg)
    } else if (grepl("^--n-top=", arg)) {
      out$n_top <- as.integer(sub("^--n-top=", "", arg))
    } else if (grepl("^--save-dir=", arg)) {
      out$save_dir <- sub("^--save-dir=", "", arg)
    }
  }
  out
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  comparisons <- extract_comparisons()
  if (length(comparisons) == 0) {
    stop("No step-4 comparison result files found in ", ISOFORM_DIR)
  }

  comparison <- args$comparison
  if (is.null(comparison) || comparison == "") {
    comparison <- comparisons[[1]]
  }

  if (!(comparison %in% comparisons)) {
    stop("Comparison not found: ", comparison, ". Available: ", paste(comparisons, collapse = ", "))
  }

  res <- load_comparison_data(comparison)

  plot_list <- list(
    isoform_volcano = plot_isoform_volcano(res$isoform),
    isoform_top = plot_top_isoforms(res$isoform, n_top = args$n_top),
    gene_top = plot_top_genes(res$gene, n_top = args$n_top),
    novel_top = plot_novel_candidates(res$novel, n_top = args$n_top),
    novel_category_counts = plot_novel_category_counts(),
    tss_volcano = plot_site_volcano(res$tss, "TSS"),
    tss_top = plot_site_usage(res$tss, "TSS", n_top = args$n_top),
    polya_volcano = plot_site_volcano(res$polya, "polyA"),
    polya_top = plot_site_usage(res$polya, "polyA", n_top = args$n_top)
  )

  for (nm in names(plot_list)) {
    if (!is.null(plot_list[[nm]])) {
      print(plot_list[[nm]])
    }
  }

  save_plot_list(plot_list, args$save_dir)

  message("Comparison plotted: ", comparison)
  message("Plots created: ", sum(vapply(plot_list, function(x) !is.null(x), logical(1))))
  if (!is.null(args$save_dir) && args$save_dir != "") {
    message("Saved plots to: ", args$save_dir)
  }

  invisible(plot_list)
}

main()