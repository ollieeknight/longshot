#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ===========================================================================
# INTERACTIVE VAF COMPARISON SCRIPT
# ===========================================================================
# USAGE:
#   In R session:
#   > source("09_plot_population_variants.R")
#   > vaf_data <- load_all_vcfs()
#   > explore_vaf(vaf_data)
#   > plot_comparison(vaf_data, donor = "HC01")
# ===========================================================================

# ---- CONFIG ----
WORK_DIR <- "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk"
JOINT_DIR <- file.path(WORK_DIR, "downstream", "07_population_variants", "joint")
OUT_DIR <- file.path(JOINT_DIR, "exploration")

# ---- FILTERS (edit these to change) ----
QUAL_THRESHOLD <- 20
DEPTH_THRESHOLD <- 10
VAF_THRESHOLD <- 0.05

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- HELPER: Extract genotype and calculate VAF ----
parse_genotype <- function(gt_string) {
  # Parse GT string (e.g., "0/1", "1|1", ".", etc.)
  if (is.na(gt_string) || gt_string == "." || gt_string == "") {
    return(list(alleles = NA, vaf = NA))
  }
  
  # Replace / and | with space for easier splitting
  gt_clean <- gsub("[|/]", " ", gt_string)
  alleles <- as.numeric(strsplit(gt_clean, " ")[[1]])
  alleles <- alleles[!is.na(alleles)]
  
  if (length(alleles) == 0) {
    return(list(alleles = NA, vaf = NA))
  }
  
  alt_count <- sum(alleles > 0)
  total_alleles <- length(alleles)
  vaf <- alt_count / total_alleles
  
  return(list(alleles = alleles, vaf = vaf))
}

# ---- HELPER: Read single VCF directly (no external tools) ----
read_vcf_vaf <- function(vcf_path, donor_name) {
  if (!file.exists(vcf_path)) {
    message("  File not found: ", vcf_path)
    return(data.table())
  }
  
  # Open gzipped VCF and read all lines
  con <- gzfile(vcf_path, "r")
  on.exit(close(con))
  vcf_lines <- readLines(con)
  
  # Find header line (starts with #CHROM)
  header_idx <- which(grepl("^#CHROM", vcf_lines))
  
  if (length(header_idx) == 0) {
    warning("No #CHROM header in ", vcf_path)
    return(data.table())
  }
  
  header_line <- vcf_lines[header_idx]
  header <- strsplit(header_line, "\t")[[1]]
  sample_names <- header[10:length(header)]
  
  # Get data lines (not comments)
  data_lines <- vcf_lines[-c(1:header_idx)]
  data_lines <- data_lines[!grepl("^#", data_lines)]
  
  if (length(data_lines) == 0) {
    message("  No variants in ", basename(vcf_path))
    return(data.table())
  }
  
  # Parse each line
  variant_list <- list()
  
  for (line in data_lines) {
    fields <- strsplit(line, "\t")[[1]]
    
    # Standard VCF columns
    chrom <- fields[1]
    pos <- as.numeric(fields[2])
    id <- fields[3]
    ref <- fields[4]
    alt <- fields[5]
    qual <- as.numeric(fields[6])
    filter <- fields[7]
    info <- fields[8]
    
    # Apply quality filters
    if (is.na(qual) || qual < QUAL_THRESHOLD) next
    
    # Extract DP from INFO field (try DP=XXX or extract from FORMAT)
    depth <- NA
    if (grepl("DP=", info)) {
      dp_match <- regmatches(info, regexpr("DP=[0-9]+", info))
      if (length(dp_match) > 0) {
        depth <- as.numeric(sub("DP=", "", dp_match))
      }
    }
    
    # If DP not in INFO, try to get from AD (if available)
    # For now, we'll accept variants and use NA for missing depth
    # (this is a fallback; ideally depth would be in INFO)
    
    if (!is.na(depth) && depth < DEPTH_THRESHOLD) next
    
    # Genotype columns start at position 10
    genotype_cols <- fields[10:length(fields)]
    
    for (j in seq_along(sample_names)) {
      if (j > length(genotype_cols)) break
      
      sample_col <- sample_names[j]
      gt_string <- genotype_cols[j]
      
      # Parse genotype
      gt_parsed <- parse_genotype(gt_string)
      vaf <- gt_parsed$vaf
      
      # Apply VAF filter
      if (!is.na(vaf) && vaf >= VAF_THRESHOLD) {
        # Extract population from sample name (e.g., "HC01_adaptive_nk" -> "adaptive_nk")
        population <- sub("^[^_]+_", "", sample_col)
        
        variant_list[[paste(chrom, pos, sample_col, sep="_")]] <- data.table(
          chrom = chrom,
          pos = pos,
          ref = ref,
          alt = alt,
          donor = donor_name,
          population = population,
          sample = sample_col,
          genotype = gt_string,
          depth = depth,
          qual = qual,
          vaf = vaf
        )
      }
    }
  }
  
  if (length(variant_list) > 0) {
    rbindlist(variant_list)
  } else {
    data.table()
  }
}

# ---- MAIN FUNCTION: Load all VCFs ----
load_all_vcfs <- function() {
  message("Finding master VCF files...")
  vcf_files <- list.files(JOINT_DIR, pattern = "_master.vcf.gz$", full.names = TRUE, recursive = TRUE)
  
  if (length(vcf_files) == 0) {
    stop("No master VCF files found in ", JOINT_DIR)
  }
  
  donors <- sub("_master.vcf.gz$", "", basename(vcf_files))
  message("Found ", length(donors), " donors: ", paste(donors, collapse = ", "))
  message("")
  
  all_vaf <- list()
  
  for (i in seq_along(vcf_files)) {
    message("Reading ", donors[i], "...")
    vaf_dt <- read_vcf_vaf(vcf_files[i], donors[i])
    if (nrow(vaf_dt) > 0) {
      all_vaf[[i]] <- vaf_dt
      message("  ✓ Loaded ", nrow(vaf_dt), " variants")
    }
  }
  
  if (length(all_vaf) == 0) {
    stop("No variants passed filters!")
  }
  
  vaf_data <- rbindlist(all_vaf, fill = TRUE)
  
  message("\n✓ Total: ", nrow(vaf_data), " variants across ",
          length(unique(vaf_data$donor)), " donors and ",
          length(unique(vaf_data$population)), " populations")
  
  return(vaf_data)
}

# ---- HELPER: Explore data ----
explore_vaf <- function(vaf_data) {
  message("\n=== VARIANT SUMMARY ===")
  
  summary <- vaf_data[, .(
    n_variants = .N,
    median_vaf = median(vaf, na.rm = TRUE),
    mean_vaf = mean(vaf, na.rm = TRUE),
    mean_depth = mean(depth, na.rm = TRUE)
  ), by = .(donor, population)]
  
  setorder(summary, donor, -n_variants)
  print(summary)
  
  message("\n=== UNIQUE VALUES ===")
  message("Donors: ", paste(unique(vaf_data$donor), collapse = ", "))
  message("Populations: ", paste(unique(vaf_data$population), collapse = ", "))
  message("VAF range: ", round(min(vaf_data$vaf, na.rm = TRUE), 3), " - ", 
          round(max(vaf_data$vaf, na.rm = TRUE), 3))
  
  invisible(summary)
}

# ---- HELPER: Create comparison plot ----
plot_comparison <- function(vaf_data, donor = NULL, pop1 = NULL, pop2 = NULL) {
  subset_data <- vaf_data
  
  if (!is.null(donor)) {
    subset_data <- subset_data[donor == donor]
  }
  
  if (!is.null(pop1) && !is.null(pop2)) {
    subset_data <- subset_data[population %in% c(pop1, pop2)]
  }
  
  p <- ggplot(subset_data, aes(x = population, y = vaf, fill = population)) +
    geom_boxplot(show.legend = FALSE, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
    facet_wrap(~donor) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Population", y = "Variant Allele Frequency")
  
  return(p)
}

# ---- HELPER: Save outputs ----
save_outputs <- function(vaf_data, vaf_summary) {
  fwrite(vaf_data, file.path(OUT_DIR, "all_variants_vaf.tsv.gz"), sep = "\t")
  fwrite(vaf_summary, file.path(OUT_DIR, "vaf_summary_by_population.tsv.gz"), sep = "\t")
  message("\nSaved to: ", OUT_DIR)
}

# ---- AUTO-RUN (if script is executed, not sourced) ----
if (!interactive()) {
  message("Loading VAF data...")
  vaf_data <- load_all_vcfs()
  
  message("\nGenerating summary...")
  vaf_summary <- vaf_data[, .(
    n_variants = .N,
    median_vaf = median(vaf, na.rm = TRUE),
    mean_vaf = mean(vaf, na.rm = TRUE),
    mean_depth = mean(depth, na.rm = TRUE)
  ), by = .(donor, population)]
  
  explore_vaf(vaf_data)
  save_outputs(vaf_data, vaf_summary)
  
  message("\n✓ Done. Use in interactive R:")
  message("  > source('09_plot_population_variants.R')")
  message("  > vaf_data <- load_all_vcfs()")
  message("  > p <- plot_comparison(vaf_data, donor = 'HC01')")
}

# ---- INTERACTIVE STEP-BY-STEP COMMANDS (paste these in R one at a time) ----
# 
# Step 1: Load just one VCF to test
# vcf_path <- "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk/downstream/07_population_variants/joint/HC01/HC01_master.vcf.gz"
# hc01_vaf <- read_vcf_vaf(vcf_path, "HC01")
# dim(hc01_vaf)
# head(hc01_vaf)
#
# Step 2: Load all three donors
# vaf_data <- load_all_vcfs()
# head(vaf_data, 20)
#
# Step 3: Explore the data
# explore_vaf(vaf_data)
#
# Step 4: Look at unique populations
# unique(vaf_data$population)
# table(vaf_data$donor, vaf_data$population)
#
# Step 5: Create a plot for one donor
# p <- plot_comparison(vaf_data, donor = "HC01")
# print(p)
#
# Step 6: Compare two specific populations
# p2 <- plot_comparison(vaf_data, donor = "HC01", pop1 = "adaptive_nk", pop2 = "cd56dim_nk")
# print(p2)
#
# Step 7: Save outputs
# vaf_summary <- vaf_data[, .(
#   n_variants = .N,
#   median_vaf = median(vaf, na.rm = TRUE),
#   mean_vaf = mean(vaf, na.rm = TRUE),
#   mean_depth = mean(depth, na.rm = TRUE)
# ), by = .(donor, population)]
# save_outputs(vaf_data, vaf_summary)

