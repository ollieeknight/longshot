repo <- '~/work/scripts/nk_xcl1/'

source(paste0(repo, 'utilities/environment.R'))
options(future.globals.maxSize = 199 * 1024^3)
source(paste0(repo, 'utilities/helper_functions.R'))
source(paste0(repo, 'utilities/themes.R'))

setwd('~/work/data')

## Import and quality control

import_sqanti_seurat <- function(
    base_dir = NULL, isoquant_dir = NULL, sqanti_dir = NULL, sample_mapping_csv = NULL,
    gtf_file = NULL
) {
  
  # 1. Define paths
  gene_mtx      <- file.path(isoquant_dir, "isoquant_out", "joint", "joint.gene_grouped_barcode_counts.matrix.mtx")
  gene_barcodes <- file.path(isoquant_dir, "isoquant_out", "joint", "joint.gene_grouped_barcode_counts.barcodes.tsv")
  gene_features <- file.path(isoquant_dir, "isoquant_out", "joint", "joint.gene_grouped_barcode_counts.features.tsv")
  
  iso_mtx       <- file.path(isoquant_dir, "isoquant_out", "joint", "joint.transcript_grouped_barcode_counts.matrix.mtx")
  iso_barcodes  <- file.path(isoquant_dir, "isoquant_out", "joint", "joint.transcript_grouped_barcode_counts.barcodes.tsv")
  iso_features  <- file.path(isoquant_dir, "isoquant_out", "joint", "joint.transcript_grouped_barcode_counts.features.tsv")
  
  sqanti_filter_file <- file.path(sqanti_dir, "sqanti_filter", "joint_filtered_RulesFilter_result_classification.txt")
  sqanti_qc_file     <- file.path(sqanti_dir, "sqanti_qc", "joint_classification.txt")
  
  # 2. Parse SQANTI filtering
  # We bypass the aggressive SQANTI3 filter file and use the raw QC file
  # Keep only high-confidence categories manually (FSM, ISM, NIC, NNC)
  sqanti_class <- fread(sqanti_qc_file)
  
  # Optional: Filter for reliable structural categories to drop intergenic/antisense noise
  sqanti_class <- sqanti_class[structural_category %in% c("full-splice_match", 
                                                          "novel_in_catalog", 
                                                          "novel_not_in_catalog", 
                                                          "incomplete-splice_match")]
  
  # Predict Seurat's upcoming forced conversion to dashes
  sqanti_class$isoform <- gsub("_", "-", sqanti_class$isoform)
  
  passing_transcripts <- sqanti_class$isoform
  message(paste("Found", length(passing_transcripts), "transcripts passing manual structural category filtering."))
  
  # 3. Load Matrix & Subsetting
  message("Loading IsoQuant sparse matrices...")
  counts_sparse <- readMM(gene_mtx)
  barcodes      <- fread(gene_barcodes, header = FALSE)[[1]]
  features_dt   <- fread(gene_features, header = FALSE)
  
  raw_ensg_ids  <- features_dt[[1]]
  
  if (!is.null(gtf_file) && file.exists(gtf_file)) {
    message("Parsing GTF to map Ensembl IDs to Gene Symbols...")
    # Extremely fast read filtering for 'gene' entries directly via awk
    gtf <- fread(cmd = paste0("awk -F'\\t' '$3 == \"gene\"' ", gtf_file), header=FALSE, sep="\t")
    # Extract attributes using regex
    gtf[, ensg_id := sub('.*gene_id "([^"]+)".*', '\\1', V9)]
    gtf[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', V9)]
    
    map_dict <- setNames(gtf$gene_name, gtf$ensg_id)
    
    # Match against features
    mapped_names <- map_dict[raw_ensg_ids]
    
    # Handle internal Ensembl versioning discrepancies
    missing_idx <- is.na(mapped_names)
    if (any(missing_idx)) {
      base_ensg_raw <- sub("\\.[0-9]+$", "", raw_ensg_ids[missing_idx])
      base_ensg_gtf <- sub("\\.[0-9]+$", "", gtf$ensg_id)
      base_map_dict <- setNames(gtf$gene_name, base_ensg_gtf)
      # Only update NAs if they exist in base GTF dict
      mapped_names[missing_idx] <- ifelse(is.na(base_map_dict[base_ensg_raw]), 
                                          raw_ensg_ids[missing_idx], 
                                          base_map_dict[base_ensg_raw])
    }
    
    # Fallback to ENSG if gene_name was absolutely completely missing from GTF
    mapped_names[is.na(mapped_names)] <- raw_ensg_ids[is.na(mapped_names)]
    
    message("Collapsing sparse matrix by Gene Symbol...")
    group_f <- factor(mapped_names)
    row_map <- sparseMatrix(
      i = as.integer(group_f),
      j = seq_along(group_f),
      x = 1,
      dims = c(nlevels(group_f), length(group_f)),
      dimnames = list(levels(group_f), NULL)
    )
    counts_sparse <- row_map %*% counts_sparse
    
  } else {
    # No collapsing, just unique rownames
    rownames(counts_sparse) <- make.unique(as.character(raw_ensg_ids))
  }
  
  colnames(counts_sparse) <- barcodes
  
  iso_sparse      <- readMM(iso_mtx)
  iso_barcodes_v  <- fread(iso_barcodes, header = FALSE)[[1]]
  iso_features_dt <- fread(iso_features, header = FALSE)
  
  # Inject dashes into raw counts matrix exactly matching Seurat
  iso_features_dt[[1]] <- gsub("_", "-", iso_features_dt[[1]])
  
  rownames(iso_sparse) <- iso_features_dt[[1]]
  colnames(iso_sparse) <- iso_barcodes_v
  
  iso_sparse <- iso_sparse[rownames(iso_sparse) %in% passing_transcripts, ]
  
  # 4. Construct Dual-Assay Seurat Object
  message("Creating Seurat Object tracking both Genes and Isoforms...")
  seurat_obj <- CreateSeuratObject(
    counts = counts_sparse, 
    project = "MSK_LongRead", 
    assay = "RNA",
    min.cells = 0,  
    min.features = 10 
  )
  
  valid_cells <- colnames(seurat_obj)
  iso_sparse_valid <- iso_sparse[, colnames(iso_sparse) %in% valid_cells]
  
  iso_assay <- CreateAssay5Object(counts = iso_sparse_valid)
  seurat_obj[["ISO"]] <- iso_assay
  
  # 5. Metadata Processing (Sample IDs via CSV)
  message("Adding sample metadata...")
  cell_names <- colnames(seurat_obj)
  seurat_obj$cell_barcode <- cell_names
  
  if (file.exists(sample_mapping_csv)) {
    # Extract the suffix from the cell barcode (e.g. "_01" from "AACCAAGGAAACCCCA_01")
    meta_df <- data.table(cell_barcode = cell_names)
    meta_df[, suffix := paste0("_", sub(".*_", "", cell_barcode))]
    
    sample_map <- fread(sample_mapping_csv)
    
    merged <- merge(meta_df, sample_map, by="suffix", all.x=TRUE)
    merged <- merged[match(cell_names, merged$cell_barcode)]
    seurat_obj$sample_id <- merged$sample_id
  }
  
  # 6. Attach SQANTI Feature Labels
  sqanti_meta <- as.data.frame(sqanti_class)
  rownames(sqanti_meta) <- sqanti_meta$isoform
  valid_iso_meta <- sqanti_meta[rownames(seurat_obj[["ISO"]]), ]
  
  if (exists("base_map_dict")) {
    meta_base_ensg <- sub("\\.[0-9]+$", "", valid_iso_meta$associated_gene)
    valid_iso_meta$associated_gene_symbol <- ifelse(
      !is.na(base_map_dict[meta_base_ensg]), 
      base_map_dict[meta_base_ensg], 
      valid_iso_meta$associated_gene
    )
  } else {
    valid_iso_meta$associated_gene_symbol <- valid_iso_meta$associated_gene
  }
  
  seurat_obj[["ISO"]]@meta.data <- valid_iso_meta
  
  return(seurat_obj)
}

alldata <- import_sqanti_seurat(
  base_dir = "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk.",
  isoquant_dir = "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk/workspace/isoquant/",
  sqanti_dir = "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk/workspace/sqanti3/",
  sample_mapping_csv = "/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/msk/workspace/pacbio/sample_mapping.csv",
  gtf_file = "~/group/work/ref/hs/pacbio/reference/gencode.v39.annotation.gtf"
)

head(rownames(alldata[['RNA']]))
head(rownames(alldata[['ISO']]))

alldata$sample_id <- factor(alldata$sample_id, levels = sort(unique(alldata$sample_id)))
libraries <- SplitObject(alldata, split.by = 'sample_id')

libraries[[1]]$donor_id <- 'HC02'
libraries[[2]]$donor_id <- 'HC01'
libraries[[3]]$donor_id <- 'HC03'

libraries <- add_qc_metrics(libraries, 
                            gtf_path = '~/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/genes/genes.gtf.gz')

for (i in 1:length(libraries)) {
  libraries[[i]] <- subset(libraries[[i]], percent_mitochondrial > 0)
  libraries[[i]] <- subset(libraries[[i]], percent_ribosomal > 0)
}

plot_metrics_x_y(libraries[[1]], x_metric = 'nCount_RNA', y_metric = 'percent_ribosomal', x_threshold = 1500, y_threshold = 6.0)
selected <- WhichCells(libraries[[1]], expression = nCount_RNA >= 1500 & percent_ribosomal >= 6.0)
libraries[[1]] <- subset(libraries[[1]], cells = selected)

plot_metrics_x_y(libraries[[2]], x_metric = 'nCount_RNA', y_metric = 'percent_ribosomal', x_threshold = 1500, y_threshold = 6.0)
selected <- WhichCells(libraries[[2]], expression = nCount_RNA >= 1500 & percent_ribosomal >= 6.0)
libraries[[2]] <- subset(libraries[[2]], cells = selected)

plot_metrics_x_y(libraries[[3]], x_metric = 'nCount_RNA', y_metric = 'percent_ribosomal', x_threshold = 1500, y_threshold = 6.0)
selected <- WhichCells(libraries[[3]], expression = nCount_RNA >= 1500 & percent_ribosomal >= 6.0)
libraries[[3]] <- subset(libraries[[3]], cells = selected)

plot_metric_violin(libraries[[1]], metric = 'percent_mitochondrial', threshold = 2.5)
libraries[[1]] <- subset(libraries[[1]], percent_mitochondrial <= 2.5)

plot_metric_violin(libraries[[2]], metric = 'percent_mitochondrial', threshold = 2.5)
libraries[[2]] <- subset(libraries[[2]], percent_mitochondrial <= 2.5)

plot_metric_violin(libraries[[3]], metric = 'percent_mitochondrial', threshold = 2.5)
libraries[[3]] <- subset(libraries[[3]], percent_mitochondrial <= 2.5)

## Normalise data

# Normalise RNA matrix according to https://doi.org/10.1101/2022.05.06.490859
libraries <- proportional_normalise_seurat_list(libraries)

## DoubletFinder removal

libraries <- run_doubletfinder_seurat_list(libraries)

# # Save processed libraries
qs::qsave(libraries, 'adaptive_nk/objects/raw/Kinnex_10x_NK_libraries.qs')

## Dimensionality reduction and annotations

plan(multisession, workers = 32)

total <- Reduce(merge, libraries)
total <- join_assay_layers(total)

plot_metric_violin(total, metric = 'percent_HSP', threshold = 2.5) + scale_y_log10()
total <- subset(total, percent_HSP < 2.5)

plot_metric_violin(total, metric = 'percent_HB', threshold = 0.001) + scale_y_log10()
total <- subset(total, percent_HB < 0.001)

total <- subset(total, rna_PPBP < 0.001)

total <- map_to_azimuth_pbmc_reference(total, assay = 'RNA')

gc()

total <- integrate_seurat_by_harmony(total, 'donor_id',
                                     nPCs = 30,
                                     cluster_resolution = 0.5,
                                     n_threads = 32)

p1 <- plotalot(total, metadata_column = 'azimuth_pbmc_reference_celltype_l2')
p2 <- plotalot(total, metadata_column = 'seurat_clusters')

p1 | p2

total <- subset(total, idents = c(9, 10), invert = T)

donordata <- SplitObject(total, split.by = 'donor_id')

for (i in 1:length(donordata)) {
  donordata[[i]] <- join_assay_layers(donordata[[i]])
  
  DefaultAssay(donordata[[i]]) <- 'RNA'
  donordata[[i]] <- DietSeurat(donordata[[i]]) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%
    FindNeighbors(reduction = 'pca', 
                  dims = 1:30, 
                  k.param = 30) %>%
    FindClusters(graph.name = 'RNA_snn', 
                 resolution = 0.5, 
                 algorithm = 4,
                 leiden_method = 'igraph',
                 group.singletons = TRUE, 
                 random.seed = 123)
  
  pca_emb <- Embeddings(donordata[[i]], 'pca')[, 1:30]
  umap_model <- uwot::umap2(pca_emb,
                            n_neighbors = 30, 
                            n_components = 2,
                            metric = 'cosine', 
                            min_dist = 0.3, 
                            ret_model = TRUE,
                            n_threads = 32, 
                            verbose = FALSE)
  
  donordata[[i]][['umap_rna']] <- CreateDimReducObject(
    embeddings = umap_model$embedding,
    key = 'pca_',
    assay = 'RNA'
  )
  donordata[[i]][['umap_rna']]@misc$model <- umap_model
  
  donordata[[i]] <- join_assay_layers(donordata[[i]])
}

p1 <- plotalot(donordata[[1]], metadata_column = 'azimuth_pbmc_reference_celltype_l2')
p2 <- plotalot(donordata[[1]], metadata_column = 'seurat_clusters')

p1 | p2

total$celltype_l1 <- 'Immune'
total$celltype_l2 <- 'NK/ILC'
total$celltype_l3 <- as.character('NK')
total$celltype_l4 <- as.character(total$celltype_l3)
total$celltype_l5 <- as.character(total$celltype_l4)

selected <- colnames(donordata[[1]])
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'CD56dim NK', total$celltype_l4)
selected <- WhichCells(donordata[[1]], idents = c(3))
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'Adaptive NK', total$celltype_l4)
total$celltype_l5 <- ifelse(colnames(total) %in% selected, 'Adaptive NK 1', total$celltype_l5)
selected <- WhichCells(donordata[[1]], idents = c(4))
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'Adaptive NK', total$celltype_l4)
total$celltype_l5 <- ifelse(colnames(total) %in% selected, 'Adaptive NK 2', total$celltype_l5)

plotalot(total, 'celltype_l4')
plotalot(total, 'celltype_l5')

p1 <- plotalot(donordata[[2]], metadata_column = 'azimuth_pbmc_reference_celltype_l2')
p2 <- plotalot(donordata[[2]], metadata_column = 'seurat_clusters')

p1 | p2

donordata[[2]] <- SetIdent(donordata[[2]], value = 'seurat_clusters')
donordata[[2]] <- FindSubCluster(donordata[[2]], cluster = 5, graph.name = 'RNA_snn', resolution = 0.2, algorithm = 4)
donordata[[2]] <- SetIdent(donordata[[2]], value = 'sub.cluster')

p2 <- plotalot(donordata[[2]], metadata_column = 'sub.cluster')

p1 | p2

selected <- colnames(donordata[[2]])
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'CD56dim NK', total$celltype_l4)
selected <- WhichCells(donordata[[2]], idents = c(4, '5_2'))
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'Adaptive NK', total$celltype_l4)
total$celltype_l5 <- ifelse(colnames(total) %in% selected, 'Adaptive NK 3', total$celltype_l5)

plotalot(total, 'celltype_l4')
plotalot(total, 'celltype_l5')

p1 <- plotalot(donordata[[3]], metadata_column = 'azimuth_pbmc_reference_celltype_l2')
p2 <- plotalot(donordata[[3]], metadata_column = 'seurat_clusters')

p1 | p2

donordata[[3]] <- SetIdent(donordata[[3]], value = 'seurat_clusters')
donordata[[3]] <- FindSubCluster(donordata[[3]], cluster = 3, graph.name = 'RNA_snn', resolution = 0.1, algorithm = 4)
donordata[[3]] <- SetIdent(donordata[[3]], value = 'sub.cluster')

p2 <- plotalot(donordata[[3]], metadata_column = 'sub.cluster')

p1 | p2

selected <- colnames(donordata[[3]])
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'CD56dim NK', total$celltype_l4)
selected <- WhichCells(donordata[[3]], idents = c(2, '3_1'))
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'Adaptive NK', total$celltype_l4)
total$celltype_l5 <- ifelse(colnames(total) %in% selected, 'Adaptive NK 4', total$celltype_l5)
selected <- WhichCells(donordata[[3]], idents = c('3_2'))
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'Adaptive NK', total$celltype_l4)
total$celltype_l5 <- ifelse(colnames(total) %in% selected, 'Adaptive NK 5', total$celltype_l5)
selected <- WhichCells(donordata[[3]], idents = c(4))
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'Adaptive NK', total$celltype_l4)
total$celltype_l5 <- ifelse(colnames(total) %in% selected, 'Adaptive NK 6', total$celltype_l5)

plotalot(total, 'celltype_l4')

total <- integrate_seurat_by_harmony(total, 'donor_id',
                                     nPCs = 30,
                                     cluster_resolution = 0.5,
                                     n_threads = 32)

p1 <- plotalot(total, metadata_column = 'azimuth_pbmc_reference_celltype_l2')
p2 <- plotalot(total, metadata_column = 'seurat_clusters')

p1 | p2

total <- SetIdent(total, value = 'seurat_clusters')
total <- FindSubCluster(total, cluster = 4, graph.name = 'RNA_snn', resolution = 0.5, algorithm = 4)
total <- SetIdent(total, value = 'sub.cluster')
total <- FindSubCluster(total, cluster = '4_2', graph.name = 'RNA_snn', resolution = 0.5, algorithm = 4)
total <- SetIdent(total, value = 'sub.cluster')

p2 <- plotalot(total, metadata_column = 'sub.cluster')

p1 | p2

selected <- WhichCells(total, idents = c(7, 8))
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'Cycling NK', total$celltype_l4)
total$celltype_l5 <- ifelse(colnames(total) %in% selected, 'Cycling NK', total$celltype_l5)
selected <- WhichCells(total, idents = c('4_1', '4_2_1', '4_3'))
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'CD56bright NK', total$celltype_l4)
total$celltype_l5 <- ifelse(colnames(total) %in% selected, 'CD56bright NK', total$celltype_l5)
selected <- WhichCells(total, idents = c('4_2_2'))
total$celltype_l3 <- ifelse(colnames(total) %in% selected, 'ILC', total$celltype_l3)
total$celltype_l4 <- ifelse(colnames(total) %in% selected, 'Naïve-like ILC', total$celltype_l4)
total$celltype_l5 <- ifelse(colnames(total) %in% selected, 'Naïve-like ILC', total$celltype_l5)

total$celltype_l5[total$celltype_l4 == 'CD56dim NK'] <- 'CD56dim NK'

plotalot(total, metadata_column = 'celltype_l3')
plotalot(total, metadata_column = 'celltype_l4')
plotalot(total, metadata_column = 'celltype_l5')

## Export for some long read bam stuff to try

barcode_output_dir <- '~/work/data/adaptive_nk/objects/longread/barcodes'

clean_name <- function(x) {
  x |>
    gsub('ï', 'i', x = _) |>
    tolower() |>
    gsub('[^a-z0-9]+', '_', x = _) |>
    gsub('^_|_$', '', x = _)
}

for (donor in unique(total$donor_id)) {
  
  donor_cells <- total[, total$donor_id == donor]
  barcodes <- sub('-.*', '', colnames(donor_cells))
  cells_df <- data.frame(barcode = barcodes, label = clean_name(donor_cells$celltype_l4))
    write.table(cells_df, barcode_output_dir, paste0(donor, '.tsv'), sep = '\t', quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

total <- join_assay_layers(total)

total$UMAP_1 <- total@reductions$umap_rna@cell.embeddings[, 1]
total$UMAP_2 <- -total@reductions$umap_rna@cell.embeddings[, 2]

total$celltype_l4 <- factor(total$celltype_l4, levels =
                              c('CD56bright NK', 'CD56dim NK', 'Adaptive NK', 'Cycling NK', 'Naïve-like ILC'))

qs::qsave(total, 'adaptive_nk/objects/Kinnex_10x_PBMC_NK.qs')

ggplot(shufflepoints(total@meta.data), 
       aes(x = UMAP_1, y = UMAP_2, color = celltype_l4)) +
  geom_point(size = 0.1) + 
  scale_colour_tableau() +
  themes$UMAP +
  add_umap_coords(total@meta.data) +
  theme(
    legend.position = 'bottom',
    legend.direction = 'horizontal')  +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

df <- total@meta.data %>%
  group_by(donor_id, celltype_l4) %>%
  summarise(
    n_cells = n(), 
    .groups = 'drop_last'
  ) %>%
  mutate(total_cells_donor = sum(n_cells)) %>%
  mutate(percentage = (n_cells / total_cells_donor) * 100) %>%
  ungroup() %>%
  filter(n_cells > 1)

df$celltype_l4 <- factor(df$celltype_l4, levels =
                           c('CD56bright NK', 'CD56dim NK', 'Adaptive NK', 'Cycling NK', 'Naïve-like ILC'))

df_summary <- df %>%
  group_by(celltype_l4) %>%
  summarise(
    mean_val = mean(percentage, na.rm = TRUE),
    n = n(),
    sd = sd(percentage, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = 'drop'
  )

ggplot(df, aes(x = celltype_l4, y = percentage, colour = celltype_l4)) +
  geom_col(data = df_summary, aes(y = mean_val, fill = celltype_l4), 
           alpha = 0.6, width = 0.6, colour = NA) +
  geom_errorbar(data = df_summary, aes(y = mean_val, ymin = mean_val - sd, ymax = mean_val + sd),
                width = 0.15, colour = 'black', linewidth = 0.5) +
  geom_jitter(data = df, aes(fill = celltype_l4), 
              height = 0, width = 0.2, size = 3, alpha = 1,
              stroke = 0.5, colour = 'black', shape = 21) +
  labs(
    y = 'Subset (%)',
    x = NULL,
    fill = NULL,
    colour = NULL
  ) +
  scale_fill_tableau() +
  scale_colour_tableau() +
  theme(
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 8)
    ) +
  scale_y_continuous(expand = c(0, 0))

a <- plot_umap_feature(total, 'KIT', modality = 'rna', colour_scale = scale$rna, cutoff = 0.999)
b <- plot_umap_feature(total, 'KLRC2', modality = 'rna', colour_scale = scale$rna, cutoff = 0.999)
c <- plot_umap_feature(total, 'FCGR3A', modality = 'rna', colour_scale = scale$rna, cutoff = 0.999)
d <- plot_umap_feature(total, 'MKI67', modality = 'rna', colour_scale = scale$rna, cutoff = 0.999)

a + b + c + d + plot_layout(ncol = 2, nrow = 2, guides = 'collect')

alldata <- subset(total, donor_id == 'HC03')

DefaultAssay(alldata) <- 'RNA'
alldata <- DietSeurat(alldata) %>%
  FindVariableFeatures(verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(npcs = 50, verbose = F) %>%
  RunUMAP(reduction = 'pca', dims = 1:30, n.components = 2, 
          umap.method = 'uwot', verbose = F, reduction.name = 'umap_rna', 
          return.model = T)

alldata$UMAP_1 <- alldata@reductions[['umap_rna']]@cell.embeddings[, 1]
alldata$UMAP_2 <- alldata@reductions[['umap_rna']]@cell.embeddings[, 2]

ggplot(shufflepoints(alldata@meta.data), 
       aes(x = UMAP_1, y = UMAP_2, color = celltype_l4)) +
  geom_point(size = 0.1) + 
  scale_colour_tableau() +
  themes$UMAP +
  add_umap_coords(total@meta.data) +
  theme(
    legend.position = 'bottom',
    legend.direction = 'horizontal')  +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

# DE - Improved Timo's version to tease out the adaptive signature

alldata$analysis_group <- case_when(
  alldata$celltype_l4 == 'Adaptive NK'~ 'Adaptive NK',
  alldata$celltype_l4 == 'CD56dim NK'~ 'CD56dim NK',
  alldata$celltype_l4 == 'CD56bright NK'  ~ 'CD56bright NK',
  TRUE ~ 'Exclude'
)

alldata$pseudobulk_id <- ifelse(alldata$analysis_group != 'Exclude',
                                paste(alldata$analysis_group, alldata$donor_id, sep = '_'),
                                'Exclude')

alldata_filtered <- subset(alldata, pseudobulk_id != 'Exclude')

comparisons <- list(
  list('bright_v_dim', 'CD56bright NK', 'CD56dim NK'),
  list('bright_v_adaptive', 'CD56bright NK', 'Adaptive NK'),
  list('adaptive_v_dim', 'Adaptive NK', 'CD56dim NK')
)

run_comparison <- function(comp, obj) {
  message('Running ', comp[[1]])
  groups <- c(comp[[2]], comp[[3]])
  
  sub_counts <- AggregateExpression(obj,
                                    group.by = 'pseudobulk_id',
                                    assays = 'RNA',
                                    slot = 'counts',
                                    return.seurat = FALSE)$RNA
  
  sub_meta <- obj@meta.data %>%
    dplyr::select(donor_id, analysis_group, pseudobulk_id) %>%
    tidyr::drop_na() %>%
    dplyr::distinct(pseudobulk_id, .keep_all = TRUE) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames('pseudobulk_id')
  
  rownames(sub_meta) <- gsub('[^[:alnum:]]', '', rownames(sub_meta))
  colnames(sub_counts) <- gsub('[^[:alnum:]]', '', colnames(sub_counts))
  
  common_ids <- intersect(rownames(sub_meta), colnames(sub_counts))
  sub_meta <- sub_meta[common_ids, ]
  sub_counts <- sub_counts[, common_ids]
  
  sub_meta$donor_id <- as.factor(sub_meta$donor_id)
  sub_meta$analysis_group <- factor(sub_meta$analysis_group, levels = groups)
  
  dds_sub <- DESeqDataSetFromMatrix(countData = sub_counts,
                                    colData = sub_meta,
                                    design = ~ donor_id + analysis_group)
  
  dds_sub <- dds_sub[rowSums(counts(dds_sub) >= 10) >= 3, ]
  dds_sub$analysis_group <- relevel(dds_sub$analysis_group, ref = groups[[2]])
  dds_sub <- DESeq(dds_sub)
  
  coef_name <- grep('^analysis_group_', resultsNames(dds_sub), value = TRUE)
  
  res <- results(dds_sub, name = coef_name, alpha = 0.01)
  
  as.data.frame(lfcShrink(dds_sub,
                          coef = coef_name,
                          res  = res,
                          type = 'apeglm'))
}

nkmarkers <- setNames(
  lapply(comparisons, function(comp) {
    groups <- c(comp[[2]], comp[[3]])
      cells_use <- alldata_filtered@meta.data %>%
        tibble::rownames_to_column('cell_barcode') %>%
        filter(analysis_group %in% groups) %>%
        pull(cell_barcode)
    
    obj <- subset(alldata_filtered, cells = cells_use)
    run_comparison(comp, obj)
  }),
  sapply(comparisons, `[[`, 1)
)

get_sig <- function(res, dir) {
  res <- res[!is.na(res$padj) & res$padj < 0.05, ]
  if (dir == 'up') {
    rownames(res[res$log2FoldChange > 1, ])
  } else {
    rownames(res[res$log2FoldChange < -1, ])
  }
}

bright_genes <- intersect(
  get_sig(nkmarkers$bright_v_dim, 'up'),
  get_sig(nkmarkers$bright_v_adaptive, 'up')
)

dim_genes <- intersect(
  get_sig(nkmarkers$bright_v_dim, 'down'),
  get_sig(nkmarkers$adaptive_v_dim, 'down')
)

adaptive_genes <- intersect(
  get_sig(nkmarkers$bright_v_adaptive, 'down'),
  get_sig(nkmarkers$adaptive_v_dim, 'up')
)

plot_genes <- unique(c(bright_genes, dim_genes, adaptive_genes))

pseudobulked_counts <- Seurat::AggregateExpression(alldata_filtered, 
                                         group.by = 'pseudobulk_id', 
                                         return.seurat = FALSE)$RNA

log_mat <- log2(sweep(pseudobulked_counts, 2, colSums(pseudobulked_counts), '/') * 10000 + 1)

pheatmap::pheatmap(log_mat[plot_genes, ], scale = 'row', color = scale$heatmap, show_rownames = F)

library(DEXSeq)

total$analysis_group <- case_when(
  total$celltype_l4 == 'Adaptive NK'~ 'Adaptive NK',
  total$celltype_l4 == 'CD56dim NK'~ 'CD56dim NK',
  total$celltype_l4 == 'CD56bright NK'  ~ 'CD56bright NK',
  TRUE ~ 'Exclude'
)

total$pseudobulk_id <- ifelse(total$analysis_group != 'Exclude',
                              paste(total$analysis_group, total$donor_id, sep = '_'),
                              'Exclude')

total_filtered <- subset(total, pseudobulk_id != 'Exclude')

comparisons <- list(
  # list('bright_v_dim', 'CD56bright NK', 'CD56dim NK'),
  # list('bright_v_adaptive', 'CD56bright NK', 'Adaptive NK'),
  list('adaptive_v_dim', 'Adaptive NK', 'CD56dim NK')
)

library(DEXSeq)

run_dexseq_comparison <- function(comp, obj) {
  message('Running DEXSeq DIU for ', comp[[1]])
  groups <- c(comp[[2]], comp[[3]])
  
  # 1. Grab isoform counts
  sub_counts <- AggregateExpression(obj,
                                    group.by = 'pseudobulk_id',
                                    assays = 'ISO',
                                    slot = 'counts',
                                    return.seurat = FALSE)$ISO
  
  # Ensure sub_counts is a standard dense matrix, since DEXSeq strictly requires it
  sub_counts <- as.matrix(sub_counts)
  
  # 2. Format sample metadata
  sub_meta <- obj@meta.data %>%
    dplyr::select(donor_id, analysis_group, pseudobulk_id) %>%
    tidyr::drop_na() %>%
    dplyr::distinct(pseudobulk_id, .keep_all = TRUE) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames('pseudobulk_id')
  
  rownames(sub_meta) <- gsub('[^[:alnum:]]', '', rownames(sub_meta))
  colnames(sub_counts) <- gsub('[^[:alnum:]]', '', colnames(sub_counts))
  
  common_ids <- intersect(rownames(sub_meta), colnames(sub_counts))
  sub_meta <- sub_meta[common_ids, ]
  sub_counts <- sub_counts[, common_ids]
  
  sub_meta$donor_id <- as.factor(sub_meta$donor_id)
  sub_meta$analysis_group <- factor(sub_meta$analysis_group, levels = groups)
  
  # 3. Get transcript -> gene groupings from our SQANTI metadata
  iso_meta <- obj[["ISO"]][[]]
  iso_meta <- iso_meta[rownames(sub_counts), ]
  
  valid_features <- rownames(iso_meta)[!is.na(iso_meta$associated_gene)]
  
  # Crucially round numeric matrices explicitly 
  sub_counts <- round(sub_counts[valid_features, ])
  iso_meta <- iso_meta[valid_features, ]
  
  feature_ids <- rownames(sub_counts)
  group_ids <- iso_meta$associated_gene
  
  # 4. Construct DEXSeq core objects
  dxd <- DEXSeqDataSet(countData = sub_counts,
                       sampleData = sub_meta,
                       design = ~ sample + exon + donor_id:exon + analysis_group:exon,
                       featureID = feature_ids,
                       groupID = group_ids)
  
  # Filter 1: Expressed in at least 3 samples
  dxd <- dxd[rowSums(featureCounts(dxd) >= 10) >= 3, ]
  
  # Filter 2: Remove genes with only a single isoform (no alternatives to use)
  genes_with_multiple_isoforms <- names(which(table(groupIDs(dxd)) > 1))
  dxd <- dxd[groupIDs(dxd) %in% genes_with_multiple_isoforms, ]
  
  # Filter 3: Fallback check
  if (nrow(dxd) == 0) {
    warning("No valid multiple-isoform genes survived filtering for ", comp[[1]])
    return(NULL)
  }
  
  # 5. Run standard DEXSeq pipeline
  message("   Estimating size factors & dispersion...")
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet = TRUE)
  
  message("   Testing for differential usage...")
  dxd <- testForDEU(dxd, reducedModel = ~ sample + exon + donor_id:exon)
  
  message("   Estimating Exon (Isoform) fold changes...")
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "analysis_group")
  
  res <- as.data.frame(DEXSeqResults(dxd))
  res <- res[order(res$padj), ]
  
  # Optional: Append the associated gene symbol directly into the results if you mapped it
  if("associated_gene_symbol" %in% colnames(iso_meta)) {
    symb_vector <- iso_meta$associated_gene_symbol
    names(symb_vector) <- rownames(iso_meta)
    res$gene_symbol <- symb_vector[res$featureID]
  }
  
  return(res)
}

dexseq_markers <- setNames(
  lapply(comparisons, function(comp) {
    groups <- c(comp[[2]], comp[[3]])
    
    cells_use <- total_filtered@meta.data %>%
      filter(analysis_group %in% groups) %>%
      pull(cell_barcode)
    
    obj <- subset(total_filtered, cells = cells_use)
    run_dexseq_comparison(comp, obj)
  }),
  sapply(comparisons, `[[`, 1)
)


markers <- dexseq_markers$adaptive_v_dim %>%
  mutate(expression = case_when(
    log2fold_CD56dim.NK_Adaptive.NK > 0.5 & padj < 0.05 ~ 'Adaptive NK',
    log2fold_CD56dim.NK_Adaptive.NK < -0.5 & padj < 0.05 ~ 'CD56dim NK',
    TRUE ~ 'Not significant'
  )) %>%
  mutate(gene = rownames(dexseq_markers$adaptive_v_dim)) %>%
  mutate(expression = factor(expression, levels = c('Adaptive NK', 'CD56dim NK', 'Not significant'))) %>%
  filter(!is.na(padj))

ggplot(markers, aes(x = log2fold_CD56dim.NK_Adaptive.NK, y = -log10(padj), colour = expression)) +
  geom_vline(xintercept = 0, linetype = 'dashed', colour = 'black') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', colour = 'grey50', alpha = 0.5) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dashed', colour = 'grey50', alpha = 0.5) +
  geom_point(alpha = 1, size = 1) +
  theme(
    legend.position = 'none',
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
  ) +
  scale_colour_manual(values = c('red', 'blue', 'black')) +
  labs(
    title = NULL,
    x = expression("Adaptive NK v CD56dim NK (log"[2]*"(fold change))"),
    y = expression("-log"[10]*"(FDR)"),
    colour = NULL
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# 1. Helper Function: Map a single Gene Symbol to ENSG ID using the GTF
map_symbol_to_ensg <- function(symbol, gtf_path) {
  # Fast extraction using awk to find exactly the gene_name match
  cmd <- sprintf("awk -F'\t' '$3 == \"gene\" && $9 ~ /gene_name \"%s\"/' %s", symbol, gtf_path)
  
  gtf <- tryCatch(
    fread(cmd = cmd, header = FALSE, sep = "\t"),
    error = function(e) return(NULL)
  )
  
  if (nrow(gtf) == 0) stop(paste("Gene symbol", symbol, "not found in GTF."))
  
  # Extract gene_id from the first match
  ensg <- sub('.*gene_id "([^"]+)".*', '\\1', gtf$V9[1])
  return(ensg)
}

# 2. Main Plotting Function
plot_isoform_distribution <- function(seurat_obj, gene_symbol, gtf_path, group_by = "seurat_clusters") {
  
  # Map symbol to ENSG ID
  message(paste("Mapping", gene_symbol, "to Ensembl ID..."))
  ensg_id <- map_symbol_to_ensg(gene_symbol, gtf_path)
  ensg_base <- sub("\\.[0-9]+$", "", ensg_id)
  message(paste("Found ENSG:", ensg_id))
  
  # Find transcripts for this gene in the Seurat object's ISO assay
  iso_meta <- seurat_obj[["ISO"]][[]]
  
  # Strip versions for safe matching
  iso_meta$associated_gene_base <- sub("\\.[0-9]+$", "", iso_meta$associated_gene)
  target_transcripts <- rownames(iso_meta)[iso_meta$associated_gene_base == ensg_base]
  
  if (length(target_transcripts) == 0) {
    stop(paste("No isoforms found for gene", gene_symbol, "in the Seurat object."))
  }
  
  message(paste("Found", length(target_transcripts), "isoforms for", gene_symbol))
  
  # Aggregate counts across the specified populations
  pb_counts <- AggregateExpression(seurat_obj, 
                                   group.by = group_by, 
                                   assays = 'ISO', 
                                   features = target_transcripts,
                                   slot = 'counts',
                                   return.seurat = FALSE)$ISO
  
  # Convert matrix to a long-format data.frame for ggplot2
  pb_df <- as.data.frame(as.matrix(pb_counts))
  pb_df$Transcript <- rownames(pb_df)
  
  df_long <- pivot_longer(pb_df, 
                          cols = -Transcript, 
                          names_to = "Population", 
                          values_to = "Counts")
  
  # Convert Population to factor to keep ordering natural rather than alphabetical (if numeric-like)
  # df_long$Population <- factor(df_long$Population, levels = unique(df_long$Population))
  
  # Generate ggplot (No theme() or titles, per your preference for easy downstream customization)
  p <- ggplot(df_long, aes(x = Population, y = Counts, fill = Transcript)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = scales::percent_format()) +
    ylab("Isoform Proportion") +
    xlab(NULL)
  
  return(p)
}

total$celltype_l5 <- factor(total$celltype_l5, levels = rev(c(
  'Naïve-like ILC', 'CD56bright NK', 'CD56dim NK', 'Cycling NK',
  'Adaptive NK 1', 'Adaptive NK 2', 'Adaptive NK 3', 'Adaptive NK 4', 'Adaptive NK 5', 'Adaptive NK 6'
)))

plot_isoform_distribution(seurat_obj = total, 
                          gene_symbol = "KIR2DL1", 
                          gtf_path = gtf_file,
                          group_by = "celltype_l5")
