#!/usr/bin/env Rscript

#################################################################
# Description: Spatial trajectory analysis of KRT15+ cells
#################################################################

# Following script performs spatial trajectory analysis along the
# hair follicle axis using KRT15+ cells from Visium HD data.
# Genes are ordered by peak expression position along the trajectory.

#################################################################
# LIBRARY LOADING
#################################################################

library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(viridis)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#################################################################
# SETUP PROJECT PARAMETERS
#################################################################

project <- "ubels26_haircycle"
main_folder <- "./"
output_dir <- paste0(main_folder, "spatial_trajectory/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

vis_obj <- readRDS(paste0(main_folder, "Spatial_scalp_S2_final.rds"))
vis_obj <- UpdateSeuratObject(vis_obj)

#################################################################
# PRIORITY GENES FOR LABELING
#################################################################

priority_genes <- c("CD34", "LGR6", "SOX9", "LHX2", "LGR5", "NFATC1")

#################################################################
# ANALYSIS PARAMETERS
#################################################################

krt15_threshold <- 0.5
n_position_bins <- 50
n_variable_features <- 3000
cor_threshold <- 0.3
fc_threshold <- 2
zscore_cap <- 2

#################################################################
# SUBSET TO KRT15+ CELLS
#################################################################

blacklist_genes <- get_blacklist_genes(vis_obj)

krt15_expr <- FetchData(vis_obj, vars = "KRT15")
krt15_positive <- rownames(krt15_expr)[krt15_expr$KRT15 > krt15_threshold]
subset_obj <- subset(vis_obj, cells = krt15_positive)

message(paste("KRT15+ cells:", length(krt15_positive)))

#################################################################
# ORDER CELLS BY SPATIAL COORDINATE (TRAJECTORY)
#################################################################

coords <- GetTissueCoordinates(subset_obj)
coords <- coords %>% arrange(desc(x))
ordered_cells <- coords$cell

#################################################################
# IDENTIFY VARIABLE GENES AND REMOVE BLACKLIST
#################################################################

subset_obj <- FindVariableFeatures(subset_obj, nfeatures = n_variable_features)
all_genes <- VariableFeatures(subset_obj)
all_genes <- all_genes[!all_genes %in% blacklist_genes]

message(paste("Variable genes after blacklist removal:", length(all_genes)))

#################################################################
# BUILD SMOOTHED EXPRESSION MATRIX
#################################################################

expr_matrix <- FetchData(subset_obj, vars = all_genes)
expr_matrix <- expr_matrix[ordered_cells, ]

cell_bins <- cut(1:nrow(expr_matrix), breaks = n_position_bins, labels = FALSE)

expr_smoothed <- sapply(1:n_position_bins, function(b) {
  colMeans(expr_matrix[cell_bins == b, , drop = FALSE], na.rm = TRUE)
})
expr_smoothed <- t(expr_smoothed)
colnames(expr_smoothed) <- all_genes

#################################################################
# FILTER GENES BY SPATIAL PATTERN
#################################################################
position_vec <- 1:n_position_bins

spatial_cor <- apply(expr_smoothed, 2, function(x) {
  abs(cor(x, position_vec, method = "spearman"))
})

fold_change <- apply(expr_smoothed, 2, function(x) {
  x_pos <- x[x > 0]
  if (length(x_pos) < 2) return(0)
  max(x_pos) / (min(x_pos) + 0.001)
})

gene_stats <- data.frame(
  gene = colnames(expr_smoothed),
  cor = spatial_cor,
  fc = fold_change
)

variable_genes <- gene_stats$gene[gene_stats$cor > cor_threshold | gene_stats$fc > fc_threshold]
variable_genes <- variable_genes[variable_genes %in% colnames(expr_smoothed)]

message(paste("Genes passing spatial pattern filter:", length(variable_genes)))

#################################################################
# Z-SCALE AND CAP EXPRESSION VALUES
#################################################################

expr_filtered <- expr_smoothed[, variable_genes, drop = FALSE]
expr_scaled <- t(scale(t(expr_filtered)))

expr_scaled[expr_scaled > zscore_cap] <- zscore_cap
expr_scaled[expr_scaled < -zscore_cap] <- -zscore_cap

#################################################################
# ORDER GENES BY PEAK POSITION
#################################################################

peak_position <- apply(expr_scaled, 2, which.max)
gene_order <- names(sort(peak_position))
expr_scaled <- expr_scaled[, gene_order]

#################################################################
# ASSIGN GENES TO SPATIAL BINS
#################################################################

priority_present <- priority_genes[priority_genes %in% gene_order]
message("Priority genes found in data:")
print(priority_present)

bin_labels <- rep("", n_position_bins)
used_genes <- c()

# Get peak bin for each gene
gene_peak_bins <- apply(expr_scaled, 2, which.max)

# Assign priority genes first to their peak bins
for (gene in priority_present) {
  peak_bin <- gene_peak_bins[gene]
  if (bin_labels[peak_bin] == "") {
    bin_labels[peak_bin] <- gene
    used_genes <- c(used_genes, gene)
  }
}

# Fill remaining bins - ONLY consider genes that actually PEAK at that bin
for (bin in 1:n_position_bins) {
  if (bin_labels[bin] == "") {
    # Get genes that peak at this specific bin
    genes_peaking_here <- names(gene_peak_bins)[gene_peak_bins == bin]
    # Exclude already used genes
    available <- genes_peaking_here[!genes_peaking_here %in% used_genes]
    
    if (length(available) > 0) {
      # Among genes peaking here, pick highest expression at peak
      peak_expr <- sapply(available, function(g) expr_scaled[bin, g])
      top_gene <- available[which.max(peak_expr)]  # Use index into available, not names()
      
      if (expr_scaled[bin, top_gene] > 0.5) {
        bin_labels[bin] <- top_gene
        used_genes <- c(used_genes, top_gene)
      }
    }
  }
}

# Store labeled genes
labeled_genes <- unique(bin_labels[bin_labels != ""])

message(paste("Genes labeled on both axes:", length(labeled_genes)))
print(labeled_genes)

#################################################################
# BUILD ANNOTATIONS FOR HEATMAP
#################################################################

# Column annotation (top): genes at their spatial bin position
col_anno <- HeatmapAnnotation(
  gene = anno_mark(
    at = which(bin_labels != ""),
    labels = bin_labels[bin_labels != ""],
    side = "top",
    labels_gp = gpar(fontsize = 7),
    link_gp = gpar(lwd = 0.5),
    link_height = unit(3, "mm")
  )
)

# Row annotation (right): SAME genes at their gene order position
label_indices_row <- which(gene_order %in% labeled_genes)

row_anno <- rowAnnotation(
  link = anno_mark(
    at = label_indices_row,
    labels = gene_order[label_indices_row],
    labels_gp = gpar(fontsize = 7),
    link_width = unit(5, "mm")
  )
)

#################################################################
# GENERATE HEATMAP
#################################################################

ht <- Heatmap(
  t(expr_scaled),
  name = "Z-score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = inferno(100),
  top_annotation = col_anno,
  right_annotation = row_anno,
  column_title = "Position along trajectory →",
  row_title = paste0("Genes (n=", length(gene_order), ")"),
  use_raster = TRUE
)

#################################################################
# SAVE OUTPUT
#################################################################

pdf(paste0(output_dir, "spatial_trajectory_heatmap.pdf"), width = 6, height = 8)
draw(ht)
dev.off()

png(paste0(output_dir, "spatial_trajectory_heatmap.png"), width = 6, height = 8, units = "in", res = 300)
draw(ht)
dev.off()

message(paste("Heatmap saved to:", output_dir))
message("Analysis complete.")

# Set expression below threshold to NA for visualization
vis_obj$KRT15_filtered <- ifelse(
  FetchData(vis_obj, vars = "KRT15")$KRT15 > 0.5,
  FetchData(vis_obj, vars = "KRT15")$KRT15,
  NA
)
SpatialFeaturePlot(vis_obj, features = "KRT15_filtered", pt.size.factor = 5)
SpatialFeaturePlot(vis_obj, features = "KRT79", pt.size.factor = 5)
SpatialFeaturePlot(vis_obj, features = "AQP1", pt.size.factor = 5)

