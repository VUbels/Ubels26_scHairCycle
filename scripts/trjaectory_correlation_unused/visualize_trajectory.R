#!/usr/bin/env Rscript

#################################################################
# GRANULAR TRAJECTORY HEATMAPS
# Side-by-side: scRNA | Chromatin Accessibility | TF Activity
# Using rolling mean smoothing like the spatial trajectory script
#################################################################

library(ArchR)
library(monocle3)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
library(Matrix)
library(zoo)
library(circlize)

set.seed(42)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

#################################################################
# PARAMETERS
#################################################################

main_folder <- "./"
scatac_dir <- paste0(main_folder, "scatac_files/")
spatial_dir <- paste0(main_folder, "spatial_trajectory/")
output_dir <- paste0(main_folder, "integrated_trajectory/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

proj_save_folder <- paste0(scatac_dir, "archr_proj/")
krt_clusters <- c("C1", "C2", "C4", "C5", "C6")

# Smoothing parameters (same as spatial script)
smooth_window <- 50  # Rolling mean window
min_expr <- 0.1      # Minimum expression filter
zscore_cap <- 3      # Cap z-scores

# Priority genes to label
priority_genes <- c("KRT15", "KRT14", "KRT5", "SOX9", "LHX2", "CD34", 
                    "COMP", "LGR5", "KRT19", "DKK3", "NFATC1")

#################################################################
# LOAD DATA
#################################################################

message("=== LOADING DATA ===")

skin_archr <- loadArchRProject(path = proj_save_folder)
cds_visium <- readRDS(paste0(spatial_dir, "monocle3_cds.rds"))
gene_modules <- read.csv(paste0(spatial_dir, "gene_modules.csv"))

message(paste("scATAC cells:", nCells(skin_archr)))
message(paste("Visium cells:", ncol(cds_visium)))

#################################################################
# PART 1: PREPARE VISIUM DATA (scRNA)
#################################################################

message("\n=== PART 1: PREPARING VISIUM scRNA DATA ===")

# Get expression matrix and pseudotime
visium_expr <- as.matrix(exprs(cds_visium))
visium_gene_names <- rowData(cds_visium)$gene_short_name
if (!is.null(visium_gene_names)) {
  rownames(visium_expr) <- visium_gene_names
}

visium_pt <- pseudotime(cds_visium)

# Order cells by pseudotime
visium_pt_order <- order(visium_pt)
visium_expr_ordered <- visium_expr[, visium_pt_order]

message(paste("Visium cells ordered:", ncol(visium_expr_ordered)))

#################################################################
# PART 2: PREPARE scATAC DATA
#################################################################

message("\n=== PART 2: PREPARING scATAC DATA ===")

# Subset to keratinocytes
clusters <- getCellColData(skin_archr)$Clusters
krt_cells <- getCellNames(skin_archr)[clusters %in% krt_clusters]
skin_archr_krt <- subsetCells(skin_archr, cellNames = krt_cells)

message(paste("Keratinocyte cells:", nCells(skin_archr_krt)))

# Get gene score matrix
gsm_se <- getMatrixFromProject(skin_archr_krt, useMatrix = "GeneScoreMatrix")
gs_matrix <- assay(gsm_se)
gene_names_atac <- gsm_se@elementMetadata@listData[["name"]]
rownames(gs_matrix) <- gene_names_atac

# Get motif matrix
motif_se <- getMatrixFromProject(skin_archr_krt, useMatrix = "MotifMatrix")
motif_matrix <- assay(motif_se)

message(paste("Gene score matrix:", nrow(gs_matrix), "x", ncol(gs_matrix)))
message(paste("Motif matrix:", nrow(motif_matrix), "x", ncol(motif_matrix)))

#################################################################
# PART 3: COMPUTE MODULE-BASED PSEUDOTIME FOR scATAC
#################################################################

message("\n=== PART 3: COMPUTING MODULE-BASED PSEUDOTIME ===")

# Build module gene lists
module_list <- split(gene_modules$id, gene_modules$module)
names(module_list) <- paste0("Module_", names(module_list))

# Compute Visium module peak positions
visium_module_scores <- sapply(names(module_list), function(mod) {
  genes <- module_list[[mod]]
  genes_present <- intersect(genes, rownames(visium_expr))
  if (length(genes_present) > 0) {
    colMeans(visium_expr[genes_present, , drop = FALSE], na.rm = TRUE)
  } else {
    rep(NA, ncol(visium_expr))
  }
})

module_peak_pt <- sapply(colnames(visium_module_scores), function(mod) {
  scores <- visium_module_scores[, mod]
  if (all(is.na(scores))) return(NA)
  weights <- scores - min(scores, na.rm = TRUE)
  weights <- weights / sum(weights, na.rm = TRUE)
  sum(weights * visium_pt, na.rm = TRUE)
})

# Compute scATAC module scores
atac_module_scores <- sapply(names(module_list), function(mod) {
  genes <- module_list[[mod]]
  genes_present <- intersect(genes, rownames(gs_matrix))
  if (length(genes_present) > 0) {
    colMeans(gs_matrix[genes_present, , drop = FALSE], na.rm = TRUE)
  } else {
    rep(NA, ncol(gs_matrix))
  }
})
rownames(atac_module_scores) <- colnames(gs_matrix)

# Assign each scATAC cell to best module
atac_module_scaled <- scale(atac_module_scores)
best_module <- apply(atac_module_scaled, 1, function(x) {
  if (all(is.na(x))) return(NA)
  colnames(atac_module_scaled)[which.max(x)]
})

# Derive pseudotime from module assignment
atac_pseudotime <- module_peak_pt[best_module]
names(atac_pseudotime) <- names(best_module)

# Remove NA cells
valid_atac_cells <- names(atac_pseudotime)[!is.na(atac_pseudotime)]
atac_pseudotime <- atac_pseudotime[valid_atac_cells]

message(paste("scATAC cells with valid pseudotime:", length(valid_atac_cells)))

# Order scATAC cells by pseudotime
atac_pt_order <- names(sort(atac_pseudotime))
gs_matrix_ordered <- gs_matrix[, atac_pt_order]
motif_matrix_ordered <- motif_matrix[, atac_pt_order]

message(paste("scATAC cells ordered:", ncol(gs_matrix_ordered)))

#################################################################
# PART 4: GET COMMON TRAJECTORY GENES
#################################################################

message("\n=== PART 4: IDENTIFYING TRAJECTORY GENES ===")

# Use genes from the Monocle modules (trajectory-associated)
trajectory_genes <- unique(gene_modules$id)

# Find genes present in both datasets
common_genes <- intersect(trajectory_genes, rownames(visium_expr))
common_genes <- intersect(common_genes, rownames(gs_matrix))

message(paste("Trajectory genes:", length(trajectory_genes)))
message(paste("Common genes in both datasets:", length(common_genes)))

# Add priority genes if present
priority_in_data <- priority_genes[priority_genes %in% common_genes]
message(paste("Priority genes found:", paste(priority_in_data, collapse = ", ")))

#################################################################
# PART 5: SMOOTH AND SCALE - VISIUM scRNA
#################################################################

message("\n=== PART 5: SMOOTHING VISIUM scRNA ===")

# Subset to common genes
visium_traj <- visium_expr_ordered[common_genes, ]

# Apply rolling mean smoothing
smooth_window_vis <- min(smooth_window, floor(ncol(visium_traj) / 10))
message(paste("Visium smoothing window:", smooth_window_vis))

visium_smoothed <- t(apply(visium_traj, 1, function(x) {
  zoo::rollmean(x, k = smooth_window_vis, fill = NA, align = "center")
}))

# Remove edge NAs
valid_cols_vis <- complete.cases(t(visium_smoothed))
visium_smoothed <- visium_smoothed[, valid_cols_vis]

# Filter lowly expressed genes (but keep priority genes)
gene_means_vis <- rowMeans(visium_smoothed, na.rm = TRUE)
genes_keep_vis <- gene_means_vis > min_expr | rownames(visium_smoothed) %in% priority_genes
visium_smoothed <- visium_smoothed[genes_keep_vis, ]

# Z-scale
visium_scaled <- t(scale(t(visium_smoothed)))
visium_scaled[is.na(visium_scaled)] <- 0
visium_scaled[visium_scaled > zscore_cap] <- zscore_cap
visium_scaled[visium_scaled < -zscore_cap] <- -zscore_cap

# Order by peak position
peak_pos_vis <- apply(visium_scaled, 1, which.max)
gene_order_vis <- names(sort(peak_pos_vis))
visium_final <- visium_scaled[gene_order_vis, ]

message(paste("Visium genes after filtering:", nrow(visium_final)))
message(paste("Visium pseudotime bins:", ncol(visium_final)))

#################################################################
# PART 6: SMOOTH AND SCALE - scATAC GENE SCORES
#################################################################

message("\n=== PART 6: SMOOTHING scATAC CHROMATIN ===")

# Subset to genes in Visium final (for matching rows)
genes_for_atac <- intersect(gene_order_vis, rownames(gs_matrix_ordered))
atac_traj <- gs_matrix_ordered[genes_for_atac, ]

# Apply rolling mean smoothing
smooth_window_atac <- min(smooth_window, floor(ncol(atac_traj) / 10))
message(paste("scATAC smoothing window:", smooth_window_atac))

atac_smoothed <- t(apply(atac_traj, 1, function(x) {
  zoo::rollmean(x, k = smooth_window_atac, fill = NA, align = "center")
}))

# Remove edge NAs
valid_cols_atac <- complete.cases(t(atac_smoothed))
atac_smoothed <- atac_smoothed[, valid_cols_atac]

# Z-scale
atac_scaled <- t(scale(t(atac_smoothed)))
atac_scaled[is.na(atac_scaled)] <- 0
atac_scaled[atac_scaled > zscore_cap] <- zscore_cap
atac_scaled[atac_scaled < -zscore_cap] <- -zscore_cap

# Order genes to match Visium (same gene order)
genes_both <- intersect(gene_order_vis, rownames(atac_scaled))
atac_final <- atac_scaled[genes_both, ]

# Also subset Visium to matching genes
visium_final <- visium_final[genes_both, ]

message(paste("Matched genes:", nrow(atac_final)))
message(paste("scATAC pseudotime bins:", ncol(atac_final)))

#################################################################
# PART 7: SMOOTH AND SCALE - TF ACTIVITY
#################################################################

message("\n=== PART 7: SMOOTHING TF ACTIVITY ===")

# Apply rolling mean smoothing to motif matrix
motif_smoothed <- t(apply(motif_matrix_ordered, 1, function(x) {
  zoo::rollmean(x, k = smooth_window_atac, fill = NA, align = "center")
}))

# Remove edge NAs
valid_cols_motif <- complete.cases(t(motif_smoothed))
motif_smoothed <- motif_smoothed[, valid_cols_motif]

# Filter to variable TFs
motif_var <- apply(motif_smoothed, 1, var, na.rm = TRUE)
top_tfs <- names(sort(motif_var, decreasing = TRUE))[1:50]
motif_smoothed <- motif_smoothed[top_tfs, ]

# Z-scale
motif_scaled <- t(scale(t(motif_smoothed)))
motif_scaled[is.na(motif_scaled)] <- 0
motif_scaled[motif_scaled > zscore_cap] <- zscore_cap
motif_scaled[motif_scaled < -zscore_cap] <- -zscore_cap

# Order by peak position
peak_pos_motif <- apply(motif_scaled, 1, which.max)
motif_order <- names(sort(peak_pos_motif))
motif_final <- motif_scaled[motif_order, ]

# Clean TF names for display
clean_tf_names <- gsub("_.*", "", rownames(motif_final))
# Handle duplicates from cleaning
clean_tf_names <- make.unique(clean_tf_names)

message(paste("TFs in heatmap:", nrow(motif_final)))
message(paste("TF pseudotime bins:", ncol(motif_final)))

#################################################################
# PART 8: CREATE GENE ANNOTATIONS
#################################################################

message("\n=== PART 8: CREATING ANNOTATIONS ===")

# Find priority genes in final matrix
priority_in_final <- priority_genes[priority_genes %in% rownames(visium_final)]
priority_positions <- match(priority_in_final, rownames(visium_final))

# Also add top genes at different pseudotime positions
n_label_bins <- 10
gene_peak_bins <- apply(visium_final, 1, which.max)
n_cols_vis <- ncol(visium_final)
col_bins <- cut(seq_len(n_cols_vis), breaks = n_label_bins, labels = FALSE)

bin_top_genes <- character(0)
bin_top_positions <- integer(0)

for (b in seq_len(n_label_bins)) {
  bin_cols <- which(col_bins == b)
  genes_peaking_here <- names(gene_peak_bins)[gene_peak_bins %in% bin_cols]
  genes_peaking_here <- setdiff(genes_peaking_here, c(priority_in_final, bin_top_genes))
  
  if (length(genes_peaking_here) > 0) {
    peak_vals <- sapply(genes_peaking_here, function(g) {
      max(visium_final[g, bin_cols])
    })
    top_gene <- genes_peaking_here[which.max(peak_vals)]
    bin_top_genes <- c(bin_top_genes, top_gene)
    bin_top_positions <- c(bin_top_positions, match(top_gene, rownames(visium_final)))
  }
}

# Combine labels
all_label_genes <- c(priority_in_final, bin_top_genes)
all_label_positions <- c(priority_positions, bin_top_positions)

# Remove NAs
valid_labels <- !is.na(all_label_positions)
all_label_genes <- all_label_genes[valid_labels]
all_label_positions <- all_label_positions[valid_labels]

message(paste("Genes to label:", length(all_label_genes)))

#################################################################
# PART 9: BUILD HEATMAPS
#################################################################

message("\n=== PART 9: BUILDING HEATMAPS ===")

# Gene annotation (between scRNA and scATAC)
gene_anno <- rowAnnotation(
  link = anno_mark(
    at = all_label_positions,
    labels = all_label_genes,
    labels_gp = gpar(fontsize = 6),
    link_width = unit(3, "mm")
  )
)

# TF annotation
tf_label_idx <- seq(1, nrow(motif_final), length.out = min(20, nrow(motif_final)))
tf_label_idx <- round(tf_label_idx)

tf_anno <- rowAnnotation(
  link = anno_mark(
    at = tf_label_idx,
    labels = clean_tf_names[tf_label_idx],
    labels_gp = gpar(fontsize = 6),
    link_width = unit(3, "mm")
  )
)

# Heatmap 1: Visium scRNA
ht_rna <- Heatmap(
  visium_final,
  name = "scRNA\nZ-score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = colorRamp2(c(-zscore_cap, 0, zscore_cap), c("purple", "black", "yellow")),
  column_title = "scRNA Expression\n(Visium Pseudotime →)",
  row_title = paste0("Genes (n=", nrow(visium_final), ")"),
  use_raster = TRUE,
  raster_quality = 2,
  width = unit(5, "cm"),
  right_annotation = gene_anno
)

# Heatmap 2: scATAC Chromatin
ht_atac <- Heatmap(
  atac_final,
  name = "Chromatin\nZ-score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = colorRamp2(c(-zscore_cap, 0, zscore_cap), c("purple", "black", "yellow")),
  column_title = "Chromatin Accessibility\n(scATAC Pseudotime →)",
  use_raster = TRUE,
  raster_quality = 2,
  width = unit(5, "cm")
)

# Heatmap 3: TF Activity
ht_tf <- Heatmap(
  motif_final,
  name = "TF\nZ-score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = colorRamp2(c(-zscore_cap, 0, zscore_cap), c("blue", "white", "red")),
  column_title = "TF Motif Activity\n(scATAC Pseudotime →)",
  row_title = paste0("TFs (n=", nrow(motif_final), ")"),
  use_raster = TRUE,
  raster_quality = 2,
  width = unit(5, "cm"),
  right_annotation = tf_anno
)

#################################################################
# PART 10: SAVE COMBINED HEATMAP
#################################################################

message("\n=== PART 10: SAVING HEATMAPS ===")

# Combined: RNA + Chromatin side by side, TF below
library(grid)

pdf(file.path(output_dir, "granular_trajectory_heatmaps.pdf"), width = 14, height = 16)

pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(1.5, 1), "null"))))

# Top: scRNA + Chromatin
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_rna + ht_atac,
     column_title = "Gene Expression & Chromatin Accessibility Along HF Differentiation",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     ht_gap = unit(0.5, "cm"),
     newpage = FALSE)
popViewport()

# Bottom: TF activity
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht_tf,
     column_title = "TF Motif Activity Along Trajectory",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     newpage = FALSE)
popViewport()

popViewport()
dev.off()

png(file.path(output_dir, "granular_trajectory_heatmaps.png"), 
    width = 14, height = 16, units = "in", res = 300)

pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(1.5, 1), "null"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_rna + ht_atac,
     column_title = "Gene Expression & Chromatin Accessibility Along HF Differentiation",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     ht_gap = unit(0.5, "cm"),
     newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht_tf,
     column_title = "TF Motif Activity Along Trajectory",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     newpage = FALSE)
popViewport()

popViewport()
dev.off()

message("Granular heatmaps saved!")

#################################################################
# PART 11: COMPUTE GENE-LEVEL CORRELATIONS
#################################################################

message("\n=== PART 11: COMPUTING CORRELATIONS ===")

# For each gene, correlate scRNA pattern with scATAC pattern
# Need to resample to same number of columns first
n_bins <- 100

# Resample Visium to n_bins
visium_resampled <- t(apply(visium_final, 1, function(x) {
  approx(x, n = n_bins)$y
}))

# Resample scATAC to n_bins
atac_resampled <- t(apply(atac_final, 1, function(x) {
  approx(x, n = n_bins)$y
}))

# Compute correlation per gene
gene_cors <- sapply(rownames(visium_resampled), function(gene) {
  if (gene %in% rownames(atac_resampled)) {
    cor(visium_resampled[gene, ], atac_resampled[gene, ], use = "complete.obs")
  } else {
    NA
  }
})

cor_summary <- data.frame(
  gene = names(gene_cors),
  rna_atac_correlation = gene_cors
) %>% 
  filter(!is.na(rna_atac_correlation)) %>%
  arrange(desc(rna_atac_correlation))

write.csv(cor_summary, file.path(output_dir, "gene_rna_atac_correlations.csv"), row.names = FALSE)

message("\nGene-level RNA-ATAC correlations:")
message(paste("Mean:", round(mean(gene_cors, na.rm = TRUE), 3)))
message(paste("Median:", round(median(gene_cors, na.rm = TRUE), 3)))
message(paste("Positive correlations:", sum(gene_cors > 0, na.rm = TRUE), "/", sum(!is.na(gene_cors))))

# Top correlated genes
message("\nTop 10 correlated genes:")
print(head(cor_summary, 10))

# Priority gene correlations
message("\nPriority gene correlations:")
priority_cors <- cor_summary %>% filter(gene %in% priority_genes)
print(priority_cors)

#################################################################
# PART 12: SAVE DATA
#################################################################

message("\n=== PART 12: SAVING DATA ===")

saveRDS(list(
  visium_smoothed = visium_final,
  atac_smoothed = atac_final,
  motif_smoothed = motif_final,
  gene_order = rownames(visium_final),
  tf_order = rownames(motif_final),
  gene_correlations = cor_summary
), file.path(output_dir, "granular_trajectory_data.rds"))

#################################################################
# SUMMARY
#################################################################

message("\n")
message("================================================================")
message("           GRANULAR TRAJECTORY ANALYSIS COMPLETE                ")
message("================================================================")
message(paste("Genes in heatmap:", nrow(visium_final)))
message(paste("TFs in heatmap:", nrow(motif_final)))
message(paste("Visium pseudotime bins:", ncol(visium_final)))
message(paste("scATAC pseudotime bins:", ncol(atac_final)))
message(paste("Mean RNA-ATAC correlation:", round(mean(gene_cors, na.rm = TRUE), 3)))
message("")
message("Key outputs:")
message(paste("  -", file.path(output_dir, "granular_trajectory_heatmaps.png")))
message(paste("  -", file.path(output_dir, "gene_rna_atac_correlations.csv")))
message("")
message("INTERPRETATION:")
message("  - Left heatmap: scRNA expression along Visium pseudotime")
message("  - Right heatmap: Chromatin accessibility along scATAC pseudotime")
message("  - Genes are in SAME ORDER (by scRNA peak position)")
message("  - If patterns match → chromatin dynamics follow expression")
message("  - Diagonal pattern in both = successful integration")
message("================================================================")