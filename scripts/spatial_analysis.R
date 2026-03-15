#!/usr/bin/env Rscript

#################################################################
# Spatial trajectory with Monocle3 pseudotime
# Using graph_test + find_gene_modules
#################################################################

library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
library(leidenbase)

source("./scripts/helper_functions.R")

#################################################################
# SETUP
#################################################################

main_folder <- "./"
output_dir <- paste0(main_folder, "spatial_trajectory/")
marker_dir <- paste0(output_dir, "marker_genes/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)

vis_obj <- readRDS(paste0(main_folder, "Spatial_scalp_S2_final.rds"))
vis_obj <- UpdateSeuratObject(vis_obj)

#################################################################
# PARAMETERS
#################################################################

krt15_threshold <- 0.5
min_expr <- 0.1
zscore_cap <- 3
module_resolution <- 1e-2
priority_genes <- c("CD34", "LGR5", "CD200", "COMP", "KRT19", "DKK3", "TGFB2", "CD71", "DIO2", "ANGPTL2", "WIF1", "LGR6", "SOX9", "LHX2", "LGR5", "NFATC1", "KRT15")

#################################################################
# SELECT LINEAGE CELLS
#################################################################

selected_cells <- spatial_lasso_selector(vis_obj, "KRT15")

vis_obj <- add_lineage_label(
  obj = vis_obj,
  selected_cells = selected_cells,
  feature = "KRT15",
  threshold = krt15_threshold
)

vis_obj$KRT15_Lineage <- factor(
  vis_obj$KRT15_Lineage,
  levels = c(FALSE, TRUE),
  labels = c("Other", "KRT15_Lineage")
)

#################################################################
# VISUALIZE KRT15 EXPRESSION AND LINEAGE SELECTION
#################################################################

# 1. SpatialFeaturePlot: All KRT15+ cells (before filtering)
# Show raw KRT15 expression across entire tissue
p_krt15_all <- SpatialFeaturePlot(
  vis_obj, 
  features = "KRT15", 
  pt.size.factor = 3,
  image.alpha = 0
) + 
  scale_fill_viridis(option = "B") +
  ggtitle("KRT15 Expression (All Cells)")

print(p_krt15_all)
ggsave(
  paste0(output_dir, "KRT15_expression_all_cells.png"), 
  p_krt15_all, 
  width = 10, 
  height = 8
)

# 2. SpatialFeaturePlot: KRT15 expression filtered (above threshold only)
vis_obj$KRT15_filtered <- ifelse(
  FetchData(vis_obj, vars = "KRT15")$KRT15 > krt15_threshold,
  FetchData(vis_obj, vars = "KRT15")$KRT15,
  NA
)

p_krt15_filtered <- SpatialFeaturePlot(
  vis_obj, 
  features = "KRT15_filtered", 
  pt.size.factor = 3,
  image.alpha = 0
) + 
  scale_fill_viridis(option = "B", na.value = "grey40") +
  ggtitle(paste0("KRT15 Expression (>", krt15_threshold, ")"))

print(p_krt15_filtered)
ggsave(
  paste0(output_dir, "KRT15_expression_filtered.png"), 
  p_krt15_filtered, 
  width = 10, 
  height = 8
)

# 3. SpatialDimPlot: Selected lineage cells
p_lineage <- SpatialDimPlot(
  vis_obj, 
  group.by = "KRT15_Lineage", 
  cols = c("Other" = "grey40", "KRT15_Lineage" = "red"),
  pt.size.factor = 3,
  image.alpha = 0
) + 
  ggtitle(paste0("KRT15 Lineage Selection (n=", sum(vis_obj$KRT15_Lineage == "KRT15_Lineage"), " cells)"))

print(p_lineage)
ggsave(
  paste0(output_dir, "KRT15_lineage_selection.png"), 
  p_lineage, 
  width = 10, 
  height = 8
)

# 4. Combined panel for comparison
p_combined_selection <- patchwork::wrap_plots(
  p_krt15_all + theme(legend.position = "bottom"),
  p_krt15_filtered + theme(legend.position = "bottom"),
  p_lineage + theme(legend.position = "bottom"),
  ncol = 3,
  image.alpha = 0
)

ggsave(
  paste0(output_dir, "KRT15_selection_overview.png"),
  p_combined_selection,
  width = 18,
  height = 6
)

message(paste("Total cells:", ncol(vis_obj)))
message(paste("KRT15+ cells (>", krt15_threshold, "):", sum(FetchData(vis_obj, vars = "KRT15")$KRT15 > krt15_threshold)))
message(paste("Selected lineage cells:", sum(vis_obj$KRT15_Lineage == "KRT15_Lineage")))

lineage_cells <- colnames(vis_obj)[vis_obj$KRT15_Lineage == "KRT15_Lineage"]
subset_seurat <- subset(vis_obj, cells = lineage_cells)

message(paste("Lineage cells:", length(lineage_cells)))

#################################################################
# CONVERT SEURAT TO MONOCLE3 CDS
#################################################################

expression_matrix <- GetAssayData(subset_seurat, layer = "counts")
cell_metadata <- subset_seurat@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

#################################################################
# PREPROCESS AND USE SPATIAL COORDINATES
#################################################################

cds <- preprocess_cds(cds, num_dim = 30)

# Use spatial coordinates as embedding
spatial_coords <- GetTissueCoordinates(subset_seurat)
spatial_matrix <- as.matrix(spatial_coords[, c("x", "y")])
rownames(spatial_matrix) <- rownames(spatial_coords)
spatial_matrix <- spatial_matrix[colnames(cds), ]

reducedDims(cds)[["UMAP"]] <- spatial_matrix

#################################################################
# CLUSTER AND LEARN GRAPH
#################################################################

cds <- cluster_cells(cds, cluster_method = "louvain")
cds <- learn_graph(cds, use_partition = FALSE)

p_graph <- plot_cells(
  cds,
  color_cells_by = "cluster",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)
print(p_graph)
ggsave(paste0(output_dir, "monocle3_trajectory_graph.png"), p_graph, width = 8, height = 6)

#################################################################
# ORDER CELLS - RIGHTMOST AS ROOT
#################################################################

coords <- reducedDims(cds)[["UMAP"]]
rightmost_cell <- rownames(coords)[which.max(coords[, 1])]

cell_to_node <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
root_node_idx <- cell_to_node[rightmost_cell, 1]
graph_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name
root_node_name <- graph_nodes[root_node_idx]

message(paste("Root node:", root_node_name))

cds <- order_cells(cds, root_pr_nodes = root_node_name)

#################################################################
# EXTRACT PSEUDOTIME
#################################################################

pseudotime_values <- pseudotime(cds)
subset_seurat$pseudotime <- pseudotime_values[colnames(subset_seurat)]

# Add pseudotime bins for visualization
n_pt_bins <- 20
colData(cds)$pseudotime_bin <- cut(
  pseudotime_values,
  breaks = n_pt_bins,
  labels = paste0("PT_", seq_len(n_pt_bins))
)

p_pt <- SpatialFeaturePlot(
  subset_seurat,
  features = "pseudotime",
  pt.size.factor = 3
) + scale_fill_viridis(option = "B")

print(p_pt)
ggsave(paste0(output_dir, "spatial_pseudotime.png"), p_pt, width = 8, height = 6)

#################################################################
# GRAPH TEST - IDENTIFY TRAJECTORY GENES
#################################################################

message("Running graph_test on principal graph...")

trajectory_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

# Filter significant genes and rank by effect size
sig_genes_df <- trajectory_genes %>%
  filter(q_value < 0.05) %>%
  arrange(desc(morans_I))

sig_gene_ids <- sig_genes_df$gene_short_name

# FORCE INCLUDE PRIORITY GENES regardless of significance
priority_in_data <- priority_genes[priority_genes %in% rowData(cds)$gene_short_name]
sig_gene_ids <- unique(c(sig_gene_ids, priority_in_data))

message(paste("Significant trajectory genes:", length(sig_gene_ids)))
message(paste("Priority genes included:", sum(priority_in_data %in% sig_gene_ids), "/", length(priority_in_data)))
#################################################################
# FIND GENE MODULES
#################################################################

message("Finding gene modules...")

# Subset to significant genes
cds_sig <- cds[rowData(cds)$gene_short_name %in% sig_gene_ids, ]

gene_module_df <- find_gene_modules(cds_sig, resolution = module_resolution)

message(paste("Number of modules:", length(unique(gene_module_df$module))))

# Save module assignments
write.csv(gene_module_df, paste0(output_dir, "gene_modules.csv"), row.names = FALSE)

#################################################################
# AGGREGATE MODULE EXPRESSION BY PSEUDOTIME BIN
#################################################################

cell_group_df <- tibble::tibble(
  cell = colnames(cds),
  cell_group = as.character(colData(cds)$pseudotime_bin)
)

agg_mat <- aggregate_gene_expression(cds_sig, gene_module_df, cell_group_df)
row.names(agg_mat) <- paste0("Module_", row.names(agg_mat))

# Order columns by pseudotime
pt_order <- paste0("PT_", seq_len(n_pt_bins))
agg_mat <- agg_mat[, pt_order[pt_order %in% colnames(agg_mat)]]

# Order modules by peak pseudotime position
module_peak <- apply(agg_mat, 1, which.max)
module_order <- names(sort(module_peak))
agg_mat <- agg_mat[module_order, ]

#################################################################
# GENE-LEVEL HEATMAP ALONG PSEUDOTIME
#################################################################

message("Building gene-level pseudotime heatmap...")
library(zoo)

# Order cells by pseudotime
pt_order <- order(pseudotime(cds_sig))
expr_mat <- as.matrix(exprs(cds_sig)[, pt_order])

# Smooth expression (rolling mean across cells)
smooth_window <- min(50, floor(ncol(expr_mat) / 10))

expr_smoothed <- t(apply(expr_mat, 1, function(x) {
  zoo::rollmean(x, k = smooth_window, fill = NA, align = "center")
}))

# Remove edge NAs
valid_cols <- complete.cases(t(expr_smoothed))
expr_smoothed <- expr_smoothed[, valid_cols]

# Filter lowly expressed genes BUT KEEP PRIORITY GENES
gene_means <- rowMeans(expr_smoothed, na.rm = TRUE)
genes_to_keep <- gene_means > min_expr | rownames(expr_smoothed) %in% priority_genes
expr_smoothed <- expr_smoothed[genes_to_keep, ]

message(paste("Genes after expression filter:", nrow(expr_smoothed)))

# Z-scale per gene
expr_scaled <- t(scale(t(expr_smoothed)))
expr_scaled[is.na(expr_scaled)] <- 0
expr_scaled[expr_scaled > zscore_cap] <- zscore_cap
expr_scaled[expr_scaled < -zscore_cap] <- -zscore_cap

# Order genes by peak pseudotime position
peak_position <- apply(expr_scaled, 1, which.max)
gene_order <- names(sort(peak_position))
expr_final <- expr_scaled[gene_order, ]

message(paste("Genes in heatmap:", nrow(expr_final)))

#################################################################
# LABEL KEY GENES ON HEATMAP
#################################################################

gene_peak_bins <- apply(expr_final, 1, which.max)
n_cols <- ncol(expr_final)

# Priority genes - check which are actually in the heatmap
priority_in_heatmap <- priority_genes[priority_genes %in% rownames(expr_final)]
priority_missing <- priority_genes[!priority_genes %in% rownames(expr_final)]

if (length(priority_missing) > 0) {
  message(paste("Priority genes NOT in heatmap:", paste(priority_missing, collapse = ", ")))
}
message(paste("Priority genes in heatmap:", paste(priority_in_heatmap, collapse = ", ")))

priority_positions <- match(priority_in_heatmap, rownames(expr_final))

# Remove any NA positions (shouldn't happen now but safety check)
valid_priority <- !is.na(priority_positions)
priority_in_heatmap <- priority_in_heatmap[valid_priority]
priority_positions <- priority_positions[valid_priority]

# Top gene per ~10% of pseudotime (excluding priority genes)
n_label_bins <- 20
col_bins <- cut(seq_len(n_cols), breaks = n_label_bins, labels = FALSE)

bin_top_genes <- character(0)
bin_top_positions <- integer(0)

for (b in seq_len(n_label_bins)) {
  bin_cols <- which(col_bins == b)
  genes_peaking_here <- names(gene_peak_bins)[gene_peak_bins %in% bin_cols]
  genes_peaking_here <- setdiff(genes_peaking_here, c(priority_in_heatmap, bin_top_genes))
  
  if (length(genes_peaking_here) > 0) {
    peak_vals <- sapply(genes_peaking_here, function(g) {
      max(expr_final[g, bin_cols])
    })
    top_gene <- genes_peaking_here[which.max(peak_vals)]
    bin_top_genes <- c(bin_top_genes, top_gene)
    bin_top_positions <- c(bin_top_positions, match(top_gene, rownames(expr_final)))
  }
}

# Combine labels - priority genes first, then bin top genes
all_label_genes <- c(priority_in_heatmap, bin_top_genes)
all_label_positions <- c(priority_positions, bin_top_positions)

# Final safety: remove any NAs
valid_labels <- !is.na(all_label_positions)
all_label_genes <- all_label_genes[valid_labels]
all_label_positions <- all_label_positions[valid_labels]

message(paste("Total genes labeled:", length(all_label_genes)))
message(paste("Labels:", paste(all_label_genes, collapse = ", ")))

# Row annotation
row_anno <- rowAnnotation(
  link = anno_mark(
    at = all_label_positions,
    labels = all_label_genes,
    labels_gp = gpar(fontsize = 7),
    link_width = unit(5, "mm")
  )
)

#################################################################
# DRAW GENE HEATMAP
#################################################################

ht_genes <- Heatmap(
  expr_final,
  name = "Z-score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = inferno(100),
  right_annotation = row_anno,
  column_title = "KRT15 Differentiation Trajectory DE Genes ->",
  row_title = paste0("Genes (n = ", nrow(expr_final), ")"),
  use_raster = TRUE,
  raster_quality = 2
)

pdf(paste0(output_dir, "gene_pseudotime_heatmap.pdf"), width = 8, height = 10)
draw(ht_genes)
dev.off()

png(paste0(output_dir, "gene_pseudotime_heatmap.png"), width = 8, height = 10, units = "in", res = 300)
draw(ht_genes)
dev.off()

message("Gene-level heatmap saved.")


#################################################################
# MODULE HEATMAP
#################################################################

ht_modules <- Heatmap(
  as.matrix(agg_mat),
  name = "Expression",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = inferno(100),
  column_title = "KRT15 Differentiation Trajectory",
  row_title = "Gene Modules",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)

pdf(paste0(output_dir, "module_pseudotime_heatmap.pdf"), width = 8, height = 6)
draw(ht_modules)
dev.off()

png(paste0(output_dir, "module_pseudotime_heatmap.png"), width = 8, height = 6, units = "in", res = 300)
draw(ht_modules)
dev.off()

#################################################################
# VISUALIZE MODULES ON SPATIAL EMBEDDING
#################################################################

# Select top modules by variance
module_var <- apply(agg_mat, 1, var)
top_modules <- names(sort(module_var, decreasing = TRUE))[1:min(6, nrow(agg_mat))]
top_module_nums <- gsub("Module_", "", top_modules)

p_modules <- plot_cells(
  cds_sig,
  genes = gene_module_df %>% filter(module %in% top_module_nums),
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE
)

ggsave(paste0(output_dir, "top_modules_spatial.png"), p_modules, width = 12, height = 8)

#################################################################
# PRIORITY GENES: PSEUDOTIME PLOTS
#################################################################

message("Generating priority gene plots...")

# Filter to priority genes present in data
priority_present <- priority_genes[priority_genes %in% rowData(cds)$gene_short_name]
message(paste("Priority genes found:", paste(priority_present, collapse = ", ")))

if (length(priority_present) > 0) {
  
  cds_priority <- cds[rowData(cds)$gene_short_name %in% priority_present, ]
  
  # Plot genes in pseudotime
  p_pseudotime <- plot_genes_in_pseudotime(
    cds_priority,
    color_cells_by = "pseudotime_bin",
    min_expr = 1,
    ncol = 3
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    paste0(marker_dir, "priority_genes_pseudotime.png"),
    p_pseudotime,
    width = 12,
    height = 8
  )
  
  ggsave(
    paste0(marker_dir, "priority_genes_pseudotime.pdf"),
    p_pseudotime,
    width = 12,
    height = 8
  )
  
  # Plot genes violin/hybrid by pseudotime bin
  # Note: plot_genes_hybrid may only be in monocle3 dev branch
  # Fall back to violin if not available
  tryCatch({
    p_hybrid <- plot_genes_hybrid(
      cds_priority,
      group_cells_by = "pseudotime_bin",
      ncol = 3
    ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      paste0(marker_dir, "priority_genes_hybrid.png"),
      p_hybrid,
      width = 12,
      height = 8
    )
  }, error = function(e) {
    message("plot_genes_hybrid not available, using violin plot instead")
    
    p_violin <- plot_genes_violin(
      cds_priority,
      group_cells_by = "pseudotime_bin",
      ncol = 3
    ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      paste0(marker_dir, "priority_genes_violin.png"),
      p_violin,
      width = 12,
      height = 8
    )
  })
  
  # Individual gene plots
  for (gene in priority_present) {
    cds_gene <- cds[rowData(cds)$gene_short_name == gene, ]
    
    p_gene <- plot_genes_in_pseudotime(
      cds_gene,
      color_cells_by = "pseudotime_bin",
      min_expr = 0.5
    ) + 
      ggtitle(gene) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      paste0(marker_dir, gene, "_pseudotime.png"),
      p_gene,
      width = 6,
      height = 4
    )
  }
}

#################################################################
# GENES PER MODULE TABLE
#################################################################

# Create summary of top genes per module
module_summary <- gene_module_df %>%
  left_join(sig_genes_df, by = c("id" = "id")) %>%
  group_by(module) %>%
  arrange(desc(morans_I)) %>%
  slice_head(n = 10) %>%
  summarise(
    n_genes = n(),
    top_genes = paste(gene_short_name, collapse = ", ")
  )

write.csv(module_summary, paste0(output_dir, "module_summary.csv"), row.names = FALSE)

#################################################################
# SAVE OBJECTS
#################################################################

saveRDS(cds, paste0(output_dir, "monocle3_cds.rds"))
saveRDS(subset_seurat, paste0(output_dir, "subset_seurat_with_pseudotime.rds"))

message("Analysis complete.")
message(paste("Main outputs:", output_dir))
message(paste("Marker gene plots:", marker_dir))

