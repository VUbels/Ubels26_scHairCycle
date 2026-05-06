#!/usr/bin/env Rscript

#################################################################
# Pre-Integrated Follicle Trajectory Analysis
# 
# CORRECT WORKFLOW ORDER:
#   1. Integrate Visium + scATAC (FULL datasets)
#   2. Preprocess integrated object
#   3. THEN select cells with spatial lasso (on integrated data)
#   4. Run Monocle3 on selected subset
#   5. Visualize RNA/Chromatin/TF along unified pseudotime
#################################################################

library(Seurat)
library(Signac)
library(ArchR)
library(SummarizedExperiment)  # For rowData() on ArchR matrices
library(monocle3)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
library(Matrix)
library(zoo)
library(circlize)
library(shiny)
library(plotly)

set.seed(42)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

#################################################################
# HELPER FUNCTIONS
#################################################################

spatial_lasso_selector <- function(obj, feature, image_name = NULL) {
  
  if (is.null(image_name)) {
    image_name <- Images(obj)[1]
  }
  
  coords <- GetTissueCoordinates(obj, image = image_name)
  expr <- FetchData(obj, vars = feature)[, 1]
  
  df <- data.frame(
    x = coords[[1]],
    y = coords[[2]],
    expr = expr,
    barcode = rownames(coords),
    stringsAsFactors = FALSE
  )
  
  selected_barcodes <- NULL
  
  ui <- fluidPage(
    titlePanel(paste("Select ROI -", feature)),
    fluidRow(
      column(9, plotlyOutput("spatial_plot", height = "600px")),
      column(3,
             h4("Selected spots:"),
             verbatimTextOutput("n_selected"),
             actionButton("done", "Confirm Selection", class = "btn-success"),
             hr(),
             helpText("Use lasso or box select tool in plot")
      )
    )
  )
  
  server <- function(input, output, session) {
    output$spatial_plot <- renderPlotly({
      plot_ly(
        df,
        x = ~x,
        y = ~y,
        color = ~expr,
        text = ~barcode,
        customdata = ~barcode,
        type = "scatter",
        mode = "markers",
        marker = list(size = 5),
        colors = viridis::viridis(100),
        source = "spatial"
      ) %>%
        layout(
          dragmode = "lasso",
          yaxis = list(scaleanchor = "x", scaleratio = 1, autorange = "reversed")
        )
    })
    
    output$n_selected <- renderPrint({
      d <- event_data("plotly_selected", source = "spatial")
      if (is.null(d)) {
        cat("0 spots selected")
      } else {
        cat(nrow(d), "spots selected")
      }
    })
    
    observeEvent(input$done, {
      d <- event_data("plotly_selected", source = "spatial")
      if (!is.null(d)) {
        selected_barcodes <<- d$customdata
      }
      stopApp(selected_barcodes)
    })
  }
  
  result <- runApp(shinyApp(ui, server))
  return(result)
}

add_lineage_label <- function(obj, 
                              selected_cells, 
                              feature, 
                              threshold, 
                              label_name = NULL,
                              assay = NULL) {
  
  if (is.null(label_name)) {
    label_name <- paste0(feature, "_Lineage")
  }
  
  # Direct data access is much faster than FetchData
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }
  
  # Check if assay exists
  if (!assay %in% names(obj@assays)) {
    stop(paste("Assay", assay, "not found. Available:", 
               paste(names(obj@assays), collapse = ", ")))
  }
  
  # Get expression directly from the data layer
  data_mat <- GetAssayData(obj, assay = assay, layer = "data")
  
  if (!feature %in% rownames(data_mat)) {
    stop(paste("Feature", feature, "not found in assay", assay))
  }
  
  expr <- data_mat[feature, ]
  
  in_region <- colnames(obj) %in% selected_cells
  above_threshold <- expr >= threshold
  
  # Handle NA values (e.g., scATAC cells won't have Spatial assay data)
  above_threshold[is.na(above_threshold)] <- FALSE
  
  obj@meta.data[[label_name]] <- in_region & above_threshold
  
  message(paste0(
    "Added '", label_name, "': ",
    sum(obj@meta.data[[label_name]]), " TRUE / ",
    sum(!obj@meta.data[[label_name]]), " FALSE"
  ))
  
  return(obj)
}

#################################################################
# SETUP
#################################################################

main_folder <- "./"
output_dir <- paste0(main_folder, "integrated_follicle_trajectory/")
marker_dir <- paste0(output_dir, "marker_genes/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)

#################################################################
# PARAMETERS
#################################################################

# Data paths
visium_path <- paste0(main_folder, "Spatial_scalp_S2_final.rds")
archr_path <- "./scatac_files/archr_proj"

# Integration
n_integration_dims <- 30
n_anchor_features <- 3000

# Lineage selection (applied AFTER integration)
krt15_threshold <- 0.5
lineage_marker <- "KRT15"

# Heatmap
smooth_window <- 50
min_expr <- 0.1
zscore_cap <- 3
n_pt_bins <- 100

# Priority genes
priority_genes <- c("CD34", "LGR5", "CD200", "COMP", "KRT19", "DKK3", 
                    "TGFB2", "DIO2", "ANGPTL2", "WIF1", "LGR6", 
                    "SOX9", "LHX2", "KRT15", "NFATC1", "KRT14", "KRT5")

#################################################################
# STEP 1: LOAD ALL DATA
#################################################################

message("================================================================")
message("   STEP 1: LOADING DATA                                        ")
message("================================================================")

# Load Visium
vis_obj <- readRDS(visium_path)
vis_obj <- UpdateSeuratObject(vis_obj)
vis_obj$modality <- "Visium"
message(paste("Visium spots:", ncol(vis_obj)))

# Verify Visium has gene names
vis_test_genes <- tryCatch(rownames(vis_obj), error = function(e) NULL)
if (is.null(vis_test_genes) || length(vis_test_genes) == 0) {
  # Try to get from assay
  default_assay <- DefaultAssay(vis_obj)
  vis_test_genes <- tryCatch(rownames(vis_obj[[default_assay]]), error = function(e) NULL)
}
message(paste("Visium genes:", length(vis_test_genes)))

# Load ArchR project (no copying)
skin_archr <- loadArchRProject(path = archr_path)
message(paste("scATAC cells:", nCells(skin_archr)))

# Identify healthy control cells (avoid AA chromatin bias)
cell_meta <- getCellColData(skin_archr)
healthy_samples <- grep("^C_", unique(cell_meta$Sample), value = TRUE)
healthy_cells <- rownames(cell_meta)[cell_meta$Sample %in% healthy_samples]
message(paste("Healthy scATAC cells:", length(healthy_cells)))

rm(cell_meta)
gc()

#################################################################
# STEP 2: PREPARE scATAC AS SEURAT OBJECT (LogNormalize)
# NOTE: SCT is NOT recommended for scATAC gene activity matrices
#################################################################

message("================================================================")
message("   STEP 2: PREPARING scATAC FOR INTEGRATION                    ")
message("================================================================")

# Extract gene scores for ALL cells, then subset to healthy
gsm <- getMatrixFromProject(skin_archr, useMatrix = "GeneScoreMatrix")
gs_counts <- assay(gsm)

# Get gene names from rowData (ArchR stores them there)
gene_names <- rowData(gsm)$name
if (is.null(gene_names)) {
  # Fallback: try rownames directly
  gene_names <- rownames(gsm)
}
if (is.null(gene_names) || length(gene_names) != nrow(gs_counts)) {
  stop("Could not extract gene names from ArchR GeneScoreMatrix")
}

# Set rownames on the matrix
rownames(gs_counts) <- gene_names
message(paste("Gene names extracted:", length(gene_names)))

# Subset to healthy cells only
healthy_in_matrix <- colnames(gs_counts)[colnames(gs_counts) %in% healthy_cells]
gs_counts <- gs_counts[, healthy_in_matrix, drop = FALSE]
message(paste("Gene score matrix:", nrow(gs_counts), "genes x", ncol(gs_counts), "healthy cells"))

# Create Seurat object with ACTIVITY assay
atac_seurat <- CreateSeuratObject(
  counts = gs_counts,
  assay = "ACTIVITY",
  meta.data = as.data.frame(getCellColData(skin_archr))[healthy_in_matrix, , drop = FALSE]
)

# Clean up large matrices no longer needed
rm(gs_counts, gsm, gene_names, healthy_in_matrix)
gc()

# Verify gene names transferred
message(paste("scATAC genes in Seurat:", nrow(atac_seurat)))

# LogNormalize (NOT SCT - SCT is inappropriate for gene activity scores)
atac_seurat <- NormalizeData(
  atac_seurat, 
  assay = "ACTIVITY",
  normalization.method = "LogNormalize",
  scale.factor = median(atac_seurat$nCount_ACTIVITY)
)
atac_seurat <- FindVariableFeatures(atac_seurat, selection.method = "vst", nfeatures = n_anchor_features)
atac_seurat <- ScaleData(atac_seurat)

# Add LSI embedding from ArchR
lsi_embed <- getReducedDims(skin_archr, reducedDims = "IterativeLSI")
if ("Harmony" %in% names(skin_archr@reducedDims)) {
  lsi_embed <- getReducedDims(skin_archr, reducedDims = "Harmony")
}

# Subset LSI to healthy cells
lsi_embed <- lsi_embed[colnames(atac_seurat), , drop = FALSE]

# Detect available LSI dimensions
n_lsi_dims <- ncol(lsi_embed)
message(paste("LSI dimensions available:", n_lsi_dims))

atac_seurat[["lsi"]] <- CreateDimReducObject(
  embeddings = lsi_embed,
  key = "LSI_",
  assay = "ACTIVITY"
)

atac_seurat$modality <- "scATAC"
message(paste("scATAC Seurat object:", ncol(atac_seurat), "cells"))

#################################################################
# STEP 3: PREPARE VISIUM (LogNormalize to match scATAC)
#################################################################

message("================================================================")
message("   STEP 3: PREPARING VISIUM                                    ")
message("================================================================")

# Determine which assay to use (prefer Spatial, fall back to RNA)
available_assays <- names(vis_obj@assays)
vis_assay <- if ("Spatial" %in% available_assays) "Spatial" else "RNA"
DefaultAssay(vis_obj) <- vis_assay
message(paste("Using Visium assay:", vis_assay))

# LogNormalize for cross-modality integration (SCT cannot be used with scATAC ACTIVITY)
# Check if data layer exists; if so, check if normalization is needed
layers_available <- Layers(vis_obj[[vis_assay]])

if ("data" %in% layers_available) {
  # Check if already normalized by comparing counts vs data
  counts_sum <- sum(GetAssayData(vis_obj, assay = vis_assay, layer = "counts")[1:10, 1:min(10, ncol(vis_obj))])
  data_sum <- sum(GetAssayData(vis_obj, assay = vis_assay, layer = "data")[1:10, 1:min(10, ncol(vis_obj))])
  
  if (abs(counts_sum - data_sum) < 0.001) {
    message("Normalizing Visium with LogNormalize...")
    vis_obj <- NormalizeData(vis_obj, normalization.method = "LogNormalize")
  } else {
    message("Visium already has normalized data layer")
  }
} else {
  message("Normalizing Visium with LogNormalize...")
  vis_obj <- NormalizeData(vis_obj, normalization.method = "LogNormalize")
}

# Variable features and scaling
vis_obj <- FindVariableFeatures(vis_obj, selection.method = "vst", nfeatures = n_anchor_features)
vis_obj <- ScaleData(vis_obj, verbose = FALSE)

# PCA if not already computed
if (!"pca" %in% names(vis_obj@reductions)) {
  vis_obj <- RunPCA(vis_obj, verbose = FALSE)
}

message(paste("Visium ready:", ncol(vis_obj), "spots"))

#################################################################
# SPATIAL LASSO SELECTION (do early while vis_obj still in memory)
# This saves memory by allowing us to remove vis_obj before integration
#################################################################

message("Opening spatial lasso selector...")
message("Select the follicle region of interest, then click 'Confirm Selection'")

selected_visium_cells <- spatial_lasso_selector(vis_obj, lineage_marker)
message(paste("Selected Visium spots:", length(selected_visium_cells)))

gc()

#################################################################
# STEP 4: INTEGRATE VISIUM + scATAC (FULL DATASETS)
#################################################################

message("================================================================")
message("   STEP 4: INTEGRATING FULL VISIUM + scATAC                    ")
message("================================================================")

# Find shared variable features (with diagnostics)
message("--- Checking features ---")

# Get feature names robustly (Seurat v5 compatible) - multiple fallback methods
get_features <- function(obj, assay_name) {
  # Method 1: rownames of GetAssayData
  feats <- tryCatch({
    rownames(GetAssayData(obj, assay = assay_name, layer = "counts"))
  }, error = function(e) NULL)
  
  if (!is.null(feats) && length(feats) > 0) return(feats)
  
  # Method 2: rownames of the assay object
  feats <- tryCatch({
    rownames(obj[[assay_name]])
  }, error = function(e) NULL)
  
  if (!is.null(feats) && length(feats) > 0) return(feats)
  
  # Method 3: Features() function
  feats <- tryCatch({
    Features(obj[[assay_name]])
  }, error = function(e) NULL)
  
  if (!is.null(feats) && length(feats) > 0) return(feats)
  
  # Method 4: rownames of object
  feats <- tryCatch({
    rownames(obj)
  }, error = function(e) NULL)
  
  if (!is.null(feats) && length(feats) > 0) return(feats)
  
  # Method 5: names from assay@features slot (Seurat v5 specific)
  feats <- tryCatch({
    names(obj[[assay_name]]@features)
  }, error = function(e) NULL)
  
  return(feats)
}

vis_features <- get_features(vis_obj, vis_assay)
atac_features <- get_features(atac_seurat, "ACTIVITY")

message(paste("Visium features:", length(vis_features)))
message(paste("scATAC features:", length(atac_features)))

if (is.null(vis_features) || length(vis_features) == 0) {
  stop("Could not extract feature names from Visium object")
}
if (is.null(atac_features) || length(atac_features) == 0) {
  stop("Could not extract feature names from scATAC object")
}

# Get variable features (may be empty if not computed)
vis_var <- tryCatch(VariableFeatures(vis_obj), error = function(e) character(0))
atac_var <- tryCatch(VariableFeatures(atac_seurat), error = function(e) character(0))

message(paste("Visium variable features:", length(vis_var)))
message(paste("scATAC variable features:", length(atac_var)))

# Find shared genes
all_shared <- intersect(vis_features, atac_features)
message(paste("Total shared features:", length(all_shared)))

# Determine which features to use for integration
if (length(vis_var) > 0 && length(atac_var) > 0) {
  shared_genes <- intersect(vis_var, atac_var)
  message(paste("Shared variable features:", length(shared_genes)))
  
  if (length(shared_genes) < 500) {
    # Use union of variable features that exist in both
    union_var <- union(vis_var, atac_var)
    shared_genes <- intersect(union_var, all_shared)
    message(paste("Using union of variable features:", length(shared_genes)))
  }
} else {
  shared_genes <- character(0)
}

# Final fallback: use all shared genes ranked by variance
if (length(shared_genes) < 500) {
  message("Insufficient variable features. Computing variance-based ranking...")
  
  # Get expression data
  vis_data <- GetAssayData(vis_obj, assay = vis_assay, layer = "data")
  vis_data <- vis_data[all_shared, , drop = FALSE]
  
  # Calculate variance for each gene (sparse matrix compatible)
  gene_means <- Matrix::rowMeans(vis_data)
  gene_vars <- Matrix::rowMeans(vis_data^2) - gene_means^2
  names(gene_vars) <- all_shared
  
  # Take top genes by variance
  shared_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(n_anchor_features, length(all_shared))]
  message(paste("Using top variance genes:", length(shared_genes)))
  
  rm(vis_data, gene_means, gene_vars)
}

# Clean up feature vectors
rm(vis_features, atac_features, vis_var, atac_var, all_shared)
gc()

stopifnot("No shared genes found between Visium and scATAC" = length(shared_genes) >= 100)

# Determine dimensions to use based on available LSI dimensions
# LSI dim 1 correlates with sequencing depth (technical), so we skip it for weight.reduction
# See: https://stuartlab.org/signac/articles/pbmc_vignette.html
n_dims_for_cca <- min(n_integration_dims, n_lsi_dims)
n_dims_for_transfer <- min(n_integration_dims, n_lsi_dims)  # Will use 2:n_dims_for_transfer

message(paste("Using CCA dims: 1:", n_dims_for_cca, sep = ""))
message(paste("Using LSI dims for transfer: 2:", n_dims_for_transfer, " (skipping dim 1 = sequencing depth)", sep = ""))

# Find transfer anchors (LogNormalize for both)
transfer_anchors <- FindTransferAnchors(
  reference = vis_obj,
  query = atac_seurat,
  features = shared_genes,
  reference.assay = vis_assay,
  query.assay = "ACTIVITY",
  normalization.method = "LogNormalize",
  reduction = "cca",
  dims = 1:n_dims_for_cca,
  verbose = TRUE
)

# Impute RNA expression into scATAC cells
# dims = 2:n because LSI dim 1 captures sequencing depth, not biological signal
imputed_data <- TransferData(
  anchorset = transfer_anchors,
  refdata = GetAssayData(vis_obj, assay = vis_assay, layer = "data"),
  weight.reduction = atac_seurat[["lsi"]],
  dims = 2:n_dims_for_transfer
)

# TransferData returns predicted values (data layer only, no counts)
# For Seurat v5, we need to properly add this assay
if (inherits(imputed_data, "Assay5") || inherits(imputed_data, "Assay")) {
  atac_seurat[["RNA_imputed"]] <- imputed_data
} else {
  # If it returns a matrix, create assay manually
  atac_seurat[["RNA_imputed"]] <- CreateAssayObject(data = imputed_data)
}
message("RNA imputed into scATAC")

# Clean up transfer objects
rm(transfer_anchors, imputed_data, shared_genes)
gc()

# Prepare for integration
DefaultAssay(atac_seurat) <- "RNA_imputed"

# Create merged object for integration
vis_obj$orig_barcode <- colnames(vis_obj)
atac_seurat$orig_barcode <- colnames(atac_seurat)

# Select integration features from both
integration_features <- SelectIntegrationFeatures(
  object.list = list(vis_obj, atac_seurat),
  nfeatures = n_anchor_features
)

# Prepare objects for RPCA integration
object_list <- list(Visium = vis_obj, scATAC = atac_seurat)

object_list <- lapply(object_list, function(x) {
  # Use appropriate assay based on modality
  if (x$modality[1] == "Visium") {
    current_assay <- if ("Spatial" %in% names(x@assays)) "Spatial" else "RNA"
  } else {
    current_assay <- "RNA_imputed"
  }
  DefaultAssay(x) <- current_assay
  x <- ScaleData(x, features = integration_features, verbose = FALSE)
  # Use standard "pca" name so FindIntegrationAnchors can find it
  x <- RunPCA(x, features = integration_features, verbose = FALSE, reduction.name = "pca")
  return(x)
})

# Determine integration dimensions based on available PCA dimensions
n_int_pca <- min(sapply(object_list, function(x) ncol(Embeddings(x, "pca"))))
n_int_dims <- min(n_integration_dims, n_int_pca)
message(paste("Using", n_int_dims, "dimensions for RPCA integration"))

# Find integration anchors
int_anchors <- FindIntegrationAnchors(
  object.list = object_list,
  anchor.features = integration_features,
  normalization.method = "LogNormalize",
  reduction = "rpca",
  dims = 1:n_int_dims
)

# Integrate data
integrated_obj <- IntegrateData(
  anchorset = int_anchors,
  normalization.method = "LogNormalize",
  dims = 1:n_int_dims
)

message(paste("Integrated object:", ncol(integrated_obj), "cells/spots"))

# Clean up integration intermediates - these are large
rm(object_list, int_anchors, integration_features)
# atac_seurat and vis_obj no longer needed - we have integrated_obj and skin_archr for ATAC data
rm(atac_seurat, vis_obj)
gc()

#################################################################
# STEP 5: PREPROCESS INTEGRATED OBJECT
#################################################################

message("================================================================")
message("   STEP 5: PREPROCESSING INTEGRATED OBJECT                     ")
message("================================================================")

DefaultAssay(integrated_obj) <- "integrated"
integrated_obj <- ScaleData(integrated_obj, verbose = FALSE)
integrated_obj <- RunPCA(integrated_obj, verbose = FALSE)

# Use available PCA dimensions
n_pca_avail <- ncol(Embeddings(integrated_obj, "pca"))
n_dims_use <- min(n_integration_dims, n_pca_avail)
message(paste("Using", n_dims_use, "PCA dimensions for UMAP/clustering"))

integrated_obj <- RunUMAP(integrated_obj, dims = 1:n_dims_use, verbose = FALSE)
integrated_obj <- FindNeighbors(integrated_obj, dims = 1:n_dims_use, verbose = FALSE)
integrated_obj <- FindClusters(integrated_obj, resolution = 0.5, verbose = FALSE)

# Visualize full integration
p_int <- DimPlot(integrated_obj, group.by = "modality",
                 cols = c("Visium" = "steelblue", "scATAC" = "coral")) +
  ggtitle("Integrated Visium + scATAC (All Cells)")
print(p_int)
ggsave(paste0(output_dir, "integration_full_umap.png"), p_int, width = 8, height = 6)

# Save integrated object before selection
saveRDS(integrated_obj, paste0(output_dir, "integrated_full_object.rds"))
message("Full integrated object saved")

# Clean up - we still need vis_obj for spatial selection but can slim down atac_seurat
# Keep atac_seurat minimal (will need for LSI later but remove unneeded assays)
gc()

#################################################################
# STEP 6: SELECT LINEAGE CELLS (AFTER INTEGRATION)
# Spatial selection was done in STEP 1 - now apply to integrated object
#################################################################

message("================================================================")
message("   STEP 6: SELECT LINEAGE CELLS (POST-INTEGRATION)             ")
message("================================================================")

# Apply lineage label to integrated object (Visium cells)
# selected_visium_cells was set in STEP 3 via spatial lasso
# Also apply threshold filter
# Use vis_assay (Spatial or RNA) since KRT15 is in original Visium assay, not "integrated"
integrated_obj <- add_lineage_label(
  obj = integrated_obj,
  selected_cells = selected_visium_cells,
  feature = lineage_marker,
  threshold = krt15_threshold,
  assay = vis_assay
)

# For scATAC cells: find those that cluster with selected Visium cells
# Use nearest neighbor in integrated space
visium_lineage_cells <- colnames(integrated_obj)[
  integrated_obj$modality == "Visium" & 
    integrated_obj[[paste0(lineage_marker, "_Lineage")]][,1] == TRUE
]

message(paste("Visium lineage cells:", length(visium_lineage_cells)))

# Find scATAC cells that are nearest neighbors to selected Visium cells
# in the integrated embedding
pca_embed <- Embeddings(integrated_obj, "pca")
n_pca_dims <- min(n_integration_dims, ncol(pca_embed))
pca_embed <- pca_embed[, 1:n_pca_dims, drop = FALSE]

visium_lineage_embed <- pca_embed[visium_lineage_cells, , drop = FALSE]
atac_cells <- colnames(integrated_obj)[integrated_obj$modality == "scATAC"]
atac_embed <- pca_embed[atac_cells, , drop = FALSE]

# For each scATAC cell, find distance to nearest Visium lineage cell
if (length(visium_lineage_cells) > 0 && length(atac_cells) > 0) {
  
  # Use FNN for fast nearest neighbor
  library(FNN)
  
  nn_result <- get.knnx(
    data = visium_lineage_embed,
    query = atac_embed,
    k = 1
  )
  
  # scATAC cells close to Visium lineage cells (distance threshold)
  dist_threshold <- quantile(nn_result$nn.dist, 0.5)  # Top 50% closest
  atac_lineage_cells <- atac_cells[nn_result$nn.dist[, 1] <= dist_threshold]
  
  message(paste("scATAC lineage cells (by proximity):", length(atac_lineage_cells)))
  
  # Mark in metadata
  integrated_obj$scATAC_Lineage <- colnames(integrated_obj) %in% atac_lineage_cells
  
  # Combined lineage
  all_lineage_cells <- c(visium_lineage_cells, atac_lineage_cells)
  integrated_obj$Combined_Lineage <- colnames(integrated_obj) %in% all_lineage_cells
  
} else {
  warning("No lineage cells found. Check selection.")
  all_lineage_cells <- visium_lineage_cells
  integrated_obj$Combined_Lineage <- colnames(integrated_obj) %in% all_lineage_cells
}

message(paste("Total lineage cells:", length(all_lineage_cells)))

# Clean up selection intermediates
rm(selected_visium_cells)
rm(pca_embed, visium_lineage_embed, atac_embed, atac_cells)
if (exists("nn_result")) rm(nn_result)
gc()

#################################################################
# STEP 7: SUBSET TO LINEAGE AND RUN MONOCLE3
#################################################################

message("================================================================")
message("   STEP 7: TRAJECTORY ON SELECTED LINEAGE                      ")
message("================================================================")

# Subset to lineage cells
# First slim down the object by dropping non-essential assays to speed up subset
message("Slimming integrated object before subset...")
assays_to_keep <- c("integrated", vis_assay)  # Keep integrated + Spatial/RNA for images
assays_in_obj <- names(integrated_obj@assays)
for (a in assays_in_obj) {
  if (!a %in% assays_to_keep) {
    integrated_obj[[a]] <- NULL
  }
}
gc()

message("Subsetting to lineage cells...")
lineage_obj <- subset(integrated_obj, cells = all_lineage_cells)
message(paste("Lineage subset:", ncol(lineage_obj), "cells"))
message(paste("  - Visium:", sum(lineage_obj$modality == "Visium")))
message(paste("  - scATAC:", sum(lineage_obj$modality == "scATAC")))

# Store cell count for summary before removing
n_integrated_cells <- ncol(integrated_obj)

# integrated_obj has been saved to RDS - can remove to save memory
rm(integrated_obj, all_lineage_cells, visium_lineage_cells)
if (exists("atac_lineage_cells")) rm(atac_lineage_cells)
gc()

# Re-cluster lineage subset (use same dimension count as full integration)
lineage_obj <- FindNeighbors(lineage_obj, dims = 1:n_dims_use, verbose = FALSE)
lineage_obj <- FindClusters(lineage_obj, resolution = 0.5, verbose = FALSE)
lineage_obj <- RunUMAP(lineage_obj, dims = 1:n_dims_use, verbose = FALSE)

# Monocle3 trajectory
message("Converting to Monocle3 cell_data_set...")

# Create cell_data_set from Seurat
# Get expression matrix from integrated assay (ensure sparse format)
expr_mat_cds <- GetAssayData(lineage_obj, assay = "integrated", layer = "data")
if (!inherits(expr_mat_cds, "dgCMatrix")) {
  expr_mat_cds <- as(expr_mat_cds, "dgCMatrix")
}

# Create gene metadata (required by monocle3)
gene_metadata <- data.frame(
  gene_short_name = rownames(expr_mat_cds),
  row.names = rownames(expr_mat_cds)
)

# Create cell_data_set
cds <- new_cell_data_set(
  expression_data = expr_mat_cds,
  cell_metadata = lineage_obj@meta.data,
  gene_metadata = gene_metadata
)

# Set size factors (required for monocle3 since we skip preprocess_cds)
size_factors(cds) <- rep(1, ncol(cds))

# Transfer UMAP embedding from Seurat (skip monocle preprocessing)
umap_coords <- Embeddings(lineage_obj, "umap")
reducedDims(cds)[["UMAP"]] <- umap_coords

# Transfer PCA from Seurat
pca_coords <- Embeddings(lineage_obj, "pca")
reducedDims(cds)[["PCA"]] <- pca_coords

# Use Seurat clusters for Monocle partitions
cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3)

# Learn the trajectory graph
message("Learning trajectory graph...")
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

# Find root cell based on rightmost spatial coordinate (bulge region)
visium_cells_in_lineage <- colnames(lineage_obj)[lineage_obj$modality == "Visium"]
coords <- GetTissueCoordinates(lineage_obj, image = Images(lineage_obj)[1])
coords_lineage <- coords[rownames(coords) %in% visium_cells_in_lineage, ]

# Find cell with maximum x-coordinate as root
root_cell <- rownames(coords_lineage)[which.max(coords_lineage[[1]])]
message(paste("Root cell (rightmost x-coordinate):", root_cell))

# Get the principal graph node closest to root cell
cds <- order_cells(cds, root_cells = root_cell)

# Extract pseudotime
pt <- pseudotime(cds)
pt[is.infinite(pt)] <- max(pt[is.finite(pt)], na.rm = TRUE)

lineage_obj$pseudotime <- pt[colnames(lineage_obj)]

# Create bins with unique breaks
pt_breaks <- unique(quantile(lineage_obj$pseudotime, probs = seq(0, 1, length.out = n_pt_bins + 1), na.rm = TRUE))
n_actual_bins <- length(pt_breaks) - 1
message(paste("Pseudotime bins:", n_actual_bins))

lineage_obj$pseudotime_bin <- cut(
  lineage_obj$pseudotime,
  breaks = pt_breaks,
  labels = paste0("PT_", seq_len(n_actual_bins)),
  include.lowest = TRUE
)

# Visualize trajectory
p_pt <- FeaturePlot(lineage_obj, features = "pseudotime", pt.size = 0.5) +
  scale_color_viridis(option = "B") +
  ggtitle("Monocle3 Pseudotime (Integrated Lineage)")
ggsave(paste0(output_dir, "pseudotime_umap.png"), p_pt, width = 8, height = 6)

p_pt_split <- FeaturePlot(lineage_obj, features = "pseudotime",
                          split.by = "modality", pt.size = 0.3) +
  scale_color_viridis(option = "B")
ggsave(paste0(output_dir, "pseudotime_by_modality.png"), p_pt_split, width = 12, height = 5)

# Monocle trajectory plot
p_monocle <- plot_cells(cds, 
                        color_cells_by = "pseudotime",
                        label_cell_groups = FALSE,
                        label_leaves = FALSE,
                        label_branch_points = FALSE,
                        graph_label_size = 3) +
  scale_color_viridis(option = "B")
ggsave(paste0(output_dir, "monocle3_trajectory.png"), p_monocle, width = 8, height = 6)

# Clean up Monocle intermediates
rm(expr_mat_cds, gene_metadata, umap_coords, pca_coords, pt)
rm(visium_cells_in_lineage, coords, coords_lineage, root_cell)
rm(p_pt, p_pt_split, p_monocle)
gc()

#################################################################
# STEP 8: BUILD HEATMAPS ALONG UNIFIED PSEUDOTIME
#################################################################

message("================================================================")
message("   STEP 8: GENERATING TRAJECTORY HEATMAPS                      ")
message("================================================================")

# Order all cells by pseudotime
pt_order <- order(lineage_obj$pseudotime)
cell_order <- colnames(lineage_obj)[pt_order]

# Get expression from integrated assay
expr_mat <- as.matrix(GetAssayData(lineage_obj, assay = "integrated", layer = "data"))
expr_mat <- expr_mat[, cell_order]

#--- 8A: Gene Expression Heatmap ---
message("Building gene expression heatmap...")

expr_smoothed <- t(apply(expr_mat, 1, function(x) {
  zoo::rollmean(x, k = min(smooth_window, length(x) %/% 3), fill = NA, align = "center")
}))

valid_cols <- complete.cases(t(expr_smoothed))
expr_smoothed <- expr_smoothed[, valid_cols]

gene_means <- rowMeans(expr_smoothed, na.rm = TRUE)
genes_keep <- gene_means > min_expr | rownames(expr_smoothed) %in% priority_genes
expr_smoothed <- expr_smoothed[genes_keep, ]

expr_scaled <- t(scale(t(expr_smoothed)))
expr_scaled[is.na(expr_scaled)] <- 0
expr_scaled[expr_scaled > zscore_cap] <- zscore_cap
expr_scaled[expr_scaled < -zscore_cap] <- -zscore_cap

peak_pos <- apply(expr_scaled, 1, which.max)
gene_order <- names(sort(peak_pos))
expr_final <- expr_scaled[gene_order, ]

message(paste("Genes in heatmap:", nrow(expr_final)))

# Clean up expression intermediates
rm(expr_mat, expr_smoothed, expr_scaled, gene_means, genes_keep, peak_pos, gene_order, valid_cols)
gc()

#--- 8B: Chromatin Accessibility Heatmap ---
message("Building chromatin accessibility heatmap...")

atac_cells_lineage <- colnames(lineage_obj)[lineage_obj$modality == "scATAC"]
atac_cells_archr <- intersect(atac_cells_lineage, getCellNames(skin_archr))
# Also filter to healthy cells
atac_cells_archr <- intersect(atac_cells_archr, healthy_cells)

if (length(atac_cells_archr) > 0) {
  gsm_lineage <- getMatrixFromProject(skin_archr, useMatrix = "GeneScoreMatrix")
  atac_gs <- assay(gsm_lineage)[, atac_cells_archr]
  
  atac_pt <- lineage_obj$pseudotime[atac_cells_archr]
  atac_order <- order(atac_pt)
  atac_gs <- atac_gs[, atac_order]
  
  atac_smoothed <- t(apply(as.matrix(atac_gs), 1, function(x) {
    zoo::rollmean(x, k = min(smooth_window, length(x) %/% 3), fill = NA, align = "center")
  }))
  
  valid_atac <- complete.cases(t(atac_smoothed))
  atac_smoothed <- atac_smoothed[, valid_atac]
  
  atac_means <- rowMeans(atac_smoothed, na.rm = TRUE)
  atac_keep <- atac_means > min_expr | rownames(atac_smoothed) %in% priority_genes
  atac_smoothed <- atac_smoothed[atac_keep, ]
  
  atac_scaled <- t(scale(t(atac_smoothed)))
  atac_scaled[is.na(atac_scaled)] <- 0
  atac_scaled[atac_scaled > zscore_cap] <- zscore_cap
  atac_scaled[atac_scaled < -zscore_cap] <- -zscore_cap
  
  common_genes <- intersect(rownames(expr_final), rownames(atac_scaled))
  message(paste("Common genes:", length(common_genes)))
  
  rna_matched <- expr_final[common_genes, ]
  atac_matched <- atac_scaled[common_genes, ]
  
  # Clean up ATAC intermediates
  rm(gsm_lineage, atac_gs, atac_smoothed, atac_scaled, atac_means, atac_keep, valid_atac, atac_order, atac_pt)
  gc()
  
  has_atac <- TRUE
} else {
  message("No scATAC cells in lineage")
  has_atac <- FALSE
}

#--- 8C: TF Motif Activity Heatmap ---
message("Building TF motif heatmap...")

available_matrices <- getAvailableMatrices(skin_archr)

if ("MotifMatrix" %in% available_matrices && has_atac) {
  motif_mat <- getMatrixFromProject(skin_archr, useMatrix = "MotifMatrix")
  motif_dev <- assay(motif_mat, "deviations")
  
  motif_cells <- intersect(colnames(motif_dev), atac_cells_archr)
  motif_dev <- motif_dev[, motif_cells]
  
  motif_pt <- lineage_obj$pseudotime[motif_cells]
  motif_order <- order(motif_pt)
  motif_dev <- motif_dev[, motif_order]
  
  motif_smoothed <- t(apply(as.matrix(motif_dev), 1, function(x) {
    zoo::rollmean(x, k = min(smooth_window, length(x) %/% 3), fill = NA, align = "center")
  }))
  
  valid_motif <- complete.cases(t(motif_smoothed))
  motif_smoothed <- motif_smoothed[, valid_motif]
  
  motif_scaled <- t(scale(t(motif_smoothed)))
  motif_scaled[is.na(motif_scaled)] <- 0
  motif_scaled[motif_scaled > zscore_cap] <- zscore_cap
  motif_scaled[motif_scaled < -zscore_cap] <- -zscore_cap
  
  motif_peak <- apply(motif_scaled, 1, which.max)
  motif_order_final <- names(sort(motif_peak))
  motif_final <- motif_scaled[motif_order_final, ]
  
  # Clean up motif intermediates
  rm(motif_mat, motif_dev, motif_smoothed, motif_scaled, motif_peak, motif_order_final, valid_motif, motif_cells, motif_order, motif_pt)
  gc()
  
  has_motif <- TRUE
  message(paste("TF motifs:", nrow(motif_final)))
} else {
  has_motif <- FALSE
}

#################################################################
# STEP 9: DRAW HEATMAPS
#################################################################

message("================================================================")
message("   STEP 9: DRAWING HEATMAPS                                    ")
message("================================================================")

if (has_atac) {
  # Gene annotations
  priority_in_heatmap <- priority_genes[priority_genes %in% rownames(rna_matched)]
  priority_positions <- match(priority_in_heatmap, rownames(rna_matched))
  valid_idx <- !is.na(priority_positions)
  priority_positions <- priority_positions[valid_idx]
  priority_in_heatmap <- priority_in_heatmap[valid_idx]
  
  row_anno <- rowAnnotation(
    link = anno_mark(
      at = priority_positions,
      labels = priority_in_heatmap,
      labels_gp = gpar(fontsize = 7),
      link_width = unit(5, "mm")
    )
  )
  
  # RNA heatmap
  ht_rna <- Heatmap(
    rna_matched,
    name = "RNA\nZ-score",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    col = inferno(100),
    column_title = paste("Gene Expression\n(n=", ncol(rna_matched), ")"),
    row_title = paste0("Genes (n=", nrow(rna_matched), ")"),
    use_raster = TRUE,
    raster_quality = 2,
    width = unit(5, "cm")
  )
  
  # Chromatin heatmap
  ht_atac <- Heatmap(
    atac_matched,
    name = "Chromatin\nZ-score",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    col = viridis(100),
    column_title = paste("Chromatin Accessibility\n(n=", ncol(atac_matched), ")"),
    use_raster = TRUE,
    raster_quality = 2,
    right_annotation = row_anno,
    width = unit(5, "cm")
  )
  
  ht_combined <- ht_rna + ht_atac
  
  pdf(paste0(output_dir, "trajectory_rna_chromatin_heatmaps.pdf"), width = 12, height = 10)
  draw(ht_combined,
       column_title = paste(lineage_marker, "Differentiation: RNA | Chromatin"),
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()
  
  png(paste0(output_dir, "trajectory_rna_chromatin_heatmaps.png"),
      width = 12, height = 10, units = "in", res = 300)
  draw(ht_combined,
       column_title = paste(lineage_marker, "Differentiation: RNA | Chromatin"),
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()
  
  message("RNA + Chromatin heatmaps saved.")
  
  # Full heatmap with TF motifs
  if (has_motif) {
    ht_motif <- Heatmap(
      motif_final,
      name = "Motif\nDeviation",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
      column_title = paste("TF Activity\n(n=", ncol(motif_final), ")"),
      row_title = paste0("TFs (n=", nrow(motif_final), ")"),
      use_raster = TRUE,
      raster_quality = 2,
      width = unit(5, "cm")
    )
    
    ht_full <- ht_rna + ht_atac + ht_motif
    
    pdf(paste0(output_dir, "trajectory_full_heatmaps.pdf"), width = 16, height = 10)
    draw(ht_full,
         column_title = paste(lineage_marker, "Trajectory: RNA | Chromatin | TF Activity"),
         column_title_gp = gpar(fontsize = 14, fontface = "bold"))
    dev.off()
    
    png(paste0(output_dir, "trajectory_full_heatmaps.png"),
        width = 16, height = 10, units = "in", res = 300)
    draw(ht_full,
         column_title = paste(lineage_marker, "Trajectory: RNA | Chromatin | TF Activity"),
         column_title_gp = gpar(fontsize = 14, fontface = "bold"))
    dev.off()
    
    message("Full heatmaps saved.")
  }
  
  # Clean up heatmap objects
  rm(ht_rna, ht_atac, ht_combined, row_anno)
  if (has_motif) rm(ht_motif, ht_full)
  gc()
}

#################################################################
# STEP 10: CORRELATION ANALYSIS
#################################################################

if (has_atac) {
  message("\n=== RNA-CHROMATIN CORRELATION ===")
  
  n_bins <- min(ncol(rna_matched), ncol(atac_matched), 50)
  
  rna_binned <- t(apply(rna_matched, 1, function(x) {
    tapply(x, cut(seq_along(x), n_bins), mean, na.rm = TRUE)
  }))
  
  atac_binned <- t(apply(atac_matched, 1, function(x) {
    tapply(x, cut(seq_along(x), n_bins), mean, na.rm = TRUE)
  }))
  
  gene_cors <- sapply(common_genes, function(g) {
    cor(rna_binned[g, ], atac_binned[g, ], use = "pairwise.complete.obs")
  })
  
  cor_df <- data.frame(gene = names(gene_cors), correlation = gene_cors) %>%
    arrange(desc(correlation))
  
  write.csv(cor_df, paste0(output_dir, "rna_chromatin_correlations.csv"), row.names = FALSE)
  
  p_cor <- ggplot(cor_df, aes(x = correlation)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = mean(gene_cors, na.rm = TRUE), color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = "RNA-Chromatin Correlation", x = "Pearson Correlation", y = "Genes")
  ggsave(paste0(output_dir, "rna_chromatin_correlation.png"), p_cor, width = 8, height = 5)
  
  message(paste("Mean correlation:", round(mean(gene_cors, na.rm = TRUE), 3)))
  
  # Clean up correlation intermediates
  rm(rna_binned, atac_binned, p_cor)
  gc()
}

#################################################################
# SAVE OUTPUTS
#################################################################

message("\n=== SAVING OUTPUTS ===")

saveRDS(lineage_obj, paste0(output_dir, "lineage_trajectory_object.rds"))
saveRDS(cds, paste0(output_dir, "monocle3_cds.rds"))

traj_meta <- data.frame(
  cell = colnames(lineage_obj),
  pseudotime = lineage_obj$pseudotime,
  modality = lineage_obj$modality,
  cluster = lineage_obj$seurat_clusters
)
write.csv(traj_meta, paste0(output_dir, "trajectory_metadata.csv"), row.names = FALSE)

#################################################################
# SUMMARY
#################################################################

message("\n================================================================")
message("   COMPLETE                                                    ")
message("================================================================")
message(paste("Integrated cells:", n_integrated_cells))
message(paste("Lineage cells:", ncol(lineage_obj)))
message(paste("  - Visium:", sum(lineage_obj$modality == "Visium")))
message(paste("  - scATAC:", sum(lineage_obj$modality == "scATAC")))
if (has_atac) {
  message(paste("Genes in heatmap:", nrow(rna_matched)))
  message(paste("Mean RNA-Chromatin correlation:", round(mean(gene_cors, na.rm = TRUE), 3)))
}
if (has_motif) {
  message(paste("TF motifs:", nrow(motif_final)))
}
message("================================================================")