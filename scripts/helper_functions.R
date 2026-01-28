###################################################
# SIMPLE FUNCTION TO PLOT QC METRICS FOR ALL DATA
###################################################

show_qc_metrics <- function(input_folder) {
  
  post_aRNA_removal_dirs <- list.dirs(
    file.path(input_folder, "preprocessing"),
    recursive = FALSE
  )
  
  h5_files <- list.files(
                        post_aRNA_removal_dirs,
                        pattern = "_filtered_seurat.h5",
                        full.names = TRUE,
                        recursive = TRUE
                        )
  
  object.list <- list()
  
  for (i in seq_along(h5_files)) {
    
    object <- h5_files[[i]]  
    stage <- dataset_names[[i]]
    
    # Load CellBender-corrected data
    data.arna_corrected <- Read10X_h5(object, use.names = TRUE)  
    obj <- CreateSeuratObject(counts = data.arna_corrected, project = stage)
    obj$orig.ident <- stage
    
    # Calculate QC metrics
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    
    # QC visualization
    gg1 <- print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
    output_folder = paste0(input_folder, "/preprocessing/qc_metrics/")
    dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    
    file_name = paste0(output_folder, "qc_VlnPlot_", stage, "_.png")
    ggsave(file_name)
    
    plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    print(plot1 + plot2)
    
    file_name = paste0(output_folder, "qc_Scatterplot_", stage, "_.png")
    ggsave(file_name, width = 12)
    
    object.list[[i]] <- obj
    
    rm(data.arna_corrected)
    rm(obj)
  }
  
  return(object.list)
  
}

###################################################
# FILTER DATA AFTER QC METRIC VISUALIZATION
###################################################

filter_by_qc <- function(input_folder, project_names, min_feature = NULL, max_feature = NULL, min_count = NULL, max_count = NULL, max_percent_mt = NULL) {
  
  post_aRNA_removal_dirs <- list.dirs(
    file.path(input_folder, "preprocessing"),
    recursive = FALSE
  )
  
  h5_files <- list.files(
    post_aRNA_removal_dirs,
    pattern = "_filtered_seurat\\.h5$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  # DEBUG: Show file order
  cat("Files found (in order):\n")
  print(basename(h5_files))
  cat("\nProject names (in order):\n")
  print(project_names)
  cat("\n")
  
  object.list <- list()
  
  for (i in seq_along(h5_files)) {
    
    cat("Processing file:", basename(h5_files[i]), "with name:", project_names[i], "\n")
    
    counts <- Seurat::Read10X_h5(h5_files[i])
    
    obj <- Seurat::CreateSeuratObject(
      counts = counts,
      project = as.character(project_names[i])
    )
    
    obj$orig.ident <- project_names[i]
    cat("Set orig.ident to:", unique(obj$orig.ident), "\n")
    
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-|^MT-")
    
    obj <- subset(
      obj, 
      subset = nFeature_RNA > min_feature & 
        nFeature_RNA < max_feature & 
        nCount_RNA > min_count &  
        nCount_RNA < max_count & 
        percent.mt < max_percent_mt
    )
    
    # Rename cells BEFORE adding to list, using project name
    obj <- RenameCells(obj, add.cell.id = project_names[i])
    
    object.list[[i]] <- obj
    
    cat("Remaining cells after QC for", project_names[i], "is", ncol(obj), "cells\n\n")
  }
  
  return(object.list)
}

###############################################################
# SIMPLE CLUSTERIZATION FOR SUBCLUSTERS
###############################################################

cluster_subcluster <- function(subset_obj, output_dir) {
  
  library(clustree)
  subset_obj <- PercentageFeatureSet(subset_obj, pattern = "^MT-", col.name = "percent.mt")
  subset_obj <- SCTransform(subset_obj, vars.to.regress = "percent.mt")
  subset_obj <- RunPCA(subset_obj)
  subset_obj <- RunUMAP(subset_obj, dims = 1:30)
  subset_obj <- FindNeighbors(subset_obj, dims = 1:30)
  subset_obj <- FindClusters(subset_obj, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), algorithm = 1)
  
  
  png(
    filename = file.path(output_dir, ("clustree.png")),
    width = 5,
    height = 5,
    units = "in",
    res = 300
  )
  
  p <- clustree(subset_obj, prefix = "SCT_snn_res.")
  print(p)
  
  dev.off()
  
  return(subset_obj)
  
}

###############################################################
# INTEGRATE SCRNA DATA USING SCANORAMA (PYTHON SCRIPT)
###############################################################

scrna_integrate <- function(object.list, output_folder = "./", dataset_names, 
                            python_script_path = "./integrate_scanorama.py",
                            python_path = NULL) { 
  
  cat("\n##################################\n")
  cat("Normalizing and saving filtered objects...\n")
  cat("##################################\n")
  
  sc <- reticulate::import("scanpy")
  temp_dir <- file.path(output_folder, "preprocessing", "filtered_for_integration")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (i in seq_along(object.list)) {
    obj <- object.list[[i]]
    cat("Processing", dataset_names[i], "...\n")
    
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- JoinLayers(obj)
    data_matrix <- LayerData(obj, assay = "RNA", layer = "data")
    
    gene_names <- rownames(data_matrix)
    
    h5_path <- file.path(temp_dir, paste0(dataset_names[i], "_qc_filtered.h5"))
    DropletUtils::write10xCounts(
      path = h5_path, 
      x = data_matrix, 
      gene.id = gene_names,
      gene.symbol = gene_names,
      version = "3", 
      overwrite = TRUE
    )
    
    cat("  Saved:", ncol(obj), "cells x", nrow(obj), "genes\n")
  }
  
  cat("\n##################################\n")
  cat("Running Scanorama via Python...\n")
  cat("##################################\n")
  
  if (is.null(python_path)) {
    python_cmd <- "python"
  } else {
    python_cmd <- python_path
  }
  
  cmd <- sprintf(
    "%s %s %s %s %s",
    shQuote(python_cmd),
    shQuote(python_script_path),
    shQuote(temp_dir),
    shQuote(file.path(output_folder, "preprocessing")),
    paste(shQuote(dataset_names), collapse = " ")
  )
  
  cat("Running:\n", cmd, "\n\n")
  exit_code <- system(cmd)
  
  if (exit_code != 0) {
    stop("Python script failed with exit code ", exit_code)
  }
  
  cat("\n##################################\n")
  cat("Loading integrated data into R...\n")
  cat("##################################\n")
  
  h5ad_path <- file.path(output_folder, "preprocessing", "integrated_scanorama.h5ad")
  adata <- sc$read_h5ad(h5ad_path)
  
  # Convert expression matrix (stays sparse)
  counts <- reticulate::py_to_r(adata$X)
  if (inherits(counts, "scipy.sparse.base.spmatrix") || 
      grepl("sparse", class(counts)[1], ignore.case = TRUE)) {
    counts <- as(counts, "CsparseMatrix")
  }
  counts <- Matrix::t(counts)
  
  gene_names <- adata$var_names$to_list()
  cell_names <- adata$obs_names$to_list()
  
  rownames(counts) <- gene_names
  colnames(counts) <- cell_names
  
  cat("First 5 genes:", paste(head(gene_names, 5), collapse = ", "), "\n")
  
  # Get metadata with leiden clusters
  metadata <- reticulate::py_to_r(adata$obs)
  rownames(metadata) <- cell_names
  
  if ("leiden" %in% colnames(metadata)) {
    metadata$seurat_clusters <- factor(metadata$leiden)
  }
  
  # Create Seurat object
  integrated_seurat <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata
  )
  
  # Transfer UMAP
  umap_coords <- reticulate::py_to_r(adata$obsm$get("X_umap"))
  colnames(umap_coords) <- paste0("UMAP_", 1:2)
  rownames(umap_coords) <- colnames(integrated_seurat)
  integrated_seurat[["umap"]] <- CreateDimReducObject(
    embeddings = umap_coords,
    key = "UMAP_",
    assay = "RNA"
  )
  
  # Transfer scanorama embedding (use as "pca" equivalent for downstream)
  scanorama_coords <- reticulate::py_to_r(adata$obsm$get("X_scanorama"))
  colnames(scanorama_coords) <- paste0("PC_", 1:ncol(scanorama_coords))
  rownames(scanorama_coords) <- colnames(integrated_seurat)
  integrated_seurat[["pca"]] <- CreateDimReducObject(
    embeddings = scanorama_coords,
    key = "PC_",
    assay = "RNA"
  )
  
  # Transfer neighbor graph
  connectivities <- reticulate::py_to_r(adata$obsp$get("connectivities"))
  if (!is.null(connectivities)) {
    connectivities <- as(connectivities, "CsparseMatrix")
    rownames(connectivities) <- colnames(integrated_seurat)
    colnames(connectivities) <- colnames(integrated_seurat)
    integrated_seurat@graphs$RNA_snn <- as.Graph(connectivities)
  }
  
  Idents(integrated_seurat) <- "seurat_clusters"
  
  # Normalize for downstream (FeaturePlot etc.)
  integrated_seurat <- NormalizeData(integrated_seurat, verbose = FALSE)
  
  cat("\n##################################\n")
  cat("Integration complete!\n")
  cat(sprintf("Cells: %d | Genes: %d | Clusters: %d\n", 
              ncol(integrated_seurat), 
              nrow(integrated_seurat),
              length(unique(integrated_seurat$seurat_clusters))))
  cat("##################################\n\n")
  
  return(integrated_seurat)
}
######################################################
# RUN NEBULOSA ON GENES FOR VISUALIZATION - FIX FOR S5
######################################################

plot_marker_genes <- function(obj, 
                                  genes, 
                                  cluster_col = "seurat_clusters",
                                  reduction = "umap", 
                                  output_dir = "./marker_genes", 
                                  pt_size = 0.3,
                                  outline_size = 0.6,
                                  concavity = 2,
                                  show_labels = FALSE,
                                  eps = 1.5,
                                  min_pts = 5,
                                  outlier_percentile = 0.99) {
  
  library(Nebulosa)
  library(ggplot2)
  library(viridis)
  library(dbscan)
  library(dplyr)
  library(concaveman)
  library(shadowtext)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  available_genes <- rownames(obj)
  
  # Custom color scale for density
  inferno_cols <- viridis::inferno(100)
  custom_cols <- c("grey85", inferno_cols[0:100])
  
  # Get UMAP coordinates and cluster labels
  umap_coords <- as.data.frame(Embeddings(obj, reduction = reduction))
  colnames(umap_coords) <- c("x", "y")
  umap_coords$cluster <- as.factor(obj@meta.data[[cluster_col]])
  
  # Compute boundaries for each cluster
  message("Computing cluster boundaries...")
  boundary_list <- list()
  
  for (clust in levels(umap_coords$cluster)) {
    
    clust_data <- umap_coords %>% filter(cluster == clust)
    
    # Remove outliers based on distance from cluster centroid
    centroid_x <- median(clust_data$x)
    centroid_y <- median(clust_data$y)
    clust_data$dist <- sqrt((clust_data$x - centroid_x)^2 + (clust_data$y - centroid_y)^2)
    
    dist_threshold <- quantile(clust_data$dist, outlier_percentile)
    clust_data <- clust_data %>% filter(dist <= dist_threshold)
    
    clust_matrix <- as.matrix(clust_data[, c("x", "y")])
    
    # Use DBSCAN to find spatially contiguous regions
    db <- dbscan(clust_matrix, eps = eps, minPts = min_pts)
    
    for (sub_clust in unique(db$cluster)) {
      if (sub_clust == 0) next
      
      sub_pts <- clust_matrix[db$cluster == sub_clust, , drop = FALSE]
      
      if (nrow(sub_pts) < 3) next
      
      hull <- concaveman(sub_pts, concavity = concavity)
      
      boundary_list[[length(boundary_list) + 1]] <- data.frame(
        x = hull[, 1],
        y = hull[, 2],
        cluster = clust
      )
    }
  }
  
  boundaries <- bind_rows(boundary_list, .id = "group_id")
  
  # Plot each gene
  for (gene in genes) {
    
    if (!gene %in% available_genes) {
      message(sprintf("Gene not found: %s", gene))
      next
    }
    
    tryCatch({
      
      p <- plot_density(obj, gene, reduction = reduction, size = pt_size) +
        scale_color_gradientn(colors = inferno_cols) +
        ggtitle(gene) +
        theme(
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)
        )
      
      # Add boundaries as dashed white lines
      for (gid in unique(boundaries$group_id)) {
        b <- boundaries %>% filter(group_id == gid)
        
        p <- p + 
          geom_path(
            data = b,
            aes(x = x, y = y),
            color = "white",
            linewidth = outline_size,
            linetype = "dashed",
            inherit.aes = FALSE
          )
      }
      
      if (show_labels) {
        label_coords <- umap_coords %>%
          group_by(cluster) %>%
          summarize(x = median(x), y = median(y), .groups = "drop")
        
        p <- p + 
          geom_shadowtext(
            data = label_coords,
            aes(x = x, y = y, label = cluster),
            color = "white",
            fontface = "bold",
            size = 4,
            inherit.aes = FALSE
          )
      }
      
      ggsave(
        filename = file.path(output_dir, paste0(gene, "_density.png")),
        plot = p,
        width = 5,
        height = 5,
        dpi = 300
      )
      
      message(sprintf("Saved: %s", gene))
      
    }, error = function(e) {
      message(sprintf("Failed to plot %s: %s", gene, e$message))
    })
  }
}

###################################################
# CELL EXTRACTION FUNCTION
###################################################

extract_DA_cells <- function(milo_obj, da_results, alpha = 0.05, 
                             direction = "both", use_pvalue = TRUE) {
  # direction: "up" (logFC > 0), "down" (logFC < 0), or "both"
  
  if(use_pvalue) {
    sig_col <- "PValue"
  } else {
    sig_col <- "SpatialFDR"
  }
  
  # Filter significant neighborhoods
  if(direction == "up") {
    sig_nhoods <- da_results$Nhood[da_results[[sig_col]] < alpha & da_results$logFC > 0]
  } else if(direction == "down") {
    sig_nhoods <- da_results$Nhood[da_results[[sig_col]] < alpha & da_results$logFC < 0]
  } else {
    sig_nhoods <- da_results$Nhood[da_results[[sig_col]] < alpha]
  }
  
  # Extract cells from significant neighborhoods
  cell_barcodes <- c()
  for(i in sig_nhoods) {
    nhood_cells <- colnames(milo_obj)[nhoods(milo_obj)[, i] == 1]
    cell_barcodes <- c(cell_barcodes, nhood_cells)
  }
  
  return(unique(cell_barcodes))
}

##############################################################################
# Simple percentage overview for clusters per hair cycle phase
##############################################################################

get_percentage_clusters <- function(seurat_obj, clusters, phases) {
  data_frame <- data.frame(
    cluster = seurat_obj@meta.data[[clusters]],
    orig.ident = seurat_obj@meta.data[[phases]]
  ) |>
    count(cluster, orig.ident) |>
    group_by(cluster) |>
    mutate(percentage = n / sum(n) * 100) |>
    ungroup()
  return(data_frame)

}

visualize_percentage_clusters <- function(seurat_obj, clusters, phases, output_dir) {
  
  phase_df <- get_percentage_clusters(seurat_obj = seurat_obj, clusters = clusters, phases = phases)
  phase_matrix <- xtabs(percentage ~ orig.ident + cluster, data = phase_df)
  
  colors <- c("#272E6A", "#D51F26", "#1A7D3F")  # Dark blue, Red, Complementary green
  
  png(
    filename = file.path(output_dir, paste0(clusters, "_percentiles.png")),
    width = 10,
    height = 5,
    units = "in",
    res = 300
  )
  
  par(mar = c(7, 4, 4, 8), xpd = TRUE)
  
  bp <- barplot(phase_matrix, 
                beside = FALSE,
                legend = FALSE,
                las = 2,
                col = colors,
                xaxt = "n")
  
  text(x = bp, 
       y = -2,
       labels = colnames(phase_matrix), 
       srt = 45, 
       adj = 1, 
       xpd = TRUE)
  
  legend("topright", 
         inset = c(-0.15, 0),
         legend = rownames(phase_matrix),
         fill = colors,
         title = "orig.ident",
         bty = "n")
  
  dev.off()
}

##############################################################################
# GRAB FEATURE GENES FOR DIFFERENTIAL GENE EXPRESSION ANALYSIS
##############################################################################

get_blacklist_genes <- function(seurat_obj) {
  
  library(GenomicFeatures)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  allGenes <- rownames(seurat_obj)
  
  # Mitochondrial
  mt.genes <- grep(pattern = "^MT-", x = allGenes, value = TRUE)
  
  # Ribosomal
  RPS.genes <- grep(pattern = "^RPS", x = allGenes, value = TRUE)
  RPL.genes <- grep(pattern = "^RPL", x = allGenes, value = TRUE)
  
  # X/Y chromosome genes
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  geneGR <- GenomicFeatures::genes(txdb)
  sexGenesGR <- geneGR[seqnames(geneGR) %in% c("chrY", "chrX")]
  matchedGeneSymbols <- AnnotationDbi::select(org.Hs.eg.db,
                                              keys = sexGenesGR$gene_id,
                                              columns = c("ENTREZID", "SYMBOL"),
                                              keytype = "ENTREZID")
  sexChr.genes <- matchedGeneSymbols$SYMBOL
  
  # Cell cycle genes
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  blacklist.genes <- unique(c(
    mt.genes,
    sexChr.genes,
    s.genes,
    g2m.genes,
    RPS.genes,
    RPL.genes
  ))
  
  # Only return genes present in the data
  blacklist.genes <- blacklist.genes[blacklist.genes %in% allGenes]
  
  message(paste("Blacklisted", length(blacklist.genes), "genes"))
  return(blacklist.genes)
}


##############################################################################
# ENHANCED SPATIAL OVERLAP WITH EXPRESSION THRESHOLDS
##############################################################################

plot_spatial_overlap_thresholded <- function(spatial_obj, gene1, gene2,
                                             sigma = 50,
                                             assay = "SCT",
                                             overlap_method = "geometric_mean",
                                             # NEW THRESHOLD PARAMETERS
                                             min_expr_percentile = 20,  # Only consider top 80% of expression
                                             min_expr_absolute = NULL,   # Optional: absolute cutoff (e.g., 0.5)
                                             min_overlap_percentile = 10, # Only show top 90% of overlap
                                             require_both_genes = TRUE,   # Both genes must pass threshold
                                             output_path = NULL,
                                             pt.size.factor = 1.5,
                                             show_diagnostics = TRUE) {
  
  
  # Get expression data
  assay_data <- GetAssayData(spatial_obj, assay = assay, layer = "data")
  
  # Check genes exist
  if(!gene1 %in% rownames(assay_data)) {
    stop(sprintf("Gene '%s' not found in assay '%s'", gene1, assay))
  }
  if(!gene2 %in% rownames(assay_data)) {
    stop(sprintf("Gene '%s' not found in assay '%s'", gene2, assay))
  }
  
  # Apply Gaussian smoothing
  smoothed <- gaussian_smooth_spatial(spatial_obj, c(gene1, gene2), 
                                      sigma = sigma, assay = assay)
  
  gene1_expr <- as.numeric(smoothed[gene1, ])
  gene2_expr <- as.numeric(smoothed[gene2, ])
  
  # ========================================================================
  # APPLY EXPRESSION THRESHOLDS
  # ========================================================================
  
  # Calculate percentile thresholds for each gene
  gene1_threshold <- quantile(gene1_expr[gene1_expr > 0], 
                              probs = min_expr_percentile/100, 
                              na.rm = TRUE)
  gene2_threshold <- quantile(gene2_expr[gene2_expr > 0], 
                              probs = min_expr_percentile/100, 
                              na.rm = TRUE)
  
  # Override with absolute threshold if provided
  if(!is.null(min_expr_absolute)) {
    gene1_threshold <- max(gene1_threshold, min_expr_absolute)
    gene2_threshold <- max(gene2_threshold, min_expr_absolute)
  }
  
  if(show_diagnostics) {
    message("\n=== EXPRESSION THRESHOLDS ===")
    message(sprintf("%s threshold: %.4f (raw range: %.4f - %.4f)", 
                    gene1, gene1_threshold, min(gene1_expr), max(gene1_expr)))
    message(sprintf("%s threshold: %.4f (raw range: %.4f - %.4f)", 
                    gene2, gene2_threshold, min(gene2_expr), max(gene2_expr)))
    message(sprintf("%s spots above threshold: %d / %d (%.1f%%)", 
                    gene1, sum(gene1_expr > gene1_threshold), 
                    length(gene1_expr),
                    100 * sum(gene1_expr > gene1_threshold) / length(gene1_expr)))
    message(sprintf("%s spots above threshold: %d / %d (%.1f%%)", 
                    gene2, sum(gene2_expr > gene2_threshold), 
                    length(gene2_expr),
                    100 * sum(gene2_expr > gene2_threshold) / length(gene2_expr)))
  }
  
  # Create thresholded expression vectors
  gene1_thresh <- gene1_expr
  gene2_thresh <- gene2_expr
  
  gene1_thresh[gene1_expr < gene1_threshold] <- 0
  gene2_thresh[gene2_expr < gene2_threshold] <- 0
  
  # Calculate overlap using thresholded expression
  if(require_both_genes) {
    # Both genes must be positive for overlap to exist
    overlap_raw <- ifelse(gene1_thresh > 0 & gene2_thresh > 0,
                          switch(overlap_method,
                                 "geometric_mean" = sqrt(gene1_thresh * gene2_thresh),
                                 "minimum" = pmin(gene1_thresh, gene2_thresh),
                                 "product" = gene1_thresh * gene2_thresh,
                                 "average" = (gene1_thresh + gene2_thresh) / 2,
                                 stop("Invalid overlap_method")),
                          0)
  } else {
    # Calculate overlap, then threshold
    overlap_raw <- switch(overlap_method,
                          "geometric_mean" = sqrt(gene1_thresh * gene2_thresh),
                          "minimum" = pmin(gene1_thresh, gene2_thresh),
                          "product" = gene1_thresh * gene2_thresh,
                          "average" = (gene1_thresh + gene2_thresh) / 2,
                          stop("Invalid overlap_method")
    )
  }
  
  # Apply overlap percentile threshold
  nonzero_overlap <- overlap_raw[overlap_raw > 0]
  if(length(nonzero_overlap) > 0) {
    overlap_threshold <- quantile(nonzero_overlap, 
                                  probs = min_overlap_percentile/100, 
                                  na.rm = TRUE)
    overlap_raw[overlap_raw < overlap_threshold] <- 0
  } else {
    overlap_threshold <- 0
  }
  
  # Align with Seurat object
  names(overlap_raw) <- colnames(smoothed)
  overlap_aligned <- overlap_raw[Cells(spatial_obj)]
  overlap_aligned[is.na(overlap_aligned)] <- 0
  
  # Calculate quantiles for color scaling (only from positive values)
  nonzero_vals <- overlap_aligned[overlap_aligned > 0]
  
  if(length(nonzero_vals) == 0) {
    warning("No overlap detected after thresholding. Try lower thresholds.")
    q01 <- 0
    q50 <- 0
    q99 <- max(overlap_aligned)
  } else {
    q01 <- quantile(nonzero_vals, 0.01, na.rm = TRUE)
    q50 <- quantile(nonzero_vals, 0.50, na.rm = TRUE)
    q99 <- quantile(nonzero_vals, 0.99, na.rm = TRUE)
  }
  
  # Diagnostics
  if(show_diagnostics) {
    message("\n=== OVERLAP DIAGNOSTICS ===")
    message(sprintf("Overlap threshold: %.4f", overlap_threshold))
    message(sprintf("Overlap range: %.4f to %.4f", 
                    min(overlap_aligned), max(overlap_aligned)))
    message(sprintf("Spots with overlap > 0: %d / %d (%.1f%%)", 
                    sum(overlap_aligned > 0), 
                    length(overlap_aligned),
                    100 * sum(overlap_aligned > 0) / length(overlap_aligned)))
    message(sprintf("Quantiles - 1%%: %.4f, 50%%: %.4f, 99%%: %.4f", 
                    q01, q50, q99))
    message("============================\n")
  }
  
  # Add to metadata
  spatial_obj$temp_overlap <- overlap_aligned
  
  # Create spatial plot
  # Use discrete colors to emphasize true overlap regions
  p <- SpatialFeaturePlot(
    spatial_obj, 
    features = "temp_overlap",
    image.alpha = 0,
    pt.size.factor = pt.size.factor
  ) +
    scale_fill_gradientn(
      colors = c("grey95", "grey90", "#FFCCCC", "#FF9999", 
                 "#FF5555", "red", "#CC0000", "darkred"),
      values = scales::rescale(c(0, 0.001, q01, q01*1.5, q50*0.8, 
                                 q50*1.2, q99*0.9, q99)),
      limits = c(0, q99),
      na.value = "white",
      name = "Overlap\nScore"
    ) +
    ggtitle(paste0(gene1, " & ", gene2, " spatial overlap"),
            subtitle = sprintf("σ=%dμm | Gene thresh=%d%% | Overlap thresh=%d%%", 
                               sigma, min_expr_percentile, min_overlap_percentile)) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
  
  # Clean up
  spatial_obj$temp_overlap <- NULL
  
  # Save if requested
  if(!is.null(output_path)) {
    ggsave(output_path, p, width = 10, height = 10, dpi = 300, bg = "white")
    message(sprintf("Plot saved to: %s", output_path))
  }
  
  return(list(
    plot = p,
    overlap_scores = overlap_aligned,
    gene1_smoothed = gene1_expr,
    gene2_smoothed = gene2_expr,
    gene1_thresholded = gene1_thresh,
    gene2_thresholded = gene2_thresh,
    thresholds = list(
      gene1 = gene1_threshold,
      gene2 = gene2_threshold,
      overlap = overlap_threshold
    ),
    sigma = sigma,
    overlap_method = overlap_method
  ))
}

# ============================================================================
# SIDE-BY-SIDE COMPARISON: RAW vs THRESHOLDED
# ============================================================================

compare_thresholding <- function(spatial_obj, gene1, gene2,
                                 sigma = 50,
                                 threshold_levels = c(0, 10, 20, 30)) {

  library(patchwork)
  
  plots <- list()
  
  for(thresh in threshold_levels) {
    message(sprintf("\n--- Testing threshold = %d%% ---", thresh))
    
    result <- plot_spatial_overlap_thresholded(
      spatial_obj = spatial_obj,
      gene1 = gene1,
      gene2 = gene2,
      sigma = sigma,
      min_expr_percentile = thresh,
      min_overlap_percentile = max(5, thresh/2),  # Scale overlap threshold too
      show_diagnostics = TRUE,
      pt.size.factor = 1.2
    )
    
    label <- ifelse(thresh == 0, "No threshold", sprintf("%d%% threshold", thresh))
    plots[[as.character(thresh)]] <- result$plot + 
      ggtitle(label)
  }
  
  combined <- wrap_plots(plots, ncol = 2)
  
  return(combined)
}

# ============================================================================
# ADAPTIVE THRESHOLDING BASED ON DISTRIBUTION
# ============================================================================

plot_spatial_overlap_adaptive <- function(spatial_obj, gene1, gene2,
                                          sigma = 50,
                                          assay = "SCT",
                                          overlap_method = "geometric_mean",
                                          sensitivity = "medium",  # "low", "medium", "high"
                                          output_path = NULL,
                                          pt.size.factor = 1.5) {

  # Map sensitivity to threshold values
  threshold_map <- list(
    "low" = list(expr = 30, overlap = 15),
    "medium" = list(expr = 20, overlap = 10),
    "high" = list(expr = 10, overlap = 5)
  )
  
  if(!sensitivity %in% names(threshold_map)) {
    stop("sensitivity must be 'low', 'medium', or 'high'")
  }
  
  thresholds <- threshold_map[[sensitivity]]
  
  message(sprintf("Using %s sensitivity (expr: %d%%, overlap: %d%%)", 
                  sensitivity, thresholds$expr, thresholds$overlap))
  
  result <- plot_spatial_overlap_thresholded(
    spatial_obj = spatial_obj,
    gene1 = gene1,
    gene2 = gene2,
    sigma = sigma,
    assay = assay,
    overlap_method = overlap_method,
    min_expr_percentile = thresholds$expr,
    min_overlap_percentile = thresholds$overlap,
    require_both_genes = TRUE,
    output_path = output_path,
    pt.size.factor = pt.size.factor,
    show_diagnostics = TRUE
  )
  
  return(result)
}


##################################################
# OLD FUNCITONS NO LONGER USED KEPT FOR POSTERITY
##################################################

plot_magic_genes <- function(seurat_obj, 
                             genes,
                             output_folder = "./marker_genes",
                             knn = 10,
                             t = 3,
                             decay = 40,
                             reduction = "umap",
                             pt.size = 0.5,
                             n_pca = 50) {
  
  library(Matrix)
  library(FNN)
  library(ggplot2)
  library(Seurat)
  
  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
  
  seurat_obj <- JoinLayers(seurat_obj)
  
  # Run PCA if not present
  if (!"pca" %in% Reductions(seurat_obj)) {
    message("PCA not found, running PCA...")
    seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = n_pca, verbose = FALSE)
  }
  
  pca <- Embeddings(seurat_obj, "pca")
  if (ncol(pca) < n_pca) {
    n_pca <- ncol(pca)
  }
  pca <- pca[, 1:n_pca]
  
  # Get expression
  expr <- t(GetAssayData(seurat_obj, layer = "data"))
  expr <- expr[, colSums(expr) > 0]
  
  valid_genes <- genes[genes %in% colnames(expr)]
  if (length(valid_genes) < length(genes)) {
    missing <- setdiff(genes, valid_genes)
    message(sprintf("Warning: %d genes not found: %s", 
                    length(missing), paste(missing, collapse = ", ")))
  }
  
  if (length(valid_genes) == 0) {
    stop("No valid genes found in expression matrix")
  }
  
  n_cells <- nrow(pca)
  cat(sprintf("Cells: %d, PCs: %d, Target genes: %d\n", 
              n_cells, n_pca, length(valid_genes)))
  
  # Build kNN graph (k+1 to include self)
  cat("Building kNN graph...\n")
  nn <- get.knn(pca, k = knn)
  
  # Adaptive bandwidth: distance to kth neighbor
  bandwidth <- nn$nn.dist[, knn]
  # Prevent division by zero
  bandwidth[bandwidth == 0] <- min(bandwidth[bandwidth > 0])
  
  cat("Computing affinity matrix...\n")
  
  # Build sparse affinity matrix with MAGIC's alpha-decay kernel
  # Include self-connections (diagonal = 1)
  i_idx <- c(1:n_cells, rep(1:n_cells, each = knn))
  j_idx <- c(1:n_cells, as.vector(nn$nn.index))
  
  # Self-affinities = 1, neighbor affinities use decay kernel
  dists <- as.vector(nn$nn.dist)
  bw_i <- bandwidth[rep(1:n_cells, each = knn)]
  
  # MAGIC kernel: exp(-(d/sigma)^decay) but decay=40 is essentially binary
  # Use simpler Gaussian: exp(-d^2 / sigma^2)
  neighbor_affinities <- exp(-dists^2 / bw_i^2)
  
  all_affinities <- c(rep(1, n_cells), neighbor_affinities)
  
  A <- sparseMatrix(i = i_idx, j = j_idx, x = all_affinities, 
                    dims = c(n_cells, n_cells))
  
  # Symmetrize: use max to preserve connections
  A <- pmax(A, t(A))
  
  # Row-normalize to Markov matrix
  cat("Building Markov matrix...\n")
  row_sums <- rowSums(A)
  D_inv <- Diagonal(x = 1 / row_sums)
  M <- D_inv %*% A
  
  # Verify row sums = 1
  cat(sprintf("Markov matrix row sums: min=%.4f, max=%.4f\n", 
              min(rowSums(M)), max(rowSums(M))))
  
  # Power the matrix
  cat(sprintf("Diffusing (t=%d)...\n", t))
  M_t <- M
  for (i in seq_len(t - 1)) {
    M_t <- M_t %*% M
  }
  
  # Get expression subset
  expr_subset <- as.matrix(expr[, valid_genes, drop = FALSE])
  
  cat(sprintf("Expression range before: min=%.4f, max=%.4f, mean=%.4f\n",
              min(expr_subset), max(expr_subset), mean(expr_subset)))
  
  # Apply diffusion
  expr_smooth <- as.matrix(M_t %*% expr_subset)
  colnames(expr_smooth) <- valid_genes
  
  cat(sprintf("Expression range after:  min=%.4f, max=%.4f, mean=%.4f\n",
              min(expr_smooth), max(expr_smooth), mean(expr_smooth)))
  
  # Get UMAP coordinates
  umap_coords <- Embeddings(seurat_obj, reduction = reduction)
  
  # Plot each gene
  for (gene in valid_genes) {
    
    plot_df <- data.frame(
      UMAP_1 = umap_coords[, 1],
      UMAP_2 = umap_coords[, 2],
      Expression = expr_smooth[, gene]
    )
    
    upper_lim <- quantile(plot_df$Expression, probs = 0.95)
    plot_df$Expression[plot_df$Expression >= upper_lim] <- upper_lim
    
    p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Expression)) +
      geom_point(size = pt.size) +
      scale_color_viridis_c() +
      ggtitle(paste0(gene, " (MAGIC smoothed)")) +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"))
    
    pdf(
      file.path(output_folder, paste0(gene, "_marker.pdf")),
      width = 8,
      height = 6
    )
    print(p)
    dev.off()
  }
  
  cat(sprintf("Plots saved to: %s\n", output_folder))
}


#############################################################################
# SPATIAL DATA SELECTOR
#############################################################################

library(shiny)
library(plotly)
library(Seurat)

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
                              label_name = NULL) {
  
  if (is.null(label_name)) {
    label_name <- paste0(feature, "_Lineage")
  }
  
  expr <- FetchData(obj, vars = feature)[, 1]
  names(expr) <- colnames(obj)
  
  in_region <- colnames(obj) %in% selected_cells
  above_threshold <- expr >= threshold
  
  obj@meta.data[[label_name]] <- in_region & above_threshold
  
  message(paste0(
    "Added '", label_name, "': ",
    sum(obj@meta.data[[label_name]]), " TRUE / ",
    sum(!obj@meta.data[[label_name]]), " FALSE"
  ))
  
  return(obj)
}

#################################################################
# MCLUST FILTERING + QC PLOT
#################################################################
mclust_filter_cells <- function(arrow_file, min_tss = 5, min_frags = 1000, plot_dir = NULL) {
  
  sample_name <- gsub("\\.arrow$", "", basename(arrow_file))
  message(paste0("Processing: ", sample_name))
  
  # Read QC data from Arrow
  cell_names <- h5read(arrow_file, "Metadata/CellNames")
  nFrags <- h5read(arrow_file, "Metadata/nFrags")
  TSSEnrichment <- h5read(arrow_file, "Metadata/TSSEnrichment")
  
  df <- data.frame(
    cellNames = cell_names,
    nFrags = nFrags,
    TSSEnrichment = TSSEnrichment,
    log10_nFrags = log10(nFrags + 1),
    stringsAsFactors = FALSE
  )
  
  # Fit 2-4 Gaussians
  mclust_input <- df[, c("log10_nFrags", "TSSEnrichment")]
  fit <- Mclust(mclust_input, G = 2:4, modelNames = "VVV")
  df$cluster <- fit$classification
  
  # Identify best cluster (highest mean TSS)
  cluster_means <- df %>%
    group_by(cluster) %>%
    summarise(mean_TSS = mean(TSSEnrichment), .groups = "drop")
  best_cluster <- cluster_means$cluster[which.max(cluster_means$mean_TSS)]
  
  # Label cells
  df$passed <- df$cluster == best_cluster & 
    df$TSSEnrichment >= min_tss & 
    df$nFrags >= min_frags
  
  n_passed <- sum(df$passed)
  
  # ADD SAMPLE PREFIX HERE - inside the function
  valid_cells <- paste0(sample_name, "#", df$cellNames[df$passed])
  
  message(paste0("  Retained: ", n_passed, " / ", nrow(df)))
  
  # Generate QC plot
  if (!is.null(plot_dir)) {
    
    p <- ggplot(df, aes(x = log10_nFrags, y = TSSEnrichment)) +
      geom_point(data = df[!df$passed, ], color = "grey80", size = 0.5, alpha = 0.5) +
      stat_density_2d(
        data = df[df$passed, ],
        aes(fill = after_stat(density)),
        geom = "raster",
        contour = FALSE
      ) +
      scale_fill_viridis_c(option = "H", name = "Density") +
      geom_hline(yintercept = min_tss, linetype = "dashed", color = "black") +
      geom_vline(xintercept = log10(min_frags), linetype = "dashed", color = "black") +
      labs(
        title = paste0(sample_name, " ", n_passed, " Cells"),
        x = expression(Log[10] ~ "Unique Fragments"),
        y = "TSS Enrichment"
      ) +
      xlim(3, 5.5) +
      ylim(0, 20) +
      theme_minimal() +
      theme(
        plot.title = element_text(color = "forestgreen", face = "bold", size = 12),
        legend.position = "none"
      )
    
    ggsave(
      filename = file.path(plot_dir, paste0(sample_name, "_mclust_QC.png")),
      plot = p,
      width = 4,
      height = 4,
      dpi = 150
    )
  }
  
  return(valid_cells)
}
