library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(grid)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(DESeq2)
library(tidyr)
library(speckle)


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
                        pattern = "postcellbender_filtered_seurat.h5",
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
    pattern = "postcellbender_filtered_seurat\\.h5$",
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
# TOP EXPRESSED GENES PER CLUSTER
###############################################################

top_expressed_per_cluster <- function(obj,
                                      cluster_col = NULL,
                                      n_genes = 10,
                                      layer = "data",
                                      assay = NULL) {
  
  if (is.null(assay)) assay <- DefaultAssay(obj)
  
  expr <- GetAssayData(obj, assay = assay, layer = layer)
  
  # Exclude mitochondrial and ribosomal genes
  exclude <- grepl("^MT-|^mt-|^RPL|^RPS", rownames(expr))
  expr <- expr[!exclude, , drop = FALSE]
  
  if (is.null(cluster_col)) {
    groups <- Idents(obj)
  } else {
    groups <- obj@meta.data[[cluster_col]]
  }
  names(groups) <- colnames(obj)
  
  results <- lapply(sort(unique(groups)), function(cl) {
    
    cells <- names(groups)[groups == cl]
    mat <- expr[, cells, drop = FALSE]
    
    gene_mean <- Matrix::rowMeans(mat)
    gene_pct  <- Matrix::rowSums(mat > 0) / length(cells) * 100
    
    top_idx <- order(gene_mean, decreasing = TRUE)[seq_len(min(n_genes, length(gene_mean)))]
    
    data.frame(
      cluster   = cl,
      gene      = rownames(mat)[top_idx],
      mean_expr = round(gene_mean[top_idx], 4),
      pct_expr  = round(gene_pct[top_idx], 2),
      rank      = seq_along(top_idx),
      row.names = NULL
    )
  })
  
  do.call(rbind, results)
}

###############################################################
# SIMPLE CLUSTERIZATION FOR SUBCLUSTERS
###############################################################

cluster_subcluster <- function(obj, output_dir, ndims, n_genes = 2000, v_features = 2000, ncells = 5000, kparam = 20) {
  
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(clustree)
  
  set.seed(123)
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
  obj <- PercentageFeatureSet(obj, pattern = "^RP[LS]", col.name = "percent.ribo")
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  obj$CC.Difference <- obj$S.Score - obj$G2M.Score

  obj <- SCTransform(
    obj,
    n_genes = n_genes,
    vars.to.regress = c("percent.mt", "percent.ribo", "CC.Difference"), 
    variable.features.n = v_features,
    ncells = ncells
  )
  
  # Remove blacklist genes from variable features so they don't
  # drive PCA/clustering. Genes remain in the object for downstream use.
  blacklist <- get_blacklist_genes(obj)
  var_feats <- VariableFeatures(obj)
  VariableFeatures(obj) <- setdiff(var_feats, blacklist)
  message(paste("Variable features:", length(var_feats), "->", 
                length(VariableFeatures(obj)), 
                "(removed", length(intersect(var_feats, blacklist)), "blacklisted)"))
  
  obj <- RunPCA(obj)
  obj <- RunTSNE(obj, dims = 1:ndims)
  obj <- RunUMAP(obj, dims = 1:ndims)
  obj <- FindNeighbors(obj, dims = 1:ndims, k.param = kparam)
  obj <- FindClusters(obj, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), algorithm = 4, random.seed = 123)
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  png(
    filename = file.path(output_dir, ("clustree.png")),
    width = 5,
    height = 5,
    units = "in",
    res = 300
  )
  
  p <- clustree(obj, prefix = "SCT_snn_res.")
  print(p)
  
  dev.off()
  
  return(obj)
  
}

subcluster_and_markers <- function(
    obj,
    cluster_name,
    cluster_col = "fine_clust",
    resolution = 0.5,
    dims = 1:20,
    npcs = 30,
    min_pct = 0.1,
    logfc_threshold = 0.5
) {
  
  message("Subsetting cluster: ", cluster_name)
  
  # identify cells
  cells_use <- colnames(obj)[obj@meta.data[[cluster_col]] == cluster_name]
  
  if (length(cells_use) == 0) {
    stop(paste("No cells found for", cluster_name))
  }
  
  # subset object
  sub_obj <- subset(obj, cells = cells_use)
  
  # recluster
  sub_obj <- RunPCA(sub_obj, npcs = npcs, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, dims = dims)
  sub_obj <- FindClusters(sub_obj, resolution = resolution)
  sub_obj <- RunUMAP(sub_obj, dims = dims)
  
  # markers
  sub_markers <- FindAllMarkers(
    sub_obj,
    only.pos = TRUE,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold
  )
  
  return(list(
    sub_obj = sub_obj,
    sub_markers = sub_markers
  ))
}

insert_subclusters <- function(
    obj,
    sub_obj,
    sub_map,
    new_col = "fine_clust"
) {
  
  # get subcluster IDs
  sub_ids <- as.character(Idents(sub_obj))
  
  # convert to biological labels
  new_labels <- sub_map[sub_ids]
  
  if (any(is.na(new_labels))) {
    stop("Some subclusters are missing from sub_map")
  }
  
  # insert back into original object
  obj[[new_col]][colnames(sub_obj), 1] <- new_labels
  
  return(obj)
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
# GROUPED DOTPLOT FOR INITIAL IDENTIFICATION
###################################################

grouped_dotplot <- function(obj,
                            group_def,
                            gene_groups,
                            ident_col   = NULL,
                            cluster_map = NULL,
                            bar_colors  = NULL,
                            col_low     = "steelblue",
                            col_mid     = "white",
                            col_high    = "darkred",
                            dot_range   = c(0.5, 6),
                            base_size   = 11,
                            bar_height  = 0.7) {
  
  ## packages ──────────────────────────────────────────────────────────────
  suppressPackageStartupMessages({
    require(Seurat)
    require(ggplot2)
    require(dplyr)
  })
  
  ## 1. Resolve identities ────────────────────────────────────────────────
  if (!is.null(cluster_map) && is.null(ident_col)) {
    mapping <- unlist(cluster_map)
    obj$grouped_id <- mapping[as.character(Idents(obj))]
    ident_col <- "grouped_id"
  }
  if (is.null(ident_col)) {
    stop("Supply either `ident_col` or `cluster_map`.")
  }
  
  cluster_order <- unlist(group_def, use.names = FALSE)
  obj@meta.data[[ident_col]] <- factor(obj@meta.data[[ident_col]],
                                       levels = cluster_order)
  Idents(obj) <- ident_col
  
  ## 2. Gene order (reversed for y-axis) ──────────────────────────────────
  all_genes <- unlist(gene_groups, use.names = FALSE)
  n_genes   <- length(all_genes)
  
  ## 3. Extract data via DotPlot ──────────────────────────────────────────
  p_tmp  <- DotPlot(obj, features = all_genes) + RotatedAxis()
  dot_df <- p_tmp$data
  dot_df <- dot_df %>%
    dplyr::rename(gene = features.plot, cluster = id)
  
  ## mappings
  gene_to_group <- stack(gene_groups) %>%
    dplyr::rename(gene = values, gene_group = ind)
  cluster_to_group <- stack(group_def) %>%
    dplyr::rename(cluster = values, cluster_group = ind)
  
  dot_df <- dot_df %>%
    dplyr::left_join(gene_to_group,    by = "gene") %>%
    dplyr::left_join(cluster_to_group, by = "cluster")
  
  dot_df$gene    <- factor(dot_df$gene,    levels = rev(all_genes))
  dot_df$cluster <- factor(dot_df$cluster, levels = cluster_order)
  dot_df$gene_group    <- factor(dot_df$gene_group,
                                 levels = names(gene_groups))
  dot_df$cluster_group <- factor(dot_df$cluster_group,
                                 levels = names(group_def))
  
  ## 4. Header bar annotation frame ───────────────────────────────────────
  if (is.null(bar_colors)) {
    palette <- c("#E8A0BF","#A0C4E8","#A8E6CF","#FFD3B6","#D5AAFF",
                 "#FFE0AC","#B5EAD7","#C7CEEA","#FFDAC1","#E2F0CB")
    bar_colors <- setNames(palette[seq_along(group_def)],
                           names(group_def))
  }
  
  x_bar <- data.frame(
    group = names(group_def),
    xmin  = sapply(group_def, function(x)
      which(cluster_order == x[1]) - 0.5),
    xmax  = sapply(group_def, function(x)
      which(cluster_order == x[length(x)]) + 0.5),
    fill  = bar_colors[names(group_def)],
    stringsAsFactors = FALSE
  )
  
  ## 5. Separator positions ───────────────────────────────────────────────
  ## vertical: between cluster groups
  v_breaks <- cumsum(sapply(group_def, length))
  v_breaks <- v_breaks[-length(v_breaks)] + 0.5
  
  ## horizontal: between gene groups
  gene_group_sizes <- sapply(gene_groups, length)
  h_breaks <- cumsum(rev(gene_group_sizes))
  h_breaks <- h_breaks[-length(h_breaks)] + 0.5
  
  ## 6. Build plot ────────────────────────────────────────────────────────
  bar_y_lo <- n_genes + 0.6
  bar_y_hi <- bar_y_lo + bar_height
  
  p <- ggplot(dot_df, aes(x = cluster, y = gene)) +
    geom_point(aes(size = pct.exp, fill = avg.exp.scaled),
               shape = 21, color = "grey30", stroke = 0.3) +
    scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                         midpoint = 0, name = "Scaled\nExpression") +
    scale_size_continuous(range = dot_range, name = "% Expressed",
                          breaks = c(10, 25, 50, 75, 100)) +
    ## separators
    {if (length(v_breaks) > 0)
      geom_vline(xintercept = v_breaks,
                 linewidth = 0.6, color = "grey20")} +
    {if (length(h_breaks) > 0)
      geom_hline(yintercept = h_breaks,
                 linewidth = 0.6, color = "grey20")} +
    ## header bars
    annotate("rect",
             xmin = x_bar$xmin, xmax = x_bar$xmax,
             ymin = bar_y_lo,   ymax = bar_y_hi,
             fill = x_bar$fill, color = "grey30", linewidth = 0.3) +
    annotate("text",
             x = (x_bar$xmin + x_bar$xmax) / 2,
             y = (bar_y_lo + bar_y_hi) / 2,
             label = x_bar$group,
             fontface = "bold", size = 3.5) +
    ## axes & theme
    coord_cartesian(ylim = c(0.5, bar_y_hi + 0.15), clip = "off") +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y  = element_text(face = "italic", size = 8.5),
      axis.title   = element_blank(),
      panel.grid.major = element_line(linewidth = 0.15,
                                      color = "grey85"),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      legend.key.height = unit(0.35, "cm"),
      legend.key.width  = unit(0.35, "cm"),
      legend.title = element_text(size = 8.5),
      legend.text  = element_text(size = 7.5),
      plot.margin  = margin(t = 25, r = 10, b = 5, l = 5)
    ) +
    guides(
      fill = guide_colorbar(order = 1, barheight = unit(3, "cm")),
      size = guide_legend(order = 2,
                          override.aes = list(fill = "grey70"))
    )
  
  return(p)
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
    dplyr::count(cluster, orig.ident) |>
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
  
  allGenes <- rownames(seurat_obj)
  
  # Mitochondrial
  mt.genes <- grep(pattern = "^MT-", x = allGenes, value = TRUE)
  
  # Ribosomal
  RPS.genes <- grep(pattern = "^RPS", x = allGenes, value = TRUE)
  RPL.genes <- grep(pattern = "^RPL", x = allGenes, value = TRUE)
  
  # X/Y chromosome genes - attempt DB lookup, fall back to pattern matching
  sexChr.genes <- tryCatch({
    library(GenomicFeatures)
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    geneGR <- GenomicFeatures::genes(txdb)
    sexGenesGR <- geneGR[seqnames(geneGR) %in% c("chrY", "chrX")]
    matchedGeneSymbols <- AnnotationDbi::select(org.Hs.eg.db,
                                                keys = sexGenesGR$gene_id,
                                                columns = c("ENTREZID", "SYMBOL"),
                                                keytype = "ENTREZID")
    message("Sex chromosome genes retrieved from TxDb")
    matchedGeneSymbols$SYMBOL
  }, error = function(e) {
    message("TxDb unavailable, using pattern-based sex chromosome filtering")
    grep(pattern = "^XIST$|^TSIX$|^RPS4Y|^DDX3Y|^USP9Y|^UTY$|^KDM5D|^EIF1AY|^ZFY$|^SRY$|^NLGN4Y",
         x = allGenes, value = TRUE)
  })
  
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
  
  message(paste("Blacklisted", length(blacklist.genes), "genes",
                "(MT:", length(mt.genes),
                "RPL:", length(RPL.genes),
                "RPS:", length(RPS.genes),
                "Sex:", length(intersect(sexChr.genes, allGenes)),
                "CC:", length(intersect(c(s.genes, g2m.genes), allGenes)), ")"))
  
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

###################################################
# BINOMIAL ENRICHMENT TEST PER CLUSTER
###################################################

#' Tests whether each cluster is statistically enriched in a focal condition
#' versus all remaining conditions, using a binomial distribution adapted
#' from Joost et al., Cell Stem Cell 2020.
#'
#' Supports two experimental layouts:
#'
#'   Two conditions (e.g., sCD83 vs PBS):
#'     Set condition_a and condition_b explicitly, or leave NULL to
#'     auto-detect from the two levels in treatment_col.
#'
#'   Three+ conditions with one-vs-rest (e.g., Anagen vs Catagen/Telogen):
#'     Set focal_condition to the condition of interest. All other
#'     levels are collapsed into a single reference group labelled
#'     "non-<focal_condition>" (or the string in rest_label).
#'     condition_a/condition_b are ignored when focal_condition is set.
#'
#' Under the null, cells from each group distribute across clusters
#' proportionally to their global representation.
#'
#' Two modes for computing expected fraction (p):
#'
#'   use_global_fraction = TRUE (default):
#'     p = total_focal_cells / total_cells
#'     Appropriate when samples were loaded at roughly equal density.
#'
#'   use_global_fraction = FALSE:
#'     Uses the Joost et al. %_a formula which cross-references
#'     sequenced counts with external cell isolation counts (cell_counts).
#'     Required when conditions differ in tissue cellularity.
#'
#' Outlier clusters that are known to be condition-specific and large
#' enough to skew the global baseline can be excluded from p_expected
#' via exclude_clusters. These clusters are still tested and reported.
#'
#' @param obj                 Seurat v5 object
#' @param treatment_col       Column in meta.data with condition labels
#' @param cluster_col         Column in meta.data with cluster identities
#' @param sample_col          Column in meta.data with sample/replicate IDs
#' @param focal_condition     (One-vs-rest mode) The condition of interest.
#'                            All other levels are collapsed into a single
#'                            reference group. When set, condition_a and
#'                            condition_b are ignored.
#' @param rest_label          Label for the collapsed reference group
#'                            (default: "non-<focal_condition>")
#' @param condition_a         (Two-condition mode) Label for condition A.
#'                            Ignored when focal_condition is set.
#' @param condition_b         (Two-condition mode) Label for condition B.
#'                            Ignored when focal_condition is set.
#' @param exclude_clusters    Character vector of cluster names to exclude
#'                            from p_expected computation.
#' @param use_global_fraction Logical (default TRUE). If TRUE, p is the
#'                            simple global fraction of focal cells.
#' @param cell_counts         Named numeric vector of total isolated cells
#'                            per sample (only for use_global_fraction=FALSE).
#' @param alpha               Significance threshold (default 0.001)
#'
#' @return data.frame with per-cluster test results. Attributes store
#'         condition labels and alpha for downstream use by plot_enrichment().

test_cluster_enrichment <- function(obj,
                                    treatment_col       = "treatment",
                                    cluster_col         = "fine_clust",
                                    sample_col          = "orig.ident",
                                    focal_condition     = NULL,
                                    rest_label          = NULL,
                                    condition_a         = NULL,
                                    condition_b         = NULL,
                                    exclude_clusters    = NULL,
                                    use_global_fraction = TRUE,
                                    cell_counts         = NULL,
                                    alpha               = 0.001) {
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(dplyr)
  })
  
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  md <- obj@meta.data
  
  for (col in c(treatment_col, cluster_col, sample_col)) {
    if (!col %in% colnames(md)) {
      stop(paste0("Column '", col, "' not found in meta.data."))
    }
  }
  
  # ------------------------------------------------------------------
  # Resolve condition labels: one-vs-rest OR two-condition mode
  # ------------------------------------------------------------------
  if (!is.null(focal_condition)) {
    
    all_levels <- if (is.factor(md[[treatment_col]])) {
      levels(md[[treatment_col]])
    } else {
      sort(unique(md[[treatment_col]]))
    }
    
    if (!focal_condition %in% all_levels) {
      stop("focal_condition '", focal_condition, "' not found in '",
           treatment_col, "'. Available levels: ",
           paste(all_levels, collapse = ", "))
    }
    
    rest_levels <- setdiff(all_levels, focal_condition)
    if (length(rest_levels) == 0) {
      stop("Only one level found in '", treatment_col,
           "'. Need at least two conditions.")
    }
    
    if (is.null(rest_label)) {
      rest_label <- paste0("non-", focal_condition)
    }
    
    cat("One-vs-rest mode: '", focal_condition, "' vs '", rest_label,
        "' (collapsing: ", paste(rest_levels, collapse = ", "), ")\n",
        sep = "")
    
    # Binarize in meta.data (temporary column)
    md$.binary_condition <- ifelse(
      md[[treatment_col]] == focal_condition,
      focal_condition,
      rest_label
    )
    
    condition_a <- focal_condition
    condition_b <- rest_label
    working_col <- ".binary_condition"
    
  } else {
    
    # Standard two-condition mode
    conditions <- if (is.factor(md[[treatment_col]])) {
      levels(md[[treatment_col]])
    } else {
      sort(unique(md[[treatment_col]]))
    }
    
    if (length(conditions) != 2) {
      stop("Two-condition mode requires exactly 2 levels in '",
           treatment_col, "', found ", length(conditions),
           ": ", paste(conditions, collapse = ", "),
           ". Use focal_condition for one-vs-rest with 3+ levels.")
    }
    
    if (is.null(condition_a)) condition_a <- conditions[1]
    if (is.null(condition_b)) condition_b <- conditions[2]
    
    working_col <- treatment_col
  }
  
  # ------------------------------------------------------------------
  # Baseline: optionally exclude clusters from p_expected computation
  # ------------------------------------------------------------------
  if (!is.null(exclude_clusters)) {
    baseline_idx <- !md[[cluster_col]] %in% exclude_clusters
    md_baseline  <- md[baseline_idx, ]
    n_excluded   <- sum(!baseline_idx)
    cat("Excluding ", length(exclude_clusters), " cluster(s) from baseline: ",
        paste(exclude_clusters, collapse = ", "), "\n", sep = "")
    cat("  -> ", n_excluded, " cells excluded from p_expected computation\n",
        sep = "")
  } else {
    md_baseline <- md
  }
  
  # ------------------------------------------------------------------
  # Compute p_expected
  # ------------------------------------------------------------------
  # Per-sample sequenced cell counts (baseline cells only)
  sample_tab <- md_baseline %>%
    group_by(across(all_of(c(sample_col, working_col)))) %>%
    summarise(n_seq = n(), .groups = "drop")
  
  samples_a <- sample_tab %>% filter(.data[[working_col]] == condition_a)
  samples_b <- sample_tab %>% filter(.data[[working_col]] == condition_b)
  
  S_a <- samples_a$n_seq
  S_b <- samples_b$n_seq
  
  n_total_baseline <- nrow(md_baseline)
  n_a_baseline     <- sum(md_baseline[[working_col]] == condition_a)
  
  if (use_global_fraction) {
    p_expected <- n_a_baseline / n_total_baseline
    cat("Mode: global fraction",
        if (!is.null(exclude_clusters)) "(baseline-corrected)" else "",
        "\n")
  } else {
    if (!is.null(cell_counts)) {
      C_a <- cell_counts[samples_a[[sample_col]]]
      C_b <- cell_counts[samples_b[[sample_col]]]
      if (any(is.na(C_a)) || any(is.na(C_b))) {
        stop("cell_counts names must match all sample IDs in '",
             sample_col, "'.")
      }
    } else {
      C_a <- S_a
      C_b <- S_b
    }
    p_expected <- (sum(S_a) * sum(C_b)) /
      (sum(S_a) * sum(C_b) + sum(S_b) * sum(C_a))
    cat("Mode: Joost et al. cross-sample correction",
        if (!is.null(exclude_clusters)) "(baseline-corrected)" else "",
        "\n")
  }
  
  cat("Condition A: ", condition_a, " | Condition B: ", condition_b, "\n",
      sep = "")
  cat("Baseline cells — ",
      condition_a, ": ", n_a_baseline, " | ",
      condition_b, ": ", n_total_baseline - n_a_baseline,
      " (total: ", n_total_baseline, ")\n", sep = "")
  cat("Per-sample baseline cells — ",
      condition_a, ": ", paste(S_a, collapse = ", "), " | ",
      condition_b, ": ", paste(S_b, collapse = ", "), "\n", sep = "")
  cat("p_expected (fraction ", condition_a, "): ",
      round(p_expected, 4), "\n", sep = "")
  
  # ------------------------------------------------------------------
  # Per-cluster binomial test (ALL clusters, including excluded ones)
  # ------------------------------------------------------------------
  clusters <- sort(unique(md[[cluster_col]]))
  
  results <- do.call(rbind, lapply(clusters, function(cl) {
    
    idx <- md[[cluster_col]] == cl
    n   <- sum(idx)
    n_a <- sum(md[[working_col]][idx] == condition_a)
    n_b <- n - n_a
    
    # P(X >= n_a): if small, cluster enriched in condition_a
    pval_a <- pbinom(n_a - 1, size = n, prob = p_expected,
                     lower.tail = FALSE)
    # P(X <= n_a): if small, cluster enriched in condition_b
    pval_b <- pbinom(n_a, size = n, prob = p_expected)
    
    enrichment <- "none"
    if (pval_a < alpha) enrichment <- condition_a
    if (pval_b < alpha) enrichment <- condition_b
    if (pval_a < alpha && pval_b < alpha) {
      enrichment <- ifelse(pval_a < pval_b, condition_a, condition_b)
    }
    
    excluded <- cl %in% exclude_clusters
    
    data.frame(
      cluster             = cl,
      n_total             = n,
      n_condA             = n_a,
      n_condB             = n_b,
      frac_condA          = n_a / n,
      p_expected          = p_expected,
      pval_enriched_condA = pval_a,
      pval_enriched_condB = pval_b,
      enrichment          = enrichment,
      excluded_from_baseline = excluded,
      stringsAsFactors    = FALSE
    )
  }))
  
  rownames(results) <- NULL
  attr(results, "condition_a")      <- condition_a
  attr(results, "condition_b")      <- condition_b
  attr(results, "alpha")            <- alpha
  attr(results, "focal_condition")  <- focal_condition
  attr(results, "rest_label")       <- if (!is.null(focal_condition)) rest_label else NULL
  
  n_a_enr <- sum(results$enrichment == condition_a)
  n_b_enr <- sum(results$enrichment == condition_b)
  n_none  <- sum(results$enrichment == "none")
  cat(n_a_enr, " cluster(s) enriched in ", condition_a, " | ",
      n_b_enr, " cluster(s) enriched in ", condition_b, " | ",
      n_none, " not enriched (alpha = ", alpha, ")\n", sep = "")
  
  if (!is.null(exclude_clusters)) {
    excl_res <- results %>% filter(excluded_from_baseline)
    cat("Excluded cluster results:\n")
    for (i in seq_len(nrow(excl_res))) {
      cat("  ", excl_res$cluster[i], ": ",
          excl_res$n_condA[i], " ", condition_a, " / ",
          excl_res$n_condB[i], " ", condition_b,
          " -> enriched in ", excl_res$enrichment[i], "\n", sep = "")
    }
  }
  
  return(results)
}

###################################################
# PLOT ENRICHMENT RESULTS PER COMPARTMENT
###################################################

#' Takes the output of test_cluster_enrichment() (run once on the
#' full object), attaches mapping_cell_type from the Seurat object,
#' and produces one separate plot per compartment.
#'
#' @param results      Output from test_cluster_enrichment()
#' @param obj          Seurat object (to look up mapping_cell_type per cluster)
#' @param cluster_col  Column in obj meta.data matching results$cluster
#' @param group_col    Column in obj meta.data for compartment grouping
#' @param colors       Named vector of bar colors per compartment
#' @param alpha        Significance threshold. NULL = pulled from results.
#' @param max_score    Y-axis cap in -log10 scale
#' @param focal_color  Shading color for focal enrichment
#' @param ref_color    Shading color for reference enrichment
#' @param base_size    Base font size
#' @param point_size   Lollipop point size
#' @param output_dir   Directory to save plots. NULL = no saving.
#' @param width        Plot width. NULL = auto-scaled to cluster count.
#' @param height       Plot height
#'
#' @return Named list of ggplot objects (one per compartment), invisibly.

plot_enrichment_by_group <- function(results,
                                     obj,
                                     cluster_col  = "fine_clust",
                                     group_col    = "mapping_cell_type",
                                     colors       = NULL,
                                     alpha        = NULL,
                                     max_score    = 20,
                                     focal_color  = "#4CAF50",
                                     ref_color    = "#E53935",
                                     base_size    = 11,
                                     point_size   = 2.5,
                                     output_dir   = NULL,
                                     width        = NULL,
                                     height       = 5) {
  
  suppressPackageStartupMessages({
    require(ggplot2)
    require(dplyr)
  })
  
  if (is.null(colors)) {
    colors <- c(
      "FOL"  = "#208A42",
      "EPI"  = "#D51F26",
      "FIB"  = "#272E6A",
      "ENDO" = "#89288F",
      "IMM"  = "#b35900",
      "NC"   = "#808080"
    )
  }
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Map each cluster to its compartment via majority vote
  md <- obj@meta.data
  cluster_to_group <- md %>%
    dplyr::group_by(.data[[cluster_col]], .data[[group_col]]) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(.data[[cluster_col]]) %>%
    dplyr::slice_max(n_cells, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(all_of(c(cluster_col, group_col)))
  colnames(cluster_to_group) <- c("cluster", "group")
  
  results <- results %>%
    left_join(cluster_to_group, by = "cluster")
  
  # Order compartments by colors vector
  ordered <- intersect(names(colors), unique(results$group))
  remaining <- setdiff(unique(results$group), ordered)
  compartments <- c(ordered, remaining)
  compartments <- compartments[!is.na(compartments)]
  
  plots <- list()
  
  for (comp in compartments) {
    
    df <- results %>% filter(group == comp)
    if (nrow(df) == 0) next
    
    p <- plot_enrichment_single(
      df,
      max_score   = max_score,
      focal_color = focal_color,
      ref_color   = ref_color,
      base_size   = base_size,
      point_size  = point_size,
      alpha       = alpha,
      condition_a = attr(results, "condition_a"),
      condition_b = attr(results, "condition_b")
    )
    
    plots[[comp]] <- p
    
    if (!is.null(output_dir)) {
      auto_width <- if (is.null(width)) max(3, nrow(df) * 0.35 + 1.5) else width
      fpath <- file.path(output_dir, paste0("enrichment_", comp, ".png"))
      ggsave(fpath, p, width = auto_width, height = height, dpi = 300)
      cat("Saved:", fpath, "\n")
    }
  }
  
  invisible(plots)
}


###################################################
# SINGLE-COMPARTMENT ENRICHMENT PLOT
###################################################

plot_enrichment_single <- function(results,
                                   max_score   = 20,
                                   focal_color = "#4CAF50",
                                   ref_color   = "#E53935",
                                   base_size   = 11,
                                   point_size  = 2.5,
                                   alpha       = NULL,
                                   condition_a = NULL,
                                   condition_b = NULL) {
  
  suppressPackageStartupMessages({
    require(ggplot2)
    require(dplyr)
  })
  
  if (is.null(condition_a)) condition_a <- attr(results, "condition_a")
  if (is.null(condition_b)) condition_b <- attr(results, "condition_b")
  if (is.null(alpha)) alpha <- attr(results, "alpha") %||% 0.001
  thresh <- -log10(alpha)
  
  results <- results %>%
    mutate(
      score = case_when(
        frac_condA >= p_expected ~ -log10(pmax(pval_enriched_condA, 1e-300)),
        TRUE                    ~  log10(pmax(pval_enriched_condB, 1e-300))
      ),
      score_capped = pmax(pmin(score, max_score), -max_score),
      point_color = case_when(
        enrichment == condition_a ~ "focal",
        enrichment == condition_b ~ "ref",
        TRUE                     ~ "ns"
      )
    ) %>%
    arrange(cluster)
  
  results$cluster <- factor(results$cluster, levels = unique(results$cluster))
  n_clust <- nrow(results)
  
  y_breaks <- sort(unique(c(-max_score, -15, -10, -thresh,
                            0, thresh, 10, 15, max_score)))
  y_breaks <- y_breaks[y_breaks >= -max_score & y_breaks <= max_score]
  
  y_labels <- sapply(y_breaks, function(v) {
    if (v == 0) return(expression(1))
    as.expression(bquote(10^{-.(abs(round(v)))}))
  })
  
  ggplot(results, aes(x = cluster, y = score_capped)) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = 0, ymax = thresh,
             fill = focal_color, alpha = 0.06) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = thresh, ymax = max_score,
             fill = focal_color, alpha = 0.14) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -thresh, ymax = 0,
             fill = ref_color, alpha = 0.06) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -max_score, ymax = -thresh,
             fill = ref_color, alpha = 0.14) +
    geom_hline(yintercept = thresh, linetype = "dashed",
               color = "grey40", linewidth = 0.35) +
    geom_hline(yintercept = -thresh, linetype = "dashed",
               color = "grey40", linewidth = 0.35) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.25) +
    geom_segment(aes(x = cluster, xend = cluster,
                     y = 0, yend = score_capped),
                 color = "grey50", linewidth = 0.3) +
    geom_point(aes(color = point_color), size = point_size) +
    scale_color_manual(
      values = c("focal" = focal_color, "ref" = ref_color, "ns" = "black"),
      guide  = "none"
    ) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels,
                       expand = expansion(mult = 0.03)) +
    scale_x_discrete(expand = expansion(add = 0.15)) +
    coord_cartesian(ylim = c(-max_score, max_score)) +
    labs(x = NULL, y = "P-value [binomial]") +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                        size = 9, face = "bold"),
      axis.text.y        = element_text(size = 8),
      axis.title.y       = element_text(size = 9, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.15),
      panel.border       = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      panel.background   = element_blank(),
      plot.margin        = margin(t = 10, r = 10, b = 10, l = 10)
    )
}

###############################################################
# SUBCLUSTER FOR MILO PIPELINE (LIGHTWEIGHT)
# Use when input is already annotated with SCT layer
###############################################################

subcluster_for_milo <- function(obj, output_dir, dims = 1:30, 
                                resolution = 0.5, 
                                n_genes = 2000,
                                v_features = 2000, 
                                ncells = 5000) {
  
  library(Seurat)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(glmGamPoi)
  
  # Cell cycle genes
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # QC metrics
  obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
  obj <- PercentageFeatureSet(obj, pattern = "^RP[LS]", col.name = "percent.ribo")
  
  # Normalize and score cell cycle
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  obj$CC.Difference <- obj$S.Score - obj$G2M.Score
  
  # SCTransform with regression
  obj <- SCTransform(
    obj,
    n_genes = n_genes,
    variable.features.n = v_features,
    ncells = ncells,
    vars.to.regress = c("percent.mt", "percent.ribo", "CC.Difference"),
    method = "glmGamPoi",
    seed.use = 123,
    verbose = FALSE
  )
  
  # Remove blacklist genes from variable features
  if (exists("get_blacklist_genes", mode = "function")) {
    blacklist <- get_blacklist_genes(obj)
    var_feats <- VariableFeatures(obj)
    VariableFeatures(obj) <- setdiff(var_feats, blacklist)
    message(paste("Variable features:", length(var_feats), "->",
                  length(VariableFeatures(obj)),
                  "(removed", length(intersect(var_feats, blacklist)), "blacklisted)"))
  }
  
  # PCA, UMAP, clustering
  obj <- RunPCA(obj, features = VariableFeatures(obj), seed.use = 123, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims, seed.use = 123, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, algorithm = 1, verbose = FALSE)
  return(obj)
}
###############################################################
# ADAPTIVE MILO PARAMETERS BASED ON DATASET SIZE
###############################################################

get_milo_params <- function(n_cells, n_samples) {
  
  # Select k, prop for Milo based on dataset size.
  
  # Returns a named list with k and prop.
  # Note: d (PCA dims) is set separately via subcluster_params in main script.
  
  
  if (n_cells < 5000) {
    list(k = 20, prop = 0.2)
  } else if (n_cells < 10000) {
    list(k = 30, prop = 0.1)
  } else if (n_cells < 20000) {
    list(k = 30, prop = 0.1)
  } else {
    list(k = 40, prop = 0.1)
  }
}

###############################################################
# RUN MILO PIPELINE
###############################################################

run_milo_pipeline <- function(obj_sub,
                              k             = 50,
                              d             = 50,
                              prop          = 0.1,
                              sample_col    = "orig.ident",
                              treatment_col = "treatment",
                              pca_name      = "pca",
                              umap_name     = "umap") {
  
  DefaultAssay(obj_sub) <- "RNA"
  obj.sce <- as.SingleCellExperiment(obj_sub)
  
  # MiloR needs a logcounts assay. The Seurat→SCE conversion should
  # produce one from the RNA "data" layer, but if it's missing, add it.
  if (!"logcounts" %in% assayNames(obj.sce)) {
    cat("  logcounts missing — running logNormCounts...\n")
    obj.sce <- scuttle::logNormCounts(obj.sce)
  }
  
  # Inject Seurat reductions
  seurat_pca  <- Embeddings(obj_sub, reduction = pca_name)
  seurat_umap <- Embeddings(obj_sub, reduction = umap_name)
  
  # Cap d to available PCA components
  if (ncol(seurat_pca) < d) {
    d <- ncol(seurat_pca)
    cat("  PCA has", ncol(seurat_pca), "components — using d =", d, "\n")
  }
  
  reducedDim(obj.sce, "PCA")  <- seurat_pca[, 1:d]
  reducedDim(obj.sce, "UMAP") <- seurat_umap
  
  # Build Milo
  cat("  Building Milo object...\n")
  obj.milo <- Milo(obj.sce)
  obj.milo <- buildGraph(obj.milo, k = k, d = d, reduced.dim = "PCA")
  obj.milo <- makeNhoods(obj.milo, prop = prop, k = k, d = d,
                         refined = TRUE, reduced_dims = "PCA")
  
  cat("  Neighborhoods:", ncol(nhoods(obj.milo)),
      " | Mean size:", round(mean(colSums(nhoods(obj.milo))), 1), "\n")
  
  # Count + distances
  obj.milo <- countCells(obj.milo,
                         meta.data = data.frame(colData(obj.milo)),
                         samples   = sample_col)
  obj.milo <- calcNhoodDistance(obj.milo, d = d, reduced.dim = "PCA")
  
  # Design matrix
  sample_meta <- colData(obj.milo) %>%
    as.data.frame() %>%
    dplyr::select(all_of(c(sample_col, treatment_col))) %>%
    distinct()
  
  design_matrix <- data.frame(
    Sample    = sample_meta[[sample_col]],
    Treatment = sample_meta[[treatment_col]]
  )
  rownames(design_matrix) <- design_matrix$Sample
  design_matrix <- design_matrix[colnames(nhoodCounts(obj.milo)), ]
  
  cat("  Design matrix:\n")
  print(design_matrix)
  
  # DA testing
  cat("  Testing for DA...\n")
  da_results <- testNhoods(obj.milo,
                           norm.method = "RLE",
                           design    = ~ Treatment,
                           design.df = design_matrix)
  
  cat("  Sig neighborhoods (P < 0.05):",
      sum(da_results$PValue < 0.05, na.rm = TRUE),
      "| Up:", sum(da_results$PValue < 0.05 & da_results$logFC > 0, na.rm = TRUE),
      "| Down:", sum(da_results$PValue < 0.05 & da_results$logFC < 0, na.rm = TRUE), "\n")
  
  # Nhood graph + annotate
  obj.milo <- buildNhoodGraph(obj.milo)
  
  if ("fine_clust" %in% colnames(colData(obj.milo))) {
    da_results <- annotateNhoods(obj.milo, da_results,
                                 coldata_col = "fine_clust")
  }
  
  return(list(milo = obj.milo, da_results = da_results))
}

################################################################################
# MILO DA NEIGHBOURHOOD GRAPH
################################################################################

#' layout="UMAP" pulls from the same reducedDim injected from Seurat.
plot_milo_da <- function(milo_obj, da_results, alpha = 0.05, title = NULL) {
  
  p <- plotNhoodGraphDA(milo_obj, da_results, layout = "UMAP", alpha = alpha)
  if (!is.null(title)) p <- p + ggtitle(title)
  
  p + theme_void(base_size = 11) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
}


################################################################################
# PLOT SPLIT UMAP FOR MILO PIPELINE
################################################################################
#'
#' Uses the same UMAP reduction that was injected into the Milo object.
#'
#' @param obj_sub   Seurat object (already has UMAP from cluster_subcluster)
#' @param split_by  Facet column (default "treatment")
#' @param color_by  Color column (default "fine_clust")
#' @param reduction Reduction name (default "umap")
#' @param pt_size   Point size (default 0.4)
#' @param title     Plot title
#' @return ggplot
plot_split_umap <- function(obj_sub,
                            split_by  = "treatment",
                            color_by  = "fine_clust",
                            reduction = "umap",
                            pt_size   = 0.4,
                            title     = NULL) {
  
  emb <- Embeddings(obj_sub, reduction = reduction)
  df  <- data.frame(
    UMAP_1    = emb[, 1],
    UMAP_2    = emb[, 2],
    split_var = obj_sub@meta.data[[split_by]],
    color_var = obj_sub@meta.data[[color_by]],
    stringsAsFactors = FALSE
  )
  
  set.seed(42)
  df <- df[sample(nrow(df)), ]
  
  if (is.null(title)) title <- paste0(ncol(obj_sub), " cells")
  
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_var)) +
    geom_point(size = pt_size, stroke = 0) +
    facet_wrap(~ split_var, ncol = 2) +
    labs(title = title, color = color_by) +
    theme_void(base_size = 11) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
      strip.text      = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.title    = element_text(size = 9),
      legend.text     = element_text(size = 8),
      legend.key.size = unit(0.35, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))
}

############################################################################
# PLOT ORDERED DA BEESWARM (CLEANED UP)
############################################################################

plot_da_beeswarm_ordered <- function(da_results, 
                                     group_by = "fine_clust",
                                     title = NULL,
                                     alpha = 0.05,
                                     use_pvalue_as_fdr = FALSE,
                                     pt_size = 1.5) {
  
  library(ggplot2)
  library(ggbeeswarm)
  
  df <- as.data.frame(da_results)
  
  # Determine significance column
  if (use_pvalue_as_fdr) {
    df$sig <- df$PValue < alpha
  } else {
    df$sig <- df$SpatialFDR < alpha
  }
  
  # Reorder: non-significant first, significant last (plotted on top)
  df <- df[order(df$sig, decreasing = FALSE), ]
  
  # Handle missing group column
  if (!group_by %in% colnames(df)) {
    warning(paste0("Column '", group_by, "' not found. Using NhoodGroup if available."))
    if ("NhoodGroup" %in% colnames(df)) {
      group_by <- "NhoodGroup"
    } else {
      df$group <- "All"
      group_by <- "group"
    }
  }
  
  p <- ggplot(df, aes(x = .data[[group_by]], y = logFC, color = sig)) +
    geom_quasirandom(size = pt_size, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "grey70"),
      labels = c("TRUE" = paste0("p < ", alpha), "FALSE" = "NS"),
      name = NULL
    ) +
    labs(
      x = NULL,
      y = "Log Fold Change",
      title = title
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x    = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y    = element_text(size = 10, face = "bold"),
      axis.title.y   = element_text(size = 11, face = "bold"),
      plot.title     = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "top",
      legend.text    = element_text(size = 10, face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  return(p)
}

################################################################################
# PLOT SPLIT UMAP FOR MILO PIPELINE (PRETTIER VERSION)
################################################################################

plot_split_dimred <- function(obj_sub,
                              split_by  = "treatment",
                              color_by  = "fine_clust",
                              reduction = "umap",
                              pt_size   = 0.5,
                              title     = NULL,
                              cmap      = NULL,
                              useRaster = TRUE) {
  
  library(ggplot2)
  library(ggrastr)
  
  # Default stallion palette
  if (is.null(cmap)) {
    cmap <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B",
              "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC","10"="#90D5E4",
              "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8",
              "16"="#6E4B9E","17"="#0C727C","18"="#7E1416","9"="#D8A767","20"="#3D3D3D")
  }
  
  emb <- Embeddings(obj_sub, reduction = reduction)
  df  <- data.frame(
    UMAP_1    = emb[, 1],
    UMAP_2    = emb[, 2],
    split_var = obj_sub@meta.data[[split_by]],
    color_var = obj_sub@meta.data[[color_by]],
    stringsAsFactors = FALSE
  )
  
  # Shuffle for random plotting order
  set.seed(42)
  df <- df[sample(nrow(df)), ]
  
  # Expand color palette if needed
  n_colors <- length(unique(df$color_var))
  if (n_colors > length(cmap)) {
    cmap <- colorRampPalette(cmap)(n_colors)
  } else {
    names(cmap) <- NULL
    cmap <- cmap[1:n_colors]
  }
  
  # Build plot
  if (useRaster) {
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_var)) +
      geom_point_rast(size = pt_size)
  } else {
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_var)) +
      geom_point(size = pt_size)
  }
  
  p <- p +
    facet_wrap(~ split_var, ncol = 2) +
    scale_color_manual(values = cmap, name = NULL, na.value = "grey35") +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_classic(base_size = 11) +
    theme(
      axis.ticks       = element_blank(),
      axis.text        = element_blank(),
      aspect.ratio     = 1,
      strip.text       = element_text(face = "bold", size = 12),
      strip.background = element_blank(),
      legend.position  = "right",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 10, face = "bold"),
      legend.key.size  = unit(0.4, "cm"),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  # Title: default to cell count, or custom, or NULL
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

################################################################################
# PLOT SINGLE UMAP (PRETTIER VERSION)
################################################################################

plot_single_dimred <- function(obj_sub,
                               color_by  = "fine_clust",
                               reduction = "umap",
                               pt_size   = 0.5,
                               title     = NULL,
                               cmap      = NULL,
                               useRaster = TRUE) {
  
  library(ggplot2)
  library(ggrastr)
  
  # Default stallion palette
  if (is.null(cmap)) {
    cmap <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B",
              "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","9"="#E6C2DC","10"="#90D5E4",
              "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8",
              "16"="#6E4B9E","17"="#0C727C","18"="#7E1416","19"="#D8A767","20"="#3D3D3D")
  }
  
  emb <- Embeddings(obj_sub, reduction = reduction)
  df  <- data.frame(
    UMAP_1    = emb[, 1],
    UMAP_2    = emb[, 2],
    color_var = obj_sub@meta.data[[color_by]],
    stringsAsFactors = FALSE
  )
  
  # Shuffle for random plotting order
  set.seed(42)
  df <- df[sample(nrow(df)), ]
  
  # Expand color palette if needed
  n_colors <- length(unique(df$color_var))
  if (n_colors > length(cmap)) {
    cmap <- colorRampPalette(cmap)(n_colors)
  } else {
    names(cmap) <- NULL
    cmap <- cmap[1:n_colors]
  }
  
  # Build plot
  if (useRaster) {
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_var)) +
      geom_point_rast(size = pt_size)
  } else {
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_var)) +
      geom_point(size = pt_size)
  }
  
  p <- p +
    scale_color_manual(values = cmap, name = NULL, na.value = "grey35") +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_classic(base_size = 11) +
    theme(
      axis.ticks       = element_blank(),
      axis.text        = element_blank(),
      aspect.ratio     = 1,
      legend.position  = "right",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 10, face = "bold"),
      legend.key.size  = unit(0.4, "cm"),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

################################################################################
# MILO WRAPPER PLOT FUNCTIONS
################################################################################

plot_nhood_size_hist <- function(milo_obj, title = "Neighbourhood Sizes") {
  plotNhoodSizeHist(milo_obj) +
    ggtitle(title) +
    theme_classic()
}

plot_pvalue_hist <- function(da_results, title = "P-value Distribution") {
  ggplot(da_results, aes(PValue)) +
    geom_histogram(bins = 50, fill = "#272E6A", color = "white") +
    theme_classic() +
    labs(title = title, x = "P-value", y = "Count")
}

plot_volcano <- function(da_results, pval_thresh = 0.05, title = "Volcano Plot") {
  ggplot(da_results, aes(logFC, -log10(PValue))) +
    geom_point(aes(color = PValue < pval_thresh), size = 2) +
    scale_color_manual(values = c("grey60", "#D51F26"),
                       labels = c("NS", paste0("P < ", pval_thresh)),
                       name = "Significance") +
    geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", color = "red") +
    theme_classic() +
    labs(title = title, x = "Log Fold Change", y = "-log10(P-value)")
}

summarise_da_by_cluster <- function(da_results, cluster_col = "fine_clust") {
  if (!cluster_col %in% colnames(da_results)) return(NULL)
  
  da_results %>%
    group_by(.data[[cluster_col]]) %>%
    summarise(
      n_nhoods     = n(),
      n_sig_p0.05  = sum(PValue < 0.05, na.rm = TRUE),
      n_up         = sum(PValue < 0.05 & logFC > 0, na.rm = TRUE),
      n_down       = sum(PValue < 0.05 & logFC < 0, na.rm = TRUE),
      mean_logFC   = mean(logFC, na.rm = TRUE),
      median_logFC = median(logFC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sig_p0.05))
}

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
################################################################################
# PLOT NHOOD GRAPH WITH P-VALUE FILTERING (ADJUSTABLE NODE SIZE)
################################################################################

plot_nhood_umap <- function(milo_obj, 
                            da_results,
                            alpha = 0.05,
                            use_pvalue = TRUE,
                            layout = "UMAP",
                            node_size_range = c(0.5, 4),
                            edge_width = 0.2,
                            title = NULL) {
  
  library(ggplot2)
  library(igraph)
  library(ggraph)
  
  # Get significance column
  if (use_pvalue) {
    da_results$sig <- da_results$PValue < alpha
  } else {
    da_results$sig <- da_results$SpatialFDR < alpha
  }
  
  # Set non-significant logFC to 0 for coloring
  da_results$logFC_plot <- ifelse(da_results$sig, da_results$logFC, 0)
  
  # Get nhood graph
  nh_graph <- nhoodGraph(milo_obj)
  
  # Get layout coordinates from UMAP
  umap_coords <- reducedDim(milo_obj, layout)
  
  # Calculate nhood centroids
  nh_pos <- do.call(rbind, lapply(seq_len(ncol(nhoods(milo_obj))), function(i) {
    cells_in_nhood <- which(nhoods(milo_obj)[, i] == 1)
    colMeans(umap_coords[cells_in_nhood, , drop = FALSE])
  }))
  
  # Build graph layout
  layout_df <- data.frame(
    x = nh_pos[, 1],
    y = nh_pos[, 2]
  )
  
  # Node data
  node_data <- data.frame(
    nhood = seq_len(ncol(nhoods(milo_obj))),
    logFC = da_results$logFC_plot,
    sig = da_results$sig,
    size = colSums(nhoods(milo_obj))
  )
  
  # Create ggraph
  g <- ggraph(nh_graph, layout = as.matrix(layout_df)) +
    geom_edge_link0(edge_colour = "grey80", edge_width = edge_width) +
    geom_node_point(aes(fill = node_data$logFC, size = node_data$size), 
                    shape = 21, stroke = 0.1, color = "grey30") +
    scale_fill_gradient2(low = "#F8766C", mid = "white", high = "#00BEC4",
                         midpoint = 0, name = "logFC",
                         limits = c(-max(abs(da_results$logFC), na.rm = TRUE),
                                    max(abs(da_results$logFC), na.rm = TRUE))) +
    scale_size_continuous(range = node_size_range, name = "Nhood size") +
    theme_void(base_size = 11) +
    theme(
      legend.position  = "right",
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 9, face = "bold"),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13)
    )
  
  if (!is.null(title)) {
    g <- g + ggtitle(title)
  }
  
  return(g)
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
  
 # Add sample prefix
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


#################################################################
# BOR_THEME FOR PLOTS
#################################################################

# Clean theme for plots

theme_BOR <- function(border = TRUE) {
  theme_classic() +
    theme(
      panel.border = if(border) element_rect(color = "black", fill = NA, size = 0.5) else element_blank(),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

#################################################################
# CALCULTING AVERAGE EXPRESSION AND PERCENT EXPRESSION PER GROUP
#################################################################

#' For each gene and each group (cluster), computes:
#' - avgExpr: mean expression across cells in that group (optionally normalized)
#' - pctExpr: percentage of cells with expression > 0
#'
#' @param mat Expression/score matrix (genes x cells)
#' @param groups Vector of group assignments per cell
#' @param feature_normalize Scale each gene to 0-1 range across groups
#' @param min_pct Minimum percent expressing to include
#' 
avg_prct_expression <- function(mat, groups, feature_normalize = TRUE, min_pct = 0) {
  
  # Ensure groups align with matrix columns
  groups <- as.character(groups)
  unique_groups <- unique(groups) %>% sort()
  
  results <- lapply(unique_groups, function(grp) {
    
    # Cells in this group
    cells_in_grp <- which(groups == grp)
    
    if (length(cells_in_grp) == 0) return(NULL)
    
    # Subset matrix to this group
    sub_mat <- mat[, cells_in_grp, drop = FALSE]
    
    # Average expression per gene
    avg_expr <- rowMeans(sub_mat, na.rm = TRUE)
    
    # Percent of cells expressing (value > 0)
    pct_expr <- rowSums(sub_mat > 0, na.rm = TRUE) / ncol(sub_mat) * 100
    
    data.frame(
      feature = rownames(mat),
      grp = grp,
      avgExpr = avg_expr,
      pctExpr = pct_expr,
      stringsAsFactors = FALSE
    )
  })
  
  df <- do.call(rbind, results)
  
  # Feature normalization: scale each gene's avgExpr to 0-1 range
  # This allows comparison across genes with different expression magnitudes
  if (feature_normalize) {
    df <- df %>%
      group_by(feature) %>%
      mutate(avgExpr = (avgExpr - min(avgExpr)) / (max(avgExpr) - min(avgExpr) + 1e-9)) %>%
      ungroup()
  }
  
  # Filter by minimum percent expressing
  if (min_pct > 0) {
    df <- df %>% filter(pctExpr >= min_pct)
  }
  
  return(as.data.frame(df))
}

#################################################################
# DOTPLOT TO VISUALIZE GENE MARKER EXPRESSION
#################################################################

plot_dotplot <- function(df, xcol, ycol, color_col, size_col, 
                    xorder = NULL, yorder = NULL, cmap = NULL, 
                    color_label = NULL, size_label = NULL, 
                    aspectRatio = NULL, sizeLims = NULL, colorLims = NULL) {
  
  # Set order of axes (important for interpretability)
  if (is.null(xorder)) xorder <- unique(df[, xcol]) %>% sort()
  if (is.null(yorder)) yorder <- unique(df[, ycol]) %>% sort()
  if (is.null(aspectRatio)) aspectRatio <- length(yorder) / length(xorder)
  
  df[, xcol] <- factor(df[, xcol], levels = xorder)
  df[, ycol] <- factor(df[, ycol], levels = yorder)
  df <- df[order(df[, xcol], df[, ycol]), ]
  
  # Build plot
  p <- ggplot(df, aes(
    x = .data[[xcol]], 
    y = .data[[ycol]], 
    color = .data[[color_col]], 
    size = ifelse(.data[[size_col]] > 0, .data[[size_col]], NA)
  )) +
    geom_point() +
    xlab(xcol) +
    ylab(ycol) +
    theme_BOR(border = TRUE) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.25, 0.5, 0.25, 1), "cm"),
      aspect.ratio = aspectRatio,
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    guides(
      colour = guide_colorbar(title = color_label),
      size = guide_legend(title = size_label)
    )
  
  # Color scale
  if (!is.null(cmap)) {
    if (!is.null(colorLims)) {
      p <- p + scale_color_gradientn(colors = cmap, limits = colorLims, oob = scales::squish)
    } else {
      p <- p + scale_color_gradientn(colors = cmap)
    }
  } else {
    p <- p + scale_color_viridis_c(option = "magma")
  }
  
  # Size scale
  if (!is.null(sizeLims)) {
    p <- p + scale_size_continuous(limits = sizeLims, range = c(0.5, 5))
  } else {
    p <- p + scale_size_continuous(range = c(0.5, 5))
  }
  
  return(p)
}

################################################################
# JACCARD SIMILARITY BETWEEN DATASETS
################################################################

run_metaneighbor <- function(
    scrna_obj,
    vis_obj,
    scrna_labels = "fine_clust",
    vis_labels = "celltype_final",
    assay_scrna = NULL,
    assay_vis = NULL,
    use_variable_genes = TRUE,
    n_var_genes = 2000,
    plot = TRUE,
    n_cores = 1
){
  library(Seurat)
  library(Matrix)
  library(MetaNeighbor)
  library(pheatmap)
  library(BiocParallel)
  
  if (is.null(assay_scrna)) assay_scrna <- DefaultAssay(scrna_obj)
  if (is.null(assay_vis)) assay_vis <- DefaultAssay(vis_obj)
  
  # Step 1: Select variable genes
  shared_genes <- intersect(rownames(scrna_obj[[assay_scrna]]),
                            rownames(vis_obj[[assay_vis]]))
  
  if (use_variable_genes) {
    var_genes <- unique(c(VariableFeatures(scrna_obj),
                          VariableFeatures(vis_obj)))
    var_genes <- intersect(var_genes, shared_genes)
    var_genes <- var_genes[seq_len(min(length(var_genes), n_var_genes))]
  } else {
    var_genes <- shared_genes
  }
  
  # Step 2: Extract sparse matrices (variable genes only)
  mat_list <- list(
    scRNA = GetAssayData(scrna_obj, assay = assay_scrna, layer = "data")[var_genes, ],
    Visium = GetAssayData(vis_obj, assay = assay_vis, layer = "data")[var_genes, ]
  )
  
  # Step 3: Prepare metadata
  dataset_labels <- unlist(lapply(names(mat_list), function(n) rep(n, ncol(mat_list[[n]]))))
  celltype_labels <- c(as.vector(scrna_obj[[scrna_labels]]),
                       as.vector(vis_obj[[vis_labels]]))
  
  metadata <- data.frame(dataset = dataset_labels, celltype = celltype_labels)
  rownames(metadata) <- unlist(lapply(mat_list, colnames))
  
  # Step 4: Create SummarizedExperiment (no combined cbind2)
  # We just store the list in assays, MetaNeighbor works on one assay at a time
  se <- SummarizedExperiment(
    assays = SimpleList(data = mat_list$scRNA)  # main assay, can be any dataset
  )
  
  # Step 5: Run MetaNeighborUS on variable genes with parallel
  result <- MetaNeighborUS(
    var_genes = var_genes,
    dat = se,
    study_id = metadata$dataset,
    cell_type = metadata$celltype,
    BPPARAM = MulticoreParam(workers = n_cores)
  )
  
  scrna_types <- unique(metadata$celltype[metadata$dataset == "scRNA"])
  vis_types   <- unique(metadata$celltype[metadata$dataset == "Visium"])
  
  similarity_matrix <- result[as.character(scrna_types), as.character(vis_types)]
  
  # Step 6: Plot
  if (plot) {
    pal <- colorRampPalette(c("#272E6A","white","#D51F26"))(100)
    pheatmap(similarity_matrix,
             color = pal,
             main = "MetaNeighbor similarity (AUROC)",
             border_color = NA)
  }
  
  return(similarity_matrix)
}

# ==============================================================================
# Leiden Clustering with Silhouette Scoring for Seurat v5
# ==============================================================================
#
# Dependencies:
#   install.packages("cluster")    # for silhouette()
#
# Usage:
#   source("leiden_silhouette.R")
#
#   # Run after standard Seurat preprocessing (NormalizeData -> SCTransform -> RunPCA)
#   result <- leiden_silhouette(obj, ndims = 30, output_dir = "results")
#
#   # Access the optimal resolution and object
#   result$optimal_resolution
#   result$obj  # Seurat object with all resolution columns + "leiden_optimal"
#   result$scores  # data.frame of resolution vs silhouette width
#
#   # Or run with custom resolution range
#   result <- leiden_silhouette(obj, ndims = 30,
#                               resolutions = seq(0.1, 2.0, by = 0.1),
#                               output_dir = "results")
# ==============================================================================


#' Sweep Leiden resolutions and score each with silhouette width
#'
#' @param obj            A Seurat v5 object with PCA computed.
#' @param ndims          Number of PCA dimensions for neighbor graph and
#'                       silhouette distance computation (default: 30).
#' @param resolutions    Numeric vector of resolutions to sweep
#'                       (default: seq(0.1, 2.0, by = 0.1)).
#' @param algorithm      Clustering algorithm for FindClusters. 4 = Leiden
#'                       (default). 1 = Louvain. 2 = Louvain with multilevel
#'                       refinement. 3 = SLM.
#' @param reduction      Dimensional reduction to use (default: "pca").
#' @param max_cells_sil  Maximum cells for silhouette computation. At >15k
#'                       cells the distance matrix becomes expensive, so a
#'                       random subset is used for silhouette scoring. Cluster
#'                       labels are still computed on all cells. NULL = use all
#'                       cells (default: 10000).
#' @param output_dir     Directory for saving plots. NULL = plot to device.
#' @param verbose        Print progress (default: TRUE).
#' @param seed           Random seed (default: 42).
#'
#' @return A list with:
#'   - obj:                 Seurat object with leiden columns for each resolution
#'                          and a "leiden_optimal" column for the best resolution
#'   - scores:              data.frame with columns: resolution, n_clusters,
#'                          avg_silhouette, min_cluster_size
#'   - optimal_resolution:  the resolution with highest average silhouette
#'   - silhouette_details:  per-cell silhouette values at optimal resolution
#'
leiden_silhouette <- function(obj,
                              ndims          = 30,
                              resolutions    = seq(0.1, 2.0, by = 0.1),
                              algorithm      = 4,
                              reduction      = "pca",
                              max_cells_sil  = 10000,
                              output_dir     = NULL,
                              verbose        = TRUE,
                              seed           = 42) {
  
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required. Install with: install.packages('cluster')")
  }
  library(cluster)
  
  set.seed(seed)
  
  # --- Ensure neighbors are computed ----------------------------------------
  if (!"nn" %in% names(obj@neighbors) &&
      !any(grepl("snn", names(obj@graphs), ignore.case = TRUE))) {
    if (verbose) message("Running FindNeighbors...")
    obj <- FindNeighbors(obj, dims = 1:ndims, reduction = reduction, verbose = FALSE)
  }
  
  # --- Get PCA embeddings for silhouette distances --------------------------
  pca_embeddings <- Embeddings(obj, reduction = reduction)
  available_dims <- ncol(pca_embeddings)
  if (ndims > available_dims) {
    warning(paste0("Requested ndims=", ndims, " but only ", available_dims,
                   " available. Using ", available_dims, "."))
    ndims <- available_dims
  }
  pca_embeddings <- pca_embeddings[, 1:ndims]
  n_cells <- nrow(pca_embeddings)
  
  # --- Subsample for silhouette if needed -----------------------------------
  if (!is.null(max_cells_sil) && n_cells > max_cells_sil) {
    sil_idx <- sort(sample(n_cells, max_cells_sil))
    sil_cells <- rownames(pca_embeddings)[sil_idx]
    pca_sil <- pca_embeddings[sil_idx, ]
    if (verbose) message(paste0("Using ", max_cells_sil, " / ", n_cells,
                                " cells for silhouette scoring."))
  } else {
    sil_cells <- rownames(pca_embeddings)
    pca_sil <- pca_embeddings
  }
  
  # Precompute distance matrix for silhouette (on subset)
  if (verbose) message("Computing distance matrix for silhouette scoring...")
  dist_sil <- dist(pca_sil, method = "euclidean")
  
  # --- Sweep resolutions ----------------------------------------------------
  scores <- data.frame(
    resolution       = numeric(),
    n_clusters       = integer(),
    avg_silhouette   = numeric(),
    min_cluster_size = integer(),
    stringsAsFactors = FALSE
  )
  
  if (verbose) message(paste0("\nSweeping ", length(resolutions), " resolutions ",
                              "(algorithm = ", algorithm, ")...\n"))
  
  for (res in resolutions) {
    # Cluster
    obj <- FindClusters(obj, resolution = res, algorithm = algorithm,
                        verbose = FALSE)
    
    # Get labels for silhouette subset
    col_name <- paste0("leiden_res.", res)
    
    # Seurat stores the latest clustering in seurat_clusters and in a
    # resolution-specific column. Find the right column.
    # The SCT-based column name format varies, so we grab active idents.
    labels_all <- as.integer(Idents(obj))
    names(labels_all) <- colnames(obj)
    
    # Store with a clean column name
    obj[[col_name]] <- Idents(obj)
    
    labels_sub <- labels_all[sil_cells]
    n_clust <- length(unique(labels_sub))
    
    if (n_clust < 2) {
      if (verbose) message(paste0("  res=", sprintf("%.2f", res),
                                  " -> ", n_clust, " cluster (skipped, need >= 2)"))
      next
    }
    
    # Silhouette
    sil <- silhouette(labels_sub, dist_sil)
    avg_sil <- mean(sil[, "sil_width"])
    min_size <- min(table(labels_sub))
    
    scores <- rbind(scores, data.frame(
      resolution       = res,
      n_clusters       = n_clust,
      avg_silhouette   = avg_sil,
      min_cluster_size = min_size
    ))
    
    if (verbose) {
      message(paste0("  res=", sprintf("%.2f", res),
                     "  k=", sprintf("%2d", n_clust),
                     "  avg_sil=", sprintf("%.4f", avg_sil),
                     "  min_size=", min_size))
    }
  }
  
  # --- Find optimal resolution ----------------------------------------------
  best_idx <- which.max(scores$avg_silhouette)
  best_res <- scores$resolution[best_idx]
  best_k   <- scores$n_clusters[best_idx]
  best_sil <- scores$avg_silhouette[best_idx]
  
  if (verbose) {
    message(paste0("\nOptimal resolution: ", best_res,
                   " (k=", best_k,
                   ", avg silhouette=", round(best_sil, 4), ")"))
  }
  
  # Set optimal as a named column
  best_col <- paste0("leiden_res.", best_res)
  obj[["leiden_optimal"]] <- obj[[best_col]]
  Idents(obj) <- "leiden_optimal"
  
  # Compute full silhouette at optimal for diagnostics
  labels_opt <- as.integer(obj[[best_col, drop = TRUE]])[match(sil_cells, colnames(obj))]
  sil_opt <- silhouette(labels_opt, dist_sil)
  
  # --- Plot -----------------------------------------------------------------
  .plot_silhouette_curve(scores, best_res, output_dir, verbose)
  .plot_silhouette_profile(sil_opt, best_res, best_k, output_dir, verbose)
  
  return(list(
    obj                = obj,
    scores             = scores,
    optimal_resolution = best_res,
    silhouette_details = sil_opt
  ))
}


# ==============================================================================
# Internal: silhouette vs resolution curve
# ==============================================================================

.plot_silhouette_curve <- function(scores, best_res, output_dir, verbose) {
  
  .do_plot <- function() {
    par(family = "sans", mar = c(5, 5, 4, 5))
    
    # Primary axis: silhouette
    plot(scores$resolution, scores$avg_silhouette,
         type = "b", pch = 19, col = "#2c7bb6", lwd = 2, cex = 1.2,
         xlab = "Resolution", ylab = "Average Silhouette Width",
         main = "Leiden Resolution Selection",
         cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.4)
    
    # Mark optimum
    best_idx <- which(scores$resolution == best_res)
    points(best_res, scores$avg_silhouette[best_idx],
           pch = 18, col = "#d7191c", cex = 2.5)
    text(best_res, scores$avg_silhouette[best_idx],
         labels = paste0("res=", best_res, "\nk=", scores$n_clusters[best_idx]),
         pos = 3, col = "#d7191c", cex = 0.9, font = 2)
    
    # Secondary axis: number of clusters
    par(new = TRUE)
    plot(scores$resolution, scores$n_clusters,
         type = "b", pch = 17, col = "#fdae61", lwd = 1.5, cex = 0.9,
         axes = FALSE, xlab = "", ylab = "")
    axis(4, col.axis = "#b35900", cex.axis = 1.1)
    mtext("Number of Clusters", side = 4, line = 3, col = "#b35900", cex = 1.1)
    
    legend("topright",
           legend = c("Avg Silhouette", "Num Clusters", "Optimal"),
           col = c("#2c7bb6", "#fdae61", "#d7191c"),
           pch = c(19, 17, 18),
           lwd = c(2, 1.5, NA),
           pt.cex = c(1.2, 0.9, 2),
           cex = 0.9, bg = "white")
  }
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    filepath <- file.path(output_dir, "silhouette_vs_resolution.png")
    png(filepath, width = 8, height = 6, units = "in", res = 300)
    .do_plot()
    dev.off()
    if (verbose) message(paste0("Resolution curve saved to: ", filepath))
  } else {
    .do_plot()
  }
}


# ==============================================================================
# Internal: silhouette profile at optimal resolution
# ==============================================================================

.plot_silhouette_profile <- function(sil_obj, best_res, best_k, output_dir, verbose) {
  
  .do_plot <- function() {
    # Use cluster's built-in silhouette plot
    plot(sil_obj, col = rainbow(best_k), border = NA,
         main = paste0("Silhouette Profile (res=", best_res, ", k=", best_k, ")"),
         cex.names = 0.8)
  }
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    filepath <- file.path(output_dir, "silhouette_profile_optimal.png")
    # Height scales with cluster count
    plot_height <- max(6, best_k * 0.8)
    png(filepath, width = 8, height = plot_height, units = "in", res = 300)
    .do_plot()
    dev.off()
    if (verbose) message(paste0("Silhouette profile saved to: ", filepath))
  } else {
    .do_plot()
  }
}

# ==============================================================================
# Cell-Cell Distance Heatmap for Seurat v5
# ==============================================================================
#
# Replicates the Python cell-cell distance heatmap from the 10X Skin notebook:
# cells ordered by cluster on both axes, colored sidebars for cluster identity,
# viridis (reversed) color scale for Euclidean distance in PCA space.
#
# Dependencies:
#   BiocManager::install("ComplexHeatmap")
#   install.packages("circlize")
#   install.packages("viridisLite")
#
# Usage:
#   source("plot_cell_distance.R")
#
#   # With default active idents
#   plot_cell_distance(obj, ndims = 30, output_dir = "results")
#
#   # With a specific metadata column
#   plot_cell_distance(obj, cluster_col = "leiden_res.0.4", ndims = 30,
#                      distance_cap = 25, output_dir = "results")
#
#   # Custom cluster colors (named vector)
#   my_colors <- c("0" = "#66c2a5", "1" = "#fc8d62", "2" = "#8da0cb")
#   plot_cell_distance(obj, cluster_col = "seurat_clusters",
#                      cluster_colors = my_colors, output_dir = "results")
# ==============================================================================


#' Plot cell-cell Euclidean distance heatmap ordered by cluster
#'
#' @param obj             A Seurat v5 object with PCA computed.
#' @param cluster_col     Metadata column containing cluster labels. If NULL
#'                        (default), uses active Idents.
#' @param ndims           Number of PCA dimensions for distance computation
#'                        (default: 30).
#' @param reduction       Dimensional reduction to use (default: "pca").
#' @param distance_cap    Upper limit for the distance color scale. Distances
#'                        above this are clamped. NULL = data maximum. Setting
#'                        this to ~25 matches the reference notebook and makes
#'                        block structure more visible.
#' @param cluster_colors  Named vector of colors for clusters. Names must match
#'                        cluster labels. NULL = auto-generate.
#' @param output_dir      Directory for saving the heatmap. NULL = plot to
#'                        current device.
#' @param filename        Output filename (default: "cell_distance_heatmap.png").
#' @param plot_width      Width in inches (default: 10).
#' @param plot_height     Height in inches (default: 9).
#' @param raster_quality  Raster quality for ComplexHeatmap (default: 3).
#'                        Increase for publication, decrease for speed.
#' @param verbose         Print progress (default: TRUE).
#'
#' @return Invisibly returns the distance matrix (ordered by cluster), for
#'         downstream use if needed.
#'
plot_cell_distance <- function(obj,
                               cluster_col    = NULL,
                               ndims          = 30,
                               reduction      = "pca",
                               distance_cap   = NULL,
                               cluster_colors = NULL,
                               output_dir     = NULL,
                               filename       = "cell_distance_heatmap.png",
                               plot_width     = 10,
                               plot_height    = 9,
                               raster_quality = 3,
                               verbose        = TRUE) {
  
  for (pkg in c("ComplexHeatmap", "circlize", "viridisLite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required but not installed."))
    }
  }
  
  library(ComplexHeatmap)
  library(circlize)
  library(viridisLite)
  
  # --- Get cluster labels ---------------------------------------------------
  if (is.null(cluster_col)) {
    cluster_labels <- Idents(obj)
    label_name <- "Active Idents"
  } else {
    if (!cluster_col %in% colnames(obj@meta.data)) {
      stop(paste0("Column '", cluster_col, "' not found in metadata. ",
                  "Available: ", paste(colnames(obj@meta.data), collapse = ", ")))
    }
    cluster_labels <- obj[[cluster_col, drop = TRUE]]
    label_name <- cluster_col
  }
  
  cluster_labels <- droplevels(as.factor(cluster_labels))
  cluster_levels <- levels(cluster_labels)
  n_clusters <- length(cluster_levels)
  n_cells <- length(cluster_labels)
  
  if (verbose) message(paste0("Cluster source: ", label_name))
  if (verbose) message(paste0("Cells: ", n_cells, ", Clusters: ", n_clusters))
  
  # --- Extract PCA and compute distances ------------------------------------
  if (!(reduction %in% Reductions(obj))) {
    stop(paste0("Reduction '", reduction, "' not found. Run RunPCA() first."))
  }
  
  pca_embeddings <- Embeddings(obj, reduction = reduction)
  available_dims <- ncol(pca_embeddings)
  
  if (ndims > available_dims) {
    warning(paste0("Requested ndims=", ndims, " but only ", available_dims,
                   " available. Using ", available_dims, "."))
    ndims <- available_dims
  }
  
  pca_embeddings <- pca_embeddings[, 1:ndims]
  
  if (verbose) message(paste0("Computing Euclidean distances (", ndims, " PCs)..."))
  dist_matrix <- as.matrix(dist(pca_embeddings, method = "euclidean"))
  
  # --- Order cells by cluster -----------------------------------------------
  cell_order <- order(cluster_labels)
  dist_ordered <- dist_matrix[cell_order, cell_order]
  labels_ordered <- cluster_labels[cell_order]
  
  # --- Distance cap ---------------------------------------------------------
  if (is.null(distance_cap)) {
    distance_cap <- max(dist_ordered)
  }
  dist_ordered[dist_ordered > distance_cap] <- distance_cap
  
  if (verbose) message(paste0("Distance scale: 0 - ", round(distance_cap, 1)))
  
  # --- Cluster colors -------------------------------------------------------
  if (is.null(cluster_colors)) {
    base_colors <- c(
      "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854",
      "#ffd92f", "#e5c494", "#b3b3b3", "#1b9e77", "#d95f02",
      "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d",
      "#666666", "#1f78b4", "#33a02c", "#fb9a99", "#cab2d6"
    )
    if (n_clusters > length(base_colors)) {
      base_colors <- rep(base_colors, length.out = n_clusters)
    }
    cluster_colors <- base_colors[seq_len(n_clusters)]
    names(cluster_colors) <- cluster_levels
  } else {
    # Validate that all levels have a color
    missing <- setdiff(cluster_levels, names(cluster_colors))
    if (length(missing) > 0) {
      stop(paste0("Missing colors for clusters: ",
                  paste(missing, collapse = ", ")))
    }
  }
  
  # --- Build heatmap --------------------------------------------------------
  viridis_cols <- viridisLite::viridis(256, direction = -1)
  col_fun <- colorRamp2(seq(0, distance_cap, length.out = 256), viridis_cols)
  
  row_anno <- rowAnnotation(
    Cluster = labels_ordered,
    col = list(Cluster = cluster_colors),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = unit(5, "mm")
  )
  
  col_anno <- HeatmapAnnotation(
    Cluster = labels_ordered,
    col = list(Cluster = cluster_colors),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    height = unit(5, "mm")
  )
  
  ht <- Heatmap(
    dist_ordered,
    name                 = "Euclidean\ndistance\n(PCA-space)",
    col                  = col_fun,
    cluster_rows         = FALSE,
    cluster_columns      = FALSE,
    show_row_names       = FALSE,
    show_column_names    = FALSE,
    show_row_dend        = FALSE,
    show_column_dend     = FALSE,
    left_annotation      = row_anno,
    bottom_annotation    = col_anno,
    row_title            = paste0("Cells (#0 - ", n_cells, ")"),
    column_title         = paste0("Cells (#0 - ", n_cells, ")"),
    row_title_gp         = gpar(fontsize = 16, fontfamily = "sans"),
    column_title_gp      = gpar(fontsize = 16, fontfamily = "sans"),
    column_title_side    = "bottom",
    use_raster           = TRUE,
    raster_quality       = raster_quality,
    heatmap_legend_param = list(
      title         = "Euclidean distance\n(PCA-space)",
      title_gp      = gpar(fontsize = 12, fontfamily = "sans"),
      labels_gp     = gpar(fontsize = 10, fontfamily = "sans"),
      legend_height = unit(6, "cm"),
      at            = round(seq(0, distance_cap, length.out = 6), 1),
      labels        = round(seq(0, distance_cap, length.out = 6), 1)
    )
  )
  
  cluster_legend <- Legend(
    labels    = cluster_levels,
    title     = "Cluster",
    legend_gp = gpar(fill = cluster_colors[cluster_levels]),
    title_gp  = gpar(fontsize = 12, fontfamily = "sans"),
    labels_gp = gpar(fontsize = 10, fontfamily = "sans")
  )
  
  # --- Render ---------------------------------------------------------------
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    filepath <- file.path(output_dir, filename)
    
    png(filepath,
        width  = plot_width,
        height = plot_height,
        units  = "in",
        res    = 300)
    
    draw(ht,
         annotation_legend_list = list(cluster_legend),
         merge_legend           = TRUE,
         heatmap_legend_side    = "right")
    dev.off()
    
    if (verbose) message(paste0("Heatmap saved to: ", filepath))
  } else {
    draw(ht,
         annotation_legend_list = list(cluster_legend),
         merge_legend           = TRUE,
         heatmap_legend_side    = "right")
  }
  
  invisible(dist_ordered)
}

################################################################################
# PLOT BROAD CLUSTER AND LABELED UMAP FOR PUBLISHING
################################################################################

plot_umap_hierarchical <- function(
    seurat_obj,
    fine_cluster_col = "fine_clust",
    broad_cluster_col = "mapping_cell_type",
    reduction = "umap",
    broad_labels = NULL,
    colors = NULL,
    point_size = 0.3,
    point_alpha = 0.8,
    legend_title_size = 10,
    legend_text_size = 7,
    legend_width = 0.22
) {
  
  if (is.null(colors)) {
    colors <- c(
      "FOL" = "#208A42",
      "EPI" = "#D51F26",
      "FIB" = "#272E6A",
      "ENDO" = "#89288F",
      "IMM" = "#b35900",
      "NC" = "#808080"
    )
  }
  
  if (is.null(broad_labels)) {
    broad_labels <- c(
      "Keratinocytes" = "KRT - Keratinocytes",
      "Fibroblasts" = "FIB - Fibroblasts",
      "Immune" = "IMM - Immune cells",
      "Endothelial" = "ENDO - Endothelial cells",
      "Remaining" = "NC - Neural-crest cells"
    )
  }
  
  # Extract embeddings
  embed_df <- as.data.frame(Embeddings(seurat_obj, reduction))
  colnames(embed_df) <- c("Dim1", "Dim2")
  
  # Pull metadata as plain character to avoid factor-level interference
  embed_df$fine_clust  <- as.character(seurat_obj@meta.data[[fine_cluster_col]])
  embed_df$broad_clust <- as.character(seurat_obj@meta.data[[broad_cluster_col]])
  
  # -----------------------------------------------------------------
  # Compute all ordering BEFORE converting anything to factors
  # -----------------------------------------------------------------
  
  # Broad order: by total cell count descending
  broad_order <- embed_df %>%
    dplyr::count(broad_clust, name = "n") %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::pull(broad_clust)
  
  fine_order_by_broad <- embed_df %>%
    dplyr::count(broad_clust, fine_clust, name = "n_cells") %>%
    dplyr::arrange(broad_clust, fine_clust) %>%   # alphabetical
    dplyr::group_by(broad_clust) %>%
    dplyr::summarise(
      fine_order = list(fine_clust),
      counts     = list(n_cells),
      .groups    = "drop"
    ) %>%
    { setNames(
      purrr::map2(.$fine_order, .$counts, ~ tibble(fine_clust = .x, n_cells = .y)),
      .$broad_clust
    )
    }
  
  # Now set factors with correct levels
  embed_df$broad_clust <- factor(embed_df$broad_clust, levels = broad_order)
  
  # -----------------------------------------------------------------
  # Main UMAP plot
  # -----------------------------------------------------------------
  p <- ggplot(embed_df, aes(x = Dim1, y = Dim2, color = broad_clust)) +
    geom_point(size = point_size, alpha = point_alpha, stroke = 0) +
    scale_color_manual(values = colors) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 0)
    ) +
    coord_fixed(ratio = 1)
  
  # -----------------------------------------------------------------
  # Legend (grid-based) — uses broad_order and fine_order_by_broad
  # directly, no factor/arrange ambiguity
  # -----------------------------------------------------------------
  create_legend <- function() {
    
    # Count total lines needed for vertical spacing
    n_fine_counts <- vapply(broad_order, function(bc) nrow(fine_order_by_broad[[bc]]), integer(1))
    total_lines   <- length(broad_order) +            # one header per broad group
      sum(n_fine_counts) +              # one line per fine cluster
      (length(broad_order) - 1) * 0.5  # small gap between groups
    
    line_height  <- 0.9 / total_lines
    bar_width    <- 0.02
    left_margin  <- 0.03
    text_start   <- left_margin + bar_width + 0.02
    
    grobs     <- gList()
    y_current <- 0.97
    
    for (bc in broad_order) {
      
      bc_color  <- colors[bc]
      bc_label  <- if (!is.null(broad_labels) && bc %in% names(broad_labels)) broad_labels[bc] else bc
      bc_fine   <- fine_order_by_broad[[bc]]   # already sorted by desc(n_cells)
      n_fine    <- nrow(bc_fine)
      
      section_height <- (n_fine + 1) * line_height
      
      # Colored vertical bar for this broad group
      grobs <- gList(grobs, rectGrob(
        x      = unit(left_margin, "npc"),
        y      = unit(y_current - section_height / 2, "npc"),
        width  = unit(bar_width, "npc"),
        height = unit(section_height, "npc"),
        hjust  = 0,
        vjust  = 0.5,
        gp     = gpar(fill = bc_color, col = NA)
      ))
      
      # Broad cluster header
      grobs <- gList(grobs, textGrob(
        label = bc_label,
        x     = unit(text_start, "npc"),
        y     = unit(y_current, "npc"),
        hjust = 0,
        vjust = 0.5,
        gp    = gpar(fontsize = legend_title_size, fontface = "bold", col = bc_color)
      ))
      
      y_current <- y_current - line_height
      
      # Fine cluster entries — order is guaranteed by fine_order_by_broad
      for (i in seq_len(n_fine)) {
        fine_label <- paste0(bc_fine$fine_clust[i], " (", bc_fine$n_cells[i], ")")
        
        grobs <- gList(grobs, textGrob(
          label = fine_label,
          x     = unit(text_start + 0.02, "npc"),
          y     = unit(y_current, "npc"),
          hjust = 0,
          vjust = 0.5,
          gp    = gpar(fontsize = legend_text_size, fontface = "bold", col = "black")
        ))
        
        y_current <- y_current - line_height
      }
      
      # Inter-group gap
      y_current <- y_current - line_height * 0.3
    }
    
    gTree(children = grobs)
  }
  
  legend_grob  <- create_legend()
  
  legend_panel <- ggplot() +
    annotation_custom(legend_grob, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  final_plot <- plot_grid(
    legend_panel, p,
    rel_widths = c(legend_width, 1 - legend_width),
    nrow  = 1,
    align = "h"
  )
  
  return(final_plot)
}