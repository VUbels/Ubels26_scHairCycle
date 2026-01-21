###############################################################
# DOUBLET REMOVAL USING SCDBLFINDER 
###############################################################

# We prefer scDblFinder over DoubletFinder due to significantly
# reduced false positive doublet results for highly heteroskedastic data
# as found in hair follicle scRNA data.

doublet_identification <- function(object.list, assumed_doublet_rate = 0.004) {

  for (i in seq_along(object.list)) {
    obj <- object.list[[i]]
    
    cat("\n###################################################\n")
    cat("Running scDblFinder on", unique(obj$orig.ident), "dataset ...\n")
    cat("###################################################\n")
    
    # Convert Seurat object to SingleCellExperiment
    sce <- as.SingleCellExperiment(obj)
    
    # Run scDblFinder
    # clusters=FALSE uses the random approach (recommended for complex datasets)
    # Set seed for reproducibility
    set.seed(123)
    sce <- scDblFinder(sce, clusters = FALSE, dbr.per1k	= assumed_doublet_rate)
    
    # Extract doublet information and add to Seurat object
    obj$scDblFinder.class <- sce$scDblFinder.class
    obj$scDblFinder.score <- sce$scDblFinder.score
    
    # Print doublet statistics
    doublet_table <- table(sce$scDblFinder.class)
    cat("Doublets detected:", doublet_table["doublet"], 
        "(", round(doublet_table["doublet"]/ncol(obj)*100, 2), "%)\n")
    cat("Singlets:", doublet_table["singlet"], 
        "(", round(doublet_table["singlet"]/ncol(obj)*100, 2), "%)\n")
    
    # Update object in list
    object.list[[i]] <- obj
    
    rm(sce)
  }
  
  return(object.list)
  
}
###############################################################
# VISUALIZATION OF DOUBLET DETECTION
###############################################################

doublet_visualization <- function(object.list, output_folder) {
    
    dir.create(paste0(output_folder, "/doublet_detection/"), showWarnings = FALSE)
    
    for (i in seq_along(object.list)) {
    obj <- object.list[[i]]
    stage <- unique(obj$orig.ident)
    
    # Preprocessing for visualization
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
    
    # Plot doublets on UMAP
    p1 <- DimPlot(obj, group.by = "scDblFinder.class", pt.size = 0.1) + 
      ggtitle(paste(stage, "- Doublet Classification"))
    
    p2 <- FeaturePlot(obj, features = "scDblFinder.score", pt.size = 0.1) + 
      ggtitle(paste(stage, "- Doublet Score"))
    
    # Save plots
    combined_plot <- p1 + p2
    ggsave(filename = paste0(output_folder, "/doublet_detection/", stage, "_doublet_detection.png"),
           plot = combined_plot, width = 12, height = 5)
    
    print(combined_plot)
    
    object.list[[i]] <- obj
  }
  
}
###############################################################
# FILTER DOUBLETS
###############################################################

filter_doublets <- function(object.list) {

  # Filter singlets only
  object.list_filtered <- lapply(object.list, function(obj) {
    stage <- unique(obj$orig.ident)
    obj_filtered <- subset(obj, subset = scDblFinder.class == "singlet")
    cat("Remaining cells after doublet removal for", stage, ":", ncol(obj_filtered), "cells\n")
    return(obj_filtered)
  })

}