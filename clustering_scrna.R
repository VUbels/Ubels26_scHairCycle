#!/usr/bin/env Rscript

#################################################################
# Description: Clustering and evaluation of subclusters using CDI
#################################################################

# Following script contains the full clusteringpipeline for the Ubels26_HairCycle
# publication. Renv/Conda environments are dynamically set. Pytorch compatability
# has to be set by user and will not be supported.

#################################################################
# LIBRARY LOADING
#################################################################

library(Seurat)
library(SeuratDisk)

library(Nebulosa)
library(CDI)
library(reticulate)

#################################################################
# SETUP PROJECT PARAMETERS
#################################################################

project <- "ubels26_haircycle"
main_folder <- "./"
obj <- readRDS(paste0(main_folder, "post_filter_integrated_objects.RDS"))
vis_obj <- readRDS(paste0(main_folder, "Spatial_scalp_S2_final.rds"))

#################################################################
# SETUP PY ENVIRONMENT
#################################################################

# Please note that for GPU support you need to manually change
# parameters in setup_py_env.R Due to this being highly user 
# dependent, questions regarding setting up appropriate pytorch
# compatibility will not be supported. CellBender can run 
# without GPU support but this will take a very long time.

source("./scripts/helper_functions.R")
source("./scripts/setup_py_env.R")
source("./scripts/gene_lists.R")

py_location <- "/home/uvictor/miniconda3/bin/conda"
conda_info_env <- setup_py_env(project, py_location)

#################################################################
# RUNNING BROAD MARKER GENES FOR INITIAL CLUSTERIZATION
#################################################################

# Easy visualization through Nebulosa to get a better overview of gene expression particularly
# for when low cell count has high gene expression in a particular cluster

broad_markers <- FindAllMarkers(obj, min.pct = 0.1, logfc.threshold = 0.3)

plot_marker_genes(obj = obj, 
                              genes = broad_gene_list, 
                              cluster_col = "seurat_clusters",
                              reduction = "umap", 
                              output_dir = "./marker_genes/broad_markers", 
                              pt_size = 1,
                              outline_size = 0.25,
                              concavity = 5,
                              show_labels = TRUE,
                              eps = 2,
                              min_pts = 25,
                              outlier_percentile = 0.98)

#################################################################
# ASSIGN BROAD MARKER IDENTIFICATION TO CLUSTERS
#################################################################

broad_cluster_identification <- list(
 `0` = "Lower.Follicle", #PDZRN3/SERPINA3/SOX9/KRT75
 `1` = "Lower.Follicle", #PDZRN3/SOX9
 `2` = "Central.Follicle", #WNT5A/KRT19
 `3` = "Central.Follicle", #WNT5A/
 `4` = "Matrix", #KRT35/KRT85/WNT10B
 `5` = "Upper.Follicle", #KRT1/KRT10/KRT79/AQP3
 `6` = "Upper.Follicle", #LGR6/SOX7/WNT10A
 `7` = "Upper.Follicle", #KRT1/KRT10/KRT79/WNT3
 `8` = "Endothelial", #PLVAP/PECAM1/VWF/SOX18/VIM
 `9` = "Fibroblasts", #RGS5 but also EDN3/VIM/CD34
 `10` = "Immune", #CD3D/CD53/CD69
 `11` = "Fibroblasts", #FGF7/TWIST2/PDGFRA
 `12` = "Melanocytes", #PMEL/MLANA/SOX10
 `13` = "Fibroblasts", #TALGN/VIM/ACTA2
 `14` = "Central.Follicle", #EGF/FGF14
 `15` = "Lower.Follicle", #Likely differentiating cells? #KR17
 `16` = "Neural.Progenitors" #CTNNA2/CDH9
)
obj$broad_cluster <- unname(unlist(broad_cluster_identification[as.character(obj$seurat_clusters)]))
visualize_percentage_clusters(seurat_obj = obj, clusters = "broad_cluster", phases = "orig.ident", output_dir = paste0(main_folder, "marker_genes"))

#################################################################
# SPECIFYING FINE CLUSTER IDENTIFICATION FOR LARGER COHORTS
#################################################################

subset_obj <- subset(obj, subset = obj$broad_cluster == c("Lower.Follicle", "Central.Follicle"))
subset_obj <- cluster_subcluster(subset_obj, output_dir = "./marker_genes/")
Idents(subset_obj) <- "SCT_snn_res.0.5"

subset_cluster_markers <- FindAllMarkers(subset_obj, min.pct = .1, only.pos = TRUE)
plot_marker_genes(obj = subset_obj, 
                  genes = unique(unlist(gene_list$keratinocytes)), 
                  cluster_col = "SCT_snn_res.0.6",
                  reduction = "umap", 
                  output_dir = "./marker_genes/lower_follicle", 
                  pt_size = 1,
                  outline_size = 0.25,
                  concavity = 5,
                  show_labels = TRUE,
                  eps = 2,
                  min_pts = 25,
                  outlier_percentile = 0.98)


subset_obj <- subset(obj, subset = obj$broad_cluster == c("Upper.Follicle"))
subset_obj <- cluster_subcluster(subset_obj, output_dir = "./marker_genes/")
Idents(subset_obj) <- "SCT_snn_res.0.6"

plot_marker_genes(obj = subset_obj, 
                  genes = unique(unlist(gene_list$keratinocytes)), 
                  cluster_col = "SCT_snn_res.0.6",
                  reduction = "umap", 
                  output_dir = "./marker_genes/upper_follicle", 
                  pt_size = 1,
                  outline_size = 0.25,
                  concavity = 5,
                  show_labels = TRUE,
                  eps = 2,
                  min_pts = 25,
                  outlier_percentile = 0.98)
