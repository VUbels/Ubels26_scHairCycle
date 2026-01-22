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

gene_list = c("PDGFRA", "FGF7", "VIM", "EDN3", "WNT5A", "KRT1", "KRT10", "KRT6A",
              "KRT17", "KRT75", "KRT35", "KRT85", "PECAM1", "VWF", "TAGLN", "DSP",
              "MLANA", "PMEL")

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
source("./scripts/ambient_rna_removal.R")
source("./scripts/doublet_removal.R")

py_location <- "/home/uvictor/miniconda3/bin/conda"
conda_info_env <- setup_py_env(project, py_location)

#################################################################
# RUNNING BROAD MARKER GENES FOR INITIAL CLUSTERIZATION
#################################################################

# Easy visualization through Nebulosa to get a better overview of gene expression particularly
# for when low cell count has high gene expression in a particular cluster

plot_marker_genes(obj, genes = gene_list, pt_size = 0.4)
