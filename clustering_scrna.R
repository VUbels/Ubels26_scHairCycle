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
library(CDI)

#################################################################
# SETUP PROJECT PARAMETERS
#################################################################

project <- "ubels26_haircycle"
main_folder <- "./"

obj <- readRDS(paste0(main_folder, "post_filter_integrated_objects.RDS"))

#################################################################
# SETUP PY ENVIRONMENT
#################################################################

# Please note that for GPU support you need to manually change
# parameters in setup_py_env.R Due to this being highly user 
# dependent, questions regarding setting up appropriate pytorch
# compatibility will not be supported. CellBender can run 
# without GPU support but this will take a very long time.

source("./helper_functions.R")
source("./setup_py_env.R")
source("./ambient_rna_removal.R")
source("./doublet_removal.R")

py_location <- "/home/uvictor/miniconda3/bin/conda"
conda_info_env <- setup_py_env(project, py_location)

#