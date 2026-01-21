#!/usr/bin/env Rscript

#################################################################
# Description: Preproccesing data using CellBender and scDblFinder
#################################################################

# Following script contains the full preprocessing pipeline for the Ubels26_HairCycle
# publication. Renv/Conda environments are dynamically set. Pytorch compatability
# has to be set by user and will not be supported.

#################################################################
# LIBRARY LOADING
#################################################################

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)

#################################################################
# SETUP PROJECT PARAMETERS
#################################################################

project <- "ubels26_haircycle"
input_folder <- "/mnt/d/scrna_datasets/ubels26_scrna_dataset"
main_folder <- "./"
scrna_data <- list.files(input_folder, recursive = FALSE, include.dirs = FALSE, pattern = ".h5")

# Variables are pure an identification list for file output,
# list has to identical in length to datasets
dataset_names <- c("Anagen", "Catagen", "Telogen")

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
# RUN CELLBENDER - AMBIENT RNA REMOVAL
#################################################################

remove_ambient_rna(input_folder = input_folder, cellbender_learning_rate = 0.00005)
make_seurat_compatible(output_folder = main_folder) 

#################################################################
# DATA LOADING AND INITIAL QC
#################################################################

show_qc_metrics(input_folder = main_folder)

#################################################################
# INITIAL QC FILTERING AND SETTING UNIQUE BARCODES PER SAMPLE
#################################################################

object.list <- filter_by_qc(
  input_folder = main_folder,
  project_names = dataset_names,
  min_feature = 150, 
  max_feature = 7000, 
  min_count = 200,
  max_count = 50000, 
  max_percent_mt = 35
)

#################################################################
# SCDBLFINDER DOUBLET DETECTION
#################################################################

# Per 10X guidelines on V4 chips we assume 0.004 per 1000 cell doublet rate. Change to 0.008 for 10k V3 chips.

object.list <- doublet_identification(object.list = object.list, assumed_doublet_rate = 0.004)
doublet_visualization(object.list = object.list, output_folder = paste0(main_folder, "preprocessing"))
object.list <- filter_doublets(object.list)

#################################################################
# INTEGRATE AND SAVE FILTERED OBJECT
#################################################################

integrated_obj <- scrna_integrate(
  object.list = object.list,
  output_folder = main_folder,
  dataset_names = dataset_names,
  python_script_path = "./scripts/integrate_scanorama.py",
  python_path = conda_info_env[["python_path"]]
)

saveRDS(integrated_obj, file = paste0(main_folder, "post_filter_integrated_objects.RDS"))
