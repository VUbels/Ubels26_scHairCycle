#!/usr/bin/env Rscript

#################################################################
# Description: Clustering and evaluation of subclusters
#################################################################

# Following script contains the full clusteringpipeline for the Ubels26_HairCycle
# publication. Renv/Conda environments are dynamically set. Pytorch compatability
# has to be set by user and will not be supported.

#################################################################
# LIBRARY LOADING
#################################################################

library(Seurat)
library(SeuratDisk)
library(reticulate)
library(XML)
library(GenomicFeatures)
library(SeuratWrappers)
library(ggplot2)
library(scCustomize)

options(future.globals.maxSize = 5 * 1024^3)
set.seed(123)

source("./scripts/helper_functions.R")

#################################################################
# SETUP PROJECT PARAMETERS
#################################################################
project <- "ubels26_haircycle"
main_folder <- "./"
obj <- readRDS(paste0(main_folder, "annotated_obj.rds"))
obj$mapping_cell_type <- sub("\\..*", "", obj$mapping_cell_type)

###################################################
# EXAMPLE USAGE
###################################################

# One-vs-rest: Anagen vs non-Anagen (Catagen + Telogen)
res <- test_cluster_enrichment(
  obj,
  treatment_col   = "orig.ident",
  cluster_col     = "fine_clust",
  sample_col      = "orig.ident",
  focal_condition = "Anagen"
)

# Faceted plot — needs obj to look up mapping_cell_type per cluster
plots <- plot_enrichment_by_group(
  res,
  obj         = obj,
  group_col   = "mapping_cell_type",
  cluster_col = "fine_clust",
  output_dir  = "./enrichment_plots/"
)
