#!/usr/bin/env Rscript

#################################################################
# scATAC preparation and GRN generation
#################################################################

library(ArchR)
library(parallel)
library(mclust)
library(dplyr)
library(rhdf5)

set.seed(42)
ArchR::addArchRThreads(threads = 6) 
ArchR::addArchRGenome("hg38")

source("./scripts/helper_functions.R")

#################################################################
# SETUP
#################################################################

scatac_folder <- "/mnt/d/scatac_input/greenleaf_23"
scatac_files <- list.files(scatac_folder, full.names = TRUE) %>% 
  .[grepl("\\.tsv\\.gz$", .)]
scatac_samples <- c("C_PB1", "C_PB2", "C_PB3", "C_SD1", "C_SD2", "C_SD3")


main_folder <- "./"
output_dir <- paste0(main_folder, "scatac_files/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#################################################################
# CREATE ARROW FILES (loose filtering)
#################################################################
ArrowFiles <- createArrowFiles(
  inputFiles = scatac_files,
  sampleNames = scatac_samples,
  outputNames = file.path(output_dir, scatac_samples),
  QCDir = file.path(output_dir, "QualityControl"),
  minTSS = 0,
  minFrags = 500,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)

#################################################################
# MCLUST-BASED CELL FILTERING
#################################################################

# Apply with plotting
valid_barcodes <- lapply(
  ArrowFiles, 
  mclust_filter_cells, 
  min_tss = 5, 
  min_frags = 1000,
  plot_dir = output_dir
)

all_valid_cells <- unlist(valid_barcodes)

#################################################################
# CREATE ARCHR PROJECT AND SUBSET FIRST
#################################################################
skin_archr <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = file.path(output_dir, "skin_scatac"),
  copyArrows = FALSE
)

# Subset to mclust-filtered cells BEFORE doublet detection
skin_archr <- subsetCells(skin_archr, cellNames = all_valid_cells)

#################################################################
# IDENTIFY DOUBLETS (on filtered cells only)
#################################################################
skin_archr <- addDoubletScores(
  input = skin_archr,
  k = 10,
  UMAPParams = list(n_neighbors = 50, min_dist = 0.4, metric = "cosine", verbose = FALSE),
  LSIParams = list(dimsToUse = 1:50, varFeatures = 50000),
  knnMethod = "UMAP",
  LSIMethod = 1
)

# Filter doublets
skin_archr <- filterDoublets(skin_archr, filterRatio = 1)
message(paste0("Final cell count: ", nCells(skin_archr)))