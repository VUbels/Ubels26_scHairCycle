#!/usr/bin/env Rscript

#################################################################
# scATAC preparation and GRN generation
#################################################################

library(ArchR)
library(parallel)
library(mclust)
library(dplyr)
library(rhdf5)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(42)
ArchR::addArchRThreads(threads = 8) 
ArchR::addArchRGenome("hg38")
getArchRGenome()

source("./scripts/helper_functions.R")
source("./scripts/setup_py_env.R")

py_location <- "/home/uvictor/miniconda3/bin/conda"
conda_info_env <- setup_py_env(project, py_location)

macs2_path <- "/home/uvictor/miniconda3/envs/macs2_env/bin/macs2"
system2(macs2_path, "--version") 


#################################################################
# PARAMETERS
#################################################################

scatac_folder <- "/mnt/d/scatac_input/greenleaf_23"
scatac_files <- list.files(scatac_folder, full.names = TRUE) %>% 
  .[grepl("\\.tsv\\.gz$", .)]
scatac_samples <- c("C_PB1", "C_PB2", "C_PB3", "C_SD1", "C_SD2", "C_SD3")


main_folder <- "./"
output_dir <- paste0(main_folder, "scatac_files/")
proj_save_folder <- output_dir
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
  outputDirectory = output_dir,
  copyArrows = FALSE
)

# Subset to mclust-filtered cells BEFORE doublet detection
skin_archr <- subsetCells(skin_archr, cellNames = all_valid_cells)

#################################################################
# IDENTIFY DOUBLETS
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

# Filtering 2703 cells from ArchRProject
# C_PB1 : 1633 of 12781 (12.8%)
# C_SD1 : 97 of 3115 (3.1%)
# C_SD3 : 220 of 4699 (4.7%)
# C_PB3 : 574 of 7580 (7.6%)
# C_PB2 : 179 of 4233 (4.2%)
# C_SD2 : 0 of 2459 (0%)
# Final cell count: 32164
# Logging for posterity

# Save after doublet filtering, overwrite to same output_dir
skin_archr <- saveArchRProject(ArchRProj = skin_archr, outputDirectory = proj_save_folder, overwrite = TRUE, load = TRUE)
# To resume from this point:
# skin_archr <- loadArchRProject(path = proj_save_folder)

#################################################################
# DIMENSIONALITY REDUCTION (LSI)
#################################################################
skin_archr <- addIterativeLSI(
  ArchRProj = skin_archr,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(resolution = c(0.2), n.start = 10),
  varFeatures = 50000,
  dimsToUse = 1:25,
  force = TRUE
)

#################################################################
# BATCH CORRECTION (if needed across samples)
#################################################################
skin_archr <- addHarmony(
  ArchRProj = skin_archr,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

#################################################################
# CLUSTERING
#################################################################
skin_archr <- addClusters(
  input = skin_archr,
  reducedDims = "Harmony",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

#################################################################
# UMAP EMBEDDING
#################################################################
skin_archr <- addUMAP(
  ArchRProj = skin_archr,
  reducedDims = "Harmony",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

# Visualize
p1 <- plotEmbedding(skin_archr, colorBy = "cellColData", name = "Sample")
p2 <- plotEmbedding(skin_archr, colorBy = "cellColData", name = "Clusters")
plotPDF(p1, p2, name = "UMAP_Clusters", addDOC = FALSE)

p_krt15 <- plotEmbedding(
  skin_archr,
  colorBy = "GeneScoreMatrix",
  name = "KRT15",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(skin_archr)
)

# Also check other keratinocyte markers
marker_genes <- c("KRT15", "KRT14", "KRT5", "KRT1", "KRT10", "ITGA6", "TP63")
p_markers <- plotEmbedding(
  skin_archr,
  colorBy = "GeneScoreMatrix",
  name = marker_genes,
  embedding = "UMAP",
  imputeWeights = getImputeWeights(skin_archr)
)

plotPDF(p_markers, name = "Keratinocyte_Markers", addDOC = FALSE, width = 8, height = 8)

#################################################################
# CALL PEAKS (required for downstream analysis)
#################################################################
#Threats need to be set to one when using wsl2 
ArchR::addArchRThreads(threads = 1) 

# First create pseudo-bulk replicates
skin_archr <- addGroupCoverages(
  ArchRProj = skin_archr,
  groupBy = "Clusters",
  force = FALSE,
  threads = 1
)

# Call peaks with MACS2
skin_archr <- addReproduciblePeakSet(
  ArchRProj = skin_archr,
  groupBy = "Clusters",
  pathToMacs2 = macs2_path,
  force = TRUE
)

# Add peak matrix
skin_archr <- addPeakMatrix(skin_archr, force = TRUE)
saveArchRProject(ArchRProj = skin_archr, outputDirectory = proj_save_folder, overwrite = TRUE)

#################################################################
# MOTIF ANNOTATIONS + CHROMVAR DEVIATIONS
#################################################################
skin_archr <- addMotifAnnotations(
  ArchRProj = skin_archr,
  motifSet = "cisbp",
  name = "Motif",
  force = TRUE
)

skin_archr <- addBgdPeaks(skin_archr, force = TRUE)

skin_archr <- addDeviationsMatrix(
  ArchRProj = skin_archr,
  peakAnnotation = "Motif",
  force = TRUE
)

#################################################################
# MARKER FEATURES PER CLUSTER
#################################################################
markersGS <- getMarkerFeatures(
  ArchRProj = skin_archr,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)")
)

markersPeaks <- getMarkerFeatures(
  ArchRProj = skin_archr,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)")
)

markersMotif <- getMarkerFeatures(
  ArchRProj = skin_archr,
  useMatrix = "MotifMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)")
)

#################################################################
# VISUALIZE MARKER HEATMAPS
#################################################################
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  labelMarkers = marker_genes,
  transpose = TRUE
)

heatmapMotifs <- plotMarkerHeatmap(
  seMarker = markersMotif,
  cutOff = "FDR <= 0.1 & abs(MeanDiff) >= 0.1",
  log2Norm = FALSE,  # Important for deviation matrix
  transpose = TRUE
)

plotPDF(heatmapGS, name = "GeneScores_Marker_Heatmap", width = 8, height = 6, addDOC = FALSE)
plotPDF(heatmapMotifs, name = "Motifs_Marker_Heatmap", width = 8, height = 6, addDOC = FALSE)

#################################################################
# PEAK-TO-GENE LINKS (GRN foundation)
#################################################################
skin_archr <- addPeak2GeneLinks(
  ArchRProj = skin_archr,
  reducedDims = "Harmony",
  useMatrix = "GeneScoreMatrix"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = skin_archr,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#################################################################
# POSITIVE TF REGULATORS (global)
#################################################################
corGS_MM <- correlateMatrices(
  ArchRProj = skin_archr,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "Harmony"
)

sum(is.na(corGS_MM$cor))
sum(is.na(corGS_MM$padj))

# Filter with NA handling
positiveTFs <- corGS_MM[which(corGS_MM$cor > 0.3 & corGS_MM$padj < 0.01), ]

write.csv(positiveTFs, file.path(output_dir, "positive_TF_regulators_global.csv"), row.names = FALSE)

#################################################################
# EXPORT KEY DATA FOR INTEGRATION
#################################################################

# Gene score matrix
gsm <- getMatrixFromProject(skin_archr, useMatrix = "GeneScoreMatrix")
saveRDS(gsm, file.path(output_dir, "GeneScoreMatrix.rds"))

# Motif matrix
motifm <- getMatrixFromProject(skin_archr, useMatrix = "MotifMatrix")
saveRDS(motifm, file.path(output_dir, "MotifMatrix.rds"))

# Peak matrix
peakm <- getMatrixFromProject(skin_archr, useMatrix = "PeakMatrix")
saveRDS(peakm, file.path(output_dir, "PeakMatrix.rds"))

# Cell metadata
cellmeta <- getCellColData(skin_archr)
write.csv(as.data.frame(cellmeta), file.path(output_dir, "cell_metadata.csv"))

# Peak-to-gene links
saveRDS(p2g, file.path(output_dir, "peak2gene_links.rds"))

# Embeddings
umap_coords <- getEmbedding(skin_archr, embedding = "UMAP")
harmony_coords <- getReducedDims(skin_archr, reducedDims = "Harmony")
saveRDS(umap_coords, file.path(output_dir, "UMAP_coords.rds"))
saveRDS(harmony_coords, file.path(output_dir, "Harmony_coords.rds"))

# Peak set
peakset <- getPeakSet(skin_archr)
saveRDS(peakset, file.path(output_dir, "peak_set.rds"))

#################################################################
# FINAL SAVE
#################################################################
saveArchRProject(ArchRProj = skin_archr, outputDirectory = proj_save_folder, overwrite = TRUE)

message("scATAC preprocessing complete.")
message(paste("Cells:", nCells(skin_archr)))
message(paste("Peaks:", length(getPeakSet(skin_archr))))
message(paste("Output:", proj_save_folder))