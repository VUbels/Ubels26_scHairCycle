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
library(XML)
library(GenomicFeatures)
library(Banksy)
library(SeuratWrappers)
library(ggplot2)

options(future.globals.maxSize = 5 * 1024^3)
set.seed(43)

# There may be an incompatability caused by XML due to libxml2
# In this case simply force install of libxml2 to the latest version
# and run renv::install("XML", type = "source", rebuild = TRUE) to
# force reinstall into the renv.

#################################################################
# SETUP PROJECT PARAMETERS
#################################################################
project <- "ubels26_haircycle"
main_folder <- "./"

obj <- readRDS(paste0(main_folder, "post_filter_integrated_objects.RDS"))

#################################################################
# SETTING UP VISIUM DATA
#################################################################
#vis_obj <- readRDS(paste0(main_folder, "Spatial_scalp_S2_final.rds"))
#vis_obj <- UpdateSeuratObject(vis_obj)

vis_obj <- readRDS(paste0(main_folder, "annotated_vis_obj.rds"))

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
#conda_info_env <- setup_py_env(project, py_location)

#################################################################
# RUNNING BROAD MARKER GENES FOR INITIAL CLUSTERIZATION
#################################################################

# Easy visualization through Nebulosa to get a better overview of gene expression particularly
# for when low cell count has high gene expression in a particular cluster
obj <- cluster_subcluster(obj, output_dir = "./")

Idents(obj) <- "SCT_snn_res.0.8"
broad_markers <- FindAllMarkers(obj, min.pct = 0.1, logfc.threshold = 1, only.pos = TRUE)
write.csv(broad_markers, file = "./marker_genes/broad_markers.csv")

plot_marker_genes(obj = obj, 
                  genes = broad_gene_list, 
                  cluster_col = "SCT_snn_res.0.8",
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
# RUNNING BROAD MARKER GENES FOR INITIAL CLUSTERIZATION SPATIAL
#################################################################

vis_markers <- FindAllMarkers(vis_obj, group.by = "broad_cluster", only.pos = TRUE)

Idents(vis_obj) <- "SCT_snn_res.0.8"
cells <- CellsByIdentities(vis_obj)
clusters_to_show <- unique(vis_obj@active.ident)

p <- SpatialDimPlot(
  vis_obj,
  pt.size.factor = 4,
  cells.highlight = cells[clusters_to_show],
  cols.highlight = c("red","grey30"),
  facet.highlight = TRUE,
  combine = TRUE,
  image.alpha = 0
) + NoLegend()

p

#################################################################
# ASSIGN BROAD MARKER IDENTIFICATION TO CLUSTERS
#################################################################

vis_cluster_identification <- list(
  `0` = "Basal.KCs",              #KRT5/KRT14/COL7A1
  `1` = "Center.Bulge",           #KRT15/COMP/LGR5
  `2` = "ORS.Suprabasal.1",       #PDZRN3/KRT6/EFNA5/ACTN1
  `3` = "Endothelial",            #PECAM1/VWF/PLVAP
  `4` = "Dermal.Sheath",          #DCN/COL6A1/COL3A1
  `5` = "FBs.Residual",           #VIM/FN1
  `6` = "Epid.Spinous",           #KRT1/KRT10/KRTDAP
  `7` = "Sebaceous.Gland",        #KRT79/GAL/MUC1
  `8` = "Junctional.Zone",        #S100A8/S100A9/CALML5
  `9` = "ORS.Suprabasal.1",       #PDZRN3/KRT6A/KRT6C
  `10` = "Epid.Spinous",          #KRT1/KRT10/KRTDAP
  `11` = "IRS.Huxley",            #KRT71/TCHH/KRT74+
  `12` = "ORS.Suprabasal.2",      #CST6/S100A1/SHF/SPRR4
  `13` = "Eccrine.Gland",         #KRT19/MUCL1/AQP5
  `14` = "Cortex",                #KRT31/KRT83/KRT86
  `15` = "Sebaceous.Gland",       #KRT79/GAL/MUC1
  `16` = "Ishtmus",               #KLK7/KLK6
  `17` = "IRS.Henley",            #KRT71/TCHH/KRT74-
  `18` = "Epid.Tran",             #LOR/KRT1/KRTDAP
  `19` = "IRS.Huxley",            #KRT71/TCHH/KRT74+
  `20` = "Dermal.Sheath",         #TAGLN/ACTA2/VIM/COL6A2
  `21` = "Epid.Granulosum",       #KRT1/KRT10/KRTDAP/LOR
  `22` = "Matrix.Melanocyte",     #KRT35/KTR85/MLANA
  `23` = "ORS.Basal",             #KRT5//COMP/FGF18
  `24` = "SG.Sebocytes",          #KRT79/GAL/MUC1/PPARG
  `25` = "Cuticle",               #KRT32/KRT82
  `26` = "Companion.Layer",       #KRT75/TCHH/WNT11
  `27` = "Upper.Matrix",          #KRT35/KRT85
  `28` = "ORS.Suprabasal.2",      #CST6/S100A1/SHF/SPRR4
  `29` = "Cortex",                #KRT31/KRT83/KRT86
  `30` = "Dermal.Papilla",        #VCAN/PDGFRA/DKK2/RSPO2/RSPO4
  `31` = "Eccrine.Sebocytes",     #SCGB1D2/PPARG/AQP5
  `32` = "Germinal.Matrix",       #TOP2A/CDK1/DLX3
  `33` = "FBs.Perifollicular",    #COL3A1/FN1/COL6A3
  `34` = "Upper.Bulge"            #COL17A1/KRT15/DKK3
)

vis_obj$broad_cluster <- unname(unlist(vis_cluster_identification[as.character(vis_obj$SCT_snn_res.0.8)]))
Idents(vis_obj) <- "broad_cluster"

p <- SpatialDimPlot(vis_obj, label = TRUE, pt.size.factor = 4, image.alpha = 0, repel = TRUE, label.size = 2)
ggsave(filename = "./marker_genes/annotated_spatial_follicle.png", p, width = 15, height = 10)

saveRDS(vis_obj, file = "./annotated_vis_obj.rds")

#################################################################
# ASSIGN BROAD MARKER IDENTIFICATION TO SCRNA CLUSTERS
#################################################################

broad_cluster_identification <- list(
  `0` = "Germinal.Centre",           #GPX2/MOXD1/BMERB1
  `1` = "ORS.Suprabasal.1",          #PDZRN3/KRT6/EFNA5/ACTN1
  `2` = "Upper.Bulge.1",             #DKK3/DIO2
  `3` = "Epid.Spinous.1",            #KRT1/KRT10/KRTDAP
  `4` = "Endothelial.1",             #PECAM1/VWF/PLVAP
  `5` = "Lower.Bulge.1",             #THBS2/LGR5/COMP
  `6` = "Center.Bulge.2",
  `7` = "Dermal.Papilla",            #VCAN/PDGFRA/TWIST2/
  `8` = "",
  `9` = "Unspecified",               #Further subclustering necessary
  `10` = "Upper.Bulge.2",
  `11` = "Center.Bulge.1",
  `12` = "Junctional.Zone",          #TMEM45A/KLK6/SPINK5
  `13` = "Epid.Tran",                #LAMB4/SERPINB2/CDHR1
  `14` = "Epid.Spinous.2",
  `15` = "Dermal.Sheath",            #TAGLN/ACTA2/VIM/COL6A2
  `16` = "",
  `17` = "Neural.Progenitor",        #CNTNAP2/CTNNA2/FMN2
  `18` = "",
  `19` = "Lower.Bulge.2",            #LGR/COMP/EDN2
  `20` = "Matrix",                   #KRT35/KRT85
  `21` = "",
  `22` = "ORS.Suprabasal.2",         #SERPINA3/FABP5/CD24
  `23` = "T.cells",                  #CD3D/TRBC1
  `24` = "Endothelial.2",            #PECAM1/VWF/PLVAP
  `25` = "Melanocytes",              #CDH19/MLANA/SOX10
  `26` = "Macrophages",              #HLA-DRA/CD74
  `27` = "Pericytes"                 #RGS5
)

sub_markers <- FindMarkers(obj, ident.1 = "2", ident.2 = "10")
DimPlot(obj, label = TRUE, split.by = "dataset")

SpatialDimPlot(vis_obj, label = TRUE, repel = TRUE, image.alpha = 0, label.size = 3, pt.size.factor = 5)

gene = "DIO2"
g1 <- DimPlot(obj, label = TRUE)
g2 <- FeaturePlot(obj, features = gene, label = TRUE, order = TRUE)
g3 <- ImageFeaturePlot(vis_obj, features = gene, size = 1)
plot(g1+g2+g3)


obj$broad_cluster <- unname(unlist(broad_cluster_identification[as.character(obj$SCT_snn_res.0.8)]))
visualize_percentage_clusters(seurat_obj = obj, clusters = "broad_cluster", phases = "orig.ident", output_dir = paste0(main_folder, "marker_genes"))
top_genes <- top_expressed_per_cluster(obj, cluster_col = "seurat_clusters", n_genes = 50)


#################################################################
# INITIALIZE FINE CLUSTER COLUMN
#################################################################

DimPlot(obj, label = TRUE, group.by = "SCT_snn_res.0.8")
FeaturePlot(obj, feature = "PIP", order = TRUE)

SpatialDimPlot(vis_obj, group.by = "celltype_final", image.alpha = 0, pt.size.factor = 4, label = TRUE, label.size = 3)
SpatialFeaturePlot(vis_obj, feature = "SERPINA3", image.alpha = 0, pt.size.factor = 4, shape = 21)

# Start with broad labels; subclustering below will overwrite
# cells belonging to each cohort with finer annotations.
obj$fine_clust <- obj$broad_cluster


##########################################################
# FINE GRAINED CLUSTERING FOR NON-STATIAL MATRIX
##########################################################
res <- subcluster_and_markers(
  obj,
  cluster_name = "Matrix.Cortex",
  cluster_col = "fine_clust",
  resolution = 0.2
)

sub_obj <- res$sub_obj
sub_markers <- res$sub_markers

sub_map <- c(
  "0" = "Matrix.Cortex", #KRT35/KRT85
  "1" = "Cuticle",       #KRT25/KRT71/KRT74-
  "2" = "Matrix.Cortex", #KRT35/KRT85
  "3" = "Eccrine.Gland"  #KRT19/PIP/MUCL1
)

obj <- insert_subclusters(
  obj,
  sub_obj,
  sub_map
)

##########################################################
# FINE GRAINED CLUSTERING FOR MATRIX
##########################################################

similarity_matrix <- run_metaneighbor(
  scrna_obj = obj,
  vis_obj = vis_obj,
  scrna_labels = "fine_clust",
  vis_labels = "celltype_final"
)

DimPlot(obj, label = TRUE)


fine_markers <- FindAllMarkers(obj, group.by = "fine_clust")
saveRDS(obj, file = paste0(main_folder, "annotated.RDS"))
