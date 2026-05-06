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
library(Nebulosa)
library(CDI)
library(reticulate)
library(XML)
library(GenomicFeatures)
library(Banksy)
library(SeuratWrappers)
library(ggplot2)
library(scCustomize)

options(future.globals.maxSize = 5 * 1024^3)
set.seed(123)

# There may be an incompatability caused by XML due to libxml2
# In this case simply force install of libxml2 to the latest version
# and run renv::install("XML", type = "source", rebuild = TRUE) to
# force reinstall into the renv.

#################################################################
# SETUP PROJECT PARAMETERS
#################################################################
project <- "ubels26_haircycle"
main_folder <- "./"
obj <- readRDS(paste0(main_folder, "post_filter_nonintegrated_objects.RDS"))

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
obj <- cluster_subcluster(obj, output_dir = "./", ndims = 30)
Idents(obj) <- "SCT_snn_res.0.6"

obj <- PrepSCTFindMarkers(obj)
broad_markers <- FindAllMarkers(obj, min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)


plot_marker_genes(obj = obj, 
                  genes = broad_gene_list, 
                  cluster_col = "SCT_snn_res.0.6",
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

vis_markers <- FindAllMarkers(vis_obj, group.by = "SCT_snn_res.0.8", only.pos = TRUE)

Idents(vis_obj) <- "SCT_snn_res.0.6"
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

ImageFeaturePlot(vis_obj, features = c("IVL", "KRT75", "KRT71", "GATA3"))


#################################################################
# ASSIGN MAPPING CELL TYPE MARKER IDENTIFICATION TO SCRNA CLUSTERS
#################################################################

mapping_cluster_identification <- list(
  `1` = "FOL.1",             #PDZRN3/KRT6B
  `2` = "FOL.2",             #DIO2/KRT15/DKK3
  `3` = "EPI.1",             #MOXD1/TRPC6/PTN
  `4` = "FOL.3",             #FGF14/GPC5/DIO2
  `5` = "FOL.4",             #LGR5/COMP/KRT15
  `6` = "EPI.2",             #KRT1/KRT10/KRTDAP
  `7` = "ENDO.1",            #PECAM1/VWF/TALAM1
  `8` = "FOL.5",             #SERPINA3/KRT16
  `9` = "EPI.3",             #SPINK5/S100A8/S100A9
  `10` = "FIB.1",            #VCAN/LAMA2/TWIST2
  `11` = "FIB.2",            #PRKG1/TAGLN/ACTA2
  `12` = "EPI.4",            #LAMB4/SOX5/KRT5
  `13` = "EPI.5",            #CYP3A5/KRT14
  `14` = "FOL.6",            #KRT85/KRT35/KRT32
  `15` = "NC.1",             #MLANA/PAX3/PMEL
  `16` = "IMM.1",            #PTPRC/CD3D/RUNX3
  `17` = "IMM.2"             #HLA-DRA/CD74/CD86
)

obj$mapping_cell_type <- unname(unlist(mapping_cluster_identification[as.character(obj$SCT_snn_res.0.6)]))
Idents(obj) <- "mapping_cell_type"

DimPlot(obj, label = TRUE)

# SHOW DISTRIBUTION OF CLUSTER SAMPLES PRE FINE GRAINED ANNOTATION
visualize_percentage_clusters(seurat_obj = obj, clusters = "mapping_cell_type", phases = "orig.ident", output_dir = paste0(main_folder, "marker_genes"))
top_genes <- top_expressed_per_cluster(obj, cluster_col = "SCT_snn_res.0.6", n_genes = 50)
write.csv(broad_markers, file = "./marker_genes/mapping_cell_type_markers.csv")
saveRDS(obj, file = "./mapping_annotated_obj.rds")

#################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO FOLLICLE CLUSTERS
#################################################################

obj$broad_cluster <- as.character(obj$mapping_cell_type)
fol_obj <- subset(obj, idents = grep("^FOL", levels(Idents(obj)), value = TRUE))
fol_obj <- cluster_subcluster(fol_obj, output_dir = "./marker_genes/fol/", ndims = 20, n_genes = 3000, v_features = 2000, ncells = 5000)

Idents(fol_obj) <- "SCT_snn_res.0.9"
fol_markers <- FindAllMarkers(fol_obj, min.pct = 0.1, logfc.threshold = 1, only.pos = TRUE)

fol_cluster_identification <- list(
  `1` = "uHF.Bulge.1",           #KRT15/CDH12/LGR6
  `2` = "aHF.ORS.1",             #CDH1/NOTCH2/SPINK5
  `3` = "aHF.Bulge.1",           #KRT15/COL11A1/WNT5A
  `4` = "aHF.ORS.2",             #SERPINA3/KRT6B/KRT16
  `5` = "uHF.Bulge.2",           #KRT15/COL17A1/DIO2
  `6` = "uHF.Bulge.3",           #DKK3/KRT15/COL17A1
  `7` = "aHF.SBulge.1",          #CD34/COMP/NGFR
  `8` = "uHF.Bulge.4",           #KRT15/NFATC2
  `9` = "aHF.ORS.3",             #SNX31/INTS6/TRIM9
  `10` = "uHF.Bulge.4",          #KRT15/CUX2
  `11` = "aHF.ORS.4",            #NOTCH3/CMYA5/VTCN1
  `12` = "aHF.Matrix",           #KRT35/KRT85
  `13` = "aHF.SBulge.2",         #LGR5/COMP/DPP10
  `14` = "aHF.ORS.5",            #EGR3/OVOL1/CST6
  `15` = "Ecc.Gland"             #KRT19/MUCL1/SOX10
)

new_ids <- unlist(fol_cluster_identification)
fol_obj$broad_cluster <- plyr::mapvalues(
  x = as.character(Idents(fol_obj)),  
  from = names(new_ids),
  to = new_ids
)
Idents(fol_obj) <- "broad_cluster"
fol_obj$fine_clust <- as.character(fol_obj$broad_cluster)

################################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO LOWER FOLLICLE KERATINOCYTES
################################################################################

unique(obj$broad_cluster)
DimPlot(obj, group.by = "broad_cluster", label = TRUE)
Idents(obj) <- "broad_cluster"

ors_obj <- subset(fol_obj, idents = grep("^aHF.ORS", levels(Idents(fol_obj)), value = TRUE))
ors_obj <- cluster_subcluster(ors_obj, output_dir = "./marker_genes/fol/", ndims = 20, kparam = 15, n_genes = 2000, v_features = 2000, ncells = 3000)

Idents(ors_obj) <- "SCT_snn_res.0.4"
ors_markers <- FindAllMarkers(ors_obj, min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE, min.diff.pct = 0.3)

top_markers <- ors_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 100) %>%
  ungroup()

write.csv(top_genes, file = "./marker_genes/fol/cell_type_markers.csv")
top_genes <- top_expressed_per_cluster(ors_obj, cluster_col = "SCT_snn_res.0.4", n_genes = 50)

ors_cluster_identification <- list(
  `1` = "aHF.ORS.B",               #THSD4/CD44
  `2` = "aHF.ORS.SB.1",            #S100A1/KRT6C
  `3` = "aHF.ORS.SB.2",            #SERPINA3/KRT8 
  `4` = "aHF.LPC",                 #MEG3/EZH2/ITGA6
  `5` = "aHF.ORS.Reg",             #VTCN1/MYO9A
  `6` = "aHF.UM",                  #OVOL1/SCARB1
  `7` = "aHF.CL"                   #SPINK5/WNT11/KRT75
)

new_ids <- unlist(ors_cluster_identification)
ors_obj$fine_clust <- plyr::mapvalues(
  x = as.character(Idents(ors_obj)),  
  from = names(new_ids),
  to = new_ids
)
non_ors_cells <- Cells(ors_obj)[ors_obj$fine_clust != "aHF.ORS"]
fol_obj$fine_clust[non_ors_cells] <- ors_obj$fine_clust[non_ors_cells]

################################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO MATRIX
################################################################################
Idents(obj) <- "fine_clust"
mx_obj <- subset(obj, idents = "aHF.Matrix")
mx_obj <- cluster_subcluster(mx_obj, output_dir = "./marker_genes/mx/", ndims = 5, kparam = 5, n_genes = 500, v_features = 500, ncells = 1000)
Idents(mx_obj) <- "SCT_snn_res.0.1"
mx_markers <- FindAllMarkers(mx_obj, min.pct = 0.1, logfc.threshold = 1, only.pos = TRUE)

mx_cluster_identification <- list(
  `1` = "aHF.Matrix",           #KRT35/KRT85
  `2` = "aHF.IRS",              #GATA3/KRT71
  `3` = "uHF.Bulge.3",          #DKK3/KRT15/COL17A1
  `4` = "aHF.Bulge.1"           #LGR5/COMP/COL4A1
)

new_ids <- unlist(mx_cluster_identification)
mx_obj$fine_clust <- plyr::mapvalues(
  x = as.character(Idents(mx_obj)),  
  from = names(new_ids),
  to = new_ids
)

obj$fine_clust[Cells(mx_obj)] <- mx_obj$fine_clust
DimPlot(obj, label = TRUE, group.by = "fine_clust")

################################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO EPIDERMAL LAYERS
################################################################################

epi_obj <- subset(obj, idents = grep("^EPI", levels(Idents(obj)), value = TRUE))
epi_obj <- cluster_subcluster(epi_obj, output_dir = "./marker_genes/epi/", ndims = 20, n_genes = 3000, v_features = 2000, ncells = 5000)

Idents(epi_obj) <- "SCT_snn_res.0.5"
epi_markers <- FindAllMarkers(epi_obj, min.pct = 0.1, logfc.threshold = 1, only.pos = TRUE)

top_markers <- epi_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 100) %>%
  ungroup()

write.csv(top_markers, file = "./marker_genes/epi/cell_type_markers.csv")

epi_cluster_identification <- list(
  `1` = "pHF.Isth",       # MOXD1/TNC/CDH13        — Isthmus keratinocytes; MOXD1 isthmus marker, no IFE basal markers (COL17A1/KRT14 absent), ORS-adjacent signature
  `2` = "IFE.Spi",        # KRT1/KRT10/DSG1        — Spinous layer; canonical suprabasal keratins + desmosomal maturation
  `3` = "pHF.Inf",        # KRT6B/SLC15A1/SPINK5   — Infundibulum; KRT6A/B trichilemmal keratins, SLC15A1 highly specific, SOX9+
  `4` = "IFE.SB",         # RCOR1/NDRG1/PLA2G4E    — Suprabasal transitional; lost basal markers, not yet spinous, chromatin remodeling active
  `5` = "IFE.TA",         # FTL/FTH1/IFITM3        — Transit-amplifying; 10/30 top markers ribosomal, ferritin+, immune-primed, KRT14 retained
  `6` = "IFE.BSC",        # COL17A1/LAMB4/IGFBP3   — Basal stem cells; COL17A1-high, full integrin complement, WNT3/WNT9A active
  `7` = "IFE.Cyc",        # MIR17HG/NABP1/RCC1     — Cycling S/G2; MIR17HG drives G1-S, NABP1 S-phase repair, RCC1 mitotic condensation
  `8` = "IFE.Gra"         # S100A7/SERPINB4/CSTA   — Granular layer; terminal differentiation antimicrobials + cornified envelope assembly
)

new_ids <- unlist(epi_cluster_identification)
epi_obj$broad_cluster <- plyr::mapvalues(
  x = as.character(Idents(epi_obj)),  
  from = names(new_ids),
  to = new_ids
)
Idents(epi_obj) <- "broad_cluster"
epi_obj$fine_clust <- as.character(epi_obj$broad_cluster)
obj$fine_clust[Cells(epi_obj)] <- epi_obj$fine_clust

################################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO IMMUNE CLUSTERS
################################################################################

im_obj <- subset(obj, idents = grep("^IMM", levels(Idents(obj)), value = TRUE))
im_obj <- cluster_subcluster(im_obj, output_dir = "./marker_genes/immune/", ndims = 10, n_genes = 1000, v_features = 1000, ncells = 1000, kparam = 10)

Idents(im_obj) <- "SCT_snn_res.0.3"
imm_markers <- FindAllMarkers(im_obj, min.pct = 0.1, logfc.threshold = 1, only.pos = TRUE)

top_markers <- imm_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 100) %>%
  ungroup()

write.csv(top_markers, file = "./marker_genes/immune/cell_type_markers.csv")

imm_cluster_identification <- list(
  `1` = "IMM.CD4.TRM",    # SKAP1/BCL11B/CD69      — CD4+ tissue-resident T helper cells
  `2` = "IMM.Mac",        # CD68/LYZ/C1QA          — Dermal macrophages (C1Q+ tissue-resident)
  `3` = "IMM.LC",         # CD207/CD1A/CDH1         — Langerhans cells
  `4` = "IMM.Th17",       # KLRB1/CCR6/IL26         — Th17 skin-homing T cells
  `5` = "IMM.Treg"        # FOXP3/CTLA4/IL2RA       — Regulatory T cells (CCR8+ tissue-resident)
)

new_ids <- unlist(im_cluster_identification)
im_obj$broad_cluster <- plyr::mapvalues(
  x = as.character(Idents(im_obj)),  
  from = names(new_ids),
  to = new_ids
)
Idents(im_obj) <- "broad_cluster"
im_obj$fine_clust <- as.character(im_obj$broad_cluster)
obj$fine_clust[Cells(im_obj)] <- im_obj$fine_clust

################################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO ENDOTHELIAL CLUSTERS
################################################################################

end_obj <- subset(obj, idents = grep("^ENDO", levels(Idents(obj)), value = TRUE))
end_obj <- cluster_subcluster(end_obj, output_dir = "./marker_genes/endo/", ndims = 20, n_genes = 1000, v_features = 1000, ncells = 2000, kparam = 10)

Idents(end_obj) <- "SCT_snn_res.0.2"
end_markers <- FindAllMarkers(end_obj, min.pct = 0.1, logfc.threshold = 1, only.pos = TRUE)

top_markers <- end_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 100) %>%
  ungroup()

write.csv(top_markers, file = "./marker_genes/endo/cell_type_markers.csv")

endo_cluster_identification <- list(
  `1` = "ENDO.Art",        # SOX17/SEMA3G/GJA4        — Arterial endothelium
  `2` = "ENDO.PCV",        # SELE/ACKR1/VCAM1         — Postcapillary venule
  `3` = "ENDO.Cap",        # EMCN/EGFL7/NOSTRIN        — Capillary endothelium
  `4` = "VSMCs"            # PDGFRB/RGS5/NOTCH3        — vSMC (non-endothelial)
)

new_ids <- unlist(endo_cluster_identification)
end_obj$broad_cluster <- plyr::mapvalues(
  x = as.character(Idents(end_obj)),  
  from = names(new_ids),
  to = new_ids
)
Idents(end_obj) <- "broad_cluster"
end_obj$fine_clust <- as.character(end_obj$broad_cluster)
obj$fine_clust[Cells(end_obj)] <- end_obj$fine_clust

################################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO FIBROBLAST CLUSTERS
################################################################################

fib_obj <- subset(obj, idents = grep("^FIB", levels(Idents(obj)), value = TRUE))
fib_obj <- cluster_subcluster(fib_obj, output_dir = "./marker_genes/fibroblast/", ndims = 10, n_genes = 1000, v_features = 1000, ncells = 2000, kparam = 10)

Idents(fib_obj) <- "SCT_snn_res.0.4"
fib_markers <- FindAllMarkers(fib_obj, min.pct = 0.1, logfc.threshold = 1, only.pos = TRUE)

top_markers <- fib_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 100) %>%
  ungroup()

write.csv(top_markers, file = "./marker_genes/fibroblast/cell_type_markers.csv")

fib_cluster_identification <- list(
  `1` = "FIB.Pericytes",    # RGS5/ABCC9/PDGFRB        — Pericytes
  `2` = "DP",               # C3/PTGDS/CFH             — Telogen dermal papilla (complement+, quiescent)
  `3` = "DP",               # ADAM12/MMP2/MMP14        — Catagen dermal papilla (ECM remodeling)
  `4` = "FIB.DS",           # MYH11/CNN1/ACTA2         — Dermal sheath (phasic contractile)
  `5` = "DP",               # CORIN/DKK2/FOXP2         — Early anagen dermal papilla (activation/transition)
  `6` = "FIB.DS",           # SOX6/RGS6/NFASC          — Early dermal sheath
  `7` = "DP"                # SOX2/HAPLN1/FGF7         — Anagen dermal papilla (peak inductive)
)

new_ids <- unlist(fib_cluster_identification)
fib_obj$fine_clust <- plyr::mapvalues(
  x = as.character(Idents(fib_obj)),  
  from = names(new_ids),
  to = new_ids
)

fib_obj$fine_clust[Cells(fib_obj)] <- fib_obj$fine_clust

################################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO DP CLUSTERS
################################################################################

Idents(fib_obj) <- "fine_clust"
dp_obj <- subset(fib_obj, idents = "DP")
dp_obj <- cluster_subcluster(dp_obj, output_dir = "./marker_genes/dp/", ndims = 5, kparam = 5, n_genes = 500, v_features = 500, ncells = 1000)
Idents(dp_obj) <- "SCT_snn_res.0.3"

dp_markers <- FindAllMarkers(dp_obj, min.pct = 0.1, logfc.threshold = 1, only.pos = TRUE)

top_markers <- dp_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 100) %>%
  ungroup()

write.csv(top_markers, file = "./marker_genes/dp/cell_type_markers.csv")

dp_cluster_identification <- list(
  `1` = "DP.Remod",       # ADAM12/MMP14/BMP6        — ECM degradation and notch signaling
  `2` = "DP.Qui",         # C3/CFH/CHRDL1            — Dermal papilla (complement+, quiescent)
  `3` = "DP.Str",         # CORIN/THSD4/VCAN         — Dermal papilla (No specific markers)
  `4` = "DP.Ind",         # PTGDS/APCDD1/RSPO4       — Anagen inductive dermal papilla
  `5` = "DP.Cup"          # ACTA2/VCAN/TAGLN         — DP/DS combined markers
)

new_ids <- unlist(dp_cluster_identification)
dp_obj$fine_clust <- plyr::mapvalues(
  x = as.character(Idents(dp_obj)),  
  from = names(new_ids),
  to = new_ids
)

fib_obj$fine_clust[Cells(dp_obj)] <- dp_obj$fine_clust
obj$fine_clust[Cells(fib_obj)] <- fib_obj$fine_clust

################################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO NEURAL CREST CLUSTERS
################################################################################

nc_obj <- subset(obj, idents = grep("^NC", levels(Idents(obj)), value = TRUE))
nc_obj <- cluster_subcluster(nc_obj, output_dir = "./marker_genes/neural_crest/", ndims = 10, n_genes = 1000, v_features = 1000, ncells = 1000, kparam = 10)

Idents(nc_obj) <- "SCT_snn_res.0.5"
nc_markers <- FindAllMarkers(nc_obj, min.pct = 0.1, logfc.threshold = 1, only.pos = TRUE)

top_markers <- nc_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 100) %>%
  ungroup()

write.csv(top_markers, file = "./marker_genes/neural_crest/cell_type_markers.csv")

nc_cluster_identification <- list(
  `1` = "MEL.1",        # RBMS3 (only sig. DEG)     
  `2` = "MEL.1",        # EDNRB/MITF/SOX10          
  `3` = "MEL.NCP",         # S100B/VCAN/ZEB1           — Neural crest progenitors
  `4` = "MEL.Diff"         # CRABP1/OCA2/SLC24A5       — Differentiated melanocytes
)

new_ids <- unlist(nc_cluster_identification)
nc_obj$broad_cluster <- plyr::mapvalues(
  x = as.character(Idents(nc_obj)),  
  from = names(new_ids),
  to = new_ids
)
Idents(nc_obj) <- "broad_cluster"
nc_obj$fine_clust <- as.character(nc_obj$broad_cluster)
obj$fine_clust[Cells(nc_obj)] <- nc_obj$fine_clust

################################################################################
# PLOTTING AND SAVING ANNOTATED OBJECT
################################################################################

obj$plotting_cell_type <- sub("\\..*", "", obj$mapping_cell_type)
unique(obj$plotting_cell_type)

broad_labels <- c(
  "FOL" = "FOL - Anagen Follicle",
  "EPI" = "EPI - Permanent Follicle Epidermis",
  "FIB" = "FIB - Fibroblasts",
  "IMM" = "IMM - Immune cells",
  "ENDO" = "ENDO - Endothelial cells",
  "NC" = "NC - Neural-crest cells"
)

p <- plot_umap_hierarchical(
  obj,
  fine_cluster_col = "fine_clust",
  broad_cluster_col = "plotting_cell_type",
  broad_labels = broad_labels,
  point_size = 0.5,
  legend_width = 0.22
)

ggsave("umap_hierarchical.svg", p, width = 9, height = 7.5)


saveRDS(obj, file = "./annotated_obj.rds")
#obj <- readRDS(file = "./annotated_obj.rds")
#################################################################
# SETUP PY ENVIRONMENT
#################################################################

# Please note that for GPU support you need to manually change
# parameters in setup_py_env.R Due to this being highly user 
# dependent, questions regarding setting up appropriate pytorch
# compatibility will not be supported. CellBender can run 
# without GPU support but this will take a very long time.

py_location <- "/home/uvictor/miniconda3/bin/conda"
conda_info_env <- setup_py_env(project, py_location)

reticulate::py_install(
  packages = c("anndata", "scipy", "h5py"),
  pip = TRUE
)

reticulate::py_module_available("anndata")

as.anndata(x = obj, file_path = "./", file_name = "obj_anndata.h5ad")
