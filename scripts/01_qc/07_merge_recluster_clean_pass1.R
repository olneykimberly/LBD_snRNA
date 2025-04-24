## ----setup, include=FALSE--------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

# Libraris, paths, colors
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender"
color.panel <- dittoColors()
set.seed(28)

## ----read_object-----------------------------------------------------------------------------------------------------
# read object
dataObject.astrocyte <- readRDS(paste0("../rObjects/", projectID, "_astrocyte.rds"))
dataObject.oligodendrocyte <- readRDS(paste0("../rObjects/",projectID,"_oligodendrocyte.rds"))
dataObject.opc <- readRDS(paste0("../rObjects/",projectID,"_opc.rds"))
dataObject.endothelial <- readRDS(paste0("../rObjects/",projectID,"_endothelial.rds"))
dataObject.mural<- readRDS(paste0("../rObjects/",projectID,"_mural.rds"))
dataObject.microglia <- readRDS(paste0("../rObjects/",projectID,"_microglia.rds"))
dataObject.neuron <- readRDS(paste0("../rObjects/",projectID,"_neuron.rds"))

## ----doublet_checks-----------------------------------------------------------------------------------------------------
# FeaturePlot(dataObject.neuron,
#             reduction = "umap",
#             features = c("TYROBP", "MOG", "AQP4", "RBFOX3"))
# DimPlot(dataObject.neuron,
#         group.by = "group")
# DimPlot(dataObject.neuron,
#         group.by = "Sample_ID")
## ----clean_objects-----------------------------------------------------------------------------------------------------
dataObject.astrocyte.clean <- subset(dataObject.astrocyte, cells = WhichCells(dataObject.astrocyte, idents = setdiff(unique(dataObject.astrocyte$seurat_clusters), c(2, 11, 13, 17, 18, 19, 21, 22, 26))))
dataObject.astrocyte.dirty <- subset(dataObject.astrocyte, cells = WhichCells(dataObject.astrocyte, idents = c(2, 11, 13, 17, 18, 19, 21, 22, 26)))
rm(dataObject.astrocyte)

dataObject.oligodendrocyte.clean <- subset(dataObject.oligodendrocyte, cells = WhichCells(dataObject.oligodendrocyte, idents = setdiff(unique(dataObject.oligodendrocyte$seurat_clusters), c(13, 19))))
dataObject.oligodendrocyte.dirty <- subset(dataObject.oligodendrocyte, cells = WhichCells(dataObject.oligodendrocyte, idents = c(13, 19)))
rm(dataObject.oligodendrocyte)

dataObject.opc.clean <- subset(dataObject.opc, cells = WhichCells(dataObject.opc, idents = setdiff(unique(dataObject.opc$seurat_clusters), c(8, 15, 19, 20))))
dataObject.opc.dirty <- subset(dataObject.opc, cells = WhichCells(dataObject.opc, idents = c(8, 15, 19, 20)))
rm(dataObject.opc)

dataObject.endothelial.clean <- subset(dataObject.endothelial, cells = WhichCells(dataObject.endothelial, idents = setdiff(unique(dataObject.endothelial$seurat_clusters), c(3, 5 , 6, 11))))
dataObject.endothelial.dirty <- subset(dataObject.endothelial, cells = WhichCells(dataObject.endothelial, idents = c(3, 5 , 6, 11)))
rm(dataObject.endothelial)

dataObject.mural.clean <- subset(dataObject.mural, cells = WhichCells(dataObject.mural, idents = setdiff(unique(dataObject.mural$seurat_clusters), c(8, 13, 14, 15))))
dataObject.mural.dirty <- subset(dataObject.mural, cells = WhichCells(dataObject.mural, idents = c(8, 13, 14, 15)))
rm(dataObject.mural)

dataObject.microglia.clean <- subset(dataObject.microglia, cells = WhichCells(dataObject.microglia, idents = setdiff(unique(dataObject.microglia$seurat_clusters), c(5, 7, 13, 14, 15, 19))))
dataObject.microglia.dirty <- subset(dataObject.microglia, cells = WhichCells(dataObject.microglia, idents =c(5, 7, 13, 14, 15, 19)))
rm(dataObject.microglia)

dataObject.neuron.clean <- subset(dataObject.neuron, cells = WhichCells(dataObject.neuron, idents = setdiff(unique(dataObject.neuron$seurat_clusters), c(27, 25, 24, 19, 18, 15, 14))))
dataObject.neuron.dirty <- subset(dataObject.neuron, cells = WhichCells(dataObject.neuron, idents = c(27, 25, 24, 19, 18, 15, 14)))
rm(dataObject.neuron)

## ----merge_objects_and_save-----------------------------------------------------------------------------------------------------
# clean
dataObject.clean <- merge(x = dataObject.astrocyte.clean, y = c(dataObject.oligodendrocyte.clean, dataObject.opc.clean, dataObject.endothelial.clean, dataObject.mural.clean, dataObject.microglia.clean, dataObject.neuron.clean))
saveRDS(dataObject.clean, paste0("../rObjects/",projectID,"_clean.rds"))
# dirty
dataObject.dirty <- merge(x = dataObject.astrocyte.dirty, y = c(dataObject.oligodendrocyte.dirty, dataObject.opc.dirty, dataObject.endothelial.dirty, dataObject.mural.dirty, dataObject.microglia.dirty, dataObject.neuron.dirty))
saveRDS(dataObject.dirty, paste0("../rObjects/",projectID,"_dirty.rds"))

## ----read_in_clean_data-----------------------------------------------------------------------------------------------------
# Read in clean data
dataObject.clean <- readRDS(paste0("../rObjects/",projectID,"_clean.rds"))

# Inspect clean
dataObject.clean

# Visualize the number of cell counts per sample
data <- as.data.frame(table(dataObject.clean$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")

ncells <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) +
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency),
            position=position_dodge(width=0.9),
            #vjust=-0.25,
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,30000, by = 5000), limits = c(0,30000)) +
  ggtitle("Filtered: nuclei per sample") +
  theme(legend.position =  "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
path <- paste0("../results/nuclei_count/",projectID,
               "_cells_per_sample_after_recluster_before_integration")
ncells
saveToPDF(paste0(path, ".pdf"), width = 12, height = 5)

# mean cell count
mean(data$frequency)

# median cell count
median(data$frequency)

# Inspect 
dataObject.clean
# Inspect layers
Layers(dataObject.clean[["RNA"]])
# Join layers
dataObject.clean[["RNA"]] <- JoinLayers(dataObject.clean[["RNA"]])

# Inspect after join layers
Layers(dataObject.clean[["RNA"]])
# Set indents to sample ID
Idents(dataObject.clean) <- "Sample_ID"
dataObject.clean # inspect

## ----split_by_sample-----------------------------------------------------------------------------------------------------
# Split
dataObject.clean[["RNA"]] <- split(dataObject.clean[["RNA"]], f = dataObject.clean$Sample_ID)
# Note that since the data is split into layers, normalization and variable feature identification is performed for each sample independently 
# (a consensus set of variable features is automatically identified).

# Inspect
dataObject.clean
# Layers
Layers(dataObject.clean[["RNA"]])

# Re-normalizing data via SCTransform
dataObject.clean <- SCTransform(dataObject.clean, verbose = TRUE, conserve.memory = TRUE)
all.genes <- rownames(dataObject.clean)
# Inspect 
dataObject.clean

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dataObject.clean), 10)
top10

# run PCA on the merged object
dataObject.clean <- RunPCA(object = dataObject.clean)
Idents(dataObject.clean) 

# Determine the K-nearest neighbor graph
dataObject.clean <- FindNeighbors(object = dataObject.clean, 
                                  assay = "SCT", 
                                  reduction = "pca",
                                  dims = 1:15)
# Run UMAP
dataObject.clean <- RunUMAP(dataObject.clean,
                            dims = 1:15,
                            reduction = "pca",
                            n.components = 3) 

# Determine the clusters for various resolutions
dataObject.clean <- FindClusters(object = dataObject.clean,
                                 algorithm = 1, # 1= Louvain
                                 resolution = 0.9)

# Set indents to SCT resolution 
Idents(dataObject.clean) <- dataObject.clean$SCT_snn_res.0.9
dataObject.clean$seurat_clusters <- dataObject.clean$SCT_snn_res.0.9

# Save 
saveRDS(dataObject.clean, paste0("../rObjects/",projectID,"_after_recluster_before_integration_SCTransform.rds"))
#dataObject.clean <- readRDS(paste0("../rObjects/",projectID,"_merged_reclusters_clean_SCTransform.rds"))


#----- plot UMAPS
# group 
pdf(paste0("../results/UMAP/unannotated/", projectID,
           "_after_recluster_before_integration_group.pdf"), width = 9, height = 7)
DimPlot(dataObject.clean,
        group.by = "group")
dev.off()

# Sample ID
pdf(paste0("../results/UMAP/unannotated/", projectID,
    "_after_recluster_before_integration_Sample_ID.pdf"), width = 9, height = 7)
DimPlot(dataObject.clean,
        group.by = "Sample_ID")
dev.off()

markers.to.plot <-
  c(
    "CLU", 
    "GFAP", 
    "AQP4", 
    "GJA1",
    "CLDN5",
    "ADGRF5",
    "FLT1",
    "COL1A1",
    "COL1A2",
    "DCN",
    "HEXB",
    "C1QA",
    "C1QB",
    "C1QC",
    "ITGAM",
    "TYROBP",
    "P2RY12",
    "AIF1",
    "RBFOX1",
    "RBFOX3", 
    "SNAP25",
    "SYT1",
    "GAD1",
    "GAD2",
    "PLP1",
    "MBP", 
    "MOG", 
    "OLIG1",
    "PDGFRA",
    "VCAN",
    "TNR",
    "ACTA2",
    "VTN"
  )

dot_ind_celltype <- DotPlot(dataObject.clean,
                            features = markers.to.plot, 
                            cluster.idents = FALSE,
                            dot.scale = 8) + RotatedAxis()
dot_ind_celltype
pdf(
  paste0(
    "../results/dot_plot/",
    projectID,
    "_after_recluster_before_integration.pdf"
  ),width = 10, height = 12)
dot_ind_celltype
dev.off()


## ----IntegrateLayers_via_Harmony-----------------------------------------------------------------------------------------------------
dataObject.integrated <- IntegrateLayers(
  object = dataObject.clean, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
dataObject.integrated
# Determine the K-nearest neighbor graph
dataObject.integrated <- FindNeighbors(object = dataObject.integrated, 
                                       reduction = "harmony", 
                                       dims = 1:30)

# Determine the clusters for various resolutions
dataObject.integrated <- FindClusters(object = dataObject.integrated, resolution = 0.9)
dataObject.integrated <- RunUMAP(dataObject.integrated, reduction = "harmony", dims = 1:30)
Idents(dataObject.integrated) 

# Set indents to SCT resolution 
Idents(dataObject.integrated) <- dataObject.integrated$SCT_snn_res.0.9
dataObject.integrated$seurat_clusters <- dataObject.integrated$SCT_snn_res.0.9

## ----save_integrated-----------------------------------------------------------------------------------------------------
saveRDS(dataObject.integrated, paste0("../rObjects/",projectID,"_after_recluster_harm_int.rds"))
dataObject.integrated <- readRDS(paste0("../rObjects/",projectID,"_after_recluster_harm_int.rds"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Visualize the number of cell counts per sample
data <- as.data.frame(table(dataObject.integrated$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")

ncells <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) +
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency),
            position=position_dodge(width=0.9),
            #vjust=-0.25,
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,30000, by = 5000), limits = c(0,30000)) +
  ggtitle("Filtered: nuclei per sample") +
  theme(legend.position =  "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
path <- paste0("../results/nuclei_count/",projectID,
               "_cells_per_sample_after_recluster_harm_int")
ncells
saveToPDF(paste0(path, ".pdf"), width = 12, height = 5)



#----- plot UMAPS
# group 
pdf(paste0("../results/UMAP/unannotated/", projectID,
           "_after_recluster_harm_int_group.pdf"), width = 9, height = 7)
DimPlot(dataObject.integrated,
        group.by = "group")
dev.off()

# Sample ID
pdf(paste0("../results/UMAP/unannotated/", projectID,
           "_after_recluster_harm_int_Sample_ID.pdf"), width = 9, height = 7)
DimPlot(dataObject.integrated,
        group.by = "Sample_ID")
dev.off()

## ----Azimuth_predictions-----------------------------------------------------------------------------------------------------
dataObject.integrated[["RNA"]] <- JoinLayers(dataObject.integrated[["RNA"]])

# Inspect after join layers
Layers(dataObject.integrated[["RNA"]])

dataObject.integrated <- RunAzimuth(dataObject.integrated, reference = "humancortexref")
# currently, the object has two layers in the RNA assay: counts, and data
# DimPlot(dataObject, group.by = "predicted.subclass", label = TRUE, label.size = 3) + NoLegend()
ditto_umap <- dittoDimPlot(object = dataObject.integrated,
                           var="predicted.subclass",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
path <- paste0("../results/UMAP/annotated/",projectID,
               "_harmony_integrated_annoated_Azimuth_after_recluster")
ditto_umap
saveToPDF(paste0(path, ".pdf"), width = 12, height = 8)
## ----End-----------------------------------------------------------------------------------------------------
