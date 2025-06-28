## ----setup-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

## ----source------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender"
color.panel <- dittoColors()

## ----read_object--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read object
dataObject <- readRDS(file = paste0("../rObjects/", projectID, "_doublets.rds"))
dataObject # inspect
## ----reprocess----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Join layers
dataObject <- NormalizeData(dataObject)
dataObject <- FindVariableFeatures(dataObject)
dataObject <- ScaleData(dataObject)
dataObject <- RunPCA(dataObject, features = VariableFeatures(dataObject))
Layers(dataObject)
# Type
table(dataObject$group) # doublets by group

dataObject.integrated <- IntegrateLayers(
  object = dataObject, method = HarmonyIntegration,
  normalization.method = "RNA",
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
# Determine the K-nearest neighbor graph
dataObject.integrated <- FindNeighbors(object = dataObject.integrated, 
                                       reduction = "harmony", # pca, harmony 
                                       dims = 1:15)

# Determine the clusters for various resolutions
dataObject.integrated <- FindClusters(object = dataObject.integrated, resolution = 0.2)
dataObject.integrated <- RunUMAP(dataObject.integrated, reduction = "harmony", dims = 1:15)

# inspect 
dataObject.integrated@reductions
dataObject.integrated@assays
dataObject.integrated[["RNA"]] <- JoinLayers(dataObject.integrated[["RNA"]])
# Azimuth annotations
dataObject.integrated <- RunAzimuth(dataObject.integrated, reference = "humancortexref")

p1 <- DimPlot(
  dataObject.integrated,
  reduction = "harmony",
  group.by = c("Sample_ID"),
  combine = FALSE, label.size = 2
)
p1
ditto_umap <- dittoDimPlot(object = dataObject.integrated,
                           var = "predicted.subclass",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)

ditto_umap


pdf(
  paste0(
    "../results/UMAP/doublets/",
    projectID,
    "_doublets_only.pdf"
  ),
  width = 7,
  height = 5
)
ditto_umap
dev.off()

dot_ind <- DotPlot(dataObject.integrated,
                   features = markers.to.plot, 
                   cluster.idents = TRUE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0(
    "../results/dot_plot/",
    projectID,
    "_doublets_only.pdf"
  ),
  width = 14,
  height = 10
)
dot_ind
dev.off()

## ----save_object,echo=FALSE,eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(dataObject.integrated, paste0("../rObjects/",projectID,"_unannotated_harmony_int_doublets_only.rds"), compress = FALSE)
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------