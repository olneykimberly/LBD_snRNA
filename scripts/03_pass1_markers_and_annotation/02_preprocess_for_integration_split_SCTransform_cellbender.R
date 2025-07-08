## ----setup-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

#https://satijalab.org/seurat/articles/integration_introduction
#https://satijalab.org/seurat/articles/seurat5_integration 

## ----source------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender_joinlayers_azimuth"
color.panel <- dittoColors()

## ----read_object--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read object
dataObject <- readRDS(file = paste0("../rObjects/", projectID, ".rds"))
dataObject # inspect

all_meta_data_columns <- colnames(dataObject@meta.data)
columns_to_remove <- all_meta_data_columns[startsWith(all_meta_data_columns, c("pANN", "DF.classifications"))]

if (length(columns_to_remove) > 0) {
  dataObject[[columns_to_remove]] <- NULL
  message(paste("Removed the following metadata columns:", paste(columns_to_remove, collapse = ", ")))
} else {
  message("No metadata columns starting with 'pANN' found to remove.")
}
all_meta_data_columns <- colnames(dataObject@meta.data)
all_meta_data_columns

# clean up
dataObject$RNA_snn_res.0.2 <- NULL
dataObject$pANN_0.25_0.01_103 <- NULL
dataObject$DF.classifications_0.25_0.01_103 <- NULL

## ----Preprocess_for_integration----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dataObject[["RNA"]] <- split(dataObject[["RNA"]], f = dataObject$Sample_ID)

# since the data is split into layers, normalization and variable feature identification is performed for each sample independently (a consensus set of variable features is automatically identified).
dataObject # inspect
# Layers
Layers(dataObject[["RNA"]])

dataObject <- SCTransform(dataObject, verbose = TRUE, conserve.memory = TRUE) # this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
DefaultAssay(dataObject) # inspect
saveRDS(dataObject, paste0("../rObjects/",projectID,"_SCTransform_only.rds"), compress = FALSE)

print("Identify the 10 most highly variable genes:")
top10 <- head(VariableFeatures(dataObject), 10)
top10

dataObject <- RunPCA(dataObject, npcs = 30, verbose = F)
dataObject <- RunUMAP(dataObject, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
dataObject <- FindNeighbors(dataObject, dims = 1:30, reduction = "pca", assay = "SCT")
dataObject <- FindClusters(dataObject, resolution = 0.6, cluster.name = "unintegrated_clusters")
dataObject <- RunUMAP(dataObject, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
projectID <- "CWOW_cellbender_SCTransform"
saveRDS(dataObject, paste0("../rObjects/",projectID,".rds"), compress = FALSE)
