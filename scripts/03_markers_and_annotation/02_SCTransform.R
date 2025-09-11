## ----setup-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

#https://satijalab.org/seurat/articles/integration_introduction
#https://satijalab.org/seurat/articles/seurat5_integration 

## ----source------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
project_ID <- "CWOW_cellbender_singlets" # CWOW_cellbender_joinlayers_azimuth
color.panel <- dittoColors()

## ----read_object--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read object
dataObject <- readRDS(file = paste0("../rObjects/", project_ID, ".rds"))
dataObject # inspect

# Notes
## un-normalized counts (layer='counts')
## normalized data (layer='data')
## z-scored/variance-stabilized data (layer='scale.data')

# Reorder samples
barcodes <- colnames(dataObject)
sample <- str_match(barcodes, "(.+)_[ACGT]+")[,2]
dataObject$sample <- factor(sample, levels = order_samples)
table(dataObject$sample)  # check
Idents(dataObject) <- dataObject$sample
rm(sample, barcodes)

# Add metadata
rownames(metadata) <- metadata$sampleID
# Reorder the metadata data frame to match the column order of the Seurat object
barcodes <- colnames(dataObject)
sample <- str_match(barcodes, "(.+)_[ACGT]+")[,2]
metadata_reordered <- metadata[sample, ]
# Add the entire metadata dataframe to the Seurat object
# Add the reordered metadata to the Seurat object
dataObject <- AddMetaData(object = dataObject, metadata = metadata_reordered)
dataObject$group <- factor(dataObject$TYPE, levels = c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS"))

## ----Preprocess_for_integration----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#dataObject[["RNA"]] <- split(dataObject[["RNA"]], f = dataObject$Sample_ID)

# since the data is split into layers, normalization and variable feature identification is performed for each sample independently (a consensus set of variable features is automatically identified).
dataObject # inspect
# Layers
Layers(dataObject[["RNA"]])

dataObject <- SCTransform(dataObject, verbose = TRUE, conserve.memory = TRUE) # this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
DefaultAssay(dataObject) # inspect
saveRDS(dataObject, paste0("../rObjects/",project_ID,"_SCTransform_only.rds"), compress = FALSE)

print("Identify the 10 most highly variable genes:")
top10 <- head(VariableFeatures(dataObject), 10)
top10

dataObject <- RunPCA(dataObject, npcs = 30, verbose = F)
dataObject <- RunUMAP(dataObject, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
dataObject <- FindNeighbors(dataObject, dims = 1:30, reduction = "pca", assay = "SCT")
dataObject <- FindClusters(dataObject, resolution = 0.6, cluster.name = "unintegrated_clusters")
dataObject <- RunUMAP(dataObject, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
project_ID <- "CWOW_cellbender_SCTransform"
saveRDS(dataObject, paste0("../rObjects/",project_ID,".rds"), compress = FALSE)
