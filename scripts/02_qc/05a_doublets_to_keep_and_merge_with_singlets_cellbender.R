## ----setup-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

## ----source------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender"
color.panel <- dittoColors()

## ----read_object--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read object
dataObject.doublets <- readRDS(file = paste0("../rObjects/", projectID, "_unannotated_harmony_int_doublets_only.rds"))
dataObject.singlets <- readRDS(file = paste0("../rObjects/", projectID, "_singlets.rds"))

## ----doublets_to_keep--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ditto_umap <- dittoDimPlot(object = dataObject.doublets,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)

ditto_umap
dataObject.doublets.keep <- subset(dataObject.doublets, cells = WhichCells(dataObject.doublets, idents = setdiff(unique(dataObject.doublets$seurat_clusters), c(11,12, 3, 19, 14))))

## ----merge_objects_and_save-----------------------------------------------------------------------------------------------------
# clean
dataObject <- merge(x = dataObject.singlets, y = c(dataObject.doublets.keep))
saveRDS(dataObject, paste0("../rObjects/",projectID,"_mereged_singlets_with_kept_doublets.rds"), compress = FALSE)
dataObject # inspect