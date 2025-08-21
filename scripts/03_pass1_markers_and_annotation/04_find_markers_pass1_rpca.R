## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender_RPCAIntegration"
color.panel <- dittoColors()

# read object
dataObject  <- readRDS(paste0("../rObjects/",projectID,".rds"))

dataObject$seurat_clusters <- dataObject$rpca_clusters
Idents(dataObject) <- "seurat_clusters"
# inspect
dataObject

# Given a merged object with multiple SCT models, this function uses minimum of the median UMI (calculated using the raw UMI counts) of individual objects to reverse the individual SCT regression model using minimum of median UMI as the sequencing depth covariate. 
# The counts slot of the SCT assay is replaced with recorrected counts and the data slot is replaced with log1p of recorrected counts.
dataObject.prepSCT <- PrepSCTFindMarkers(dataObject, assay = "SCT", verbose = TRUE)

## Markers per cluster
markers <- SeuratWrappers::RunPrestoAll(
  object = dataObject.prepSCT,
  assay = "SCT",
  slot = "counts",
  only.pos = FALSE
)
write.table(markers, 
            paste0("../results/markers/", projectID, "_markers.tsv"),
            quote = FALSE,
            row.names = FALSE)
saveRDS(markers, paste0("../rObjects/", projectID, "_markers.rds"))
markers <- readRDS(paste0("../rObjects/", projectID,"_markers.rds"))

# rearrange to order by cluster & filter to only include log2FC > 1 & FDR < 0.05
all.markers.strict <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05)

saveRDS(all.markers.strict, paste0("../rObjects/", projectID,"_markers_log2FC1_q0.01.rds"))
#all.markers.strict <- readRDS(paste0("../rObjects/", projectID,"_markers_log2FC1_q0.01.rds"))