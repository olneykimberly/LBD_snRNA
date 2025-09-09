## ----setup-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

#https://satijalab.org/seurat/articles/integration_introduction
#https://satijalab.org/seurat/articles/seurat5_integration 

## ----source------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender_mereged_singlets_with_kept_doublets"
color.panel <- dittoColors()

## ----read_object--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read object
dataObject <- readRDS(file = paste0("../rObjects/", projectID, ".rds"))
dataObject # inspect
# Notes
## un-normalized counts (layer='counts')
## normalized data (layer='data')
## z-scored/variance-stabilized data (layer='scale.data')

## ----JoinLayers_AZIMUTH--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# update projectID
projectID <- "CWOW_cellbender_joinlayers_azimuth"
## We will obtain predicted cell annotations via AZIMUTH for comparing integration methods 
dataObject[["RNA"]] <- JoinLayers(dataObject[["RNA"]]) # first join layers
# Azimuth annotations
dataObject <- RunAzimuth(dataObject, reference = "humancortexref")
saveRDS(dataObject, paste0("../rObjects/",projectID,".rds"), compress = FALSE)
dataObject # inspect
Idents(dataObject) <- dataObject$predicted.subclass # set idents

# annotation
ditto_umap <- dittoDimPlot(object = dataObject, var = "predicted.subclass", reduction.use = "ref.umap", do.label = TRUE, labels.highlight = TRUE)
path <- paste0("../results/integration_exploration/",projectID,"_UMAP_annotations")
ditto_umap
saveToPDF(paste0(path, ".pdf"), width = 8, height = 7)

# doublet status
UMAP_doublets <- DimPlot(dataObject, reduction = "ref.umap", group.by = c("CellTypes_DF"), cols=c("black", "#66C2A5"))
path <- paste0("../results/integration_exploration/",projectID,"_UMAP_doubletStatus")
UMAP_doublets
saveToPDF(paste0(path, ".pdf"), width = 8, height = 7)

# feature plot 
UMAP_feature <- FeaturePlot(dataObject,  reduction = "ref.umap", features = c("AQP4", "PLP1", "RBFOX3", "GAD1"))
path <- paste0("../results/integration_exploration/",projectID, "_UMAP_feature_celltype_markers")
UMAP_feature
saveToPDF(paste0(path, ".pdf"), width = 12, height = 9)

# sample
UMAP_sample <- DimPlot(dataObject, reduction = "ref.umap", group.by = c("Sample_ID"))
path <- paste0("../results/integration_exploration/",projectID,"_UMAP_sample")
UMAP_sample
saveToPDF(paste0(path, ".pdf"), width = 9, height = 7)

# group
UMAP_group <- DimPlot(dataObject, reduction = "ref.umap", group.by = c("group"))
path <- paste0("../results/integration_exploration/",projectID,"_UMAP_group")
UMAP_group
saveToPDF(paste0(path, ".pdf"), width = 9, height = 7)

# violin features
v1 <- VlnPlot(dataObject, features = c("AQP4", "PLP1", "FLT1", "P2RY12","RBFOX1", "GAD1"), pt.size = 0.01, stack = TRUE, flip = TRUE,
  group.by = "predicted.subclass") + NoLegend() + ggtitle(paste0(projectID))
pdf(paste0("../results/integration_exploration/", projectID,"_violin_celltype_markers.pdf"),width = 12,height = 8)
v1
dev.off()

# DotPlot annotations
dot_anno <- DotPlot(dataObject, features = markers.to.plot,cluster.idents = TRUE)+ RotatedAxis()
pdf(paste0("../results/integration_exploration/", projectID,"_DotPlot_annotations.pdf"),width = 12,height = 6)
dot_anno
dev.off()

# DotPlot split by group
dot_group <- DotPlot(dataObject, features = markers.to.plot,cluster.idents = TRUE, split.by = "group", cols = ATSColors)+ RotatedAxis()
pdf(paste0("../results/integration_exploration/", projectID,"_DotPlot_group.pdf"),width = 12,height = 20)
dot_group
dev.off()

# Nuclei count per cluster
count_per_cluster <- FetchData(dataObject, vars = c("ident", "predicted.subclass")) %>%
  dplyr::count(ident, predicted.subclass) %>%
  tidyr::spread(ident, n)
count_per_cluster

count_melt <- reshape2::melt(count_per_cluster)
colnames(count_melt) <- c("ident", "cluster", "number of nuclei")
count_max <- count_melt[which.max(count_melt$`number of nuclei`), ]
count_max_value <- count_max$`number of nuclei`
cellmax <- count_max_value + 25000 # so that the figure doesn't cut off the text
count_bar <- ggplot(count_melt, aes(x = factor(cluster), y = `number of nuclei`, fill = `cluster`)) +
  geom_bar(
    stat = "identity",
    colour = "black",
    width = 1,
    position = position_dodge(width = 0.8)
  ) +
  geom_text(
    aes(label = `number of nuclei`),
    position = position_dodge(width = 0.9),
    vjust = -0.25,
    angle = 45,
    hjust = -.01
  ) +
  theme_classic() + scale_fill_manual(values = color.panel) +
  ggtitle("Number of nuclei per cluster") +  xlab("cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, cellmax))
count_bar
pdf(paste0("../results/integration_exploration/", projectID,"_nuclei_count_per_cluster.pdf"),width = 12,height = 6)
count_bar
dev.off()

# clean up
rm(gene_info, UMAP_doublets, UMAP_group, UMAP_sample, UMAP_feature, v1, count_bar, count_melt, count_max_value, count_max, count_per_cluster, dot_anno, dot_group)