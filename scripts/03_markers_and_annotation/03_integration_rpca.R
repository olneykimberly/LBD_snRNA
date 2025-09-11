# Integration 
# Identify cell subpopulations that are present among all samples
# Obtain cell type markers that are conserved among groups
# Compare the samples to find cell-type specific responses to disease group

# Aim to integrate the disease conditions together, so that we can jointly identify cell subpopulations across disease group.
# Then explore how each group differs across conditions

## ----setup-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

#https://satijalab.org/seurat/articles/integration_introduction
#https://satijalab.org/seurat/articles/seurat5_integration 

## ----source------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
project_ID <- "CWOW_cellbender"
color_panel <- dittoColors()

## ----read_object--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read object
dataObject <- readRDS(file = paste0("../rObjects/", project_ID, "_SCTransform.rds"))
dataObject # inspect

## ----Perform_integration----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# We now aim to integrate data so that cells from the same cell type/subpopulation will cluster together.
# We often refer to this procedure as intergration/alignment. When aligning two genome sequences together, identification of shared/homologous regions can help to interpret differences between the sequences as well. 
# Similarly for scRNA-seq integration, our goal is not to remove biological differences across conditions, but to learn shared cell types/states in an initial step - specifically because that will enable us to compare among disease profiles for these individual cell types.

# The Seurat v5 integration procedure aims to return a single dimensional reduction that captures the shared sources of variance across multiple layers, so that cells in a similar biological state will cluster. 
# The method returns a dimensional reduction (i.e. integrated.cca) which can be used for visualization and unsupervised clustering analysis. 
# For evaluating performance, we can use cell type labels from the AZIMUTH analysis that we ran in the prior script metadata column.

#---- RPCAIntegration
print("RPCAIntegration")
dataObject <- RunPCA(dataObject, npcs = 30, verbose = F)
features <- dataObject@assays$SCT@var.features
write.table(features, paste0("../results/",project_ID,"_SCT_variable_features.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## add assay name 
dataObject <- IntegrateLayers(object = dataObject, orig.reduction = "pca", method = RPCAIntegration, normalization.method = "SCT", verbose = F, new.reduction = "integrated.rpca")
dataObject <- FindNeighbors(dataObject, reduction = "integrated.rpca", dims = 1:30)
dataObject <- FindClusters(dataObject, resolution = 0.6, cluster.name = "rpca_clusters")
dataObject <- RunUMAP(dataObject, dims = 1:30, reduction = "integrated.rpca", reduction.name = "integrated.rpca")

project_ID <- "CWOW_cellbender_RPCAIntegration"

## ----save_object,echo=FALSE,eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(dataObject, paste0("../rObjects/",project_ID,".rds"), compress = FALSE)
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# pdf(paste0("../results/UMAP/",project_ID,"_UMAP_annotations_and_group.pdf"), width = 14, height = 7)
# DimPlot(dataObject, reduction = "integrated.rpca", group.by = c("group", "predicted.subclass")) # disease and azimuth annotations
# dev.off()
# 
# ditto_umap <- dittoDimPlot(object = dataObject, var = "predicted.subclass", reduction.use = "integrated.rpca", do.label = TRUE, labels.highlight = TRUE)
# pdf(paste0("../results/UMAP/",project_ID,"_UMAP_annotations.pdf"), width = 8, height = 7)
# ditto_umap
# dev.off()

# doublet status
#UMAP_doublets <- DimPlot(dataObject, reduction = "integrated.rpca", group.by = c("CellTypes_DF"), cols=c("black", "#66C2A5"))
#pdf(paste0("../results/UMAP/",project_ID,"_UMAP_doubletStatus.pdf"), width = 8, height = 7)
#UMAP_doublets
#dev.off()

# feature plot 
UMAP_feature <- FeaturePlot(dataObject,  reduction = "integrated.rpca", features = c("AQP4", "PLP1", "RBFOX3", "GAD1"))
pdf(paste0("../results/UMAP/",project_ID, "_UMAP_feature_celltype_markers.pdf"), width = 12, height = 9)
UMAP_feature
dev.off()

# sample
UMAP_sample <- DimPlot(dataObject, reduction = "integrated.rpca", group.by = c("Sample_ID"))
pdf(paste0("../results/UMAP/",project_ID,"_UMAP_sample.pdf"), width = 9, height = 7)
UMAP_sample
dev.off()

# group
UMAP_group <- DimPlot(dataObject, reduction = "integrated.rpca", group.by = c("group"))
pdf(paste0("../results/UMAP/",project_ID,"_UMAP_group.pdf"), width = 9, height = 7)
UMAP_group
dev.off()

# violin features predicted.subclass
# v1 <- VlnPlot(dataObject, features = c("AQP4", "PLP1", "FLT1", "P2RY12","RBFOX1", "GAD1"), pt.size = 0.01, stack = TRUE, flip = TRUE,
#               group.by = "predicted.subclass") + NoLegend() + ggtitle(paste0(project_ID))
# pdf(paste0("../results/violin/", project_ID,"_violin_celltype_markers_azimuth.pdf"),width = 12,height = 8)
# v1
# dev.off()

# violin features clusters
v1 <- VlnPlot(dataObject, features = c("AQP4", "PLP1", "FLT1", "P2RY12","RBFOX1", "GAD1"), pt.size = 0.01, stack = TRUE, flip = TRUE,
              group.by = "rpca_clusters") + NoLegend() + ggtitle(paste0(project_ID))
pdf(paste0("../results/violin/", project_ID,"_violin_celltype_markers_by_cluster.pdf"),width = 12,height = 8)
v1
dev.off()

# DotPlot clusters
Idents(dataObject) <- dataObject$rpca_clusters
dot_clusters <- DotPlot(dataObject, features = genes_markers, cluster.idents = TRUE)+ RotatedAxis()
pdf(paste0("../results/dot_plot/", project_ID,"_DotPlot_clusters.pdf"),width = 12,height = 15)
dot_clusters
dev.off()

# Idents(dataObject) <- dataObject$predicted.subclass
# # DotPlot annotations
# dot_anno <- DotPlot(dataObject, features = genes_markers, cluster.idents = TRUE)+ RotatedAxis()
# pdf(paste0("../results/dot_plot/", project_ID,"_DotPlot_annotations.pdf"),width = 12,height = 6)
# dot_anno
# dev.off()

# Nuclei count per cluster
count_per_cluster <- FetchData(dataObject, vars = c("group", "rpca_clusters")) %>%
  dplyr::count(group, rpca_clusters) %>%
  tidyr::spread(group, n)
count_per_cluster

count_melt <- reshape2::melt(count_per_cluster)
colnames(count_melt) <- c("group", "cluster", "number of nuclei")
count_max <- count_melt[which.max(count_melt$`number of nuclei`), ]
count_max_value <- count_max$`number of nuclei`
cellmax <- count_max_value + 250 # so that the figure doesn't cut off the text
count_bar <- ggplot(count_melt, aes(x = factor(group), y = `number of nuclei`, fill = `cluster`)) +
  geom_bar(
    stat = "identity",
    colour = "black",
    width = 1,
    position = position_dodge(width = 0.8)
  ) +
  theme_classic() + scale_fill_manual(values = color_panel) +
  ggtitle("Number of nuclei per cluster and group") +  xlab("cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, cellmax))
count_bar
pdf(paste0("../results/nuclei_count/", project_ID,"_nuclei_count_per_cluster_group.pdf"),width = 12,height = 6)
count_bar
dev.off()

count_per_cluster <- FetchData(dataObject, vars = c("rpca_clusters")) %>%
  dplyr::count(rpca_clusters)
colnames(count_per_cluster) <- c("rpca_clusters", "number of nuclei")
count_max <- count_per_cluster[which.max(count_per_cluster$`number of nuclei`), ]
count_max_value <- count_max$`number of nuclei`
cellmax <- count_max_value + 2500 # so that the figure doesn't cut off the text
count_bar <- ggplot(count_per_cluster, aes(x = factor(rpca_clusters), y = `number of nuclei`, fill = `rpca_clusters`)) +
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
  theme_classic() + scale_fill_manual(values = color_panel) +
  ggtitle("Number of nuclei per cluster") +  xlab("cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, cellmax))
count_bar
pdf(paste0("../results/nuclei_count/", project_ID,"_nuclei_count_per_cluster.pdf"),width = 12,height = 6)
count_bar
dev.off()
