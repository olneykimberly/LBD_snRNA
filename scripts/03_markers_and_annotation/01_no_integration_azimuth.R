## ----setup-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

#https://satijalab.org/seurat/articles/integration_introduction
#https://satijalab.org/seurat/articles/seurat5_integration 

## ----source------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
project_ID <- "CWOW_cellbender_singlets" #CWOW_cellbender_mereged_singlets_with_kept_doublets
color_panel <- dittoColors()

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


## ----JoinLayers_AZIMUTH--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# update project_ID
project_ID <- "CWOW_cellbender_joinlayers_azimuth"
## We will obtain predicted cell annotations via AZIMUTH for comparing integration methods 
dataObject[["RNA"]] <- JoinLayers(dataObject[["RNA"]]) # first join layers
# Azimuth annotations
dataObject <- RunAzimuth(dataObject, reference = "humancortexref")
saveRDS(dataObject, paste0("../rObjects/",project_ID,".rds"), compress = FALSE)
dataObject # inspect
Idents(dataObject) <- dataObject$predicted.subclass # set idents

# annotation
ditto_umap <- dittoDimPlot(object = dataObject, var = "predicted.subclass", reduction.use = "ref.umap", do.label = TRUE, labels.highlight = TRUE)
pdf(paste0("../results/unintegrated_exploration/",project_ID,"_UMAP_annotations.pdf"), width = 8, height = 7)
ditto_umap
dev.off()

# doublet status
#UMAP_doublets <- DimPlot(dataObject, reduction = "ref.umap", group.by = c("CellTypes_DF"), cols=c("black", "#66C2A5"))
#path <- paste0("../results/unintegrated_exploration/",project_ID,"_UMAP_doubletStatus")
#UMAP_doublets
#saveToPDF(paste0(path, ".pdf"), width = 8, height = 7)

# feature plot 
UMAP_feature <- FeaturePlot(dataObject,  reduction = "ref.umap", features = c("AQP4", "PLP1", "RBFOX3", "GAD1"))
pdf(paste0("../results/unintegrated_exploration/",project_ID, "_UMAP_feature_celltype_markers.pdf"), width = 12, height = 9)
UMAP_feature
dev.off()

# sample
UMAP_sample <- DimPlot(dataObject, reduction = "ref.umap", group.by = c("Sample_ID"))
pdf(paste0("../results/unintegrated_exploration/",project_ID,"_UMAP_sample.pdf"), width = 9, height = 7)
UMAP_sample
dev.off()

# group
UMAP_group <- DimPlot(dataObject, reduction = "ref.umap", group.by = c("group"))
pdf(paste0("../results/unintegrated_exploration/",project_ID,"_UMAP_group.pdf"),width = 9, height = 7)
UMAP_group
dev.off()

# violin features
v1 <- VlnPlot(dataObject, features = c("AQP4", "PLP1", "FLT1", "P2RY12","RBFOX1", "GAD1"), pt.size = 0.01, stack = TRUE, flip = TRUE,
  group.by = "predicted.subclass") + NoLegend() + ggtitle(paste0(project_ID))
pdf(paste0("../results/unintegrated_exploration/", project_ID,"_violin_celltype_markers.pdf"),width = 12,height = 8)
v1
dev.off()

# DotPlot annotations
dot_anno <- DotPlot(dataObject, features = genes_markers,cluster.idents = TRUE)+ RotatedAxis()
pdf(paste0("../results/unintegrated_exploration/", project_ID,"_DotPlot_annotations.pdf"),width = 12,height = 6)
dot_anno
dev.off()

# DotPlot split by group
dot_group <- DotPlot(dataObject, features = genes_markers,cluster.idents = TRUE, split.by = "group", cols = color_ATS)+ RotatedAxis()
pdf(paste0("../results/unintegrated_exploration/", project_ID,"_DotPlot_group.pdf"),width = 12,height = 20)
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
  theme_classic() + scale_fill_manual(values = color_panel) +
  ggtitle("Number of nuclei per cluster") +  xlab("cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, cellmax))
count_bar
pdf(paste0("../results/unintegrated_exploration/", project_ID,"_nuclei_count_per_cluster.pdf"),width = 12,height = 6)
count_bar
dev.off()

# clean up
rm(gene_info, UMAP_doublets, UMAP_group, UMAP_sample, UMAP_feature, v1, count_bar, count_melt, count_max_value, count_max, count_per_cluster, dot_anno, dot_group)