## ----setup, include=FALSE--------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

# Libraris, paths, colors
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender_RPCAIntegration_dirty"
color.panel <- dittoColors()

# read object
dataObject  <- readRDS(paste0("../rObjects/",projectID,"_RNA.rds"))
dataObject[["RNA"]] <- JoinLayers(dataObject[["RNA"]])

## ----metadata-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
metadata <- subset(metadata, Sample_ID != "LBD_AS_F4")
metadata$sampleID <- factor(metadata$Sample_ID, levels = c(metadata$Sample_ID))
samples <- metadata$sampleID 
sex_order <- factor(metadata$sex_inferred, levels = unique(metadata$sex_inferred))
disease_order <- factor(metadata$TYPE, levels = c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS"))

metadata <- metadata %>%
  mutate(sampleID = gsub(".*_(\\d+)_.*_(BR_Nuclei).*", "\\2_\\1", Lane.Name))
samples <- metadata$sampleID 

# sampleID with disease_order
order <- metadata %>%
  arrange(disease_order) %>%
  dplyr::select(TYPE, sampleID, Sample_ID)
samples <- order$sampleID
disease_order <- order$TYPE
sample_order <- factor(order$Sample_ID, levels = order$Sample_ID)

seurat_sample_order <- as.character(dataObject$sample)
matched_metadata <- metadata[match(seurat_sample_order, as.character(metadata$sampleID)), ]

dataObject$Sample_ID <- factor(matched_metadata$Sample_ID, levels = sample_order)
table(dataObject$Sample_ID )

genes <- readRDS("../rObjects/annotation.rds")
mt.genes.df <- subset(genes, seqnames == "chrM")
mt.genes <- mt.genes.df$gene_name
## ----QC_metrics-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("Summary nCount")
summary(dataObject$nCount_RNA)
print("Summary nFeature")
summary(dataObject$nFeature_RNA)
# cell.complexity
dataObject$cell.complexity <- log10(dataObject$nFeature_RNA) / log10(dataObject$nCount_RNA)

# Chromosome M
gene.names <- rownames(dataObject)
dataObject$percent.mt <- PercentageFeatureSet(dataObject, features = mt.genes)
summary(dataObject$percent.mt)

## ----cells_per_sample-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Visualize the number of cell counts per sample
data <- as.data.frame(table(dataObject$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")

ncells1 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  #scale_fill_manual(values = sample_colors) + 
  scale_y_continuous(breaks = seq(0,30000, by = 2000), limits = c(0,30000)) +
  ggtitle("Nuclei per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
pdf(paste0("../results/dirty/", projectID,"_nuclei_count_per_sample.pdf"),width = 12,height = 4.5)
ncells1
dev.off()

mean_counts <- mean(data$frequency)
median(data$frequency)
sd_counts <- sd(data$frequency)

upper_threshold <- mean_counts + 2 * sd_counts
lower_threshold <- mean_counts - 2 * sd_counts

data$is_outlier <- ifelse(data$frequency > upper_threshold | data$frequency < lower_threshold, TRUE, FALSE)

# To view the outlier samples:
outlier_samples <- data[data$is_outlier == TRUE, ]
print(outlier_samples)

## ----Density-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df <- dataObject@meta.data
density_plot <- ggplot(df, aes(x = nCount_RNA, color = Sample_ID, fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  xlab("nCount_RNA") +
  ylab("Density") +
  ggtitle("Density nCount") +
  theme(axis.title.x = element_blank(), # Remove x-axis title for dirtyer alignment
        axis.text.x = element_blank(),   # Remove x-axis text as boxplot will provide labels
        axis.ticks.x = element_blank())  # Remove x-axis ticks

# Calculate summary statistics per Sample_ID
summary_stats <- df %>%
  group_by(Sample_ID) %>%
  summarise(
    ymin = min(nCount_RNA),
    ymax = max(nCount_RNA),
    Q1 = quantile(nCount_RNA, 0.25),
    Median = quantile(nCount_RNA, 0.5),
    Q3 = quantile(nCount_RNA, 0.75),
    Mean = mean(nCount_RNA)
  )

# Convert to log10 if your x-axis in the plot is log10
summary_stats_log10 <- summary_stats %>%
  mutate(across(c(ymin, ymax, Q1, Median, Q3, Mean), ~log10(.x)))
# For multiple boxplots:
box_plot_multiple <- ggplot(df, aes(x = nCount_RNA, y = Sample_ID, color = Sample_ID, fill = Sample_ID)) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA) +
  theme_classic() +
  scale_x_log10() +
  ggtitle("Boxplot nCount") +
  xlab("") + # No x-label here
  ylab("") + # No y-label here, as it will be provided by the density plot
  theme(legend.position = "none",
        axis.text.y = element_blank(), # Hide y-axis text to dirty up
        axis.ticks.y = element_blank(), # Hide y-axis ticks
        axis.title.y = element_blank(), # Hide y-axis title
        plot.margin = margin(b = 0)) # Reduce bottom margin
# For individual boxplots for each Sample_ID, arranged on top:
combined_plot <- box_plot_multiple / density_plot + plot_layout(heights = c(2, 3)) # Adjust heights as needed
pdf(paste0("../results/dirty/density/", projectID,"_nCount.pdf"),width = 8,height = 8)
print(combined_plot)
dev.off()

density_plot <- ggplot(df, aes(x = nFeature_RNA, color = Sample_ID, fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  xlab("nFeature_RNA") +
  ylab("Density") +
  ggtitle("Density nFeature") +
  theme(axis.title.x = element_blank(), # Remove x-axis title for dirtyer alignment
        axis.text.x = element_blank(),   # Remove x-axis text as boxplot will provide labels
        axis.ticks.x = element_blank())  # Remove x-axis ticks

# Calculate summary statistics per Sample_ID
summary_stats <- df %>%
  group_by(Sample_ID) %>%
  summarise(
    ymin = min(nFeature_RNA),
    ymax = max(nFeature_RNA),
    Q1 = quantile(nFeature_RNA, 0.25),
    Median = quantile(nFeature_RNA, 0.5),
    Q3 = quantile(nFeature_RNA, 0.75),
    Mean = mean(nFeature_RNA)
  )

# Convert to log10 if your x-axis in the plot is log10
summary_stats_log10 <- summary_stats %>%
  mutate(across(c(ymin, ymax, Q1, Median, Q3, Mean), ~log10(.x)))
# For multiple boxplots:
box_plot_multiple <- ggplot(df, aes(x = nFeature_RNA, y = Sample_ID, color = Sample_ID, fill = Sample_ID)) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA) +
  theme_classic() +
  scale_x_log10() +
  ggtitle("Boxplot nFeature") +
  xlab("") + # No x-label here
  ylab("") + # No y-label here, as it will be provided by the density plot
  theme(legend.position = "none",
        axis.text.y = element_blank(), # Hide y-axis text to dirty up
        axis.ticks.y = element_blank(), # Hide y-axis ticks
        axis.title.y = element_blank(), # Hide y-axis title
        plot.margin = margin(b = 0)) # Reduce bottom margin
# For individual boxplots for each Sample_ID, arranged on top:
combined_plot <- box_plot_multiple / density_plot + plot_layout(heights = c(2, 3)) # Adjust heights as needed
pdf(paste0("../results/dirty/density/", projectID,"_nFeature.pdf"),width = 8,height = 8)
print(combined_plot)
dev.off()

density_plot <- ggplot(df, aes(x = percent.mt, color = Sample_ID, fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  xlab("percent.mt") +
  ylab("Density") +
  ggtitle("Density percent.mt") +
  theme(axis.title.x = element_blank(), # Remove x-axis title for dirtyer alignment
        axis.text.x = element_blank(),   # Remove x-axis text as boxplot will provide labels
        axis.ticks.x = element_blank())  # Remove x-axis ticks

# Calculate summary statistics per Sample_ID
summary_stats <- df %>%
  group_by(Sample_ID) %>%
  summarise(
    ymin = min(percent.mt),
    ymax = max(percent.mt),
    Q1 = quantile(percent.mt, 0.25),
    Median = quantile(percent.mt, 0.5),
    Q3 = quantile(percent.mt, 0.75),
    Mean = mean(percent.mt)
  )
# Convert to log10 if your x-axis in the plot is log10
summary_stats_log10 <- summary_stats %>%
  mutate(across(c(ymin, ymax, Q1, Median, Q3, Mean),))
# For multiple boxplots:
box_plot_multiple <- ggplot(df, aes(x = percent.mt, y = Sample_ID, color = Sample_ID, fill = Sample_ID)) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA) +
  theme_classic() +
  ggtitle("Boxplot percent.mt") +
  xlab("") + # No x-label here
  ylab("") + # No y-label here, as it will be provided by the density plot
  theme(legend.position = "none",
        axis.text.y = element_blank(), # Hide y-axis text to dirty up
        axis.ticks.y = element_blank(), # Hide y-axis ticks
        axis.title.y = element_blank(), # Hide y-axis title
        plot.margin = margin(b = 0)) # Reduce bottom margin
# For individual boxplots for each Sample_ID, arranged on top:
combined_plot <- box_plot_multiple / density_plot + plot_layout(heights = c(2, 3)) # Adjust heights as needed
pdf(paste0("../results/dirty/density/", projectID,"_mt.pdf"),width = 8,height = 8)
print(combined_plot)
dev.off()

## ----Preprocess_for_integration----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dataObject[["RNA"]] <- split(dataObject[["RNA"]], f = dataObject$Sample_ID)

# since the data is split into layers, normalization and variable feature identification is performed for each sample independently (a consensus set of variable features is automatically identified).
dataObject # inspect
# Layers
print("Inspect Layers")
Layers(dataObject[["RNA"]])

dataObject <- SCTransform(dataObject, verbose = TRUE, conserve.memory = TRUE) # this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
print("Insepect SCT")
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
projectID <- "CWOW_cellbender_SCTransform_dirty"
saveRDS(dataObject, paste0("../rObjects/",projectID,".rds"), compress = FALSE)

#--------------------
print("RPCAIntegration")
projectID <- "CWOW_cellbender_RPCAIntegration_dirty"

dataObject <- RunPCA(dataObject, npcs = 30, verbose = F)
features <- dataObject@assays$SCT@var.features
write.table(features, paste0("../results/dirty/",projectID,"_SCT_variable_features.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

dataObject <- IntegrateLayers(object = dataObject, orig.reduction = "pca", method = RPCAIntegration, normalization.method = "SCT", verbose = F, new.reduction = "integrated.rpca")
dataObject <- FindNeighbors(dataObject, reduction = "integrated.rpca", dims = 1:30)
dataObject <- FindClusters(dataObject, resolution = 0.6, cluster.name = "rpca_clusters")
dataObject <- RunUMAP(dataObject, dims = 1:30, reduction = "integrated.rpca", reduction.name = "integrated.rpca")

projectID <- "CWOW_cellbender_RPCAIntegration_dirty"
## ----save_object,echo=FALSE,eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(dataObject, paste0("../rObjects/",projectID,"_SCT.rds"), compress = FALSE)
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

pdf(paste0("../results/dirty/",projectID,"_UMAP_annotations_and_group.pdf"), width = 14, height = 7)
DimPlot(dataObject, reduction = "integrated.rpca", group.by = c("group", "predicted.subclass")) # disease and azimuth annotations
dev.off()

ditto_umap <- dittoDimPlot(object = dataObject, var = "predicted.subclass", reduction.use = "integrated.rpca", do.label = TRUE, labels.highlight = TRUE)
pdf(paste0("../results/dirty/",projectID,"_UMAP_annotations.pdf"), width = 8, height = 7)
ditto_umap
dev.off()

# doublet status
UMAP_doublets <- DimPlot(dataObject, reduction = "integrated.rpca", group.by = c("CellTypes_DF"), cols=c("black", "#66C2A5"))
pdf(paste0("../results/dirty/",projectID,"_UMAP_doubletStatus.pdf"), width = 8, height = 7)
UMAP_doublets
dev.off()

# feature plot 
UMAP_feature <- FeaturePlot(dataObject,  reduction = "integrated.rpca", features = c("AQP4", "PLP1", "RBFOX3", "GAD1"))
pdf(paste0("../results/dirty/",projectID, "_UMAP_feature_celltype_markers.pdf"), width = 12, height = 9)
UMAP_feature
dev.off()

# sample
UMAP_sample <- DimPlot(dataObject, reduction = "integrated.rpca", group.by = c("Sample_ID"))
pdf(paste0("../results/dirty/",projectID,"_UMAP_sample.pdf"), width = 9, height = 7)
UMAP_sample
dev.off()

# group
UMAP_group <- DimPlot(dataObject, reduction = "integrated.rpca", group.by = c("group"))
pdf(paste0("../results/dirty/",projectID,"_UMAP_group.pdf"), width = 9, height = 7)
UMAP_group
dev.off()

# violin features predicted.subclass
v1 <- VlnPlot(dataObject, features = c("AQP4", "PLP1", "FLT1", "P2RY12","RBFOX1", "GAD1"), pt.size = 0.01, stack = TRUE, flip = TRUE,
              group.by = "predicted.subclass") + NoLegend() + ggtitle(paste0(projectID))
pdf(paste0("../results/dirty/", projectID,"_violin_celltype_markers_azimuth.pdf"),width = 12,height = 8)
v1
dev.off()

# violin features clusters
v1 <- VlnPlot(dataObject, features = c("AQP4", "PLP1", "FLT1", "P2RY12","RBFOX1", "GAD1"), pt.size = 0.01, stack = TRUE, flip = TRUE,
              group.by = "rpca_clusters") + NoLegend() + ggtitle(paste0(projectID))
pdf(paste0("../results/dirty/", projectID,"_violin_celltype_markers_by_cluster.pdf"),width = 12,height = 8)
v1
dev.off()

# DotPlot clusters
Idents(dataObject) <- dataObject$rpca_clusters
dot_clusters <- DotPlot(dataObject, features = markers.to.plot, cluster.idents = TRUE)+ RotatedAxis()
pdf(paste0("../results/dirty/", projectID,"_DotPlot_clusters.pdf"),width = 12,height = 15)
dot_clusters
dev.off()

Idents(dataObject) <- dataObject$predicted.subclass
# DotPlot annotations
dot_anno <- DotPlot(dataObject, features = markers.to.plot, cluster.idents = TRUE)+ RotatedAxis()
pdf(paste0("../results/dirty/", projectID,"_DotPlot_annotations.pdf"),width = 12,height = 6)
dot_anno
dev.off()

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
  theme_classic() + scale_fill_manual(values = color.panel) +
  ggtitle("Number of nuclei per cluster and group") +  xlab("cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, cellmax))
count_bar
pdf(paste0("../results/dirty/", projectID,"_nuclei_count_per_cluster_group.pdf"),width = 12,height = 6)
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
  theme_classic() + scale_fill_manual(values = color.panel) +
  ggtitle("Number of nuclei per cluster") +  xlab("cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, cellmax))
count_bar
pdf(paste0("../results/dirty/", projectID,"_nuclei_count_per_cluster.pdf"),width = 12,height = 6)
count_bar
dev.off()
## ----End-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------