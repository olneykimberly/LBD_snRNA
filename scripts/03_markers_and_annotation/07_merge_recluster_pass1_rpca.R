## ----setup, include=FALSE--------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

# Libraris, paths, colors
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender_RPCAIntegration"
color.panel <- dittoColors()
set.seed(28)

#---- get markers 
# markers <- readRDS(paste0("../rObjects/", projectID,"_markers.rds"))
# all.markers.strict <- readRDS(paste0("../rObjects/", projectID,"_markers_log2FC1_q0.01.rds"))
# # Get markers for each cluster
# # unique clusters variable
# unique_clusters <- unique(markers$cluster)
# # empty list to store individual cluster data frames
# cluster_list <- list()
# # loop through each cluster and create a data frame
# for (i in unique_clusters) {
#   cluster_name <- paste0("cluster", i)
#   cluster_data <- all.markers.strict[all.markers.strict$cluster == i, ]
#   assign(cluster_name, cluster_data)
#   cluster_list[[cluster_name]] <- cluster_data
# }

## ----read_object-----------------------------------------------------------------------------------------------------
# read object
dataObject.astrocyte <- readRDS(paste0("../rObjects/", projectID, "_astrocyte.rds"))
dataObject.oligodendrocyte <- readRDS(paste0("../rObjects/",projectID,"_oligodendrocyte.rds"))
dataObject.opc <- readRDS(paste0("../rObjects/",projectID,"_polydencrocyte.rds"))
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

dir_path <- "rpca_pass1_removals_before_merge/"

#------------------------------ astrocyte
celltype <- "astrocyte"

print("astrocyte nuclei total count before filter")
sum(summary(dataObject.astrocyte$seurat_clusters))

data <- as.data.frame(table(dataObject.astrocyte$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells1 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,3500, by = 500), limits = c(0,3500)) +
  ggtitle(paste0(celltype, " nuclei count per sample before percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)
# remove dirty clusters 
dataObject.astrocyte.clean <- subset(dataObject.astrocyte, cells = WhichCells(dataObject.astrocyte, idents = setdiff(unique(dataObject.astrocyte$seurat_clusters), c(19, 17, 2, 5, 12))))
dataObject.astrocyte.dirty <- subset(dataObject.astrocyte, cells = WhichCells(dataObject.astrocyte, idents = c(19, 17, 2, 5, 12)))

v1 <- VlnPlot(dataObject.astrocyte.clean,
        features = c("AQP4", "GJA1", "PLP1", "MBP", "RBFOX1"),
        stack = TRUE,
        flip = TRUE) +ggtitle("Before percentage filter")

pdf(paste0("../results/", dir_path,  celltype, "_violin_before_percent_filter.pdf"), width = 7, height = 5)
v1
dev.off()

astrocyte.genes <- c("AQP4", "GJA1")
contaminate.genes <- c("RBFOX1", "MBP", "PLP1")
dataObject.astrocyte.clean$percent.astrocyte <- PercentageFeatureSet(dataObject.astrocyte.clean, features = astrocyte.genes)
dataObject.astrocyte.clean$percent.contaminate <- PercentageFeatureSet(dataObject.astrocyte.clean, features = contaminate.genes)

print("summary of percent astrocyte genes")
summary(dataObject.astrocyte.clean$percent.astrocyte)
print("summary of percent contaminate genes")
summary(dataObject.astrocyte.clean$percent.contaminate)

prev_nuclei <- ncol(dataObject.astrocyte.clean)
dataObject.astrocyte.clean <- subset(dataObject.astrocyte.clean, subset = percent.astrocyte > 0 & percent.contaminate == 0)
removed_count <- prev_nuclei - ncol(dataObject.astrocyte.clean)
print(paste0("removed nuclei based on percentage: ", removed_count))

v2 <- VlnPlot(dataObject.astrocyte.clean,
              features = c("AQP4", "GJA1"),
              stack = TRUE,
              flip = TRUE) +ggtitle("After percentage filter")

pdf(paste0("../results/", dir_path, celltype, "_violin_after_percent_filter.pdf"), width = 7, height = 3.5)
v2
dev.off()

data <- as.data.frame(table(dataObject.astrocyte.clean$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells2 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,3500, by = 500), limits = c(0,3500)) +
  ggtitle(paste0(celltype, " nuclei count per sample after percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", dir_path, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 12, height = 7)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("astrocyte nuclei total count after filter")
sum(summary(dataObject.astrocyte.clean$seurat_clusters))

# clean up 
rm(dataObject.astrocyte, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2)

#------------------------------ oligodendrocyte
celltype <- "oligodendrocyte"

print("oligodendrocyte nuclei total count before filter")
sum(summary(dataObject.oligodendrocyte$seurat_clusters))

data <- as.data.frame(table(dataObject.oligodendrocyte$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells1 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,20000, by = 1500), limits = c(0,20000)) +
  ggtitle(paste0(celltype, " nuclei count per sample before percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)
# remove dirty clusters 
dataObject.oligodendrocyte.clean <- subset(dataObject.oligodendrocyte, cells = WhichCells(dataObject.oligodendrocyte, idents = setdiff(unique(dataObject.oligodendrocyte$seurat_clusters), c(19, 20, 22, 4))))
dataObject.oligodendrocyte.dirty <- subset(dataObject.oligodendrocyte, cells = WhichCells(dataObject.oligodendrocyte, idents = c(19, 20, 22, 4)))

v1 <- VlnPlot(dataObject.oligodendrocyte.clean,
              features = c("PLP1", "MBP", "RBFOX1", "SNAP25"),
              stack = TRUE,
              flip = TRUE) +ggtitle("Before percentage filter")

pdf(paste0("../results/", dir_path,  celltype, "_violin_before_percent_filter.pdf"), width = 7, height = 5)
v1
dev.off()

oligodendrocyte.genes <- c("PLP1", "MBP")
contaminate.genes <- c("RBFOX1", "SNAP25")
dataObject.oligodendrocyte.clean$percent.oligodendrocyte <- PercentageFeatureSet(dataObject.oligodendrocyte.clean, features = oligodendrocyte.genes)
dataObject.oligodendrocyte.clean$percent.contaminate <- PercentageFeatureSet(dataObject.oligodendrocyte.clean, features = contaminate.genes)

print("summary of percent oligodendrocyte genes")
summary(dataObject.oligodendrocyte.clean$percent.oligodendrocyte)
print("summary of percent contaminate genes")
summary(dataObject.oligodendrocyte.clean$percent.contaminate)

prev_nuclei <- ncol(dataObject.oligodendrocyte.clean)
dataObject.oligodendrocyte.clean <- subset(dataObject.oligodendrocyte.clean, subset = percent.oligodendrocyte > 0 & percent.contaminate == 0)
removed_count <- prev_nuclei - ncol(dataObject.oligodendrocyte.clean)
print(paste0("removed nuclei based on percentage: ", removed_count))

v2 <- VlnPlot(dataObject.oligodendrocyte.clean,
              features = c("PLP1", "MBP"),
              stack = TRUE,
              flip = TRUE) +ggtitle("After percentage filter")

pdf(paste0("../results/", dir_path, celltype, "_violin_after_percent_filter.pdf"), width = 7, height = 3.5)
v2
dev.off()

data <- as.data.frame(table(dataObject.oligodendrocyte.clean$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells2 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,20000, by = 1500), limits = c(0,20000)) +
  ggtitle(paste0(celltype, " nuclei count per sample after percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", dir_path, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 12, height = 7)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("oligodendrocyte nuclei total count after filter")
sum(summary(dataObject.oligodendrocyte.clean$seurat_clusters))

# clean up 
rm(dataObject.oligodendrocyte, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2)

#------------------------------ microglia
celltype <- "microglia"

print("microglia nuclei total count before filter")
sum(summary(dataObject.microglia$seurat_clusters))

data <- as.data.frame(table(dataObject.microglia$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells1 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,3500, by = 500), limits = c(0,3500)) +
  ggtitle(paste0(celltype, " nuclei count per sample before percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)
# remove dirty clusters 
dataObject.microglia.clean <- subset(dataObject.microglia, cells = WhichCells(dataObject.microglia, idents = setdiff(unique(dataObject.microglia$seurat_clusters), c(13, 5, 4))))
dataObject.microglia.dirty <- subset(dataObject.microglia, cells = WhichCells(dataObject.microglia, idents = c(13, 5, 4)))

v1 <- VlnPlot(dataObject.microglia.clean,
              features = c("DOCK8", "P2RY12", "PLP1", "MBP", "RBFOX1"),
              stack = TRUE,
              flip = TRUE) +ggtitle("Before percentage filter")

pdf(paste0("../results/", dir_path,  celltype, "_violin_before_percent_filter.pdf"), width = 7, height = 5)
v1
dev.off()

microglia.genes <- c("DOCK8", "P2RY12")
contaminate.genes <- c("RBFOX1", "PLP1")
dataObject.microglia.clean$percent.microglia <- PercentageFeatureSet(dataObject.microglia.clean, features = microglia.genes)
dataObject.microglia.clean$percent.contaminate <- PercentageFeatureSet(dataObject.microglia.clean, features = contaminate.genes)

print("summary of percent microglia genes")
summary(dataObject.microglia.clean$percent.microglia)
print("summary of percent contaminate genes")
summary(dataObject.microglia.clean$percent.contaminate)

prev_nuclei <- ncol(dataObject.microglia.clean)
dataObject.microglia.clean <- subset(dataObject.microglia.clean, subset = percent.microglia > 0 & percent.contaminate == 0)
removed_count <- prev_nuclei - ncol(dataObject.microglia.clean)
print(paste0("removed nuclei based on percentage: ", removed_count))

v2 <- VlnPlot(dataObject.microglia.clean,
              features = c("DOCK8", "P2RY12"),
              stack = TRUE,
              flip = TRUE) +ggtitle("After percentage filter")

pdf(paste0("../results/", dir_path, celltype, "_violin_after_percent_filter.pdf"), width = 7, height = 3.5)
v2
dev.off()

data <- as.data.frame(table(dataObject.microglia.clean$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells2 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,3500, by = 500), limits = c(0,3500)) +
  ggtitle(paste0(celltype, " nuclei count per sample after percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", dir_path, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 12, height = 7)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("microglia nuclei total count after filter")
sum(summary(dataObject.microglia.clean$seurat_clusters))

# clean up 
rm(dataObject.microglia, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2)
#------------------------------ mural
celltype <- "mural"

print("mural nuclei total count before filter")
sum(summary(dataObject.mural$seurat_clusters))

data <- as.data.frame(table(dataObject.mural$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells1 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,500, by = 100), limits = c(0,500)) +
  ggtitle(paste0(celltype, " nuclei count per sample before percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)
# remove dirty clusters 
dataObject.mural.clean <- subset(dataObject.mural, cells = WhichCells(dataObject.mural, idents = setdiff(unique(dataObject.mural$seurat_clusters), c(2,5))))
dataObject.mural.dirty <- subset(dataObject.mural, cells = WhichCells(dataObject.mural, idents = c(2,5)))

v1 <- VlnPlot(dataObject.mural.clean,
              features = c("EBF1", "RBPMS", "PLP1", "MBP", "RBFOX1"),
              stack = TRUE,
              flip = TRUE) +ggtitle("Before percentage filter")

pdf(paste0("../results/", dir_path,  celltype, "_violin_before_percent_filter.pdf"), width = 7, height = 5)
v1
dev.off()

mural.genes <- c("EBF1", "RBPMS")
contaminate.genes <- c("RBFOX1", "PLP1", "MBP")
dataObject.mural.clean$percent.mural <- PercentageFeatureSet(dataObject.mural.clean, features = mural.genes)
dataObject.mural.clean$percent.contaminate <- PercentageFeatureSet(dataObject.mural.clean, features = contaminate.genes)

print("summary of percent mural genes")
summary(dataObject.mural.clean$percent.mural)
print("summary of percent contaminate genes")
summary(dataObject.mural.clean$percent.contaminate)

prev_nuclei <- ncol(dataObject.mural.clean)
dataObject.mural.clean <- subset(dataObject.mural.clean, subset = percent.mural > 0 & percent.contaminate == 0)
removed_count <- prev_nuclei - ncol(dataObject.mural.clean)
print(paste0("removed nuclei based on percentage: ", removed_count))

v2 <- VlnPlot(dataObject.mural.clean,
              features = c("EBF1", "RBPMS"),
              stack = TRUE,
              flip = TRUE) +ggtitle("After percentage filter")

pdf(paste0("../results/", dir_path, celltype, "_violin_after_percent_filter.pdf"), width = 7, height = 3.5)
v2
dev.off()

data <- as.data.frame(table(dataObject.mural.clean$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells2 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,500, by = 100), limits = c(0,500)) +
  ggtitle(paste0(celltype, " nuclei count per sample after percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", dir_path, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 12, height = 7)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("mural nuclei total count after filter")
sum(summary(dataObject.mural.clean$seurat_clusters))

# clean up 
rm(dataObject.mural, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2)
#------------------------------ opc
celltype <- "opc"

print("opc nuclei total count before filter")
sum(summary(dataObject.opc$seurat_clusters))

data <- as.data.frame(table(dataObject.opc$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells1 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,2500, by = 500), limits = c(0,2500)) +
  ggtitle(paste0(celltype, " nuclei count per sample before percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)
# remove dirty clusters 
dataObject.opc.clean <- subset(dataObject.opc, cells = WhichCells(dataObject.opc, idents = setdiff(unique(dataObject.opc$seurat_clusters), c(13, 14, 15, 17, 4))))
dataObject.opc.dirty <- subset(dataObject.opc, cells = WhichCells(dataObject.opc, idents = c(13, 14, 15, 17, 4)))

v1 <- VlnPlot(dataObject.opc.clean,
              features = c("OLIG1", "VCAN", "PLP1", "MBP", "RBFOX1", "SNAP25"),
              stack = TRUE,
              flip = TRUE) +ggtitle("Before percentage filter")

pdf(paste0("../results/", dir_path,  celltype, "_violin_before_percent_filter.pdf"), width = 7, height = 6)
v1
dev.off()

opc.genes <- c("OLIG1", "VCAN")
contaminate.genes <- c("PLP1")
dataObject.opc.clean$percent.opc <- PercentageFeatureSet(dataObject.opc.clean, features = opc.genes)
dataObject.opc.clean$percent.contaminate <- PercentageFeatureSet(dataObject.opc.clean, features = contaminate.genes)

print("summary of percent opc genes")
summary(dataObject.opc.clean$percent.opc)
print("summary of percent contaminate genes")
summary(dataObject.opc.clean$percent.contaminate)

prev_nuclei <- ncol(dataObject.opc.clean)
dataObject.opc.clean <- subset(dataObject.opc.clean, subset = percent.opc > 0 & percent.contaminate == 0)
removed_count <- prev_nuclei - ncol(dataObject.opc.clean)
print(paste0("removed nuclei based on percentage: ", removed_count))

v2 <- VlnPlot(dataObject.opc.clean,
              features = c("OLIG1", "VCAN", "MBP", "RBFOX1", "SNAP25"),
              stack = TRUE,
              flip = TRUE) +ggtitle("After percentage filter")

pdf(paste0("../results/", dir_path, celltype, "_violin_after_percent_filter.pdf"), width = 7, height = 5)
v2
dev.off()

data <- as.data.frame(table(dataObject.opc.clean$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells2 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,2500, by = 1500), limits = c(0,2500)) +
  ggtitle(paste0(celltype, " nuclei count per sample after percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", dir_path, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 12, height = 7)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("opc nuclei total count after filter")
sum(summary(dataObject.opc.clean$seurat_clusters))

# clean up 
rm(dataObject.opc, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2)

#------------------------------ neuron
celltype <- "neuron"

print("neuron nuclei total count before filter")
sum(summary(dataObject.neuron$seurat_clusters))

data <- as.data.frame(table(dataObject.neuron$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells1 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,10500, by = 1500), limits = c(0,10500)) +
  ggtitle(paste0(celltype, " nuclei count per sample before percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)

DotPlot(dataObject.neuron, features = markers.to.plot, cluster.idents = FALSE, dot.scale = 8) + RotatedAxis()
# remove dirty clusters 
dataObject.neuron.clean <- subset(dataObject.neuron, cells = WhichCells(dataObject.neuron, idents = setdiff(unique(dataObject.neuron$seurat_clusters), c(7, 10, 14, 21, 26, 33, 35, 38))))
dataObject.neuron.dirty <- subset(dataObject.neuron, cells = WhichCells(dataObject.neuron, idents = c(7, 10, 14, 21, 26, 33, 35, 38)))

v1 <- VlnPlot(dataObject.neuron.clean,
              features = c("RBFOX1", "SNAP25", "GAD1", "GAD2", "PLP1", "MBP"),
              stack = TRUE,
              flip = TRUE) +ggtitle("Before percentage filter")

pdf(paste0("../results/", dir_path,  celltype, "_violin_before_percent_filter.pdf"), width = 10, height = 7)
v1
dev.off()

neuron.genes <- c("RBFOX1", "SNAP25")
contaminate.genes <- c("PLP1")
dataObject.neuron.clean$percent.neuron <- PercentageFeatureSet(dataObject.neuron.clean, features = neuron.genes)
dataObject.neuron.clean$percent.contaminate <- PercentageFeatureSet(dataObject.neuron.clean, features = contaminate.genes)

print("summary of percent neuron genes")
summary(dataObject.neuron.clean$percent.neuron)
print("summary of percent contaminate genes")
summary(dataObject.neuron.clean$percent.contaminate)

prev_nuclei <- ncol(dataObject.neuron.clean)
dataObject.neuron.clean <- subset(dataObject.neuron.clean, subset = percent.neuron > 0 & percent.contaminate == 0)
removed_count <- prev_nuclei - ncol(dataObject.neuron.clean)
print(paste0("removed nuclei based on percentage: ", removed_count))

v2 <- VlnPlot(dataObject.neuron.clean,
              features = c("RBFOX1", "SNAP25", "GAD1", "GAD2", "MBP"),
              stack = TRUE,
              flip = TRUE) +ggtitle("After percentage filter")

pdf(paste0("../results/", dir_path, celltype, "_violin_after_percent_filter.pdf"), width = 10, height = 7)
v2
dev.off()

data <- as.data.frame(table(dataObject.neuron.clean$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells2 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,10500, by = 1500), limits = c(0,10500)) +
  ggtitle(paste0(celltype, " nuclei count per sample after percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", dir_path, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 12, height = 7)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("neuron nuclei total count after filter")
sum(summary(dataObject.neuron.clean$seurat_clusters))

# clean up 
rm(dataObject.neuron, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2)

#------------------------------ endothelial
celltype <- "endothelial"

print("endothelial nuclei total count before filter")
sum(summary(dataObject.endothelial$seurat_clusters))

data <- as.data.frame(table(dataObject.endothelial$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells1 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,200, by = 50), limits = c(0,200)) +
  ggtitle(paste0(celltype, " nuclei count per sample before percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)
# remove dirty clusters 
dataObject.endothelial.clean <- subset(dataObject.endothelial, cells = WhichCells(dataObject.endothelial, idents = setdiff(unique(dataObject.endothelial$seurat_clusters), c(5, 4))))
dataObject.endothelial.dirty <- subset(dataObject.endothelial, cells = WhichCells(dataObject.endothelial, idents = c(5, 4)))

v1 <- VlnPlot(dataObject.endothelial.clean,
              features = c("FLT1", "ADGRF5", "PLP1", "MBP", "RBFOX1"),
              stack = TRUE,
              flip = TRUE) +ggtitle("Before percentage filter")

pdf(paste0("../results/", dir_path,  celltype, "_violin_before_percent_filter.pdf"), width = 7, height = 5)
v1
dev.off()

endothelial.genes <- c("FLT1", "ADGRF5")
contaminate.genes <- c("RBFOX1", "PLP1", "MBP")
dataObject.endothelial.clean$percent.endothelial <- PercentageFeatureSet(dataObject.endothelial.clean, features = endothelial.genes)
dataObject.endothelial.clean$percent.contaminate <- PercentageFeatureSet(dataObject.endothelial.clean, features = contaminate.genes)

print("summary of percent endothelial genes")
summary(dataObject.endothelial.clean$percent.endothelial)
print("summary of percent contaminate genes")
summary(dataObject.endothelial.clean$percent.contaminate)

prev_nuclei <- ncol(dataObject.endothelial.clean)
dataObject.endothelial.clean <- subset(dataObject.endothelial.clean, subset = percent.endothelial > 0 & percent.contaminate == 0)
removed_count <- prev_nuclei - ncol(dataObject.endothelial.clean)
print(paste0("removed nuclei based on percentage: ", removed_count))

v2 <- VlnPlot(dataObject.endothelial.clean,
              features = c("FLT1", "ADGRF5"),
              stack = TRUE,
              flip = TRUE) +ggtitle("After percentage filter")

pdf(paste0("../results/", dir_path, celltype, "_violin_after_percent_filter.pdf"), width = 7, height = 3.5)
v2
dev.off()

data <- as.data.frame(table(dataObject.endothelial.clean$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells2 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  scale_y_continuous(breaks = seq(0,200, by = 50), limits = c(0,200)) +
  ggtitle(paste0(celltype, " nuclei count per sample after percent filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", dir_path, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 12, height = 7)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("endothelial nuclei total count after filter")
sum(summary(dataObject.endothelial.clean$seurat_clusters))

# clean up 
rm(dataObject.endothelial, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2)


## ----merge_objects_and_save-----------------------------------------------------------------------------------------------------
# clean
DefaultAssay(dataObject.astrocyte.clean) <- "RNA"
dataObject.astrocyte.clean@assays$SCT <- NULL
dataObject.astrocyte.clean <- JoinLayers(dataObject.astrocyte.clean)

DefaultAssay(dataObject.oligodendrocyte.clean) <- "RNA"
dataObject.oligodendrocyte.clean@assays$SCT <- NULL
dataObject.oligodendrocyte.clean <- JoinLayers(dataObject.oligodendrocyte.clean)

DefaultAssay(dataObject.opc.clean) <- "RNA"
dataObject.opc.clean@assays$SCT <- NULL
dataObject.opc.clean <- JoinLayers(dataObject.opc.clean)

DefaultAssay(dataObject.endothelial.clean) <- "RNA"
dataObject.endothelial.clean@assays$SCT <- NULL
dataObject.endothelial.clean <- JoinLayers(dataObject.endothelial.clean)

DefaultAssay(dataObject.mural.clean) <- "RNA"
dataObject.mural.clean@assays$SCT <- NULL
dataObject.mural.clean <- JoinLayers(dataObject.mural.clean)

DefaultAssay(dataObject.microglia.clean) <- "RNA"
dataObject.microglia.clean@assays$SCT <- NULL
dataObject.microglia.clean <- JoinLayers(dataObject.microglia.clean)

DefaultAssay(dataObject.neuron.clean) <- "RNA"
dataObject.neuron.clean@assays$SCT <- NULL
dataObject.neuron.clean <- JoinLayers(dataObject.neuron.clean)

print("merge clean")
dataObject.clean <- merge(x = dataObject.astrocyte.clean, y = c(dataObject.oligodendrocyte.clean, dataObject.opc.clean, dataObject.endothelial.clean, dataObject.mural.clean, dataObject.microglia.clean, dataObject.neuron.clean))
saveRDS(dataObject.clean, paste0("../rObjects/",projectID,"_clean_RNA.rds"), compress = FALSE)

#-----------------------------
# dirty
DefaultAssay(dataObject.astrocyte.dirty) <- "RNA"
dataObject.astrocyte.dirty@assays$SCT <- NULL
dataObject.astrocyte.dirty <- JoinLayers(dataObject.astrocyte.dirty)

DefaultAssay(dataObject.oligodendrocyte.dirty) <- "RNA"
dataObject.oligodendrocyte.dirty@assays$SCT <- NULL
dataObject.oligodendrocyte.dirty <- JoinLayers(dataObject.oligodendrocyte.dirty)

DefaultAssay(dataObject.opc.dirty) <- "RNA"
dataObject.opc.dirty@assays$SCT <- NULL
dataObject.opc.dirty <- JoinLayers(dataObject.opc.dirty)

DefaultAssay(dataObject.endothelial.dirty) <- "RNA"
dataObject.endothelial.dirty@assays$SCT <- NULL
dataObject.endothelial.dirty <- JoinLayers(dataObject.endothelial.dirty)

DefaultAssay(dataObject.mural.dirty) <- "RNA"
dataObject.mural.dirty@assays$SCT <- NULL
dataObject.mural.dirty <- JoinLayers(dataObject.mural.dirty)

DefaultAssay(dataObject.microglia.dirty) <- "RNA"
dataObject.microglia.dirty@assays$SCT <- NULL
dataObject.microglia.dirty <- JoinLayers(dataObject.microglia.dirty)

DefaultAssay(dataObject.neuron.dirty) <- "RNA"
dataObject.neuron.dirty@assays$SCT <- NULL
dataObject.neuron.dirty <- JoinLayers(dataObject.neuron.dirty)

print("merge dirty")
dataObject.dirty <- merge(x = dataObject.astrocyte.dirty, y = c(dataObject.oligodendrocyte.dirty, dataObject.opc.dirty, dataObject.endothelial.dirty, dataObject.mural.dirty, dataObject.microglia.dirty, dataObject.neuron.dirty))
saveRDS(dataObject.dirty, paste0("../rObjects/",projectID,"_dirty_RNA.rds"), compress = FALSE)