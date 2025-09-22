## ----setup, include=FALSE--------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

# Libraris, paths, colors
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
project_ID <- "CWOW_cellbender_RPCAIntegration"
path_dir <- "recluster_pass1/"

## ----read_object-----------------------------------------------------------------------------------------------------
# read object
dataObject.astrocyte <- readRDS(paste0("../rObjects/", project_ID, "_astrocyte.rds"))
dataObject.oligodendrocyte <- readRDS(paste0("../rObjects/",project_ID,"_oligodendrocyte.rds"))
dataObject.opc <- readRDS(paste0("../rObjects/",project_ID,"_polydencrocyte.rds"))
dataObject.endothelial <- readRDS(paste0("../rObjects/",project_ID,"_endothelial.rds"))
dataObject.mural<- readRDS(paste0("../rObjects/",project_ID,"_mural.rds"))
dataObject.microglia <- readRDS(paste0("../rObjects/",project_ID,"_microglia.rds"))
dataObject.neuron <- readRDS(paste0("../rObjects/",project_ID,"_neuron.rds"))

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
  scale_y_continuous(breaks = seq(0,4000, by = 500), limits = c(0,4000)) +
  ggtitle(paste0(celltype, " nuclei count per sample before extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1

v1 <- VlnPlot(dataObject.astrocyte,
              features = c("AQP4", "GJA1", "PLP1", "MBP", "RBFOX1", "SNAP25"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("Before extra filtering")

pdf(paste0("../results/", path_dir,  celltype, "_violin_before_extra_filtering.pdf"), width = 7, height = 5)
v1
dev.off()

# remove dirty clusters 
dataObject.astrocyte.clean <- subset(dataObject.astrocyte, cells = WhichCells(dataObject.astrocyte, idents = setdiff(unique(dataObject.astrocyte$seurat_clusters), c(23,13,12,11,10,7,4))))
dataObject.astrocyte.dirty <- subset(dataObject.astrocyte, cells = WhichCells(dataObject.astrocyte, idents = c(23,13,12,11,10,7,4)))

# percentage filter
astrocyte.genes <- c("AQP4", "GJA1")
contaminate.genes <- c("SNAP25", "PLP1")
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

# Keep cells where expression is greater than 0.1
prev_nuclei <- ncol(dataObject.astrocyte.clean)
cells_to_keep <- WhichCells(dataObject.astrocyte.clean, expression = AQP4 > 0.5)
dataObject.astrocyte.clean <- subset(dataObject.astrocyte.clean, cells = cells_to_keep)
removed_count <- prev_nuclei - ncol(dataObject.astrocyte.clean)
print(paste0("removed nuclei based on min expression: ", removed_count))

v2 <- VlnPlot(dataObject.astrocyte.clean,
              features = c("AQP4", "GJA1", "PLP1", "MBP", "RBFOX1", "SNAP25"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("After extra filtering")
pdf(paste0("../results/", path_dir, celltype, "_violin_after_percent_filter.pdf"), width = 7, height = 5)
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
  ggtitle(paste0(celltype, " nuclei count per sample after extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", path_dir, celltype, "_cells_per_sample_before_and_after_extra_filter.pdf"), width = 10, height = 6)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("astrocyte nuclei total count after filter")
sum(summary(dataObject.astrocyte.clean$seurat_clusters))

# clean up 
rm(dataObject.astrocyte, astrocyte.genes, cells_to_keep, celltype, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2, contaminate.genes, prev_nuclei, removed_count)

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
  scale_y_continuous(breaks = seq(0,15000, by = 1500), limits = c(0,15000)) +
  ggtitle(paste0(celltype, " nuclei count per sample before extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)

v1 <- VlnPlot(dataObject.oligodendrocyte,
              features = c("PLP1", "MBP", "RBFOX1", "SNAP25"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("Before extra filtering")

pdf(paste0("../results/", path_dir,  celltype, "_violin_before_extra_filtering.pdf"), width = 7, height = 5)
v1
dev.off()

# remove dirty clusters 
dataObject.oligodendrocyte.clean <- subset(dataObject.oligodendrocyte, cells = WhichCells(dataObject.oligodendrocyte, idents = setdiff(unique(dataObject.oligodendrocyte$seurat_clusters), c(21))))
dataObject.oligodendrocyte.dirty <- subset(dataObject.oligodendrocyte, cells = WhichCells(dataObject.oligodendrocyte, idents = c(21)))

# percentage filter
oligodendrocyte.genes <- c("PLP1", "MBP")
contaminate.genes <- c("SNAP25")
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

# Keep cells where expression is greater than 0.1
prev_nuclei <- ncol(dataObject.oligodendrocyte.clean)
cells_to_keep <- WhichCells(dataObject.oligodendrocyte.clean, expression = PLP1 > 0.5)
dataObject.oligodendrocyte.clean <- subset(dataObject.oligodendrocyte.clean, cells = cells_to_keep)
removed_count <- prev_nuclei - ncol(dataObject.oligodendrocyte.clean)
print(paste0("removed nuclei based on min expression: ", removed_count))

v1 <- VlnPlot(dataObject.oligodendrocyte.clean,
              features = c("PLP1", "MBP", "RBFOX1", "SNAP25"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("After extra filtering")

pdf(paste0("../results/", path_dir,  celltype, "_violin_after_extra_filtering.pdf"), width = 7, height = 5)
v1
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
  scale_y_continuous(breaks = seq(0,15000, by = 1500), limits = c(0,15000)) +
  ggtitle(paste0(celltype, " nuclei count per sample after extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", path_dir, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 10, height = 6)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("oligodendrocyte nuclei total count after filter")
sum(summary(dataObject.oligodendrocyte.clean$seurat_clusters))

# clean up 
rm(dataObject.oligodendrocyte, oligodendrocyte.genes, cells_to_keep, celltype, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2, contaminate.genes, prev_nuclei, removed_count)

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
  scale_y_continuous(breaks = seq(0,2500, by = 500), limits = c(0,2500)) +
  ggtitle(paste0(celltype, " nuclei count per sample before extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)

v1 <- VlnPlot(dataObject.microglia,
              features = c("P2RY12", "DOCK8", "PLP1", "MBP", "RBFOX1", "SNAP25"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("Before extra filtering")
pdf(paste0("../results/", path_dir, celltype, "_violin_before_extra_filter.pdf"), width = 7, height = 5)
v1
dev.off()

# remove dirty clusters 
dataObject.microglia.clean <- subset(dataObject.microglia, cells = WhichCells(dataObject.microglia, idents = setdiff(unique(dataObject.microglia$seurat_clusters), c(4,3,12,15))))
dataObject.microglia.dirty <- subset(dataObject.microglia, cells = WhichCells(dataObject.microglia, idents = c(4,3,12,15)))

microglia.genes <- c("DOCK8", "P2RY12")
contaminate.genes <- c("SNAP25", "PLP1")
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

prev_nuclei <- ncol(dataObject.microglia.clean)
cells_to_keep <- WhichCells(dataObject.microglia.clean, expression = DOCK8 > 0.5)
dataObject.microglia.clean <- subset(dataObject.microglia.clean, cells = cells_to_keep)
removed_count <- prev_nuclei - ncol(dataObject.microglia.clean)
print(paste0("removed nuclei based on min expression: ", removed_count))

v2 <- VlnPlot(dataObject.microglia.clean,
              features = c("P2RY12", "DOCK8", "PLP1", "MBP", "RBFOX1", "SNAP25"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("After extra filtering")
pdf(paste0("../results/", path_dir, celltype, "_violin_after_extra_filter.pdf"), width = 7, height = 5)
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
  scale_y_continuous(breaks = seq(0,2500, by = 500), limits = c(0,2500)) +
  ggtitle(paste0(celltype, " nuclei count per sample after extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", path_dir, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 10, height = 6)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("microglia nuclei total count after filter")
sum(summary(dataObject.microglia.clean$seurat_clusters))

# clean up 
rm(dataObject.microglia, microglia.genes, cells_to_keep, celltype, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2, contaminate.genes, prev_nuclei, removed_count)

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
  scale_y_continuous(breaks = seq(0,1500, by = 100), limits = c(0,1500)) +
  ggtitle(paste0(celltype, " nuclei count per sample before extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)

v1 <- VlnPlot(dataObject.mural,
              features = c("DCN", "COL1A2", "PLP1", "MBP", "RBFOX1", "SNAP25", "EBF1", "RBPMS"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("Before extra filtering")
pdf(paste0("../results/", path_dir, celltype, "_violin_before_extra_filter.pdf"), width = 7, height = 5)
v1
dev.off()

# remove dirty clusters 
dataObject.mural.clean <- subset(dataObject.mural, cells = WhichCells(dataObject.mural, idents = setdiff(unique(dataObject.mural$seurat_clusters), c(3,6,7,8,11,12))))
dataObject.mural.dirty <- subset(dataObject.mural, cells = WhichCells(dataObject.mural, idents = c(3,6,7,8,11,12)))


mural.genes <- c("EBF1", "RBPMS")
contaminate.genes <- c("SNAP25", "PLP1")
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

# Example: Keep cells where expression is greater than 0.1
prev_nuclei <- ncol(dataObject.mural.clean)
cells_to_keep <- WhichCells(dataObject.mural.clean, expression = EBF1 > 0.5)
dataObject.mural.clean <- subset(dataObject.mural.clean, cells = cells_to_keep)
removed_count <- prev_nuclei - ncol(dataObject.mural.clean)
print(paste0("removed nuclei based on min expression: ", removed_count))


v2 <- VlnPlot(dataObject.mural.clean,
              features = c("DCN", "COL1A2", "PLP1", "MBP", "RBFOX1", "SNAP25", "EBF1", "RBPMS"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("After extra filter")

pdf(paste0("../results/", path_dir, celltype, "_violin_after_percent_filter.pdf"), width = 7, height = 5)
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
  scale_y_continuous(breaks = seq(0,1500, by = 100), limits = c(0,1500)) +
  ggtitle(paste0(celltype, " nuclei count per sample after extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", path_dir, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 10, height = 6)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("mural nuclei total count after filter")
sum(summary(dataObject.mural.clean$seurat_clusters))

# clean up 
rm(dataObject.mural, mural.genes, cells_to_keep, celltype, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2, contaminate.genes, prev_nuclei, removed_count)

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
  ggtitle(paste0(celltype, " nuclei count per sample before extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)

v1 <- VlnPlot(dataObject.opc,
              features = c("OLIG1", "GJA1", "PLP1", "MBP", "RBFOX1", "SNAP25", "VCAN"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("Before extra filtering")

pdf(paste0("../results/", path_dir,  celltype, "_violin_before_extra_filtering.pdf"), width = 7, height = 5)
v1
dev.off()

# remove dirty clusters 
dataObject.opc.clean <- subset(dataObject.opc, cells = WhichCells(dataObject.opc, idents = setdiff(unique(dataObject.opc$seurat_clusters), c(20,5,11))))
dataObject.opc.dirty <- subset(dataObject.opc, cells = WhichCells(dataObject.opc, idents = c(20,5,11)))

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


prev_nuclei <- ncol(dataObject.opc.clean)
cells_to_keep <- WhichCells(dataObject.opc.clean, expression = OLIG1 > 0.5)
dataObject.opc.clean <- subset(dataObject.opc.clean, cells = cells_to_keep)
removed_count <- prev_nuclei - ncol(dataObject.opc.clean)
print(paste0("removed nuclei based on min expression: ", removed_count))

v2 <- VlnPlot(dataObject.opc.clean,
              features = c("OLIG1", "GJA1", "PLP1", "MBP", "RBFOX1", "SNAP25", "VCAN"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("After extra filter")

pdf(paste0("../results/", path_dir, celltype, "_violin_after_percent_filter.pdf"), width = 7, height = 5)
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
  ggtitle(paste0(celltype, " nuclei count per sample after extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", path_dir, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 10, height = 6)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("opc nuclei total count after filter")
sum(summary(dataObject.opc.clean$seurat_clusters))

# clean up 
rm(dataObject.opc, opc.genes, cells_to_keep, celltype, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2, contaminate.genes, prev_nuclei, removed_count)

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
  ggtitle(paste0(celltype, " nuclei count per sample before extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)

v1 <- VlnPlot(dataObject.neuron,
              features = c("AQP4", "PLP1","RBFOX1", "SNAP25", "GAD1", "GAD2"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("Before extra filtering")

pdf(paste0("../results/", path_dir,  celltype, "_violin_before_extra_filtering.pdf"), width = 7, height = 5)
v1
dev.off()

# remove dirty clusters 
dataObject.neuron.clean <- subset(dataObject.neuron, cells = WhichCells(dataObject.neuron, idents = setdiff(unique(dataObject.neuron$seurat_clusters), c(9,18,32,33,38))))
dataObject.neuron.dirty <- subset(dataObject.neuron, cells = WhichCells(dataObject.neuron, idents = c(9,18,32,33,38)))

neuron.genes <- c("SNAP25")
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

prev_nuclei <- ncol(dataObject.neuron.clean)
cells_to_keep <- WhichCells(dataObject.neuron.clean, expression = SNAP25 > 0.5)
dataObject.neuron.clean <- subset(dataObject.neuron.clean, cells = cells_to_keep)
removed_count <- prev_nuclei - ncol(dataObject.neuron.clean)
print(paste0("removed nuclei based on min expression: ", removed_count))

v2 <- VlnPlot(dataObject.neuron.clean,
              features = c("AQP4", "PLP1","RBFOX1", "SNAP25", "GAD1", "GAD2"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("After extra filter")

pdf(paste0("../results/", path_dir, celltype, "_violin_after_percent_filter.pdf"), width = 10, height = 7)
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
  ggtitle(paste0(celltype, " nuclei count per sample after extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", path_dir, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 10, height = 6)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("neuron nuclei total count after filter")
sum(summary(dataObject.neuron.clean$seurat_clusters))

# clean up 
rm(dataObject.neuron, neuron.genes, cells_to_keep, celltype, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2, contaminate.genes, prev_nuclei, removed_count)

#------------------------------ endothelial
celltype <- "endothelial"

print("endothelial nuclei total count before filter")
sum(summary(dataObject.endothelial$seurat_clusters))

dataObject.endothelial <- subset(dataObject.endothelial, subset = Sample_ID != "AD_AT_M4")
dataObject.endothelial <- subset(dataObject.endothelial, subset = Sample_ID != "Ctr_M4")
dataObject.endothelial <- subset(dataObject.endothelial, subset = Sample_ID != "LBD_ATS_F1")

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
  ggtitle(paste0(celltype, " nuclei count per sample before extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
rm(data)

v1 <- VlnPlot(dataObject.endothelial,
              features = c("AQP4", "FLT1", "ADGRF5", "PLP1", "MBP", "RBFOX1", "SNAP25"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("Before extra filtering")

pdf(paste0("../results/", path_dir,  celltype, "_violin_before_extra_filtering.pdf"), width = 7, height = 5)
v1
dev.off()

# remove dirty clusters 
dataObject.endothelial.clean <- subset(dataObject.endothelial, cells = WhichCells(dataObject.endothelial, idents = setdiff(unique(dataObject.endothelial$seurat_clusters), c(3,4))))
dataObject.endothelial.dirty <- subset(dataObject.endothelial, cells = WhichCells(dataObject.endothelial, idents = c(3,4)))


endothelial.genes <- c("FLT1")
contaminate.genes <- c("SNAP25", "PLP1")
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

# Example: Keep cells where expression is greater than 0.1
prev_nuclei <- ncol(dataObject.endothelial.clean)
cells_to_keep <- WhichCells(dataObject.endothelial.clean, expression = FLT1 > 0.5)
dataObject.endothelial.clean <- subset(dataObject.endothelial.clean, cells = cells_to_keep)
removed_count <- prev_nuclei - ncol(dataObject.endothelial.clean)
print(paste0("removed nuclei based on min expression: ", removed_count))


v2 <- VlnPlot(dataObject.endothelial.clean,
              features = c("AQP4", "FLT1", "ADGRF5", "PLP1", "MBP", "RBFOX1", "SNAP25"),
              stack = TRUE,
              same.y.lims = TRUE,
              flip = TRUE) +ggtitle("After extra filter")

pdf(paste0("../results/", path_dir, celltype, "_violin_after_percent_filter.pdf"), width = 7, height = 5)
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
  ggtitle(paste0(celltype, " nuclei count per sample after extra filter")) +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

pdf(paste0("../results/", path_dir, celltype, "_cells_per_sample_before_and_after_percent_filter.pdf"), width = 10, height = 6)
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
dev.off()

print("endothelial nuclei total count after filter")
sum(summary(dataObject.endothelial.clean$seurat_clusters))

# clean up 
rm(dataObject.endothelial, endothelial.genes, cells_to_keep, celltype, data, ncells1, ncells2, v1, v2, grid2, layout2, plots2, contaminate.genes, prev_nuclei, removed_count)

## ----merge_objects_and_save-----------------------------------------------------------------------------------------------------
# print("merge clean start")
# # clean
# DefaultAssay(dataObject.astrocyte.clean) <- "RNA"
# dataObject.astrocyte.clean@assays$SCT <- NULL
# dataObject.astrocyte.clean <- JoinLayers(dataObject.astrocyte.clean)
# 
# DefaultAssay(dataObject.oligodendrocyte.clean) <- "RNA"
# dataObject.oligodendrocyte.clean@assays$SCT <- NULL
# dataObject.oligodendrocyte.clean <- JoinLayers(dataObject.oligodendrocyte.clean)
# 
# DefaultAssay(dataObject.opc.clean) <- "RNA"
# dataObject.opc.clean@assays$SCT <- NULL
# dataObject.opc.clean <- JoinLayers(dataObject.opc.clean)
# 
# # remove samples with less than 3 cells, else join layers won't work 
# DefaultAssay(dataObject.endothelial.clean) <- "RNA"
# dataObject.endothelial.clean@assays$SCT <- NULL
# dataObject.endothelial.clean <- JoinLayers(dataObject.endothelial.clean)
# 
# DefaultAssay(dataObject.mural.clean) <- "RNA"
# dataObject.mural.clean@assays$SCT <- NULL
# dataObject.mural.clean <- JoinLayers(dataObject.mural.clean)
# 
# DefaultAssay(dataObject.microglia.clean) <- "RNA"
# dataObject.microglia.clean@assays$SCT <- NULL
# dataObject.microglia.clean <- JoinLayers(dataObject.microglia.clean)
# 
# DefaultAssay(dataObject.neuron.clean) <- "RNA"
# dataObject.neuron.clean@assays$SCT <- NULL
# dataObject.neuron.clean <- JoinLayers(dataObject.neuron.clean)
# 
# dataObject.clean <- merge(x = dataObject.astrocyte.clean, y = c(dataObject.oligodendrocyte.clean, dataObject.opc.clean, dataObject.endothelial.clean, dataObject.mural.clean, dataObject.microglia.clean, dataObject.neuron.clean))
# saveRDS(dataObject.clean, paste0("../rObjects/",project_ID,"_clean_RNA.rds"), compress = FALSE)

#------- Dirty 
print("merge dirty start")
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

print("merge dirty save")
dataObject.dirty <- merge(x = dataObject.astrocyte.dirty, y = c(dataObject.oligodendrocyte.dirty, dataObject.opc.dirty, dataObject.endothelial.dirty, dataObject.mural.dirty, dataObject.microglia.dirty, dataObject.neuron.dirty))
saveRDS(dataObject.dirty, paste0("../rObjects/",project_ID,"_dirty_RNA.rds"), compress = FALSE)