## ----setup, include=FALSE--------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

# Libraris, paths, colors
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
project_ID <- "CWOW_cellbender_RPCAIntegration"
color_panel <- dittoColors()


## ----read_object-----------------------------------------------------------------------------------------------------
# read object
dataObject<- readRDS(paste0("../rObjects/",project_ID,"_annotated_before_recluster.rds"))

# Split by cell type 
object_list <- SplitObject(dataObject, split.by = "individual_clusters")
## cell types
# neuron
# opc
# astrocyte
# oligodendrocyte
# endothelial
# mural
# microglia  

## ----astrocyte-------------------------------------------------------------------------------------------------------
# transform
dataObject.astrocyte <- SCTransform(object_list$astrocyte, verbose = FALSE)

# run PCA on the merged object
dataObject.astrocyte <- RunPCA(object = dataObject.astrocyte)
Idents(dataObject.astrocyte) <- "Sample_ID"

# Determine the K-nearest neighbor graph
dataObject.astrocyte <- FindNeighbors(object = dataObject.astrocyte,
                                      assay = "SCT",
                                      reduction = "pca",
                                      dims = 1:30)
# Run UMAP
dataObject.astrocyte <- RunUMAP(dataObject.astrocyte,
                                dims = 1:30,
                                reduction = "pca",
                                n.components = 3)

# Determine the clusters for various resolutions
dataObject.astrocyte <- FindClusters(object = dataObject.astrocyte,
                                     algorithm = 1,
                                     resolution = 0.7)

Idents(dataObject.astrocyte) <- dataObject.astrocyte$SCT_snn_res.0.7
dataObject.astrocyte$seurat_clusters <- dataObject.astrocyte$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.astrocyte,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/recluster_pass1/UMAP_recluster_astrocyte.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()


dot_ind <- DotPlot(dataObject.astrocyte,
                   features = genes_markers,
                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
pdf(
  paste0("../results/recluster_pass1/DotPlot_recluster_astrocyte.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.astrocyte, paste0("../rObjects/",project_ID,"_astrocyte.rds"), compress = FALSE)
rm(dataObject.astrocyte)

## ----oligodendrocyte-------------------------------------------------------------------------------------------------
# transform
dataObject.oligodendrocyte <- SCTransform(object_list$oligodendrocyte, verbose = FALSE)

# run PCA on the merged object
dataObject.oligodendrocyte <- RunPCA(object = dataObject.oligodendrocyte)
Idents(dataObject.oligodendrocyte) <- "Sample_ID"

# Determine the K-nearest neighbor graph
dataObject.oligodendrocyte <- FindNeighbors(object = dataObject.oligodendrocyte,
                                            assay = "SCT",
                                            reduction = "pca",
                                            dims = 1:30)
# Run UMAP
dataObject.oligodendrocyte <- RunUMAP(dataObject.oligodendrocyte,
                                      dims = 1:30,
                                      reduction = "pca",
                                      n.components = 3)

# Determine the clusters for various resolutions
dataObject.oligodendrocyte <- FindClusters(object = dataObject.oligodendrocyte,
                                           algorithm = 1,
                                           resolution = 0.7)

Idents(dataObject.oligodendrocyte) <- dataObject.oligodendrocyte$SCT_snn_res.0.7
dataObject.oligodendrocyte$seurat_clusters <- dataObject.oligodendrocyte$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.oligodendrocyte,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/recluster_pass1/UMAP_recluster_oligodendrocyte.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()

dot_ind <- DotPlot(dataObject.oligodendrocyte,
                   features = genes_markers,

                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/recluster_pass1/DotPlot_recluster_oligodendrocyte.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.oligodendrocyte, paste0("../rObjects/",project_ID,"_oligodendrocyte.rds"), compress = FALSE)
rm(dataObject.oligodendrocyte)

## ----polydencrocyte--------------------------------------------------------------------------------------------------
# transform
dataObject.polydencrocyte <- SCTransform(object_list$polydencrocyte, verbose = FALSE)

# run PCA on the merged object
dataObject.polydencrocyte <- RunPCA(object = dataObject.polydencrocyte)
Idents(dataObject.polydencrocyte) <- "Sample_ID"

# Determine the K-nearest neighbor graph
dataObject.polydencrocyte <- FindNeighbors(object = dataObject.polydencrocyte,
                                           assay = "SCT",
                                           reduction = "pca",
                                           dims = 1:30)
# Run UMAP
dataObject.polydencrocyte <- RunUMAP(dataObject.polydencrocyte,
                                     dims = 1:30,
                                     reduction = "pca",
                                     n.components = 3)

# Determine the clusters for various resolutions
dataObject.polydencrocyte <- FindClusters(object = dataObject.polydencrocyte,
                                          algorithm = 1,
                                          resolution = 0.7)

Idents(dataObject.polydencrocyte) <- dataObject.polydencrocyte$SCT_snn_res.0.7
dataObject.polydencrocyte$seurat_clusters <- dataObject.polydencrocyte$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.polydencrocyte,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/recluster_pass1/UMAP_recluster_polydencrocyte.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()


dot_ind <- DotPlot(dataObject.polydencrocyte,
                   features = genes_markers,

                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/recluster_pass1/DotPlot_recluster_polydencrocyte.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.polydencrocyte, paste0("../rObjects/",project_ID,"_polydencrocyte.rds"), compress = FALSE)
rm(dataObject.polydencrocyte)


## ----mural-----------------------------------------------------------------------------------------------------------
# transform
dataObject.mural <- SCTransform(object_list$mural, verbose = FALSE)

# run PCA on the merged object
dataObject.mural <- RunPCA(object = dataObject.mural)
Idents(dataObject.mural) <- "Sample_ID"

# Determine the K-nearest neighbor graph
dataObject.mural <- FindNeighbors(object = dataObject.mural,
                                  assay = "SCT",
                                  reduction = "pca",
                                  dims = 1:30)
# Run UMAP
dataObject.mural <- RunUMAP(dataObject.mural,
                            dims = 1:30,
                            reduction = "pca",
                            n.components = 3)

# Determine the clusters for various resolutions
dataObject.mural <- FindClusters(object = dataObject.mural,
                                 algorithm = 1,
                                 resolution = 0.7)

Idents(dataObject.mural) <- dataObject.mural$SCT_snn_res.0.7
dataObject.mural$seurat_clusters <- dataObject.mural$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.mural,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/recluster_pass1/UMAP_recluster_mural.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()


dot_ind <- DotPlot(dataObject.mural,
                   features = genes_markers,
                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/recluster_pass1/DotPlot_recluster_mural.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.mural, paste0("../rObjects/",project_ID,"_mural.rds"), compress = FALSE)
rm(dataObject.mural)

## ----microglia-------------------------------------------------------------------------------------------------------
# transform
dataObject.microglia <- SCTransform(object_list$microglia, verbose = FALSE)

# run PCA on the merged object
dataObject.microglia <- RunPCA(object = dataObject.microglia)
Idents(dataObject.microglia) <- "Sample_ID"

# Determine the K-nearest neighbor graph
dataObject.microglia <- FindNeighbors(object = dataObject.microglia,
                                      assay = "SCT",
                                      reduction = "pca",
                                      dims = 1:30)
# Run UMAP
dataObject.microglia <- RunUMAP(dataObject.microglia,
                                dims = 1:30,
                                reduction = "pca",
                                n.components = 3)

# Determine the clusters for various resolutions
dataObject.microglia <- FindClusters(object = dataObject.microglia,
                                     algorithm = 1,
                                     resolution = 0.7)

Idents(dataObject.microglia) <- dataObject.microglia$SCT_snn_res.0.7
dataObject.microglia$seurat_clusters <- dataObject.microglia$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.microglia,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/recluster_pass1/UMAP_recluster_microglia.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()


dot_ind <- DotPlot(dataObject.microglia,
                   features = genes_markers,

                   cluster.idents = TRUE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/recluster_pass1/DotPlot_recluster_microglia.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.microglia, paste0("../rObjects/",project_ID,"_microglia.rds"))
rm(dataObject.microglia)

## ----neuron----------------------------------------------------------------------------------------------------------
# transform
dataObject.neuron <- SCTransform(object_list$neuron, verbose = FALSE)

# run PCA on the merged object
dataObject.neuron <- RunPCA(object = dataObject.neuron)
Idents(dataObject.neuron) <- "Sample_ID"

# Determine the K-nearest neighbor graph
dataObject.neuron <- FindNeighbors(object = dataObject.neuron,
                                   assay = "SCT",
                                   reduction = "pca",
                                   dims = 1:30)
# Run UMAP
dataObject.neuron <- RunUMAP(dataObject.neuron,
                             dims = 1:30,
                             reduction = "pca",
                             n.components = 3)

# Determine the clusters for various resolutions
dataObject.neuron <- FindClusters(object = dataObject.neuron,
                                  algorithm = 1,
                                  resolution = 0.7)

Idents(dataObject.neuron) <- dataObject.neuron$SCT_snn_res.0.7
dataObject.neuron$seurat_clusters <- dataObject.neuron$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.neuron,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/recluster_pass1/UMAP_recluster_neuron.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()


dot_ind <- DotPlot(dataObject.neuron,
                   features = genes_markers,
                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/recluster_pass1/DotPlot_recluster_neuron.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.neuron, paste0("../rObjects/",project_ID,"_neuron.rds"), compress = FALSE)
rm(dataObject.neuron)

## ----endothelial-----------------------------------------------------------------------------------------------------
# transform

data <- as.data.frame(table(object_list$endothelial$Sample_ID))
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
  scale_y_continuous(breaks = seq(0,200, by = 50), limits = c(0,200)) +
  ggtitle("Nuclei per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1

# LBD_ATS_F1 only has 3 cells and thus can't have SCTransform applied 
# remove this sample from the endothelial cluster
dataObject.endothelial <- subset(object_list$endothelial, subset = Sample_ID != "LBD_ATS_F1")
DefaultAssay(dataObject.endothelial) <- "RNA"
dataObject.endothelial <- SCTransform(dataObject.endothelial, verbose = TRUE)

# run PCA on the merged object
dataObject.endothelial <- RunPCA(object = dataObject.endothelial)
Idents(dataObject.endothelial) <- "Sample_ID"

# Determine the K-nearest neighbor graph
dataObject.endothelial <- FindNeighbors(object = dataObject.endothelial, 
                                        assay = "SCT", 
                                        reduction = "pca",
                                        dims = 1:15)
# Run UMAP
dataObject.endothelial <- RunUMAP(dataObject.endothelial,
                                  dims = 1:15,
                                  reduction = "pca",
                                  n.components = 3) 

# Determine the clusters for various resolutions
dataObject.endothelial <- FindClusters(object = dataObject.endothelial,
                                       algorithm = 1, 
                                       resolution = 0.7)

Idents(dataObject.endothelial) <- dataObject.endothelial$SCT_snn_res.0.7
dataObject.endothelial$seurat_clusters <- dataObject.endothelial$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.endothelial,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/recluster_pass1/UMAP_recluster_endothelial.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()

dot_ind <- DotPlot(dataObject.endothelial,
                   features = genes_markers, 
                   cluster.idents = TRUE,
                   dot.scale = 8) + RotatedAxis()
pdf(
  paste0("../results/recluster_pass1/DotPlot_recluster_endothelial.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.endothelial, paste0("../rObjects/",project_ID,"_endothelial.rds"), compress = FALSE)
rm(dataObject.endothelial)
