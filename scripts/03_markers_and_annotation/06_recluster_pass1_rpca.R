## ----setup, include=FALSE--------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

# Libraris, paths, colors
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender_RPCAIntegration"
color.panel <- dittoColors()


## ----read_object-----------------------------------------------------------------------------------------------------
# read object
dataObject<- readRDS(paste0("../rObjects/",projectID,"_annotated_before_recluster.rds"))
# inspect
#table(dataObject$individual_clusters)
#DimPlot(dataObject, reduction = "integrated.rpca", group.by = "individual_clusters") 
#dataObject$individual_clusters
# Split by cell type 
object_list <- SplitObject(dataObject, split.by = "individual_clusters")
# cell types
# neuron
# opc
# astrocyte
# oligodendrocyte
# endothelial
# mural
# microglia  

# ## ----astrocyte-------------------------------------------------------------------------------------------------------
# # transform
# dataObject.astrocyte <- SCTransform(object_list$astrocyte, verbose = FALSE)
# 
# # run PCA on the merged object
# dataObject.astrocyte <- RunPCA(object = dataObject.astrocyte)
# Idents(dataObject.astrocyte) <- "SampleID"
# 
# # Determine the K-nearest neighbor graph
# dataObject.astrocyte <- FindNeighbors(object = dataObject.astrocyte,
#                                       assay = "SCT",
#                                       reduction = "pca",
#                                       dims = 1:30)
# # Run UMAP
# dataObject.astrocyte <- RunUMAP(dataObject.astrocyte,
#                                 dims = 1:30,
#                                 reduction = "pca",
#                                 n.components = 3)
# 
# # Determine the clusters for various resolutions
# dataObject.astrocyte <- FindClusters(object = dataObject.astrocyte,
#                                      algorithm = 1, 
#                                      resolution = 0.7)
# 
# Idents(dataObject.astrocyte) <- dataObject.astrocyte$SCT_snn_res.0.7
# dataObject.astrocyte$seurat_clusters <- dataObject.astrocyte$SCT_snn_res.0.7
# ditto_umap <- dittoDimPlot(object = dataObject.astrocyte,
#                            var = "seurat_clusters",
#                            reduction.use = "umap",
#                            do.label = TRUE,
#                            labels.highlight = TRUE)
# ditto_umap
# pdf(
#   paste0("../results/UMAP/reclusters/",
#          projectID, "_recluster_astrocyte.pdf"
#   ),
#   width = 7,
#   height = 6
# )
# ditto_umap
# dev.off()
# 
# 
# dot_ind <- DotPlot(dataObject.astrocyte,
#                    features = markers.to.plot,
#                    cluster.idents = FALSE,
#                    dot.scale = 8) + RotatedAxis()
# pdf(
#   paste0("../results/dot_plot/reclusters/",
#          projectID, "_recluster_astrocyte.pdf"
#   ),
#   width = 14,
#   height = 7
# )
# dot_ind
# dev.off()
# 
# saveRDS(dataObject.astrocyte, paste0("../rObjects/",projectID,"_astrocyte.rds"), compress = FALSE)
# rm(dataObject.astrocyte)
# 
# ## ----oligodendrocyte-------------------------------------------------------------------------------------------------
# # transform
# dataObject.oligodendrocyte <- SCTransform(object_list$oligodendrocyte, verbose = FALSE)
# 
# # run PCA on the merged object
# dataObject.oligodendrocyte <- RunPCA(object = dataObject.oligodendrocyte)
# Idents(dataObject.oligodendrocyte) <- "SampleID"
# 
# # Determine the K-nearest neighbor graph
# dataObject.oligodendrocyte <- FindNeighbors(object = dataObject.oligodendrocyte,
#                                             assay = "SCT",
#                                             reduction = "pca",
#                                             dims = 1:30)
# # Run UMAP
# dataObject.oligodendrocyte <- RunUMAP(dataObject.oligodendrocyte,
#                                       dims = 1:30,
#                                       reduction = "pca",
#                                       n.components = 3)
# 
# # Determine the clusters for various resolutions
# dataObject.oligodendrocyte <- FindClusters(object = dataObject.oligodendrocyte,
#                                            algorithm = 1, 
#                                            resolution = 0.7)
# 
# Idents(dataObject.oligodendrocyte) <- dataObject.oligodendrocyte$SCT_snn_res.0.7
# dataObject.oligodendrocyte$seurat_clusters <- dataObject.oligodendrocyte$SCT_snn_res.0.7
# ditto_umap <- dittoDimPlot(object = dataObject.oligodendrocyte,
#                            var = "seurat_clusters",
#                            reduction.use = "umap",
#                            do.label = TRUE,
#                            labels.highlight = TRUE)
# ditto_umap
# pdf(
#   paste0("../results/UMAP/reclusters/",
#          projectID, "_recluster_oligodendrocyte.pdf"
#   ),
#   width = 7,
#   height = 6
# )
# ditto_umap
# dev.off()
# 
# dot_ind <- DotPlot(dataObject.oligodendrocyte,
#                    features = markers.to.plot,
# 
#                    cluster.idents = FALSE,
#                    dot.scale = 8) + RotatedAxis()
# dot_ind
# pdf(
#   paste0("../results/dot_plot/reclusters/",
#          projectID, "_recluster_oligodendrocyte.pdf"
#   ),
#   width = 14,
#   height = 7
# )
# dot_ind
# dev.off()
# 
# saveRDS(dataObject.oligodendrocyte, paste0("../rObjects/",projectID,"_oligodendrocyte.rds"), compress = FALSE)
# rm(dataObject.oligodendrocyte)
# 
# ## ----polydencrocyte--------------------------------------------------------------------------------------------------
# # transform
# dataObject.polydencrocyte <- SCTransform(object_list$polydencrocyte, verbose = FALSE)
# 
# # run PCA on the merged object
# dataObject.polydencrocyte <- RunPCA(object = dataObject.polydencrocyte)
# Idents(dataObject.polydencrocyte) <- "SampleID"
# 
# # Determine the K-nearest neighbor graph
# dataObject.polydencrocyte <- FindNeighbors(object = dataObject.polydencrocyte,
#                                            assay = "SCT",
#                                            reduction = "pca",
#                                            dims = 1:30)
# # Run UMAP
# dataObject.polydencrocyte <- RunUMAP(dataObject.polydencrocyte,
#                                      dims = 1:30,
#                                      reduction = "pca",
#                                      n.components = 3)
# 
# # Determine the clusters for various resolutions
# dataObject.polydencrocyte <- FindClusters(object = dataObject.polydencrocyte,
#                                           algorithm = 1, 
#                                           resolution = 0.7)
# 
# Idents(dataObject.polydencrocyte) <- dataObject.polydencrocyte$SCT_snn_res.0.7
# dataObject.polydencrocyte$seurat_clusters <- dataObject.polydencrocyte$SCT_snn_res.0.7
# ditto_umap <- dittoDimPlot(object = dataObject.polydencrocyte,
#                            var = "seurat_clusters",
#                            reduction.use = "umap",
#                            do.label = TRUE,
#                            labels.highlight = TRUE)
# ditto_umap
# pdf(
#   paste0("../results/UMAP/reclusters/",
#          projectID, "_recluster_polydencrocyte.pdf"
#   ),
#   width = 7,
#   height = 6
# )
# ditto_umap
# dev.off()
# 
# 
# dot_ind <- DotPlot(dataObject.polydencrocyte,
#                    features = markers.to.plot,
# 
#                    cluster.idents = FALSE,
#                    dot.scale = 8) + RotatedAxis()
# dot_ind
# pdf(
#   paste0("../results/dot_plot/reclusters/",
#          projectID, "_recluster_polydencrocyte.pdf"
#   ),
#   width = 14,
#   height = 7
# )
# dot_ind
# dev.off()
# 
# saveRDS(dataObject.polydencrocyte, paste0("../rObjects/",projectID,"_polydencrocyte.rds"), compress = FALSE)
# rm(dataObject.polydencrocyte)
# 
# 
# ## ----mural-----------------------------------------------------------------------------------------------------------
# # transform
# dataObject.mural <- SCTransform(object_list$mural, verbose = FALSE)
# 
# # run PCA on the merged object
# dataObject.mural <- RunPCA(object = dataObject.mural)
# Idents(dataObject.mural) <- "SampleID"
# 
# # Determine the K-nearest neighbor graph
# dataObject.mural <- FindNeighbors(object = dataObject.mural, 
#                                   assay = "SCT", 
#                                   reduction = "pca",
#                                   dims = 1:30)
# # Run UMAP
# dataObject.mural <- RunUMAP(dataObject.mural,
#                             dims = 1:30,
#                             reduction = "pca",
#                             n.components = 3) 
# 
# # Determine the clusters for various resolutions
# dataObject.mural <- FindClusters(object = dataObject.mural,
#                                  algorithm = 1, 
#                                  resolution = 0.7)
# 
# Idents(dataObject.mural) <- dataObject.mural$SCT_snn_res.0.7
# dataObject.mural$seurat_clusters <- dataObject.mural$SCT_snn_res.0.7
# ditto_umap <- dittoDimPlot(object = dataObject.mural,
#                            var = "seurat_clusters",
#                            reduction.use = "umap",
#                            do.label = TRUE,
#                            labels.highlight = TRUE)
# ditto_umap
# pdf(
#   paste0("../results/UMAP/reclusters/",
#          projectID, "_recluster_mural.pdf"
#   ),
#   width = 7,
#   height = 6
# )
# ditto_umap
# dev.off()
# 
# 
# dot_ind <- DotPlot(dataObject.mural,
#                    features = markers.to.plot, 
#                    cluster.idents = FALSE,
#                    dot.scale = 8) + RotatedAxis()
# dot_ind
# pdf(
#   paste0("../results/dot_plot/reclusters/",
#          projectID, "_recluster_mural.pdf"
#   ),
#   width = 14,
#   height = 7
# )
# dot_ind
# dev.off()
# 
# saveRDS(dataObject.mural, paste0("../rObjects/",projectID,"_mural.rds"), compress = FALSE)
# rm(dataObject.mural)
# 
# ## ----microglia-------------------------------------------------------------------------------------------------------
# # transform
# dataObject.microglia <- SCTransform(object_list$microglia, verbose = FALSE)
# 
# # run PCA on the merged object
# dataObject.microglia <- RunPCA(object = dataObject.microglia)
# Idents(dataObject.microglia) <- "SampleID"
# 
# # Determine the K-nearest neighbor graph
# dataObject.microglia <- FindNeighbors(object = dataObject.microglia, 
#                                       assay = "SCT", 
#                                       reduction = "pca",
#                                       dims = 1:30)
# # Run UMAP
# dataObject.microglia <- RunUMAP(dataObject.microglia,
#                                 dims = 1:30,
#                                 reduction = "pca",
#                                 n.components = 3) 
# 
# # Determine the clusters for various resolutions
# dataObject.microglia <- FindClusters(object = dataObject.microglia,
#                                      algorithm = 1, 
#                                      resolution = 0.7)
# 
# Idents(dataObject.microglia) <- dataObject.microglia$SCT_snn_res.0.7
# dataObject.microglia$seurat_clusters <- dataObject.microglia$SCT_snn_res.0.7
# ditto_umap <- dittoDimPlot(object = dataObject.microglia,
#                            var = "seurat_clusters",
#                            reduction.use = "umap",
#                            do.label = TRUE,
#                            labels.highlight = TRUE)
# ditto_umap
# pdf(
#   paste0("../results/UMAP/reclusters/",
#          projectID, "_recluster_microglia.pdf"
#   ),
#   width = 7,
#   height = 6
# )
# ditto_umap
# dev.off()
# 
# 
# dot_ind <- DotPlot(dataObject.microglia,
#                    features = markers.to.plot, 
#                    
#                    cluster.idents = TRUE,
#                    dot.scale = 8) + RotatedAxis()
# dot_ind
# pdf(
#   paste0("../results/dot_plot/reclusters/",
#          projectID, "_recluster_microglia.pdf"
#   ),
#   width = 14,
#   height = 7
# )
# dot_ind
# dev.off()
# 
# saveRDS(dataObject.microglia, paste0("../rObjects/",projectID,"_microglia.rds"))
# rm(dataObject.microglia)
# 
# ## ----neuron----------------------------------------------------------------------------------------------------------
# # transform
# dataObject.neuron <- SCTransform(object_list$neuron, verbose = FALSE)
# 
# # run PCA on the merged object
# dataObject.neuron <- RunPCA(object = dataObject.neuron)
# Idents(dataObject.neuron) <- "SampleID"
# 
# # Determine the K-nearest neighbor graph
# dataObject.neuron <- FindNeighbors(object = dataObject.neuron, 
#                                    assay = "SCT", 
#                                    reduction = "pca",
#                                    dims = 1:30)
# # Run UMAP
# dataObject.neuron <- RunUMAP(dataObject.neuron,
#                              dims = 1:30,
#                              reduction = "pca",
#                              n.components = 3) 
# 
# # Determine the clusters for various resolutions
# dataObject.neuron <- FindClusters(object = dataObject.neuron,
#                                   algorithm = 1, 
#                                   resolution = 0.7)
# 
# Idents(dataObject.neuron) <- dataObject.neuron$SCT_snn_res.0.7
# dataObject.neuron$seurat_clusters <- dataObject.neuron$SCT_snn_res.0.7
# ditto_umap <- dittoDimPlot(object = dataObject.neuron,
#                            var = "seurat_clusters",
#                            reduction.use = "umap",
#                            do.label = TRUE,
#                            labels.highlight = TRUE)
# ditto_umap
# pdf(
#   paste0("../results/UMAP/reclusters/",
#          projectID, "_recluster_neuron.pdf"
#   ),
#   width = 7,
#   height = 6
# )
# ditto_umap
# dev.off()
# 
# 
# dot_ind <- DotPlot(dataObject.neuron,
#                    features = markers.to.plot, 
#                    cluster.idents = FALSE,
#                    dot.scale = 8) + RotatedAxis()
# dot_ind
# pdf(
#   paste0("../results/dot_plot/reclusters/",
#          projectID, "_recluster_neuron.pdf"
#   ),
#   width = 14,
#   height = 7
# )
# dot_ind
# dev.off()
# 
# saveRDS(dataObject.neuron, paste0("../rObjects/",projectID,"_neuron.rds"), compress = FALSE)
# rm(dataObject.neuron)
# 
# ## ----noise----------------------------------------------------------------------------------------------------------
# # transform
# dataObject.noise <- SCTransform(object_list$noise, verbose = FALSE)
# 
# # run PCA on the merged object
# dataObject.noise <- RunPCA(object = dataObject.noise)
# Idents(dataObject.noise) <- "SampleID"
# 
# # Determine the K-nearest neighbor graph
# dataObject.noise <- FindNeighbors(object = dataObject.noise, 
#                                   assay = "SCT", 
#                                   reduction = "pca",
#                                   dims = 1:30)
# # Run UMAP
# dataObject.noise <- RunUMAP(dataObject.noise,
#                             dims = 1:30,
#                             reduction = "pca",
#                             n.components = 3) 
# 
# # Determine the clusters for various resolutions
# dataObject.noise <- FindClusters(object = dataObject.noise,
#                                  algorithm = 1, 
#                                  resolution = 0.7)
# 
# Idents(dataObject.noise) <- dataObject.noise$SCT_snn_res.0.7
# dataObject.noise$seurat_clusters <- dataObject.noise$SCT_snn_res.0.7
# ditto_umap <- dittoDimPlot(object = dataObject.noise,
#                            var = "seurat_clusters",
#                            reduction.use = "umap",
#                            do.label = TRUE,
#                            labels.highlight = TRUE)
# ditto_umap
# pdf(
#   paste0("../results/UMAP/reclusters/",
#          projectID, "_recluster_noise.pdf"
#   ),
#   width = 7,
#   height = 6
# )
# ditto_umap
# dev.off()
# 
# 
# dot_ind <- DotPlot(dataObject.noise,
#                    features = markers.to.plot, 
#                    cluster.idents = FALSE,
#                    dot.scale = 8) + RotatedAxis()
# dot_ind
# pdf(
#   paste0("../results/dot_plot/reclusters/",
#          projectID, "_recluster_noise.pdf"
#   ),
#   width = 14,
#   height = 7
# )
# dot_ind
# dev.off()
# 
# saveRDS(dataObject.noise, paste0("../rObjects/",projectID,"_noise.rds"), compress = FALSE)
# rm(dataObject.noise)

## ----endothelial-----------------------------------------------------------------------------------------------------
# transform
print("endothelial")
# Assuming 'orig.ident' is the metadata column containing sample IDs
object_list$endothelial_subset <- subset(object_list$endothelial, subset = Sample_ID != "LBD_ATS_F1")
object_list$endothelial_subset <- subset(object_list$endothelial_subset, subset = Sample_ID != "LBD_ATS_F4")

dataObject.endothelial <- SCTransform(object_list$endothelial_subset, verbose = TRUE)

dataObject.endothelial
# run PCA on the merged object
dataObject.endothelial <- RunPCA(object = dataObject.endothelial)

print("dataObject.endothelial")
dataObject.endothelial
Idents(dataObject.endothelial) <- "SampleID"

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
  paste0("../results/UMAP/reclusters/",
         projectID, "_recluster_endothelial.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()

dot_ind <- DotPlot(dataObject.endothelial,
                   features = markers.to.plot, 
                   cluster.idents = TRUE,
                   dot.scale = 8) + RotatedAxis()
pdf(
  paste0("../results/dot_plot/reclusters/",
         projectID, "_recluster_endothelial.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.endothelial, paste0("../rObjects/",projectID,"_endothelial.rds"), compress = FALSE)
rm(dataObject.endothelial)
