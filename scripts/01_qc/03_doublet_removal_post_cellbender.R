## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

## ----echo=FALSE, message=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender"
color.panel <- dittoColors()


## ----read_object--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read object
dataObject <- readRDS(file = paste0("../rObjects/", projectID, "_filtered.rds"))


## ----split_object-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# split object by sample
dataObject.split <- SplitObject(dataObject, split.by = "sample") 


## ----doubletFinder, message=FALSE, warning=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for (i in 1:length(dataObject.split)) {
  # normalize and find PCs
  print(i)
  obj_sample <- NormalizeData(dataObject.split[[i]])
  sampleID <- levels(droplevels(obj_sample@meta.data$sample))
  obj_sample <- FindVariableFeatures(obj_sample, selection.method = "vst", nfeatures = 2000)
  obj_sample <- ScaleData(obj_sample)
  obj_sample <- RunPCA(obj_sample)
  
  # get significant PCs
  stdv <- obj_sample[["pca"]]@stdev
  sum.stdv <- sum(obj_sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # run umap
  obj_sample <- RunUMAP(obj_sample, dims = 1:min.pc, reduction = "pca")
  
  # cluster
  obj_sample <- FindNeighbors(object = obj_sample, dims = 1:min.pc)                           
  obj_sample <- FindClusters(object = obj_sample, resolution = 0.2)
  
  # Assign identity of clusters
  Idents(object = obj_sample) <- "seurat_clusters"

  # number of cells in each cluster
  n_cells <- FetchData(obj_sample, vars = c("ident")) %>% dplyr::count(ident) %>%tidyr::spread(ident, n)
  
  ## pK Identification (no ground-truth) 
  sweep.res.list <- paramSweep(obj_sample, PCs = 1:min.pc, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK for any scRNA-seq data can be manually discerned as maxima in BCmvn distributions
  bcmvn_max <- bcmvn[which.max(bcmvn$BCmetric),]
  pK_value <- bcmvn_max$pK
  pK_value <- as.numeric(levels(pK_value))[pK_value]
  
  # Homotypic Doublet Proportion Estimate 
  annotations <- obj_sample@meta.data$individual_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp_poi <- round(pK_value*nrow(obj_sample@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  # Run DoubletFinder with varying classification  
  obj_sample <- doubletFinder(obj_sample, PCs = 1:min.pc, 
                pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, 
                reuse.pANN = FALSE, sct = FALSE)
  
  # set DF class for calling doublets
  DF_class <- obj_sample@meta.data[, grep("DF.classifications",colnames(obj_sample@meta.data)),]
  DF_class[which(DF_class == "Doublet")] <- "Doublet"
  table(DF_class)
  
  # table showing the number of doublets and singlets
  write.table(table(DF_class), paste0("../results/DoubletFinder/",projectID, "_doubletFinder_table_",sampleID), sep = "\t", 
              row.names = FALSE, quote = FALSE)
  obj_sample@meta.data[,"CellTypes_DF"] <- DF_class
  
  # plot
  d2 <- DimPlot(obj_sample, group.by="CellTypes_DF", reduction="umap",
          order=c("Coll.Duct.TC","Doublet"), 
          cols=c("#66C2A5","black"))
  path <- paste0("../results/DoubletFinder/",projectID,
               "_doubletFinder_UMAP_",sampleID)
  pdf(paste0(path, ".pdf"), width = 5,height = 4)
  print(d2)
  dev.off()
  
  # plot
  f1 <- FeaturePlot(obj_sample, 
            reduction = "umap", 
            features = c("nFeature_RNA", "nCount_RNA", 
                         "cell.complexity", "percent.mt"),
            pt.size = 0.4, 
            order = TRUE,
            label = TRUE)
  path <- paste0("../results/DoubletFinder/",projectID,
               "_doubletFinder_FeaturePlot_",sampleID)
  pdf(paste0(path, ".pdf"), width = 10, height = 7)
  print(f1)
  dev.off()
  
  #only keep singlets
  obj_sample_singlets <- subset(obj_sample, subset = CellTypes_DF == "Singlet")
  

  # number of cells in each cluster per and post removing doublets
  n_cells_singlets <- FetchData(obj_sample_singlets, vars = c("ident")) %>% dplyr::count(ident) %>% tidyr::spread(ident, n)
  n_cells_singlets
  # Find the column names that are in n_cells_singlets but not in n_cells
  missing_cols <- setdiff(colnames(n_cells), colnames(n_cells_singlets))

  # Add the missing column(s) to n_cells, setting their values to 0
  for (col in missing_cols) {
  n_cells_singlets[[col]] <- 0
  }

  # Now combine the two data frames
  ncells_per_cluster <- rbind(n_cells, n_cells_singlets)

  # View the combined result
  ncells_per_cluster

  row.names(ncells_per_cluster) <- c("Doublets and singlets", "Singlets only")
  ncells_per_cluster
  difference <- diff(as.matrix(ncells_per_cluster))
  difference <- as.data.frame(difference)
  row.names(difference) <- c("difference")
  cbind(difference, ncells_per_cluster)
  write.table(ncells_per_cluster, paste0(
    "../results/DoubletFinder/",projectID,
    "_doubletFinder_table_ncells_per_cluster_",sampleID, ".txt"), sep = "\t", 
    row.names = FALSE, quote = FALSE)
  # plot the number of cells in each cluster per and post doubletFinder
  ncell_matrix <- as.matrix(ncells_per_cluster)
  ncells_melt <- melt(ncell_matrix)
  colnames(ncells_melt) <- c("doublet type","cluster","number of cells")
  ncell_max <- ncells_melt[which.max(ncells_melt$`number of cells`),]
  ncell_max_value <- ncell_max$`number of cells`
  cellmax <- ncell_max_value + 800 # so that the figure doesn't cut off the text
  b1 <- ggplot(ncells_melt, aes(x = factor(cluster), y = `number of cells`,
                          fill = `doublet type`)) + 
    geom_bar(stat="identity", colour="black", width=1, position = position_dodge(width=0.8)) +
    geom_text(aes(label = `number of cells`), 
              position=position_dodge(width=0.9), vjust=-0.25, angle = 45, hjust=-.01) + 
    theme_classic() + scale_fill_manual(values = c("gray", "#66C2A5")) +
    ggtitle("Number of cells per cluster") +  xlab("cluster") +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    scale_y_continuous(limits = c(0,cellmax))
  path <- paste0("../results/DoubletFinder/",projectID,
               "_doubletFinder_barplot_ncells_per_cluster_",sampleID)
  pdf(paste0(path, ".pdf"), width = 7,height = 5)
  print(b1)
  dev.off()
  f2 <- FeaturePlot(obj_sample_singlets, 
            reduction = "umap", 
            features = c("nFeature_RNA", "nCount_RNA", 
                         "cell.complexity", "percent.mt"),
            pt.size = 0.4, 
            order = TRUE,
            label = TRUE)
  path <- paste0("../results/DoubletFinder/",projectID,
               "_doubletFinder_FeaturePlot_singlets_",sampleID)
  pdf(paste0(path, ".pdf"), width = 10,height = 7)
  print(f2)
  dev.off()
  
  # put the PMI together again
  dataObject.split[[i]] <- obj_sample_singlets
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# converge obj.split
dataObject.singlets <- merge(x = dataObject.split[[1]],
                             y = c(dataObject.split[[2]], dataObject.split[[3]], dataObject.split[[4]], dataObject.split[[5]], dataObject.split[[6]],
                                   dataObject.split[[7]], dataObject.split[[8]], dataObject.split[[9]], dataObject.split[[10]], dataObject.split[[11]],
                                   dataObject.split[[12]], dataObject.split[[13]], dataObject.split[[14]], dataObject.split[[15]], dataObject.split[[16]],
                                   dataObject.split[[17]], dataObject.split[[18]], dataObject.split[[19]], dataObject.split[[20]], dataObject.split[[21]], 
                                   dataObject.split[[22]], dataObject.split[[23]], dataObject.split[[24]], dataObject.split[[25]], dataObject.split[[26]],
                                   dataObject.split[[27]], dataObject.split[[28]], dataObject.split[[29]], dataObject.split[[30]], dataObject.split[[31]],
                                   dataObject.split[[32]], dataObject.split[[33]], dataObject.split[[34]], dataObject.split[[35]], dataObject.split[[36]],
                                   dataObject.split[[37]], dataObject.split[[38]], dataObject.split[[39]]),
                             project = paste0("cwow"))

# print how many cells removed
print(paste0(dim(dataObject)[2] - dim(dataObject.singlets)[2]," nuclie removed"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(dataObject.singlets, paste0("../rObjects/",projectID,"_singlets.rds"))
#dataObject.singlets <- readRDS(paste0("../rObjects/",projectID,"_singlets.rds"))
# Insepect
dataObject.singlets
dataObject <- dataObject.singlets


## ----reprocess----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# transform
dataObject <- SCTransform(dataObject, verbose = FALSE)

# run PCA on the merged object
dataObject <- RunPCA(object = dataObject)
Idents(dataObject) <- "Sample_ID"

# re-join layers
dataObject[["RNA"]] <- JoinLayers(dataObject[["RNA"]])

# Determine the K-nearest neighbor graph
dataObject <- FindNeighbors(object = dataObject, 
                                     assay = "SCT", # set as default after SCTransform
                                     reduction = "pca", # pca, harmony 
                                     dims = 1:15)
# Run UMAP
dataObject <- RunUMAP(dataObject,
                               dims = 1:15,
                               reduction = "pca",
                               n.components = 3) # set to 3 to use with VR

# Determine the clusters for various resolutions
dataObject <- FindClusters(object = dataObject,
                                 algorithm = 1, # 1= Louvain
                                 resolution = seq(0.2,1,by=0.2))

saveRDS(dataObject, paste0("../rObjects/",projectID,"_unannotated_doublets_removed.rds"))

## ----umap_noDoublets----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Idents(dataObject) <- dataObject$SCT_snn_res.0.6
dataObject$seurat_clusters <- dataObject$SCT_snn_res.0.6
ditto_umap <- dittoDimPlot(object = dataObject,
             var = "seurat_clusters",
             reduction.use = "umap",
             do.label = TRUE,
             labels.highlight = TRUE)

path <- paste0("../results/UMAP/unannotated/",projectID,
               "_UMAP_unannotated_doublets_removed")
ditto_umap
saveToPDF(paste0(path, ".pdf"), width = 7, height = 6.6)

## ----harmony_int----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(45)
dataObject.integrated <- IntegrateLayers(
  object = dataObject, method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
# Determine the K-nearest neighbor graph
dataObject.integrated <- FindNeighbors(object = dataObject.integrated, 
                                       reduction = "harmony", # pca, harmony 
                                       dims = 1:30)

# Determine the clusters for various resolutions
dataObject.integrated <- FindClusters(object = dataObject.integrated, resolution = 0.2)
dataObject.integrated <- RunUMAP(dataObject.integrated, reduction = "harmony", dims = 1:30)

p1 <- DimPlot(
  dataObject.integrated,
  reduction = "harmony",
  group.by = c("Sample_ID"),
  combine = FALSE, label.size = 2
)
p1
## ----nuclei_per_cluster_noDoublets--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Idents(dataObject) <- dataObject$seurat_clusters
sample_ncells <- FetchData(dataObject, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident,sample) %>%
  tidyr::spread(ident, n)
write.table(sample_ncells, 
            paste0("../results/nuclei_count/",
                   projectID, 
                   "_nuclei_per_cluster_doublets_removed.txt"),
            quote = FALSE, sep = "\t")
sample_ncells


## ----dot_individual-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
markers.to.plot <-
  c(
"CLU", 
"GFAP", 
"AQP4", 
"GJA1",
"CLDN5",
"ADGRF5",
"FLT1",
"COL1A1",
"COL1A2",
"DCN",
"HEXB",
"C1QA",
"C1QB",
"C1QC",
"TMEM119",
"ITGAM",
"TYROBP",
"P2RY12",
"AIF1",
"RBFOX1",
"RBFOX3", 
"SNAP25",
"SYT1",
"GAD1",
"GAD2",
"PLP1",
"MBP", 
"MOG", 
"OLIG1",
"PDGFRA",
"VCAN",
"TNR",
"ACTA2",
"VTN"
  )

dot_ind <- DotPlot(dataObject,
                   features = markers.to.plot, 
                   cluster.idents = TRUE,
                   dot.scale = 8) + RotatedAxis()

## ----save_dot_individual, echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf(
  paste0(
    "../results/dot_plot/",
    projectID,
    "_clusters_DotPlot_no_integration_doublets_removed.pdf"
  ),
  width = 14,
  height = 10
)
dot_ind
dev.off()


## ----save_object,echo=FALSE,eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(dataObject, paste0("../rObjects/",projectID,"_unannotated_doublets_removed_harmony_int.rds"))
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------