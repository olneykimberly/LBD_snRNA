## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender_after_recluster_harm_int_noise_removed_after_annotation"
color.panel <- dittoColors()

# read object
dataObject  <- readRDS(paste0("../rObjects/",projectID,".rds"))
DefaultAssay(object = dataObject) <- "RNA"
dataObject[["RNA"]] <- JoinLayers(dataObject[["RNA"]])
dataObject <- NormalizeData(dataObject)

dataObject$celltype.stim <- paste(dataObject$individual_clusters, dataObject$group, sep = "_")
Idents(dataObject) <- "celltype.stim"
VlnPlot(dataObject, features <- c('DPYD'), idents = c("microglia_AD_AT", "microglia_CONTROL"), group.by = "Sample_ID", ncol = 1) 
VlnPlot(dataObject, features <- c('ACSL1'), idents = c("microglia_AD_AT", "microglia_CONTROL"), group.by = "Sample_ID", ncol = 1) 

dataObject.pseudo <- AggregateExpression(
  dataObject, 
  assays = "RNA", # DESeq works with raw counts
  #  features = protein_coding_genes$gene_name,
  return.seurat = TRUE, # If return.seurat = TRUE, aggregated values are placed in the 'counts' layer of the returned object. 
  # The data is then normalized by running NormalizeData on the aggregated counts. 
  # ScaleData is then run on the default assay before returning the object.
  group.by = c("sample", "group", "individual_clusters")
)
# group comparison with covariates
dataObject.pseudo$allcells.group <- paste(dataObject.pseudo$individual_clusters, dataObject.pseudo$group, sep = "_")
Idents(dataObject.pseudo) <- "allcells.group"
VlnPlot(dataObject.pseudo, features <- c('DPYD'), idents = c("microglia_AD-AT", "microglia_CONTROL"),  ncol = 1) 
VlnPlot(dataObject.pseudo, features <- c('ACSL1'), idents = c("microglia_AD-AT", "microglia_CONTROL"),  ncol = 1) 

cell_types <- c(unique(dataObject$individual_clusters))
group_types <- c(unique(dataObject$group))


#------------------
#dataObject <- PrepSCTFindMarkers(dataObject)


DE <- FindMarkers(dataObject, ident.1 = "microglia_AD_AT", ident.2 = "microglia_CONTROL", verbose = FALSE, latent.vars = c("sex"), test.use = "DESeq2")
DE$gene <- rownames(DE)
write.table(DE, 
            paste0("../results/DEGs/AD_vs_CONTROL_MAST_markers.tsv"),
            quote = FALSE,
            row.names = FALSE)
sig <- subset(DE, p_val_adj < 0.05)

