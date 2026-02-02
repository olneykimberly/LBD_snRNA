knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

#source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender_RPCAIntegration_clean_RNA_annotation"

.libPaths(c("/tgen_labs/jfryer/kolney/R/rstudio-with_modules-4.4.0-3.sif", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

.libPaths()
unloadNamespace("RSpectra")
unloadNamespace("rtracklayer")
unloadNamespace("GenomicAlignments")
unloadNamespace("SummarizedExperiment")
unloadNamespace("DelayedArray")
unloadNamespace("SparseArray")
unloadNamespace("S4Arrays")
library(Matrix, lib.loc = "/usr/local/lib/R/site-library")
library(SeuratObject) 
#library(Signac)
library(Seurat) 
library(stringr)
library(ggplot2)
#library(harmony)
#library(remaCor)
library(gridExtra)
library(grid)
library(lattice)
library(R.utils)


library(ggdendro)
library(ShinyCell)

dataObject <- readRDS(paste0("../rObjects/",projectID,".rds"))
#DimPlot(dataObject,
#        group.by = "individual_clusters", reduction = "integrated.rpca")

# inspect 
#dataObject[["RNA"]] <- JoinLayers(dataObject[["RNA"]])
#dataObject.PSCT <- PrepSCTFindMarkers(dataObject, assay = "SCT", verbose = TRUE)

metadata <- colnames(dataObject@meta.data)
df <- as.data.frame(metadata)
table(dataObject@meta.data$TYPE)
#metadata <- metadata[c(138,1:5, 22)]
metadata <- metadata[c(138,1:5, 22,41,23,38,26,27,31,42,43,44,48,50,51,52,53)]

sc.config <- createConfig(obj = dataObject, meta.to.include = metadata)

makeShinyApp(obj = dataObject,
             scConf = sc.config, 
             gex.assay = "RNA", 
           #  gene.mapping = TRUE,
             defPtSiz = 0.75,
            # default.gene1 = "SNAP25",
            # default.gene2 = "PLP1",
             default.dimred = "integrated.rpca",
             shiny.dir = paste0("../shiny_apps/LBD_CWOW_snRNA_annotation_RNA"),
             shiny.title = "LBD CWOW snRNA; n = 35")
