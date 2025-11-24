knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender_RPCAIntegration_clean_SCT_annotation"
color.panel <- dittoColors()

install.packages("ShinyCell")
BiocManager::install("ShinyCell", lib="/tgen_labs/jfryer/kolney/R/rstudio-with_modules-4.4.0-3.sif")

reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
           "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
           "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}
devtools::install_github("SGDDNB/ShinyCell")

# If you are using h5ad file as input, run the code below as well
# reticulate::py_install("anndata")


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
metadata <- metadata[c(138,1:5, 22)]
sc.config <- createConfig(obj = dataObject, meta.to.include = metadata)

makeShinyApp(obj = dataObject,
             scConf = sc.config, 
             gex.assay = "SCT", 
           #  gene.mapping = TRUE,
             defPtSiz = 0.75,
            # default.gene1 = "SNAP25",
            # default.gene2 = "PLP1",
             default.dimred = "integrated.rpca",
             shiny.dir = paste0("../shiny_apps/LBD_CWOW_snRNA_annotation"),
             shiny.title = "LBD CWOW snRNA; n = 35")
