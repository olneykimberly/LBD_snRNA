knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender_RPCAIntegration_clean_SCT_annotation"
color.panel <- dittoColors()

library(ggdendro)
library(ShinyCell)

dataObject.annotated <- readRDS(paste0("../rObjects/",projectID,".rds"))
# inspect 
#dataObject.annotated[["RNA"]] <- JoinLayers(dataObject.annotated[["RNA"]])

metadata <- colnames(dataObject.annotated@meta.data)
df <- as.data.frame(metadata)
metadata <- metadata[c(37,1:17)]
sc.config <- createConfig(obj = dataObject.annotated, meta.to.include = metadata)

makeShinyApp(obj = dataObject.annotated,
             scConf = sc.config, 
             gex.assay = "SCT", 
             gene.mapping = TRUE,
             defPtSiz = 0.75,
             default.gene1 = "RBFOX1",
             default.gene2 = "PLP1",
             default.dimred = "integrated.rpca",
             shiny.dir = paste0("../shiny_apps/rpca_annotations"),
             shiny.title = "CWOW rpca annotations")
