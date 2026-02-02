#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Build ShinyCell app (Seurat v5 -> v4-style bridge) and limit expression matrix
# to only the top 8,000 VariableFeatures to reduce sc1gexpr.h5 size.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Setup & Object Cleaning
# ------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ShinyCell)
  library(dittoSeq)
})

# ------------------------------------------------------------------------------
# load object + join layers 
# ------------------------------------------------------------------------------
dataObject <- readRDS("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/rObjects/subclusters/CWOW_cellbender_interneuron_subclustered.rds")
Idents(dataObject) <- dataObject$subcluster

# bad_subclusters <- c(
#   "neuron_9" # Basically a single sample in this cluster
# )
# 
# # how many cells will drop
# sum(dataObject$subcluster %in% bad_subclusters)         
# 
# # ensures subtype name is in object
# bad_present <- intersect(bad_subclusters, unique(dataObject$subcluster))
# cells_keep <- colnames(dataObject)[!(dataObject$subcluster %in% bad_present)]
# obj_final <- subset(dataObject, cells = cells_keep) # subset to remove bad subclusters
obj_final <- dataObject

UMAP_plot <- DimPlot(obj_final, reduction = "umap.ct", label = TRUE, label.size = 3) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(
    axis.title.x = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.title = element_text(size = 8),
    plot.margin = margin(0.5, 0.025, 0.25, 0.025, "cm")
  )
UMAP_plot
DefaultAssay(obj_final)

# ------------------------------------------------------------------------------
# Metadata 
# ------------------------------------------------------------------------------
meta_keep <- c(
  "orig.ident", "cell_type", "subcluster", "seurat_clusters", "group", "Sample_ID",
  "nCount_RNA", "nFeature_RNA", "Braak.NFT", "Thal.amyloid", "Cing.LB",
  "Age", "sex_inferred"
)

# Only include columns that exist
meta_exists <- meta_keep[meta_keep %in% colnames(obj_final@meta.data)]
cat("Including metadata columns:", paste(meta_exists, collapse = ", "), "\n")

sc.config <- createConfig(obj_final, meta.to.include = meta_exists)

# ------------------------------------------------------------------------------
# Build the App
# ------------------------------------------------------------------------------
shiny_path <- "../shiny_apps/LBD_CWOW_snRNAseq_Inhibitory_Neurons_human"
if (dir.exists(shiny_path)) {
  unlink(shiny_path, recursive = TRUE)
}

makeShinyApp(
  obj = obj_final,
  scConf = sc.config,
  gex.assay = "RNA",
  gex.slot = "data",         # ShinyCell expects slot-like semantics
  default.dimred = "umap.ct",
  default.gene1 = "SNAP25",
  default.gene2 = "RBFOX3",
  shiny.dir = shiny_path,
  shiny.title = "LBD CWOW human snRNA Inhibitory Neurons; n = 35"
)

message("Done")