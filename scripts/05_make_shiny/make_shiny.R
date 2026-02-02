#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Build ShinyCell app (Seurat v5 -> v4-style bridge) and limit expression matrix
# to only the top 8,000 VariableFeatures to reduce sc1gexpr.h5 size.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1. Setup & Object Cleaning
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
# 2. Load object + join layers (Seurat v5)
# ------------------------------------------------------------------------------
dataObject <- readRDS("../rObjects/CWOW_cellbender_filtered_subclusters_pass2.rds")
dataObject <- JoinLayers(dataObject)

# Standardize Factor Levels
order_cell_type <- c("neuron", "interneuron", "oligodendrocyte", "opc",
                     "astrocyte", "microglia", "endothelial", "fibroblast", "mural")
dataObject$cell_type <- factor(as.character(dataObject$cell_type), levels = order_cell_type)

# ------------------------------------------------------------------------------
# 3. Create Clean V4-style Object using layer-aware getters/setters
#    (bypasses Seurat v5 Assay5 quirks in ShinyCell)
# ------------------------------------------------------------------------------
# Pull counts + log-normalized data from the v5 object
rna_counts <- GetAssayData(dataObject, assay = "RNA", layer = "counts")
rna_data   <- GetAssayData(dataObject, assay = "RNA", layer = "data")

# Build clean object
cleanObj <- CreateSeuratObject(
  counts = rna_counts,
  meta.data = dataObject@meta.data
)

# Ensure RNA assay exists and set both counts + data layers explicitly
cleanObj[["RNA"]] <- CreateAssayObject(counts = rna_counts)
cleanObj <- SetAssayData(
  cleanObj,
  assay = "RNA",
  layer = "data",
  new.data = rna_data
)

# ------------------------------------------------------------------------------
# 4. Limit genes to ONLY the top 8,000 VariableFeatures
#    This is the key step to shrinking sc1gexpr.h5.
# ------------------------------------------------------------------------------
# Prefer variable features from the original object if present; otherwise compute.
vf <- VariableFeatures(dataObject)

# Compute variable features if missing/empty
if (length(vf) == 0) {
  message("No VariableFeatures found in dataObject; computing VariableFeatures (nfeatures=8000) on cleanObj ...")
  DefaultAssay(cleanObj) <- "RNA"
  # Ensure data layer exists; if not, create it (safety)
  if (nrow(GetAssayData(cleanObj, assay = "RNA", layer = "data")) == 0) {
    cleanObj <- NormalizeData(cleanObj, verbose = FALSE)
  }
  cleanObj <- FindVariableFeatures(cleanObj, nfeatures = 8000, verbose = FALSE)
  vf <- VariableFeatures(cleanObj)
} else {
  # Keep only genes that actually exist in cleanObj (safety)
  vf <- intersect(vf, rownames(cleanObj))
  if (length(vf) < 8000) {
    message("VariableFeatures(dataObject) has <8000 genes present in cleanObj (n=", length(vf),
            "). Will use available genes.")
  } else {
    vf <- vf[seq_len(8000)]
  }
  VariableFeatures(cleanObj) <- vf
}

# Subset the Seurat object to only those 8,000 variable genes.
# IMPORTANT: This reduces both counts + data layers to those genes.
message("Subsetting cleanObj to VariableFeatures only (n=", length(vf), ") ...")
cleanObj <- subset(cleanObj, features = vf)

# ------------------------------------------------------------------------------
# 5. Fix UMAP Coordinates (The Viewport fix)
# ------------------------------------------------------------------------------
umap_coords <- Embeddings(dataObject, "umap.rpca")
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# Ensure rownames match the cell names in cleanObj
rownames(umap_coords) <- colnames(dataObject)
umap_coords <- umap_coords[colnames(cleanObj), , drop = FALSE]

cleanObj[["umap"]] <- CreateDimReducObject(
  embeddings = umap_coords,
  key = "UMAP_",
  assay = "RNA"
)

# ------------------------------------------------------------------------------
# 6. Manual Config Reconstruction (The "4 vs 3 items" fix)
# ------------------------------------------------------------------------------
meta_keep <- c(
  "orig.ident", "cell_type", "cell_type_subcluster", "group", "Sample_ID",
  "nCount_RNA", "nFeature_RNA", "Braak.NFT", "Thal.amyloid", "Cing.LB",
  "Age", "sex_inferred"
)

# Only include columns that exist
meta_exists <- meta_keep[meta_keep %in% colnames(cleanObj@meta.data)]
cat("Including metadata columns:", paste(meta_exists, collapse = ", "), "\n")

sc.config <- createConfig(cleanObj, meta.to.include = meta_exists)

# 1) Force reduction name (string)
sc.config$dimred <- "umap"
# Force cell_type to be the default coloring on UMAP
if ("cell_type" %in% sc.config$meta$ID) {
  sc.config$meta$default <- "cell_type"
}

# 2) Build color string for cell_type
actual_levels <- levels(droplevels(cleanObj$cell_type))
color_string <- paste(dittoColors()[seq_len(length(actual_levels))], collapse = "|")

# 3) Inject colors safely
if ("cell_type" %in% sc.config$meta$ID) {
  sc.config$meta$cols[sc.config$meta$ID == "cell_type"] <- color_string
}

# Optional: quick UMAP check
Idents(cleanObj) <- cleanObj$cell_type
UMAP_plot <- DimPlot(cleanObj, reduction = "umap", label = FALSE, label.size = 3) +
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
print(UMAP_plot)

# ------------------------------------------------------------------------------
# 7. Build the App
# ------------------------------------------------------------------------------
shiny_path <- "../shiny_apps/LBD_CWOW_snRNAseq_human"
if (dir.exists(shiny_path)) {
  unlink(shiny_path, recursive = TRUE)
}

makeShinyApp(
  obj = cleanObj,
  scConf = sc.config,
  gex.assay = "RNA",
  gex.slot = "data",         # ShinyCell expects slot-like semantics
  default.dimred = "umap",
  default.gene1 = "P2RY12",
  default.gene2 = "APOE",
  shiny.dir = shiny_path,
  shiny.title = "LBD CWOW human snRNA; n = 35"
)

message("Done. Shiny app written to: ", normalizePath(shiny_path))
message("NOTE: sc1gexpr.h5 size should now be substantially smaller (8,000 genes only).")
