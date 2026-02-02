#!/usr/bin/env Rscript
# =============================================================================
# Cell-type score filtering diagnostics
# - Uses existing module scores in obj meta.data:
#     Score_<celltype>_1, Contam_<type>_1
# - Computes per-cell max_contam, per-cell-type thresholds, KEEP/DROP labels
# - Saves violin plots + summary tables per cell type
#
# Output:
#   ../results/subclustering_cell_type/<cell_type>/
#     - *_Vln_targetScore.pdf
#     - *_Vln_maxContam.pdf
#     - *_Vln_markers.pdf
#     - *_Vln_contam_markers.pdf
#   ../results/subclustering_cell_type/score_filter_summary_by_celltype.csv
#   ../rObjects/<project_ID>_with_score_filter_labels.rds
# =============================================================================

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

project_ID <- "CWOW_cellbender"
obj <- readRDS(paste0("../rObjects/", project_ID, "_annotated_with_celltype_subclusters.rds"))

DimPlot(obj, reduction = "umap.rpca", group.by = "cell_type",
        label = TRUE, repel = TRUE, raster = TRUE)
# Subclusters to remove (edit this)
bad_subclusters <- c(
  "astrocyte_11",
  "astrocyte_21", 
  "endothelial_4", 
  "fibroblast_7", 
  "fibroblast_4", 
  "fibroblast_2", 
  "interneuron_13", 
  "interneuron_12", 
  "microglia_6", 
  "microglia_12", 
  "microglia_8", 
  "microglia_14", 
  "mural_4", 
  "mural_7", 
  "oligodendrocyte_16",
  "opc_15",
  "opc_16",
  "opc_14",
  "opc_13"
)

# how many cells will drop
sum(obj$cell_type_subcluster %in% bad_subclusters)         

# ensures subtype name is in object
bad_present <- intersect(bad_subclusters, unique(obj$cell_type_subcluster))
cells_keep <- colnames(obj)[!(obj$cell_type_subcluster %in% bad_present)]
obj_final <- subset(obj, cells = cells_keep) # subset to remove bad subclusters
p_umap <- DimPlot(obj_final, reduction = "umap.rpca", group.by = "cell_type",
        label = TRUE, repel = TRUE, raster = TRUE) +
  ggtitle("UMAP by cell type")
ggsave(file.path(paste0("../results/UMAP/", project_ID, "_UMAP_by_cell_type_post_pass1_cleaning.pdf")),
       p_umap, width = 11, height = 7, units = "in")

# Save
out_rds <- paste0("../rObjects/", project_ID, "_filtered_subclusters.rds")
saveRDS(obj_final, out_rds, compress = FALSE)

message("Saved: ", out_rds)
obj_final