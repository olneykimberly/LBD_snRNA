#!/usr/bin/env Rscript
# ==============================================================================
# Subcluster markers within each cell type (Seurat v5 + Presto)
#   - For each cell_type: Find markers for cell_type_subcluster within that cell_type
#   - Save outputs by cell_type (astrocyte/, microglia/, neuron/, etc.)
# ==============================================================================

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(SeuratWrappers)
  library(ggplot2)
})

# ---- threads (avoid oversubscription) ----
cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
if (is.na(cpus) || cpus < 1) cpus <- 8
Sys.setenv(
  OMP_NUM_THREADS = as.character(cpus),
  OPENBLAS_NUM_THREADS = as.character(cpus),
  MKL_NUM_THREADS = as.character(cpus),
  VECLIB_MAXIMUM_THREADS = as.character(cpus),
  NUMEXPR_NUM_THREADS = as.character(cpus)
)

# ---- inputs ----
project_ID <- "CWOW_cellbender"
in_rds_path <- paste0("../rObjects/", project_ID, "_filtered_subclusters_pass2.rds")

# ---- outputs ----
out_root <- "../results/markers/cell_type_subclusters"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ---- load ----
message(Sys.time(), " | Loading Seurat object: ", in_rds_path)
dataObject <- readRDS(in_rds_path)

# ---- sanity checks ----
required_cols <- c("cell_type", "cell_type_subcluster")
missing_cols <- setdiff(required_cols, colnames(dataObject[[]]))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

DefaultAssay(dataObject) <- "RNA"

# Optional: quick UMAP by subcluster (saved once)
try({
  p_umap <- DimPlot(
    dataObject,
    reduction = "umap.rpca",
    group.by = "cell_type_subcluster",
    label = TRUE, repel = TRUE, raster = TRUE
  ) + ggtitle("UMAP by cell_type_subcluster")
  ggsave(
    filename = file.path(out_root, paste0(project_ID, "_UMAP_cell_type_subcluster.pdf")),
    plot = p_umap, width = 20, height = 10, limitsize = FALSE
  )
}, silent = TRUE)

# ---- iterate over cell types ----
cell_types <- sort(unique(as.character(dataObject$cell_type)))
message(Sys.time(), " | Found ", length(cell_types), " cell types.")

# thresholds (edit if desired)
min_pct <- 0.25
logfc_thr <- 0.25
only_pos <- TRUE

# "strict" filter settings
strict_log2fc <- 1
strict_q <- 0.05

for (ct in cell_types) {
  message(Sys.time(), " | Cell type: ", ct)
  
  # subset to this cell type
  obj_ct <- subset(dataObject, subset = cell_type == ct)
  
  # ensure there are multiple subclusters; otherwise skip
  subclust_levels <- sort(unique(as.character(obj_ct$cell_type_subcluster)))
  if (length(subclust_levels) < 2) {
    message(Sys.time(), " |  - Skipping ", ct, " (only ", length(subclust_levels), " subcluster).")
    next
  }
  
  # output dirs/files for this cell type
  ct_dir <- file.path(out_root, ct)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_markers_tsv <- file.path(ct_dir, paste0(project_ID, "_", ct, "_markers_presto.tsv"))
  out_markers_rds <- file.path(ct_dir, paste0(project_ID, "_", ct, "_markers_presto.rds"))
  
  out_strict_tsv <- file.path(ct_dir, paste0(project_ID, "_", ct, "_markers_presto_log2FC",
                                             strict_log2fc, "_q", strict_q, ".tsv"))
  out_strict_rds <- file.path(ct_dir, paste0(project_ID, "_", ct, "_markers_presto_log2FC",
                                             strict_log2fc, "_q", strict_q, ".rds"))
  
  # run presto within this cell type, grouped by subcluster
  # (this compares astrocyte_0 vs other astrocyte_* only, etc.)
  message(Sys.time(), " |  - RunPrestoAll within ", ct, " (n_subclusters=", length(subclust_levels), ") ...")
  markers_ct <- SeuratWrappers::RunPrestoAll(
    object = obj_ct,
    assay = "RNA",
    slot = "data",
    group.by = "cell_type_subcluster",
    min.pct = min_pct,
    logfc.threshold = logfc_thr,
    only.pos = only_pos
  )
  
  # Add helpful columns (optional)
  markers_ct <- markers_ct %>%
    dplyr::mutate(cell_type = ct)
  
  # save full
  saveRDS(markers_ct, out_markers_rds)
  write.table(markers_ct, out_markers_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # strict filter (optional, same logic you had)
  # NOTE: RunPrestoAll typically returns columns like: cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2, etc.
  if (all(c("cluster", "avg_log2FC", "p_val_adj") %in% colnames(markers_ct))) {
    markers_ct_strict <- markers_ct %>%
      dplyr::group_by(cluster) %>%
      dplyr::filter(avg_log2FC > strict_log2fc, p_val_adj < strict_q) %>%
      dplyr::arrange(cluster, dplyr::desc(avg_log2FC)) %>%
      dplyr::ungroup()
    
    saveRDS(markers_ct_strict, out_strict_rds)
    write.table(markers_ct_strict, out_strict_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    message(Sys.time(), " |  - Strict filter skipped (expected columns missing).")
  }
  
  message(Sys.time(), " |  - Saved outputs to: ", ct_dir)
}

message(Sys.time(), " | Done.")
