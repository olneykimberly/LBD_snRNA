#!/usr/bin/env Rscript
# ==============================================================================
# Cluster markers for annotation (Seurat v5)
# 1) Load FINAL object (reproClustersUMAP)
# 2) Make RNA "data" exist by JoinLayers + NormalizeData ONCE (checkpoint)
# 3) Run fast markers with SeuratWrappers::RunPrestoAll on RNA log-normalized data
# ==============================================================================

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

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

# ---- paths ----
project_ID <- "CWOW_cellbender"
in_rds <- paste0("../rObjects/", project_ID, "_FULL_projected_integratedrpca_reproClustersUMAP.rds")

out_rds_join_norm <- paste0("../rObjects/", project_ID, "_FULL_RNA_joined_normalized.rds")

dir.create("../results/markers", recursive = TRUE, showWarnings = FALSE)
out_markers_tsv <- paste0("../results/markers/", project_ID, "_markers_presto.tsv")
out_markers_rds <- paste0("../rObjects/", project_ID, "_markers_presto.rds")

out_strict_rds <- paste0("../rObjects/", project_ID, "_markers_presto_log2FC1_q0.05.rds")
out_strict_tsv <- paste0("../results/markers/", project_ID, "_markers_presto_log2FC1_q0.05.tsv")

# ---- load ----
message(Sys.time(), " | Loading: ", in_rds)
dataObject <- readRDS(in_rds)

# ---- build RNA data once (JoinLayers + NormalizeData) ----
if (file.exists(out_rds_join_norm)) {
  message(Sys.time(), " | Loading joined+normalized checkpoint: ", out_rds_join_norm)
  dataObject <- readRDS(out_rds_join_norm)
} else {
  message(Sys.time(), " | JoinLayers + NormalizeData (one-time) ...")
  DefaultAssay(dataObject) <- "RNA"
  dataObject <- JoinLayers(dataObject)
  dataObject <- NormalizeData(dataObject, normalization.method = "LogNormalize", verbose = FALSE)
  
  message(Sys.time(), " | Saving joined+normalized checkpoint: ", out_rds_join_norm)
  saveRDS(dataObject, out_rds_join_norm, compress = FALSE)
}

# ---- markers ----
DefaultAssay(dataObject) <- "RNA"
Idents(dataObject) <- "seurat_clusters"

message(Sys.time(), " | RunPrestoAll (RNA, slot=data) ...")
markers <- SeuratWrappers::RunPrestoAll(
  object = dataObject,
  assay = "RNA",
  slot = "data",
  group.by = "seurat_clusters",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  only.pos = TRUE
)

message(Sys.time(), " | Saving markers: ", out_markers_rds)
saveRDS(markers, out_markers_rds)

message(Sys.time(), " | Writing markers: ", out_markers_tsv)
write.table(markers, out_markers_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# ---- strict filter (optional, kept because you had it) ----
all.markers.strict <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  dplyr::arrange(cluster, dplyr::desc(avg_log2FC))

saveRDS(all.markers.strict, out_strict_rds)
write.table(all.markers.strict, out_strict_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

message(Sys.time(), " | Done.")
