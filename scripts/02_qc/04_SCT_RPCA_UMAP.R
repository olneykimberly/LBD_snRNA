#!/usr/bin/env Rscript
# ==============================================================================
# CWOW Seurat v5.2.1 integration + fully reproducible clustering/UMAP 
# What this script does (and why):
#  1) Builds a balanced "sketch" (a representative subset of cells per sample)
#     using SketchData(). This makes integration feasible at your scale while
#     preserving sample diversity (critical for correct cell type structure).
#  2) Runs SCTransform+PCA+IntegrateLayers (RPCAIntegration) on the sketch assay
#     (atomic cells) to learn an integration space robust to batch/sample effects.
#  3) Projects ("unsketches") the learned integrated embedding to ALL cells
#     using ProjectIntegration(), producing an integrated embedding for the full
#     dataset (integrated.rpca) while keeping layered counts intact.
#  4) Reruns FindNeighbors/FindClusters/RunUMAP with deterministic settings to
#     guarantee identical results across reruns (no future RNG warnings).
#  5) Writes the same style of diagnostics as your standalone reproducible script:
#        - UMAP PDFs (by cluster, by Sample_ID)
#        - cluster count tables (overall and by sample)
#        - a text summary
#     and saves a final RDS that includes the reproducible clusters + UMAP.
#
# Reproducibility principles used here:
#  - Heavy steps (integration + projection) are checkpointed to avoid re-running.
#  - Neighbors/clusters/UMAP are forced to a clean slate before recomputing.
#  - UMAP is forced to single-thread (uwot n.threads=1) and fixed seeds are used.
#  - Future plan is set to sequential during stochastic steps to avoid
#    future.apply RNG misuse warnings.
# ==============================================================================

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(ggplot2)
})

# ----------------------------- SLURM resources --------------------------------
# Expected SLURM: 8 cores, 250GB. We cap BLAS/OpenMP threads to cpus.
cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
cpus <- ifelse(is.na(cpus) || cpus < 1, 8, cpus)

Sys.setenv(
  OMP_NUM_THREADS = as.character(cpus),
  OPENBLAS_NUM_THREADS = as.character(cpus),
  MKL_NUM_THREADS = as.character(cpus),
  VECLIB_MAXIMUM_THREADS = as.character(cpus),
  NUMEXPR_NUM_THREADS = as.character(cpus)
)

# futures memory cap
options(future.globals.maxSize = 500 * 1024^3)  # 500 GB

# --------------------------- User parameters ----------------------------------
projectID <- "CWOW_cellbender"
in_rds <- paste0("../rObjects/", projectID, "_singlets_scDblFinder_exprate.rds")

split_var <- "Sample_ID"

# sketching
cells_per_sample_sketch <- 5000
nfeatures <- 3000
npcs <- 50

# integration / clustering / UMAP
dims_use <- 1:30
resolution_use <- 0.3
reduction_rpca <- "integrated.rpca"
umap_name <- "umap.rpca"
graph_nn  <- "integratedrpca_nn"
graph_snn <- "integratedrpca_snn"

# reproducibility seeds
seed_global <- 12345
seed_umap   <- 4242

# ----------------------------- Output dirs ------------------------------------
out_dir <- "../results/final_full_rpca_diagnostics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- Checkpoints ------------------------------------
# Heavy steps checkpoints
out_rds_obj_sketch_integrated <- paste0("../rObjects/", projectID, "_FULL_sketchAssay_integrateddr.rds")
out_rds_full_proj_ckpt        <- paste0("../rObjects/", projectID, "_FULL_projected_integratedrpca_ckpt.rds")

# Final output object (with reproducible neighbors/clusters/umap)
out_rds_full_proj_final <- paste0("../rObjects/", projectID, "_FULL_projected_integratedrpca_reproClustersUMAP.rds")

# Diagnostics outputs (FULL, not sketch)
out_umap_by_cluster    <- file.path(out_dir, paste0(projectID, "_FULL_UMAP_by_cluster.pdf"))
out_umap_by_sample     <- file.path(out_dir, paste0(projectID, "_FULL_UMAP_by_", split_var, ".pdf"))
out_counts_overall     <- file.path(out_dir, paste0(projectID, "_FULL_cluster_counts_overall.csv"))
out_counts_by_sample   <- file.path(out_dir, paste0(projectID, "_FULL_cluster_counts_by_sample.csv"))
out_summary_txt        <- file.path(out_dir, paste0(projectID, "_FULL_cluster_summary.txt"))

# ----------------------------- Parallel plan ----------------------------------
# Use futures for heavy steps if needed (SketchData/SCTransform can benefit),
# but we will force sequential for stochastic steps later.
n_workers <- max(1, cpus - 1)
if (Sys.info()[["sysname"]] == "Linux" && Sys.getenv("RSTUDIO") == "") {
  plan(multicore, workers = n_workers)
  message(Sys.time(), " | future plan: multicore | workers=", n_workers)
} else {
  plan(multisession, workers = n_workers)
  message(Sys.time(), " | future plan: multisession | workers=", n_workers)
}

# Global RNG seed
set.seed(seed_global)

message(Sys.time(), " | SLURM_CPUS_PER_TASK=", cpus)
message(Sys.time(), " | Loading Seurat object: ", in_rds)

# ------------------------------ Load FULL object ------------------------------
obj <- readRDS(in_rds)

stopifnot(as.character(packageVersion("Seurat")) == "5.2.1")
stopifnot(inherits(obj, "Seurat"))
stopifnot("RNA" %in% names(obj@assays))

DefaultAssay(obj) <- "RNA"

meta <- obj[[]]
stopifnot(is.data.frame(meta))
stopifnot(split_var %in% colnames(meta))

message(Sys.time(), " | Cells(nuclei)=", ncol(obj),
        " | Features=", nrow(obj),
        " | Unique ", split_var, "=", length(unique(meta[[split_var]])))

rna_layers <- tryCatch(Layers(obj[["RNA"]]), error = function(e) character(0))
message(Sys.time(), " | RNA layers present: n=", length(rna_layers),
        " | example: ", paste(head(rna_layers, 5), collapse = ", "))

# ------------------------------ Build sketch cell set --------------------------
message(Sys.time(), " | Building balanced sketch cell set via meta.data ('", split_var, "') ...")
samples <- sort(unique(meta[[split_var]]))

cells.sketch <- unique(unlist(lapply(samples, function(s) {
  cells_s <- rownames(meta)[meta[[split_var]] == s]
  if (length(cells_s) == 0) return(character(0))
  if (length(cells_s) <= cells_per_sample_sketch) return(cells_s)
  sample(cells_s, cells_per_sample_sketch)
})))

message(Sys.time(), " | Sketch cells selected: ", length(cells.sketch),
        " (target ~", cells_per_sample_sketch, " x ", length(samples), ")")
if (length(cells.sketch) < 5000) stop("Too few sketch cells selected.")

# ------------------------------ SketchData ------------------------------------
# Always ensure sketch assay exists (needed for downstream steps and for re-runs).
message(Sys.time(), " | Running SketchData (creates 'sketch' assay on FULL object) ...")
obj <- SketchData(
  object = obj,
  cells = cells.sketch,
  assay = "RNA",
  verbose = FALSE
)
gc()

stopifnot("sketch" %in% names(obj@assays))
sk_layers <- tryCatch(Layers(obj[["sketch"]]), error = function(e) character(0))
message(Sys.time(), " | sketch assay layers present: n=", length(sk_layers),
        " | example: ", paste(head(sk_layers, 5), collapse = ", "))

# ------------------------------ Integrate on sketch assay ---------------------
if (file.exists(out_rds_obj_sketch_integrated)) {
  message(Sys.time(), " | Found FULL sketch-assay integration checkpoint. Loading: ",
          out_rds_obj_sketch_integrated)
  obj <- readRDS(out_rds_obj_sketch_integrated)
} else {
  message(Sys.time(), " | Running SCT+PCA on FULL object's SKETCH assay (atomic) ...")
  old_plan <- future::plan()
  future::plan(sequential)  # keep deterministic + reduce RNG/future issues in this block
  
  DefaultAssay(obj) <- "sketch"
  
  obj <- SCTransform(
    obj,
    assay = "sketch",
    method = "glmGamPoi",
    vst.flavor = "v2",
    variable.features.n = nfeatures,
    verbose = FALSE
  )
  gc()
  
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  gc()
  
  message(Sys.time(), " | IntegrateLayers on FULL object (RPCAIntegration; assay='sketch') ...")
  obj <- IntegrateLayers(
    object = obj,
    method = RPCAIntegration,
    normalization.method = "SCT",
    verbose = FALSE
  )
  gc()
  
  message(Sys.time(), " | Saving FULL sketch-assay integration checkpoint: ",
          out_rds_obj_sketch_integrated)
  saveRDS(obj, out_rds_obj_sketch_integrated, compress = FALSE)
  
  future::plan(old_plan)
  message(Sys.time(), " | Restored future plan.")
  gc()
}

stopifnot("integrated.dr" %in% names(obj@reductions))
stopifnot("SCT" %in% names(obj@assays))

# ------------------------------ Prepare features for ProjectIntegration --------
message(Sys.time(), " | Preparing features for ProjectIntegration (Seurat v5.2.1) ...")

# Variable features live in SCT assay (created by SCTransform on sketch assay).
features_target <- VariableFeatures(obj[["SCT"]])
if (length(features_target) == 0L) stop("VariableFeatures(obj[['SCT']]) is empty.")
message(Sys.time(), " | features_target (from obj[['SCT']] VariableFeatures) = ", length(features_target))

# Use RNA layers matching counts.* (layered counts per sample)
proj_layers <- Layers(obj[["RNA"]], search = "counts")
if (length(proj_layers) == 0L) stop("No RNA layers matched search='counts'.")
message(Sys.time(), " | ProjectIntegration will use RNA layers matching 'counts': n=",
        length(proj_layers), " | example: ", paste(head(proj_layers, 5), collapse = ", "))

# Ensure those same layers exist in sketch assay
sketched_layers <- proj_layers
missing_sk <- setdiff(sketched_layers, Layers(obj[["sketch"]]))
if (length(missing_sk) > 0L) {
  stop("Some RNA count layers are missing in sketch assay. Example: ",
       paste(head(missing_sk, 5), collapse = ", "))
}

# Build projection features that exist in the sketch assay layers (atomic)
features_atom <- Reduce(
  intersect,
  lapply(sketched_layers, function(lyr) Features(obj[["sketch"]], layer = lyr))
)

features.proj <- intersect(features_target, features_atom)
message(Sys.time(), " | features_atom (present in sketch assay layers) = ", length(features_atom))
message(Sys.time(), " | features.proj (used for projection) = ", length(features.proj))
if (length(features.proj) < 500L) stop("Projection feature set too small.")

# Important: set variable features on sketch assay to satisfy Seurat internal defaults
VariableFeatures(obj[["sketch"]]) <- features.proj

# ------------------------------ ProjectIntegration to FULL RNA ----------------
if (file.exists(out_rds_full_proj_ckpt)) {
  message(Sys.time(), " | Found FULL projected checkpoint. Loading: ", out_rds_full_proj_ckpt)
  obj <- readRDS(out_rds_full_proj_ckpt)
} else {
  message(Sys.time(), " | ProjectIntegration: unsketch to FULL RNA ...")
  
  old_plan2 <- future::plan()
  future::plan(sequential)
  
  obj <- ProjectIntegration(
    object = obj,
    sketched.assay = "sketch",
    assay = "RNA",
    reduction = "integrated.dr",
    reduction.name = "integrated.rpca",
    layers = "counts",
    sketched.layers = sketched_layers,
    features = features.proj,
    verbose = TRUE
  )
  gc()
  
  future::plan(old_plan2)
  message(Sys.time(), " | Restored future plan after ProjectIntegration.")
  gc()
  
  message(Sys.time(), " | Saving FULL projected checkpoint: ", out_rds_full_proj_ckpt)
  saveRDS(obj, out_rds_full_proj_ckpt, compress = FALSE)
}

# ------------------------------ Fix embedding colnames (cosmetic) --------------
# Ensure the projected embeddings have stable column names (helps downstream consistency).
if (reduction_rpca %in% names(obj@reductions)) {
  emb <- Embeddings(obj[[reduction_rpca]])
  if (is.null(colnames(emb)) || any(colnames(emb) == "")) {
    colnames(emb) <- paste0("integratedrpca_", seq_len(ncol(emb)))
    obj[[reduction_rpca]]@cell.embeddings <- emb
  }
  rm(emb)
}

message(Sys.time(), " | object information before clusters/UMAP")
print(obj)

# ==============================================================================
# FULL reproducible Neighbors + Clusters + UMAP
# This section is designed to match your standalone reproducible script EXACTLY.
# ==============================================================================
message(Sys.time(), " | Neighbors/Clusters/UMAP on FULL using ", reduction_rpca, " ...")

# Use a stable assay for graph storage/name conventions (does NOT affect reduction used)
DefaultAssay(obj) <- "RNA"

# Force clean slate so reruns are identical and don't reuse stale graphs/UMAP
if (graph_nn %in% names(obj@graphs))  obj@graphs[[graph_nn]] <- NULL
if (graph_snn %in% names(obj@graphs)) obj@graphs[[graph_snn]] <- NULL
if ("seurat_clusters" %in% colnames(obj[[]])) obj$seurat_clusters <- NULL
if (umap_name %in% names(obj@reductions)) obj@reductions[[umap_name]] <- NULL

# Also clear Idents (optional, avoids confusion)
Idents(obj) <- factor(rep("NA", ncol(obj)))

gc()

# IMPORTANT: force sequential execution here (prevents future.apply RNG warnings)
old_plan3 <- future::plan()
future::plan(sequential)

# Deterministic neighbors (FindNeighbors generally not stochastic)
message(Sys.time(), " | FindNeighbors: reduction=", reduction_rpca,
        " dims=", min(dims_use), ":", max(dims_use),
        " graph.name=", graph_nn, "/", graph_snn)

obj <- FindNeighbors(
  obj,
  reduction = reduction_rpca,
  dims = dims_use,
  graph.name = c(graph_nn, graph_snn),
  verbose = FALSE
)

message(Sys.time(), " | Graphs present: ", paste(names(obj@graphs), collapse = ", "))

# Deterministic clustering
message(Sys.time(), " | FindClusters: graph.name=", graph_snn,
        " resolution=", resolution_use,
        " seed=", seed_global)

obj <- FindClusters(
  obj,
  graph.name = graph_snn,
  resolution = resolution_use,
  algorithm = 1,          # Louvain
  n.start = 1,            # reduce nondeterminism
  n.iter = 10,
  random.seed = seed_global,
  verbose = FALSE
)

message(Sys.time(), " | Cluster count: ", length(unique(obj$seurat_clusters)))

# Deterministic UMAP: single-thread uwot + fixed seed
message(Sys.time(), " | RunUMAP: reduction=", reduction_rpca,
        " dims=", min(dims_use), ":", max(dims_use),
        " reduction.name=", umap_name,
        " seed=", seed_umap)

obj <- RunUMAP(
  obj,
  reduction = reduction_rpca,
  dims = dims_use,
  reduction.name = umap_name,
  seed.use = seed_umap,
  umap.method = "uwot",
  metric = "cosine",
  n.threads = 1,
  verbose = FALSE
)

gc()

# restore prior plan
future::plan(old_plan3)

# ------------------------------ Write diagnostics ------------------------------
message(Sys.time(), " | Writing outputs to: ", out_dir)

# UMAP by cluster
p1 <- DimPlot(
  obj,
  reduction = umap_name,
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) + ggtitle(paste0("FULL UMAP (", umap_name, ") by seurat_clusters"))

ggsave(out_umap_by_cluster, p1, width = 10, height = 7, units = "in")
message(Sys.time(), " | Saved: ", out_umap_by_cluster)

# UMAP by sample
p2 <- DimPlot(
  obj,
  reduction = umap_name,
  group.by = split_var,
  raster = TRUE
) + ggtitle(paste0("FULL UMAP (", umap_name, ") by ", split_var))

ggsave(out_umap_by_sample, p2, width = 10, height = 7, units = "in")
message(Sys.time(), " | Saved: ", out_umap_by_sample)

# Cluster counts overall
cluster_counts <- as.data.frame(table(obj$seurat_clusters), stringsAsFactors = FALSE)
colnames(cluster_counts) <- c("cluster", "n_cells")
suppressWarnings({
  if (all(!is.na(as.integer(cluster_counts$cluster)))) {
    cluster_counts <- cluster_counts[order(as.integer(cluster_counts$cluster)), , drop = FALSE]
  } else {
    cluster_counts <- cluster_counts[order(cluster_counts$n_cells, decreasing = TRUE), , drop = FALSE]
  }
})
write.csv(cluster_counts, out_counts_overall, row.names = FALSE)
message(Sys.time(), " | Saved: ", out_counts_overall)

# Cluster counts by sample
tab_by_sample <- as.data.frame(
  table(obj[[split_var]][, 1], obj$seurat_clusters),
  stringsAsFactors = FALSE
)
colnames(tab_by_sample) <- c(split_var, "cluster", "n_cells")
write.csv(tab_by_sample, out_counts_by_sample, row.names = FALSE)
message(Sys.time(), " | Saved: ", out_counts_by_sample)

# Summary TXT
sample_counts <- table(obj[[split_var]][, 1])
cluster_sizes <- table(obj$seurat_clusters)

sink(out_summary_txt)
cat("FULL OBJECT SUMMARY (REPRO CLUSTERS + UMAP)\n")
cat("=========================================\n")
cat("Project: ", projectID, "\n")
cat("Seurat: ", as.character(packageVersion("Seurat")), "\n")
cat("Input RDS: ", in_rds, "\n")
cat("Integration checkpoint: ", out_rds_obj_sketch_integrated, "\n")
cat("Projection checkpoint: ", out_rds_full_proj_ckpt, "\n")
cat("Reduction used: ", reduction_rpca, "\n")
cat("Dims: ", min(dims_use), ":", max(dims_use), "\n")
cat("Resolution: ", resolution_use, "\n")
cat("Graph names: ", graph_nn, " / ", graph_snn, "\n")
cat("UMAP reduction: ", umap_name, "\n")
cat("Seeds: global=", seed_global, " umap=", seed_umap, "\n\n")

cat("Total cells: ", ncol(obj), "\n")
cat("Unique samples: ", length(unique(obj[[split_var]][, 1])), "\n")
cat("Clusters: ", length(unique(obj$seurat_clusters)), "\n\n")

cat("Cells per sample (min/median/mean/max):\n")
cat("  min   = ", min(sample_counts), "\n")
cat("  median= ", median(as.numeric(sample_counts)), "\n")
cat("  mean  = ", mean(as.numeric(sample_counts)), "\n")
cat("  max   = ", max(sample_counts), "\n\n")

cat("Cluster sizes (min/median/mean/max):\n")
cat("  min   = ", min(cluster_sizes), "\n")
cat("  median= ", median(as.numeric(cluster_sizes)), "\n")
cat("  mean  = ", mean(as.numeric(cluster_sizes)), "\n")
cat("  max   = ", max(cluster_sizes), "\n\n")

cat("Top 20 largest clusters:\n")
print(head(sort(cluster_sizes, decreasing = TRUE), 20))

cat("\nTop 20 smallest clusters:\n")
print(head(sort(cluster_sizes, decreasing = FALSE), 20))
sink()

message(Sys.time(), " | Saved: ", out_summary_txt)

# ------------------------------ Save final object ------------------------------
message(Sys.time(), " | Saving FULL object with reproducible clusters/UMAP: ", out_rds_full_proj_final)
saveRDS(obj, out_rds_full_proj_final, compress = FALSE)

message(Sys.time(), " | Done.")
