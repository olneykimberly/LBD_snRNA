#!/usr/bin/env Rscript
# =============================================================================
# FINAL FULL object: reproducible Neighbors + Clusters + UMAP
# - Designed for SLURM (8 cores, 250GB)
# - Forces sequential execution (no future RNG warnings)
# - Explicit seeds + explicit graph names
# - Saves: updated FULL object + PDFs + CSVs + TXT summary
# =============================================================================
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

# ----------------------------- Paths / IDs -----------------------------------
projectID <- "CWOW_cellbender"

in_full_rds <- "../rObjects/CWOW_cellbender_FULL_projected_integratedrpca.rds"

out_dir <- "../results/final_full_rpca_diagnostics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save updated object (with reproducible clusters/UMAP)
out_full_rds <- file.path(dirname(in_full_rds),
                          paste0(projectID, "_FULL_projected_integratedrpca_reproClustersUMAP.rds"))

# ----------------------------- Parameters ------------------------------------
split_var <- "Sample_ID"

reduction_rpca <- "integrated.rpca"   # integrated embedding for ALL cells
umap_name <- "umap.rpca"

dims_use <- 1:30
resolution_use <- 0.3

# explicit graph names (prevents missing-graph errors)
graph_nn  <- "integratedrpca_nn"
graph_snn <- "integratedrpca_snn"

# reproducibility
seed_global <- 12345
seed_umap   <- 4242

# ----------------------------- Outputs ---------------------------------------
out_umap_by_cluster <- file.path(out_dir, paste0(projectID, "_FULL_UMAP_by_cluster.pdf"))
out_umap_by_sample  <- file.path(out_dir, paste0(projectID, "_FULL_UMAP_by_", split_var, ".pdf"))
out_counts_overall  <- file.path(out_dir, paste0(projectID, "_FULL_cluster_counts_overall.csv"))
out_counts_by_sample <- file.path(out_dir, paste0(projectID, "_FULL_cluster_counts_by_sample.csv"))
out_summary_txt     <- file.path(out_dir, paste0(projectID, "_FULL_cluster_summary.txt"))

# ----------------------------- CPU / threads ---------------------------------
# SLURM provides SLURM_CPUS_PER_TASK. Use it to cap BLAS/OpenMP thread usage.
cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
cpus <- ifelse(is.na(cpus) || cpus < 1, 8, cpus)

Sys.setenv(
  OMP_NUM_THREADS = as.character(cpus),
  OPENBLAS_NUM_THREADS = as.character(cpus),
  MKL_NUM_THREADS = as.character(cpus),
  VECLIB_MAXIMUM_THREADS = as.character(cpus),
  NUMEXPR_NUM_THREADS = as.character(cpus)
)

# Make randomness deterministic across the entire run
set.seed(seed_global)

message(Sys.time(), " | SLURM_CPUS_PER_TASK=", cpus)
message(Sys.time(), " | Loading FULL object: ", in_full_rds)

obj <- readRDS(in_full_rds)
stopifnot(inherits(obj, "Seurat"))

if (!(split_var %in% colnames(obj[[]]))) {
  stop("'", split_var, "' not found in meta.data. First columns: ",
       paste(head(colnames(obj[[]]), 50), collapse = ", "))
}

if (!(reduction_rpca %in% names(obj@reductions))) {
  stop("Missing reduction '", reduction_rpca, "'. Available: ",
       paste(names(obj@reductions), collapse = ", "))
}

message(Sys.time(), " | Cells=", ncol(obj), " | Features=", nrow(obj),
        " | Unique ", split_var, "=", length(unique(obj[[split_var]][, 1])))

# Use a stable assay for graphs/metadata storage (does NOT affect reduction used)
DefaultAssay(obj) <- "RNA"

# -----------------------------------------------------------------------------
# Force clean slate so reruns are identical and don't reuse stale graphs/UMAP
# -----------------------------------------------------------------------------
if (graph_nn %in% names(obj@graphs))  obj@graphs[[graph_nn]] <- NULL
if (graph_snn %in% names(obj@graphs)) obj@graphs[[graph_snn]] <- NULL
if ("seurat_clusters" %in% colnames(obj[[]])) obj$seurat_clusters <- NULL
if (umap_name %in% names(obj@reductions)) obj@reductions[[umap_name]] <- NULL

# Also clear Idents (optional, avoids confusion)
Idents(obj) <- factor(rep("NA", ncol(obj)))

gc()

# -----------------------------------------------------------------------------
# Neighbors (deterministic)
# -----------------------------------------------------------------------------
message(Sys.time(), " | FindNeighbors: reduction=", reduction_rpca,
        " dims=", min(dims_use), ":", max(dims_use),
        " graph.name=", graph_nn, "/", graph_snn)

# Determinism: FindNeighbors isn't supposed to be stochastic, but we keep seeds set.
obj <- FindNeighbors(
  obj,
  reduction = reduction_rpca,
  dims = dims_use,
  graph.name = c(graph_nn, graph_snn),
  verbose = FALSE
)

message(Sys.time(), " | Graphs present: ", paste(names(obj@graphs), collapse = ", "))

# -----------------------------------------------------------------------------
# Clustering (deterministic)
# -----------------------------------------------------------------------------
message(Sys.time(), " | FindClusters: graph.name=", graph_snn,
        " resolution=", resolution_use,
        " seed=", seed_global)

# Seurat's clustering uses randomness for some algorithms; seed.use makes it explicit.
obj <- FindClusters(
  obj,
  graph.name = graph_snn,
  resolution = resolution_use,
  algorithm = 1,          # Louvain (stable + common); change only if you prefer Leiden
  n.start = 1,            # reduce nondeterminism
  n.iter = 10,
  random.seed = seed_global,
  verbose = FALSE
)

message(Sys.time(), " | Cluster count: ", length(unique(obj$seurat_clusters)))

# -----------------------------------------------------------------------------
# UMAP (deterministic)
# -----------------------------------------------------------------------------
message(Sys.time(), " | RunUMAP: reduction=", reduction_rpca,
        " dims=", min(dims_use), ":", max(dims_use),
        " reduction.name=", umap_name,
        " seed=", seed_umap)

# Deterministic UMAP settings:
# - uwot, fixed seed.use, and force single-thread for uwot to avoid nondeterminism.
#   (uwot can be parallel; parallel can introduce slight nondeterminism.)
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

# -----------------------------------------------------------------------------
# Write diagnostics
# -----------------------------------------------------------------------------
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
cat("Input RDS: ", in_full_rds, "\n")
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

# -----------------------------------------------------------------------------
# Save updated object
# -----------------------------------------------------------------------------
message(Sys.time(), " | Saving FULL object with reproducible clusters/UMAP: ", out_full_rds)
saveRDS(obj, out_full_rds, compress = FALSE)

message(Sys.time(), " | Done.")
