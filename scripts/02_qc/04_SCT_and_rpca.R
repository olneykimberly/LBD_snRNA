#!/usr/bin/env Rscript
# ==============================================================================
# Seurat v5.2.1
# SketchData -> SCTransform+PCA+IntegrateLayers (on FULL obj, assay="sketch") ->
# ProjectIntegration to FULL RNA (layered counts) ->
# FindNeighbors/FindClusters/RunUMAP on projected embedding
#
# Fixes:
#  - Remove future.apply RNG misuse warning by running Neighbors/Clusters/UMAP sequential
#  - Keep heavy integration/projection steps checkpointed
# ==============================================================================

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(ggplot2)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

options(future.globals.maxSize = 500 * 1024^3)

# --------------------------- User parameters ----------------------------------
projectID <- "CWOW_cellbender"
in_rds <- paste0("../rObjects/", projectID, "_singlets_scDblFinder_exprate.rds")

out_dir <- file.path("../results", "integration_rpca_sct_sketch")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Checkpoints
out_rds_obj_sketch_integrated <- paste0("../rObjects/", projectID, "_FULL_sketchAssay_integrateddr.rds")
out_rds_full_proj_ckpt        <- paste0("../rObjects/", projectID, "_FULL_projected_integratedrpca_ckpt.rds")
out_rds_full_proj_final       <- paste0("../rObjects/", projectID, "_FULL_projected_integratedrpca.rds")

split_var <- "Sample_ID"
dims_use <- 1:30
npcs <- 50
resolution_use <- 0.3
nfeatures <- 3000

cells_per_sample_sketch <- 5000
sketch_seed <- 12345

# Parallelism for heavy steps
n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
n_workers <- max(1, n_workers - 1)

if (Sys.info()[["sysname"]] == "Linux" && Sys.getenv("RSTUDIO") == "") {
  plan(multicore, workers = n_workers)
  message(Sys.time(), " | future plan: multicore | workers=", n_workers)
} else {
  plan(multisession, workers = n_workers)
  message(Sys.time(), " | future plan: multisession | workers=", n_workers)
}

set.seed(sketch_seed)
options(future.rng.onMisuse = "warning")

# ------------------------------ Load FULL object ------------------------------
message(Sys.time(), " | Loading Seurat object: ", in_rds)
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

# ------------------------------ SketchData -----------------------------------
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
  future::plan(sequential)
  
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

features_target <- VariableFeatures(obj[["SCT"]])
if (length(features_target) == 0L) {
  stop("VariableFeatures(obj[['SCT']]) is empty.")
}
message(Sys.time(), " | features_target (from obj[['SCT']] VariableFeatures) = ", length(features_target))

proj_layers <- Layers(obj[["RNA"]], search = "counts")
if (length(proj_layers) == 0L) stop("No RNA layers matched search='counts'.")
message(Sys.time(), " | ProjectIntegration will use RNA layers matching 'counts': n=",
        length(proj_layers), " | example: ", paste(head(proj_layers, 5), collapse = ", "))

sketched_layers <- proj_layers
missing_sk <- setdiff(sketched_layers, Layers(obj[["sketch"]]))
if (length(missing_sk) > 0L) stop("Some RNA count layers are missing in sketch assay.")

features_atom <- Reduce(
  intersect,
  lapply(sketched_layers, function(lyr) Features(obj[["sketch"]], layer = lyr))
)

features.proj <- intersect(features_target, features_atom)
message(Sys.time(), " | features_atom (present in sketch assay layers) = ", length(features_atom))
message(Sys.time(), " | features.proj (used for projection) = ", length(features.proj))
if (length(features.proj) < 500L) stop("Projection feature set too small.")

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
if ("integrated.rpca" %in% names(obj@reductions)) {
  emb <- Embeddings(obj[["integrated.rpca"]])
  if (is.null(colnames(emb)) || any(colnames(emb) == "")) {
    colnames(emb) <- paste0("integratedrpca_", seq_len(ncol(emb)))
    obj[["integrated.rpca"]]@cell.embeddings <- emb
  }
  rm(emb)
}

print("object information before clusters/UMAP")
print(obj)

# ------------------------------ Neighbors/Clusters/UMAP ------------------------
# IMPORTANT: force sequential here to remove future.apply RNG misuse warnings
message(Sys.time(), " | Neighbors/Clusters/UMAP on FULL using integrated.rpca ...")

old_plan3 <- future::plan()
future::plan(sequential)

# If anything tries to RNG in parallel, this would now be deterministic anyway
options(future.rng.onMisuse = "error")

DefaultAssay(obj) <- "RNA"

graph_nn  <- "integratedrpca_nn"
graph_snn <- "integratedrpca_snn"

obj <- FindNeighbors(
  obj,
  reduction = "integrated.rpca",
  dims = dims_use,
  graph.name = c(graph_nn, graph_snn),
  verbose = FALSE
)

message(Sys.time(), " | Graphs present: ", paste(names(obj@graphs), collapse = ", "))

obj <- FindClusters(
  obj,
  graph.name = graph_snn,
  resolution = resolution_use,
  verbose = FALSE
)

obj <- RunUMAP(
  obj,
  reduction = "integrated.rpca",
  dims = dims_use,
  reduction.name = "umap.rpca",
  umap.method = "uwot",
  metric = "cosine",
  seed.use = 42,
  verbose = FALSE
)
gc()

# restore prior plan
future::plan(old_plan3)
options(future.rng.onMisuse = "warning")

message(Sys.time(), " | Cluster count: ", length(unique(obj$seurat_clusters)))
message(Sys.time(), " | Cells per cluster (top 10):")
print(head(sort(table(obj$seurat_clusters), decreasing = TRUE), 10))

# ------------------------------ Save final ------------------------------------
message(Sys.time(), " | Saving FULL projected object: ", out_rds_full_proj_final)
saveRDS(obj, out_rds_full_proj_final, compress = FALSE)

message(Sys.time(), " | Done.")
