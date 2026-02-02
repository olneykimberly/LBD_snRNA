#!/usr/bin/env Rscript
# =============================================================================
# Per-cell-type subclustering + contamination scoring + diagnostics
# (OPTION 1: write subcluster labels + scores BACK onto the FULL object)
#
# What this does:
#  1) Loads your annotated FULL object (with cell_type already assigned)
#  2) For each cell type in `order_cell_type`:
#       - subsets to that cell type
#       - runs HVG -> Scale -> PCA -> Neighbors -> Clusters -> UMAP (within-cell-type)
#       - computes module scores for identity + contamination (optional but useful)
#       - saves per-cell-type PDFs/CSVs/RDS for review
#       - writes the resulting subcluster labels + scores back into the FULL object
#  3) Saves a new FULL object that contains per-cell-type subcluster metadata so you
#     can later remove unwanted subclusters directly from the FULL object.
#
# Outputs:
#   ../results/subclustering_cell_type/<celltype>/
#   ../rObjects/<project_ID>_annotated_with_celltype_subclusters.rds
# =============================================================================

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

# ----------------------------- Inputs / Outputs ------------------------------
project_ID <- "CWOW_cellbender"
in_rds <- paste0("../rObjects/", project_ID, "_annotated_before_subcluster.rds")

out_base <- "../results/subclustering_cell_type"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Save FULL object with subclusters written into meta.data
out_full_with_subclusters <- paste0("../rObjects/", project_ID, "_annotated_with_celltype_subclusters.rds")

# ----------------------------- Cell type order -------------------------------
order_cell_type <- c(
  "neuron", "interneuron", "oligodendrocyte", "opc",
  "astrocyte", "microglia", "mural", "fibroblast", "endothelial"
)

# ----------------------------- Parameters ------------------------------------
dims_use <- 1:30
npcs <- 50
nfeatures <- 3000
resolution_use <- 0.2

seed_global <- 12345
seed_umap   <- 4242
umap_nthreads <- 1  # deterministic UMAP

# ----------------------------- Load object -----------------------------------
dataObject <- readRDS(in_rds)
message("Loaded FULL object:")
print(dataObject)

celltype_col <- "cell_type"
if (!(celltype_col %in% colnames(dataObject[[]]))) {
  stop("Column '", celltype_col, "' not found in meta.data. Available columns include: ",
       paste(head(colnames(dataObject[[]]), 50), collapse = ", "))
}
message(Sys.time(), " | Using cell type column: ", celltype_col)

# Ensure we're using RNA for marker visualization / module scoring on log-normalized data
DefaultAssay(dataObject) <- "RNA"

# ----------------------------- Marker sets -----------------------------------
markers_list <- list(
  neuron          = c("SNAP25","SYT1","RBFOX3","SLC17A7"),
  interneuron     = c("GAD1","GAD2","SLC6A1","DLX1","DLX2"),
  oligodendrocyte = c("PLP1","MBP","MOG","MAG"),
  opc             = c("PDGFRA","CSPG4","VCAN","PTPRZ1"),
  astrocyte       = c("AQP4","SLC1A3","ALDH1L1","GFAP"),
  microglia       = c("P2RY12","C1QA","TYROBP","SPI1"),
  endothelial     = c("FLT1","KDR","PECAM1","VWF"),
  mural           = c("RGS5","PDGFRB","MCAM","CSPG4"),
  fibroblast      = c("DCN","COL1A1","COL1A2","LUM")
)

score_sets <- list(
  neuron          = c("SNAP25","SYT1","RBFOX3","SLC17A7"),
  interneuron     = c("GAD1","GAD2","SLC6A1","DLX1"),
  oligodendrocyte = c("PLP1","MBP","MOG","MAG"),
  opc             = c("PDGFRA","CSPG4","VCAN","PTPRZ1"),
  astrocyte       = c("AQP4","SLC1A3","ALDH1L1","GFAP"),
  microglia       = c("P2RY12","C1QA","TYROBP","SPI1"),
  endothelial     = c("FLT1","KDR","PECAM1","VWF"),
  mural           = c("RGS5","PDGFRB","MCAM","RERGL"),
  fibroblast      = c("DCN","COL1A1","COL1A2","LUM")
)

contam_map <- list(
  neuron          = c("interneuron","oligodendrocyte","opc","astrocyte","microglia","endothelial","mural","fibroblast"),
  interneuron     = c("neuron","oligodendrocyte","opc","astrocyte","microglia","endothelial","mural","fibroblast"),
  oligodendrocyte = c("opc","astrocyte","microglia","neuron","interneuron","endothelial","mural","fibroblast"),
  opc             = c("oligodendrocyte","astrocyte","microglia","neuron","interneuron","endothelial","mural","fibroblast"),
  astrocyte       = c("microglia","oligodendrocyte","opc","neuron","interneuron","endothelial","mural","fibroblast"),
  microglia       = c("astrocyte","oligodendrocyte","opc","neuron","interneuron","endothelial","mural","fibroblast"),
  endothelial     = c("mural","fibroblast","astrocyte","microglia","oligodendrocyte","opc","neuron","interneuron"),
  mural           = c("endothelial","fibroblast","astrocyte","microglia","oligodendrocyte","opc","neuron","interneuron"),
  fibroblast      = c("endothelial","mural","astrocyte","microglia","oligodendrocyte","opc","neuron","interneuron")
)

# ----------------------------- OPTION 1: initialize FULL columns --------------
# These will be filled for cells belonging to each cell type.
dataObject$cell_type_subcluster     <- NA_character_
dataObject$cell_type_subcluster_id  <- NA_character_

# Score columns are created on-demand (only once) as we discover them
# (we'll add them to FULL object as numeric columns when computed)
# Example names created by AddModuleScore will end with "_1"
#   Score_<celltype>_1
#   Contam_<othertype>_1

# ----------------------------- Main loop --------------------------------------
set.seed(seed_global)

for (ct_name in order_cell_type) {
  message("\n", "============================================================")
  message(Sys.time(), " | Cell type: ", ct_name)
  message("============================================================")
  
  out_dir <- file.path(out_base, ct_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Subset to this cell type
  ct <- subset(dataObject, subset = cell_type == ct_name)
  if (ncol(ct) < 200) {
    message(Sys.time(), " | Skipping ", ct_name, " (too few cells: ", ncol(ct), ")")
    next
  }
  message(Sys.time(), " | Cells in ", ct_name, ": ", ncol(ct))
  
  DefaultAssay(ct) <- "RNA"
  
  # ------------------- Preprocess within cell type ---------------------------
  set.seed(seed_global)
  
  ct <- FindVariableFeatures(ct, nfeatures = nfeatures, verbose = FALSE)
  ct <- ScaleData(ct, features = VariableFeatures(ct), verbose = FALSE)
  ct <- RunPCA(ct, features = VariableFeatures(ct), npcs = npcs, verbose = FALSE)
  
  ct <- FindNeighbors(ct, reduction = "pca", dims = dims_use,
                      graph.name = c("ct_nn","ct_snn"), verbose = FALSE)
  
  ct <- FindClusters(ct, graph.name = "ct_snn",
                     resolution = resolution_use,
                     algorithm = 1, n.start = 1, n.iter = 10,
                     random.seed = seed_global, verbose = FALSE)
  
  ct <- RunUMAP(ct, reduction = "pca", dims = dims_use,
                reduction.name = "umap.ct",
                umap.method = "uwot", metric = "cosine",
                n.threads = umap_nthreads,
                seed.use = seed_umap, verbose = FALSE)
  
  # Stable per-cell-type subcluster labels
  ct$subcluster <- paste0(ct_name, "_", ct$seurat_clusters)
  
  # ------------------- Module scores (target + contamination) ----------------
  score_cols <- c()
  
  # Target score
  target_genes <- score_sets[[ct_name]]
  if (!is.null(target_genes)) {
    target_genes <- intersect(target_genes, rownames(ct))
    if (length(target_genes) >= 2) {
      ct <- AddModuleScore(ct, features = list(target_genes), name = paste0("Score_", ct_name, "_"))
      sc <- paste0("Score_", ct_name, "_1")
      if (sc %in% colnames(ct[[]])) score_cols <- c(score_cols, sc)
    }
  }
  
  # Contamination scores
  contam_types <- contam_map[[ct_name]]
  for (cname in contam_types) {
    genes <- score_sets[[cname]]
    if (is.null(genes)) next
    genes <- intersect(genes, rownames(ct))
    if (length(genes) < 2) next
    ct <- AddModuleScore(ct, features = list(genes), name = paste0("Contam_", cname, "_"))
    sc <- paste0("Contam_", cname, "_1")
    if (sc %in% colnames(ct[[]])) score_cols <- c(score_cols, sc)
  }
  
  # ------------------- OPTION 1: write subclusters + scores back to FULL ------
  cells_ct <- colnames(ct)
  
  # Subcluster labels
  dataObject$cell_type_subcluster[cells_ct]    <- ct$subcluster[cells_ct]
  dataObject$cell_type_subcluster_id[cells_ct] <- as.character(ct$seurat_clusters[cells_ct])
  
  # Score columns (create on FULL if missing, then fill values for these cells)
  if (length(score_cols) > 0) {
    for (sc in score_cols) {
      if (!(sc %in% colnames(dataObject[[]]))) {
        dataObject[[sc]] <- NA_real_
      }
      dataObject@meta.data[cells_ct, sc] <- as.numeric(ct@meta.data[cells_ct, sc])
    }
  }
  
  # ------------------- Outputs: UMAPs + DotPlot + tables ----------------------
  message(Sys.time(), " | Writing outputs: ", out_dir)
  
  # UMAP by subcluster
  p_umap <- DimPlot(ct, reduction = "umap.ct", group.by = "subcluster",
                    label = TRUE, repel = TRUE, raster = TRUE) +
    ggtitle(paste0(ct_name, ": UMAP by subcluster (RNA PCA)"))
  ggsave(file.path(out_dir, paste0(project_ID, "_", ct_name, "_UMAP_by_subcluster.pdf")),
         p_umap, width = 11, height = 7, units = "in")
  
  # UMAP colored by scores
  if (length(score_cols) > 0) {
    p_scores <- FeaturePlot(ct, reduction = "umap.ct", features = score_cols,
                            raster = TRUE, keep.scale = "feature") +
      ggtitle(paste0(ct_name, ": module scores"))
    ggsave(file.path(out_dir, paste0(project_ID, "_", ct_name, "_UMAP_scores.pdf")),
           p_scores, width = 14, height = 9, units = "in")
  }
  
  # DotPlot of canonical + contamination markers
  genes_dot <- unique(c(
    markers_list[[ct_name]],
    unlist(markers_list[contam_types], use.names = FALSE)
  ))
  genes_dot <- genes_dot[!is.na(genes_dot)]
  genes_dot <- intersect(genes_dot, rownames(ct))
  
  if (length(genes_dot) > 0) {
    Idents(ct) <- ct$subcluster
    p_dot <- DotPlot(ct, features = genes_dot, assay = "RNA") +
      RotatedAxis() +
      ggtitle(paste0(ct_name, ": DotPlot (canonical + contamination markers)"))
    ggsave(file.path(out_dir, paste0(project_ID, "_", ct_name, "_DotPlot_markers.pdf")),
           p_dot, width = 14, height = 6 + 0.08 * length(unique(ct$subcluster)), units = "in")
  }
  
  # Per-cell table with scores + subcluster
  keep_cols <- c(celltype_col, "Sample_ID", "group", "seurat_clusters", "subcluster", score_cols)
  keep_cols <- keep_cols[keep_cols %in% colnames(ct[[]])]
  
  md <- ct[[]][, keep_cols, drop = FALSE]
  md$barcode <- rownames(md)
  write.csv(md, file.path(out_dir, paste0(project_ID, "_", ct_name, "_subcluster_scores_percell.csv")),
            row.names = FALSE)
  
  # Save per-cell-type object for interactive inspection
  saveRDS(ct, file.path(out_dir, paste0(project_ID, "_", ct_name, "_subclustered.rds")),
          compress = FALSE)
  
  message(Sys.time(), " | Done cell type: ", ct_name,
          " | cells=", ncol(ct),
          " | subclusters=", length(unique(ct$subcluster)))
  rm(ct)
  gc()
}

# ----------------------------- Save FULL object with subclusters --------------
message("\n", Sys.time(), " | Writing FULL object with per-cell-type subclusters to: ", out_full_with_subclusters)
saveRDS(dataObject, out_full_with_subclusters, compress = FALSE)

message("\n", Sys.time(), " | ALL DONE. Outputs in: ", out_base)
message(Sys.time(), " | FULL object saved: ", out_full_with_subclusters)
