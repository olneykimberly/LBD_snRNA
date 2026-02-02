#!/usr/bin/env Rscript
# =============================================================================
# Per-cell-type subclustering + Harmony Integration + Contamination Scoring
# =============================================================================
# Set working directory and load project-specific colors/paths
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

# ----------------------------- Inputs / Outputs ------------------------------
project_ID <- "CWOW_cellbender"
in_rds <- paste0("../rObjects/", project_ID, "_filtered_subclusters.rds")
out_base <- "../results/subclustering_cell_type_pass2"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Final object with all integrated metadata
out_full_with_subclusters <- paste0("../rObjects/", project_ID, "_annotated_with_celltype_subclusters_pass2_harmony.rds")

# ----------------------------- Parameters ------------------------------------
order_cell_type <- c("neuron", "interneuron", "oligodendrocyte", "opc",
                     "astrocyte", "microglia", "mural", "fibroblast", "endothelial")

dims_use <- 1:30
npcs <- 50
nfeatures <- 3000
resolution_use <- 0.2
seed_global <- 12345
seed_umap <- 4242
umap_nthreads <- 1

# ----------------------------- Marker & Score Sets ---------------------------
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

# ----------------------------- Load & Clean Object ---------------------------
dataObject <- readRDS(in_rds)
DefaultAssay(dataObject) <- "RNA"

# Initialize/Reset metadata columns
dataObject$cell_type_subcluster <- NA_character_
dataObject$cell_type_subcluster_id <- NA_character_

# Remove old scoring columns to avoid duplication
meta <- dataObject[[]]
drop_idx <- grepl("^Contam_\\.", colnames(meta)) | grepl("^Score_\\.", colnames(meta))
dataObject[[]] <- meta[, !drop_idx, drop = FALSE]

# ----------------------------- Main Loop --------------------------------------
set.seed(seed_global)

for (ct_name in order_cell_type) {
  message("\n", "============================================================")
  message(Sys.time(), " | Processing Cell type: ", ct_name)
  message("============================================================")
  
  out_dir <- file.path(out_base, ct_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Subset to current cell type
  ct <- subset(dataObject, subset = cell_type == ct_name)
  if (ncol(ct) < 200) {
    message("Skipping ", ct_name, ": insufficient cell count.")
    next
  }
  
  # 1. Harmony Integration Within Cell Type
  DefaultAssay(ct) <- "RNA"
  ct <- FindVariableFeatures(ct, nfeatures = nfeatures, verbose = FALSE)
  ct <- ScaleData(ct, features = VariableFeatures(ct), verbose = FALSE)
  ct <- RunPCA(ct, features = VariableFeatures(ct), npcs = npcs, verbose = FALSE)
  
  ct <- RunHarmony(ct, group.by.vars = "Sample_ID", reduction = "pca", 
                   assay.use = "RNA", reduction.save = "harmony", verbose = FALSE)
  
  ct <- FindNeighbors(ct, reduction = "harmony", dims = dims_use, 
                      graph.name = c("ct_nn","ct_snn"), verbose = FALSE)
  
  ct <- FindClusters(ct, graph.name = "ct_snn", resolution = resolution_use, 
                     random.seed = seed_global, verbose = FALSE)
  
  ct <- RunUMAP(ct, reduction = "harmony", dims = dims_use, reduction.name = "umap.ct",
                umap.method = "uwot", metric = "cosine", n.threads = umap_nthreads,
                seed.use = seed_umap, verbose = FALSE)
  
  ct$subcluster <- paste0(ct_name, "_", ct$seurat_clusters)
  
  # 2. Module Scoring (Identity & Contamination)
  score_cols <- c()
  # Target Score
  target_genes <- intersect(markers_list[[ct_name]], rownames(ct))
  if (length(target_genes) >= 2) {
    ct <- AddModuleScore(ct, features = list(target_genes), name = paste0("Score_", ct_name, "_"))
    score_cols <- c(score_cols, paste0("Score_", ct_name, "_1"))
  }
  # Contamination Scores
  other_types <- setdiff(names(markers_list), ct_name)
  for (ot in other_types) {
    ot_genes <- intersect(markers_list[[ot]], rownames(ct))
    if (length(ot_genes) >= 2) {
      ct <- AddModuleScore(ct, features = list(ot_genes), name = paste0("Contam_", ot, "_"))
      score_cols <- c(score_cols, paste0("Contam_", ot, "_1"))
    }
  }
  
  # 3. Write Metadata Back to dataObject
  cells_ct <- colnames(ct)
  dataObject$cell_type_subcluster[cells_ct] <- ct$subcluster
  dataObject$cell_type_subcluster_id[cells_ct] <- as.character(ct$seurat_clusters)
  
  for (sc in score_cols) {
    if (!(sc %in% colnames(dataObject[[]]))) { dataObject[[sc]] <- NA_real_ }
    dataObject@meta.data[cells_ct, sc] <- as.numeric(ct@meta.data[cells_ct, sc])
  }
  
  # 4. Visualization & Outputs
  # Diagnostic Side-by-Side UMAP
  common_theme <- theme(axis.title=element_text(size=8), axis.text=element_text(size=8),
                        legend.position="none", plot.title=element_text(size=9))
  
  p1 <- DimPlot(ct, reduction = "umap.ct", group.by = "Sample_ID") + 
    ggtitle(paste(ct_name, "- Harmony by Sample")) + common_theme + scale_color_manual(values = color_panel)
  
  p2 <- DimPlot(ct, reduction = "umap.ct", group.by = "subcluster", label = TRUE, label.size = 3) + 
    ggtitle(paste(ct_name, "- Integrated Clusters")) + common_theme + scale_color_manual(values = color_panel)
  
  pdf(file.path(out_dir, paste0(project_ID, "_", ct_name, "_UMAP_Comparison.pdf")), width = 12, height = 6)
  print(p1 + p2); dev.off()
  
  # DotPlot of markers
  genes_dot <- unique(unlist(markers_list))
  genes_dot <- intersect(genes_dot, rownames(ct))
  pdf(file.path(out_dir, paste0(project_ID, "_", ct_name, "_DotPlot_Markers.pdf")), width = 14, height = 8)
  print(DotPlot(ct, features = genes_dot, assay = "RNA") + RotatedAxis() + ggtitle(paste(ct_name, "Markers"))); dev.off()
  
  # Save per-cell-type RDS
  saveRDS(ct, file.path(out_dir, paste0(project_ID, "_", ct_name, "_subclustered.rds")), compress = FALSE)
  
  rm(ct); gc()
}

# ----------------------------- Final Save -------------------------------------
message(Sys.time(), " | Saving Full Object...")
saveRDS(dataObject, out_full_with_subclusters, compress = FALSE)
message("Pipeline Complete.")