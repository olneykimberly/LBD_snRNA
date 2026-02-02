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

# ----------------------------- Inputs -----------------------------------------
project_ID <- "CWOW_cellbender"
in_rds <- paste0("../rObjects/", project_ID, "_annotated_with_celltype_subclusters.rds")

obj <- readRDS(in_rds)
message("Loaded object:")
print(obj)

stopifnot("cell_type" %in% colnames(obj[[]]))

# Where your per-cell-type subclustering outputs live
out_base <- "../results/subclustering_cell_type"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- Knobs ------------------------------------------
min_target_q <- 0.10   # drop bottom 10% identity score within each cell type
max_contam_q <- 0.95   # drop top 5% contamination within each cell type

min_cells_ct <- 200    # skip tiny cell types

# ----------------------------- Marker gene sets (edit freely) -----------------
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

# ----------------------------- Helpers ----------------------------------------
target_col <- function(ct) paste0("Score_", ct, "_1")
contam_cols_all <- grep("^Contam_.*_1$", colnames(obj[[]]), value = TRUE)

safe_quant <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

# ----------------------------- Compute max_contam globally --------------------
# Treat missing contamination scores as -Inf (so they do not inflate max_contam)
if (length(contam_cols_all) == 0L) {
  warning("No contamination columns (Contam_*_1) found in meta.data. max_contam will be NA.")
  obj$max_contam <- NA_real_
} else {
  contam_mat <- obj[[]][, contam_cols_all, drop = FALSE]
  contam_mat <- as.matrix(contam_mat)
  
  # treat NAs as -Inf so they don't inflate max_contam
  contam_mat[is.na(contam_mat)] <- -Inf
  
  obj$max_contam <- apply(contam_mat, 1, max)
  rm(contam_mat)
}

# Prepare columns to store thresholds + keep/drop labels
obj$thr_target <- NA_real_
obj$thr_contam <- NA_real_
obj$score_filter <- NA_character_

cells_keep <- rep(TRUE, ncol(obj))
names(cells_keep) <- colnames(obj)

summary_rows <- list()

# Use RNA for gene expression violins
DefaultAssay(obj) <- "RNA"

# ----------------------------- Main per-cell-type loop ------------------------
cts <- sort(unique(as.character(obj$cell_type)))
cts <- cts[!is.na(cts) & cts != ""]

for (ct in cts) {
  ct_cells <- colnames(obj)[which(obj$cell_type == ct)]
  if (length(ct_cells) < min_cells_ct) {
    message("Skipping ", ct, " (too few cells: ", length(ct_cells), ")")
    next
  }
  
  tc <- target_col(ct)
  if (!tc %in% colnames(obj[[]])) {
    message("Skipping target filter for ", ct, " (missing ", tc, ")")
    # Still label (no thresholds)
    obj$score_filter[ct_cells] <- "KEEP"  # default keep if no score
    summary_rows[[ct]] <- data.frame(
      cell_type = ct,
      n_cells = length(ct_cells),
      n_keep = length(ct_cells),
      n_drop = 0,
      keep_frac = 1,
      thr_target = NA_real_,
      thr_contam = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }
  
  # Identity score threshold (within cell type)
  target_vals <- obj[[tc]][ct_cells, 1]
  thr_target <- safe_quant(target_vals, min_target_q)
  
  # Contamination threshold (within cell type)
  max_contam_vals <- obj$max_contam[ct_cells]
  thr_contam <- safe_quant(max_contam_vals, max_contam_q)
  
  keep_ct <- rep(TRUE, length(ct_cells))
  names(keep_ct) <- ct_cells
  
  if (!is.na(thr_target)) keep_ct <- keep_ct & (target_vals >= thr_target)
  if (!is.na(thr_contam)) keep_ct <- keep_ct & (max_contam_vals <= thr_contam)
  
  cells_keep[ct_cells] <- keep_ct
  
  # Write thresholds + KEEP/DROP labels into obj metadata
  obj$thr_target[ct_cells] <- thr_target
  obj$thr_contam[ct_cells] <- thr_contam
  obj$score_filter[ct_cells] <- ifelse(keep_ct, "KEEP", "DROP")
  
  # Summaries
  n_keep <- sum(keep_ct, na.rm = TRUE)
  n_drop <- length(ct_cells) - n_keep
  
  summary_rows[[ct]] <- data.frame(
    cell_type = ct,
    n_cells = length(ct_cells),
    n_keep = n_keep,
    n_drop = n_drop,
    keep_frac = n_keep / length(ct_cells),
    thr_target = thr_target,
    thr_contam = thr_contam,
    stringsAsFactors = FALSE
  )
  
  message(
    ct, ": keep ", n_keep, "/", length(ct_cells),
    " | target>=Q", min_target_q, " (", round(thr_target, 3), ")",
    " | maxContam<=Q", max_contam_q, " (", round(thr_contam, 3), ")"
  )
  
  # --------------------------- Plot outputs per cell type ---------------------
  out_dir <- file.path(out_base, ct)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ct_obj <- subset(obj, cells = ct_cells)
  
  # Make score_filter an ordered factor for nicer violins
  ct_obj$score_filter <- factor(ct_obj$score_filter, levels = c("KEEP", "DROP"))
  
  # 1) Target score violin
  p_target <- VlnPlot(
    ct_obj,
    features = tc,
    group.by = "score_filter",
    pt.size = 0.01
  ) + ggtitle(paste0(ct, " | ", tc, " | KEEP vs DROP"))
  
  ggsave(
    filename = file.path(out_dir, paste0(project_ID, "_", ct, "_Vln_targetScore.pdf")),
    plot = p_target, width = 7, height = 5, units = "in"
  )
  
  # 2) Max contamination violin
  p_contam <- VlnPlot(
    ct_obj,
    features = "max_contam",
    group.by = "score_filter",
    pt.size = 0.01
  ) + ggtitle(paste0(ct, " | max_contam | KEEP vs DROP"))
  
  ggsave(
    filename = file.path(out_dir, paste0(project_ID, "_", ct, "_Vln_maxContam.pdf")),
    plot = p_contam, width = 7, height = 5, units = "in"
  )
  
  # 3) Marker genes for this cell type
  genes_target <- markers_list[[ct]]
  genes_target <- intersect(genes_target, rownames(ct_obj))
  
  if (length(genes_target) > 0) {
    p_markers <- VlnPlot(
      ct_obj,
      features = genes_target,
      group.by = "score_filter",
      pt.size = 0
    ) + ggtitle(paste0(ct, " | marker genes | KEEP vs DROP"))
    
    ggsave(
      filename = file.path(out_dir, paste0(project_ID, "_", ct, "_Vln_markers.pdf")),
      plot = p_markers, width = 12, height = max(5, 1.2 * length(genes_target)), units = "in"
    )
  }
  
  # 4) Contaminant marker genes (all other lineages)
  contam_types <- contam_map[[ct]]
  contam_genes <- unique(unlist(markers_list[contam_types], use.names = FALSE))
  contam_genes <- intersect(contam_genes, rownames(ct_obj))
  
  if (length(contam_genes) > 0) {
    p_contam_genes <- VlnPlot(
      ct_obj,
      features = contam_genes,
      group.by = "score_filter",
      pt.size = 0
    ) + ggtitle(paste0(ct, " | contaminant marker genes | KEEP vs DROP"))
    
    ggsave(
      filename = file.path(out_dir, paste0(project_ID, "_", ct, "_Vln_contam_markers.pdf")),
      plot = p_contam_genes, width = 14, height = max(5, 0.9 * length(contam_genes)), units = "in"
    )
  }
  
  # Optional: write per-cell table for THIS cell type
  md_ct <- ct_obj[[]] %>%
    dplyr::select(any_of(c("Sample_ID", "group", "cell_type", "cell_type_subcluster",
                           "score_filter", "thr_target", "thr_contam", "max_contam", tc))) %>%
    dplyr::mutate(barcode = rownames(ct_obj[[]]))
  
  write.csv(
    md_ct,
    file = file.path(out_dir, paste0(project_ID, "_", ct, "_score_filter_percell.csv")),
    row.names = FALSE
  )
  
  rm(ct_obj)
  gc()
}

# ----------------------------- Save summary tables ----------------------------
summary_df <- bind_rows(summary_rows) %>%
  arrange(factor(cell_type, levels = unique(cell_type)))

summary_out <- file.path(out_base, "score_filter_summary_by_celltype.csv")
write.csv(summary_df, summary_out, row.names = FALSE)
message("Saved summary table: ", summary_out)

# Save an updated object with the new metadata columns (no filtering applied yet)
out_obj <- paste0("../rObjects/", project_ID, "_with_score_filter_labels.rds")
saveRDS(obj, out_obj, compress = FALSE)
message("Saved object with KEEP/DROP labels: ", out_obj)

# ----------------------------- Apply filter (optional) -------------------------
obj_clean_scores <- subset(obj, cells = names(cells_keep)[cells_keep])
out_clean <- paste0("../rObjects/", project_ID, "_clean_by_score_filters.rds")
saveRDS(obj_clean_scores, out_clean, compress = FALSE)
message("Saved filtered object: ", out_clean)

message("DONE.")
