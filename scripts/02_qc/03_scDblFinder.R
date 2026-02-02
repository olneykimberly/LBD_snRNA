#!/usr/bin/env Rscript

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

options(future.globals.maxSize = 250 * 1024^3)
projectID <- "CWOW_cellbender"
out_dir <- file.path("../results", "scDblFinder_exprate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(Sys.time(), " | Loading filtered Seurat object...")
dataObject <- readRDS(file = paste0("../rObjects/", projectID, "_filtered.rds"))

# -------------------------------------------------------------------------
# Build verified technical-to-descriptive map
# -------------------------------------------------------------------------
all_layers <- Layers(dataObject, assay = "RNA", search = "counts")

mapping_table <- dataObject@meta.data %>%
  select(Sample_BR, Sample_ID) %>%
  distinct() %>%
  mutate(layer_name = paste0("counts.", Sample_BR)) %>%
  filter(layer_name %in% all_layers)

if (nrow(mapping_table) == 0) {
  stop("No samples found after layer filtering. Check your layer naming (counts.<Sample_BR>).")
}

# -------------------------------------------------------------------------
# Read Cell Ranger metrics to get Estimated Number of Cells per sample
# Path pattern: /.../cellranger/<sample_br>/outs/metrics_summary.csv
# -------------------------------------------------------------------------
metrics_root <- "../cellranger"

read_est_cells <- function(sample_br, root_dir) {
  f <- file.path(root_dir, sample_br, "outs", "metrics_summary.csv")
  if (!file.exists(f)) {
    stop("Missing metrics_summary.csv for ", sample_br, " at: ", f)
  }
  
  # metrics_summary.csv is usually 2-row wide table: header row + 1 values row
  m <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
  
  # Robust lookup: column name can vary slightly across cellranger versions
  # Prefer exact, then fallback to partial match
  if ("Estimated Number of Cells" %in% names(m)) {
    x <- m[["Estimated Number of Cells"]][1]
  } else {
    idx <- grep("Estimated.*Number.*Cells|Estimated.*Cells", names(m), ignore.case = TRUE)
    if (length(idx) == 0) stop("Could not find Estimated Number of Cells column in: ", f)
    x <- m[[ idx[1] ]][1]
  }
  
  # Value is often like "13,953" as a string -> remove commas
  x <- as.numeric(gsub(",", "", as.character(x)))
  if (is.na(x) || x <= 0) stop("Bad Estimated Cells value in ", f, ": ", x)
  return(x)
}

# Add estimated cells + dbr_use to mapping table
mapping_table <- mapping_table %>%
  rowwise() %>%
  mutate(
    estimated_cells = read_est_cells(Sample_BR, metrics_root),
    # your linear heuristic, but anchored to the CellRanger metric
    dbr_use = (estimated_cells / 1000) * 0.008
  ) %>%
  ungroup()

message(Sys.time(), " | CellRanger estimated cells summary:")
print(mapping_table %>% select(Sample_BR, Sample_ID, estimated_cells, dbr_use))

doublet_results_list <- list()

# -------------------------------------------------------------------------
# Loop
# -------------------------------------------------------------------------
bp_param <- BiocParallel::MulticoreParam(workers = 16)

message(Sys.time(), " | Starting doublet detection for ", nrow(mapping_table), " samples...")

for (i in 1:nrow(mapping_table)) {
  s_br <- mapping_table$Sample_BR[i]
  s_id <- mapping_table$Sample_ID[i]
  l_nm <- mapping_table$layer_name[i]
  dbr_use <- mapping_table$dbr_use[i]
  
  message(Sys.time(), " | Sample ", i, "/", nrow(mapping_table),
          ": ", s_id, " (", s_br, ") | dbr=", round(dbr_use, 4),
          " | est_cells=", mapping_table$estimated_cells[i])
  
  tryCatch({
    # 1) Subsetting
    cells_to_keep <- colnames(dataObject)[dataObject$Sample_BR == s_br]
    
    # sanity: ensure the layer actually contains those cells
    mat_cols <- colnames(GetAssayData(dataObject, assay = "RNA", layer = l_nm))
    if (!all(cells_to_keep %in% mat_cols)) {
      stop("Layer ", l_nm, " does not contain all cells for Sample_BR=", s_br,
           " (mapping mismatch).")
    }
    
    raw_counts <- GetAssayData(dataObject, assay = "RNA", layer = l_nm)[, cells_to_keep, drop = FALSE]
    sub_meta   <- dataObject@meta.data[cells_to_keep, , drop = FALSE]
    
    # 2) Local Processing
    obj_sample <- CreateSeuratObject(counts = raw_counts, meta.data = sub_meta)
    obj_sample <- NormalizeData(obj_sample, verbose = FALSE)
    obj_sample <- FindVariableFeatures(obj_sample, nfeatures = 2000, verbose = FALSE)
    obj_sample <- ScaleData(obj_sample, verbose = FALSE)
    obj_sample <- RunPCA(obj_sample, npcs = 30, verbose = FALSE, approx = TRUE)
    obj_sample <- FindNeighbors(obj_sample, dims = 1:20, verbose = FALSE)
    obj_sample <- FindClusters(obj_sample, resolution = 0.3, verbose = FALSE)
    
    # group safety (avoid unique() returning multiple values)
    grp <- unique(na.omit(obj_sample$group))
    if (length(grp) != 1) stop("Sample ", s_id, " has multiple/zero groups: ", paste(grp, collapse = ","))
    grp <- grp[1]
    
    # 3) scDblFinder
    sce <- as.SingleCellExperiment(obj_sample, assay = "RNA")
    
    sce <- scDblFinder(
      sce,
      clusters = "seurat_clusters",   # <-- updated
      dbr = dbr_use,
      BPPARAM = bp_param,
      verbose = FALSE
    )
    
    # 4) Store result
    doublet_results_list[[s_br]] <- data.frame(
      barcode = colnames(obj_sample),
      scDblFinder_class = sce$scDblFinder.class,
      scDblFinder_score = sce$scDblFinder.score,
      Sample_ID = s_id,
      Sample_BR = s_br,
      group = grp,
      estimated_cells = mapping_table$estimated_cells[i],
      dbr_used = dbr_use,
      stringsAsFactors = FALSE
    )
    
    rm(obj_sample, sce, raw_counts, sub_meta); gc()
    
  }, error = function(e) {
    message("Error processing ", s_id, " (", s_br, "): ", e$message)
  })
}

# Combine
all_results <- do.call(rbind, doublet_results_list)

# -------------------------------------------------------------------------
# Post-combine checks you requested
# -------------------------------------------------------------------------
if (anyDuplicated(all_results$barcode)) stop("Duplicate barcodes in all_results.")
match_idx <- match(colnames(dataObject), all_results$barcode)
if (anyNA(match_idx)) stop("Some cells missing from all_results.")

# -------------------------------------------------------------------------
# Map back + save
# -------------------------------------------------------------------------
message(Sys.time(), " | Mapping results back to main object...")
dataObject$scDblFinder_class <- all_results$scDblFinder_class[match_idx]
dataObject$scDblFinder_score <- all_results$scDblFinder_score[match_idx]

saveRDS(dataObject, paste0("../rObjects/", projectID, "_with_doublet_scores_exprate.rds"), compress = FALSE)

dataObject.singlets <- subset(dataObject, subset = scDblFinder_class == "singlet")
saveRDS(dataObject.singlets, paste0("../rObjects/", projectID, "_singlets_scDblFinder_exprate.rds"), compress = FALSE)

## ----visualizations-------------------------------------------------------
## ----visualizations-------------------------------------------------------
# Make sure plyr isn't hijacking summarize()
if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE, character.only = TRUE)
}

# [Stacked Bar Plot]
combined_doublet_data <- all_results %>%
  dplyr::group_by(Sample_ID, scDblFinder_class) %>%
  dplyr::tally(name = "Freq") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    scDblFinder_class = factor(scDblFinder_class, levels = c("doublet", "singlet"))
  )

label_data <- all_results %>%
  dplyr::group_by(Sample_ID) %>%
  dplyr::summarise(
    Total_Freq = dplyr::n(),
    pct = 100 * mean(scDblFinder_class == "doublet"),
    .groups = "drop"
  )

stacked_bar <- ggplot2::ggplot(combined_doublet_data, ggplot2::aes(x = Sample_ID, y = Freq)) +
  ggplot2::geom_bar(stat = "identity", position = "stack", ggplot2::aes(fill = scDblFinder_class)) +
  ggplot2::scale_fill_manual(values = c("singlet" = "#66C2A5", "doublet" = "black")) +
  ggplot2::geom_text(
    data = label_data,
    ggplot2::aes(x = Sample_ID, y = Total_Freq, label = paste0(round(pct, 1), "%")),
    vjust = -0.5, size = 3
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(title = "Doublet Classification by Sample")

ggplot2::ggsave(
  filename = file.path(out_dir, paste0(projectID, "_stacked_barplot.pdf")),
  plot = stacked_bar, width = 15, height = 6, units = "in"
)

# [Pathology Boxplot]
stats_df <- all_results %>%
  dplyr::group_by(Sample_ID, group) %>%
  dplyr::summarise(
    percent_doublets = 100 * mean(scDblFinder_class == "doublet"),
    .groups = "drop"
  )

stats_df$group <- factor(stats_df$group, levels = c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS"))

boxplot_path <- ggplot2::ggplot(stats_df, ggplot2::aes(x = group, y = percent_doublets, fill = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  ggplot2::geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c(
    "CONTROL" = "#4682B4", "AD_AT" = "#B4464B",
    "LBD_S" = "gray", "LBD_AS" = "gray35", "LBD_ATS" = "grey65"
  )) +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none") +
  ggplot2::labs(title = "Doublet Rates across Pathology Groups", y = "Doublet Percentage (%)")

ggplot2::ggsave(
  filename = file.path(out_dir, paste0(projectID, "_Doublet_Rate_by_Pathology.pdf")),
  plot = boxplot_path, width = 8, height = 6, units = "in"
)

message("--- End of Script ---")
message(Sys.time(), " | SLURM Job Complete.")
