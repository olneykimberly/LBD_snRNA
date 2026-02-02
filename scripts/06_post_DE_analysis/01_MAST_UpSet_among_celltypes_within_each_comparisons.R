#!/usr/bin/env Rscript
# make_upset_plots_complexupset_final.R
# ComplexUpset UpSet plotting with set sizes shown, neuron on top, sets descending by size, theme text size 8

library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(here)

# ---- working dir / project setup ----
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
project_ID <- "CWOW_cellbender"

# UpSet and table of shared DEGs among cell types within each comparison
main_output_dir <- "../results/DEGs_RNA_pct0.25"
out_plots_dir <- "../results/DEG_MAST_upset_plots"
out_shared_dir <- "../results/DEG_MAST_shared_with_stats"
dir.create(out_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_shared_dir, recursive = TRUE, showWarnings = FALSE)

# Fixed cell type order
fixed_cell_order <- c("neuron", "interneuron", "oligodendrocyte", "opc", "astrocyte", "microglia", "mural", "endothelial", "fibroblast")
cell_types <- fixed_cell_order

comparisons <- list(
  c("AD_AT", "CONTROL"),
  c("LBD_S", "CONTROL"),
  c("LBD_AS", "CONTROL"),
  c("LBD_ATS", "CONTROL"),
  c("LBD_S", "AD_AT"),
  c("LBD_AS", "AD_AT"),
  c("LBD_ATS", "AD_AT"),
  c("LBD_AS", "LBD_S"),
  c("LBD_ATS", "LBD_S"),
  c("LBD_ATS", "LBD_AS")
)

# thresholds
up_logfc_cutoff <- 0.25
down_logfc_cutoff <- -0.25
padj_cutoff   <- 0.05
min_shared    <- 2
min_intersection_size <- 10

# base theme size 8
base_theme_8 <- theme_grey(base_size = 8)

# debug helper
dump_mem_filtered_debug <- function(mem_filtered, comp_name, dir_label) {
  message("---- DEBUG: structure of mem_filtered for ", comp_name, " / ", dir_label, " ----")
  try(str(mem_filtered, max.level = 2), silent = TRUE)
  message("---- END DEBUG ----")
}

# ---- main loop ----
for (comp in comparisons) {
  pathA <- comp[1]; pathB <- comp[2]
  comp_name <- paste0(pathA, "_vs_", pathB)
  message("== ", comp_name, " ==")
  comp_plots_dir <- file.path(out_plots_dir, comp_name); dir.create(comp_plots_dir, recursive = TRUE, showWarnings = FALSE)
  comp_shared_dir <- file.path(out_shared_dir, comp_name); dir.create(comp_shared_dir, recursive = TRUE, showWarnings = FALSE)
  
  # read deg tables for cell types
  deg_tables <- setNames(vector("list", length(cell_types)), cell_types)
  up_sets <- setNames(vector("list", length(cell_types)), cell_types)
  down_sets <- setNames(vector("list", length(cell_types)), cell_types)
  
  for (ct in cell_types) {
    df <- read_deg(ct, pathA, pathB)
    if (is.null(df)) next
    deg_tables[[ct]] <- df
    up_sets[[ct]] <- unique(df$gene[!is.na(df$avg_log2FC) & !is.na(df$p_val_adj) & df$avg_log2FC > up_logfc_cutoff & df$p_val_adj < padj_cutoff])
    down_sets[[ct]] <- unique(df$gene[!is.na(df$avg_log2FC) & !is.na(df$p_val_adj) & df$avg_log2FC < down_logfc_cutoff & df$p_val_adj < padj_cutoff])
    message(" ", ct, ": up=", length(up_sets[[ct]]), " down=", length(down_sets[[ct]]))
  }
  
  # ---- UpSet plots (up / down) ----
  for (dir_label in c("up", "down")) {
    sets_list <- if (dir_label == "up") up_sets else down_sets
    mem <- membership_df(sets_list)
    mem_filtered <- filter_members_by_combo_size(mem, min_intersection_size)
    if (is.null(mem_filtered) || nrow(mem_filtered) == 0) {
      message(" No ", dir_label, " intersections >= ", min_intersection_size)
      next
    }
    
    # defensive types
    mem_filtered$gene <- as.character(mem_filtered$gene)
    set_cols <- setdiff(names(mem_filtered), "gene")
    for (sc in set_cols) {
      mem_filtered[[sc]] <- as.logical(mem_filtered[[sc]])
      mem_filtered[[sc]][is.na(mem_filtered[[sc]])] <- FALSE
    }
    
    # total per set, sort descending
    set_totals <- sapply(set_cols, function(x) sum(as.integer(mem_filtered[[x]])))
    ordered_sets <- names(sort(set_totals, decreasing = TRUE))
    intersect_factor_order <- rev(ordered_sets) # largest on top
    
    # title
    plot_title <- sprintf("%s — %s (avg_log2FC %s, p_val_adj < %s)",
                          comp_name,
                          ifelse(dir_label=="up","Upregulated","Downregulated"),
                          ifelse(dir_label=="up",paste0(">",up_logfc_cutoff),paste0("<",down_logfc_cutoff)),
                          padj_cutoff)
    
    theme_set(base_theme_8)
    
    my_themes <- upset_modify_themes(
      list(
        'intersections_matrix' = theme(axis.text=element_text(size=8), axis.title=element_text(size=8), plot.margin=margin(0,0,0,0,"cm"), text=element_text(size=8)),
        'overall_sizes' = theme(axis.text.x=element_text(size=8), axis.title=element_text(size=8), plot.margin=margin(0,0,0,0,"cm"), text=element_text(size=8)),
        'set_sizes' = theme(axis.text=element_text(size=8), axis.title=element_text(size=8), plot.margin=margin(0,0,0,0,"cm"), text=element_text(size=8))
      )
    )
    
    # intersection size annotation
    intersection_annotation <- intersection_size(
      mapping = aes(), 
      width = 0.7,
      text = list(size = 2)
    ) +
      scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
      theme(axis.text.y=element_text(size=8),
            axis.title.y=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.margin=margin(0,0,0,0.1,"cm"))
    
    # matrix
    matrix_geom <- intersection_matrix(
      geom = geom_point(shape=16, size=1.5, stroke=0.4)
    ) + theme(axis.text.y=element_text(size=8), axis.title.y=element_text(size=8), plot.margin=margin(0,0,0,0,"cm"))
    
    # set sizes
    set_size_annotation <- upset_set_size()
    
    # plot
    p <- ComplexUpset::upset(
      data = mem_filtered,
      intersect = intersect_factor_order,
      base_annotations = list('Intersection size' = intersection_annotation),
      matrix = matrix_geom,
      set_sizes = set_size_annotation,
      themes = my_themes,
      name = toupper(dir_label),
      width_ratio = 0.18,
      sort_sets = FALSE,
      sort_intersections = "descending"
    ) +
      ggtitle(plot_title) +
      theme(plot.title=element_text(size=8,hjust=0.5), strip.text=element_text(size=8), legend.text=element_text(size=8))
    
    # save PDF
    out_pdf <- file.path(comp_plots_dir, sprintf("UpSet_%s_%s.pdf", dir_label, comp_name))
    ggsave(out_pdf, plot=p, width=7, height=4.2, units="in", device=cairo_pdf)
    message(" Saved UpSet plot: ", out_pdf)
  }
  
  # ---- shared tables ----
  up_shared <- build_shared_table(up_sets, deg_tables, min_shared, "up")
  down_shared <- build_shared_table(down_sets, deg_tables, min_shared, "down")
  
  if (!is.null(up_shared)) {
    up_fn <- file.path(comp_shared_dir, paste0("Shared_DEGs_up_2plus_with_stats_", comp_name, ".tsv"))
    write.table(up_shared, up_fn, sep="\t", row.names=FALSE, quote=FALSE, na="")
    message(" Wrote UP shared table: ", up_fn, " (n=", nrow(up_shared), ")")
  } else message(" No UP shared genes >= ", min_shared)
  
  if (!is.null(down_shared)) {
    down_fn <- file.path(comp_shared_dir, paste0("Shared_DEGs_down_2plus_with_stats_", comp_name, ".tsv"))
    write.table(down_shared, down_fn, sep="\t", row.names=FALSE, quote=FALSE, na="")
    message(" Wrote DOWN shared table: ", down_fn, " (n=", nrow(down_shared), ")")
  } else message(" No DOWN shared genes >= ", min_shared)
  
  combined <- rbind(
    if (!is.null(up_shared)) up_shared else NULL,
    if (!is.null(down_shared)) down_shared else NULL
  )
  if (!is.null(combined) && nrow(combined) > 0) {
    combined_fn <- file.path(comp_shared_dir, paste0("Shared_DEGs_updown_2plus_with_stats_", comp_name, ".tsv"))
    combined <- combined[order(combined$direction, -combined$count, combined$gene), ]
    write.table(combined, combined_fn, sep="\t", row.names=FALSE, quote=FALSE, na="")
    message(" Wrote combined table: ", combined_fn, " (n=", nrow(combined), ")")
  }
  
} # comparisons loop

message("Done. Upset plots -> ", normalizePath(out_plots_dir), " ; shared tables -> ", normalizePath(out_shared_dir))
