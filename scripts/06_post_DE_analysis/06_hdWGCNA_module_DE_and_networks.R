#!/usr/bin/env Rscript
#================================================================================
# hdWGCNA per cell type (metacells built per cell type)
#
# For each cell type:
#   A) DME (Module Eigengenes): FindDMEs()  + lollipop per comparison (side-by-side PDF)
#   B) DMS (Module Scores):     ModuleExprScore() + FindDMEs(features="ModuleScores") + lollipop (side-by-side PDF)
#
# PLUS: Subcluster looping (cell_type_subcluster) to reproduce tutorial-style heatmaps
#   - For each comparison, loop across subclusters and run FindDMEs in each subcluster
#   - Make a heatmap (cluster x module) of avg_log2FC with significance stars
#
# Outputs (per cell type):
#   ../results/hdWGCNA/
#     lollipop_ME_<ct>_comparisons.pdf
#     lollipop_MS_<ct>_comparisons.pdf
#     heatmap_ME_<ct>_<comparison>.pdf
#     heatmap_MS_<ct>_<comparison>.pdf
#     tables_ME_<ct>_subclusters.tsv
#     tables_MS_<ct>_subclusters.tsv
#================================================================================

knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source("file_paths_and_colours.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(qs)
  library(cowplot)
  library(patchwork)
  library(dplyr)
  library(ggplot2)
  library(gtools)   # stars.pval
})

theme_set(theme_cowplot())
set.seed(12345)

#--------------------------------------------------------------------------------
# Cell types to run
#--------------------------------------------------------------------------------
cell_types <- c(
  "microglia",
  "astrocyte"
)

#--------------------------------------------------------------------------------
# Output
#--------------------------------------------------------------------------------
out_dir <- "../results/hdWGCNA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#--------------------------------------------------------------------------------
# Comparisons
#--------------------------------------------------------------------------------
comparisons_to_read <- list(
  c("AD_AT",   "CONTROL"),
  c("LBD_S",   "CONTROL"),
  c("LBD_AS",  "CONTROL"),
  c("LBD_ATS", "CONTROL")
)

#--------------------------------------------------------------------------------
# Module score parameters (DMS)
#--------------------------------------------------------------------------------
ms_n_genes <- 25
ms_method  <- "UCell"

#--------------------------------------------------------------------------------
# Helpers
#--------------------------------------------------------------------------------
get_barcodes <- function(obj, ct, grp, subcluster = NULL) {
  md <- obj@meta.data
  if (!is.null(subcluster)) {
    md |>
      dplyr::filter(.data$cell_type == ct,
                    .data$group == grp,
                    .data$cell_type_subcluster == subcluster) |>
      rownames()
  } else {
    md |>
      dplyr::filter(.data$cell_type == ct,
                    .data$group == grp) |>
      rownames()
  }
}

# Extract a ggplot from PlotDMEsLollipop output (handles list-returning versions)
as_ggplot_from_lollipop <- function(x) {
  if (inherits(x, "ggplot")) return(x)
  
  if (is.list(x)) {
    for (i in seq_along(x)) {
      if (inherits(x[[i]], "ggplot")) return(x[[i]])
    }
    for (i in seq_along(x)) {
      if (is.list(x[[i]])) {
        for (j in seq_along(x[[i]])) {
          if (inherits(x[[i]][[j]], "ggplot")) return(x[[i]][[j]])
        }
      }
    }
  }
  stop("Could not extract a ggplot from PlotDMEsLollipop() return object.")
}

infer_n_modules <- function(df) {
  if (!is.null(df) && nrow(df) > 0) {
    if ("module" %in% colnames(df)) return(length(unique(df$module)))
    if ("module_name" %in% colnames(df)) return(length(unique(df$module_name)))
    cn <- colnames(df)
    mod_col <- cn[grepl("module", cn, ignore.case = TRUE)][1]
    if (!is.na(mod_col) && nzchar(mod_col)) return(length(unique(df[[mod_col]])))
  }
  30L
}

# get module order from hdWGCNA modules table
get_module_levels <- function(seurat_obj, wgcna_name = "CWOW") {
  modules <- GetModules(seurat_obj, wgcna_name = wgcna_name)
  mods <- levels(modules$module)
  mods[mods != "grey"]
}

# build a tutorial-style heatmap for ONE comparison from a DME/DMS subcluster table
make_subcluster_heatmap <- function(df, mods, title, maxval = 0.5, minval = -0.5) {
  plot_df <- df
  
  # enforce module order
  plot_df$module <- factor(as.character(plot_df$module), levels = mods)
  
  # clip effect sizes for nicer color scaling
  plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
  plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)
  
  # add significance stars
  plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)
  
  # text color for readability
  plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, "black", "white")
  
  ggplot(plot_df, aes(y = cluster, x = module, fill = avg_log2FC)) +
    geom_tile() +
    geom_text(aes(label = Significance, color = textcolor), size = 2, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
    Seurat::RotatedAxis() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 0.6),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = margin(2, 2, 2, 2),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7)
    ) +
    labs(title = title, x = NULL, y = NULL) +
    coord_fixed()
}

#================================================================================
# LOOP OVER CELL TYPES
#================================================================================
for (ct in cell_types) {
  
  message("======================================")
  message("Running DME + DMS (and subcluster heatmaps) for cell type: ", ct)
  message("======================================")
  
  rds_path <- paste0("../rObjects/hdWGCNA_", ct, "_final.rds")
  if (!file.exists(rds_path)) {
    warning("Missing file: ", rds_path, " -- skipping ", ct)
    next
  }
  dataObject.ct <- readRDS(rds_path)
  
  # sanity: require subcluster column
  if (!"cell_type_subcluster" %in% colnames(dataObject.ct@meta.data)) {
    warning("[", ct, "] meta.data lacks 'cell_type_subcluster' -- skipping subcluster heatmaps, will still do lollipops.")
  }
  
  # module order for heatmaps
  mods <- get_module_levels(dataObject.ct, wgcna_name = "CWOW")
  
  #--------------------------------------------------------------------------------
  # A) DME (Module Eigengenes): per comparison lollipops (all cells in ct)
  #--------------------------------------------------------------------------------
  message("[", ct, "] DME: all-cells lollipops...")
  
  DME_list <- list()
  for (cmp in comparisons_to_read) {
    g1 <- cmp[[1]]; g2 <- cmp[[2]]
    cmp_name <- paste0(g1, "_vs_", g2)
    
    group1 <- get_barcodes(dataObject.ct, ct = ct, grp = g1)
    group2 <- get_barcodes(dataObject.ct, ct = ct, grp = g2)
    
    message("[", ct, "][DME] ", cmp_name, " : n1=", length(group1), " n2=", length(group2))
    if (length(group1) == 0L || length(group2) == 0L) next
    
    dme <- FindDMEs(
      dataObject.ct,
      barcodes1 = group1,
      barcodes2 = group2,
      test.use  = "wilcox"
    )
    dme$comparison <- cmp_name
    dme$group1 <- g1
    dme$group2 <- g2
    DME_list[[cmp_name]] <- dme
  }
  
  if (length(DME_list) > 0L) {
    DME_all <- dplyr::bind_rows(DME_list, .id = "comparison_id")
    n_modules_DME <- infer_n_modules(DME_all)
    
    DME_plots <- list()
    for (cmp_name in names(DME_list)) {
      p_raw <- PlotDMEsLollipop(
        dataObject.ct,
        DME_list[[cmp_name]],
        wgcna_name = "CWOW",
        pvalue     = "p_val_adj"
      )
      p <- as_ggplot_from_lollipop(p_raw) +
        ggplot2::ggtitle(paste0(ct, " (ME): ", cmp_name))
      DME_plots[[cmp_name]] <- p
    }
    
    p_combined_DME <- patchwork::wrap_plots(DME_plots, nrow = 1)
    
    out_pdf_DME <- file.path(out_dir, paste0("lollipop_ME_", ct, "_comparisons.pdf"))
    pdf(out_pdf_DME,
        width  = 4.5 * length(DME_plots),
        height = 0.35 * n_modules_DME + 2.5)
    print(p_combined_DME)
    dev.off()
    message("[", ct, "][DME] Saved: ", out_pdf_DME)
  } else {
    warning("[", ct, "][DME] No all-cells DME results computed.")
  }
  
  #--------------------------------------------------------------------------------
  # A2) DME subcluster loop + heatmaps
  #--------------------------------------------------------------------------------
  if ("cell_type_subcluster" %in% colnames(dataObject.ct@meta.data)) {
    
    clusters <- sort(unique(dataObject.ct@meta.data$cell_type_subcluster))
    # keep only this cell type's subclusters (some objects might still contain others)
    clusters <- clusters[grepl(paste0("^", ct), clusters) | TRUE]  # safe: keep all if already ct-only
    
    message("[", ct, "] DME: subcluster loop across ", length(clusters), " clusters...")
    
    DME_sc <- data.frame()
    
    for (cmp in comparisons_to_read) {
      g1 <- cmp[[1]]; g2 <- cmp[[2]]
      cmp_name <- paste0(g1, "_vs_", g2)
      
      for (cur_cluster in clusters) {
        group1 <- get_barcodes(dataObject.ct, ct = ct, grp = g1, subcluster = cur_cluster)
        group2 <- get_barcodes(dataObject.ct, ct = ct, grp = g2, subcluster = cur_cluster)
        
        if (length(group1) == 0L || length(group2) == 0L) next
        
        cur_DMEs <- FindDMEs(
          dataObject.ct,
          barcodes1 = group1,
          barcodes2 = group2,
          test.use  = "wilcox"
          # pseudocount.use can be added here if desired
        )
        
        cur_DMEs$cluster <- cur_cluster
        cur_DMEs$comparison <- cmp_name
        cur_DMEs$group1 <- g1
        cur_DMEs$group2 <- g2
        
        DME_sc <- rbind(DME_sc, cur_DMEs)
      }
    }
    
    # save table
    if (nrow(DME_sc) > 0) {
      out_tsv <- file.path(out_dir, paste0("tables_ME_", ct, "_subclusters.tsv"))
      write.table(DME_sc, out_tsv, sep = "\t", quote = FALSE, row.names = TRUE)
      message("[", ct, "][DME] Saved subcluster table: ", out_tsv)
      
      # heatmap per comparison
      for (cmp in comparisons_to_read) {
        g1 <- cmp[[1]]; g2 <- cmp[[2]]
        cmp_name <- paste0(g1, "_vs_", g2)
        
        plot_df <- DME_sc |> dplyr::filter(.data$comparison == cmp_name)
        if (nrow(plot_df) == 0) next
        
        p_hm <- make_subcluster_heatmap(
          df    = plot_df,
          mods  = mods,
          title = paste0(ct, " (ME) subclusters: ", cmp_name),
          maxval = 0.5,
          minval = -0.5
        )
        
        out_pdf <- file.path(out_dir, paste0("heatmap_subcluster_ME_", ct, "_", cmp_name, ".pdf"))
        pdf(out_pdf, width = 0.28 * length(mods) + 4, height = 0.28 * length(unique(plot_df$cluster)) + 3)
        print(p_hm)
        dev.off()
        message("[", ct, "][DME] Saved heatmap: ", out_pdf)
      }
    } else {
      warning("[", ct, "][DME] No subcluster DME results computed (likely missing cells for comparisons).")
    }
  }
  
  #--------------------------------------------------------------------------------
  # B) DMS (Module Scores): compute scores once, then lollipops + subcluster heatmaps
  #--------------------------------------------------------------------------------
  message("[", ct, "] DMS: computing ModuleScores (", ms_method, ", n_genes=", ms_n_genes, ")...")
  
  dataObject.ct <- ModuleExprScore(
    dataObject.ct,
    n_genes = ms_n_genes,
    method  = ms_method
  )
  
  #--------------------------------------------------------------------------------
  # B1) DMS all-cells lollipops
  #--------------------------------------------------------------------------------
  message("[", ct, "] DMS: all-cells lollipops...")
  
  DMS_list <- list()
  for (cmp in comparisons_to_read) {
    g1 <- cmp[[1]]; g2 <- cmp[[2]]
    cmp_name <- paste0(g1, "_vs_", g2)
    
    group1 <- get_barcodes(dataObject.ct, ct = ct, grp = g1)
    group2 <- get_barcodes(dataObject.ct, ct = ct, grp = g2)
    
    message("[", ct, "][DMS] ", cmp_name, " : n1=", length(group1), " n2=", length(group2))
    if (length(group1) == 0L || length(group2) == 0L) next
    
    dms <- FindDMEs(
      dataObject.ct,
      features   = "ModuleScores",
      barcodes1  = group1,
      barcodes2  = group2,
      test.use   = "wilcox",
      wgcna_name = "CWOW"
    )
    dms$comparison <- cmp_name
    dms$group1 <- g1
    dms$group2 <- g2
    DMS_list[[cmp_name]] <- dms
  }
  
  if (length(DMS_list) > 0L) {
    DMS_all <- dplyr::bind_rows(DMS_list, .id = "comparison_id")
    n_modules_DMS <- infer_n_modules(DMS_all)
    
    DMS_plots <- list()
    for (cmp_name in names(DMS_list)) {
      p_raw <- PlotDMEsLollipop(
        dataObject.ct,
        DMS_list[[cmp_name]],
        wgcna_name = "CWOW",
        pvalue     = "p_val_adj"
      )
      p <- as_ggplot_from_lollipop(p_raw) +
        ggplot2::ggtitle(paste0(ct, " (MS): ", cmp_name))
      DMS_plots[[cmp_name]] <- p
    }
    
    p_combined_DMS <- patchwork::wrap_plots(DMS_plots, nrow = 1)
    
    out_pdf_DMS <- file.path(out_dir, paste0("lollipop_MS_", ct, "_comparisons.pdf"))
    pdf(out_pdf_DMS,
        width  = 4.5 * length(DMS_plots),
        height = 0.35 * n_modules_DMS + 2.5)
    print(p_combined_DMS)
    dev.off()
    message("[", ct, "][DMS] Saved: ", out_pdf_DMS)
  } else {
    warning("[", ct, "][DMS] No all-cells DMS results computed.")
  }
  
  #--------------------------------------------------------------------------------
  # B2) DMS subcluster loop + heatmaps
  #--------------------------------------------------------------------------------
  if ("cell_type_subcluster" %in% colnames(dataObject.ct@meta.data)) {
    
    clusters <- sort(unique(dataObject.ct@meta.data$cell_type_subcluster))
    clusters <- clusters[grepl(paste0("^", ct), clusters) | TRUE]
    
    message("[", ct, "] DMS: subcluster loop across ", length(clusters), " clusters...")
    
    DMS_sc <- data.frame()
    
    for (cmp in comparisons_to_read) {
      g1 <- cmp[[1]]; g2 <- cmp[[2]]
      cmp_name <- paste0(g1, "_vs_", g2)
      
      for (cur_cluster in clusters) {
        group1 <- get_barcodes(dataObject.ct, ct = ct, grp = g1, subcluster = cur_cluster)
        group2 <- get_barcodes(dataObject.ct, ct = ct, grp = g2, subcluster = cur_cluster)
        
        if (length(group1) == 0L || length(group2) == 0L) next
        
        cur_DMS <- FindDMEs(
          dataObject.ct,
          features   = "ModuleScores",
          barcodes1  = group1,
          barcodes2  = group2,
          test.use   = "wilcox",
          wgcna_name = "CWOW"
        )
        
        cur_DMS$cluster <- cur_cluster
        cur_DMS$comparison <- cmp_name
        cur_DMS$group1 <- g1
        cur_DMS$group2 <- g2
        
        DMS_sc <- rbind(DMS_sc, cur_DMS)
      }
    }
    
    if (nrow(DMS_sc) > 0) {
      out_tsv <- file.path(out_dir, paste0("tables_MS_", ct, "_subclusters.tsv"))
      write.table(DMS_sc, out_tsv, sep = "\t", quote = FALSE, row.names = TRUE)
      message("[", ct, "][DMS] Saved subcluster table: ", out_tsv)
      
      for (cmp in comparisons_to_read) {
        g1 <- cmp[[1]]; g2 <- cmp[[2]]
        cmp_name <- paste0(g1, "_vs_", g2)
        
        plot_df <- DMS_sc |> dplyr::filter(.data$comparison == cmp_name)
        if (nrow(plot_df) == 0) next
        
        p_hm <- make_subcluster_heatmap(
          df    = plot_df,
          mods  = mods,
          title = paste0(ct, " (MS) subclusters: ", cmp_name),
          maxval = 0.5,
          minval = -0.5
        )
        
        out_pdf <- file.path(out_dir, paste0("heatmap_subcluster_MS_", ct, "_", cmp_name, ".pdf"))
        pdf(out_pdf, width = 0.28 * length(mods) + 4, height = 0.28 * length(unique(plot_df$cluster)) + 3)
        print(p_hm)
        dev.off()
        message("[", ct, "][DMS] Saved heatmap: ", out_pdf)
      }
    } else {
      warning("[", ct, "][DMS] No subcluster DMS results computed (likely missing cells for comparisons).")
    }
  }
  
}

message("ALL CELL TYPES COMPLETED SUCCESSFULLY")
