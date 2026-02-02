#!/usr/bin/env Rscript
# ==============================================================================
# Correlation plots of log2FC between comparisons (per cell type)
#
# Updates requested:
#  1) Keep labels for concordant (shared-direction) genes in BLACK
#     and discordant (opposite-direction) genes in RED.
#     - Concordant = x,y both >0 (up/up) OR both <0 (down/down)
#     - Discordant = sign differs between x and y
#     - We label top genes by:
#         * Discordant: largest |x - y| (direction flip + big delta)
#         * Concordant: largest |x|+|y| within up/up and down/down (top 3 each)
#  2) Do NOT bold the title
#  3) Add third comparison set:
#       LBD types vs CONTROL compared to AD_AT vs CONTROL
#     with 3 plots side-by-side:
#       AD_AT vs CONTROL  vs  LBD_S vs CONTROL
#       AD_AT vs CONTROL  vs  LBD_AS vs CONTROL
#       AD_AT vs CONTROL  vs  LBD_ATS vs CONTROL
#
# Output (per cell type):
#   - corr_log2FC_<cell_type>_LBDtypes_vs_CONTROL.pdf
#   - corr_log2FC_<cell_type>_LBDtypes_vs_AD_AT.pdf
#   - corr_log2FC_<cell_type>_AD_AT_vs_CONTROL_vs_LBDtypes_vs_CONTROL.pdf
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

deg_base_dir <- "../results/DEGs_MAST_RNA_pct0.25"
out_dir <- "../results/log2FC_correlations_LBDtypes"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- comparison sets ----
set1 <- list(
  c("LBD_S",  "CONTROL"),
  c("LBD_AS", "CONTROL"),
  c("LBD_ATS","CONTROL")
)

set2 <- list(
  c("LBD_S",  "AD_AT"),
  c("LBD_AS", "AD_AT"),
  c("LBD_ATS","AD_AT")
)

comparisons_to_plot <- list(
  c("AD_AT",  "CONTROL"),
  c("LBD_S",  "CONTROL"),
  c("LBD_AS", "CONTROL"),
  c("LBD_ATS","CONTROL")
)

# For set1/set2 (3 plots each)
pair_idx <- list(c(1,2), c(1,3), c(2,3))

# For set3: AD_AT vs CONTROL compared to each LBD type vs CONTROL (3 plots)
pair_set3 <- list(
  c(1,2), # AD_AT vs CTRL vs LBD_S vs CTRL
  c(1,3), # AD_AT vs CTRL vs LBD_AS vs CTRL
  c(1,4)  # AD_AT vs CTRL vs LBD_ATS vs CTRL
)

# ---- helpers ----
pretty_comp <- function(x) gsub("_", " ", x)

deg_file <- function(cell_type, A, B) {
  file.path(deg_base_dir, cell_type, paste0("DEG_", cell_type, "_", A, "_vs_", B, ".tsv"))
}

read_deg_fc <- function(path) {
  if (!file.exists(path)) return(NULL)
  d <- readr::read_tsv(path, show_col_types = FALSE)
  
  req <- c("gene", "avg_log2FC")
  if (!all(req %in% colnames(d))) return(NULL)
  
  d <- d[, req]
  d$gene <- as.character(d$gene)
  d$avg_log2FC <- suppressWarnings(as.numeric(d$avg_log2FC))
  d <- d[!is.na(d$gene) & d$gene != "" & !is.na(d$avg_log2FC), , drop = FALSE]
  
  # de-dup gene (keep largest abs FC)
  if (any(duplicated(d$gene))) {
    ord <- order(d$gene, -abs(d$avg_log2FC))
    d <- d[ord, , drop = FALSE]
    d <- d[!duplicated(d$gene), , drop = FALSE]
  }
  d
}

merge_fc <- function(d1, d2, xlab, ylab) {
  if (is.null(d1) || is.null(d2)) return(NULL)
  colnames(d1) <- c("gene", "x")
  colnames(d2) <- c("gene", "y")
  m <- merge(d1, d2, by = "gene", all = FALSE)
  if (nrow(m) == 0) return(NULL)
  attr(m, "xlab") <- xlab
  attr(m, "ylab") <- ylab
  m
}

# Discordant labels: opposite direction + largest |x - y|
pick_discordant <- function(df, top_n = 5) {
  disc <- df[(df$x > 0 & df$y < 0) | (df$x < 0 & df$y > 0), , drop = FALSE]
  if (nrow(disc) == 0) return(disc)
  disc$score <- abs(disc$x - disc$y)
  disc <- disc[order(-disc$score), , drop = FALSE]
  disc[seq_len(min(top_n, nrow(disc))), , drop = FALSE]
}

# Concordant labels: up/up and down/down, pick top by |x|+|y|
pick_concordant <- function(df, top_each = 3) {
  df$score <- abs(df$x) + abs(df$y)
  upup <- df[df$x > 0 & df$y > 0, , drop = FALSE]
  dndn <- df[df$x < 0 & df$y < 0, , drop = FALSE]
  
  take <- function(sub, n) {
    if (nrow(sub) == 0) return(sub)
    sub <- sub[order(-sub$score), , drop = FALSE]
    sub[seq_len(min(n, nrow(sub))), , drop = FALSE]
  }
  
  rbind(take(upup, top_each), take(dndn, top_each))
}

make_corr_plot <- function(df, cell_type, xlab, ylab,
                           n_discordant = 5, n_concordant_each = 3) {
  
  # symmetric limits for display
  lim <- max(abs(c(df$x, df$y)), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  lim <- lim * 1.05
  
  r <- suppressWarnings(cor(df$x, df$y, method = "pearson", use = "pairwise.complete.obs"))
  n <- nrow(df)
  
  xlab_pretty <- pretty_comp(xlab)
  ylab_pretty <- pretty_comp(ylab)
  
  # title layout: cell type, comparison, r/n
  title_line1 <- cell_type
  title_line2 <- paste0(xlab_pretty, "  vs  ", ylab_pretty)
  title_line3 <- paste0("r = ", ifelse(is.finite(r), sprintf("%.3f", r), "NA"),
                        "   n = ", n)
  
  # labels
  disc_df <- pick_discordant(df, top_n = n_discordant)
  conc_df <- pick_concordant(df, top_each = n_concordant_each)
  
  # remove any overlap (prioritize discordant as red)
  if (nrow(disc_df) > 0 && nrow(conc_df) > 0) {
    conc_df <- conc_df[!conc_df$gene %in% disc_df$gene, , drop = FALSE]
  }
  
  do_repel <- requireNamespace("ggrepel", quietly = TRUE)
  
  p <- ggplot(df, aes(x = x, y = y)) +
    # concordant quadrants shading
    annotate("rect", xmin = 0, xmax = lim, ymin = 0, ymax = lim,
             fill = "lightpink", alpha = 0.35) +
    annotate("rect", xmin = -lim, xmax = 0, ymin = -lim, ymax = 0,
             fill = "lightblue", alpha = 0.35) +
    # diagonal
    geom_abline(intercept = 0, slope = 1, linewidth = 0.6, color = "grey40") +
    # points
    geom_point(size = 1.2, color = "black") +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(
      title = paste(title_line1, title_line2, title_line3, sep = "\n"),
      x = paste0("log2FC (", xlab_pretty, ")"),
      y = paste0("log2FC (", ylab_pretty, ")")
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "plain", lineheight = 1.05),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 8)
    )
  
  # concordant labels (black)
  if (nrow(conc_df) > 0) {
    if (do_repel) {
      p <- p + ggrepel::geom_text_repel(
        data = conc_df,
        aes(label = gene),
        size = 8 / ggplot2::.pt,
        color = "black",
        min.segment.length = 0,
        box.padding = 0.2,
        point.padding = 0.15,
        max.overlaps = 50
      )
    } else {
      p <- p + geom_text(
        data = conc_df,
        aes(label = gene),
        size = 8 / ggplot2::.pt,
        color = "black",
        vjust = -0.6
      )
    }
  }
  
  # discordant labels (red)
  if (nrow(disc_df) > 0) {
    if (do_repel) {
      p <- p + ggrepel::geom_text_repel(
        data = disc_df,
        aes(label = gene),
        size = 8 / ggplot2::.pt,
        color = "red3",
        min.segment.length = 0,
        box.padding = 0.2,
        point.padding = 0.15,
        max.overlaps = 50
      )
    } else {
      p <- p + geom_text(
        data = disc_df,
        aes(label = gene),
        size = 8 / ggplot2::.pt,
        color = "red3",
        vjust = -0.6
      )
    }
  }
  
  p
}

save_three_panel <- function(p1, p2, p3, out_pdf, width = 13.5, height = 4.2) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    comb <- p1 + p2 + p3 + patchwork::plot_layout(nrow = 1)
    ggsave(out_pdf, comb, width = width, height = height)
    return(invisible(NULL))
  }
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gr <- gridExtra::arrangeGrob(p1, p2, p3, nrow = 1)
    ggsave(out_pdf, gr, width = width, height = height)
    return(invisible(NULL))
  }
  pdf(out_pdf, width = width, height = height)
  print(p1); print(p2); print(p3)
  dev.off()
}

# ==============================================================================
# Main loop per cell type
# ==============================================================================
for (cell_type in cell_types) {
  
  message("Cell type: ", cell_type)
  
  # -------- SET 1 (LBD types vs CONTROL; pairwise among LBD types) --------
  d_set1 <- lapply(set1, function(comp) {
    read_deg_fc(deg_file(cell_type, comp[[1]], comp[[2]]))
  })
  
  plots1 <- vector("list", length(pair_idx))
  for (k in seq_along(pair_idx)) {
    i <- pair_idx[[k]][1]; j <- pair_idx[[k]][2]
    xlab <- paste0(set1[[i]][1], "_vs_", set1[[i]][2])
    ylab <- paste0(set1[[j]][1], "_vs_", set1[[j]][2])
    
    m <- merge_fc(d_set1[[i]], d_set1[[j]], xlab, ylab)
    if (is.null(m)) {
      plots1[[k]] <- ggplot() + theme_void() +
        ggtitle(paste0(cell_type, "\n", pretty_comp(xlab), " vs ", pretty_comp(ylab), "\nMissing file(s) or no overlap")) +
        theme(plot.title = element_text(size = 9))
    } else {
      plots1[[k]] <- make_corr_plot(m, cell_type, xlab, ylab,
                                    n_discordant = 5, n_concordant_each = 3)
    }
  }
  out_pdf1 <- file.path(out_dir, paste0("corr_log2FC_", cell_type, "_LBDtypes_vs_CONTROL.pdf"))
  save_three_panel(plots1[[1]], plots1[[2]], plots1[[3]], out_pdf1)
  
  # -------- SET 2 (LBD types vs AD_AT; pairwise among LBD types) --------
  d_set2 <- lapply(set2, function(comp) {
    read_deg_fc(deg_file(cell_type, comp[[1]], comp[[2]]))
  })
  
  plots2 <- vector("list", length(pair_idx))
  for (k in seq_along(pair_idx)) {
    i <- pair_idx[[k]][1]; j <- pair_idx[[k]][2]
    xlab <- paste0(set2[[i]][1], "_vs_", set2[[i]][2])
    ylab <- paste0(set2[[j]][1], "_vs_", set2[[j]][2])
    
    m <- merge_fc(d_set2[[i]], d_set2[[j]], xlab, ylab)
    if (is.null(m)) {
      plots2[[k]] <- ggplot() + theme_void() +
        ggtitle(paste0(cell_type, "\n", pretty_comp(xlab), " vs ", pretty_comp(ylab), "\nMissing file(s) or no overlap")) +
        theme(plot.title = element_text(size = 9))
    } else {
      plots2[[k]] <- make_corr_plot(m, cell_type, xlab, ylab,
                                    n_discordant = 5, n_concordant_each = 3)
    }
  }
  out_pdf2 <- file.path(out_dir, paste0("corr_log2FC_", cell_type, "_LBDtypes_vs_AD_AT.pdf"))
  save_three_panel(plots2[[1]], plots2[[2]], plots2[[3]], out_pdf2)
  
  # -------- SET 3 (AD_AT vs CONTROL compared to each LBD type vs CONTROL) --------
  d_set3 <- lapply(comparisons_to_plot, function(comp) {
    read_deg_fc(deg_file(cell_type, comp[[1]], comp[[2]]))
  })
  
  plots3 <- vector("list", length(pair_set3))
  for (k in seq_along(pair_set3)) {
    i <- pair_set3[[k]][1]; j <- pair_set3[[k]][2]
    xlab <- paste0(comparisons_to_plot[[i]][1], "_vs_", comparisons_to_plot[[i]][2]) # AD_AT vs CONTROL
    ylab <- paste0(comparisons_to_plot[[j]][1], "_vs_", comparisons_to_plot[[j]][2]) # LBD_* vs CONTROL
    
    m <- merge_fc(d_set3[[i]], d_set3[[j]], xlab, ylab)
    if (is.null(m)) {
      plots3[[k]] <- ggplot() + theme_void() +
        ggtitle(paste0(cell_type, "\n", pretty_comp(xlab), " vs ", pretty_comp(ylab), "\nMissing file(s) or no overlap")) +
        theme(plot.title = element_text(size = 9))
    } else {
      plots3[[k]] <- make_corr_plot(m, cell_type, xlab, ylab,
                                    n_discordant = 5, n_concordant_each = 3)
    }
  }
  out_pdf3 <- file.path(out_dir, paste0("corr_log2FC_", cell_type, "_AD_AT_vs_CONTROL_vs_LBDtypes_vs_CONTROL.pdf"))
  save_three_panel(plots3[[1]], plots3[[2]], plots3[[3]], out_pdf3)
  
  message("Wrote:\n  ", out_pdf1, "\n  ", out_pdf2, "\n  ", out_pdf3)
}
