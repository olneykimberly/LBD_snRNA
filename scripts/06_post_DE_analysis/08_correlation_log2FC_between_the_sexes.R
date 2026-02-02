#!/usr/bin/env Rscript
# ==============================================================================
# log2FC correlation plots from DEG tables (per cell type)
#
# Outputs:
# (A) Between-sex (female vs male) disease-vs-control (4 panels)
#   - corr_log2FC_<cell_type>_female_vs_male_disease_vs_CONTROL.pdf
#
# (B) Within-sex correlation plots across disease comparisons (3-panel PDFs),
#     produced separately for female and male:
#   - corr_log2FC_<cell_type>_<sex>_LBDtypes_vs_CONTROL.pdf
#   - corr_log2FC_<cell_type>_<sex>_LBDtypes_vs_AD_AT.pdf
#   - corr_log2FC_<cell_type>_<sex>_AD_AT_vs_CONTROL_vs_LBDtypes_vs_CONTROL.pdf
#
# DEG files:
#   ../results/DEGs_RNA_pct0.25_by_sex/<cell_type>/DEG_<cell_type>_<A>_vs_<B>.tsv
#
# Title format (4 lines):
#   1) cell type
#   2) sex
#   3) disease-vs-disease line (NO sex info; NO underscores)
#   4) r and n
#
# Labeling:
#   - Concordant (same direction): BLACK labels
#   - Discordant (opposite direction): RED labels
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

# ---- paths ----
deg_base_dir <- "../results/DEGs_RNA_pct0.25_by_sex"
out_dir <- "../results/log2FC_correlations_by_sex"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- infer cell types from directory structure ----
cell_types <- list.dirs(deg_base_dir, full.names = FALSE, recursive = FALSE)
cell_types <- cell_types[nzchar(cell_types)]
if (length(cell_types) == 0) {
  stop("No cell type subdirectories found under: ", deg_base_dir)
}

# ---- diseases ----
diseases <- c("AD_AT", "LBD_S", "LBD_AS", "LBD_ATS")
sexes <- c("female", "male")

# ==============================================================================
# Helpers
# ==============================================================================
pretty_comp <- function(x) gsub("_", " ", x)

# remove sex tokens from comparison strings (for the disease-vs-disease title line)
strip_sex <- function(x) gsub("_(female|male)", "", x)

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

merge_fc <- function(d1, d2) {
  if (is.null(d1) || is.null(d2)) return(NULL)
  colnames(d1) <- c("gene", "x")
  colnames(d2) <- c("gene", "y")
  m <- merge(d1, d2, by = "gene", all = FALSE)
  if (nrow(m) == 0) return(NULL)
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

# ---- title format is 4 lines exactly ----
# 1) cell type
# 2) sex line
# 3) disease comparison line (NO sex; no underscores)
# 4) r and n line
make_corr_plot <- function(df,
                           cell_type,
                           sex_line,
                           comp_line,
                           x_comp_label,
                           y_comp_label,
                           n_discordant = 5,
                           n_concordant_each = 3) {
  
  lim <- max(abs(c(df$x, df$y)), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  lim <- lim * 1.05
  
  r <- suppressWarnings(cor(df$x, df$y, method = "pearson", use = "pairwise.complete.obs"))
  n <- nrow(df)
  
  title_line1 <- cell_type
  title_line2 <- sex_line
  title_line3 <- comp_line
  title_line4 <- paste0(
    "r = ", ifelse(is.finite(r), sprintf("%.3f", r), "NA"),
    "   n = ", n
  )
  
  disc_df <- pick_discordant(df, top_n = n_discordant)
  conc_df <- pick_concordant(df, top_each = n_concordant_each)
  if (nrow(disc_df) > 0 && nrow(conc_df) > 0) {
    conc_df <- conc_df[!conc_df$gene %in% disc_df$gene, , drop = FALSE]
  }
  
  do_repel <- requireNamespace("ggrepel", quietly = TRUE)
  
  p <- ggplot(df, aes(x = x, y = y)) +
    annotate("rect", xmin = 0, xmax = lim, ymin = 0, ymax = lim,
             fill = "lightpink", alpha = 0.35) +
    annotate("rect", xmin = -lim, xmax = 0, ymin = -lim, ymax = 0,
             fill = "lightblue", alpha = 0.35) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.6, color = "grey40") +
    geom_point(size = 1.2, color = "black") +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(
      title = paste(title_line1, title_line2, title_line3, title_line4, sep = "\n"),
      x = paste0("log2FC (", pretty_comp(x_comp_label), ")"),
      y = paste0("log2FC (", pretty_comp(y_comp_label), ")")
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "plain", lineheight = 1.05),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 8)
    )
  
  if (nrow(conc_df) > 0) {
    if (do_repel) {
      p <- p + ggrepel::geom_text_repel(
        data = conc_df, aes(label = gene),
        size = 8 / ggplot2::.pt,
        color = "black",
        min.segment.length = 0,
        box.padding = 0.2,
        point.padding = 0.15,
        max.overlaps = 50
      )
    } else {
      p <- p + geom_text(
        data = conc_df, aes(label = gene),
        size = 8 / ggplot2::.pt,
        color = "black",
        vjust = -0.6
      )
    }
  }
  
  if (nrow(disc_df) > 0) {
    if (do_repel) {
      p <- p + ggrepel::geom_text_repel(
        data = disc_df, aes(label = gene),
        size = 8 / ggplot2::.pt,
        color = "red3",
        min.segment.length = 0,
        box.padding = 0.2,
        point.padding = 0.15,
        max.overlaps = 50
      )
    } else {
      p <- p + geom_text(
        data = disc_df, aes(label = gene),
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

save_four_panel <- function(p1, p2, p3, p4, out_pdf, width = 18, height = 4.2) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    comb <- p1 + p2 + p3 + p4 + patchwork::plot_layout(nrow = 1)
    ggsave(out_pdf, comb, width = width, height = height)
    return(invisible(NULL))
  }
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gr <- gridExtra::arrangeGrob(p1, p2, p3, p4, nrow = 1)
    ggsave(out_pdf, gr, width = width, height = height)
    return(invisible(NULL))
  }
  pdf(out_pdf, width = width, height = height)
  print(p1); print(p2); print(p3); print(p4)
  dev.off()
}

# ==============================================================================
# (A) Between-sex: disease vs CONTROL (4 panels)
# ==============================================================================
for (cell_type in cell_types) {
  
  message("Cell type (between-sex): ", cell_type)
  
  plots <- vector("list", length(diseases))
  
  for (i in seq_along(diseases)) {
    dis <- diseases[[i]]
    
    A_f <- paste0(dis, "_female")
    B_f <- "CONTROL_female"
    A_m <- paste0(dis, "_male")
    B_m <- "CONTROL_male"
    
    d_f <- read_deg_fc(deg_file(cell_type, A_f, B_f))
    d_m <- read_deg_fc(deg_file(cell_type, A_m, B_m))
    m <- merge_fc(d_f, d_m)
    
    xlab <- paste0(dis, "_female_vs_CONTROL_female")
    ylab <- paste0(dis, "_male_vs_CONTROL_male")
    
    sex_line <- "Sex: female vs male"
    # in between-sex plots, this line can include (female x, male y) since it's NOT redundant with sex_line
    comp_line <- paste0(pretty_comp(dis), " vs CONTROL (female x, male y)")
    
    if (is.null(m)) {
      msg <- paste(
        cell_type,
        sex_line,
        comp_line,
        "Missing file(s) or no overlap",
        sep = "\n"
      )
      plots[[i]] <- ggplot() + theme_void() +
        ggtitle(msg) +
        theme(plot.title = element_text(size = 9, face = "plain", lineheight = 1.05))
    } else {
      plots[[i]] <- make_corr_plot(
        m,
        cell_type = cell_type,
        sex_line = sex_line,
        comp_line = comp_line,
        x_comp_label = xlab,
        y_comp_label = ylab,
        n_discordant = 5,
        n_concordant_each = 3
      )
    }
  }
  
  out_pdf <- file.path(out_dir, paste0("corr_log2FC_", cell_type, "_female_vs_male_disease_vs_CONTROL.pdf"))
  save_four_panel(plots[[1]], plots[[2]], plots[[3]], plots[[4]], out_pdf)
  message("Wrote: ", out_pdf)
}

# ==============================================================================
# (B) Within-sex: compare disease comparisons (3-panel PDFs) for each sex
# ==============================================================================
pair_idx <- list(c(1,2), c(1,3), c(2,3))

for (cell_type in cell_types) {
  
  for (sex in sexes) {
    
    message("Cell type (within-sex): ", cell_type, "  sex: ", sex)
    
    sex_line <- paste0("Sex: ", sex)
    
    # ----- define within-sex comparison sets -----
    set1 <- list(
      c(paste0("LBD_S_", sex),   paste0("CONTROL_", sex)),
      c(paste0("LBD_AS_", sex),  paste0("CONTROL_", sex)),
      c(paste0("LBD_ATS_", sex), paste0("CONTROL_", sex))
    )
    
    set2 <- list(
      c(paste0("LBD_S_", sex),   paste0("AD_AT_", sex)),
      c(paste0("LBD_AS_", sex),  paste0("AD_AT_", sex)),
      c(paste0("LBD_ATS_", sex), paste0("AD_AT_", sex))
    )
    
    comparisons_to_plot <- list(
      c(paste0("AD_AT_", sex),   paste0("CONTROL_", sex)),
      c(paste0("LBD_S_", sex),   paste0("CONTROL_", sex)),
      c(paste0("LBD_AS_", sex),  paste0("CONTROL_", sex)),
      c(paste0("LBD_ATS_", sex), paste0("CONTROL_", sex))
    )
    
    pair_set3 <- list(
      c(1,2), # AD_AT vs CTRL  vs  LBD_S vs CTRL
      c(1,3), # AD_AT vs CTRL  vs  LBD_AS vs CTRL
      c(1,4)  # AD_AT vs CTRL  vs  LBD_ATS vs CTRL
    )
    
    # ----- SET 1: LBD types vs CONTROL (pairwise among LBD types) -----
    d_set1 <- lapply(set1, function(comp) {
      read_deg_fc(deg_file(cell_type, comp[[1]], comp[[2]]))
    })
    
    plots1 <- vector("list", length(pair_idx))
    for (k in seq_along(pair_idx)) {
      i <- pair_idx[[k]][1]; j <- pair_idx[[k]][2]
      
      x_comp <- paste0(set1[[i]][1], "_vs_", set1[[i]][2])
      y_comp <- paste0(set1[[j]][1], "_vs_", set1[[j]][2])
      
      m <- merge_fc(d_set1[[i]], d_set1[[j]])
      
      # comp line must NOT include sex (already on line above)
      comp_line <- paste0(
        pretty_comp(strip_sex(x_comp)),
        " vs ",
        pretty_comp(strip_sex(y_comp))
      )
      
      if (is.null(m)) {
        msg <- paste(
          cell_type,
          sex_line,
          comp_line,
          "Missing file(s) or no overlap",
          sep = "\n"
        )
        plots1[[k]] <- ggplot() + theme_void() +
          ggtitle(msg) +
          theme(plot.title = element_text(size = 9, face = "plain", lineheight = 1.05))
      } else {
        plots1[[k]] <- make_corr_plot(
          m,
          cell_type = cell_type,
          sex_line = sex_line,
          comp_line = comp_line,
          x_comp_label = x_comp,
          y_comp_label = y_comp,
          n_discordant = 5,
          n_concordant_each = 3
        )
      }
    }
    
    out_pdf1 <- file.path(out_dir, paste0("corr_log2FC_", cell_type, "_", sex, "_LBDtypes_vs_CONTROL.pdf"))
    save_three_panel(plots1[[1]], plots1[[2]], plots1[[3]], out_pdf1)
    
    # ----- SET 2: LBD types vs AD_AT (pairwise among LBD types) -----
    d_set2 <- lapply(set2, function(comp) {
      read_deg_fc(deg_file(cell_type, comp[[1]], comp[[2]]))
    })
    
    plots2 <- vector("list", length(pair_idx))
    for (k in seq_along(pair_idx)) {
      i <- pair_idx[[k]][1]; j <- pair_idx[[k]][2]
      
      x_comp <- paste0(set2[[i]][1], "_vs_", set2[[i]][2])
      y_comp <- paste0(set2[[j]][1], "_vs_", set2[[j]][2])
      
      m <- merge_fc(d_set2[[i]], d_set2[[j]])
      
      # comp line must NOT include sex (already on line above)
      comp_line <- paste0(
        pretty_comp(strip_sex(x_comp)),
        " vs ",
        pretty_comp(strip_sex(y_comp))
      )
      
      if (is.null(m)) {
        msg <- paste(
          cell_type,
          sex_line,
          comp_line,
          "Missing file(s) or no overlap",
          sep = "\n"
        )
        plots2[[k]] <- ggplot() + theme_void() +
          ggtitle(msg) +
          theme(plot.title = element_text(size = 9, face = "plain", lineheight = 1.05))
      } else {
        plots2[[k]] <- make_corr_plot(
          m,
          cell_type = cell_type,
          sex_line = sex_line,
          comp_line = comp_line,
          x_comp_label = x_comp,
          y_comp_label = y_comp,
          n_discordant = 5,
          n_concordant_each = 3
        )
      }
    }
    
    out_pdf2 <- file.path(out_dir, paste0("corr_log2FC_", cell_type, "_", sex, "_LBDtypes_vs_AD_AT.pdf"))
    save_three_panel(plots2[[1]], plots2[[2]], plots2[[3]], out_pdf2)
    
    # ----- SET 3: AD_AT vs CONTROL compared to each LBD type vs CONTROL -----
    d_set3 <- lapply(comparisons_to_plot, function(comp) {
      read_deg_fc(deg_file(cell_type, comp[[1]], comp[[2]]))
    })
    
    plots3 <- vector("list", length(pair_set3))
    for (k in seq_along(pair_set3)) {
      i <- pair_set3[[k]][1]; j <- pair_set3[[k]][2]
      
      x_comp <- paste0(comparisons_to_plot[[i]][1], "_vs_", comparisons_to_plot[[i]][2]) # AD_AT vs CTRL
      y_comp <- paste0(comparisons_to_plot[[j]][1], "_vs_", comparisons_to_plot[[j]][2]) # LBD_* vs CTRL
      
      m <- merge_fc(d_set3[[i]], d_set3[[j]])
      
      # comp line must NOT include sex (already on line above)
      comp_line <- paste0(
        pretty_comp(strip_sex(x_comp)),
        " vs ",
        pretty_comp(strip_sex(y_comp))
      )
      
      if (is.null(m)) {
        msg <- paste(
          cell_type,
          sex_line,
          comp_line,
          "Missing file(s) or no overlap",
          sep = "\n"
        )
        plots3[[k]] <- ggplot() + theme_void() +
          ggtitle(msg) +
          theme(plot.title = element_text(size = 9, face = "plain", lineheight = 1.05))
      } else {
        plots3[[k]] <- make_corr_plot(
          m,
          cell_type = cell_type,
          sex_line = sex_line,
          comp_line = comp_line,
          x_comp_label = x_comp,
          y_comp_label = y_comp,
          n_discordant = 5,
          n_concordant_each = 3
        )
      }
    }
    
    out_pdf3 <- file.path(out_dir, paste0("corr_log2FC_", cell_type, "_", sex, "_AD_AT_vs_CONTROL_vs_LBDtypes_vs_CONTROL.pdf"))
    save_three_panel(plots3[[1]], plots3[[2]], plots3[[3]], out_pdf3)
    
    message("Wrote:\n  ", out_pdf1, "\n  ", out_pdf2, "\n  ", out_pdf3)
  }
}
