#!/usr/bin/env Rscript
# ==============================================================================
# hdWGCNA modules: fraction of DEGs per module (up vs down; diverging bars)
# - Grey module removed
# - Fixed x-scale across facets
# - Plot height adapts tightly to number of modules (more compact for small n_mod)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(tidyr)
})

hdwgcna_dir  <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/hdWGCNA"
deg_base_dir <- "../results/DEGs_MAST_RNA_pct0.25"
out_dir      <- "../results/hdWGCNA_DEG_fraction_by_module"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

qval <- 0.05
FC   <- 0.25

comparisons <- list(
  c("AD_AT",  "CONTROL"),
  c("LBD_S",  "CONTROL"),
  c("LBD_AS", "CONTROL"),
  c("LBD_ATS","CONTROL")
)

modules_file_for_celltype <- function(cell_type) {
  file.path(hdwgcna_dir, paste0("modules_", cell_type, ".tsv"))
}

deg_file_for_celltype_comp <- function(cell_type, A, B) {
  file.path(deg_base_dir, cell_type, paste0("DEG_", cell_type, "_", A, "_vs_", B, ".tsv"))
}

read_modules <- function(cell_type) {
  f <- modules_file_for_celltype(cell_type)
  if (!file.exists(f)) {
    message("Modules file missing for ", cell_type, ": ", f)
    return(NULL)
  }
  m <- readr::read_tsv(f, show_col_types = FALSE)
  
  req <- c("gene_name", "module")
  if (!all(req %in% colnames(m))) {
    stop("Modules file ", f, " missing required columns: ",
         paste(setdiff(req, colnames(m)), collapse = ", "))
  }
  
  m %>%
    transmute(
      gene   = as.character(gene_name),
      module = as.character(module)
    ) %>%
    filter(!is.na(gene), gene != "",
           !is.na(module), module != "") %>%
    filter(module != "grey")
}

# Return a DEG dataframe (not just genes) so we can split up/down
read_deg_df <- function(cell_type, A, B, qval = 0.05, FC = 0.25) {
  f <- deg_file_for_celltype_comp(cell_type, A, B)
  if (!file.exists(f)) {
    message("DEG file missing: ", f)
    return(NULL)
  }
  
  d <- readr::read_tsv(f, show_col_types = FALSE)
  
  req <- c("gene", "p_val_adj", "avg_log2FC")
  if (!all(req %in% colnames(d))) {
    message("Skipping DEG file (missing cols): ", f, " :: ",
            paste(setdiff(req, colnames(d)), collapse = ", "))
    return(NULL)
  }
  
  d %>%
    filter(!is.na(gene), gene != "",
           !is.na(p_val_adj), !is.na(avg_log2FC)) %>%
    filter(p_val_adj < qval, abs(avg_log2FC) > FC) %>%
    transmute(
      gene = as.character(gene),
      avg_log2FC = as.numeric(avg_log2FC)
    )
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
for (cell_type in cell_types) {
  
  modules <- read_modules(cell_type)
  if (is.null(modules) || nrow(modules) == 0L) next
  
  module_sizes <- modules %>%
    group_by(module) %>%
    summarise(module_size = length(gene), .groups = "drop")
  
  per_comp <- lapply(comparisons, function(comp) {
    A <- comp[[1]]
    B <- comp[[2]]
    comp_tag <- paste0(A, "_vs_", B)
    
    deg_df <- read_deg_df(cell_type, A, B, qval, FC)
    if (is.null(deg_df) || nrow(deg_df) == 0L) {
      return(
        module_sizes %>%
          mutate(
            n_up = 0L, n_down = 0L,
            frac_up = 0, frac_down = 0,
            frac_up_plot = 0, frac_down_plot = 0,
            comparison = comp_tag
          )
      )
    }
    
    deg_mod <- modules %>%
      inner_join(deg_df, by = "gene") %>%
      mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down"))
    
    dir_counts <- deg_mod %>%
      group_by(module, direction) %>%
      summarise(n = length(gene), .groups = "drop")
    
    up_counts <- dir_counts %>%
      filter(direction == "Up") %>%
      select(module, n_up = n)
    
    down_counts <- dir_counts %>%
      filter(direction == "Down") %>%
      select(module, n_down = n)
    
    module_sizes %>%
      left_join(up_counts, by = "module") %>%
      left_join(down_counts, by = "module") %>%
      mutate(
        n_up = ifelse(is.na(n_up), 0L, n_up),
        n_down = ifelse(is.na(n_down), 0L, n_down),
        frac_up = ifelse(module_size > 0, n_up / module_size, NA_real_),
        frac_down = ifelse(module_size > 0, n_down / module_size, NA_real_),
        frac_up_plot = frac_up,
        frac_down_plot = -frac_down,
        comparison = comp_tag
      )
  }) %>% bind_rows()
  
  plot_df <- per_comp %>%
    select(module, module_size, comparison, frac_up_plot, frac_down_plot, n_up, n_down) %>%
    tidyr::pivot_longer(
      cols = c(frac_down_plot, frac_up_plot),
      names_to = "direction",
      values_to = "frac_plot"
    ) %>%
    mutate(direction = ifelse(direction == "frac_up_plot", "Up", "Down"))
  
  module_order <- plot_df %>%
    group_by(module) %>%
    summarise(mean_abs = mean(abs(frac_plot), na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_abs) %>%
    pull(module)
  
  plot_df <- plot_df %>%
    mutate(
      module = factor(module, levels = module_order),
      comparison = factor(
        comparison,
        levels = vapply(comparisons, function(x) paste0(x[[1]], "_vs_", x[[2]]), character(1))
      )
    )
  
  max_abs <- max(abs(plot_df$frac_plot), na.rm = TRUE)
  max_abs <- ifelse(is.finite(max_abs), max_abs, 0)
  max_abs <- max(max_abs, 0.01)
  
  p <- ggplot(plot_df, aes(x = frac_plot, y = module, fill = direction)) +
    geom_col(width = 0.9) +
    facet_wrap(~ comparison, nrow = 1, scales = "fixed") +
    scale_x_continuous(
      limits = c(-max_abs, max_abs),
      labels = percent_format(accuracy = 1)
    ) +
    scale_fill_manual(values = c("Up" = "red3", "Down" = "dodgerblue3")) +
    labs(
      title = paste0(cell_type, ": fraction of module genes that are DEGs (Down vs Up)"),
      x = "Fraction of module genes that are DEGs (Down left, Up right)",
      y = "Module"
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(size = 9),
      strip.text = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # ---- UPDATED HEIGHT RULE (more compact for small module counts) ----
  # This is tuned so:
  #   n_mod = 6  -> ~3.3 inches
  #   n_mod = 11 -> ~4.2 inches
  n_mod <- length(unique(plot_df$module))
  pdf_h <- 2.2 + 0.18 * n_mod
  pdf_h <- min(max(pdf_h, 3.0), 8.0)  # clamp: never too short or too tall
  pdf_w <- 12
  
  out_pdf <- file.path(out_dir, paste0("hdWGCNA_fracDEG_byModule_diverging_", cell_type, ".pdf"))
  ggsave(out_pdf, p, width = pdf_w, height = pdf_h)
  
  out_tsv <- file.path(out_dir, paste0("hdWGCNA_fracDEG_byModule_diverging_", cell_type, ".tsv"))
  write.table(per_comp, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Wrote: ", out_pdf, " (modules=", n_mod, ", height=", round(pdf_h, 2), "in)")
}
