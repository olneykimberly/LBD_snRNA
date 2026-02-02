source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))

main_output_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/DEGs_RNA_pct0.25"
output_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/correlation_foldchange"

# DEG Filtering Parameters
up_logfc_cutoff <- 0.25
down_logfc_cutoff <- -0.25
padj_cutoff <- 0.05

# List of comparison pairs
comparisons_to_read <- list(
  c("AD_AT",    "CONTROL"),
  c("LBD_S",    "CONTROL"),
  c("LBD_AS",  "CONTROL"),
  c("LBD_ATS", "CONTROL")
)
comparison_names <- sapply(comparisons_to_read, function(p) paste0(p[1], "_vs_", p[2]))

# List of cell types to process
cell_types <- c("neuron", "interneuron", "oligodendrocyte",
                "opc", "astrocyte", "microglia", "fibroblast", "endothelial", "mural")

# Define the user-requested order for vertical facets (Top to Bottom)
new_plot_labels_order <- c("LBD_ATS", "LBD_AS", "LBD_S", "AD_AT")

# Map short names to full comparison names, ensuring the order is maintained
plot_order <- sapply(new_plot_labels_order, function(short_name) {
  if (short_name == "AD_AT") {
    return("AD_AT_vs_CONTROL")
  } else if (short_name == "LBD_S") {
    return("LBD_S_vs_CONTROL")
  } else if (short_name == "LBD_AS") {
    return("LBD_AS_vs_CONTROL")
  } else if (short_name == "LBD_ATS") {
    return("LBD_ATS_vs_CONTROL")
  }
})

# Create the map and factor levels based on the new order
comparison_map <- data.frame(
  comparison = plot_order,
  plot_label = new_plot_labels_order,
  stringsAsFactors = FALSE
)
plot_labels <- new_plot_labels_order

# ------------------------------
# Main Execution: Read files into a list of dataframes 
# ------------------------------
all_degs <- list()
for (ct in cell_types) {
  all_degs[[ct]] <- list()
  for (pair in comparisons_to_read) {
    comp_nm <- paste0(pair[1], "_vs_", pair[2])
    df <- read_deg_file(ct, pair[1], pair[2], main_output_dir)
    if (!is.null(df)) {
      all_degs[[ct]][[comp_nm]] <- df
    } else {
      all_degs[[ct]][[comp_nm]] <- tibble()
    }
  }
}

# --- Correlation Plot Generation ---
plot_celltype_correlations <- function(celltype_degs, ct_name,
                                       p_adj_cut, comp_map, output_dir) {
  message("\nProcessing correlations for cell type: ", ct_name)
  filtered_degs <- lapply(celltype_degs, function(df) {
    if (nrow(df) == 0) return(NULL)
    df %>%
      filter(p_val_adj <= 1) %>%
      select(gene, avg_log2FC)
  })
  filtered_degs <- filtered_degs[!sapply(filtered_degs, is.null)]
  if (length(filtered_degs) < 2) {
    warning("Not enough comparison dataframes (need at least 2) for cell type ", ct_name, ". Skipping plot.")
    return(NULL)
  }
  #  Merge all log2FC columns for this cell type
  comparison_names_present <- names(filtered_degs)
  # Rename log2FC columns with comparison names
  data_list <- lapply(comparison_names_present, function(comp) {
    filtered_degs[[comp]] %>% rename(!!comp := avg_log2FC)
  })
  # Sequentially join all by 'gene'
  merged_df <- Reduce(function(x, y) full_join(x, y, by = "gene"), data_list)
  # Impute NA log2FC values with 0 (gene not significant/measured in that comparison)
  merged_df[is.na(merged_df)] <- 0
  
  # 3. Prepare data for pairwise plotting (Unique Combinations)
  comp_pairs <- t(combn(comparison_names_present, 2))
  plot_data_list <- list()
  for (i in 1:nrow(comp_pairs)) {
    comp_x <- comp_pairs[i, 1]
    comp_y <- comp_pairs[i, 2]
    x_label_short <- comp_map$plot_label[comp_map$comparison == comp_x]
    y_label_short <- comp_map$plot_label[comp_map$comparison == comp_y]
    facet_label <- paste0(y_label_short, " vs ", x_label_short)
    pair_df <- merged_df %>%
      select(gene, x_log2FC = !!comp_x, y_log2FC = !!comp_y) %>%
      mutate(
        x_comp = comp_x,
        y_comp = comp_y,
        facet_combination = facet_label
      )
    plot_data_list[[i]] <- pair_df
  }
  
  plot_data <- bind_rows(plot_data_list)
  facet_levels <- unique(plot_data$facet_combination)
  plot_data$facet_combination <- factor(plot_data$facet_combination, levels = facet_levels)
  max_abs_lfc <- max(abs(c(plot_data$x_log2FC, plot_data$y_log2FC)))
  max_abs_lfc_buffered <- max_abs_lfc * 1.05
  
  # Generate Plot
  cor_df <- plot_data %>%
    group_by(facet_combination) %>%
    summarise(
      cor_val = cor(x_log2FC, y_log2FC, method = "pearson"),
      cor_label = sprintf("r = %.2f", cor_val),
      .groups = "drop"
    ) %>%
    mutate(
      x_pos = -max_abs_lfc_buffered * 0.9,
      y_pos = max_abs_lfc_buffered * 0.9
    )
  
  p <- ggplot(plot_data, aes(x = x_log2FC, y = y_log2FC)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    facet_wrap(~ facet_combination, scales = "fixed") +
    geom_text(
      data = cor_df,
      mapping = aes(x = x_pos, y = y_pos, label = cor_label),
      hjust = 0, vjust = 1,
      size = 3.5 
    ) +
    coord_fixed() + 
    xlim(-max_abs_lfc_buffered, max_abs_lfc_buffered) +
    ylim(-max_abs_lfc_buffered, max_abs_lfc_buffered) +
    labs(
      title = paste0("Log2FC correlation (adjusted p < ", p_adj_cut, "): ", str_to_title(ct_name)),
      x = "Log2FC",
      y = "Log2FC"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = "gray50"),
      strip.text = element_text(face = "bold")
    )
  dir.create(output_dir, showWarnings = FALSE)
  plot_file <- file.path(output_dir, paste0(ct_name, "_log2FC_correlations_optimized.pdf"))
  ggsave(plot_file, plot = p, width = 7, height = 5) 
  message("  Optimized plot saved to: ", plot_file)
  return(p)
}

## Run the function for all cell types
all_plots <- list()
for (ct in cell_types) {
  ct_degs <- all_degs[[ct]]
  
  if (length(ct_degs) > 0) {
    p_plot <- plot_celltype_correlations(
      celltype_degs = ct_degs,
      ct_name = ct,
      p_adj_cut = padj_cutoff,
      comp_map = comparison_map,
      output_dir = output_dir
    )
    if (!is.null(p_plot)) {
      all_plots[[ct]] <- p_plot
    }
  }
}

message("\nAll optimized correlation plots have been generated and saved.")