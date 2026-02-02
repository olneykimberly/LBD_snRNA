#!/usr/bin/env Rscript
# ==============================================================================
# Sex-pattern consistency across disease-vs-control comparisons
# INCLUDING SEX-UNIQUE DEGs (simple definition; no extra padj guard)
#
# For each cell type, make FOUR heatmaps:
#   1) DEGs with shared expression between the sexes
#   2) DEGs with opposite expression between the sexes
#   3) DEGs uniquely dysregulated in females
#   4) DEGs uniquely dysregulated in males
#
# Also writes a small per-cell-type summary TSV with counts of genes
# consistent in >= 3 comparisons for each category.
#
# Consistency rule:
#   - A gene must fall into the SAME category in >= 3 disease-vs-control comparisons.
#
# Membership rules (per disease-vs-control comparison):
#   female_deg = padj<=0.05 & |log2FC|>=0.25 in FEMALE
#   male_deg   = padj<=0.05 & |log2FC|>=0.25 in MALE
#
# Categories:
#   shared_expression   : female_deg & male_deg & same direction
#   opposite_expression : female_deg & male_deg & opposite direction
#   female_unique       : female_deg & !male_deg
#   male_unique         : male_deg & !female_deg
#
# Heatmap display:
#   - Show ALL log2FC values (even when not significant)
#   - Columns: disease names only (AD AT, LBD S, LBD AS, LBD ATS)
#   - Facet band: Female then Male
#   - Colors: red=positive, blue=negative
#   - Gene order: by row-mean log2FC (largest positive at top -> most negative at bottom)
#   - Scale rules:
#       * shared_expression, female_unique, male_unique: fixed PER CELL TYPE (same scale for those 3 within ct)
#       * opposite_expression: scale computed from THAT heatmap only
#
# Outputs:
#   {output_dir}/{celltype}/
#     - {ct}_shared_expression_heatmap.pdf
#     - {ct}_opposite_expression_heatmap.pdf
#     - {ct}_female_unique_heatmap.pdf
#     - {ct}_male_unique_heatmap.pdf
#     - {ct}_sex_pattern_tables.tsv
#     - {ct}_sex_pattern_summary.tsv   <-- NEW
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(scales)
})

# ------------------------------
# Paths
# ------------------------------
main_output_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/DEGs_RNA_pct0.25_by_sex"
output_dir      <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/Heatmaps/sex_pattern_consistency_with_unique"

# ------------------------------
# Parameters
# ------------------------------
logfc_cutoff <- 0.25
padj_cutoff  <- 0.05
min_comparisons_required <- 3

cell_types <- c("neuron", "interneuron", "oligodendrocyte",
                "opc", "astrocyte", "microglia")

comparisons_base <- list(
  c("AD_AT",   "CONTROL"),
  c("LBD_S",   "CONTROL"),
  c("LBD_AS",  "CONTROL"),
  c("LBD_ATS", "CONTROL")
)

# Drop "vs CONTROL" in column labels -> use disease only
disease_order <- vapply(comparisons_base, function(p) p[1], character(1))
disease_labels_pretty <- gsub("_", " ", disease_order)

# ------------------------------
# Helpers
# ------------------------------
ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

deg_flag <- function(padj, lfc) {
  !is.na(padj) & !is.na(lfc) & (padj <= padj_cutoff) & (abs(lfc) >= logfc_cutoff)
}

deg_dir <- function(lfc) {
  dplyr::case_when(
    is.na(lfc) ~ NA_character_,
    lfc > 0 ~ "Up",
    lfc < 0 ~ "Down",
    TRUE ~ NA_character_
  )
}

read_deg_file <- function(ct, g1, g2) {
  f <- file.path(main_output_dir, ct, sprintf("DEG_%s_%s_vs_%s.tsv", ct, g1, g2))
  if (!file.exists(f)) return(NULL)
  
  df <- readr::read_tsv(f, show_col_types = FALSE)
  
  required <- c("gene", "avg_log2FC", "p_val_adj")
  if (!all(required %in% colnames(df))) {
    stop("Missing required columns in: ", f, " | need: ", paste(required, collapse = ", "))
  }
  
  df %>%
    dplyr::select(gene, avg_log2FC, p_val_adj) %>%
    dplyr::mutate(gene = as.character(gene))
}

plot_heatmap <- function(df_long, title, outfile, max_abs) {
  if (nrow(df_long) == 0) {
    message("No rows to plot: ", outfile)
    return(invisible(NULL))
  }
  if (!is.finite(max_abs) || max_abs == 0) {
    message("Non-finite/zero max_abs; skipping: ", outfile)
    return(invisible(NULL))
  }
  
  n_genes <- dplyr::n_distinct(df_long$gene)
  h <- min(max(3.5, 0.12 * n_genes + 1.8), 18)
  
  p <- ggplot(df_long, aes(disease, gene, fill = value)) +
    geom_tile() +
    facet_grid(. ~ sex, switch = "x") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      limits = c(-max_abs, max_abs),
      oob = scales::squish,
      na.value = "grey90"
    ) +
    labs(title = title, x = NULL, y = NULL, fill = "log2FC") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text.x = element_text(size = 9, face = "bold"),
      strip.placement = "outside",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
  
  pdf(outfile, width = 9, height = h)
  print(p)
  dev.off()
}

make_long_for_pattern <- function(df_ct, consistent_tbl, pattern_label, disease_order, disease_labels_pretty) {
  genes <- consistent_tbl %>%
    dplyr::filter(sex_pattern == pattern_label) %>%
    dplyr::pull(gene) %>%
    unique()
  
  if (length(genes) == 0) return(tibble())
  
  long <- df_ct %>%
    dplyr::filter(gene %in% genes) %>%
    dplyr::select(gene, disease, female_log2FC, male_log2FC) %>%
    tidyr::pivot_longer(
      cols = c(female_log2FC, male_log2FC),
      names_to = "sex_field",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      sex = dplyr::case_when(
        sex_field == "female_log2FC" ~ "Female",
        sex_field == "male_log2FC"   ~ "Male",
        TRUE ~ NA_character_
      ),
      sex = factor(sex, levels = c("Female", "Male")),
      disease = factor(disease, levels = disease_order, labels = disease_labels_pretty)
    ) %>%
    dplyr::select(gene, disease, sex, value)
  
  # Order genes by row-mean log2FC (largest -> smallest)
  gene_order <- long %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(row_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(row_mean), gene) %>%
    dplyr::pull(gene)
  
  long %>%
    dplyr::mutate(gene = factor(gene, levels = rev(gene_order)))
}

# ------------------------------
# Main
# ------------------------------
ensure_dir(output_dir)

for (ct in cell_types) {
  message("\n============================================================")
  message("Cell type: ", ct)
  message("============================================================")
  
  rows <- list()
  
  for (cmp in comparisons_base) {
    dis  <- cmp[1]
    ctrl <- cmp[2]
    
    df_f <- read_deg_file(ct, paste0(dis, "_female"), paste0(ctrl, "_female"))
    df_m <- read_deg_file(ct, paste0(dis, "_male"),   paste0(ctrl, "_male"))
    
    if (is.null(df_f) && is.null(df_m)) next
    
    if (is.null(df_f)) df_f <- tibble(gene = character(), avg_log2FC = numeric(), p_val_adj = numeric())
    if (is.null(df_m)) df_m <- tibble(gene = character(), avg_log2FC = numeric(), p_val_adj = numeric())
    
    df_f <- df_f %>% dplyr::rename(female_log2FC = avg_log2FC, female_padj = p_val_adj)
    df_m <- df_m %>% dplyr::rename(male_log2FC   = avg_log2FC, male_padj   = p_val_adj)
    
    rows[[length(rows) + 1]] <-
      dplyr::full_join(df_f, df_m, by = "gene") %>%
      dplyr::mutate(
        cell_type = ct,
        disease   = dis,
        
        female_deg = deg_flag(female_padj, female_log2FC),
        male_deg   = deg_flag(male_padj,   male_log2FC),
        
        female_dir = dplyr::if_else(female_deg, deg_dir(female_log2FC), NA_character_),
        male_dir   = dplyr::if_else(male_deg,   deg_dir(male_log2FC),   NA_character_),
        
        female_unique = female_deg & !male_deg,
        male_unique   = male_deg & !female_deg,
        
        sex_pattern = dplyr::case_when(
          female_deg & male_deg & female_dir == male_dir ~ "shared_expression",
          female_deg & male_deg & female_dir != male_dir ~ "opposite_expression",
          female_unique ~ "female_unique",
          male_unique   ~ "male_unique",
          TRUE ~ NA_character_
        )
      )
  }
  
  df_ct <- dplyr::bind_rows(rows)
  
  if (nrow(df_ct) == 0) {
    message("  No DEG tables found for: ", ct)
    next
  }
  
  ct_dir <- file.path(output_dir, ct)
  ensure_dir(ct_dir)
  
  # Save audit table
  readr::write_tsv(df_ct, file.path(ct_dir, paste0(ct, "_sex_pattern_tables.tsv")))
  
  # Require same category in >=3 comparisons
  consistent <- df_ct %>%
    dplyr::filter(!is.na(sex_pattern)) %>%
    dplyr::count(gene, sex_pattern, name = "n_comp") %>%
    dplyr::filter(n_comp >= min_comparisons_required)
  
  # ------------------------------
  # NEW: per-cell-type summary TSV
  # ------------------------------
  summary_tbl <- consistent %>%
    dplyr::group_by(sex_pattern) %>%
    dplyr::summarise(n_genes = dplyr::n_distinct(gene), .groups = "drop") %>%
    dplyr::right_join(
      tibble(sex_pattern = c("shared_expression", "opposite_expression", "female_unique", "male_unique")),
      by = "sex_pattern"
    ) %>%
    dplyr::mutate(
      n_genes = tidyr::replace_na(n_genes, 0L),
      cell_type = ct,
      min_comparisons_required = min_comparisons_required
    ) %>%
    dplyr::select(cell_type, sex_pattern, n_genes, min_comparisons_required) %>%
    dplyr::arrange(factor(sex_pattern, levels = c("shared_expression", "opposite_expression", "female_unique", "male_unique")))
  
  readr::write_tsv(summary_tbl, file.path(ct_dir, paste0(ct, "_sex_pattern_summary.tsv")))
  
  # ------------------------------
  # Build long dfs for heatmaps
  # ------------------------------
  hm_shared <- make_long_for_pattern(df_ct, consistent, "shared_expression",
                                     disease_order, disease_labels_pretty)
  hm_opp    <- make_long_for_pattern(df_ct, consistent, "opposite_expression",
                                     disease_order, disease_labels_pretty)
  hm_funiq  <- make_long_for_pattern(df_ct, consistent, "female_unique",
                                     disease_order, disease_labels_pretty)
  hm_muniq  <- make_long_for_pattern(df_ct, consistent, "male_unique",
                                     disease_order, disease_labels_pretty)
  
  # Scales:
  ct_max_abs <- max(abs(c(df_ct$female_log2FC, df_ct$male_log2FC)), na.rm = TRUE)
  opp_max_abs <- max(abs(hm_opp$value), na.rm = TRUE)
  
  # Plot
  plot_heatmap(
    hm_shared,
    paste0(ct, "\nDEGs with shared expression between the sexes"),
    file.path(ct_dir, paste0(ct, "_shared_expression_heatmap.pdf")),
    ct_max_abs
  )
  
  plot_heatmap(
    hm_opp,
    paste0(ct, "\nDEGs with opposite expression between the sexes"),
    file.path(ct_dir, paste0(ct, "_opposite_expression_heatmap.pdf")),
    opp_max_abs
  )
  
  plot_heatmap(
    hm_funiq,
    paste0(ct, "\nDEGs uniquely dysregulated in females"),
    file.path(ct_dir, paste0(ct, "_female_unique_heatmap.pdf")),
    ct_max_abs
  )
  
  plot_heatmap(
    hm_muniq,
    paste0(ct, "\nDEGs uniquely dysregulated in males"),
    file.path(ct_dir, paste0(ct, "_male_unique_heatmap.pdf")),
    ct_max_abs
  )
  
  message("  Wrote heatmaps + summary for: ", ct)
}

message("\nDONE: outputs written to ", output_dir)
