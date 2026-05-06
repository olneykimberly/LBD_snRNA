#!/usr/bin/env Rscript
# ==============================================================================
# hdWGCNA module eigengene vs DISEASE GROUP plots donor-level
# - Boxplot by disease group + donor points
# - Points shaped by sex: female = circle, male = square
# - No p-values plotted
# - Removes NA/unused groups
# - Legend removed
# - Output size: width = 2.25, height = 1.85
# - Outputs Wilcoxon p-values for CONTROL vs each disease group
# ==============================================================================

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE, character.only = TRUE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

theme_set(cowplot::theme_cowplot(font_size = 7))

# ----------------------------- paths ------------------------------------------
obj_dir      <- "../rObjects"
base_out_dir <- "../results/hdWGCNA"

out_dir <- file.path(base_out_dir, "ME_boxplot_disease_no_pvalues")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- settings ---------------------------------------
cell_types <- c("neuron",
  "interneuron",
  "oligodendrocyte",
  "opc",
  "astrocyte",
  "microglia",
  "endothelial",
  "fibroblast",
  "mural")

min_donors_total     <- 8
jitter_width_default <- 0.18

group_levels_raw <- c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS")

group_labels_map <- c(
  CONTROL = "CONTROL",
  AD_AT   = "AD (AT)",
  LBD_S   = "LBD (S)",
  LBD_AS  = "LBD (AS)",
  LBD_ATS = "LBD (ATS)"
)

group_levels_label <- unname(group_labels_map[group_levels_raw])

disease_groups <- c("AD_AT", "LBD_S", "LBD_AS", "LBD_ATS")

all_pval_tbl <- list()

# ----------------------------- helper functions -------------------------------
wrap_genes_3_per_line <- function(genes) {
  if (length(genes) == 0) return("(hub genes not found)")
  idx <- ceiling(seq_along(genes) / 3)
  paste(tapply(genes, idx, paste, collapse = "  "), collapse = "\n")
}

is_blackish <- function(col) {
  if (is.na(col) || length(col) == 0) return(FALSE)
  if (tolower(col) == "black") return(TRUE)
  
  rgb <- tryCatch(col2rgb(col), error = function(e) NULL)
  if (is.null(rgb)) return(FALSE)
  
  (0.2126 * rgb[1] + 0.7152 * rgb[2] + 0.0722 * rgb[3]) / 255 < 0.08
}

get_sex_from_sample <- function(x) {
  dplyr::case_when(
    grepl("_F", as.character(x)) ~ "Female",
    grepl("_M", as.character(x)) ~ "Male",
    TRUE ~ NA_character_
  )
}

make_control_vs_disease_pvals <- function(df_g, ct, m) {
  
  bind_rows(lapply(disease_groups, function(g) {
    
    test_dat <- df_g %>%
      filter(group_raw %in% c("CONTROL", g)) %>%
      mutate(group_raw = factor(group_raw, levels = c("CONTROL", g))) %>%
      droplevels()
    
    n_control <- sum(test_dat$group_raw == "CONTROL", na.rm = TRUE)
    n_disease <- sum(test_dat$group_raw == g, na.rm = TRUE)
    
    if (
      n_distinct(test_dat$group_raw) < 2 ||
      n_control < 2 ||
      n_disease < 2
    ) {
      return(data.frame(
        cell_type = ct,
        module = m,
        disease_group = unname(group_labels_map[g]),
        comparison = paste0(g, "_vs_CONTROL"),
        n_control = n_control,
        n_disease = n_disease,
        p_value = NA_real_
      ))
    }
    
    wt <- wilcox.test(value ~ group_raw, data = test_dat, exact = FALSE)
    
    data.frame(
      cell_type = ct,
      module = m,
      disease_group = unname(group_labels_map[g]),
      comparison = paste0(g, "_vs_CONTROL"),
      n_control = n_control,
      n_disease = n_disease,
      p_value = wt$p.value
    )
  })) %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH")
    ) %>%
    select(
      cell_type,
      module,
      disease_group,
      comparison,
      n_control,
      n_disease,
      p_value,
      p_adj_BH
    )
}

# ==============================================================================
# MAIN LOOP
# ==============================================================================
for (ct in cell_types) {
  
  message("---- ", ct, " ----")
  
  seu_path <- file.path(obj_dir, paste0("hdWGCNA_", ct, "_final.qs"))
  
  if (!file.exists(seu_path)) {
    warning("Missing object: ", seu_path)
    next
  }
  
  seu <- qread(seu_path)
  md  <- seu[[]]
  md$cell_id <- rownames(md)
  
  if (!("Sample_ID" %in% colnames(md))) stop(ct, ": missing Sample_ID")
  if (!("group" %in% colnames(md)))     stop(ct, ": missing group")
  
  donor_group <- md %>%
    mutate(group = as.character(group)) %>%
    filter(!is.na(Sample_ID)) %>%
    filter(!is.na(group)) %>%
    filter(group %in% group_levels_raw) %>%
    group_by(Sample_ID) %>%
    summarise(
      group = first(group),
      sex = get_sex_from_sample(first(Sample_ID)),
      .groups = "drop"
    )
  
  mods_raw <- GetModules(seu) %>%
    distinct(module, color) %>%
    filter(!is.na(module)) %>%
    filter(module != "grey")
  
  MEs <- as.data.frame(GetMEs(seu, harmonized = TRUE))
  MEs$cell_id <- rownames(MEs)
  
  hub_file <- file.path(base_out_dir, paste0("hubs_", ct, ".tsv"))
  hubs_tbl <- if (file.exists(hub_file)) {
    read.delim(hub_file, stringsAsFactors = FALSE)
  } else {
    NULL
  }
  
  ME_long <- MEs %>%
    pivot_longer(
      cols = -cell_id,
      names_to = "module",
      values_to = "value"
    ) %>%
    filter(!is.na(module)) %>%
    filter(module != "grey") %>%
    filter(!is.na(value)) %>%
    inner_join(md[, c("cell_id", "Sample_ID")], by = "cell_id")
  
  donor_me <- ME_long %>%
    group_by(module, Sample_ID) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  
  donor_df <- donor_me %>%
    left_join(donor_group, by = "Sample_ID") %>%
    left_join(mods_raw, by = "module") %>%
    mutate(group = as.character(group)) %>%
    filter(!is.na(group)) %>%
    filter(group %in% group_levels_raw) %>%
    filter(!is.na(value)) %>%
    filter(!is.na(module))
  
  for (m in unique(donor_df$module)) {
    
    if (m == "grey" || is.na(m)) next
    
    df_m <- donor_df %>%
      filter(module == m)
    
    mod_col <- df_m$color[which(!is.na(df_m$color))[1]]
    if (is.na(mod_col) || length(mod_col) == 0) mod_col <- "grey70"
    
    txt_col <- if (is_blackish(mod_col)) "white" else "black"
    
    hubs_text <- "(hub genes not found)"
    
    if (!is.null(hubs_tbl) && all(c("module", "gene_name") %in% colnames(hubs_tbl))) {
      hubs <- hubs_tbl %>%
        filter(module == m) %>%
        pull(gene_name) %>%
        unique() %>%
        head(9)
      
      if (length(hubs) > 0) hubs_text <- wrap_genes_3_per_line(hubs)
    }
    
    df_g <- df_m %>%
      mutate(
        group_raw = as.character(group),
        group_label = unname(group_labels_map[group_raw]),
        group_label = factor(group_label, levels = group_levels_label),
        sex = factor(sex, levels = c("Female", "Male"))
      ) %>%
      filter(!is.na(group_raw)) %>%
      filter(group_raw %in% group_levels_raw) %>%
      filter(!is.na(group_label)) %>%
      droplevels()
    
    if (nrow(df_g) < min_donors_total) next
    if (nlevels(droplevels(df_g$group_label)) < 2) next
    
    pval_tbl <- make_control_vs_disease_pvals(df_g, ct, m)
    all_pval_tbl[[paste(ct, m, sep = "_")]] <- pval_tbl
    
    lowN_df <- df_g %>%
      group_by(group_label) %>%
      summarise(
        n_bin = n(),
        med = median(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n_bin < 3)
    
    p_main <- ggplot(df_g, aes(x = group_label, y = value)) +
      geom_boxplot(
        fill = mod_col,
        outlier.shape = NA,
        width = 0.65
      ) +
      geom_crossbar(
        data = lowN_df,
        aes(x = group_label, y = med, ymin = med, ymax = med),
        inherit.aes = FALSE,
        width = 0.50,
        linewidth = 0.8,
        color = "black"
      ) +
      geom_point(
        aes(shape = sex),
        position = position_jitter(width = jitter_width_default, height = 0),
        size = 1.35,
        alpha = 0.9,
        color = "black",
        stroke = 0.4
      ) +
      scale_shape_manual(
        values = c(
          Female = 16,
          Male = 15
        ),
        na.translate = FALSE
      ) +
      scale_x_discrete(
        limits = group_levels_label,
        drop = TRUE,
        na.translate = FALSE
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0.05, 0.12))
      ) +
      labs(
        title = paste(m),
        x = NULL,
        y = "Module eigengene"
      ) +
      theme(
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        plot.margin = margin(4, 6, 4, 6)
      ) +
      coord_cartesian(clip = "off")
    
    p_bubble <- ggplot(
      data.frame(x = 1, y = 1, label = hubs_text),
      aes(x, y, label = label)
    ) +
      geom_label(
        fill = mod_col,
        color = txt_col,
        label.size = 0,
        size = 1.75,
        lineheight = 0.95,
        label.padding = unit(0.45, "lines"),
        label.r = unit(0.8, "lines")
      ) +
      theme_void() +
      theme(plot.margin = margin(0, 10, 1, 10))
    
    final <- cowplot::plot_grid(
      p_main,
      p_bubble,
      ncol = 1,
      rel_heights = c(1, 0.22)
    )
    
    ggsave(
      filename = file.path(out_dir, paste0("ME_group_", m, ".pdf")),
      plot = final,
      width = 1.8,
      height = 1.55,
      device = cairo_pdf
    )
  }
  
  rm(seu)
  invisible(gc())
}

# ==============================================================================
# Save combined p-value table
# ==============================================================================
combined_pval_tbl <- bind_rows(all_pval_tbl)

write.table(
  combined_pval_tbl,
  file = file.path(out_dir, "ME_CONTROL_vs_disease_pvalues_all_modules.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("DONE.")