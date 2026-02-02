#!/usr/bin/env Rscript
# ==============================================================================
# UpSet plots + DEG count bars for LBD types vs AD_AT, stratified by sex
#
# For each SEX (female, male) and each CELL TYPE:
#   A) DEG count bar plot (Up vs Down) for LBD_S/LBD_AS/LBD_ATS vs AD_AT
#   B) UpSet for DOWN-regulated genes
#   C) UpSet for UP-regulated genes
#
# Input DEG files expected at:
#   {main_output_dir}/{celltype}/DEG_{celltype}_{group1}_vs_{group2}.tsv
#
# Example:
#   DEG_microglia_LBD_S_female_vs_AD_AT_female.tsv
#   DEG_microglia_LBD_S_male_vs_AD_AT_male.tsv
#
# Outputs written under:
#   {output_dir}/female/
#   {output_dir}/male/
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexUpset)
  library(readr)
  library(openxlsx)
  library(cowplot)
  library(ggplot2)
})

# ------------------------------
# Directories and Parameters
# ------------------------------
main_output_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/DEGs_RNA_pct0.25_by_sex"
output_dir      <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/UpSet/LBDtypes_vs_AD_by_sex"

# DEG Filtering Parameters
up_logfc_cutoff   <- 0.25
down_logfc_cutoff <- -0.25
padj_cutoff       <- 0.05

# Cell types to process
cell_types <- c("neuron", "interneuron", "oligodendrocyte",
                "opc", "astrocyte", "microglia", "mural", "endothelial")

# Desired order (Top to Bottom) for the LBD facets
new_plot_labels_order <- c("LBD_ATS", "LBD_AS", "LBD_S")

# Sexes to run
sexes <- c("female", "male")

# ------------------------------
# Helpers
# ------------------------------
ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------
# Function to read one DEG file
# ------------------------------
read_deg_file <- function(celltype, group1, group2, main_output_dir) {
  fname <- file.path(
    main_output_dir,
    celltype,
    sprintf("DEG_%s_%s_vs_%s.tsv", celltype, group1, group2)
  )
  
  if (!file.exists(fname)) {
    warning("DEG file not found: ", fname)
    return(NULL)
  }
  
  df <- read_tsv(fname, show_col_types = FALSE)
  
  required <- c("gene", "avg_log2FC", "p_val_adj")
  if (!all(required %in% colnames(df))) {
    stop("DEG file ", fname, " missing required columns: ", paste(required, collapse = ", "))
  }
  
  df %>% select(all_of(required))
}

# ------------------------------
# Utility: fromList (binary membership matrix)
# ------------------------------
fromList <- function(input) {
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = FALSE))
  data <- data[which(rowSums(data) != 0), , drop = FALSE]
  names(data) <- names(input)
  row.names(data) <- elements
  data
}

# ------------------------------
# UpSet plot function (3-set version)
# ------------------------------
generate_upset_plot <- function(data, celltype, sex, reg_status, set_names, dot_color) {
  if (nrow(data) == 0) {
    return(
      ggplot() +
        geom_text(
          aes(x = 0.5, y = 0.5, label = paste(celltype, sex, reg_status, "\nNo DEGs found")),
          size = 4
        ) +
        theme_void()
    )
  }
  
  # Ensure columns match set_names in the provided order
  data <- data[, set_names, drop = FALSE]
  
  upset_queries <- lapply(set_names, function(set_name) {
    upset_query(set = set_name, fill = dot_color)
  })
  
  # Intersections list for exactly 3 sets (in provided order)
  sets_list <- list(
    set_names[1], set_names[2], set_names[3],
    c(set_names[1], set_names[2], set_names[3]),
    c(set_names[3], set_names[2]),
    c(set_names[1], set_names[3]),
    c(set_names[1], set_names[2])
  )
  
  plot_title <- paste(celltype, sex, reg_status)
  
  upset(
    data,
    set_names,
    set_sizes = FALSE,
    themes = upset_modify_themes(
      list(
        "intersections_matrix" = theme(
          text = element_text(size = 9),
          plot.margin = margin(0, 0, 0, 0, "cm")
        ),
        "overall_sizes" = theme(
          axis.text.x = element_text(size = 9),
          plot.margin = margin(0, 0, 0, 0, "cm")
        )
      )
    ),
    queries = upset_queries,
    intersections = sets_list,
    base_annotations = list(
      "Intersection size" = (
        intersection_size(
          size = 2,
          text = list(size = 2.5),
          text_mapping = aes(),
          bar_number_threshold = 3,
          width = 0.3,
          mapping = aes(fill = "bars_color")
        ) +
          scale_fill_manual(values = c("grey"), guide = "none") +
          scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
          theme(
            axis.text.y = element_text(size = 8),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(0, 0, 0, 0.1, "cm"),
            axis.line = element_line(colour = "black")
          )
      )
    ),
    matrix = intersection_matrix(
      geom = geom_point(
        shape = "circle filled",
        size = 2,
        stroke = 0.45,
        color = "black"
      )
    ) +
      theme(
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.margin = margin(0, 0, 0, 0, "cm")
      ),
    sort_sets = FALSE,
    sort_intersections = FALSE,
    name = plot_title
  )
}

# ------------------------------
# DEG count bar plot for ONE cell type + ONE sex
# ------------------------------
generate_deg_count_plot <- function(deg_counts_df, current_celltype, comparison_order) {
  celltype_data <- deg_counts_df %>%
    filter(cell_type == current_celltype) %>%
    mutate(plot_label = factor(plot_label, levels = comparison_order))
  
  if (nrow(celltype_data) == 0 || sum(celltype_data$n) == 0) {
    return(
      ggplot() +
        geom_text(
          aes(x = 0.5, y = 0.5, label = paste(current_celltype, "DEGs: No DEGs found")),
          size = 3
        ) +
        theme_void()
    )
  }
  
  max_count <- max(celltype_data$n, na.rm = TRUE)
  max_abs_x <- ceiling(max_count * 1.5 / 10) * 10
  
  ggplot(celltype_data, aes(x = value, y = cell_type, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(
      data = celltype_data %>% filter(n > 0),
      aes(y = cell_type, label = as.character(n), x = value),
      hjust = ifelse((celltype_data %>% filter(n > 0))$value < 0, 1.1, -0.1),
      size = 2.5,
      color = "black"
    ) +
    scale_fill_manual(values = c("Down" = "blue", "Up" = "red")) +
    scale_x_continuous(
      labels = function(x) abs(x),
      limits = c(-max_abs_x, max_abs_x),
      breaks = c(-max_abs_x, 0, max_abs_x)
    ) +
    facet_grid(plot_label ~ ., switch = "y") +
    labs(
      title = paste(current_celltype, "DEG Counts"),
      x = "Number of DEGs",
      y = NULL,
      fill = "Direction"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_text(size = 8),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5),
      panel.grid.major.x = element_line(linetype = "dashed", color = "lightgray"),
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "gray90", color = NA),
      strip.text = element_text(face = "bold", size = 7.5),
      strip.placement = "outside",
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box.margin = margin(-.2, -.2, -.2, -.2, "cm"),
      legend.key.size = unit(0.5, "lines"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      plot.margin = margin(0.2, 0.15, 0.1, 0.2, "cm")
    )
}

# ------------------------------
# Main loop: run per sex
# ------------------------------
for (sex in sexes) {
  message("\n============================================================")
  message("RUNNING SEX: ", sex)
  message("============================================================\n")
  
  sex_outdir <- file.path(output_dir, sex)
  ensure_dir(sex_outdir)
  
  # Comparisons for this sex (LBD vs AD_AT)
  comparisons_to_read <- list(
    c(paste0("LBD_S_", sex),   paste0("AD_AT_", sex)),
    c(paste0("LBD_AS_", sex),  paste0("AD_AT_", sex)),
    c(paste0("LBD_ATS_", sex), paste0("AD_AT_", sex))
  )
  comparison_names <- vapply(comparisons_to_read, function(p) paste0(p[1], "_vs_", p[2]), character(1))
  
  # The exact comparison strings used in filenames (includes sex)
  plot_order <- vapply(new_plot_labels_order, function(lbl) {
    paste0(lbl, "_", sex, "_vs_AD_AT_", sex)
  }, character(1))
  
  comparison_map <- data.frame(
    comparison = plot_order,
    plot_label = new_plot_labels_order,
    stringsAsFactors = FALSE
  )
  
  # UpSet set labels (no sex)
  set_names <- new_plot_labels_order  # c("LBD_ATS","LBD_AS","LBD_S")
  
  # ------------------------------------------------------------
  # Read DEG tables for this sex
  # ------------------------------------------------------------
  all_degs <- list()
  
  for (ct in cell_types) {
    message("Reading files for cell type: ", ct)
    all_degs[[ct]] <- list()
    
    for (pair in comparisons_to_read) {
      comp_nm <- paste0(pair[1], "_vs_", pair[2])
      df <- read_deg_file(ct, pair[1], pair[2], main_output_dir)
      
      if (!is.null(df)) {
        all_degs[[ct]][[comp_nm]] <- df
        message("  Successfully read: ", comp_nm, " (Rows: ", nrow(df), ")")
      } else {
        all_degs[[ct]][[comp_nm]] <- tibble()
      }
    }
  }
  message("\nFinished reading all DEG files for sex: ", sex)
  
  # ------------------------------------------------------------
  # Build DEG count table for this sex (only these 3 comparisons)
  # ------------------------------------------------------------
  safe_read_deg <- function(path) {
    tryCatch(
      read_tsv(path, show_col_types = FALSE),
      error = function(e) tryCatch(
        read_csv(path, show_col_types = FALSE),
        error = function(e2) read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE)
      )
    )
  }
  
  all_files <- list.files(main_output_dir,
                          pattern = "^DEG_.*_vs_.*\\.(txt|csv|tsv)$",
                          recursive = TRUE,
                          full.names = TRUE)
  
  if (length(all_files) == 0) stop("No DEG files found under main_output_dir.")
  
  file_df <- tibble(path = all_files) %>%
    mutate(bn = basename(path),
           noext = tools::file_path_sans_ext(bn),
           after_deg = str_remove(noext, "^DEG_")) %>%
    mutate(
      cell_type = word(after_deg, 1, sep = fixed("_")),
      rest = str_remove(after_deg, paste0("^", cell_type, "_")),
      comparison = rest
    ) %>%
    select(path, cell_type, comparison) %>%
    filter(cell_type %in% cell_types) %>%
    filter(comparison %in% plot_order)
  
  deg_long <- file_df %>%
    mutate(data = map(path, safe_read_deg)) %>%
    filter(map_lgl(data, ~ !is.null(.x))) %>%
    select(-path) %>%
    tidyr::unnest(data) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::mutate(
      p_val_adj  = dplyr::coalesce(.data$p_val_adj),
      avg_log2fc = dplyr::coalesce(.data$avg_log2fc)
    ) %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        !is.na(p_val_adj) & !is.na(avg_log2fc) & p_val_adj <= padj_cutoff & avg_log2fc >= up_logfc_cutoff   ~ "Up",
        !is.na(p_val_adj) & !is.na(avg_log2fc) & p_val_adj <= padj_cutoff & avg_log2fc <= down_logfc_cutoff ~ "Down",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(direction)) %>%
    dplyr::count(comparison, cell_type, direction, name = "n")
  
  all_combinations <- expand_grid(
    comparison = plot_order,
    cell_type = cell_types,
    direction = c("Down", "Up")
  )
  
  deg_counts_all <- all_combinations %>%
    left_join(deg_long, by = c("comparison", "cell_type", "direction")) %>%
    mutate(n = replace_na(n, 0)) %>%
    mutate(value = if_else(direction == "Down", -n, n)) %>%
    mutate(cell_type = factor(cell_type, levels = cell_types)) %>%
    left_join(comparison_map, by = "comparison") %>%
    mutate(plot_label = factor(plot_label, levels = new_plot_labels_order))
  
  message("\nFinished calculating DEG counts for sex: ", sex)
  
  # ------------------------------------------------------------
  # Make UpSet + Count plots per cell type and save
  # ------------------------------------------------------------
  for (ct in cell_types) {
    message("\nProcessing cell type: ", ct, " (", sex, ") ...")
    
    # Map set label -> comparison string (with sex)
    set_to_comp <- setNames(
      vapply(set_names, function(lbl) paste0(lbl, "_", sex, "_vs_AD_AT_", sex), character(1)),
      set_names
    )
    
    up_list_input <- list()
    down_list_input <- list()
    
    for (set_lab in set_names) {
      comp <- set_to_comp[[set_lab]]
      df <- all_degs[[ct]][[comp]]
      
      if (!is.null(df) && nrow(df) > 0) {
        up_list_input[[set_lab]] <- df %>%
          filter(avg_log2FC >= up_logfc_cutoff, p_val_adj <= padj_cutoff) %>%
          pull(gene) %>% as.character()
        
        down_list_input[[set_lab]] <- df %>%
          filter(avg_log2FC <= down_logfc_cutoff, p_val_adj <= padj_cutoff) %>%
          pull(gene) %>% as.character()
      } else {
        up_list_input[[set_lab]] <- character(0)
        down_list_input[[set_lab]] <- character(0)
      }
    }
    
    up_data <- fromList(up_list_input)
    down_data <- fromList(down_list_input)
    
    # Save membership tables
    write.table(
      up_data,
      file.path(sex_outdir, paste0(ct, "_", sex, "_upregulated.txt")),
      sep = "\t", quote = FALSE, col.names = NA
    )
    write.table(
      down_data,
      file.path(sex_outdir, paste0(ct, "_", sex, "_downregulated.txt")),
      sep = "\t", quote = FALSE, col.names = NA
    )
    
    # UpSet plots
    p_up   <- generate_upset_plot(up_data, ct, sex, "upregulated", set_names, "red")
    p_down <- generate_upset_plot(down_data, ct, sex, "downregulated", set_names, "blue")
    
    # DEG count plot
    p_deg_counts <- generate_deg_count_plot(deg_counts_all, ct, new_plot_labels_order)
    
    # Combine and save
    combined_plot_fname <- file.path(sex_outdir, paste0(ct, "_", sex, "_Combined_UpSet_Horizontal_Counts.pdf"))
    message("Saving combined plot: ", combined_plot_fname)
    
    p_combined <- plot_grid(
      p_deg_counts,
      p_down,
      p_up,
      ncol = 3,
      rel_widths = c(1, 1.25, 1.25),
      labels = c("A", "B", "C"),
      label_size = 10,
      align = "v"
    )
    
    pdf(combined_plot_fname, width = 7, height = 2.75)
    print(p_combined)
    dev.off()
  }
  
  message("\nAll cell types processed for sex: ", sex)
  message("Saved outputs in: ", sex_outdir)
  
  # ------------------------------------------------------------
  # Excel summary (one workbook per sex)
  # ------------------------------------------------------------
  excel_output_path <- file.path(sex_outdir, paste0("All_CellType_DEGs_Summary_Lists_LBDtypes_vs_AD_AT_", sex, ".xlsx"))
  wb <- createWorkbook()
  
  message("\nGenerating Excel Summary file: ", excel_output_path)
  
  for (ct in cell_types) {
    message("  Preparing Excel sheets for cell type: ", ct)
    
    # Map set label -> comparison string (with sex)
    set_to_comp <- setNames(
      vapply(set_names, function(lbl) paste0(lbl, "_", sex, "_vs_AD_AT_", sex), character(1)),
      set_names
    )
    
    for (reg_status in c("upregulated", "downregulated")) {
      logfc_filter_fn <- if (reg_status == "upregulated") {
        function(df) df %>% filter(avg_log2FC >= up_logfc_cutoff, p_val_adj <= padj_cutoff)
      } else {
        function(df) df %>% filter(avg_log2FC <= down_logfc_cutoff, p_val_adj <= padj_cutoff)
      }
      
      comp_gene_lists <- list()
      for (set_lab in set_names) {
        comp <- set_to_comp[[set_lab]]
        df <- all_degs[[ct]][[comp]]
        
        if (!is.null(df) && nrow(df) > 0) {
          comp_gene_lists[[set_lab]] <- logfc_filter_fn(df) %>%
            pull(gene) %>%
            as.character()
        } else {
          comp_gene_lists[[set_lab]] <- character(0)
        }
      }
      
      max_len <- max(lengths(comp_gene_lists), 0)
      if (max_len > 0) {
        padded_lists <- lapply(comp_gene_lists, function(x) {
          length(x) <- max_len
          x
        })
        
        combined_df <- as.data.frame(padded_lists, stringsAsFactors = FALSE)
        
        sheet_name <- paste0(ct, "_", ifelse(reg_status == "upregulated", "UP", "DOWN"))
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, combined_df, startRow = 1, startCol = 1)
        message("    - Added sheet: ", sheet_name, " (Max rows: ", max_len, ")")
      } else {
        message("    - Skipping ", reg_status, " for ", ct, " (No DEGs).")
      }
    }
  }
  
  tryCatch({
    saveWorkbook(wb, excel_output_path, overwrite = TRUE)
    message("\nSuccessfully saved Excel file: ", excel_output_path)
  }, error = function(e) {
    stop("Error saving Excel file: ", e$message)
  })
}

message("\nDONE. All sex-stratified LBD-type vs AD_AT UpSet + count plots and Excel summaries are complete.")
message("Base output directory: ", output_dir)
