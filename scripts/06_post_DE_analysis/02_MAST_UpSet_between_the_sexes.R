#!/usr/bin/env Rscript
# ==============================================================================
# Sex-comparison DEG counts + ONE UpSet showing Up/Down in BOTH sexes
#
# For each CELL TYPE and each DISEASE-vs-CONTROL comparison:
#   A) DEG count bar plot showing FEMALE (top) vs MALE (bottom) (Up vs Down counts)
#   B) ONE UpSet plot with 4 sets:
#        Female Up, Female Down, Male Up, Male Down
#      Key intersections:
#        - same-direction:   Female Up ∩ Male Up, Female Down ∩ Male Down
#        - opposite:         Female Up ∩ Male Down, Female Down ∩ Male Up
#
# Input DEG files expected at:
#   {main_output_dir}/{celltype}/DEG_{celltype}_{group1}_vs_{group2}.tsv
# Example:
#   DEG_neuron_LBD_ATS_female_vs_CONTROL_female.tsv
#   DEG_neuron_LBD_ATS_male_vs_CONTROL_male.tsv
#
# Outputs written under:
#   {output_dir}/{celltype}/
#     - {celltype}_{comparison}_SexCompare_Combined_CountsPlusUpDownUpSet.pdf
#     - {celltype}_{comparison}_SexCompare_UpDown4set.txt
#   Plus one Excel workbook:
#     - All_CellType_SexCompare_UpDown4set_GeneLists.xlsx
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
output_dir      <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/UpSet/sex_compare_within_disease_vs_control"

# DEG Filtering Parameters
up_logfc_cutoff   <- 0.25
down_logfc_cutoff <- -0.25
padj_cutoff       <- 0.05

# Cell types to process
cell_types <- c("neuron", "interneuron", "oligodendrocyte",
                "opc", "astrocyte", "microglia", "mural", "endothelial")

# Disease-vs-control comparisons to run (sex will be appended)
comparisons_base <- list(
  c("AD_AT",   "CONTROL"),
  c("LBD_S",   "CONTROL"),
  c("LBD_AS",  "CONTROL"),
  c("LBD_ATS", "CONTROL")
)

# ------------------------------
# Helpers
# ------------------------------
ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

make_sex_comp <- function(g1_base, g2_base, sex) {
  g1 <- paste0(g1_base, "_", sex)
  g2 <- paste0(g2_base, "_", sex)
  list(g1 = g1, g2 = g2, comp = paste0(g1, "_vs_", g2))
}

pretty_comp_label <- function(comp_label) {
  # Example: "AD_AT_vs_CONTROL" -> "AD AT vs CONTROL"
  gsub("_", " ", comp_label)
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
  
  df <- readr::read_tsv(fname, show_col_types = FALSE)
  
  required <- c("gene", "avg_log2FC", "p_val_adj")
  if (!all(required %in% colnames(df))) {
    stop("DEG file ", fname, " missing required columns: ", paste(required, collapse = ", "))
  }
  
  df %>% dplyr::select(dplyr::all_of(required))
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
# Panel A: DEG count bar plot (Female top, Male bottom)
# ------------------------------
generate_deg_count_plot_sexcompare <- function(deg_counts_df, current_celltype, comparison_label_pretty) {
  df <- deg_counts_df %>%
    dplyr::filter(cell_type == current_celltype) %>%
    dplyr::mutate(sex = factor(sex, levels = c("Female", "Male")))
  
  if (nrow(df) == 0 || sum(df$n) == 0) {
    return(
      ggplot() +
        geom_text(
          aes(x = 0.5, y = 0.5, label = paste(current_celltype, "DEGs: No DEGs found")),
          size = 3
        ) +
        theme_void()
    )
  }
  
  max_count <- max(df$n, na.rm = TRUE)
  max_abs_x <- ceiling(max_count * 1.5 / 10) * 10
  
  ggplot(df, aes(x = value, y = sex, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(
      data = df %>% dplyr::filter(n > 0),
      aes(y = sex, label = as.character(n), x = value),
      hjust = ifelse((df %>% dplyr::filter(n > 0))$value < 0, 1.1, -0.1),
      size = 2.8,
      color = "black"
    ) +
    scale_fill_manual(values = c("Down" = "blue", "Up" = "red")) +
    scale_x_continuous(
      labels = function(x) abs(x),
      limits = c(-max_abs_x, max_abs_x),
      breaks = c(-max_abs_x, 0, max_abs_x)
    ) +
    labs(
      title = paste0(current_celltype, " DEG Counts\n", comparison_label_pretty),
      x = "Number of DEGs",
      y = NULL,
      fill = "Direction"
    ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 8),
      plot.title = element_text(size = 8, hjust = 0.5),
      panel.grid.major.x = element_line(linetype = "dashed", color = "lightgray"),
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.key.size = unit(0.5, "lines"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      plot.margin = margin(0.2, 0.15, 0.1, 0.2, "cm")
    )
}

# ------------------------------
# Panel B: ONE UpSet plot (4 sets) with cleaned titles/labels
#   - Drops "_" in titles/labels
#   - Does NOT include the 4-way intersection unless it has genes
# ------------------------------
generate_upset_plot_4set <- function(data, celltype, comparison_label_pretty) {
  set_names <- c("Female Up", "Female Down", "Male Up", "Male Down")
  
  if (nrow(data) == 0) {
    return(
      ggplot() +
        geom_text(
          aes(x = 0.5, y = 0.5,
              label = paste(celltype, comparison_label_pretty, "\nNo DEGs found")),
          size = 4
        ) +
        theme_void()
    )
  }
  
  # Ensure columns exist in correct order (matching set_names)
  data <- data[, set_names, drop = FALSE]
  
  # Color matrix dots by direction (Up=red, Down=blue)
  upset_queries <- list(
    upset_query(set = "Female Up",   fill = "red"),
    upset_query(set = "Male Up",     fill = "red"),
    upset_query(set = "Female Down", fill = "blue"),
    upset_query(set = "Male Down",   fill = "blue")
  )
  
  # Add the 4-way intersection ONLY if present (rowname == gene, but intersection size is computed from membership)
  # We'll detect presence by checking if any gene has all 4 columns == 1
  has_all4 <- any(rowSums(data[, set_names, drop = FALSE]) == 4)
  
  # Curated intersections to highlight same vs opposite directions
  sets_list <- list(
    "Female Up", "Female Down", "Male Up", "Male Down",
    c("Female Up", "Male Up"),
    c("Female Down", "Male Down"),
    c("Female Up", "Male Down"),
    c("Female Down", "Male Up")
  )
  if (has_all4) {
    sets_list <- c(sets_list, list(c("Female Up", "Female Down", "Male Up", "Male Down")))
  }
  
  plot_title <- paste0(celltype, "\n", comparison_label_pretty)
  
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
          size = 2.25,
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
# Gene list helpers
# ------------------------------
get_up_genes <- function(df) {
  if (nrow(df) == 0) return(character(0))
  df %>%
    dplyr::filter(avg_log2FC >= up_logfc_cutoff, p_val_adj <= padj_cutoff) %>%
    dplyr::pull(gene) %>%
    as.character()
}

get_down_genes <- function(df) {
  if (nrow(df) == 0) return(character(0))
  df %>%
    dplyr::filter(avg_log2FC <= down_logfc_cutoff, p_val_adj <= padj_cutoff) %>%
    dplyr::pull(gene) %>%
    as.character()
}

# Collect for Excel
excel_records <- list()

ensure_dir(output_dir)

# ==============================================================================
# MAIN
# ==============================================================================
for (pair_base in comparisons_base) {
  
  g1_base <- pair_base[1]
  g2_base <- pair_base[2]
  comparison_label <- paste0(g1_base, "_vs_", g2_base)      # e.g., AD_AT_vs_CONTROL
  comparison_label_pretty <- pretty_comp_label(comparison_label)
  
  message("\n============================================================")
  message("COMPARISON: ", comparison_label, " (Female vs Male; combined Up/Down UpSet)")
  message("============================================================\n")
  
  for (ct in cell_types) {
    message("Processing: ", ct, " | ", comparison_label)
    
    ct_outdir <- file.path(output_dir, ct)
    ensure_dir(ct_outdir)
    
    fem <- make_sex_comp(g1_base, g2_base, "female")
    mal <- make_sex_comp(g1_base, g2_base, "male")
    
    df_f <- read_deg_file(ct, fem$g1, fem$g2, main_output_dir); if (is.null(df_f)) df_f <- tibble()
    df_m <- read_deg_file(ct, mal$g1, mal$g2, main_output_dir); if (is.null(df_m)) df_m <- tibble()
    
    # --------------------------
    # A) DEG counts (Up/Down) by sex
    # --------------------------
    count_one <- function(df, sex_label) {
      if (nrow(df) == 0) {
        return(tibble(sex = sex_label, direction = c("Down", "Up"), n = c(0L, 0L)))
      }
      
      df %>%
        dplyr::mutate(direction = dplyr::case_when(
          !is.na(p_val_adj) & !is.na(avg_log2FC) & p_val_adj <= padj_cutoff & avg_log2FC >= up_logfc_cutoff   ~ "Up",
          !is.na(p_val_adj) & !is.na(avg_log2FC) & p_val_adj <= padj_cutoff & avg_log2FC <= down_logfc_cutoff ~ "Down",
          TRUE ~ NA_character_
        )) %>%
        dplyr::filter(!is.na(direction)) %>%
        dplyr::count(direction, name = "n") %>%
        dplyr::mutate(sex = sex_label) %>%
        dplyr::right_join(tibble(direction = c("Down", "Up")), by = "direction") %>%
        dplyr::mutate(n = tidyr::replace_na(n, 0L)) %>%
        dplyr::select(sex, direction, n)
    }
    
    deg_counts <- dplyr::bind_rows(
      count_one(df_f, "Female"),
      count_one(df_m, "Male")
    ) %>%
      dplyr::mutate(
        value = dplyr::if_else(direction == "Down", -n, n),
        cell_type = ct
      )
    
    p_deg_counts <- generate_deg_count_plot_sexcompare(deg_counts, ct, comparison_label_pretty)
    
    # --------------------------
    # B) ONE UpSet with 4 sets (Female Up/Down, Male Up/Down)
    # --------------------------
    lists_4set <- list(
      "Female Up"   = get_up_genes(df_f),
      "Female Down" = get_down_genes(df_f),
      "Male Up"     = get_up_genes(df_m),
      "Male Down"   = get_down_genes(df_m)
    )
    
    upset_data <- fromList(lists_4set)
    
    prefix <- paste0(ct, "_", comparison_label, "_SexCompare")
    write.table(
      upset_data,
      file.path(ct_outdir, paste0(prefix, "_UpDown4set.txt")),
      sep = "\t", quote = FALSE, col.names = NA
    )
    
    p_updown <- generate_upset_plot_4set(upset_data, ct, comparison_label_pretty)
    
    # --------------------------
    # Combine and save (A + B)
    # --------------------------
    combined_plot_fname <- file.path(ct_outdir, paste0(prefix, "_Combined_CountsPlusUpDownUpSet.pdf"))
    message("  Saving combined plot: ", combined_plot_fname)
    
    p_combined <- plot_grid(
      p_deg_counts,
      p_updown,
      ncol = 2,
      rel_widths = c(1, 2.1),
      labels = c("A", "B"),
      label_size = 10,
      align = "v"
    )
    
    pdf(combined_plot_fname, width = 7, height = 2.5)
    print(p_combined)
    dev.off()
    
    # --------------------------
    # Collect for Excel
    # --------------------------
    excel_records[[length(excel_records) + 1]] <- list(
      cell_type = ct,
      comparison = comparison_label,
      Female_Up = lists_4set[["Female Up"]],
      Female_Down = lists_4set[["Female Down"]],
      Male_Up = lists_4set[["Male Up"]],
      Male_Down = lists_4set[["Male Down"]]
    )
  }
}

message("\nAll plots generated.")
message("Outputs saved under: ", output_dir)

# ==============================================================================
# Excel summary workbook (one sheet per ct+comparison)
# ==============================================================================
excel_output_path <- file.path(output_dir, "All_CellType_SexCompare_UpDown4set_GeneLists.xlsx")
wb <- createWorkbook()
message("\nGenerating Excel Summary file: ", excel_output_path)

for (rec in excel_records) {
  ct <- rec$cell_type
  comp <- rec$comparison
  
  sheet_name <- paste0(ct, "__", comp) # Excel limit = 31 chars
  sheet_name <- substr(sheet_name, 1, 31)
  
  # Ensure unique sheet names if truncation causes collisions
  base_name <- sheet_name
  k <- 1
  while (sheet_name %in% names(wb)) {
    k <- k + 1
    suffix <- paste0("_", k)
    sheet_name <- substr(base_name, 1, 31 - nchar(suffix))
    sheet_name <- paste0(sheet_name, suffix)
  }
  
  addWorksheet(wb, sheet_name)
  
  cols <- list(
    Female_Up   = rec$Female_Up,
    Female_Down = rec$Female_Down,
    Male_Up     = rec$Male_Up,
    Male_Down   = rec$Male_Down
  )
  
  max_len <- max(lengths(cols), 0)
  if (max_len == 0) {
    writeData(wb, sheet_name, data.frame(Note = "No DEGs found in either sex for this comparison/cell type."))
    next
  }
  
  padded <- lapply(cols, function(x) {
    length(x) <- max_len
    x
  })
  out_df <- as.data.frame(padded, stringsAsFactors = FALSE)
  
  writeData(wb, sheet_name, out_df, startRow = 1, startCol = 1)
}

tryCatch({
  saveWorkbook(wb, excel_output_path, overwrite = TRUE)
  message("\nSuccessfully saved Excel file: ", excel_output_path)
}, error = function(e) {
  stop("Error saving Excel file: ", e$message)
})

message("\nDONE.")
