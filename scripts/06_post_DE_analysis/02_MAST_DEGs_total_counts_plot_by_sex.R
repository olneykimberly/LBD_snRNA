source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(grid)   # unit(), margin()
})

addSmallLegend <- function(myPlot, pointSize = 1, textSize = 7, spaceLegend = .5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(
      legend.title = element_text(size = textSize),
      legend.text  = element_text(size = textSize),
      legend.key.size = unit(spaceLegend, "lines"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = -.25, unit = "cm")
    )
}

# ------------------------------------------------------------------------------
# User inputs / thresholds
# ------------------------------------------------------------------------------
main_output_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/DEGs_RNA_pct0.25_by_sex"
cell_types <- c("neuron", "interneuron", "oligodendrocyte", "opc", "astrocyte",
                "microglia", "endothelial", "fibroblast", "mural")

up_logfc_cutoff <- 0.25
down_logfc_cutoff <- -0.25
padj_cutoff <- 0.05

# ------------------------------------------------------------------------------
# Define comparisons for each sex
# ------------------------------------------------------------------------------
comparisons_female <- list(
  c("AD_AT_female", "CONTROL_female"),
  c("LBD_S_female", "CONTROL_female"),
  c("LBD_AS_female", "CONTROL_female"),
  c("LBD_ATS_female", "CONTROL_female")
)

comparisons_male <- list(
  c("AD_AT_male", "CONTROL_male"),
  c("LBD_S_male", "CONTROL_male"),
  c("LBD_AS_male", "CONTROL_male"),
  c("LBD_ATS_male", "CONTROL_male")
)

# ------------------------------------------------------------------------------
# Helpers for pretty facet labels
# ------------------------------------------------------------------------------
pretty <- function(x) {
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "\\bCONTROL\\b", "Control")
  x
}

make_comparison_map <- function(comparisons_to_plot) {
  plot_order <- comparisons_to_plot %>% purrr::map_chr(~ paste0(.x[1], "_vs_", .x[2]))
  
  plot_labels <- comparisons_to_plot %>%
    purrr::map_chr(function(x) {
      sex <- if (stringr::str_detect(x[1], "_female$")) "Female" else if (stringr::str_detect(x[1], "_male$")) "Male" else ""
      disease <- stringr::str_remove(x[1], "_(female|male)$")
      control <- stringr::str_remove(x[2], "_(female|male)$")
      paste0(sex, "\n", pretty(disease), " vs ", pretty(control))
    })
  
  list(
    plot_order = plot_order,
    plot_labels = plot_labels,
    comparison_map = tibble::tibble(comparison = plot_order, plot_label = plot_labels)
  )
}

# ------------------------------------------------------------------------------
# Safe DEG reader (expects gene + avg_log2FC + p_val_adj)
# ------------------------------------------------------------------------------
safe_read_deg <- function(path) {
  if (!file.exists(path)) return(NULL)
  
  ext <- tolower(tools::file_ext(path))
  d <- tryCatch({
    if (ext == "csv") {
      readr::read_csv(path, show_col_types = FALSE)
    } else if (ext %in% c("tsv", "txt")) {
      readr::read_tsv(path, show_col_types = FALSE)
    } else {
      return(NULL)
    }
  }, error = function(e) NULL)
  
  if (is.null(d)) return(NULL)
  
  names(d) <- tolower(names(d))
  
  # tolerate common variants
  if (!"gene" %in% names(d)) {
    if ("genes" %in% names(d)) d <- dplyr::rename(d, gene = genes)
    if ("feature" %in% names(d)) d <- dplyr::rename(d, gene = feature)
  }
  if (!"avg_log2fc" %in% names(d)) {
    if ("avg_logfc" %in% names(d)) d <- dplyr::rename(d, avg_log2fc = avg_logfc)
    if ("log2fc" %in% names(d)) d <- dplyr::rename(d, avg_log2fc = log2fc)
  }
  if (!"p_val_adj" %in% names(d)) {
    if ("padj" %in% names(d)) d <- dplyr::rename(d, p_val_adj = padj)
  }
  
  req <- c("gene", "avg_log2fc", "p_val_adj")
  if (!all(req %in% names(d))) return(NULL)
  
  d %>%
    dplyr::transmute(
      gene = as.character(.data$gene),
      avg_log2fc = suppressWarnings(as.numeric(.data$avg_log2fc)),
      p_val_adj  = suppressWarnings(as.numeric(.data$p_val_adj))
    ) %>%
    dplyr::filter(!is.na(gene) & gene != "", !is.na(avg_log2fc), !is.na(p_val_adj))
}

# ------------------------------------------------------------------------------
# Core function: build mirrored DEG count plot for a given sex + comparisons list
# ------------------------------------------------------------------------------
make_deg_count_plot <- function(comparisons_to_plot, out_pdf) {
  
  map_obj <- make_comparison_map(comparisons_to_plot)
  plot_order <- map_obj$plot_order
  plot_labels <- map_obj$plot_labels
  comparison_map <- map_obj$comparison_map
  
  # find all DEG files in the directory tree
  all_files <- list.files(
    main_output_dir,
    pattern = "^DEG_.*_vs_.*\\.(txt|csv|tsv)$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(all_files) == 0) {
    stop("No DEG files found under main_output_dir with pattern 'DEG_*_vs_*.(txt|csv|tsv)'")
  }
  
  # parse filenames: DEG_<cell_type>_<comparison>.{txt|csv|tsv}
  file_df <- tibble::tibble(path = all_files) %>%
    dplyr::mutate(
      bn = basename(path),
      noext = tools::file_path_sans_ext(bn),
      after_deg = stringr::str_remove(noext, "^DEG_"),
      cell_type = stringr::word(after_deg, 1, sep = stringr::fixed("_")),
      rest = stringr::str_remove(after_deg, paste0("^", cell_type, "_")),
      comparison = rest
    ) %>%
    dplyr::select(path, cell_type, comparison) %>%
    dplyr::filter(cell_type %in% cell_types) %>%
    dplyr::filter(comparison %in% plot_order)
  
  if (nrow(file_df) == 0) {
    stop("After filtering to requested comparisons, no DEG files matched. Check filenames under: ", main_output_dir)
  }
  
  # read and classify up/down; count
  deg_long <- file_df %>%
    dplyr::mutate(data = purrr::map(path, safe_read_deg)) %>%
    dplyr::filter(purrr::map_lgl(data, ~ !is.null(.x))) %>%
    dplyr::select(-path) %>%
    tidyr::unnest(data) %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        p_val_adj <= padj_cutoff & avg_log2fc >= up_logfc_cutoff ~ "Up",
        p_val_adj <= padj_cutoff & avg_log2fc <= down_logfc_cutoff ~ "Down",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(direction)) %>%
    # IMPORTANT: namespace to avoid masked count()
    dplyr::count(comparison, cell_type, direction, name = "n")
  
  # ensure every cell_type x comparison x direction exists (fill zeros)
  all_combinations <- tidyr::expand_grid(
    comparison = plot_order,
    cell_type = cell_types,
    direction = c("Down", "Up")
  )
  
  deg_counts <- all_combinations %>%
    dplyr::left_join(deg_long, by = c("comparison", "cell_type", "direction")) %>%
    dplyr::mutate(n = tidyr::replace_na(n, 0L)) %>%
    dplyr::mutate(value = dplyr::if_else(direction == "Down", -n, n)) %>%
    dplyr::mutate(cell_type = factor(cell_type, levels = rev(cell_types))) %>%
    dplyr::left_join(comparison_map, by = "comparison") %>%
    dplyr::mutate(plot_label = factor(plot_label, levels = plot_labels))
  
  # axis limits
  max_count <- max(deg_counts$n, na.rm = TRUE)
  max_abs_x <- ceiling(max_count * 1.5 / 10) * 10
  if (!is.finite(max_abs_x) || max_abs_x == 0) max_abs_x <- 10
  
  # build plot
  lbl_df <- deg_counts %>% dplyr::filter(n > 0)
  
  final_plot <- ggplot(deg_counts, aes(x = value, y = cell_type, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(
      data = lbl_df,
      aes(label = as.character(n), x = value),
      hjust = ifelse(lbl_df$value < 0, 1.1, -0.1),
      size = 2.5,
      color = "black"
    ) +
    scale_fill_manual(values = c("Down" = "blue", "Up" = "red")) +
    scale_x_continuous(
      labels = function(x) abs(x),
      limits = c(-max_abs_x, max_abs_x)
    ) +
    facet_grid(. ~ plot_label) +
    labs(
      title = "",
      x = "Number of DEGs\nadjusted p-value < 0.05 & |log2FC| > 0.25",
      y = NULL,
      fill = "Direction"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      plot.title   = element_text(size = 8),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(linetype = "dashed", color = "lightgray"),
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "gray90", color = NA),
      strip.text = element_text(face = "plain", size = 7.5, lineheight = 0.95),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      plot.margin = margin(-.25, 0.025, 0.05, 0.05, "cm")
    )
  
  final_plot <- addSmallLegend(final_plot)
  
  print(final_plot)
  saveToPDF(out_pdf, width = 7, height = 3)
  
  invisible(final_plot)
}

# ------------------------------------------------------------------------------
# Run: separate plots for Female and Male
# ------------------------------------------------------------------------------
out_base <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/manuscript_figures"
female_pdf <- file.path(out_base, "DEGs_total_MAST_female.pdf")
male_pdf   <- file.path(out_base, "DEGs_total_MAST_male.pdf")

female_plot <- make_deg_count_plot(comparisons_female, female_pdf)
male_plot   <- make_deg_count_plot(comparisons_male, male_pdf)
