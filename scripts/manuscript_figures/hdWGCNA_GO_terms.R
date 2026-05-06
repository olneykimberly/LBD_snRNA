#!/usr/bin/env Rscript

# ============================================================
# GO tile plot per module
# ONLY 1_Summary text = white, others = black
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

base_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/metascape_hdWGCNA_modules"

outdir <- file.path(base_dir, "GO_tileplots_by_module")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

top_n_summary <- 6

plot_width  <- 1.75
plot_height <- 1.35
base_size <- 6

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

clean_go_description <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_replace_all(regex("^HALLMARK\\s+", ignore_case = TRUE), "") %>%
    str_replace_all(regex("\\bHALLMARK\\s+", ignore_case = TRUE), "") %>%
    str_to_lower() %>%
    str_replace_all("positive regulation of", "pos. reg. of") %>%
    str_replace_all("negative regulation of", "neg. reg. of") %>%
    str_replace_all("regulation of", "reg. of") %>%
    str_replace_all("response to", "resp. to") %>%
    str_squish()
}

read_metascape_module <- function(module_dir) {

  xlsx_file <- file.path(module_dir, "metascape_result.xlsx")

  if (!file.exists(xlsx_file)) {
    message("Skipping: ", module_dir)
    return(NULL)
  }

  module_folder <- basename(module_dir)

  sheets <- excel_sheets(xlsx_file)

  if (length(sheets) < 2) {
    message("Skipping: fewer than 2 sheets in ", module_folder)
    return(NULL)
  }

  df <- read_excel(xlsx_file, sheet = 2)

  required_cols <- c("GroupID", "Description", "LogP")
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    message("Skipping ", module_folder, " missing: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }

  summary_ids <- paste0(seq_len(top_n_summary), "_Summary")

  df %>%
    mutate(
      GroupID = as.character(GroupID),
      module_folder = module_folder,
      Description_clean = clean_go_description(Description),
      LogP = as.numeric(LogP),
      neg_log10_p = abs(LogP),
      summary_rank = match(GroupID, summary_ids),

      # 🔥 ONLY first summary white
      text_color = ifelse(GroupID == "1_Summary", "white", "black")
    ) %>%
    filter(
      GroupID %in% summary_ids,
      !is.na(summary_rank),
      !is.na(Description_clean),
      Description_clean != "",
      !is.na(neg_log10_p)
    ) %>%
    arrange(summary_rank)
}

make_go_plot <- function(df_one_module) {

  df_plot <- df_one_module %>%
    arrange(desc(summary_rank)) %>%
    mutate(
      Description_clean = factor(
        Description_clean,
        levels = unique(Description_clean)
      ),
      x_position = 1
    )

  ggplot(df_plot, aes(x = x_position, y = Description_clean)) +

    geom_tile(
      aes(fill = neg_log10_p),
      width = 1,
      height = 0.95,
      color = "white",
      linewidth = 0.15
    ) +

    geom_text(
      aes(label = Description_clean, color = text_color),
      size = 6 / .pt
    ) +

    scale_color_identity() +

    theme_bw(base_size = base_size) +
    labs(
      x = NULL,
      y = NULL,
      fill = expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")
    ) +

    scale_x_continuous(
      limits = c(0.5, 1.5),
      breaks = NULL,
      expand = c(0, 0)
    ) +

    scale_fill_gradient(
      low = "gray95",
      high = "navy",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(0.7, "in"),
        barheight = unit(0.06, "in")
      )
    ) +

    theme(
      text = element_text(size = base_size),
      plot.title = element_blank(),

      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),

      panel.border = element_blank(),
      panel.grid = element_blank(),

      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(0.2, "lines"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-5, 0, 0, 0),

      plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm")
    )
}

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

module_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
module_dirs <- module_dirs[
  str_detect(basename(module_dirs), "^hdWGCNA_.+_M[0-9]+$")
]

for (module_dir in module_dirs) {

  df_one <- read_metascape_module(module_dir)
  if (is.null(df_one) || nrow(df_one) == 0) next

  module_folder <- unique(df_one$module_folder)

  p <- make_go_plot(df_one)

  ggsave(
    file.path(outdir, paste0(module_folder, "_GO_tileplot.pdf")),
    plot = p,
    width = plot_width,
    height = plot_height,
    units = "in",
    useDingbats = FALSE
  )

  ggsave(
    file.path(outdir, paste0(module_folder, "_GO_tileplot.png")),
    plot = p,
    width = plot_width,
    height = plot_height,
    units = "in",
    dpi = 300
  )

  message("Saved: ", module_folder)
}

message("Done.")