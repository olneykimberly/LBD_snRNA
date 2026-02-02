source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(forcats)
library(tidyr)
library(stringr)
library(patchwork)

addSmallLegend <- function(myPlot, pointSize = 1, textSize = 7, spaceLegend = .5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = -.25, unit = "cm"))
}
# --- New User Input for Comparison Control ---
# Define the comparisons to include and their order
comparisons_to_plot <- list(
  c("AD_AT", "CONTROL"),
  c("LBD_S", "CONTROL"),
  c("LBD_AS", "CONTROL"),
  c("LBD_ATS", "CONTROL")
)

# Convert to the expected "condition_vs_control" string format and order
plot_order <- comparisons_to_plot %>% 
  map_chr(~ paste0(.x[1], "_vs_", .x[2]))

# Create the corresponding human-readable labels
plot_labels <- comparisons_to_plot %>%
  map_chr(~ paste0(.x[1], " vs ", .x[2]))

# Create a mapping data frame
comparison_map <- tibble(
  comparison = plot_order,
  plot_label = plot_labels
)
# ---------------------------------------------


# user inputs / thresholds
main_output_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/DEGs_pseudobulk_limma_voom"
cell_types <- c("neuron", "interneuron", "oligodendrocyte", "opc", "astrocyte", "microglia", "endothelial", "fibroblast", "mural")

up_logfc_cutoff<- 0.25
down_logfc_cutoff <- -0.25
padj_cutoff <- 0.05

# find all DEG files in the directory tree that follow the pattern DEG_<cell>_..._vs_....txt
all_files <- list.files(main_output_dir, pattern = "^DEG_.*_vs_.*\\.(txt|csv|tsv)$", recursive = TRUE, full.names = TRUE)
if(length(all_files) == 0) stop("No DEG files found under main_output_dir with pattern 'DEG_*_vs_*.txt'")

# parse filename to get cell_type and comparison
# expected basename like: DEG_astrocyte_LBD_ATS_vs_LBD_AS.txt
file_df <- tibble(path = all_files) %>%
  mutate(bn = basename(path)) %>%
  mutate(noext = tools::file_path_sans_ext(bn)) %>%
  # Updated to remove the longer prefix
  mutate(after_prefix = str_remove(noext, "^DEG_pseudobulk_voom_")) %>%
  # Now "after_prefix" starts with the cell type (e.g., "astrocyte")
  mutate(cell_type = word(after_prefix, 1, sep = fixed("_")),
         comparison = str_remove(after_prefix, paste0("^", cell_type, "_"))) %>%
  select(path, cell_type, comparison)

# Keep only files where cell_type is one of your known cell types (safe-guard)
file_df <- file_df %>% filter(cell_type %in% cell_types)

file_df <- file_df %>% 
  filter(comparison %in% plot_order)

# ---------- Replace deg_long with this block (identical to previous version) ----------
deg_long <- file_df %>%
  mutate(data = map(path, safe_read_deg)) %>%
  # drop failed reads
  filter(map_lgl(data, ~ !is.null(.x))) %>%
  select(-path) %>%
  tidyr::unnest(data) %>%# expands the read tables into rows
  # standardize column names
  dplyr::rename_with(~ tolower(.x)) %>%
  # create standardized columns (use variants if present)
  dplyr::mutate(
    adj.p.val = dplyr::coalesce(
      .data$adj.p.val
    ),
    logfc = dplyr::coalesce(
      .data$logfc
    )
  ) %>%
  # classify up / down
  dplyr::mutate(
    direction = dplyr::case_when(
      !is.na(adj.p.val) & !is.na(logfc) & adj.p.val <= padj_cutoff & logfc >= up_logfc_cutoff~ "Up",
      !is.na(adj.p.val) & !is.na(logfc) & adj.p.val <= padj_cutoff & logfc <= down_logfc_cutoff ~ "Down",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(direction)) %>%
  # count per combination (safe and clear)
  dplyr::count(comparison, cell_type, direction, name = "n")
# -------------------------------------------------------


# ensure every cell_type x comparison x direction exists (fill zeros)
all_combinations <- expand_grid(
  comparison = plot_order, # Use the defined order
  cell_type = cell_types,
  direction = c("Down", "Up")
)

deg_counts <- all_combinations %>%
  left_join(deg_long, by = c("comparison","cell_type","direction")) %>%
  mutate(n = replace_na(n, 0)) %>%
  # convert to mirrored values: Down negative, Up positive
  mutate(value = if_else(direction == "Down", -n, n)) %>%
  # order cell types so they plot top->bottom as in your example
  mutate(cell_type = factor(cell_type, levels = rev(cell_types))) %>%
  # --- NEW: Order comparisons and add plot labels ---
  left_join(comparison_map, by = "comparison") %>%
  mutate(plot_label = factor(plot_label, levels = plot_labels))
# --------------------------------------------------


##Faceted Plotting

# Calculate the maximum absolute count to set uniform limits across all comparisons
max_count <- max(deg_counts$n, na.rm = TRUE)
# Set the x-axis limits (e.g., max_count + 10% for padding)
max_abs_x <- ceiling(max_count * 1.1 / 10) * 10 # Round up to the nearest 10 for cleaner limits

# Create a single, faceted plot
##Faceted Plotting

# Calculate the maximum absolute count to set uniform limits across all comparisons
max_count <- max(deg_counts$n, na.rm = TRUE)
# Set the x-axis limits (e.g., max_count + 10% for padding)
max_abs_x <- ceiling(max_count * 1.5 / 10) * 10 # Round up to the nearest 10 for cleaner limits

final_plot <- ggplot(deg_counts, aes(x = value, y = cell_type, fill = direction)) +
  geom_col(width = 0.7) +
  geom_text(
    data = deg_counts %>% filter(n > 0),
    aes(
      y = cell_type, 
      label = as.character(n), 
      x = value
    ),
    hjust = ifelse((deg_counts %>% filter(n > 0))$value < 0, 1.1, -0.1),
    size = 2.5,
    color = "black"
  ) +
  scale_fill_manual(values = c("Down" = "blue", "Up" = "red")) +
  scale_x_continuous(
    labels = function(x) abs(x),
    limits = c(-max_abs_x, max_abs_x) 
  ) +
  facet_grid(. ~ plot_label) + 
  labs(title = "",
       x = "Number of DEGs\nadjusted p-value < 0.05 & |log2FC| > 0.25", 
       y = NULL, 
       fill = "Direction") +
  theme_bw() +
  theme(
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dashed", color = "lightgray"), 
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "gray90", color = NA),
        strip.text = element_text(face = "bold", size = 7.5),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"), # Minimal margin around the legend box
        plot.margin = margin(-.25, 0.025, 0.05, 0.05, "cm")
      )
    
final_plot <- addSmallLegend(final_plot)
print(final_plot)
path <- paste0("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/manuscript_figures/DEGs_total_pseudobulk")
print(final_plot)
saveToPDF(paste0(path, ".pdf"), width = 7, height = 3)
