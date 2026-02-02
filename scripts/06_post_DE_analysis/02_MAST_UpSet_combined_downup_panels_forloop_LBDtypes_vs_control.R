# ------------------------------
# Libraries
# ------------------------------
library(tidyverse)
library(ComplexUpset)
library(readr)
library(openxlsx)
library(cowplot) 
library(ggplot2)
library(ggpubr)

# ------------------------------
# Directories and Parameters 
# ------------------------------
# NOTE: These directories must exist and contain the expected files for the script to run fully.
main_output_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/DEGs_RNA_pct0.25"
output_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/UpSet/LBDtypes_vs_control" 

# DEG Filtering Parameters
up_logfc_cutoff <- 0.25
down_logfc_cutoff <- -0.25
padj_cutoff <- 0.05

# List of comparison pairs
comparisons_to_read <- list(
  c("LBD_S",   "CONTROL"),
  c("LBD_AS",  "CONTROL"),
  c("LBD_ATS", "CONTROL")
)
comparison_names <- sapply(comparisons_to_read, function(p) paste0(p[1], "_vs_", p[2]))

# List of cell types to process
cell_types <- c("neuron", "interneuron", "oligodendrocyte",
                "opc", "astrocyte", "microglia", "mural", "endothelial")

# ------------------------------
# NEW: Setup for DEG Count Plot Order
# ------------------------------

# Define the user-requested order for vertical facets (Top to Bottom)
new_plot_labels_order <- c("LBD_ATS", "LBD_AS", "LBD_S") 

# Map short names to full comparison names, ensuring the order is maintained
plot_order <- sapply(new_plot_labels_order, function(short_name) {
 if (short_name == "LBD_S") {
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
# Function to read DEG files (ORIGINAL FUNCTION)
# ------------------------------
read_deg_file <- function(celltype, group1, group2, main_output_dir) {
  # Construct the expected file name
  fname <- file.path(main_output_dir, celltype,
                     sprintf("DEG_%s_%s_vs_%s.tsv", celltype, group1, group2))
  
  # Check if the file exists
  if (!file.exists(fname)) {
    warning("DEG file not found: ", fname)
    return(NULL)
  }
  
  # Read the tab-separated file into a data frame
  df <- read_tsv(fname, show_col_types = FALSE)
  
  # Check for required columns
  required <- c("gene", "avg_log2FC", "p_val_adj")
  if (!all(required %in% colnames(df))) {
    stop("DEG file ", fname, " missing required columns: ", paste(required, collapse = ", "))
  }
  
  # Select only the required columns and return the data frame
  return(df %>% select(all_of(required)))
}

# ------------------------------
# Main Execution: Read files into a list of dataframes (ORIGINAL READING LOOP)
# ------------------------------

# Initialize a list to store all dataframes
all_degs <- list()

for (ct in cell_types) {
  message("Reading files for cell type: ", ct)
  
  # Initialize a list for the current cell type
  all_degs[[ct]] <- list()
  
  for (pair in comparisons_to_read) {
    comp_nm <- paste0(pair[1], "_vs_", pair[2])
    
    # Read the file
    df <- read_deg_file(ct, pair[1], pair[2], main_output_dir)
    
    # Save the dataframe to the nested list using the comparison name
    if (!is.null(df)) {
      all_degs[[ct]][[comp_nm]] <- df
      message("  Successfully read: ", comp_nm, " (Rows: ", nrow(df), ")")
    } else {
      # If the file wasn't found, store an empty tibble to indicate it
      all_degs[[ct]][[comp_nm]] <- tibble()
    }
  }
}
message("\nFinished reading all DEG files.")

# ------------------------------
# Utility Function: fromList (ORIGINAL FUNCTION)
# ------------------------------
fromList <- function (input) {
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  row.names(data) <- elements
  return(data)
}

# ------------------------------
# Function: Generate Upset Plot (ORIGINAL FUNCTION)
# ------------------------------
generate_upset_plot <- function(data, celltype, reg_status, comparison_names, dot_color) {
  
  # Determine plot title
  plot_title <- paste(celltype, reg_status, "")
  
  # Convert comparison names to short names for the plot labels
  short_comp_names <- sapply(strsplit(comparison_names, "_vs_"), `[`, 1)
  
  # ComplexUpset queries for coloring the matrix dots
  upset_queries <- lapply(short_comp_names, function(set_name) {
    upset_query(set=set_name, fill=dot_color)
  })
  
  # Intersections list (using short names for sets)
  sets_list <- list(
    short_comp_names[1], short_comp_names[2], short_comp_names[3],
    c(short_comp_names[1], short_comp_names[2], short_comp_names[3]),
    c(short_comp_names[3], short_comp_names[2]),
    c(short_comp_names[1], short_comp_names[3]),
    c(short_comp_names[1], short_comp_names[2])
  )
  
  # Rename data columns to short names for Upset plot
  colnames(data) <- short_comp_names
  
  # Create the plot
  upset_gene <- upset(
    data, 
    short_comp_names,
    set_sizes = FALSE,
    themes = upset_modify_themes(
      list('intersections_matrix'=theme(text=element_text(size=9), plot.margin = margin(0, 0, 0, 0, "cm")),
           'overall_sizes' = theme(axis.text.x = element_text(size =9), plot.margin = margin(0, 0, 0, 0, "cm")))
    ),
    queries = upset_queries,
    intersections = sets_list,
    base_annotations = list(
      'Intersection size' = (
        intersection_size(
          size = 2,
          text = list(size = 2.5),
          text_mapping = aes(),
          bar_number_threshold = 3,
          width = 0.3,
          mapping = aes(fill = 'bars_color')
        )
        + scale_fill_manual(values = c('grey'), guide = 'none')
        + scale_y_continuous(expand = expansion(mult = c(0, 0.5)))
        + theme(axis.text.y = element_text(size = 8),
                axis.title.y = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.margin = margin(0, 0, 0, 0.1, "cm"),
                axis.line = element_line(colour = 'black'))
      )
    ),
    matrix = intersection_matrix(
      geom = geom_point(
        shape = 'circle filled',
        size = 2,
        stroke = 0.45,
        color = 'black' 
      )) +
      theme(axis.text.y = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            plot.margin = margin(0, 0, 0, 0, "cm")),
    sort_sets = FALSE,
    sort_intersections = FALSE,
    name = plot_title
  )
  
  return(upset_gene)
}

# ------------------------------
# UPDATED FUNCTION: Generate DEG Count Bar Plot
# ------------------------------

generate_deg_count_plot <- function(deg_data, current_celltype, comparison_order) {
  
  celltype_data <- deg_data %>% 
    filter(cell_type == current_celltype) %>%
    mutate(plot_label = factor(plot_label, levels = comparison_order))
  
  if (nrow(celltype_data) == 0 || sum(celltype_data$n) == 0) {
    return(ggplot() + 
             geom_text(aes(x=0.5, y=0.5, label=paste(current_celltype, "DEGs: No DEGs found")), size=3) + 
             theme_void())
  }
  
  # Calculate the maximum absolute count
  max_count <- max(celltype_data$n, na.rm = TRUE)
  max_abs_x <- ceiling(max_count * 1.5 / 10) * 10 
  
  p_deg_counts <- ggplot(celltype_data, aes(x = value, y = cell_type, fill = direction)) +
    geom_col(width = 0.7) +
    # Add labels for counts
    geom_text(
      data = celltype_data %>% filter(n > 0),
      aes(
        y = cell_type, 
        label = as.character(n), 
        x = value
      ),
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
    # Facet by row (vertical stacking)
    facet_grid(plot_label ~ ., switch = "y") + 
    labs(title = paste(current_celltype, "DEG Counts"),
         x = "Number of DEGs", 
         y = NULL, 
         fill = "Direction") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_text(size = 8),
      axis.ticks.x = element_blank(),
      # KEY CHANGE 1: Remove redundant Y-axis (cell type) labels
      axis.title.y = element_blank(), 
      axis.text.y = element_blank(), # Removed for visual alignment reliance
      axis.ticks.y = element_blank(),
      
      plot.title = element_text(size = 8, hjust = 0.5),
      
      panel.grid.major.x = element_line(linetype = "dashed", color = "lightgray"), 
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "gray90", color = NA),
      strip.text = element_text(face = "bold", size = 7.5),
      strip.placement = "outside",
      
      # KEY CHANGE 2: Small legend above the plot
      legend.position = "top", 
      legend.direction = "horizontal",
      legend.box.margin = margin(-.2, -.2, -.2, -.2, "cm"),
      legend.key.size = unit(0.5, "lines"), 
      legend.text = element_text(size = 7), 
      legend.title = element_text(size = 7),
      plot.margin = margin(0.2, 0.15, 0.1, 0.2, "cm")
    )
  
  return(p_deg_counts)
}

# ------------------------------
# Main Execution: Calculate all DEG counts
# ------------------------------

safe_read_deg <- function(path) {
  tryCatch(read_tsv(path, show_col_types = FALSE),
           error = function(e) tryCatch(read_csv(path, show_col_types = FALSE),
                                        error = function(e2) read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE)))
}

all_files <- list.files(main_output_dir, pattern = "^DEG_.*_vs_.*\\.(txt|csv|tsv)$", recursive = TRUE, full.names = TRUE)

if(length(all_files) == 0) stop("No DEG files found under main_output_dir with pattern 'DEG_*_vs_*.txt'")

file_df <- tibble(path = all_files) %>%
  mutate(bn = basename(path)) %>%
  mutate(noext = tools::file_path_sans_ext(bn)) %>%
  mutate(after_deg = str_remove(noext, "^DEG_")) %>%
  mutate(cell_type = word(after_deg, 1, sep = fixed("_")),
         rest = str_remove(after_deg, paste0("^", cell_type, "_")),
         comparison = rest) %>%
  select(path, cell_type, comparison) %>%
  filter(cell_type %in% cell_types) %>%
  filter(comparison %in% plot_order)

deg_long <- file_df %>%
  mutate(data = map(path, safe_read_deg)) %>%
  filter(map_lgl(data, ~ !is.null(.x))) %>%
  select(-path) %>%
  tidyr::unnest(data) %>%
  dplyr::rename_with(~ tolower(.x)) %>%
  dplyr::mutate(
    p_val_adj = dplyr::coalesce(
      .data$p_val_adj
    ),
    avg_log2fc = dplyr::coalesce(
      .data$avg_log2fc
    )
  ) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      !is.na(p_val_adj) & !is.na(avg_log2fc) & p_val_adj <= padj_cutoff & avg_log2fc >= up_logfc_cutoff~ "Up",
      !is.na(p_val_adj) & !is.na(avg_log2fc) & p_val_adj <= padj_cutoff & avg_log2fc <= down_logfc_cutoff ~ "Down",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(direction)) %>%
  dplyr::count(comparison, cell_type, direction, name = "n")

# Ensure every cell_type x comparison x direction exists (fill zeros)
all_combinations <- expand_grid(
  comparison = plot_order, 
  cell_type = cell_types,
  direction = c("Down", "Up")
)

deg_counts_all <- all_combinations %>%
  left_join(deg_long, by = c("comparison","cell_type","direction")) %>%
  mutate(n = replace_na(n, 0)) %>%
  mutate(value = if_else(direction == "Down", -n, n)) %>%
  mutate(cell_type = factor(cell_type, levels = cell_types)) %>% 
  left_join(comparison_map, by = "comparison") %>%
  mutate(plot_label = factor(plot_label, levels = new_plot_labels_order)) 

message("\nFinished calculating all DEG counts.")

# ------------------------------
# Main Execution: Filter, Prepare, Plot, and Save UpSet data (UPDATED COMBINE STEP)
# ------------------------------

for (ct in cell_types) {
  message("\nProcessing cell type: ", ct, "...")
  
  # 1. Filter, plot, and save UP-regulated genes 
  up_list_input <- list()
  for (comp in comparison_names) { 
    df <- all_degs[[ct]][[comp]]
    if (nrow(df) > 0) {
      up_genes <- df %>%
        filter(avg_log2FC >= up_logfc_cutoff & p_val_adj <= padj_cutoff) %>%
        pull(gene)
      short_name <- strsplit(comp, "_vs_")[[1]][1]
      up_list_input[[short_name]] <- up_genes
    } else {
      short_name <- strsplit(comp, "_vs_")[[1]][1]
      up_list_input[[short_name]] <- character(0)
    }
  }
  
  up_data <- fromList(up_list_input)
  write.table(
    up_data,
    file.path(output_dir, paste0(ct, "_upregulated.txt")),
    sep = "\t", quote = FALSE, col.names = NA
  )
  message("Saved UP-regulated data table for ", ct)
  
  if(nrow(up_data) > 0) {
    p_up <- generate_upset_plot(up_data, ct, "upregulated", comparison_names, "red")
  } else {
    message("Skipping UP-regulated UpSet plot object for ", ct, " (No overlapping DEGs).")
    p_up <- ggplot() + geom_text(aes(x=0.5, y=0.5, label=paste(ct, "UP-regulated: No DEGs found")), size=4) + theme_void() + ggtitle(paste(ct, "UP-regulated"))
  }
  
  # 2. Filter, plot, and save DOWN-regulated genes
  down_list_input <- list()
  for (comp in comparison_names) {
    df <- all_degs[[ct]][[comp]]
    if (nrow(df) > 0) {
      down_genes <- df %>%
        filter(avg_log2FC <= down_logfc_cutoff & p_val_adj <= padj_cutoff) %>%
        pull(gene)
      short_name <- strsplit(comp, "_vs_")[[1]][1]
      down_list_input[[short_name]] <- down_genes
    } else {
      short_name <- strsplit(comp, "_vs_")[[1]][1]
      down_list_input[[short_name]] <- character(0)
    }
  }
  
  down_data <- fromList(down_list_input)
  write.table(
    down_data,
    file.path(output_dir, paste0(ct, "_downregulated.txt")),
    sep = "\t", quote = FALSE, col.names = NA
  )
  message("Saved DOWN-regulated data table for ", ct)
  
  if(nrow(down_data) > 0) {
    p_down <- generate_upset_plot(down_data, ct, "downregulated", comparison_names, "blue")
  } else {
    message("Skipping DOWN-regulated UpSet plot object for ", ct, " (No overlapping DEGs).")
    p_down <- ggplot() + geom_text(aes(x=0.5, y=0.5, label=paste(ct, "DOWN-regulated: No DEGs found")), size=4) + theme_void() + ggtitle(paste(ct, "DOWN-regulated"))
  }
  
  # 3. GENERATE DEG COUNT PLOT
  p_deg_counts <- generate_deg_count_plot(deg_counts_all, ct, new_plot_labels_order)
  message("Generated DEG count bar plot for ", ct)
  
  # 4. Combine and Save the three plots into a single PDF (Horizontal Arrangement)
  combined_plot_fname <- file.path(output_dir, paste0(ct, "_Combined_UpSet_Horizontal_Counts.pdf"))
  message("Saving combined DEG Count and UpSet plot to: ", combined_plot_fname)
  
  # Arrange all three plots on one row (A, B, C)
  p_combined <- plot_grid(
    p_deg_counts, # Panel A
    p_down,       # Panel B
    p_up,         # Panel C
    ncol = 3,
    rel_widths = c(1, 1.25, 1.25), # Adjust widths to prioritize the UpSet complexity
    labels = c("A", "B", "C"),
    label_size = 10,
    align = 'v' # Align plots vertically (top edge alignment)
  )
  
  # Save the final combined PDF
  pdf(combined_plot_fname, width = 7, height = 2.75) 
  print(p_combined)
  dev.off()
}

message("\nAll cell types processed.")
message("Combined plots (PDF) and data tables (TXT) are saved in: ", output_dir)

# ------------------------------
# Main Execution: Generate Excel File with DEG details (UNCHANGED)
# ------------------------------

# Define the output Excel file path
excel_output_path <- file.path(output_dir, "All_CellType_DEGs_Summary_Lists.xlsx") 
wb <- createWorkbook() 

message("\nGenerating Excel Summary file: ", excel_output_path)

for (ct in cell_types) {
  message("  Preparing data for cell type: ", ct)
  
  for (reg_status in c("upregulated", "downregulated")) {
    
    # Set filtering parameters based on regulation status
    logfc_filter_fn <- if (reg_status == "upregulated") {
      function(df) filter(df, avg_log2FC >= up_logfc_cutoff & p_val_adj <= padj_cutoff)
    } else {
      function(df) filter(df, avg_log2FC <= down_logfc_cutoff & p_val_adj <= padj_cutoff)
    }
    
    # 1. Initialize a list to hold the filtered gene vectors
    comp_gene_lists <- list()
    
    # 2. Process each comparison for the current cell type and regulation status
    for (comp in comparison_names) {
      df <- all_degs[[ct]][[comp]]
      
      if (nrow(df) > 0) {
        
        # Apply the regulation filter and extract only the gene column
        gene_list <- logfc_filter_fn(df) %>%
          pull(gene) %>%
          as.character() 
        
        short_name <- strsplit(comp, "_vs_")[[1]][1]
        
        # Store the gene list using the short name
        comp_gene_lists[[short_name]] <- gene_list
      }
    }
    
    # 3. Combine gene lists into a data frame (padding with NA)
    if (length(comp_gene_lists) > 0) {
      
      max_len <- if(length(comp_gene_lists) > 0) max(sapply(comp_gene_lists, length)) else 0
      
      if (max_len > 0) {
        # Pad shorter lists with NA to match the maximum length
        padded_lists <- lapply(comp_gene_lists, function(x) {
          if (length(x) < max_len) {
            length(x) <- max_len 
          }
          return(x)
        })
        
        # Convert the list of padded vectors into a dataframe
        combined_df <- as.data.frame(padded_lists, stringsAsFactors = FALSE)
        
        # 4. Create a sheet name and add to the workbook
        sheet_name <- paste0(ct, "_", ifelse(reg_status == "upregulated", "UP", "DOWN"))
        addWorksheet(wb, sheet_name)
        
        # 5. Add the combined data to the sheet
        writeData(wb, sheet_name, combined_df, startRow = 1, startCol = 1)
        message("    - Added sheet: ", sheet_name, " (Max rows: ", max_len, ")")
      } else {
        message("    - Skipping ", reg_status, " sheet for ", ct, " (No DEGs found in any comparison).")
      }
      
    } else {
      message("    - Skipping ", reg_status, " sheet for ", ct, " (No combined DEGs found).")
    }
  }
}

# 6. Save the entire workbook
tryCatch({
  saveWorkbook(wb, excel_output_path, overwrite = TRUE)
  message("\nSuccessfully saved Excel file: ", excel_output_path)
}, error = function(e) {
  stop("Error saving Excel file: ", e$message)
})

# ------------------------------
# End of Script
# ------------------------------