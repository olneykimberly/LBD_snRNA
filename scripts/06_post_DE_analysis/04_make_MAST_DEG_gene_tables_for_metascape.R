setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")
project_ID <- "CWOW_cellbender"
order_cell_type

qval = 0.05
posFC = .25
negFC = -.25
FC = .25

comparisons <- list(
  c("AD_AT", "CONTROL"),
  c("LBD_S", "CONTROL"),
  c("LBD_AS", "CONTROL"),
  c("LBD_ATS", "CONTROL"),
  c("LBD_S", "AD_AT"),
  c("LBD_AS", "AD_AT"),
  c("LBD_ATS", "AD_AT"), 
  c("LBD_AS", "LBD_S"),
  c("LBD_ATS", "LBD_S"),
  c("LBD_ATS", "LBD_AS")
)

output_dir <- "../results/metascape/MAST_DEGs_"

for (cell_type in order_cell_type) {
  
  for (comp in comparisons) {
    
    pathology_A <- comp[1]
    pathology_B <- comp[2]
    
    # Construct the file name based on the specified format
    file_name <- paste0(cell_type,"/DEG_", cell_type, "_", pathology_A, "_vs_", pathology_B, ".tsv")
    file_path <- paste0("../results/DEGs_RNA_pct0.25/", file_name)
    plot_title <- paste0(cell_type, " (", pathology_A, " vs ", pathology_B, ")")
    output_pdf <- paste0("../results/volcanoes_MAST/", cell_type, "_", pathology_A, "_vs_", pathology_B, ".pdf")
    
    # Check if the file exists before attempting to read
    if (!file.exists(file_path)) {
      message(paste("Skipping:", file_name, "- File not found."))
      next # Move to the next comparison/cell_type
    }
    
    # 1. Read the DEG results
    data <- read.delim(file_path)
    data <- na.omit(data) # Ensure no NA rows
    
    # 2. Check for the required columns (assuming 'gene' column exists) # p_val	avg_log2FC	pct.1	pct.2	p_val_adj	gene
    if (!all(c("avg_log2FC", "p_val", "p_val_adj", "gene") %in% colnames(data))) {
      message(paste("Skipping:", file_name, "- Missing required columns (avg_log2FC, p_val, p_val_adj, or gene)."))
      next
    }
    significant_up_DEGs <- subset(data, p_val_adj < 0.05 & avg_log2FC > 0.25) 
    significant_down_DEGs <- subset(data, p_val_adj < 0.05 & avg_log2FC < -0.25)
    
    write.table(significant_up_DEGs$gene, paste0(output_dir, cell_type, "_", pathology_A, "_vs_", pathology_B, "_upregulated.txt"), row.names = FALSE, quote = FALSE)
    write.table(significant_down_DEGs$gene, paste0(output_dir, cell_type, "_", pathology_A, "_vs_", pathology_B, "_downregulated.txt"), row.names = FALSE, quote = FALSE)
  }
}

#--- multi gene list saved as CLEAN excel for Metascape (openpyxl-safe)
library(writexl)

# --- 1. Setup Paths ---
input_dir   <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/metascape"
staging_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/metascape_xlsx_staging"

if (!dir.exists(staging_dir)) dir.create(staging_dir, recursive = TRUE)

# --- 2. Define Groups ---
groups <- list(
  LBD_types_vs_CONTROL = c("LBD_S_vs_CONTROL", "LBD_AS_vs_CONTROL", "LBD_ATS_vs_CONTROL"),
  Disease_vs_CONTROL   = c("AD_AT_vs_CONTROL", "LBD_S_vs_CONTROL", "LBD_AS_vs_CONTROL", "LBD_ATS_vs_CONTROL"),
  LBD_types_vs_AD      = c("LBD_S_vs_AD_AT", "LBD_AS_vs_AD_AT", "LBD_ATS_vs_AD_AT")
)

# --- 3. Identify Cell Types ---
all_files  <- list.files(input_dir, pattern = "^MAST_DEGs_.*\\.txt$", full.names = FALSE)
cell_types <- unique(vapply(strsplit(all_files, "_"), function(x) x[3], character(1)))

# --- helper: read genes robustly ---
read_gene_file <- function(path) {
  if (!file.exists(path)) return(character(0))
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

# --- 4. Loop and Create CLEAN Excel Files ---
for (cell in cell_types) {
  for (group_name in names(groups)) {
    for (direction in c("upregulated", "downregulated")) {
      
      combined_list <- list()
      max_len <- 0L
      
      for (comp in groups[[group_name]]) {
        file_path <- file.path(input_dir, paste0("MAST_DEGs_", cell, "_", comp, "_", direction, ".txt"))
        genes <- read_gene_file(file_path)
        
        if (length(genes) > 0) {
          combined_list[[comp]] <- genes
          max_len <- max(max_len, length(genes))
        }
      }
      
      if (length(combined_list) == 0) next
      
      # pad to rectangular data.frame (Metascape multigene columns)
      padded_list <- lapply(combined_list, function(x) {
        length(x) <- max_len
        x
      })
      
      df <- as.data.frame(padded_list, stringsAsFactors = FALSE, check.names = FALSE)
      
      # Write a minimal XLSX: one worksheet, no styling, no drawings
      file_name <- paste0(cell, "_", group_name, "_", direction, ".xlsx")
      out_path  <- file.path(staging_dir, file_name)
      
      writexl::write_xlsx(x = list(GeneLists = df), path = out_path)
      
      message("Generated (clean): ", file_name)
    }
  }
}
