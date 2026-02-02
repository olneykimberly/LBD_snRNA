#!/usr/bin/env Rscript

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

library(readr)
library(stringr)
library(purrr)

# Prefer writexl for "plain" xlsx that openpyxl can read robustly
suppressPackageStartupMessages({
  ok <- requireNamespace("writexl", quietly = TRUE)
})
if (!ok) {
  stop("Package 'writexl' is not installed. Install it (once) with: install.packages('writexl')")
}

# ---- output locations ----
staging_dir <- "../results/metascape_hdWGCNA_modules"   # where your batch runner looks
dir.create(staging_dir, showWarnings = FALSE, recursive = TRUE)

# If you still want these:
txt_out_dir <- file.path(staging_dir, "hdWGCNA_module_txtlists")
dir.create(txt_out_dir, showWarnings = FALSE, recursive = TRUE)

modules_dir <- "../results/hdWGCNA"

sanitize_colname <- function(x) {
  x %>%
    str_replace_all("-", "_") %>%                 # your original intent
    str_replace_all("[^A-Za-z0-9_]+", "_") %>%    # only safe chars
    str_replace_all("^_+|_+$", "") %>%
    str_replace_all("_+", "_")
}

for (ct in order_cell_type) {
  
  input_file <- file.path(modules_dir, paste0("modules_", ct, ".tsv"))
  if (!file.exists(input_file)) {
    warning("File not found: ", input_file)
    next
  }
  
  df <- readr::read_tsv(input_file, show_col_types = FALSE)
  
  if (!all(c("gene_name", "module") %in% colnames(df))) {
    warning("Missing required columns gene_name/module in: ", input_file)
    next
  }
  
  # clean inputs
  df$gene_name <- str_trim(as.character(df$gene_name))
  df$module    <- str_trim(as.character(df$module))
  df <- df[!is.na(df$gene_name) & df$gene_name != "" & !is.na(df$module) & df$module != "", , drop = FALSE]
  
  # split genes by module
  module_list <- split(df$gene_name, df$module)
  
  # per-module cleanup: trim, drop blanks, unique
  module_list <- lapply(module_list, function(genes) {
    genes <- str_trim(as.character(genes))
    genes <- genes[!is.na(genes) & genes != ""]
    unique(genes)
  })
  module_list <- module_list[lengths(module_list) > 0]
  
  if (length(module_list) == 0) {
    warning("No genes after cleaning for cell type: ", ct)
    next
  }
  
  # clean/safe column names
  coln <- sanitize_colname(names(module_list))
  coln <- make.unique(coln, sep = "_")
  names(module_list) <- coln
  
  # pad to rectangle
  max_len <- max(lengths(module_list))
  padded_list <- lapply(module_list, function(x) { length(x) <- max_len; x })
  df_excel <- as.data.frame(padded_list, check.names = FALSE, stringsAsFactors = FALSE)
  
  # ---- write CLEAN XLSX via writexl ----
  xlsx_out <- file.path(staging_dir, paste0("hdWGCNA_", ct, "_modules.xlsx"))
  writexl::write_xlsx(list(GeneLists = df_excel), path = xlsx_out)
  
  # ---- ALSO write TSV fallback (same multi-column layout) ----
  tsv_out <- file.path(staging_dir, paste0("hdWGCNA_", ct, "_modules.tsv"))
  readr::write_tsv(df_excel, file = tsv_out, na = "")
  
  message("Created: ", xlsx_out)
  message("Created: ", tsv_out)
  
  # ---- optional txt lists per module ----
  purrr::iwalk(module_list, function(genes, module_name) {
    out_txt <- file.path(txt_out_dir, paste0("hdWGCNA_", ct, "_", module_name, ".txt"))
    readr::write_lines(genes, out_txt)
  })
}
