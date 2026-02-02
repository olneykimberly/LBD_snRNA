#!/usr/bin/env Rscript
# Nuclear fraction 

## ----working_directory----------------------------------------------------------
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")
## ----dataObject-----------------------------------------------------------------
project_ID <- "CWOW_cellbender_singlets_scDblFinder_exprate"
dataObject <- readRDS(file = paste0("../rObjects/", project_ID, ".rds"))

# output dir
out_dir <- "../results/nuclear_fraction"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ----read_nuclear_fraction------------------------------------------------------
# Seurat cell names look like: "BR_Nuclei_0373_<barcode>-1"
# So prefix must be Sample_BR (not Sample_ID)
sample_brs <- unique(dataObject$Sample_BR)

read_nf_file <- function(path_file) {
  # Try robust read: any whitespace delimiter
  df <- tryCatch(
    read.table(path_file, header = TRUE, sep = "", stringsAsFactors = FALSE, comment.char = ""),
    error = function(e) NULL
  )
  
  # If header parsing failed, try without header
  if (is.null(df) || ncol(df) < 2) {
    df <- read.table(path_file, header = FALSE, sep = "", stringsAsFactors = FALSE, comment.char = "")
  }
  
  # Normalize column names
  cn <- tolower(colnames(df))
  colnames(df) <- cn
  
  # If there is no header, assign expected names by position
  # Expect 2 columns: nuclear_fraction, barcode
  if (!("barcode" %in% colnames(df)) && ncol(df) >= 2) {
    # If the first column looks numeric and the second looks like barcodes, assume order is correct
    colnames(df)[1:2] <- c("nuclear_fraction", "barcode")
  }
  
  # Final checks
  if (!("barcode" %in% colnames(df))) {
    stop("Could not find a 'barcode' column in ", path_file,
         ". Columns seen: ", paste(colnames(df), collapse = ", "))
  }
  if (!("nuclear_fraction" %in% colnames(df))) {
    stop("Could not find a 'nuclear_fraction' column in ", path_file,
         ". Columns seen: ", paste(colnames(df), collapse = ", "))
  }
  
  df$nuclear_fraction <- as.numeric(df$nuclear_fraction)
  if (anyNA(df$nuclear_fraction)) {
    # Allow a few NAs, but warn
    warning("Some nuclear_fraction values are NA after numeric coercion in ", path_file)
  }
  
  return(df[, c("nuclear_fraction", "barcode")])
}

nuclear_fraction_list <- list()

for (s_br in sample_brs) {
  path_file <- paste0("../nuclear_fraction/", s_br, "_nuclear_fraction.tsv")
  
  if (file.exists(path_file)) {
    df_nf <- read_nf_file(path_file)
    
    # Build full Seurat-style cell names
    df_nf$barcode_full <- paste0(s_br, "_", df_nf$barcode)
    rownames(df_nf) <- df_nf$barcode_full
    
    nuclear_fraction_list[[s_br]] <- df_nf
    message("Loaded nuclear fraction: ", s_br, " | n=", nrow(df_nf))
  } else {
    message("Warning: Nuclear fraction file not found for ", s_br, " at ", path_file)
  }
}

if (length(nuclear_fraction_list) == 0) {
  stop("No nuclear fraction files were found in ../nuclear_fraction/. Check filenames and paths.")
}

## ----add_nuclear_fraction_to_object--------------------------------------------
df_nf_combined <- do.call(rbind, nuclear_fraction_list)

# Keep only rows that correspond to cells in the Seurat object
df_nf_combined <- df_nf_combined[df_nf_combined$barcode_full %in% Cells(dataObject), , drop = FALSE]

# Reorder to match Seurat cell order (will introduce NAs for missing cells)
df_nf_reordered <- df_nf_combined[Cells(dataObject), , drop = FALSE]

message("NAs in nuclear_fraction after merge (expected if some barcodes missing):")
print(table(is.na(df_nf_reordered$nuclear_fraction)))

# Add to Seurat metadata (leave NA as NA; don't set to 0)
dataObject$nuclear_fraction <- df_nf_reordered$nuclear_fraction

df_combined_qc_data <- dataObject[[]]

# cleanup
rm(df_nf, df_nf_combined, df_nf_reordered, nuclear_fraction_list, df_nf_combined)

## ----plot_nuclear_fraction------------------------------------------------------
den <- ggplot(df_combined_qc_data, aes(x = nuclear_fraction, color = Sample_ID)) +
  geom_density(na.rm = TRUE) +
  labs(
    title = "Nuclear fraction",
    x = "Nuclear fraction",
    y = "Density",
    color = "Sample ID"
  ) +
  theme_minimal()

pdf(file.path(out_dir, paste0(project_ID, "_density.pdf")), width = 8, height = 5)
print(den)
dev.off()

## ----identify_empty_drops-------------------------------------------------------
df_nf_umi <- data.frame(
  nf = df_combined_qc_data$nuclear_fraction,
  umi = df_combined_qc_data$nCount_RNA,
  Sample_ID = df_combined_qc_data$Sample_ID
)

df_cell_status <- identify_empty_drops(nf_umi = df_nf_umi)

message("Cell status overall:")
print(table(df_cell_status$cell_status))

message("Cell status by sample:")
table_cell_status_sample <- table(df_cell_status$cell_status, df_cell_status$Sample_ID)
print(table_cell_status_sample)

write.table(
  table_cell_status_sample,
  file.path(out_dir, paste0(project_ID, "_table_cell_status_by_sample.tsv")),
  sep = "\t",
  quote = FALSE
)

pdf(file.path(out_dir, paste0(project_ID, "_summary.pdf")), width = 10, height = 10)
identify_empty_drops(nf_umi = df_nf_umi, include_plot = TRUE)
dev.off()

## ----save_updated_object--------------------------------------------------------
saveRDS(
  dataObject,
  file = paste0("../rObjects/", project_ID, "_with_nuclear_fraction.rds"),
  compress = FALSE
)

message("Done. Saved updated object: ",
        paste0("../rObjects/", project_ID, "_with_nuclear_fraction.rds"))
