## ----working_directory-------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
## ----echo=FALSE, message=FALSE-----------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "path_files_and_colours.R"))

## ----dataObject---------------------------------------------------------------------------------------------------------------------------
project_ID <- "CWOW_cellbender_filtered"
dataObject <- readRDS(file = paste0("../rObjects/", project_ID, ".rds"))

## ----read_nuclear_fraction---------------------------------------------------------------------------------------------------------------------------
nuclear_fraction_list <- list()
for (sample_id in order_samples) {
  path_file <- paste0("../nuclear_fraction/", sample_id, "_nuclear_fraction.tsv")
  if (file.exists(path_file)) {
    df_nf <- read.table(path_file, header = TRUE)
    df_nf$barcode <- paste0(sample_id, "_", df_nf$barcode)
    rownames(df_nf) <- df_nf$barcode
    nuclear_fraction_list[[sample_id]] <- df_nf
    print(paste0("Nuclear fraction data for ", sample_id))
  } else {
    print(paste0("Warning: Nuclear fraction file not found for ", sample_id))
  }
}

## ----add_nuclear_fraction_to_object---------------------------------------------------------------------------------------------------------------------------
# Combine all data frames from the list into a single data frame
df_nf_combined <- do.call(rbind, nuclear_fraction_list)
df_nf_combined_filtered <- df_nf_combined[df_nf_combined$barcode %in% Cells(dataObject), ]
rownames(df_nf_combined_filtered) <- df_nf_combined_filtered$barcode

# The order of the rows in the metadata must match the order of the cells in the object
df_nf_combined_reordered <- df_nf_combined_filtered[Cells(dataObject), , drop = FALSE]
# Is nuclear fraction NA? 
table(is.na(df_nf_combined_reordered$nuclear_fraction))
# Replace NA with zero
df_nf_combined_reordered$nuclear_fraction[is.na(df_nf_combined_reordered$nuclear_fraction)] <- 0
# Add the new metadata column to the main Seurat object
dataObject$nuclear_fraction <- df_nf_combined_reordered$nuclear_fraction
df_combined_qc_data <- dataObject@meta.data # create df of metadata 

# clean up
rm(df_nf, df_nf_combined, df_nf_combined_filtered, df_nf_combined_reordered, nuclear_fraction_list)
## ----plot_nuclear_fraction---------------------------------------------------------------------------------------------------------------------------
den <- ggplot(df_combined_qc_data, aes(x = nuclear_fraction, color = Sample_ID)) +
  geom_density() +
  labs(title = "Nuclear fraction",
       x = "Nuclear fraction",
       y = "Density",
       color = "Sample ID") +
  theme_minimal()

pdf(paste0("../results/nuclear_fraction/", project_ID, "_density.pdf"), width = 7, height = 5)
den
dev.off()

## ----identify_empty_drops---------------------------------------------------------------------------------------------------------------------------
# Get data frame with the nuclear fraction and umi counts 
df_nf_umi <- data.frame(nf=df_combined_qc_data$nuclear_fraction,
                         umi=df_combined_qc_data$nCount_RNA, 
                         Sample_ID = df_combined_qc_data$Sample_ID)

# Run identify_empty_drops
df_cell_status <- identify_empty_drops(nf_umi=df_nf_umi)
# Table of cell status
table(df_cell_status$cell_status)
# Table of cell status by sample ID
table_cell_status_sample <- table(gbm.ed$cell_status, gbm.ed$Sample_ID)
table_cell_status_sample
# Save the table to a CSV file
write.table(table_cell_status_sample, paste0("../results/nuclear_fraction/", project_ID, "_table_cell_status_by_sample.tsv"), sep = "\t", quote = FALSE)

pdf(paste0("../results/nuclear_fraction/", project_ID, "_summary.pdf"), width = 10, height = 10)
identify_empty_drops(nf_umi=df_nf_umi, include_plot = TRUE)
dev.off()