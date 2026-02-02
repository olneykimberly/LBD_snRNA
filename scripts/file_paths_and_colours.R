#----------------- Libraries
set.seed(28)
.libPaths(c("/tgen_labs/jfryer/kolney/R/x86_64-pc-linux-gnu-library/4.3", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
.libPaths()

library(Matrix)
packageVersion("Matrix")
find.package("Matrix")

# suppress startup messages for ALL package loads + remove duplicates
suppressPackageStartupMessages({
 # library(Matrix, lib.loc = "/usr/local/lib/R/site-library")
 # library(Matrix)
  #library(reticulate)
  library(SeuratObject)
  library(Signac)
  library(Seurat)
  library(stringr)
  library(ggplot2)
  library(harmony)
  library(remaCor)
  library(gridExtra)
  library(grid)
  library(lattice)
  library(R.utils)
  library(SeuratWrappers)
  library(Azimuth)
  library(dittoSeq)
  library(dplyr)
  library(RColorBrewer)
  library(DESeq2) # adds matrix
  library(openxlsx)
  library(ggrepel)
  library(glmGamPoi)
  library(devtools)
  library(DoubletFinder)
  library(reshape2)
  #library(ggtree)
  library(BiocParallel)
  library(edgeR)
  library(limma)
  library(gplots)
  library(grDevices)
  library(scales)
  #library(tximport)
  library(tidyverse)
  library(GenomicFeatures)
  library(plyr)
  library(data.table)
  library(readxl)
  library(pheatmap)
  #library(NatParksPalettes)
  library(UpSetR)
  library(cowplot)
  library(ggpubr)
  library(patchwork) # For combining plots
  library(SeuratData)
  #library(batchelor)
  library(DropletQC)
  library(Rsamtools)
  library(GenomicRanges)
  library(ComplexUpset)
  library(UCell)
  library(GeneOverlap)
  #BiocManager::install(c("WGCNA", "UCell", "GenomicRanges", "GeneOverlap"))
  # co-expression network analysis packages:
  library(WGCNA)
  #devtools::install_github("immunogenomics/harmony")
  #devtools::install_github('smorabit/hdWGCNA', ref='dev')
  library(hdWGCNA)
  
  # (kept from your original; duplicates already removed above)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(tidyr)
  library(readr)
  library(scuttle)
  library(BiocSingular)
  library(future)
  library(future.apply)
  library(BiocParallel)
  
})

#--- variables
# paths, colors, shapes and more
color_control <- "#4682B4" 
color_AD <- "#B4464B" 
color_LBD <- "gray35" 
shape_control <- c(15) # square
shape_AD <- c(16) # circle
shape_PA <- c(17) # triangle
shape_LBD <- c(18) # diamond

color_type <- c("#4682B4","#B4464B", "gray35")
color_ATS <- c("#4682B4","#B4464B", "gray","gray65", "gray35")
color_sex <- c("#490092", "#D55E00")
color_blind <- dittoColors()
color_panel <- dittoColors()
sample_colors <- c("#4682B4","#4682B4","#4682B4","#4682B4","#4682B4","#4682B4","#4682B4",
                   "#B4464B","#B4464B","#B4464B","#B4464B","#B4464B","#B4464B","#B4464B","#B4464B",
                   "gray","gray","gray","gray","gray","gray","gray",
                   "gray65","gray65","gray65","gray65","gray65","gray65","gray65",
                   "gray35","gray35","gray35","gray35","gray35","gray35")
order_cell_type <- c("neuron", "interneuron", "oligodendrocyte", "opc",  "astrocyte", "microglia", "endothelial", "fibroblast", "mural")

#--- references and metadata
metadata_all <-
  read.delim(
    "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/metadata/metadata_seq_info.txt")

metadata_batch <-
  read.delim(
    "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/metadata/metadata_batch.txt")
metadata_batch <- metadata_batch[, c("Sample_ID", "prep.batch", "cdna.avg.size..bp.", "library.pcr.cycles", "library.pcr.batch", "insert.size..bp.")] # Selects multiple columns

metadata <- merge(metadata_all, metadata_batch, by = "Sample_ID")


samples_to_remove <- c("LBD_AS_F4", "Ctr_F1", "LBD_S_F4", "LBD_ATS_F3", "LBD_ATS_F4") # BR_Nuclei_0404, BR_Nuclei_0376, BR_Nuclei_0407, BR_Nuclei_0394, BR_Nuclei_0403
metadata_samples_removed <- subset(metadata, !Sample_ID %in% samples_to_remove)
metadata <- metadata_samples_removed
rm(metadata_all, metadata_samples_removed, metadata_batch)

metadata$sampleID <- factor(metadata$Sample_ID, levels = c(metadata$Sample_ID)) # keep sample order 
samples <- metadata$sampleID 
order_sex <- factor(metadata$sex_inferred, levels = unique(metadata$sex_inferred))
order_disease <- factor(metadata$TYPE, levels = c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS")) # The order in which to display the groups

# obtaining the sample ID from the fastq naming
metadata <- metadata %>%
  mutate(sampleID = gsub(".*_(\\d+)_.*_(BR_Nuclei).*", "\\2_\\1", Lane.Name))
samples <- metadata$sampleID 

# Order the sampleID by disease group, so when plotting the samples are ordered by disease group 
# sampleID with disease_order
order <- metadata %>%
  arrange(order_disease) %>%
  dplyr::select(TYPE, sampleID, Sample_ID)
write.table(order, "order.txt", quote = F, row.names = F, sep = "\t") # export the order as a text file
order_disease <- order$TYPE
order_sample <- factor(order$Sample_ID, levels = order$Sample_ID)
# Convert the column to a factor with the desired order
metadata$Sample_ID <- factor(metadata$Sample_ID, levels = order_sample)

# Sort the dataframe by the newly-ordered factor column
metadata_sorted <- metadata[order(metadata$Sample_ID), ]
order_samples <- metadata_sorted$sampleID
rm(samples, order, samples_to_remove, metadata)
metadata <- metadata_sorted
rm(metadata_sorted)

metadata_all <- read.delim("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/metadata/metadata_seq_info.txt")
metadata_batch <- read.delim("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/metadata/metadata_batch.txt")

metadata_batch <- metadata_batch[, c(
  "Sample_ID", "prep.batch", "cdna.avg.size..bp.", "library.pcr.cycles",
  "library.pcr.batch", "insert.size..bp."
)]

metadata <- merge(metadata_all, metadata_batch, by = "Sample_ID")

# remove samples (these correspond to specific BR_Nuclei IDs)

samples_to_remove <- c("LBD_AS_F4", "Ctr_F1", "LBD_S_F4", "LBD_ATS_F3", "LBD_ATS_F4")
metadata <- subset(metadata, !Sample_ID %in% samples_to_remove)

rm(metadata_all, metadata_batch)

# derive BR_Nuclei_#### name from Lane.Name (this is what matches Seurat prefixes)

metadata <- metadata %>%
  dplyr::mutate(
    sampleID = gsub(
      ".*_(BR_Nuclei)_([0-9]+)_.*",
      "\\1_\\2",
      Lane.Name
    ),
    order_disease = factor(TYPE, levels = c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS"))
  )

# order samples by disease group (and preserve within-group order as currently present)

metadata <- metadata %>%
  dplyr::arrange(order_disease)

# this is the *desired plotting order* for the short IDs (Ctr_F2 etc.)
order_samples_short <- metadata$Sample_ID

# quick sanity
stopifnot(!anyDuplicated(metadata$Sample_ID))
stopifnot(!anyDuplicated(metadata$sampleID))


columns_to_fill <- c(
  "Cing.LB",
  "Braak.NFT",
  "Thal.amyloid",
  "MF.SP",
  "MF.NFT",
  "MF.LB",
  "MF.Amyloid",
  "MF.Tau",
  "Cing.Synuclein"
)

# Loop through each column and replace NA values with 0.
# The 'matched_metadata' object is assumed to be a data frame.
for (col_name in columns_to_fill) {
  # Check if the column exists in the data frame to prevent errors.
  if (col_name %in% names(metadata)) {
    metadata[[col_name]][is.na(metadata[[col_name]])] <- 0
  }
}

# The 'matched_metadata' data frame is now updated.
# You can verify the changes by checking for NA values in the columns.
# For example:
# sum(is.na(matched_metadata$Cing.LB))

# reference
path_ref = c("/tgen_labs/jfryer/projects/references/human/GRCh38/refdata-gex-GRCh38-2024-A")
gene_info <- read.delim(paste0(path_ref, "/star/geneInfo.tab"), header = FALSE)
gene_info = gene_info[-1,]
colnames(gene_info) <- c("gene_ID", "gene_name", "type")


# cell cycle 
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
genes_s <- cc.genes$s.genes
genes_g2m <- cc.genes$g2m.genes

#--- functions 
fun_color_correlation <-
  colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

fun_from_list <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

saveToPDF <- function(...) {
  d = dev.copy(pdf,...)
  dev.off(d)
}

genes_markers <-
  c(
    "CLU", 
    "GFAP", 
    "AQP4", 
    "GJA1",
    "CLDN5",
    "ADGRF5",
    "FLT1",
    "COL1A1",
    "COL1A2",
    "DCN",
    "HEXB",
    "C1QA",
    "C1QB",
    "C1QC",
    "TMEM119",
    "ITGAM",
    "TYROBP",
    "P2RY12",
    "AIF1",
    "RBFOX1",
    "RBFOX3", 
    "SNAP25",
    "SYT1",
    "GAD1",
    "GAD2",
    "PLP1",
    "MBP", 
    "MOG", 
    "OLIG1",
    "PDGFRA",
    "VCAN",
    "TNR",
    "ACTA2",
    "VTN"
  )

# ===== helper functions (compact) =====

# read DEG file and return data.frame with standardized column names present in file
read_deg <- function(cell_type, pathA, pathB) {
  fn <- file.path(main_output_dir, cell_type, sprintf("DEG_%s_%s_vs_%s.tsv", cell_type, pathA, pathB))
  if (!file.exists(fn)) {
    message("  (missing) ", fn)
    return(NULL)
  }
  df <- tryCatch(read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""),
                 error = function(e) { message("  error reading ", fn, " : ", conditionMessage(e)); NULL })
  if (is.null(df)) return(NULL)
  # verify expected columns (simple check)
  expected <- c("gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")
  if (!all(expected %in% names(df))) {
    stop("File ", fn, " does not contain the expected columns: ", paste(expected, collapse = ", "))
  }
  # ensure gene is a character vector column
  df$gene <- as.character(df$gene)
  df
}

# given a named list of sets (cell_type -> character vector of genes), make membership data.frame
membership_df <- function(sets_list) {
  sets_list <- sets_list[!vapply(sets_list, is.null, logical(1))]
  if (length(sets_list) == 0) return(NULL)
  all_genes <- sort(unique(unlist(sets_list)))
  mat <- sapply(sets_list, function(g) as.integer(all_genes %in% g), USE.NAMES = TRUE)
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  df$gene <- all_genes
  df[, c("gene", setdiff(names(df), "gene")), drop = FALSE]
}

# filter membership_df to keep only genes that belong to combinations that have >= min_size members
filter_members_by_combo_size <- function(mem_df, min_size) {
  if (is.null(mem_df)) return(NULL)
  set_cols <- setdiff(names(mem_df), "gene")
  comb <- apply(mem_df[, set_cols, drop = FALSE], 1, paste, collapse = "")
  comb_table <- table(comb)
  keep_comb <- names(comb_table)[comb_table >= min_size]
  keep_idx <- comb %in% keep_comb
  res <- mem_df[keep_idx, , drop = FALSE]
  if (nrow(res) == 0) return(NULL)
  res
}

# build shared table that includes per-cell stats columns (compact)
build_shared_table <- function(sets_list, deg_tables, min_shared = 2, direction = "up") {
  sets_list <- sets_list[!vapply(sets_list, is.null, logical(1))]
  if (length(sets_list) == 0) return(NULL)
  all_genes <- sort(unique(unlist(sets_list)))
  presence <- sapply(sets_list, function(s) as.integer(all_genes %in% s), USE.NAMES = TRUE)
  counts <- rowSums(presence)
  keep <- which(counts >= min_shared)
  if (length(keep) == 0) return(NULL)
  genes_keep <- all_genes[keep]
  res <- data.frame(gene = genes_keep, count = counts[keep], cell_types = apply(presence[keep, , drop = FALSE], 1, function(x) paste(names(sets_list)[which(x==1)], collapse = ",")), stringsAsFactors = FALSE)
  # for each cell type, add the four columns (pct.1, pct.2, avg_log2FC, p_val_adj) where available
  for (ct in cell_types) {
    res[[paste0(ct, "_pct.1")]] <- NA_real_
    res[[paste0(ct, "_pct.2")]] <- NA_real_
    res[[paste0(ct, "_avg_log2FC")]] <- NA_real_
    res[[paste0(ct, "_p_val_adj")]] <- NA_real_
    if (!is.null(deg_tables[[ct]])) {
      tbl <- deg_tables[[ct]]
      idx <- match(genes_keep, tbl$gene)
      present <- !is.na(idx)
      if (any(present)) {
        res[[paste0(ct, "_pct.1")]][present] <- tbl$pct.1[idx[present]]
        res[[paste0(ct, "_pct.2")]][present] <- tbl$pct.2[idx[present]]
        res[[paste0(ct, "_avg_log2FC")]][present] <- tbl$avg_log2FC[idx[present]]
        res[[paste0(ct, "_p_val_adj")]][present] <- tbl$p_val_adj[idx[present]]
      }
    }
  }
  res$direction <- direction
  res
}

addSmallLegend <- function(myPlot, pointSize = 3, textSize = 8, spaceLegend = 1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

safe_read_deg <- function(path) {
  # try read_tsv then read_csv; adjust if your files use different separators
  tryCatch(read_tsv(path, show_col_types = FALSE),
           error = function(e) tryCatch(read_csv(path, show_col_types = FALSE),
                                        error = function(e2) read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE)))
}

