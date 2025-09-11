#----------------- Libraries
set.seed(28)
.libPaths(c("/tgen_labs/jfryer/kolney/R/rstudio-4.3.0-4-with_modules.sif/libs", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
.libPaths()
#unloadNamespace("RSpectra")
unloadNamespace("rtracklayer")
unloadNamespace("GenomicAlignments")
unloadNamespace("SummarizedExperiment")
unloadNamespace("DelayedArray")
unloadNamespace("SparseArray")
unloadNamespace("S4Arrays")
library(Matrix, lib.loc = "/usr/local/lib/R/site-library")
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
require(openxlsx)
library(ggrepel)
library(glmGamPoi)
library(devtools)
library(harmony)
library(DoubletFinder)
library(reshape2)
library(ggtree)
library(BiocParallel) 
library(edgeR)  
library(limma)  
library(ggrepel) 
library(ggplot2) 
library(gplots) 
library(grDevices)  
library(stringr) 
library(remaCor)
library(scales)
library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(dplyr)
library(plyr)
library(gridExtra)
library(grid)
library(lattice)
library(data.table)
library(openxlsx)
library(readxl)
library(pheatmap)
library(NatParksPalettes)
library(UpSetR)
library(cowplot)
library(ggpubr)
library(patchwork) # For combining plots
library(SeuratData)
library(batchelor)
library(DropletQC)
library(Rsamtools)
library(GenomicRanges)

#options(future.globals.maxSize = 1e9)
options(future.globals.maxSize = 150 * 1024^3)


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
color_ATS <- c("#4682B4","#B4464B", "gray35", "gray65", "gray", "gray85")
color_sex <- c("#490092", "#D55E00")
color_blind <- dittoColors()
color_panel <- dittoColors()

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
