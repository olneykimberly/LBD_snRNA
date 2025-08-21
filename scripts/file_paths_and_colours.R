#----------------- Libraries
set.seed(28)
.libPaths(c("/tgen_labs/jfryer/kolney/R/rstudio-4.3.0-4-with_modules.sif/libs", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
.libPaths()
#unloadNamespace("RSpectra")
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

#options(future.globals.maxSize = 1e9)
options(future.globals.maxSize = 150 * 1024^3)


#--- variables
# paths, colors, shapes and more
LBD <- "LBD"
AD <- "AD"
PA <- "PA"
CONTROL <- "CONTROL"
control_color <- "#4682B4" 
AD_color <- "#B4464B" 
LBD_color <- "gray35" 
control_shape <- c(15) # square
AD_shape <- c(16) # circle
PA_shape <- c(17) # triangle
LBD_shape <- c(18) # diamond

TypeColors <- c("#4682B4","#B4464B", "gray35")
ATSColors <- c("#4682B4","#B4464B", "gray35", "gray65", "gray", "gray85")

SexColors <- c("#490092", "#D55E00")
colorbindColors <- dittoColors()
correlationColors <-
  colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

color.panel <- dittoColors()

#--- references and metadata
metadata <-
  read.delim(
    "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/metadata/metadata_seq_info.txt")

pathToRef = c("/tgen_labs/jfryer/projects/references/human/GRCh38/refdata-gex-GRCh38-2024-A")
gene_info <- read.delim(paste0(pathToRef, "/star/geneInfo.tab"), header = FALSE)
gene_info = gene_info[-1,]
colnames(gene_info) <- c("gene_ID", "gene_name", "type")


# cell cycle 
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#--- functions 
saveToPDF <- function(...) {
  d = dev.copy(pdf, ...)
  dev.off(d)
}

fromList <- function (input) {
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

markers.to.plot <-
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

