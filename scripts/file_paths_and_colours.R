#----------------- Libraries
.libPaths(c("/tgen_labs/jfryer/kolney/R/rstudio-4.3.0-4-with_modules.sif/libs", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
.libPaths()
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
#library(philentropy) 
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


#--- references and metadata
metadata <-
  read.delim(
    "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/metadata/metadata_seq_info.txt")

pathToRef = c("/tgen_labs/jfryer/projects/references/human/GRCh38/refdata-gex-GRCh38-2024-A")
gene_info <- read.delim(paste0(pathToRef, "/star/geneInfo.tab"), header = FALSE)
gene_info = gene_info[-1,]
colnames(gene_info) <- c("gene_ID", "gene_name", "type")


# cell cycle 
#cell_cycle_markers <- read.delim("/research/labs/neurology/fryer/projects/references/mouse/cell_cycle_mouse.tsv")
#m.s.genes <- subset(cell_cycle_markers, phase == "S")
#m.g2m.genes <- subset(cell_cycle_markers, phase != "S")

#--- functions 
saveToPDF <- function(...) {
  d = dev.copy(pdf, ...)
  dev.off(d)
}

