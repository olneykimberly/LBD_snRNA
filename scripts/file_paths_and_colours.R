#--- libraries
.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library/4.3", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
library(Seurat)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
#library(Azimuth)
library(dittoSeq)
library(dplyr)
library(RColorBrewer)
library(DESeq2)
require(openxlsx)
library(ggrepel)
library(glmGamPoi)
library(devtools)
library(harmony)
library(DoubletFinder)
library(reshape2)
library(ggtree)
library(DoubletFinder)


#detach("package:xlsx", unload = TRUE)

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
#pathToRawData = c("/research/labs/neurology/fryer/projects/PMI/")
gene_info <- read.delim(paste0(pathToRef, "/star/geneInfo.tab"), header = FALSE)
gene_info = gene_info[-1,]
gene_info <- gene_info %>% 
  rename(
    V1 = "gene_ID",
    V2 = "gene_name", 
    V3 = "type"
  )

# cell cycle 
#cell_cycle_markers <- read.delim("/research/labs/neurology/fryer/projects/references/mouse/cell_cycle_mouse.tsv")
#m.s.genes <- subset(cell_cycle_markers, phase == "S")
#m.g2m.genes <- subset(cell_cycle_markers, phase != "S")

#--- functions 
saveToPDF <- function(...) {
  d = dev.copy(pdf, ...)
  dev.off(d)
}

