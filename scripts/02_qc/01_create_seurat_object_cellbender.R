## ----working_directory-------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

## ----echo=FALSE, message=FALSE-----------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
project_ID <- "CWOW_cellbender"

## ----gene_info---------------------------------------------------------------------------------------------------------------------------
# create an gene annotation rds file. This is helpful for downstream processes. 
if (file.exists("../rObjects/annotation.rds")) {
  genes <- readRDS("../rObjects/annotation.rds")
} else {
  gtf.file <- paste0(path_ref, "/genes/genes.gtf")
  genes <- rtracklayer::import(gtf.file)
  genes <- as.data.frame(genes)
  genes <- genes[genes$type == "gene",]
  saveRDS(genes, "../rObjects/annotation.rds")
}

# gene table is an overview of the types of genes in the annotation file 
gene_type_table <- table(genes$gene_type)
write.table(gene_type_table, "gene_type_table.tsv", row.names = F, quote = F, sep = "\t")
write.table(genes, "genes.tsv", row.names = F, quote = F, sep = "\t")

# obtain the MT genes. Will be used for filtering down stream. 
mt.genes.df <- subset(genes, seqnames == "chrM")
genes_mt <- mt.genes.df$gene_name # 13 Mito genes 

rm(mt.genes.df, gene_type_table)

## ----seurat_object-----------------------------------------------------------------------------------------------------------------------
prefix <- "../cellbender/"
suffix <- "/outs/filtered_feature_bc_matrix.h5"
# Create an empty list to store Seurat objects
seurat_list <- list()
if (file.exists(paste0("../rObjects/", project_ID, ".rds"))) {
  dataObject <- readRDS(paste0("../rObjects/", project_ID, ".rds"))
} else {
  for (i in order_samples) {
    print(i)
    sample_name <- basename(i)
    obj <- CreateSeuratObject(Read10X_h5(paste0(prefix,i,suffix)))
    assign(i, obj)
      # Store the object in the list
    assign(i, obj)
  }
  # merge objects
  dataObject <- merge(x = BR_Nuclei_0368,
                 y = c(BR_Nuclei_0381, BR_Nuclei_0376, BR_Nuclei_0373, BR_Nuclei_0384, BR_Nuclei_0389, BR_Nuclei_0400, BR_Nuclei_0392, BR_Nuclei_0377, BR_Nuclei_0382, BR_Nuclei_0369, BR_Nuclei_0390, BR_Nuclei_0385, BR_Nuclei_0374, BR_Nuclei_0401, BR_Nuclei_0393, BR_Nuclei_0380, BR_Nuclei_0372, BR_Nuclei_0391, BR_Nuclei_0388, BR_Nuclei_0396, BR_Nuclei_0407, BR_Nuclei_0397, BR_Nuclei_0406, BR_Nuclei_0379, BR_Nuclei_0371, BR_Nuclei_0383, BR_Nuclei_0387, BR_Nuclei_0395, BR_Nuclei_0399, BR_Nuclei_0405, BR_Nuclei_0386, BR_Nuclei_0370, BR_Nuclei_0398, BR_Nuclei_0378, BR_Nuclei_0375, BR_Nuclei_0394, BR_Nuclei_0402, BR_Nuclei_0403),
                 add.cell.ids = samples, 
                 project = "CWOW")
    saveRDS(dataObject, paste0("../rObjects/", project_ID, ".rds")) # remove BR_Nuclei_0404, BR_Nuclei_0376, BR_Nuclei_0407, BR_Nuclei_0394, BR_Nuclei_0403
} 
# Inspect 
dataObject

rm(BR_Nuclei_0368, BR_Nuclei_0381, BR_Nuclei_0376, BR_Nuclei_0373, BR_Nuclei_0384, BR_Nuclei_0389, BR_Nuclei_0400, BR_Nuclei_0392, BR_Nuclei_0377, BR_Nuclei_0382, BR_Nuclei_0369, BR_Nuclei_0390, BR_Nuclei_0385, BR_Nuclei_0374, BR_Nuclei_0401, BR_Nuclei_0393, BR_Nuclei_0380, BR_Nuclei_0372, BR_Nuclei_0391, BR_Nuclei_0388, BR_Nuclei_0396, BR_Nuclei_0407, BR_Nuclei_0397, BR_Nuclei_0406, BR_Nuclei_0379, BR_Nuclei_0371, BR_Nuclei_0383, BR_Nuclei_0387, BR_Nuclei_0395, BR_Nuclei_0399, BR_Nuclei_0405, BR_Nuclei_0386, BR_Nuclei_0370, BR_Nuclei_0398, BR_Nuclei_0378, BR_Nuclei_0375, BR_Nuclei_0394, BR_Nuclei_0402, BR_Nuclei_0403)

