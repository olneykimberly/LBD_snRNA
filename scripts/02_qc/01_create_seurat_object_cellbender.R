## ----working_directory-------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")


## ----echo=FALSE, message=FALSE-----------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellranger"

# Remove sample
metadata <- subset(metadata, Sample_ID != "LBD_AS_F4") # This sample had a low number of nuclei post filtering 
metadata$sampleID <- factor(metadata$Sample_ID, levels = c(metadata$Sample_ID)) # keep sample order 
samples <- metadata$sampleID 
sex_order <- factor(metadata$sex_inferred, levels = unique(metadata$sex_inferred))
disease_order <- factor(metadata$TYPE, levels = c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS")) # The order in which to display the groups

# obtaining the sample ID from the fastq naming
metadata <- metadata %>%
  mutate(sampleID = gsub(".*_(\\d+)_.*_(BR_Nuclei).*", "\\2_\\1", Lane.Name))
samples <- metadata$sampleID 

# Order the sampleID by disease group, so when plotting the samples are ordered by disease group 
# sampleID with disease_order
order <- metadata %>%
  arrange(disease_order) %>%
  dplyr::select(TYPE, sampleID, Sample_ID)
write.table(order, "order.txt", quote = F, row.names = F, sep = "\t") # export the order as a text file
samples <- order$sampleID
disease_order <- order$TYPE
sample_order <- factor(order$Sample_ID, levels = order$Sample_ID)

## ----gene_info---------------------------------------------------------------------------------------------------------------------------
# create an gene annotation rds file. This is helpful for downstream processes. 
if (file.exists("../rObjects/annotation.rds")) {
  genes <- readRDS("../rObjects/annotation.rds")
} else {
  gtf.file <- paste0(pathToRef, "/genes/genes.gtf")
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
mt.genes <- mt.genes.df$gene_name # 13 Mito genes 

## ----seurat_object-----------------------------------------------------------------------------------------------------------------------
prefix <- "../cellranger/"
suffix <- "/outs/filtered_feature_bc_matrix.h5"
# Create an empty list to store Seurat objects
seurat_list <- list()
if (file.exists(paste0("../rObjects/", projectID, ".rds"))) {
  dataObject <- readRDS(paste0("../rObjects/", projectID, ".rds"))
} else {
  for (i in samples) {
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
    saveRDS(dataObject, paste0("../rObjects/", projectID, ".rds"))
} 
# Inspect 
dataObject

rm(BR_Nuclei_0368, BR_Nuclei_0381, BR_Nuclei_0376, BR_Nuclei_0373, BR_Nuclei_0384, BR_Nuclei_0389, BR_Nuclei_0400, BR_Nuclei_0392, BR_Nuclei_0377, BR_Nuclei_0382, BR_Nuclei_0369, BR_Nuclei_0390, BR_Nuclei_0385, BR_Nuclei_0374, BR_Nuclei_0401, BR_Nuclei_0393, BR_Nuclei_0380, BR_Nuclei_0372, BR_Nuclei_0391, BR_Nuclei_0388, BR_Nuclei_0396, BR_Nuclei_0407, BR_Nuclei_0397, BR_Nuclei_0406, BR_Nuclei_0379, BR_Nuclei_0371, BR_Nuclei_0383, BR_Nuclei_0387, BR_Nuclei_0395, BR_Nuclei_0399, BR_Nuclei_0405, BR_Nuclei_0386, BR_Nuclei_0370, BR_Nuclei_0398, BR_Nuclei_0378, BR_Nuclei_0375, BR_Nuclei_0394, BR_Nuclei_0402, BR_Nuclei_0403)

