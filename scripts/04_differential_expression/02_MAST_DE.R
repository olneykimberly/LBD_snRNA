#--------------------------------------------------------------------------------
# MAST differential expression 
#--------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source("file_paths_and_colours.R")

#--------------------------------------------------------------------------------
# 1) Load object
#--------------------------------------------------------------------------------
project_ID <- "CWOW_cellbender"
dataObject <- readRDS(paste0("../rObjects/", project_ID, "_annotated_with_celltype_subclusters_pass2.rds"))
DefaultAssay(dataObject) <- "RNA"
Layers(dataObject[["RNA"]])

genes <- readRDS("../rObjects/annotation.rds")
protein_coding_genes <- subset(genes, gene_type == "protein_coding")
Idents(dataObject) <- factor(dataObject$cell_type, levels = unique(dataObject$cell_type))
table(Idents(dataObject)) 

cell_types <- levels(factor(dataObject$cell_type))
cell_types
dataObject$celltype_group <- paste(dataObject$cell_type, dataObject$group, sep = "_")
Idents(dataObject) <- "celltype_group"

genes_pc <- protein_coding_genes$gene_name
mt.genes.df <- subset(genes, seqnames == "chrM")
genes_mt <- mt.genes.df$gene_name # 13 Mito genes 
genes_pc_exclude_mt <- genes_pc[!genes_pc %in% genes_mt]

#  only protein coding genes
genes_pc_filtered <- intersect(genes_pc_exclude_mt, rownames(dataObject))

comparisons <- list(
  c("AD_AT", "CONTROL"),
  c("LBD_S", "CONTROL"),
  c("LBD_AS", "CONTROL"),
  c("LBD_ATS", "CONTROL"),
  c("LBD_S", "AD_AT"),
  c("LBD_AS", "AD_AT"),
  c("LBD_ATS", "AD_AT"), 
  c("LBD_AS", "LBD_S"),
  c("LBD_ATS", "LBD_S"),
  c("LBD_ATS", "LBD_AS")
)

main_output_dir <- "../results/DEGs_RNA_pct0.25/"
cell_type <- "microglia"

# Iterate through Cell Types
for (cell_type in cell_types) {
  message(paste("Starting differential expression for cell type:", cell_type))
  
  # Create a subdirectory for the current cell type
  cell_type_dir <- file.path(main_output_dir, cell_type)
  dir.create(cell_type_dir, showWarnings = FALSE)
  
  # Inner Loop: Iterate through Pairwise Comparisons
  for (comp in comparisons) {
    pathology_A <- comp[1]
    pathology_B <- comp[2]
    
    # Construct the full identity labels for FindMarkers
    ident_1_label <- paste0(cell_type, "_", pathology_A)
    ident_2_label <- paste0(cell_type, "_", pathology_B)
    
    message(paste("  Running DE for:", pathology_A, "vs", pathology_B))
    
    de_results <- FindMarkers(
      object = dataObject,
      ident.1 = ident_1_label,
      ident.2 = ident_2_label,
      features = genes_pc_filtered,
      test.use = "MAST",
      min.pct = 0.25, # in either group Default 0.01
      latent.vars = c("sex_inferred", "Age", "Sample_ID"), # Crucial for controlling sex/age variance 
      min.cells.group = 100,
      verbose = FALSE # Set to FALSE in a loop for cleaner output
    )
    
    de_results$gene <- rownames(de_results)
    
    # Save the results
    filename <- paste0("DEG_", cell_type, "_", pathology_A, "_vs_", pathology_B, ".txt")
    write.table(
      de_results, 
      file.path(cell_type_dir, filename), 
      sep = "\t", 
      quote = FALSE, 
      row.names = FALSE
    )
    message(paste("  --> Saved to:", filename))
  }
}
