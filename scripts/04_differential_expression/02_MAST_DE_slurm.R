#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
})

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

project_ID <- "CWOW_cellbender"
obj <- readRDS(paste0("../rObjects/", project_ID, "_filtered_subclusters_pass2.rds"))
DefaultAssay(obj) <- "RNA"

# Ensure DE uses log-normalized data layer (Seurat v5)
# If you have multiple layers, make sure "data" exists and is current
if (!"data" %in% Layers(obj[["RNA"]])) {
  stop("RNA assay has no 'data' layer. NormalizeData/JoinLayers may not have been run as expected.")
}
DefaultLayer(obj[["RNA"]]) <- "data"

genes <- readRDS("../rObjects/annotation.rds")
protein_coding <- subset(genes, gene_type == "protein_coding")$gene_name
mt_genes <- subset(genes, seqnames == "chrM")$gene_name
features_use <- intersect(setdiff(protein_coding, mt_genes), rownames(obj))

# identities: celltype_group
obj$celltype_group <- paste(obj$cell_type, obj$group, sep = "_")
Idents(obj) <- "celltype_group"

cell_types <- sort(unique(obj$cell_type))

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

out_base <- "../results/DEGs_RNA_pct0.25/"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Reproducibility
set.seed(12345)
options(stringsAsFactors = FALSE)

latent <- c("sex_inferred", "Age")

for (ct in cell_types) {
  message("\n==== Cell type: ", ct, " ====")
  
  # Subset ONCE per cell type (major speedup)
  cells_ct <- rownames(obj[[]])[obj$cell_type == ct]
  if (length(cells_ct) < 200) {
    message("Skipping ", ct, " (too few cells: ", length(cells_ct), ")")
    next
  }
  obj_ct <- subset(obj, cells = cells_ct)
  
  # ensure Idents available in subset
  Idents(obj_ct) <- "celltype_group"
  
  # output dir
  out_dir <- file.path(out_base, ct)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # run comparisons
  for (comp in comparisons) {
    A <- comp[1]; B <- comp[2]
    ident1 <- paste0(ct, "_", A)
    ident2 <- paste0(ct, "_", B)
    
    n1 <- sum(Idents(obj_ct) == ident1)
    n2 <- sum(Idents(obj_ct) == ident2)
    if (n1 < 100 || n2 < 100) {
      message("Skipping ", ident1, " vs ", ident2, " (cells: ", n1, " vs ", n2, ")")
      next
    }
    
    message("  DE: ", A, " vs ", B, " (", n1, " vs ", n2, " cells)")
    
    res <- FindMarkers(
      object = obj_ct,
      ident.1 = ident1,
      ident.2 = ident2,
      features = features_use,
      test.use = "MAST",
      min.pct = 0.25,
      latent.vars = latent,
      verbose = FALSE
    )
    
    res$gene <- rownames(res)
    fn <- paste0("DEG_", ct, "_", A, "_vs_", B, ".tsv")
    write.table(res, file.path(out_dir, fn),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  rm(obj_ct); gc()
}

message("\nDONE. Outputs in: ", out_base)
