## ----working_directory-------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

## ----setup-------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
project_ID <- "CWOW_cellbender"

library(Seurat)
library(future)
# Set up parallel processing for faster loading
plan("multicore", workers = 8)
options(future.globals.maxSize = 250 * 1024^3)  # 250 GB
options(future.rng.onMisuse = "ignore")

## ----gene_info---------------------------------------------------------------------------------
if (file.exists("../rObjects/annotation.rds")) {
  genes <- readRDS("../rObjects/annotation.rds")
} else {
  gtf.file <- paste0(path_ref, "/genes/genes.gtf")
  genes <- rtracklayer::import(gtf.file)
  genes <- as.data.frame(genes)
  genes <- genes[genes$type == "gene", ]
  saveRDS(genes, "../rObjects/annotation.rds")
}

# Extract MT genes for QC
mt.genes.df <- subset(genes, seqnames == "chrM")
genes_mt <- mt.genes.df$gene_name

## ----seurat_object_v5_complete------------------------------------------------------------------
prefix <- "../cellbender/"
suffix <- "_filtered_filtered.h5"
sample_ids <- basename(order_samples)

# Helper: Read one sample and create a Seurat object
read_cb_to_seurat <- function(sample_id, prefix, suffix) {
  h5_path <- file.path(prefix, sample_id, paste0(sample_id, suffix))
  
  if (!file.exists(h5_path)) {
    warning("File not found: ", h5_path)
    return(NULL)
  }
  
  # Read CellBender filtered H5
  mat <- Read10X_h5(h5_path)
  
  # Handle cases where Read10X_h5 returns a list (e.g., Gene Expression, Peaks, etc.)
  if (is.list(mat)) {
    mat <- if ("Gene Expression" %in% names(mat)) mat[["Gene Expression"]] else mat[[1]]
  }
  
  obj <- CreateSeuratObject(
    counts = mat, 
    project = sample_id, # Setting project as sample_id for initial ident
    min.cells = 3, 
    min.features = 200
  )
  
  # Explicitly store Sample_ID in metadata
  obj$Sample_ID <- sample_id
  return(obj)
}

rds_path <- file.path("../rObjects", paste0(project_ID, ".rds"))

if (file.exists(rds_path)) {
  dataObject <- readRDS(rds_path)
} else {
  # Load samples in parallel using the future framework
  message("Loading samples...")
  seurat_list <- lapply(sample_ids, read_cb_to_seurat, prefix = prefix, suffix = suffix)
  
  # Filter out any NULLs if files were missing
  seurat_list <- seurat_list[!sapply(seurat_list, is.null)]
  names(seurat_list) <- sapply(seurat_list, function(x) x$Sample_ID[1])
  
  # Merge all samples
  message("Merging samples...")
  dataObject <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list),
    project = "CWOW_LBD_AD"
  )
  saveRDS(dataObject, rds_path)
}

## ----quick_sanity_checks-----------------------------------------------------------------------
message("Total nuclei in object: ", ncol(dataObject))

# Confirm per-sample distribution
print(table(dataObject$Sample_ID))
