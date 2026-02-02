#!/usr/bin/env Rscript
# =============================================================================
# Pseudobulk DE per cell type using limma/voom (patient-level replicates)
# =============================================================================

knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source("file_paths_and_colours.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
})

# ----------------------------- Inputs ----------------------------------------
project_ID <- "CWOW_cellbender"
in_rds <- paste0("../rObjects/", project_ID, "_filtered_subclusters_pass2.rds")

genes <- readRDS("../rObjects/annotation.rds")
protein_coding_genes <- subset(genes, gene_type == "protein_coding")
genes_pc <- unique(protein_coding_genes$gene_name)

mt.genes.df <- subset(genes, seqnames == "chrM")
genes_mt <- unique(mt.genes.df$gene_name)

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

main_output_dir <- "../results/DEGs_pseudobulk_limma_voom/"
dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- Load object -----------------------------------
obj <- readRDS(in_rds)
DefaultAssay(obj) <- "RNA"

md <- obj[[]]

# ----------------------------- Identify patient ID column --------------------
candidate_patient_cols <- c("Individual_ID", "patient_id", "donor_id", "Sample_ID", "orig.ident")
patient_col <- candidate_patient_cols[candidate_patient_cols %in% colnames(md)][1]

if (is.na(patient_col) || length(patient_col) == 0) {
  stop("No patient ID column found. Expected one of: ",
       paste(candidate_patient_cols, collapse = ", "),
       "\nAvailable columns include: ", paste(head(colnames(md), 50), collapse = ", "))
}

message("Using patient ID column: ", patient_col)

# ----------------------------- Required metadata checks ----------------------
required_cols <- c("cell_type", "group", "sex_inferred", "Age",
                   "APOE", "Cing.LB", "Thal.amyloid", "Braak.NFT",
                   "prep.batch", "library.pcr.batch", "library.pcr.cycles")

missing_cols <- setdiff(required_cols, colnames(md))
if (length(missing_cols) > 0) {
  stop("Missing required meta.data columns: ", paste(missing_cols, collapse = ", "))
}

# ----------------------------- Gene set for DE -------------------------------
genes_pc_exclude_mt <- setdiff(genes_pc, genes_mt)
genes_use <- intersect(genes_pc_exclude_mt, rownames(obj))
message("Genes used (protein-coding, non-mt, present in object): ", length(genes_use))

# ----------------------------- Helper: safe factor/numeric -------------------
to_factor <- function(x) factor(x)
to_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ----------------------------- Covariates ------------------------------------
# NOTE: With ~8 patients/group, this is a LOT of covariates.
# If your design becomes rank-deficient, reduce covariates (see notes below).
build_design <- function(df) {
  df$group <- factor(df$group)
  df$sex_inferred <- to_factor(df$sex_inferred)
  df[[patient_col]] <- to_factor(df[[patient_col]])
  
  # Try to make numeric where appropriate (keep as factor if non-numeric)
  df$Age <- to_numeric(df$Age)
  if (all(is.na(df$Age))) stop("Age could not be parsed as numeric. Check meta.data$Age.")
  
  # If these are coded as strings like "0/1/2" or "E3/E4", keep as factor by default:
  df$APOE <- to_factor(df$APOE)
  
  # Pathology scores: try numeric, otherwise factor
  for (v in c("Cing.LB","Thal.amyloid","Braak.NFT","library.pcr.cycles")) {
    vv <- to_numeric(df[[v]])
    if (!all(is.na(vv))) df[[v]] <- vv else df[[v]] <- to_factor(df[[v]])
  }
  
  df$prep.batch <- to_factor(df$prep.batch)
  df$library.pcr.batch <- to_factor(df$library.pcr.batch)
  
  # Design without intercept so contrasts are easy: groupA - groupB
  # Keep patient out of design (patient is the replicate unit now)
  design <- model.matrix(
    ~ 0 + group + sex_inferred + Age,
    data = df
  )
  colnames(design) <- make.names(colnames(design))
  list(df = df, design = design)
}

# ----------------------------- Pseudobulk builder ----------------------------
# Updated to include an expression prevalence filter
make_pseudobulk <- function(obj, cell_type_name) {
  md <- obj[[]]
  cells <- rownames(md)[md$cell_type == cell_type_name]
  if (length(cells) == 0) return(NULL)
  
  # 1. Extract counts for this cell type
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[genes_use, cells, drop = FALSE]
  
  # 2. ADDED: Filter genes based on prevalence within this cell type
  # Requirement: Gene must be detected in at least 5% of cells in this specific cell type
  min.pct <- 0.25
  n.cells <- Matrix::rowSums(counts > 0)
  keep_pct <- n.cells >= (length(cells) * min.pct)
  
  # Also ensure a minimum absolute number of cells (e.g., 10 cells)
  keep_min <- n.cells >= 10
  
  genes_to_keep <- names(which(keep_pct & keep_min))
  
  if (length(genes_to_keep) < 10) {
    message("Too few genes pass prevalence filter for ", cell_type_name)
    return(NULL)
  }
  
  counts <- counts[genes_to_keep, , drop = FALSE]
  message("Genes passing prevalence filter (>5% cells): ", length(genes_to_keep))
  
  # 3. Aggregate by patient
  patient_ids <- md[cells, patient_col, drop = TRUE]
  patient_ids <- as.character(patient_ids)
  patients <- sort(unique(patient_ids))
  
  G <- Matrix::sparseMatrix(
    i = seq_along(patient_ids),
    j = match(patient_ids, patients),
    x = 1,
    dims = c(length(patient_ids), length(patients)),
    dimnames = list(cells, patients)
  )
  
  pb_counts <- counts %*% G 
  pb_counts <- as.matrix(pb_counts)
  
  # Metadata aggregation
  md_pat <- md[cells, c(patient_col, required_cols), drop = FALSE]
  md_pat[[patient_col]] <- as.character(md_pat[[patient_col]])
  
  pick_one <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)
    x[1]
  }
  
  df <- do.call(rbind, lapply(patients, function(pid) {
    sub <- md_pat[md_pat[[patient_col]] == pid, , drop = FALSE]
    out <- as.list(sapply(sub[1, , drop = FALSE], identity))
    for (cc in colnames(sub)) out[[cc]] <- pick_one(sub[[cc]])
    as.data.frame(out, stringsAsFactors = FALSE)
  }))
  rownames(df) <- df[[patient_col]]
  df <- df[colnames(pb_counts), , drop = FALSE]
  
  list(counts = pb_counts, meta = df)
}# ----------------------------- DE runner per cell type -----------------------
run_de_for_celltype <- function(cell_type_name) {
  message("\n========== Cell type: ", cell_type_name, " ==========")
  
  pb <- make_pseudobulk(obj, cell_type_name)
  if (is.null(pb)) {
    message("No cells found for ", cell_type_name, "; skipping.")
    return(invisible(NULL))
  }
  
  counts <- pb$counts
  df <- pb$meta
  
  # Require enough patients per group for any meaningful DE
  grp_tab <- table(df$group)
  message("Patients per group (within ", cell_type_name, "):")
  print(grp_tab)
  
  # Build DGEList + filtering
  y <- DGEList(counts = counts)
  y <- calcNormFactors(y)
  
  # Filter using the groups present in this cell type
  keep <- filterByExpr(y, group = df$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  
  message("Genes kept after filterByExpr: ", nrow(y))
  
  # Design matrix (may fail if too many covariates vs #patients)
  bd <- build_design(df)
  df2 <- bd$df
  design <- bd$design
  
  # Guard against rank deficiency (common with many covariates + small N)
  qr_rank <- qr(design)$rank
  if (qr_rank < ncol(design)) {
    stop("Design matrix is rank-deficient for cell type ", cell_type_name,
         ". Reduce covariates or collapse factors with sparse levels.\n",
         "Columns: ", paste(colnames(design), collapse = ", "))
  }
  
  # voom (use quality weights for robustness)
  v <- voomWithQualityWeights(y, design, plot = FALSE)
  
  fit <- lmFit(v, design)
  
  # Output directory for this cell type
  out_dir <- file.path(main_output_dir, cell_type_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Run pairwise contrasts
  for (comp in comparisons) {
    A <- comp[1]
    B <- comp[2]
    
    # Skip if either group absent for this cell type
    if (!(A %in% levels(df2$group)) || !(B %in% levels(df2$group))) next
    
    # Ensure enough patients in each group (recommend >= 5; you can relax)
    if (sum(df2$group == A) < 5 || sum(df2$group == B) < 5) {
      message("Skipping ", A, " vs ", B, " (too few patients in one group for ", cell_type_name, ")")
      next
    }
    
    coefA <- make.names(paste0("group", A))
    coefB <- make.names(paste0("group", B))
    
    if (!(coefA %in% colnames(design)) || !(coefB %in% colnames(design))) next
    
    contr <- makeContrasts(contrasts = paste0(coefA, " - ", coefB), levels = design)
    fit2 <- contrasts.fit(fit, contr)
    fit2 <- eBayes(fit2)
    
    tt <- topTable(fit2, number = Inf, sort.by = "P")
    tt$gene <- rownames(tt)
    
    # Save
    fn <- paste0("DEG_pseudobulk_voom_", cell_type_name, "_", A, "_vs_", B, ".tsv")
    write.table(tt, file.path(out_dir, fn),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message("Saved: ", fn)
  }
  
  invisible(NULL)
}

# ----------------------------- Run all cell types ----------------------------
cell_types <- sort(unique(obj$cell_type))

for (ct in cell_types) {
  run_de_for_celltype(ct)
}

message("\nALL DONE. Results in: ", main_output_dir)
