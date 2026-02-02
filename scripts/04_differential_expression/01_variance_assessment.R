#================================================================================
# Pseudobulk (AggregateExpression) + n_nuclei + variancePartition + PCA QC
# Robust to: dash/underscore IDs, joins that drop rownames, and Seurat meta.data
#================================================================================

knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source("file_paths_and_colours.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(edgeR)
  library(variancePartition)
  library(corrplot)
  library(cowplot)
})

#--------------------------------------------------------------------------------
# 1) Load object
#--------------------------------------------------------------------------------
project_ID <- "CWOW_cellbender"
dataObject <- readRDS(paste0("../rObjects/", project_ID, "_annotated_with_celltype_subclusters_pass2.rds"))
DefaultAssay(dataObject) <- "RNA"

stopifnot(all(c("Patient_ID", "Sample_ID", "group", "cell_type") %in% colnames(dataObject[[]])))

cat("Cells per group:\n")
print(table(dataObject$group))

#--------------------------------------------------------------------------------
# 2) Protein-coding features (safe overlap with Seurat features)
#--------------------------------------------------------------------------------
genes_annot <- readRDS("../rObjects/annotation.rds")
protein_coding_genes <- genes_annot %>% dplyr::filter(gene_type == "protein_coding")

seurat_features <- rownames(dataObject[["RNA"]])
pc_gene_symbols <- unique(protein_coding_genes$gene_name)
pc_gene_symbols <- pc_gene_symbols[!is.na(pc_gene_symbols) & pc_gene_symbols != ""]

features_use <- intersect(pc_gene_symbols, seurat_features)

cat(sprintf("\nProtein-coding overlap with Seurat RNA features: %d / %d (%.1f%%)\n",
            length(features_use), length(pc_gene_symbols),
            100 * length(features_use) / length(pc_gene_symbols)))

if (length(features_use) < 1000) {
  warning("Low overlap between annotation gene_name and Seurat RNA features. ",
          "Your object may not be using gene symbols. Consider using rownames(dataObject[['RNA']]) instead.")
}

#--------------------------------------------------------------------------------
# 3) Create harmonized keys (prevents subtle join mismatches)
#    We will aggregate using these keys, so counts + metadata + n_nuclei all match.
#--------------------------------------------------------------------------------
dataObject$pb_patient  <- trimws(as.character(dataObject$Patient_ID))
dataObject$pb_group    <- trimws(as.character(dataObject$group))
dataObject$pb_celltype <- trimws(tolower(as.character(dataObject$cell_type)))

# normalize punctuation just in case (optional but safe)
dataObject$pb_patient <- gsub("_", "-", dataObject$pb_patient, fixed = TRUE)
dataObject$pb_group   <- gsub("_", "-", dataObject$pb_group, fixed = TRUE)

grouping_vars <- c("pb_patient", "pb_group", "pb_celltype")

#--------------------------------------------------------------------------------
# 4) AggregateExpression (pseudo-bulk) using harmonized keys
#--------------------------------------------------------------------------------
dataObject.pseudo <- AggregateExpression(
  object = dataObject,
  assays = "RNA",
  features = features_use,
  group.by = grouping_vars,
  return.seurat = TRUE,
  slot = "counts"
)
DefaultAssay(dataObject.pseudo) <- "RNA"

#--------------------------------------------------------------------------------
# 5) Ensure meta.data rownames EXACTLY match Cells() (critical for subset(), joins)
#    Joins often drop rownames; Seurat requires rownames(meta.data) == Cells(object)
#--------------------------------------------------------------------------------
pseudo_cells <- Cells(dataObject.pseudo)
md <- dataObject.pseudo@meta.data

if (nrow(md) != length(pseudo_cells)) {
  stop("meta.data rows (", nrow(md), ") != number of pseudo cells (", length(pseudo_cells),
       "). Something is inconsistent right after AggregateExpression.")
}

# Force alignment by order
md <- as.data.frame(md)
rownames(md) <- pseudo_cells
dataObject.pseudo@meta.data <- md

stopifnot(identical(rownames(dataObject.pseudo@meta.data), Cells(dataObject.pseudo)))

#--------------------------------------------------------------------------------
# 6) Attach n_nuclei using the SAME harmonized keys
#--------------------------------------------------------------------------------
cell_counts <- dataObject[[]] %>%
  dplyr::count(pb_patient, pb_group, pb_celltype, name = "n_nuclei")

# Do join in a way that preserves cell order and then restore rownames
md <- dataObject.pseudo@meta.data %>% tibble::as_tibble(rownames = "pseudo_cell")
md2 <- md %>%
  dplyr::left_join(cell_counts, by = grouping_vars)

# Restore as data.frame with correct rownames
md2 <- as.data.frame(md2)
rownames(md2) <- md2$pseudo_cell
md2$pseudo_cell <- NULL

# Reorder to Cells() just to be safe
md2 <- md2[Cells(dataObject.pseudo), , drop = FALSE]

if (!"n_nuclei" %in% colnames(md2)) stop("n_nuclei column missing after join (unexpected).")
md2$n_nuclei <- suppressWarnings(as.numeric(md2$n_nuclei))

if (any(is.na(md2$n_nuclei))) {
  bad <- md2 %>%
    dplyr::as_tibble(rownames = "cell") %>%
    dplyr::filter(is.na(n_nuclei)) %>%
    dplyr::select(dplyr::all_of(grouping_vars)) %>%
    dplyr::distinct()
  stop("Failed to attach n_nuclei for some pseudo-bulk profiles. Examples:\n",
       paste(capture.output(print(head(bad, 20))), collapse = "\n"))
}

dataObject.pseudo@meta.data <- md2
stopifnot(identical(rownames(dataObject.pseudo@meta.data), Cells(dataObject.pseudo)))

# Optional: make friendly column names for downstream plotting/modeling
dataObject.pseudo@meta.data <- dataObject.pseudo@meta.data %>%
  dplyr::rename(
    Patient_ID = pb_patient,
    group = pb_group,
    cell_type = pb_celltype
  )

#--------------------------------------------------------------------------------
# 7) Filter pseudo-bulk profiles by nuclei count (robust; will not pass bad cell IDs)
#--------------------------------------------------------------------------------
min_nuclei <- 20
md <- dataObject.pseudo[[]]

cat("\nPseudo-bulk profiles total:", ncol(dataObject.pseudo), "\n")
cat("Summary n_nuclei:\n")
print(summary(md$n_nuclei))

keep_cells <- Cells(dataObject.pseudo)[md$n_nuclei >= min_nuclei]

cat("Pseudo-bulk profiles kept :", length(keep_cells), "\n")
cat("Pseudo-bulk profiles dropped:", ncol(dataObject.pseudo) - length(keep_cells), "\n")

if (length(keep_cells) == 0) {
  md_small <- md %>% dplyr::arrange(n_nuclei) %>% dplyr::select(n_nuclei, Patient_ID, group, cell_type)
  print(head(md_small, 30))
  stop("No pseudo-bulk profiles meet min_nuclei. Lower min_nuclei or fix n_nuclei.")
}

dataObject.pseudo <- subset(dataObject.pseudo, cells = keep_cells)

#--------------------------------------------------------------------------------
# 8) Gene filtering (recommended)
#--------------------------------------------------------------------------------
raw_counts <- GetAssayData(dataObject.pseudo, assay = "RNA", slot = "counts")

# Keep genes with >=10 counts in >=3 pseudo samples (tune if needed)
keep_genes <- rowSums(raw_counts >= 10) >= 3
raw_counts <- raw_counts[keep_genes, , drop = FALSE]

# safer in Seurat v5
dataObject.pseudo <- SetAssayData(
  object = dataObject.pseudo,
  assay = "RNA",
  layer = "counts",
  new.data = raw_counts
)

DefaultAssay(dataObject.pseudo) <- "RNA"

message("\nPost-filter dimensions: genes = ", nrow(raw_counts),
        ", pseudo-samples = ", ncol(raw_counts))

#--------------------------------------------------------------------------------
# 9) Merge donor-level covariates (keyed by Patient_ID, which equals your subject)
#--------------------------------------------------------------------------------
donor_metadata <- dataObject[[]] %>%
  dplyr::distinct(pb_patient, .keep_all = TRUE) %>%
  dplyr::select(
    pb_patient,
    Sex, Age, APOE, Cing.LB, Thal.amyloid, Braak.NFT,
    prep.batch, library.pcr.cycles, library.pcr.batch
  )

# Preserve rownames while joining
cell_order <- Cells(dataObject.pseudo)

md <- dataObject.pseudo@meta.data %>% tibble::as_tibble(rownames = "pseudo_cell")
md <- md %>%
  dplyr::left_join(donor_metadata, by = c("Patient_ID" = "pb_patient"))

md <- as.data.frame(md)
rownames(md) <- md$pseudo_cell
md$pseudo_cell <- NULL
md <- md[cell_order, , drop = FALSE]

dataObject.pseudo@meta.data <- md
stopifnot(identical(rownames(dataObject.pseudo@meta.data), Cells(dataObject.pseudo)))

#--------------------------------------------------------------------------------
# 10) variancePartition inputs (log2 CPM) + formatted metadata
#--------------------------------------------------------------------------------
raw_counts <- GetAssayData(dataObject.pseudo, assay = "RNA", slot = "counts")
expr_mat <- log2(edgeR::cpm(raw_counts) + 1)

varmetadata <- dataObject.pseudo[[]] %>%
  dplyr::mutate(
    disease = as.factor(group),
    sex = as.factor(Sex),
    celltype = as.factor(cell_type),
    Age = suppressWarnings(as.numeric(Age)),
    APOE = as.factor(APOE),
    Cing.LB = suppressWarnings(as.numeric(Cing.LB)),
    Thal.amyloid = suppressWarnings(as.numeric(Thal.amyloid)),
    Braak.NFT = suppressWarnings(as.numeric(Braak.NFT)),
    prep.batch = as.factor(prep.batch),
    library.pcr.batch = as.factor(library.pcr.batch),
    library.pcr.cycles = suppressWarnings(as.numeric(library.pcr.cycles))
  )

stopifnot(identical(colnames(expr_mat), rownames(varmetadata)))
message("✅ Expression matrix and metadata aligned for variancePartition.")

# Model note: random effects for 2-level factors (e.g., disease/sex) are unstable as variance components.
# Here: celltype is random; disease/sex are fixed for interpretability.
vp_formula <- ~ (1|celltype) + (1|disease) + (1|sex) + Age + (1|APOE) + Cing.LB + Thal.amyloid + Braak.NFT 

var_part <- fitExtractVarPartModel(expr_mat, vp_formula, varmetadata)
vp_sorted <- sortCols(var_part)
vp_plot <- plotVarPart(vp_sorted)
print(vp_plot)

#dir.create("../manuscript_figures/variancePartition", showWarnings = FALSE, recursive = TRUE)
pdf(paste0("../manuscript_figures/Supplemental_Figure_variancePartition.pdf"), height = 5, width = 7.5)
print(vp_plot)
dev.off()

cca_form <- ~ disease + sex + celltype + Age + APOE + Cing.LB + Thal.amyloid + Braak.NFT 
C_cov <- canCorPairs(cca_form, varmetadata)
plotCorrMatrix(C_cov)
saveToPDF(paste0("../manuscript_figures/Supplemental_Figure_CCA_celltype.pdf"), width = 6.5, height = 5)

#--------------------------------------------------------------------------------
# 11) PCA QC on pseudo-bulk (visualization only; not for DESeq2)
#--------------------------------------------------------------------------------
dataObject.pseudo <- NormalizeData(dataObject.pseudo, normalization.method = "LogNormalize", scale.factor = 1e4)
dataObject.pseudo <- FindVariableFeatures(dataObject.pseudo, selection.method = "vst", nfeatures = 2000)
dataObject.pseudo <- ScaleData(dataObject.pseudo, features = rownames(dataObject.pseudo))
dataObject.pseudo <- RunPCA(dataObject.pseudo, features = VariableFeatures(dataObject.pseudo))

dataObject.pseudo$cell_type <- factor(dataObject.pseudo$pb_celltype, levels = order_cell_type) 
dataObject.pseudo$pb_group <- factor(dataObject.pseudo$pb_group, levels = c("CONTROL", "AD-AT", "LBD-S", "LBD-AS", "LBD-ATS")) 
p1 <- DimPlot(dataObject.pseudo, reduction = "pca", group.by = "pb_celltype", cols= color_panel)
p2 <- DimPlot(dataObject.pseudo, reduction = "pca", group.by = "pb_group", cols= color_ATS)
print(p1 + p2)

PCA_celltype <- DimPlot(dataObject.pseudo, reduction = "pca", group.by = "pb_celltype", label = FALSE, label.size = 3) + 
  xlab("PCA 1") + 
  ylab("PCA 2") +
  theme(
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 8),
    plot.margin = margin(0.5, 0.025, 0.25, 0.025, "cm")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values = color_panel) +
  ggtitle("Cell type")
PCA_celltype <- addSmallLegend(PCA_celltype)
PCA_celltype


PCA_sample <- DimPlot(dataObject.pseudo, reduction = "pca", group.by = "pb_group", label = FALSE, label.size = 3) + 
  xlab("PCA 1") + 
  ylab("PCA 2") +
  theme(
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 8),
    plot.margin = margin(0.5, 0.025, 0.25, 0.025, "cm")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values = color_ATS) +
  ggtitle("Group")
PCA_sample <- addSmallLegend(PCA_sample)
PCA_sample

dir.create("../results/pca", showWarnings = FALSE, recursive = TRUE)
pdf(paste0("../manuscript_figures/Supplemental_Figure_pseudobulk_RNA_logNormalize_PCA.pdf"), width = 7.5, height = 3.5)
print(PCA_celltype + PCA_sample)
dev.off()
getwd()
