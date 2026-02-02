setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")
project_ID <- "CWOW_cellbender"
dataObject <- readRDS(paste0("../rObjects/", project_ID, "_filtered_subclusters_pass2.rds"))
library(ggpubr)

metadata$Brain_weight <- c(1060,960,889,1370,1020,1100,1320,1020,900,980,1000,1180,1120,900,1040,1480,1100,1040,1180,1220,1120,1080,1160,1140,1120,1260,1220,1240,1240,1140,1080,1140,1180,1280,1380)

# Nuclei by disease type
data <- as.data.frame(table(dataObject$group))
colnames(data) <- c("group","frequency")
ncells_group <- ggplot(data, aes(x = group, y = frequency, fill = group)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9),
            vjust=-0.25,
            size = 2.5
            # angle = 90
  ) +
  scale_fill_manual(values = color_ATS) + 
  scale_y_continuous(breaks = seq(0,100000, by = 20000), limits = c(0,100000)) +
  ggtitle("Nuclei per disease type") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 9, face = "plain"),
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust=1))
ncells_group

# Nuclei by sample 
data <- as.data.frame(table(dataObject$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")
ncells_sample <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            size = 2.5, 
            hjust = -.025,
            angle = 90) +
  scale_fill_manual(values = sample_colors) + 
  scale_y_continuous(breaks = seq(0,16000, by = 2000), limits = c(0,16000)) +
  ggtitle("Nuclei per sample") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 9, face = "plain"),
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust=1))
ncells_sample

# QC plot by group
plots <- lapply(c("nCount_RNA", "nFeature_RNA", "cell.complexity", "percent.mt"), function(feature) {
  VlnPlot(dataObject, features = feature, group.by = 'group', pt.size = 0) + 
    theme(axis.text = element_text(size = 9),
          axis.title = element_text(size = 9),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 9, face = "plain"),
          legend.position = "none") +
    scale_fill_manual(values = color_ATS)
  
})
v_group <- CombinePlots(plots, ncol = 4)
v_group


# QC plot by sample
plots <- lapply(c("nCount_RNA", "nFeature_RNA", "cell.complexity", "percent.mt"), function(feature) {
  VlnPlot(dataObject, features = feature, group.by = 'Sample_ID', pt.size = 0) + 
    theme(axis.text = element_text(size = 9),
          axis.title = element_text(size = 9),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 9, face = "plain"),
          legend.position = "none") +
    scale_fill_manual(values = sample_colors)

})
# Combine the individual plots into a single plot
v_sample <- CombinePlots(plots, nrow = 4)
v_sample

row1 <- ggarrange(
  ncells_group,
  ncells_sample,
  ncol=2,
  widths = c(1,2.35),
  labels = c("A", "B"), 
  font.label = list(size = 10)
)

row2 <- ggarrange(
  v_group,
  #v_sample,
 # nrow=2,
  #heights = c(1,3),
  labels = c("C"), 
  font.label = list(size = 10)
)
  
combined <- ggarrange(
  row1,
  row2,
  nrow=2
  #heights = c(1,3)
)
combined

path <- paste0("../manuscript_figures/Supplemental_Figure_QC")
saveToPDF(paste0(path, ".pdf"), width = 7.5, height = 7)

# Sex check 
plots <- lapply(c("XIST", "UTY"), function(feature) {
  VlnPlot(dataObject, features = feature, group.by = 'Sample_ID', pt.size = 0) + 
    theme(axis.text = element_text(size = 9),
          axis.title = element_text(size = 9),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 9, face = "plain"),
          legend.position = "none") +
    scale_fill_manual(values = sample_colors)
  
})
v_sex_check <- CombinePlots(plots, nrow = 2)
v_sex_check
path <- paste0("../manuscript_figures/Supplemental_Figure_sex_check")
saveToPDF(paste0(path, ".pdf"), width = 7.5, height = 4.5)


#---- Clinical data
metadata_continuous <-
  data.frame(
    metadata$Brain_weight,
    metadata$Age, 
    metadata$Cing.LB, 
    metadata$Thal.amyloid, 
    metadata$Braak.NFT  )

column_variables <-
  c(
    "Brain weight in grams", 
    "Age in years", 
    "Lewy bodies per um2", "
    Thal amyloid phase", 
    "Braak NFT stage"
    
  )

# 1. Label helper for the mu symbol [cite: 23, 30, 61]
get_label <- function(j) {
  if (j == "Lewy bodies per um2") {
    return(expression(paste("Lewy bodies per ", mu, m^2)))
  } else {
    return(j)
  }
}

# 2. Updated Group Plot Function
# ----------------------------------------
violin_plot_fun_group <- function(i, j) {
  plot_label <- get_label(j)
  # Calculate top position based on the specific plot's data
  top_pos <- max(i, na.rm = TRUE) * 1.20 
  
  # Special handling for Age cap [cite: 10, 14]
  if (j == "Age") { top_pos <- 100 }
  
  p <- ggplot(metadata, aes(TYPE, i, fill = TYPE)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.5) + 
    geom_jitter(aes(color = TYPE), shape = 16, 
                position = position_jitter(0.1), alpha = 0.9, size = 1.5) +
    theme_bw() + 
    theme(
      plot.title = element_text(size = 8),
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(plot_label) +
    ylab(plot_label) +
    stat_compare_means(size = 2.5, label.y = top_pos, hjust = -0.25) +
    scale_fill_manual(values = color_ATS) +
    scale_color_manual(values = color_ATS)
  
  # Apply specific Y-axis constraints [cite: 2, 46, 57]
  if (j == "Age") {
    p <- p + scale_y_continuous(limits = c(NA, 105), breaks = seq(0, 100, 10))
  } else if (j == "Thal amyloid phase") {
    p <- p + scale_y_continuous(breaks = 0:6, limits = c(0, 7.5)) 
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
  }
  
  return(p)
}

# 3. Updated Sex Facet Function
# --------------------------------------------
violin_plot_fun_sex <- function(i, j) {
  plot_label <- get_label(j)
  top_pos <- max(i, na.rm = TRUE) * 1.20
  
  if (j == "Age") { top_pos <- 100 }
  
  p <- ggplot(metadata, aes(TYPE, i, fill = TYPE)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.5) + 
    geom_jitter(aes(color = TYPE), shape = 16, 
                position = position_jitter(0.1), alpha = 0.9, size = 1.2) +
    theme_bw() + 
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      strip.text = element_text(size = 8), 
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    ylab(plot_label) +
    stat_compare_means(size = 2.5, label.y = top_pos, hjust = -0.25) +
    scale_fill_manual(values = color_ATS) +
    scale_color_manual(values = color_ATS) +
    facet_grid(. ~ sex_inferred)
  
  if (j == "Age") {
    p <- p + scale_y_continuous(limits = c(NA, 105), breaks = seq(0, 100, 10))
  } else if (j == "Thal amyloid phase") {
    p <- p + scale_y_continuous(breaks = 0:6, limits = c(0, 7.5)) 
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
  }
  
  return(p)
}

# 4. Generate, Interleave, and Save [cite: 13, 22, 39, 45]
# ------------------------------------
plots_col1 <- Map(violin_plot_fun_group, i = metadata_continuous, j = column_variables)
plots_col2 <- Map(violin_plot_fun_sex, i = metadata_continuous, j = column_variables)

all_plots <- list(
  plots_col1[[1]], plots_col2[[1]], 
  plots_col1[[2]], plots_col2[[2]], 
  plots_col1[[3]], plots_col2[[3]], 
  plots_col1[[4]], plots_col2[[4]],
  plots_col1[[5]], plots_col2[[5]]
  
)

final_plot <- ggarrange(
  plotlist = all_plots,
  ncol = 2, nrow = 5,
  widths = c(1, 1.8), 
  labels = c("A", "", "B", "", "C", "", "D", "", "E", ""), 
  font.label = list(size = 10),
  common.legend = FALSE
)
path <- "../manuscript_figures/Supplemental_Figure_pathology"
ggsave(paste0(path, ".pdf"), final_plot, width = 7.5, height = 11, device = cairo_pdf)

#--- PCA
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
