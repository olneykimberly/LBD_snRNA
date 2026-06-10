setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(rstatix)
  library(cowplot)
})

project_ID <- "CWOW_cellbender"
dataObject <- readRDS(paste0("../rObjects/", project_ID, "_filtered_subclusters_pass2.rds"))

# ==============================================================================
# Group labels and colors
# ==============================================================================

pretty_names <- c(
  "CONTROL" = "CONTROL",
  "AD_AT"   = "AD (AT)",
  "LBD_S"   = "LBD (S)",
  "LBD_AS"  = "LBD (AS)",
  "LBD_ATS" = "LBD (ATS)"
)

pretty_levels <- c(
  "CONTROL",
  "AD (AT)",
  "LBD (S)",
  "LBD (AS)",
  "LBD (ATS)"
)

group_cols <- c(
  "CONTROL"   = "#4682B4",
  "AD (AT)"   = "#B4464B",
  "LBD (S)"   = "gray",
  "LBD (AS)"  = "#88778D",
  "LBD (ATS)" = "gray10"
)

recode_group_pretty <- function(x) {
  x <- trimws(as.character(x))
  out <- dplyr::recode(x, !!!pretty_names, .default = NA_character_)
  factor(out, levels = pretty_levels)
}

metadata$Brain_weight <- c(
  1060, 960, 889, 1370, 1020, 1100, 1320, 1020, 900, 980,
  1000, 1180, 1120, 900, 1040, 1480, 1100, 1040, 1180, 1220,
  1120, 1080, 1160, 1140, 1120, 1260, 1220, 1240, 1240, 1140,
  1080, 1140, 1180, 1280, 1380
)

metadata <- metadata %>%
  mutate(
    TYPE_pretty = recode_group_pretty(TYPE),
    sex_inferred = factor(sex_inferred, levels = c("female", "male"))
  )

dataObject$group_pretty <- recode_group_pretty(dataObject$group)

cat("\nMetadata group check:\n")
print(table(metadata$TYPE, metadata$TYPE_pretty, useNA = "ifany"))

cat("\nSeurat object group check:\n")
print(table(dataObject$group, dataObject$group_pretty, useNA = "ifany"))

# ==============================================================================
# Nuclei by disease type
# ==============================================================================

data <- as.data.frame(table(dataObject$group_pretty))
colnames(data) <- c("group", "frequency")
data$group <- factor(data$group, levels = pretty_levels)

ncells_group <- ggplot(data, aes(x = group, y = frequency, fill = group)) + 
  geom_col() +
  theme_classic() +
  geom_text(
    aes(label = frequency), 
    position = position_dodge(width = 0.9),
    vjust = -0.25,
    size = 2.5
  ) +
  scale_fill_manual(values = group_cols, drop = FALSE) + 
  scale_y_continuous(breaks = seq(0, 100000, by = 20000), limits = c(0, 100000)) +
  ggtitle("Nuclei per disease type") +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 7, face = "plain"),
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ncells_group

# ==============================================================================
# Nuclei by sample
# ==============================================================================

data <- as.data.frame(table(dataObject$Sample_ID))
colnames(data) <- c("Sample_ID", "frequency")

ncells_sample <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(
    aes(label = frequency), 
    position = position_dodge(width = 0.9), 
    size = 2.5, 
    hjust = -0.025,
    angle = 90
  ) +
  scale_fill_manual(values = sample_colors) + 
  scale_y_continuous(breaks = seq(0, 16000, by = 2000), limits = c(0, 16000)) +
  ggtitle("Nuclei per sample") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 7),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 7, face = "plain"),
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ncells_sample

# ==============================================================================
# QC plot by group
# ==============================================================================

plots <- lapply(c("nCount_RNA", "nFeature_RNA", "cell.complexity", "percent.mt"), function(feature) {
  VlnPlot(dataObject, features = feature, group.by = "group_pretty", pt.size = 0) + 
    theme(
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 7, face = "plain"),
      legend.position = "none"
    ) +
    scale_fill_manual(values = group_cols, drop = FALSE)
})

v_group <- CombinePlots(plots, ncol = 4)
v_group

# ==============================================================================
# QC plot by sample
# ==============================================================================

plots <- lapply(c("nCount_RNA", "nFeature_RNA", "cell.complexity", "percent.mt"), function(feature) {
  VlnPlot(dataObject, features = feature, group.by = "Sample_ID", pt.size = 0) + 
    theme(
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 7, face = "plain"),
      legend.position = "none"
    ) +
    scale_fill_manual(values = sample_colors)
})

v_sample <- CombinePlots(plots, nrow = 4)
v_sample

row1 <- ggarrange(
  ncells_group,
  ncells_sample,
  ncol = 2,
  widths = c(1, 2.35),
  labels = c("A", "B"), 
  font.label = list(size = 10)
)

row2 <- ggarrange(
  v_group,
  labels = c("C"), 
  font.label = list(size = 10)
)

combined <- ggarrange(
  row1,
  row2,
  nrow = 2
)

combined

path <- paste0("../manuscript_figures/Supplemental_Figure_QC")
saveToPDF(paste0(path, ".pdf"), width = 6.5, height = 4.5)

# ==============================================================================
# Sex check
# ==============================================================================

plots <- lapply(c("XIST", "UTY"), function(feature) {
  VlnPlot(dataObject, features = feature, group.by = "Sample_ID", pt.size = 0) + 
    theme(
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 7, face = "plain"),
      legend.position = "none"
    ) +
    scale_fill_manual(values = sample_colors)
})

v_sex_check <- CombinePlots(plots, nrow = 2)
v_sex_check

path <- paste0("../manuscript_figures/Supplemental_Figure_sex_check")
saveToPDF(paste0(path, ".pdf"), width = 6.5, height = 4.5)

# ==============================================================================
# Clinical/pathology data: Kruskal-Wallis + Dunn post hoc BH
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
})

clinical_out_dir <- "../results/clinical_pathology_stats"
dir.create(clinical_out_dir, recursive = TRUE, showWarnings = FALSE)

pretty_names <- c(
  "CONTROL" = "CONTROL",
  "AD_AT"   = "AD (AT)",
  "LBD_S"   = "LBD (S)",
  "LBD_AS"  = "LBD (AS)",
  "LBD_ATS" = "LBD (ATS)"
)

pretty_levels <- c(
  "CONTROL",
  "AD (AT)",
  "LBD (S)",
  "LBD (AS)",
  "LBD (ATS)"
)

group_cols <- c(
  "CONTROL"   = "#4682B4",
  "AD (AT)"   = "#B4464B",
  "LBD (S)"   = "gray",
  "LBD (AS)"  = "#88778D",
  "LBD (ATS)" = "gray10"
)

recode_group_pretty <- function(x) {
  x <- trimws(as.character(x))
  out <- dplyr::recode(x, !!!pretty_names, .default = NA_character_)
  factor(out, levels = pretty_levels)
}

metadata <- metadata %>%
  dplyr::mutate(
    TYPE_pretty = recode_group_pretty(TYPE),
    sex_inferred = factor(sex_inferred, levels = c("female", "male"))
  )

cat("\nGroup recode check:\n")
print(table(metadata$TYPE, metadata$TYPE_pretty, useNA = "ifany"))

# ----------------------------- Clinical variables ------------------------------

metadata_continuous <- data.frame(
  Brain_weight  = metadata$Brain_weight,
  Age           = metadata$Age,
  Cing_LB       = metadata$Cing.LB,
  Thal_amyloid  = metadata$Thal.amyloid,
  Braak_NFT     = metadata$Braak.NFT
)

column_variables <- c(
  "Brain weight in grams",
  "Age in years",
  "Lewy bodies per mm2",
  "Thal amyloid phase",
  "Braak NFT stage"
)

get_label <- function(j) {
  if (j == "Lewy bodies per mm2") {
    return(expression(paste("Lewy bodies per ", "mm"^2)))
  } else {
    return(j)
  }
}

safe_dunn <- function(df) {
  df <- as.data.frame(df)

  df <- df %>%
    dplyr::filter(!is.na(plot_value), !is.na(TYPE_pretty))

  df$TYPE_pretty <- droplevels(df$TYPE_pretty)

  if (nrow(df) < 2 || dplyr::n_distinct(df$TYPE_pretty) < 2) {
    return(tibble::tibble())
  }

  tryCatch(
    rstatix::dunn_test(df, plot_value ~ TYPE_pretty, p.adjust.method = "BH"),
    error = function(e) tibble::tibble()
  )
}

make_stats_df <- function(i, j) {
  metadata %>%
    dplyr::mutate(
      plot_value = i,
      variable = j
    ) %>%
    dplyr::filter(!is.na(plot_value), !is.na(TYPE_pretty))
}

stats_df <- dplyr::bind_rows(
  Map(make_stats_df, i = metadata_continuous, j = column_variables)
)

dunn_table <- stats_df %>%
  dplyr::group_by(variable) %>%
  dplyr::group_modify(~ safe_dunn(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    comparison = paste(group1, "vs", group2),
    test = "Dunn post hoc BH"
  )

write.csv(
  dunn_table,
  file.path(clinical_out_dir, "clinical_dunn_posthoc_BH_by_group.csv"),
  row.names = FALSE
)

dunn_table_by_sex <- stats_df %>%
  dplyr::filter(!is.na(sex_inferred)) %>%
  dplyr::group_by(variable, sex_inferred) %>%
  dplyr::group_modify(~ safe_dunn(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    comparison = paste(group1, "vs", group2),
    test = "Dunn post hoc BH"
  )

write.csv(
  dunn_table_by_sex,
  file.path(clinical_out_dir, "clinical_dunn_posthoc_BH_by_group_within_sex.csv"),
  row.names = FALSE
)

# ----------------------------- Group plot function -----------------------------

violin_plot_fun_group <- function(i, j) {

  plot_label <- get_label(j)

  plot_df <- metadata %>%
    dplyr::mutate(plot_value = i) %>%
    dplyr::filter(!is.na(plot_value), !is.na(TYPE_pretty))

  p <- ggplot(plot_df, aes(TYPE_pretty, plot_value, fill = TYPE_pretty)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.5, linewidth = 0.25) +
    geom_jitter(
      aes(color = TYPE_pretty),
      shape = 16,
      position = position_jitter(0.1),
      alpha = 0.9,
      size = 1
    ) +
    stat_compare_means(
      method = "kruskal.test",
      aes(label = paste0("Kruskal-Wallis, p = ", after_stat(p.format))),
      label.x.npc = 0.5,
      label.y.npc = 0.98,
      hjust = 0.5,
      size = 5.5 / ggplot2::.pt
    ) +
    theme_bw(base_size = 7) +
    theme(
      text = element_text(size = 7),
      plot.title = element_text(size = 7, face = "plain"),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(plot_label) +
    ylab(plot_label) +
    scale_x_discrete(drop = FALSE, limits = pretty_levels) +
    scale_fill_manual(values = group_cols, drop = FALSE) +
    scale_color_manual(values = group_cols, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20)))

  if (j == "Age in years") {
    p <- p +
      scale_y_continuous(
        limits = c(NA, 110),
        breaks = seq(0, 110, 20),
        expand = expansion(mult = c(0.05, 0.20))
      )
  } else if (j == "Thal amyloid phase") {
    p <- p +
      scale_y_continuous(
        breaks = 0:6,
        limits = c(0, 7),
        expand = expansion(mult = c(0.05, 0.20))
      )
    }
  else if (j == "Braak NFT stage") {
    p <- p +
      scale_y_continuous(
        breaks = 0:6,
        limits = c(0, 7),
        expand = expansion(mult = c(0.05, 0.20))
      )
  }

  return(p)
}

# ----------------------------- Sex facet plot function -------------------------

violin_plot_fun_sex <- function(i, j) {

  plot_label <- get_label(j)

  plot_df <- metadata %>%
    dplyr::mutate(plot_value = i) %>%
    dplyr::filter(!is.na(plot_value), !is.na(TYPE_pretty), !is.na(sex_inferred))

  p <- ggplot(plot_df, aes(TYPE_pretty, plot_value, fill = TYPE_pretty)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.5, linewidth = 0.25) +
    geom_jitter(
      aes(color = TYPE_pretty),
      shape = 16,
      position = position_jitter(0.1),
      alpha = 0.9,
      size = 1
    ) +
    stat_compare_means(
      method = "kruskal.test",
      aes(label = paste0("Kruskal-Wallis, p = ", after_stat(p.format))),
      label.x.npc = 0.5,
      label.y.npc = 0.98,
      hjust = 0.5,
      size = 5.5 / ggplot2::.pt
    ) +
    theme_bw(base_size = 7) +
    theme(
      text = element_text(size = 7),
      plot.title = element_blank(),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      strip.text = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    ylab(plot_label) +
    scale_x_discrete(drop = FALSE, limits = pretty_levels) +
    scale_fill_manual(values = group_cols, drop = FALSE) +
    scale_color_manual(values = group_cols, drop = FALSE) +
    facet_grid(. ~ sex_inferred) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20)))

  if (j == "Age in years") {
    p <- p +
      scale_y_continuous(
        limits = c(NA, 110),
        breaks = seq(0, 110, 20),
        expand = expansion(mult = c(0.05, 0.20))
      )
  } else if (j == "Thal amyloid phase") {
    p <- p +
      scale_y_continuous(
        breaks = 0:6,
        limits = c(0, 7),
        expand = expansion(mult = c(0.05, 0.20))
      )
  }
   else if (j == "Braak NFT stage") {
    p <- p +
      scale_y_continuous(
        breaks = 0:6,
        limits = c(0, 7),
        expand = expansion(mult = c(0.05, 0.20))
      )
  }

  return(p)
}

# ----------------------------- Generate and save -------------------------------

plots_col1 <- Map(
  f = violin_plot_fun_group,
  i = metadata_continuous,
  j = column_variables
)

plots_col2 <- Map(
  f = violin_plot_fun_sex,
  i = metadata_continuous,
  j = column_variables
)

all_plots <- list(
  plots_col1[[1]], plots_col2[[1]],
  plots_col1[[2]], plots_col2[[2]],
  plots_col1[[3]], plots_col2[[3]],
  plots_col1[[4]], plots_col2[[4]],
  plots_col1[[5]], plots_col2[[5]]
)

final_plot <- ggarrange(
  plotlist = all_plots,
  ncol = 2,
  nrow = 5,
  widths = c(1, 1.8),
  labels = c("A", "", "B", "", "C", "", "D", "", "E", ""),
  font.label = list(size = 10),
  common.legend = FALSE
)

path <- "../manuscript_figures/Supplemental_Figure_pathology"

ggsave(
  paste0(path, ".pdf"),
  final_plot,
  width = 6.5,
  height = 8,
  device = cairo_pdf
)
# ==============================================================================
# PCA
# ==============================================================================
# ==============================================================================
# Pseudobulk PCA
# ==============================================================================


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

DefaultAssay(dataObject) <- "RNA"

stopifnot(all(c("Patient_ID", "Sample_ID", "group", "cell_type") %in% colnames(dataObject[[]])))

cat("Cells per group:\n")
print(table(dataObject$group))

# ----------------------------- Labels and colors -------------------------------

pretty_names <- c(
  "CONTROL" = "CONTROL",
  "AD-AT"   = "AD (AT)",
  "LBD-S"   = "LBD (S)",
  "LBD-AS"  = "LBD (AS)",
  "LBD-ATS" = "LBD (ATS)"
)

pretty_levels <- c(
  "CONTROL",
  "AD (AT)",
  "LBD (S)",
  "LBD (AS)",
  "LBD (ATS)"
)

group_cols <- c(
  "CONTROL"   = "#4682B4",
  "AD (AT)"   = "#B4464B",
  "LBD (S)"   = "gray",
  "LBD (AS)"  = "#88778D",
  "LBD (ATS)" = "gray10"
)

recode_group_pretty <- function(x) {
  x <- trimws(as.character(x))
  out <- dplyr::recode(x, !!!pretty_names, .default = NA_character_)
  factor(out, levels = pretty_levels)
}

cell_type_pretty_names <- c(
  "neuron"          = "Excitatory neuron",
  "interneuron"     = "Inhibitory neuron",
  "oligodendrocyte" = "oligodendrocyte",
  "opc"             = "opc",
  "astrocyte"       = "astrocyte",
  "microglia"       = "microglia",
  "endothelial"     = "endothelial",
  "fibroblast"      = "fibroblast",
  "mural"           = "mural"
)

cell_type_levels <- c(
  "Excitatory neuron",
  "Inhibitory neuron",
  "oligodendrocyte",
  "opc",
  "astrocyte",
  "microglia",
  "endothelial",
  "fibroblast",
  "mural"
)

cell_type_cols <- c(
  "Excitatory neuron" = "#E69F00",
  "Inhibitory neuron" = "#56B4E9",
  "oligodendrocyte"   = "#009E73",
  "opc"               = "#F0E442",
  "astrocyte"         = "#0072B2",
  "microglia"         = "#D55E00",
  "endothelial"       = "#CC79A7",
  "fibroblast"        = "gray40",
  "mural"             = "#B8860B"
)

# ----------------------------- Protein-coding genes ----------------------------

genes_annot <- readRDS("../rObjects/annotation.rds")

protein_coding_genes <- genes_annot %>%
  dplyr::filter(gene_type == "protein_coding")

seurat_features <- rownames(dataObject[["RNA"]])

pc_gene_symbols <- unique(protein_coding_genes$gene_name)
pc_gene_symbols <- pc_gene_symbols[!is.na(pc_gene_symbols) & pc_gene_symbols != ""]

features_use <- intersect(pc_gene_symbols, seurat_features)

cat(sprintf(
  "\nProtein-coding overlap with Seurat RNA features: %d / %d (%.1f%%)\n",
  length(features_use),
  length(pc_gene_symbols),
  100 * length(features_use) / length(pc_gene_symbols)
))

if (length(features_use) < 1000) {
  warning(
    "Low overlap between annotation gene_name and Seurat RNA features. ",
    "Your object may not be using gene symbols. Consider using rownames(dataObject[['RNA']]) instead."
  )
}

# ----------------------------- Pseudobulk aggregation --------------------------

dataObject$pb_patient  <- trimws(as.character(dataObject$Patient_ID))
dataObject$pb_group    <- trimws(as.character(dataObject$group))
dataObject$pb_celltype <- trimws(tolower(as.character(dataObject$cell_type)))

dataObject$pb_patient <- gsub("_", "-", dataObject$pb_patient, fixed = TRUE)

grouping_vars <- c("pb_patient", "pb_group", "pb_celltype")

dataObject.pseudo <- AggregateExpression(
  object = dataObject,
  assays = "RNA",
  features = features_use,
  group.by = grouping_vars,
  return.seurat = TRUE,
  slot = "counts"
)

dataObject.pseudo <- NormalizeData(
  dataObject.pseudo,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)

dataObject.pseudo <- FindVariableFeatures(
  dataObject.pseudo,
  selection.method = "vst",
  nfeatures = 2000
)

dataObject.pseudo <- ScaleData(
  dataObject.pseudo,
  features = rownames(dataObject.pseudo)
)

dataObject.pseudo <- RunPCA(
  dataObject.pseudo,
  features = VariableFeatures(dataObject.pseudo)
)

# ----------------------------- Pretty labels -----------------------------------

dataObject.pseudo$cell_type_pretty <- dplyr::recode(
  as.character(dataObject.pseudo$pb_celltype),
  !!!cell_type_pretty_names,
  .default = NA_character_
)

dataObject.pseudo$cell_type_pretty <- factor(
  dataObject.pseudo$cell_type_pretty,
  levels = cell_type_levels
)

dataObject.pseudo$pb_group_pretty <- recode_group_pretty(dataObject.pseudo$pb_group)

cat("\nPseudo-bulk cell type check:\n")
print(table(dataObject.pseudo$pb_celltype, dataObject.pseudo$cell_type_pretty, useNA = "ifany"))

cat("\nPseudo-bulk group check:\n")
print(table(dataObject.pseudo$pb_group, dataObject.pseudo$pb_group_pretty, useNA = "ifany"))

# ----------------------------- Extract PCA coordinates --------------------------

pca_df <- as.data.frame(Embeddings(dataObject.pseudo, reduction = "pca")) %>%
  dplyr::mutate(
    cell_type_pretty = factor(
      as.character(dataObject.pseudo$cell_type_pretty),
      levels = cell_type_levels
    ),
    pb_group_pretty = factor(
      as.character(dataObject.pseudo$pb_group_pretty),
      levels = pretty_levels
    )
  )

cat("\nPCA dataframe cell type check:\n")
print(table(pca_df$cell_type_pretty, useNA = "ifany"))

cat("\nPCA dataframe group check:\n")
print(table(pca_df$pb_group_pretty, useNA = "ifany"))

# ----------------------------- PCA theme ---------------------------------------

pca_theme <- theme_classic(base_size = 7) +
  theme(
    text = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    plot.title = element_text(size = 7, face = "plain", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.22, "cm"),
    legend.spacing.y = unit(0.01, "cm"),
    plot.margin = margin(3, 3, 3, 3)
  )

# ----------------------------- PCA plots ---------------------------------------

PCA_celltype <- ggplot(
  pca_df,
  aes(x = PC_1, y = PC_2, color = cell_type_pretty)
) +
  geom_point(size = 1, alpha = 0.9) +
  scale_color_manual(
    values = cell_type_cols,
    breaks = cell_type_levels,
    limits = cell_type_levels,
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(size = 3, alpha = 1),
      ncol = 1
    )
  ) +
  xlab("PCA 1") +
  ylab("PCA 2") +
  ggtitle("Cell type") +
  pca_theme

PCA_group <- ggplot(
  pca_df,
  aes(x = PC_1, y = PC_2, color = pb_group_pretty)
) +
  geom_point(size = 1, alpha = 0.9) +
  scale_color_manual(
    values = group_cols,
    breaks = pretty_levels,
    limits = pretty_levels,
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(size = 3, alpha = 1),
      ncol = 1
    )
  ) +
  xlab("PCA 1") +
  ylab("PCA 2") +
  ggtitle("Group") +
  pca_theme

PCA_celltype
PCA_group

# ----------------------------- Save figure -------------------------------------

dir.create("../results/pca", showWarnings = FALSE, recursive = TRUE)

combined_pca <- PCA_celltype + PCA_group +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1))

pdf(
  "../manuscript_figures/Supplemental_Figure_pseudobulk_RNA_logNormalize_PCA.pdf",
  width = 6.5,
  height = 3.5
)
print(combined_pca)
dev.off()