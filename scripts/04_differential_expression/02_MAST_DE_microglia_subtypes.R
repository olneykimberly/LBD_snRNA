#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
})

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

project_ID <- "CWOW_cellbender"
obj <- readRDS(paste0("../rObjects/subclusters/", project_ID, "_microglia_subclustered.rds"))
DefaultAssay(obj) <- "RNA"

obj$subcluster

p1 <- dittoDimPlot(object = obj,
            size = 0.5,
            var = "subcluster",
            reduction.use = "umap.ct",
            do.label = TRUE,
            labels.highlight = TRUE)


# annotation
obj.annotated <- RenameIdents(object = obj, 
                                     "0" = "homeostatic",
                                     "1" = "activated",
                                     "2" = "activated",
                                     "3" = "activated",
                                     "4" = "cellcyle",
                                     "5" = "MT"
)
obj.annotated$subtype <- factor(Idents(obj.annotated))

p2 <- dittoDimPlot(object = obj.annotated,
             size = 0.5,
             var = "subtype",
             reduction.use = "umap.ct",
             do.label = TRUE,
             labels.highlight = TRUE)

pdf_name <- paste0("../results/DEG_MAST_RNA_pct0.25_celltype_subtypes/microglia/UMAP_clusters_and_subtype_annotation.pdf")
pdf(pdf_name, width = 12, height = 6)
print(p1 + p2)
dev.off()

p3 <- dittoDimPlot(object = obj.annotated,
                   size = 0.5,
                   var = "subtype",
                   reduction.use = "umap.rpca",
                   do.label = TRUE,
                   labels.highlight = TRUE)

pdf_name <- paste0("../results/DEG_MAST_RNA_pct0.25_celltype_subtypes/microglia/UMAP.rpca_clusters_and_subtype_annotation.pdf")
pdf(pdf_name, width = 6, height = 6)
print(p3)
dev.off()


obj.annotated$subtype
count_per_cluster <- FetchData(obj.annotated,
                               vars = c("ident", "group")) %>%
  dplyr::count(ident, group) %>%
  tidyr::spread(ident, n)
count_per_cluster

count_melt <- reshape2::melt(count_per_cluster)
colnames(count_melt) <- c("group", "cluster", "number of nuclei")
count_max <- count_melt[which.max(count_melt$`number of nuclei`), ]
count_max_value <- count_max$`number of nuclei`
cellmax <- count_max_value + 250 # so that the figure doesn't cut off the text

# Re-order the melted dataframe factor
# Re-order the melted dataframe factor if needed
count_bar <- ggplot(count_melt, aes(x = cluster, y = `number of nuclei`, fill = cluster)) +
  geom_bar(
    stat = "identity",
    colour = "black",
    width = 0.8
  ) +
  geom_text(
    aes(label = `number of nuclei`),
    vjust = -0.5,
    angle = 45,
    hjust = -0.1,
    size = 3
  ) +
  # This splits the plot into multiple panels based on your group (e.g., Control/Experimental)
  facet_wrap(~group) + 
  theme_classic() + 
  scale_fill_manual(values = color_panel) + 
  # Use expand instead of a hard 'cellmax' if you want it to be dynamic per facet
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  labs(title = "Number of nuclei per subtype", x = "subtype") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none" # Hide legend since X-axis labels already show clusters
  )
count_bar
pdf_name <- paste0("../results/DEG_MAST_RNA_pct0.25_celltype_subtypes/microglia/subtype_count_per_group.pdf")
pdf(pdf_name, width = 6, height = 6)
print(count_bar)
dev.off()

# Ensure DE uses log-normalized data layer (Seurat v5)
# If you have multiple layers, make sure "data" exists and is current
if (!"data" %in% Layers(obj.annotated[["RNA"]])) {
  stop("RNA assay has no 'data' layer. NormalizeData/JoinLayers may not have been run as expected.")
}
DefaultLayer(obj.annotated[["RNA"]]) <- "data"

genes <- readRDS("../rObjects/annotation.rds")
protein_coding <- subset(genes, gene_type == "protein_coding")$gene_name
mt_genes <- subset(genes, seqnames == "chrM")$gene_name
features_use <- intersect(setdiff(protein_coding, mt_genes), rownames(obj))

# identities: celltype_group
obj.annotated$subtype_group <- paste(obj.annotated$subtype, obj.annotated$group, sep = "_")
Idents(obj.annotated) <- "subtype_group"

subtype_group <- sort(unique(obj.annotated$subtype_group))
table(obj.annotated$subtype_group)

comparisons <- list(
  c("AD_AT", "CONTROL"),
  c("LBD_S", "CONTROL"),
  c("LBD_AS", "CONTROL"),
  c("LBD_ATS", "CONTROL")
)

out_base <- "../results/DEG_MAST_RNA_pct0.25_celltype_subtypes/microglia/"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Reproducibility
set.seed(12345)
options(stringsAsFactors = FALSE)

latent <- c("sex_inferred", "Age")
subtypes <- sort(unique(obj.annotated$subtype))
for (ct in subtypes) {
  message("\n==== subtype: ", ct, " ====")
  
  # Subset ONCE per cell type (major speedup)
  cells_ct <- rownames(obj.annotated[[]])[obj.annotated$subtype == ct]
  if (length(cells_ct) < 50) {
    message("Skipping ", ct, " (too few cells: ", length(cells_ct), ")")
    next
  }
  obj.annotated_ct <- subset(obj.annotated, cells = cells_ct)
  
  # ensure Idents available in subset
  Idents(obj.annotated_ct) <- "subtype_group"
  
  # output dir
  out_dir <- file.path(out_base, ct)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # run comparisons
  for (comp in comparisons) {
    A <- comp[1]; B <- comp[2]
    ident1 <- paste0(ct, "_", A)
    ident2 <- paste0(ct, "_", B)
    
    n1 <- sum(Idents(obj.annotated_ct) == ident1)
    n2 <- sum(Idents(obj.annotated_ct) == ident2)
    if (n1 < 100 || n2 < 100) {
      message("Skipping ", ident1, " vs ", ident2, " (cells: ", n1, " vs ", n2, ")")
      next
    }
    
    message("  DE: ", A, " vs ", B, " (", n1, " vs ", n2, " cells)")
    
    res <- FindMarkers(
      object = obj.annotated_ct,
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
  
  rm(obj.annotated_ct); gc()
}

message("\nDONE. Outputs in: ", out_base)
