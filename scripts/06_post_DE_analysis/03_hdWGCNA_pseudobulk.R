#================================================================================
# hdWGCNA pseudobulk 
#================================================================================
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source("file_paths_and_colours.R")
# Load object
project_ID <- "CWOW_cellbender"
dataObject <- readRDS(paste0("../rObjects/", project_ID, "_annotated_with_celltype_subclusters_pass2.rds"))
DefaultAssay(dataObject) <- "RNA"

# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(12345)

dataObject <- FindVariableFeatures(
  dataObject,
  selection.method = "vst",
  nfeatures = 10000
)
length(VariableFeatures(dataObject))

#-- create pseduobulk counts matrix
dataObject <- SetupForWGCNA(
  dataObject,
  gene_select = "custom", 
  features = VariableFeatures(dataObject),
  wgcna_name = "pseudobulk"
)
print("inspect pseudobulk hdWGCNA")
length(GetWGCNAGenes(dataObject))
head(GetWGCNAGenes(dataObject))
dim(GetDatExpr(dataObject))

# Construt the pseudobulk expression profiles
datExpr <- ConstructPseudobulk(
  dataObject,
  group.by = 'cell_type',
  replicate_col = 'Sample_ID',
  assay = 'RNA',
  slot = 'counts', # this should always be counts!
  min_reps = 6
)

# compute log2CPM normalization
# You can substitute this with another normalization of your choosing.
cpm <- log2(
  t(t(datExpr) / colSums(datExpr)) * 1e6 + 1
)


dataObject <- SetDatExpr(
  dataObject,
  mat = cpm, 
  assay = "RNA",
  wgcna_name = "pseudobulk"
)

# Co-expression network analysis
# Now that we have our pseudobulk matrix, we can perform co-expression network analysis.
# select the soft power threshold
dataObject <- TestSoftPowers(
  dataObject,
  wgcna_name = "pseudobulk"
)

# plot the results:
plot_list <- PlotSoftPowers(dataObject)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
pdf("../results/hdWGCNA/pseudobulk_soft_power.pdf", width = 7, height = 5)
wrap_plots(plot_list, ncol=2)
dev.off()

# construct the co-expression network and identify gene modules
dataObject <- ConstructNetwork(
  dataObject,
  wgcna_name = "pseudobulk",
  tom_name = "pseudobulk",
  overwrite_tom = TRUE,
  mergeCutHeight = 0.15
)

pdf("../results/hdWGCNA/pseudobulk_dendrogram.pdf", width = 7, height = 3.5)
PlotDendrogram(dataObject, main='pseudobulk dendrogram')
dev.off()

# Unique wgcna modules
unique(dataObject@misc$pseudobulk$wgcna_modules$color)

# We can see that this analysis has resulted in 32 co-expression modules across the eight cell types.
# Next we compute the module eigengenes (MEs) and eigengene-based connectivity (kMEs) at the single-cell level, and we can plot the MEs for each module in the different cell types.
# compute the MEs and kMEs
dataObject <- ModuleEigengenes(
  dataObject,
  scale.model.use = "linear",
  wgcna_name = "pseudobulk"
)
# Modules are learned from pseudobulk
# MEs are projected back to single cells

# By default, ModuleConnectivity(dataObject) will calculate connectivity across the entire dataset.
# Instead we will specify the group.by and group_name to ensure kME is calculated within the context of the cell types.
dataObject <- ModuleConnectivity(
  dataObject,
  group.by = 'cell_type',
  group_name = unique(dataObject$cell_type),
  wgcna_name = "pseudobulk"
)

# get MEs from seurat object
MEs <- GetMEs(dataObject)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add MEs to Seurat meta-data for plotting:
meta <- dataObject@meta.data
dataObject@meta.data <- cbind(meta, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(dataObject, features=mods, group.by = 'cell_type')
p
# reset the metadata
dataObject@meta.data <- meta

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient(high='red', low='grey95') + 
  xlab('') + ylab('')

pdf("../results/hdWGCNA/pseudobulk_ME_dotPlot.pdf", width = 5, height = 3.5)
p 
dev.off()
saveRDS(dataObject, file='../rObjects/hdWGCNA_dataObject_pseudobulk.rds', compress = FALSE)

# Up Next:
# Check if modules are cell-type–specific vs shared
# Plot module–trait correlations (e.g., disease, sex, age)
# Extract hub genes per cell type
# plot genes ranked by kME for each module
plot_kMEs <- PlotKMEs(dataObject, ncol=5)
pdf("../results/hdWGCNA/pseudobulk_kME_per_module.pdf", width = 12, height = 3.5)
plot_kMEs
dev.off()

# get the module assignment table:
modules <- GetModules(dataObject) %>% subset(module != 'grey')
write.table(modules, "../results/hdWGCNA/pseudobulk_modules.tsv", sep = "\t", quote = FALSE)
# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(dataObject, n_hubs = 10)
write.table(hub_df, "../results/hdWGCNA/pseudobulk_module_hubs.tsv", sep = "\t", quote = FALSE)
head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
dataObject <- ModuleExprScore(
  dataObject,
  n_genes = 25,
  method='UCell'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  dataObject,
  features='hMEs', # plot the hMEs
  reduction = "integrated.rpca",
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
pdf(paste0("../results/hdWGCNA/pseudobulk_hME_featureplot.pdf"), 10, 7)
wrap_plots(plot_list, ncol=5)
dev.off()

dataObject@reductions
# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  dataObject,
  features='scores', # plot the hub gene scores
  reduction = "umap.rpca",
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

pdf(paste0("../results/hdWGCNA/pseudobulk_hub_gene_score_featureplot.pdf"), 14, 4)
wrap_plots(plot_list, ncol=5)
dev.off()

# radar plot
dataObject$cluster <- do.call(rbind, strsplit(as.character(dataObject$group), ' '))[,1]
dataObject$cluster <- factor(dataObject$cluster, levels = c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS"))
pdf(paste0("../results/hdWGCNA/pseudobulk_radar_all_celltypes.pdf"), 16, 12)
ModuleRadarPlot(
  dataObject,
  group.by = 'cluster',
#  barcodes = dataObject@meta.data %>% subset(individual_clusters == 'microglia') %>% rownames(),
  axis.label.size=4,
  grid.label.size=4
)
dev.off()

# Define the cell types
cell_types <- c(
  "neuron",
  "interneuron",
  "oligodendrocyte",
  "opc",
  "astrocyte",
  "microglia",
  "mural",
  "endothelial", 
  "fibroblast"
)

for (ct in cell_types) {
  print(paste("Generating radar plot for:", ct))
  current_barcodes <- rownames(dataObject@meta.data[dataObject@meta.data$cell_type == ct, ])
  pdf(paste0("../results/hdWGCNA/pseudobulk_radar_", ct, ".pdf"), 16, 12)
  try({
    p <- ModuleRadarPlot(
      dataObject,
      group.by = 'cluster',
      barcodes = current_barcodes,
      axis.label.size = 3.5,
      grid.label.size = 3.5
    )
    print(p)
  })
  dev.off()
}

#------- differential
group1 <- dataObject@meta.data %>% subset(cell_type == 'microglia' & group == "AD_AT") %>% rownames
group2 <- dataObject@meta.data %>% subset(cell_type == 'microglia' & group == "CONTROL") %>% rownames

DMEs <- FindDMEs(
  dataObject,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='pseudobulk'
)

head(DMEs)
DMEs_sub <- subset(DMEs, pct.1 > 0.5 | pct.2 > 0.5)
PlotDMEsVolcano(
  dataObject,
  DMEs_sub,
  wgcna_name = 'pseudobulk'
)


ModuleNetworkPlot(
  dataObject,
  outdir = '../results/hdWGCNA/pseudobulK_ModuleNetworks'
)

dataObject <- RunModuleUMAP(
  dataObject,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)
# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(dataObject, wgcna_name = "pseudobulk")
dataObject
# plot with ggplot
pdf(paste0("../results/hdWGCNA/pseudobulk_module_umap_no_hubs.pdf"), 8, 7)
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
dev.off()

future::plan(sequential)  # keep it sequential if you want
options(future.globals.maxSize = 250 * 1024^3)  # 250 GiB (you have 350GB)

pdf(paste0("../results/hdWGCNA/pseudobulk_module_umap.pdf"), 8, 7)
ModuleUMAPPlot(
  dataObject,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=1 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)
dev.off()

# convert sex to factor
dataObject$Sex <- as.factor(dataObject$Sex)
# convert age_death to numeric
dataObject$Age <- as.numeric(dataObject$Age)
# list of traits to correlate
cur_traits <- c('Thal.amyloid', 'Braak.NFT', 'Cing.LB', 'Age', 'Sex')

dataObject <- ModuleTraitCorrelation(
  dataObject,
  traits = cur_traits,
  group.by='cell_type'
)
# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(dataObject)

names(mt_cor)
names(mt_cor$cor)
mt_cor$cor$neuron

pdf(paste0("../results/hdWGCNA/pseudobulk_module_module_trait_correlation.pdf"), 8, 10)
PlotModuleTraitCorrelation(
  dataObject,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  plot_max = 0.2,
  combine=TRUE
)
dev.off()
