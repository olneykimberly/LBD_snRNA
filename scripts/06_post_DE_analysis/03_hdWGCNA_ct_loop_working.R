#================================================================================
# hdWGCNA per cell type (metacells built per cell type)
#================================================================================
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source("file_paths_and_colours.R")

library(Seurat)
library(hdWGCNA)
library(qs)
library(cowplot)
library(patchwork)
library(dplyr)

theme_set(theme_cowplot())
set.seed(12345)
allowWGCNAThreads(nThreads = 8)

#--------------------------------------------------------------------------------
# Load object
#--------------------------------------------------------------------------------
projectID <- "CWOW_cellbender"
dataObject <- readRDS(paste0("../rObjects/", projectID, "_filtered_subclusters_pass2.rds"))

#--------------------------------------------------------------------------------
# Global WGCNA genes ONCE (on the full object)
#--------------------------------------------------------------------------------
DefaultAssay(dataObject) <- "RNA"

dataObject <- FindVariableFeatures(
  dataObject,
  selection.method = "vst",
  nfeatures = 8000
)

wgcna_genes <- VariableFeatures(dataObject)
softpower_genes <- sample(wgcna_genes, 3000)
message("Total WGCNA genes: ", length(wgcna_genes))

#--------------------------------------------------------------------------------
# Cell types to run
#--------------------------------------------------------------------------------
cell_types <- c(
   "interneuron",
   "oligodendrocyte",
   "opc",
   "fibroblast",
   "endothelial",
   "mural",
   "microglia",
   "astrocyte",
   "neuron"
)

#--------------------------------------------------------------------------------
# Parameters
#--------------------------------------------------------------------------------
min_cells_ct <- 20       # donor-level minimum cells within this cell type
k_floor <- 5
k_cap <- 25
reduction_use <- "integrated.rpca"

out_dir <- "../results/hdWGCNA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj_dir <- "../rObjects"
dir.create(obj_dir, showWarnings = FALSE, recursive = TRUE)

#================================================================================
# LOOP OVER CELL TYPES
#================================================================================
for (ct in cell_types) {
  
  message("======================================")
  message("Running metacells + hdWGCNA for cell type: ", ct)
  message("======================================")
  
  #----------------------------------------------------------------------------
  # 1) Subset to cell type
  #----------------------------------------------------------------------------
  dataObject.ct <- subset(dataObject, subset = cell_type == ct)
  
  # sanity: skip if no cells
  if (ncol(dataObject.ct) == 0) {
    warning("No cells found for cell type: ", ct, " -- skipping.")
    next
  }
  
  #----------------------------------------------------------------------------
  # 2) Filter donors with too few cells in this cell type
  #----------------------------------------------------------------------------
  meta_ct <- dataObject.ct[[]]
  donor_n <- dplyr::count(meta_ct, Sample_ID, name = "n")
  
  keep_donors <- donor_n %>% dplyr::filter(n >= min_cells_ct) %>% dplyr::pull(Sample_ID)
  dataObject.ct <- subset(dataObject.ct, subset = Sample_ID %in% keep_donors)
  
  # if too few donors remain, skip
  donor_n2 <- dplyr::count(dataObject.ct[[]], Sample_ID, name = "n")
  if (nrow(donor_n2) < 10) {
    warning("After filtering, <10 donors remain for ", ct, " (", nrow(donor_n2), "). Skipping.")
    rm(dataObject.ct); invisible(gc()); next
  }
  
  #----------------------------------------------------------------------------
  # 3) Compute safe k for THIS cell type
  #----------------------------------------------------------------------------
  min_group_size <- min(donor_n2$n)
  safe_k <- max(k_floor, min(k_cap, min_group_size - 1))
  message(ct, " donors kept: ", nrow(donor_n2),
          " | min cells/donor: ", min_group_size,
          " | using k = ", safe_k)
  
  #----------------------------------------------------------------------------
  # 4) Setup WGCNA for THIS cell type object
  #    (important: do this on the object you will use for metacells/network)
  #----------------------------------------------------------------------------
  dataObject.ct <- SetupForWGCNA(
    dataObject.ct,
    gene_select = "custom",
    features = wgcna_genes,
    wgcna_name = "CWOW"
  )
  
  #----------------------------------------------------------------------------
  # 5) Build metacells within this cell type, grouped by Sample_ID
  #----------------------------------------------------------------------------
  dataObject.ct <- MetacellsByGroups(
    seurat_obj  = dataObject.ct,
    group.by    = "Sample_ID",
    reduction   = reduction_use,
    k           = safe_k,
    min_cells   = min_cells_ct,
    max_shared  = 10,
    ident.group = "Sample_ID"
  )
  
  dataObject.ct <- NormalizeMetacells(dataObject.ct)
  message("Metacells constructed successfully for: ", ct)
  
  #----------------------------------------------------------------------------
  # 6) Set datExpr using metacells (for soft power testing, use subset genes)
  #----------------------------------------------------------------------------
  dataObject.ct <- SetDatExpr(
    seurat_obj = dataObject.ct,
    assay = "RNA",
    layer = "data",
    use_metacells = TRUE
  )
  
  datExpr <- GetDatExpr(dataObject.ct)
  keep_genes <- intersect(colnames(datExpr), softpower_genes)
  dataObject.ct@misc$CWOW$datExpr <- datExpr[, keep_genes, drop = FALSE]
  
  #----------------------------------------------------------------------------
  # 7) Soft power testing
  #----------------------------------------------------------------------------
  dataObject.ct <- TestSoftPowers(
    dataObject.ct,
    networkType = "signed",
    powers = c(1:10, seq(12, 30, 2))
  )
  
  plot_list <- PlotSoftPowers(dataObject.ct)
  pdf(file.path(out_dir, paste0("SoftPower_", ct, ".pdf")), 10, 8)
  print(wrap_plots(plot_list, ncol = 2))
  dev.off()
  
  sft <- GetPowerTable(dataObject.ct)
  
  # Pick a soft power automatically (edit this rule as you like)
  # Rule: first power with SFT.R.sq >= 0.8; otherwise the power with max SFT.R.sq
  if ("SFT.R.sq" %in% colnames(sft)) {
    hit <- which(sft$SFT.R.sq >= 0.8)
    if (length(hit) > 0) {
      soft_power <- sft$Power[min(hit)]
    } else {
      soft_power <- sft$Power[which.max(sft$SFT.R.sq)]
    }
  } else {
    # fallback if table column naming differs
    soft_power <- sft$Power[which.max(sft[[2]])]
  }
  message("Chosen soft power for ", ct, ": ", soft_power)
  
  #----------------------------------------------------------------------------
  # 8) Restore full gene matrix for network construction (all WGCNA genes)
  #----------------------------------------------------------------------------
  dataObject.ct <- SetDatExpr(
    dataObject.ct,
    assay = "RNA",
    layer = "data",
    use_metacells = TRUE
  )
  
  #----------------------------------------------------------------------------
  # 9) Construct network
  #----------------------------------------------------------------------------
  dataObject.ct <- ConstructNetwork(
    dataObject.ct,
    tom_name = ct,
    overwrite_tom = TRUE,
    soft_power = soft_power,
    networkType = "signed"
  )
  
  pdf(file.path(out_dir, paste0("Dendrogram_", ct, ".pdf")), 7, 3.5)
  PlotDendrogram(dataObject.ct, main = paste(ct, "hdWGCNA"))
  dev.off()
  
  dataObject.ct <- ScaleData(
    dataObject.ct,
    assay = "RNA",
    features = GetWGCNAGenes(dataObject.ct),
    verbose = FALSE
  )
  
  #----------------------------------------------------------------------------
  # 10) Module eigengenes (harmonize by Sample_ID) Takes some time to run
  #----------------------------------------------------------------------------
  dataObject.ct <- ModuleEigengenes(
    dataObject.ct,
    group.by.vars = "Sample_ID"
  )
  
  #----------------------------------------------------------------------------
  # 11) kME connectivity (on this cell type object)
  #----------------------------------------------------------------------------
  dataObject.ct <- ModuleConnectivity(dataObject.ct)
  print("modules:")
  
  #----------------------------------------------------------------------------
  # 12) Rename modules
  #----------------------------------------------------------------------------
  dataObject.ct <- ResetModuleNames(
    dataObject.ct,
    new_name = paste0(ct, "-M")
  )
  
  #----------------------------------------------------------------------------
  # 13) Save modules + hubs
  #----------------------------------------------------------------------------
  modules <- GetModules(dataObject.ct) %>% dplyr::filter(module != "grey")
  write.table(
    modules,
    file.path(out_dir, paste0("modules_", ct, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  hubs <- GetHubGenes(dataObject.ct, n_hubs = 10)
  write.table(
    hubs,
    file.path(out_dir, paste0("hubs_", ct, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  #----------------------------------------------------------------------------
  # 14) Plots
  #----------------------------------------------------------------------------
  # first get module sizes, this will help with formating the figure sizes
  # gene count in modules
  modules_all <- GetModules(dataObject.ct)
  # subcluster count
  n_subclusters <- length(table(dataObject.ct$cell_type_subcluster))
  
  # count genes per module (exclude grey by default)
  module_sizes <- modules_all %>%
    dplyr::count(module, name = "n_genes") %>%
    dplyr::filter(module != "grey") %>%
    dplyr::arrange(dplyr::desc(n_genes)) %>%
    dplyr::mutate(module = factor(module, levels = module))
  
  # Plot KMEs 
  pdf(file.path(out_dir, paste0("kME_", ct, ".pdf")), 14, height = .35 * nrow(module_sizes) + 2.5)
  plot_kMEs <- PlotKMEs(dataObject.ct, ncol = 5)
  print(plot_kMEs)
  dev.off()
  
  # Feature plot of hMEs should use the CT object
  plot_list2 <- ModuleFeaturePlot(
    dataObject.ct,
    features = "hMEs",
    reduction = "umap.rpca",
    order = TRUE
  )
  pdf(file.path(out_dir, paste0("Feature_hMEs_", ct, ".pdf")), 10, height = .6 * nrow(module_sizes) + 3)
  print(wrap_plots(plot_list2, ncol = 5))
  dev.off()
  
  #-- radar plots
  # disease group
  dataObject.ct$group <- do.call(rbind, strsplit(as.character(dataObject.ct$group), ' '))[,1]
  ModuleRadarPlot_group <-  ModuleRadarPlot(
    dataObject.ct,
    group.by = 'group',
    axis.label.size=2.5,
    grid.label.size=3
  )
  pdf(file.path(out_dir, paste0("RadarPlot_", ct, "_by_group.pdf")), width = 12, height = .6 * nrow(module_sizes) + 3)
  print(ModuleRadarPlot_group)
  dev.off()
  
  # subclusters
  dataObject.ct$cell_type_subcluster <- do.call(rbind, strsplit(as.character(dataObject.ct$cell_type_subcluster), ' '))[,1]
  ModuleRadarPlot_subcluster <- ModuleRadarPlot(
    dataObject.ct,
    group.by = 'cell_type_subcluster',
    axis.label.size=2.5,
    grid.label.size=3
  )
  pdf(file.path(out_dir, paste0("RadarPlot_", ct, "_by_subcluster.pdf")), 12, height = .5 * nrow(module_sizes) + 2.5)
  print(ModuleRadarPlot_subcluster)
  dev.off()

  # get hMEs from seurat object
  MEs <- GetMEs(dataObject.ct, harmonized=TRUE)
  modules <- GetModules(dataObject.ct)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  
  # add hMEs to Seurat meta-data:
  dataObject.ct@meta.data <- cbind(dataObject.ct@meta.data, MEs)
  
  # plot with Seurat's DotPlot function
  plot_dotplot <- DotPlot(dataObject.ct, features=mods, group.by = 'cell_type_subcluster')
  # flip the x/y axes, rotate the axis labels, and change color scheme:
  plot_dotplot <- plot_dotplot +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')
  pdf(file.path(out_dir, paste0("DotPlot_", ct, "_by_subcluster.pdf")), width = .5 * nrow(module_sizes) + 2.5, height = .4 * n_subclusters + 2)
  print(plot_dotplot)
  dev.off()

  # module trait correlation
  # must be numeric
  dataObject.ct$Age <- as.numeric(dataObject.ct$Age)
  dataObject.ct$Thal.amyloid <- as.numeric(dataObject.ct$Thal.amyloid)
  dataObject.ct$Braak.NFT <- as.numeric(dataObject.ct$Braak.NFT)
  dataObject.ct$Cing.LB <- as.numeric(dataObject.ct$Cing.LB)
  dataObject.ct$Sex_female <- ifelse(dataObject.ct$sex_inferred == "female", 1, 0)
  dataObject.ct$group <- factor(
    dataObject.ct$group,
    levels = c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS")
  )
  dataObject.ct$AD_AT  <- ifelse(dataObject.ct$group == "AD_AT",  1,
                                         ifelse(dataObject.ct$group == "CONTROL", 0, NA))
  dataObject.ct$LBD_S  <- ifelse(dataObject.ct$group == "LBD_S",  1,
                                         ifelse(dataObject.ct$group == "CONTROL", 0, NA))
  dataObject.ct$LBD_AS <- ifelse(dataObject.ct$group == "LBD_AS", 1,
                                         ifelse(dataObject.ct$group == "CONTROL", 0, NA))
  dataObject.ct$LBD_ATS <- ifelse(dataObject.ct$group == "LBD_ATS", 1,
                                          ifelse(dataObject.ct$group == "CONTROL", 0, NA))
  cur_traits <- c(
    "AD_AT",
    "LBD_S",
    "LBD_AS",
    "LBD_ATS",
    "Thal.amyloid",
    "Braak.NFT",
    "Cing.LB",
    "Sex_female",     
    "Age"
  )  
  
  # By cell type 
  dataObject.ct <- ModuleTraitCorrelation(
    dataObject.ct,
    traits = cur_traits
  )
  # get the mt-correlation results
  mt_cor <- GetModuleTraitCorrelation(dataObject.ct)
  names(mt_cor)
  names(mt_cor$cor)
  # note about sex
  # positive correlation = higher in females
  # negative correlation = higher in males
  plot_MT_ct <- PlotModuleTraitCorrelation(
    dataObject.ct,
    label = 'fdr',
    label_symbol = 'stars',
    text_size = 2,
    text_digits = 2,
    text_color = "white",
    high_color = "#B2182B",   # deep red
    mid_color  = "#FFFFFF",  # white
    low_color  = "#2166AC",  # deep blue
    plot_max = 0.2,
    combine=TRUE
  )

  pdf(file.path(out_dir, paste0("ModuleTrait_", ct, "_by_cell_type.pdf")), width = .5 * nrow(module_sizes) + 3, 6)
  print(plot_MT_ct)
  dev.off()

  # By subclusters
  dataObject.ct <- ModuleTraitCorrelation(
    dataObject.ct,
    traits = cur_traits, 
    group.by = "cell_type_subcluster"
  )
  # get the mt-correlation results
  mt_cor <- GetModuleTraitCorrelation(dataObject.ct)
  names(mt_cor)
  names(mt_cor$cor)
  module_order <- colnames(mt_cor$cor$all_cells)
  
  # height scaling parameters (tweak once, reuse everywhere)
  height_per_module <- 1.5   # inches per module row
  min_height        <- 6
  max_height        <- 20

  pdf_height <- min(
    max_height,
    max(min_height, height_per_module * n_subclusters + 4)
  )
  # note about sex
  # positive correlation = higher in females
  # negative correlation = higher in males
  plot_MT_ct_subcluster <- PlotModuleTraitCorrelation(
    dataObject.ct,
    label = 'fdr',
    label_symbol = 'stars',
    text_size = 2,
    text_digits = 2,
    text_color = "white",
    high_color = "#B2182B",   # deep red
    mid_color  = "#FFFFFF",  # white
    low_color  = "#2166AC",  # deep blue
    plot_max = 0.2,
    combine=TRUE
  )
  pdf(
    file.path(out_dir, paste0("ModuleTrait_", ct, "_by_cell_type_subcluster.pdf")),
    width = .5 * nrow(module_sizes) + 3,
    height = pdf_height
  )
  print(plot_MT_ct_subcluster)
  dev.off()

  # module sizes
  df <- as.data.frame(module_sizes)
  module_order_df <- as.data.frame(module_order)
  df$module <- factor(df$module, levels = module_order_df$module_order)  # save table
  write.table(
    module_sizes,
    file.path(out_dir, paste0("ModuleSizes_", ct, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # plot
  p_module_sizes <- ggplot(df, aes(x = module, y = n_genes)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste0(ct, ": genes per module (excluding grey)"),
      x = "Module",
      y = "Number of genes"
    ) +
    theme_cowplot() +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      plot.title  = element_text(size = 10)
    )
  
  pdf(file.path(out_dir, paste0("ModuleSizes_", ct, ".pdf")), width = 7, height = 0.18 * nrow(module_sizes) + 2.5)
  print(p_module_sizes)
  dev.off()

  #----------------------------------------------------------------------------
  # 15) Save final objects
  #----------------------------------------------------------------------------
  qsave(
    dataObject.ct,
    file.path(obj_dir, paste0("hdWGCNA_", ct, "_final.qs"))
  )
  
  saveRDS(
    dataObject.ct,
    file.path(obj_dir, paste0("hdWGCNA_", ct, "_final.rds")),
    compress = FALSE
  )
  
  #----------------------------------------------------------------------------
  # 16) Cleanup
  #----------------------------------------------------------------------------
  rm(dataObject.ct)
  rm(
    donor_n,
    donor_n2,
    hubs,
    MEs,
    meta,
    meta_ct,
    metadata,
    modules,
    mt_cor,
    plot_kMEs,
    plot_list,
    plot_list2,
    plot_MT_ct,
    sft,
    ModuleRadarPlot_group,
    ModuleRadarPlot_subcluster,
    plot_dotplot,
    plot_MT_ct_subcluster,
    df
  )
  invisible(gc())
  invisible(gc())
  Sys.sleep(2)
}

message("ALL CELL TYPES COMPLETED SUCCESSFULLY")


