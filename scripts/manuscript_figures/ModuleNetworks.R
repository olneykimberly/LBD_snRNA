setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE, character.only = TRUE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

theme_set(cowplot::theme_cowplot(font_size = 7))

# ----------------------------- paths ------------------------------------------
obj_dir       <- "../rObjects"
base_out_dir  <- "../results/hdWGCNA"
network_dir   <- file.path(base_out_dir, "ModuleNetworks")

dir.create(network_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- settings ---------------------------------------
wgcna_name <- "CWOW"

cell_types <- c(
#  "neuron",
#  "interneuron",
#  "oligodendrocyte",
#  "opc"
#  "astrocyte",
#  "microglia"
  #"endothelial",
  "fibroblast",
  "mural"
)

# ----------------------------- helper -----------------------------------------
fix_yellow_to_gold <- function(seu, wgcna_name = "CWOW") {
  
  mod_df <- seu@misc[[wgcna_name]]$wgcna_modules
  
  if (!is.null(mod_df) && "color" %in% colnames(mod_df)) {
    mod_df$color <- ifelse(
      tolower(mod_df$color) == "yellow",
      "gold",
      mod_df$color
    )
    seu@misc[[wgcna_name]]$wgcna_modules <- mod_df
  }
  
  return(seu)
}

# ==============================================================================
# MAIN LOOP
# ==============================================================================
for (ct in cell_types) {
  
  message("---- ", ct, " ----")
  
  seu_path <- file.path(obj_dir, paste0("hdWGCNA_", ct, "_final.qs"))
  
  if (!file.exists(seu_path)) {
    warning("Missing object: ", seu_path)
    next
  }
  
  seu <- qread(seu_path)
  
  # Change yellow module color to gold
  seu <- fix_yellow_to_gold(seu, wgcna_name = wgcna_name)
  
  # Output directory for this cell type
  ct_outdir <- file.path(network_dir, ct)
  dir.create(ct_outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Get module names
  mod_df <- seu@misc[[wgcna_name]]$wgcna_modules
  
  modules <- mod_df %>%
    dplyr::filter(module != "grey") %>%
    dplyr::pull(module) %>%
    unique()
  
  message("Modules found: ", paste(modules, collapse = ", "))
  
  # ---------------------------------------------------------------------------
  # Make one network plot per module
  # ---------------------------------------------------------------------------
  for (mod in modules) {
    
    message("  plotting ", mod)
    
    mod_outdir <- file.path(ct_outdir, mod)
    dir.create(mod_outdir, showWarnings = FALSE, recursive = TRUE)
    
    ModuleNetworkPlot(
      seurat_obj = seu,
      modules = mod,
      outdir = mod_outdir,
      wgcna_name = wgcna_name
    )
  }
}