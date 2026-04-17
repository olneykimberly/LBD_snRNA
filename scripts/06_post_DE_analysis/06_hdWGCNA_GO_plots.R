#!/usr/bin/env Rscript

# ==============================================================================
# hdWGCNA Final GO Plot: Continuous Side-by-Side (Zero Gap)
# Features: No white space, M numbers on top, Vertical terms, Contrast fixed
# ==============================================================================

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(qs)
  library(dplyr)
  library(stringr)
  library(enrichR)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# --- Settings ---
cell_types <- c("microglia")
wgcna_name <- "CWOW"
obj_dir <- "../rObjects"
out_dir <- "../results/hdWGCNA"
requested_db <- "GO_Biological_Process_2023"
font_size <- 8

# Helper: Force white text for dark modules (Black, Blue, Brown)
get_text_contrast <- function(color_name) {
  dark_colors <- c("black", "blue", "brown", "darkred", "darkgreen", "purple", "magenta")
  if (color_name %in% dark_colors) return("white")
  return("black")
}

for (ct in cell_types) {
  message("\n--- GENERATING CONTINUOUS GO PLOT: ", ct, " ---")
  seurat_obj <- qs::qread(file.path(obj_dir, paste0("hdWGCNA_", ct, "_final.qs")))
  
  # 1. Prepare Data and Dendrogram Order
  modules <- hdWGCNA::GetModules(seurat_obj)
  dendro_obj <- seurat_obj@misc[[wgcna_name]]$wgcna_net$dendrograms[[1]]
  dendro_obj$labels <- seurat_obj@misc[[wgcna_name]]$wgcna_genes
  
  # Identify the left-to-right order of modules as they appear in the tree
  gene_order <- dendro_obj$labels[dendro_obj$order]
  
  # Determine module order by the first appearance of each module in the gene_order
  mod_order <- unique(modules$module[match(gene_order, modules$gene_name)])
  mod_order <- setdiff(mod_order, "grey")
  
  go_df <- data.frame()
  for(mod in mod_order){
    mod_genes <- modules %>% filter(module == mod) %>% pull(gene_name)
    enr <- tryCatch({ enrichr(mod_genes, requested_db) }, error = function(e) NULL)
    
    if(!is.null(enr) && nrow(enr[[requested_db]]) > 0){
      term <- enr[[requested_db]]$Term[1] %>% 
              str_replace_all(" \\(GO:.*\\)", "") %>% 
              str_trunc(40)
      
      mod_color <- modules$color[modules$module == mod][1]
      
      go_df <- rbind(go_df, data.frame(
        module = factor(mod, levels = mod_order), # Preserves dendrogram order
        m_num = paste0("M", gsub(".*-M", "", mod)),
        term = term,
        color = mod_color,
        text_col = get_text_contrast(mod_color)
      ))
    }
  }

  # 2. Create the Plot: Side-by-Side with Zero Gap
  p_go <- ggplot(go_df, aes(x = module, y = 1)) +
    # geom_tile with width=1 ensures boxes touch each other perfectly
    geom_tile(aes(fill = color), width = 1, height = 1) +
    
    # M numbers placed at the top of the tile
    geom_text(aes(label = m_num, color = text_col), 
              y = 1.4, size = font_size/.pt, fontface = "bold") +
    
    # GO terms placed vertically inside the tiles
    geom_text(aes(label = term, color = text_col),
              y = 0.85, 
              angle = 90, 
              size = (font_size - 1.5)/.pt, 
              hjust = 0.5) +
    
    scale_fill_identity() +
    scale_color_identity() +
    # Expand = c(0,0) removes the white margin padding at the edges of the boxes
    scale_x_discrete(expand = c(0, 0)) + 
    ylim(0.5, 1.5) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) # Complete removal of plot-level margins

  # 3. Save Final Version
  save_path <- file.path(out_dir, paste0("GO_Terms_", ct, "_Continuous.pdf"))
  pdf(save_path, width = 6.5, height = 3.5)
  print(p_go)
  dev.off()
  
  message("Success: Continuous side-by-side GO plot saved to ", save_path)
}