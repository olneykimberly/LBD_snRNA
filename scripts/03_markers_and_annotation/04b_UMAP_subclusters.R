library(Seurat)
library(ggplot2)
library(patchwork)

source("file_paths_and_colours.R")

# Define cell types in the specified order
order_cell_type <- c("neuron", "interneuron", "oligodendrocyte", "opc", 
                     "astrocyte", "microglia", "endothelial", "fibroblast", "mural")

# Define base path for the RDS files based on your directory structure
base_path <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/subclustering_cell_type_pass1/"

# Note: Ensure 'color_panel' is defined in your environment before running this loop
# Example: color_panel <- DiscretePalette(n = 40, palette = "polychrome")

for (cell_type in order_cell_type) {
  
  # Construct file path
  file_path <- paste0(base_path, cell_type, "/CWOW_cellbender_", cell_type, "_subclustered.rds")
  
  if (file.exists(file_path)) {
    message(paste("Processing:", cell_type))
    
    # Read the RDS file
    obj <- readRDS(file_path)
    
    # Common theme components for consistency
    common_theme <- theme(
      axis.title.x = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "none",
      plot.title = element_text(size = 8),
      plot.margin = margin(0.5, 0.025, 0.25, 0.025, "cm")
    )
    
    # UMAP 1: Grouped by Sample_ID
    p1 <- DimPlot(obj, reduction = "umap.ct", label = TRUE, label.size = 3, group.by = "Sample_ID") + 
      xlab("UMAP 1") + ylab("UMAP 2") + ggtitle(paste(cell_type, "- By Sample")) +
      common_theme +
      scale_color_manual(values = color_panel)
    
    # UMAP 2: Default (Grouped by Seurat Clusters)
    p2 <- DimPlot(obj, reduction = "umap.ct", label = TRUE, label.size = 3) + 
      xlab("UMAP 1") + ylab("UMAP 2") + ggtitle(paste(cell_type, "- Clusters")) +
      common_theme +
      scale_color_manual(values = color_panel)
    
    # Combine plots side-by-side and save
    pdf_name <- paste0(base_path, "/UMAP_clusters_", cell_type, ".pdf")
    pdf(pdf_name, width = 12, height = 6)
    print(p1 + p2)
    dev.off()
    
  } else {
    warning(paste("File not found for cell type:", cell_type))
  }
}