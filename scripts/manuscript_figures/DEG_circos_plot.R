# =========================================================
# Circos overlap plot: shared/unique UP-regulated DEGs
# Cell type: neuron / excitatory neuron
# =========================================================

library(tidyverse)
library(circlize)
library(scales)

# ------------------------------
# Directories and parameters
# ------------------------------
main_output_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/DEGs_MAST_RNA_pct0.25"
output_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/circos_DEG_overlap"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

celltype <- "neuron"

up_logfc_cutoff <- 0.25
padj_cutoff <- 0.05

comparisons_to_read <- list(
  AD  = c("AD_AT",   "CONTROL"),
  S   = c("LBD_S",   "CONTROL"),
  AS  = c("LBD_AS",  "CONTROL"),
  ATS = c("LBD_ATS", "CONTROL")
)

pretty_names <- c(
  AD  = "AD(AT)",
  S   = "LBD(S)",
  AS  = "LBD(AS)",
  ATS = "LBD(ATS)"
)

group_cols <- c(
  "AD(AT)"   = "#B4464B",
  "LBD(S)"   = "gray",
  "LBD(AS)"  = "gray65",
  "LBD(ATS)" = "gray35"
)

# ------------------------------
# Helper: read DEG file
# ------------------------------
read_deg_file <- function(celltype, group1, group2, main_output_dir) {
  
  fname <- file.path(
    main_output_dir,
    celltype,
    sprintf("DEG_%s_%s_vs_%s.tsv", celltype, group1, group2)
  )
  
  if (!file.exists(fname)) {
    stop("File not found: ", fname)
  }
  
  df <- readr::read_tsv(fname, show_col_types = FALSE)
  
  required <- c("gene", "avg_log2FC", "p_val_adj")
  if (!all(required %in% colnames(df))) {
    stop("Missing required columns in: ", fname)
  }
  
  df %>%
    dplyr::select(gene, avg_log2FC, p_val_adj)
}

# ------------------------------
# Get UP-regulated DEG lists
# ------------------------------
deg_lists <- lapply(names(comparisons_to_read), function(group_label) {
  
  pair <- comparisons_to_read[[group_label]]
  
  read_deg_file(
    celltype = celltype,
    group1 = pair[1],
    group2 = pair[2],
    main_output_dir = main_output_dir
  ) %>%
    filter(
      !is.na(gene),
      !is.na(avg_log2FC),
      !is.na(p_val_adj),
      avg_log2FC >= up_logfc_cutoff,
      p_val_adj <= padj_cutoff
    ) %>%
    pull(gene) %>%
    unique()
})

names(deg_lists) <- pretty_names[names(comparisons_to_read)]

# ------------------------------
# Save long DEG list
# ------------------------------
deg_long <- tibble(
  group = names(deg_lists),
  genes = deg_lists
) %>%
  unnest_longer(genes, values_to = "gene")

write_tsv(
  deg_long,
  file.path(output_dir, paste0(celltype, "_upregulated_DEG_lists_long.tsv"))
)

# ------------------------------
# Build pairwise overlap table
# ------------------------------
make_pairwise_links <- function(gene_lists) {
  
  group_names <- names(gene_lists)
  links <- list()
  
  for (i in seq_along(group_names)) {
    for (j in seq_along(group_names)) {
      
      if (i < j) {
        
        g1 <- group_names[i]
        g2 <- group_names[j]
        
        shared_genes <- intersect(gene_lists[[g1]], gene_lists[[g2]])
        
        if (length(shared_genes) > 0) {
          links[[paste(g1, g2, sep = "_")]] <- tibble(
            from = g1,
            to = g2,
            value = length(shared_genes),
            genes = paste(shared_genes, collapse = ";")
          )
        }
      }
    }
  }
  
  bind_rows(links)
}

links_df <- make_pairwise_links(deg_lists)

write_tsv(
  links_df,
  file.path(output_dir, paste0(celltype, "_upregulated_pairwise_overlap_links.tsv"))
)

# ------------------------------
# Build overlap matrix
# ------------------------------
overlap_mat <- matrix(
  0,
  nrow = length(deg_lists),
  ncol = length(deg_lists),
  dimnames = list(names(deg_lists), names(deg_lists))
)

for (i in names(deg_lists)) {
  for (j in names(deg_lists)) {
    if (i != j) {
      overlap_mat[i, j] <- length(intersect(deg_lists[[i]], deg_lists[[j]]))
    }
  }
}

diag(overlap_mat) <- lengths(deg_lists)

write.csv(
  overlap_mat,
  file.path(output_dir, paste0(celltype, "_upregulated_overlap_matrix.csv"))
)

# ------------------------------
# Helper: make annotation label positions
# ------------------------------
make_label_df <- function(links_df, group_names) {
  
  if (nrow(links_df) == 0) {
    return(tibble())
  }
  
  links_df %>%
    mutate(
      label = as.character(value),
      from_id = match(from, group_names),
      to_id = match(to, group_names),
      pair_id = row_number()
    )
}

label_df <- make_label_df(links_df, names(deg_lists))

# ------------------------------
# Plot function
# ------------------------------
plot_deg_circos <- function(outfile, device = c("pdf", "png")) {
  
  device <- match.arg(device)
  
  if (device == "pdf") {
    pdf(outfile, width = 7, height = 7)
  } else {
    png(outfile, width = 2400, height = 2400, res = 300)
  }
  
  circos.clear()
  
  circos.par(
    start.degree = 90,
    gap.degree = 6,
    track.margin = c(0.01, 0.01),
    canvas.xlim = c(-1.35, 1.35),
    canvas.ylim = c(-1.35, 1.35)
  )
  
  chordDiagram(
    x = overlap_mat,
    grid.col = group_cols,
    transparency = 0.35,
    directional = 0,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.12)
  )
  
  # group labels
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      
      sector_name <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      
      circos.text(
        x = mean(xlim),
        y = ylim[1] + 1.2,
        labels = sector_name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1.15
      )
    },
    bg.border = NA
  )
  
  # add overlap counts as labels near the middle of each pairwise link
  if (nrow(label_df) > 0) {
    
    for (k in seq_len(nrow(label_df))) {
      
      g1 <- label_df$from[k]
      g2 <- label_df$to[k]
      n_shared <- label_df$value[k]
      
      # approximate midpoint using the first sector position
      xlim1 <- get.cell.meta.data("xlim", sector.index = g1, track.index = 1)
      x_mid <- mean(xlim1)
      
      # stagger labels slightly to reduce overplotting
      y_pos <- 0.55 + 0.08 * ((k - 1) %% 3)
      
      circos.text(
        sector.index = g1,
        track.index = 1,
        x = x_mid,
        y = y_pos,
        labels = n_shared,
        cex = 0.75,
        col = "black",
        facing = "inside",
        niceFacing = TRUE
      )
    }
  }
  
  title(
    main = paste0(celltype, " upregulated DEG overlap"),
    cex.main = 1
  )
  
  circos.clear()
  dev.off()
}

# ------------------------------
# Save plots
# ------------------------------
plot_deg_circos(
  outfile = file.path(output_dir, paste0(celltype, "_upregulated_DEG_circos_overlap_labeled.pdf")),
  device = "pdf"
)

plot_deg_circos(
  outfile = file.path(output_dir, paste0(celltype, "_upregulated_DEG_circos_overlap_labeled.png")),
  device = "png"
)

message("Done. Circos plots and overlap tables saved to: ", output_dir)