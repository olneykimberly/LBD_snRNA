setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")

# Protein coding genes
genes <- readRDS("../rObjects/annotation.rds")
protein_coding <- subset(genes, gene_type == "protein_coding")$gene_name
mt_genes <- subset(genes, seqnames == "chrM")$gene_name
features_use <- intersect(setdiff(protein_coding, mt_genes), rownames(obj))
# ------------------------------------------------------------------------------
# Read cell-type pathology DEG tables into a named list
# ------------------------------------------------------------------------------

# thresholds (you'll use these later for "significant")
qval <- 0.05
FC   <- 0.25

# comparisons (as you defined)
comparisons <- list(
  c("AD_AT", "CONTROL"),
  c("LBD_S", "CONTROL"),
  c("LBD_AS", "CONTROL"),
  c("LBD_ATS", "CONTROL"),
  c("LBD_S", "AD_AT"),
  c("LBD_AS", "AD_AT"),
  c("LBD_ATS", "AD_AT"),
  c("LBD_AS", "LBD_S"),
  c("LBD_ATS", "LBD_S"),
  c("LBD_ATS", "LBD_AS")
)

deg_base_dir <- "../results/DEGs_MAST_RNA_pct0.25"

# storage
deg_tables <- list()   # named list of data.frames
deg_index  <- list()   # quick metadata / QC log

required_cols <- c("avg_log2FC", "p_val", "p_val_adj", "gene")

for (cell_type in cell_types) {
  for (comp in comparisons) {
    
    pathology_A <- comp[[1]]
    pathology_B <- comp[[2]]
    
    # Construct file path exactly like your format
    file_name <- paste0("DEG_", cell_type, "_", pathology_A, "_vs_", pathology_B, ".tsv")
    file_path <- file.path(deg_base_dir, cell_type, file_name)
    
    key <- paste(cell_type, pathology_A, "vs", pathology_B, sep = "__")
    
    if (!file.exists(file_path)) {
      message("Skipping (missing): ", file_path)
      deg_index[[key]] <- data.frame(
        key = key, cell_type = cell_type,
        pathology_A = pathology_A, pathology_B = pathology_B,
        file_path = file_path, status = "missing",
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Read
    dat <- tryCatch(
      read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) e
    )
    if (inherits(dat, "error")) {
      message("Skipping (read error): ", file_path, " :: ", dat$message)
      deg_index[[key]] <- data.frame(
        key = key, cell_type = cell_type,
        pathology_A = pathology_A, pathology_B = pathology_B,
        file_path = file_path, status = "read_error",
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Drop NA rows (optional; keep if you prefer)
    dat <- na.omit(dat)
    
    # Column check
    missing_cols <- setdiff(required_cols, colnames(dat))
    if (length(missing_cols) > 0L) {
      message("Skipping (missing cols): ", file_path, " :: ",
              paste(missing_cols, collapse = ", "))
      deg_index[[key]] <- data.frame(
        key = key, cell_type = cell_type,
        pathology_A = pathology_A, pathology_B = pathology_B,
        file_path = file_path, status = paste0("missing_cols: ", paste(missing_cols, collapse = ",")),
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Coerce numeric columns safely (sometimes read.delim makes them character)
    dat$avg_log2FC <- suppressWarnings(as.numeric(dat$avg_log2FC))
    dat$p_val      <- suppressWarnings(as.numeric(dat$p_val))
    dat$p_val_adj  <- suppressWarnings(as.numeric(dat$p_val_adj))
    
    # store
    deg_tables[[key]] <- dat
    
    # optional quick counts (helps sanity check)
    n_sig <- sum(dat$p_val_adj < qval & abs(dat$avg_log2FC) > FC, na.rm = TRUE)
    
    deg_index[[key]] <- data.frame(
      key = key, cell_type = cell_type,
      pathology_A = pathology_A, pathology_B = pathology_B,
      file_path = file_path, status = "ok",
      n_rows = nrow(dat),
      n_sig = n_sig,
      stringsAsFactors = FALSE
    )
  }
}

# Combine the index into a single data.frame
deg_index_df <- do.call(rbind, deg_index)

# What you now have:
# - deg_tables[[key]] is the DEG table for that cell type + comparison
# - deg_index_df summarizes which files were read / skipped and sig counts
print(deg_index_df)

#--------------------------------------
# Bulk RNAseq
# Columns expected:
# seqnames start end width gene_id gene_type gene_name logFC CI.L CI.R AveExpr t P.Value adj.P.Val B
# ------------------------------------------------------------------------------

bulk_dirs <- c(
  "/tgen_labs/jfryer/kolney/LBD_CWOW/bulkRNA/results/star/DEGs",
  "/tgen_labs/jfryer/kolney/LBD_CWOW/bulkRNA/results/star_LBD_redefined/DEGs"
)

bulk_files <- c(
  # star/DEGs
  "TYPE_ADvsControl_gene_DEGs_FDRq0.05_logFC_0.25.txt",
  "TYPE_LBDvsAD_gene_DEGs_FDRq0.05_logFC_0.25.txt",
  "TYPE_LBDvsControl_gene_DEGs_FDRq0.05_logFC_0.25.txt",
  
  # star_LBD_redefined/DEGs (S/AS/ATS)
  "LBD_redefined_LBD_S_vs_AD_gene_DEGs_FDRq0.05_logFC_0.25.txt",
  "LBD_redefined_LBD_S_vs_Control_gene_DEGs_FDRq0.05_logFC_0.25.txt",
  "LBD_redefined_LBD_AS_vs_AD_gene_DEGs_FDRq0.05_logFC_0.25.txt",
  "LBD_redefined_LBD_AS_vs_Control_gene_DEGs_FDRq0.05_logFC_0.25.txt",
  "LBD_redefined_LBD_ATS_vs_AD_gene_DEGs_FDRq0.05_logFC_0.25.txt",
  "LBD_redefined_LBD_ATS_vs_Control_gene_DEGs_FDRq0.05_logFC_0.25.txt"
)

resolve_first_existing <- function(fname, dirs) {
  paths <- file.path(dirs, fname)
  ok <- file.exists(paths)
  if (!any(ok)) return(NA_character_)
  paths[which(ok)[1]]
}

read_bulk_deg <- function(file_path) {
  # strict read
  df <- read.delim(
    file_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # required columns (tight)
  req <- c("gene_name", "logFC")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0L) {
    stop("Missing required columns: ", paste(miss, collapse = ", "),
         " in file: ", file_path)
  }
  
  # standardize
  df$gene  <- as.character(df$gene_name)
  df$logFC <- suppressWarnings(as.numeric(df$logFC))
  
  # drop junk
  df <- df[!is.na(df$gene) & df$gene != "" & !is.na(df$logFC), , drop = FALSE]
  
  # de-duplicate by gene (keep largest |logFC|)
  if (any(duplicated(df$gene))) {
    ord <- order(df$gene, -abs(df$logFC))
    df <- df[ord, , drop = FALSE]
    df <- df[!duplicated(df$gene), , drop = FALSE]
  }
  
  # keep a clean, consistent column set (you can add/remove here)
  keep <- intersect(
    c("gene", "logFC", "adj.P.Val", "P.Value", "AveExpr",
      "CI.L", "CI.R", "t", "B", "gene_id", "gene_type",
      "seqnames", "start", "end", "width"),
    colnames(df)
  )
  df <- df[, keep, drop = FALSE]
  
  df
}

# Read all files
bulk_deg_tables <- list()
bulk_index <- list()

for (fname in bulk_files) {
  key <- sub("\\.txt$", "", fname)
  fpath <- resolve_first_existing(fname, bulk_dirs)
  
  if (is.na(fpath)) {
    message("Skipping (missing): ", fname)
    bulk_index[[key]] <- data.frame(
      key = key, file = fname, file_path = NA_character_,
      status = "missing",
      stringsAsFactors = FALSE
    )
    next
  }
  
  df <- tryCatch(
    read_bulk_deg(fpath),
    error = function(e) e
  )
  
  if (inherits(df, "error")) {
    message("Skipping (parse error): ", fpath, " :: ", df$message)
    bulk_index[[key]] <- data.frame(
      key = key, file = fname, file_path = fpath,
      status = "parse_error",
      stringsAsFactors = FALSE
    )
    next
  }
  
  bulk_deg_tables[[key]] <- df
  bulk_index[[key]] <- data.frame(
    key = key, file = fname, file_path = fpath,
    status = "ok",
    n_rows = nrow(df),
    n_genes = length(unique(df$gene)),
    stringsAsFactors = FALSE
  )
}

bulk_index_df <- do.call(rbind, bulk_index)
print(bulk_index_df)

# Convenience: gene vectors for overlap (already significant)
bulk_genes <- lapply(bulk_deg_tables, function(df) sort(unique(df$gene)))

# Example:
# names(bulk_deg_tables)
# head(bulk_deg_tables[[1]])
# length(bulk_genes[["TYPE_ADvsControl_gene_DEGs_FDRq0.05_logFC_0.25"]])

#--------
# ------------------------------------------------------------------------------
# AD vs Control overlap: snRNA cell-type DEGs vs bulk RNA-seq DEGs
# Outputs:
#   1) Multi-page PDF of 2-set Venns (one per cell type)
#   2) Bar plot of fraction overlap per cell type
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Bulk vs snRNA overlap for ALL comparisons
#  - One Venn per comparison: Bulk vs UNION(sn cell-type DEGs)
#  - Two fraction plots per comparison:
#       A) overlap / sn cell-type DEGs  (frac_sn_in_bulk)
#       B) overlap / bulk DEGs          (frac_bulk_in_sn)
#  - Output table of shared genes per comparison (bulk ∩ sn_union)
#  - Output per-cell-type overlap summary tables per comparison
#
# Prereqs already in your session:
#   - deg_tables: named list of sn DEG data.frames keyed as "celltype__A__vs__B"
#                columns include: gene, p_val_adj, avg_log2FC
#   - bulk_deg_tables: named list of bulk DEG data.frames, each with column: gene, logFC
#   - cell_types: character vector in desired order (neuron ... endothelial)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Bulk vs snRNA overlap for ALL comparisons
#  - One Venn per comparison: Bulk vs UNION(sn cell-type DEGs)
#  - Two fraction plots per comparison:
#       A) overlap / sn cell-type DEGs  (frac_sn_in_bulk)
#       B) overlap / bulk DEGs          (frac_bulk_in_sn)
#  - Output table of shared genes per comparison (bulk ∩ sn_union)
#  - Output per-cell-type overlap summary tables per comparison
#
# Prereqs already in your session:
#   - deg_tables: named list of sn DEG data.frames keyed as "celltype__A__vs__B"
#                columns include: gene, p_val_adj, avg_log2FC
#   - bulk_deg_tables: named list of bulk DEG data.frames, each with column: gene, logFC
#   - cell_types: character vector in desired order (neuron ... endothelial)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Bulk vs snRNA overlap for ALL comparisons (UPDATED per your requests)
#  - One Venn per comparison: Bulk vs UNION(sn cell-type DEGs)
#     * labels: "all cell types" (pink) and "bulk RNAseq" (green)
#     * region labels show counts + percentages: "n (p%)"
#  - Three per-comparison barplots:
#       A) frac_sn_in_bulk  = overlap / sn cell-type DEGs
#       B) frac_bulk_in_sn  = overlap / bulk DEGs
#       C) n_overlap        = total overlap count per cell type
#     * bars colored by cell type using color_panel
#  - Output shared gene table per comparison (bulk ∩ sn_union)
#
# Prereqs already in your session:
#   - deg_tables: named list of sn DEG dfs keyed "celltype__A__vs__B"
#                columns: gene, p_val_adj, avg_log2FC
#   - bulk_deg_tables: named list of bulk DEG dfs with columns incl: gene, logFC
#   - cell_types: character vector in desired order (neuron ... endothelial)
#   - color_panel: named vector mapping cell_type -> color (from file_paths_and_colours.R)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---- thresholds for sn filtering (bulk already filtered) ----
qval <- 0.05
FC   <- 0.25

# ---- comparisons you want ----
comparisons <- list(
  c("AD_AT", "CONTROL"),
  c("LBD_S", "CONTROL"),
  c("LBD_AS", "CONTROL"),
  c("LBD_ATS", "CONTROL"),
  c("LBD_S", "AD_AT"),
  c("LBD_AS", "AD_AT"),
  c("LBD_ATS", "AD_AT"),
  c("LBD_AS", "LBD_S"),
  c("LBD_ATS", "LBD_S"),
  c("LBD_ATS", "LBD_AS")
)

# ---- output dir ----
out_dir <- "../results/DEG_overlap_bulk_vs_sn_allComparisons"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- enforce order exactly as cell_types; neuron should appear at TOP in plots ----
cell_types <- as.character(cell_types)
# Ensure correct length
stopifnot(length(color_panel) >= length(cell_types))

# Assign names in the correct order
color_panel <- color_panel[seq_along(cell_types)]
names(color_panel) <- cell_types

stopifnot(all(cell_types %in% names(color_panel)))  # ensure colors exist for all

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || any(is.na(x))) y else x

get_sn_deg_genes <- function(deg_tables, cell_type, A, B, qval = 0.05, FC = 0.25) {
  key <- paste(cell_type, A, "vs", B, sep = "__")
  if (!key %in% names(deg_tables)) return(character(0))
  
  df <- deg_tables[[key]]
  req <- c("gene", "p_val_adj", "avg_log2FC")
  if (!all(req %in% colnames(df))) return(character(0))
  
  df %>%
    filter(!is.na(gene), gene != "",
           !is.na(p_val_adj), !is.na(avg_log2FC)) %>%
    filter(p_val_adj < qval, abs(avg_log2FC) > FC) %>%
    pull(gene) %>%
    unique()
}

# Map requested comparisons to bulk file keys.
bulk_key_for_comp <- function(A, B) {
  # AD vs Control
  if (A == "AD_AT" && B == "CONTROL") {
    return("TYPE_ADvsControl_gene_DEGs_FDRq0.05_logFC_0.25")
  }
  
  # LBD subtype vs Control (prefer redefined if present)
  if (A %in% c("LBD_S", "LBD_AS", "LBD_ATS") && B == "CONTROL") {
    return(paste0("LBD_redefined_", A, "_vs_Control_gene_DEGs_FDRq0.05_logFC_0.25"))
  }
  
  # LBD subtype vs AD (prefer redefined if present)
  if (A %in% c("LBD_S", "LBD_AS", "LBD_ATS") && B == "AD_AT") {
    return(paste0("LBD_redefined_", A, "_vs_AD_gene_DEGs_FDRq0.05_logFC_0.25"))
  }
  
  # If you *also* want to fall back to TYPE_* keys when redefined is missing, uncomment:
  # if (A %in% c("LBD_S","LBD_AS","LBD_ATS") && B == "CONTROL") return("TYPE_LBDvsControl_gene_DEGs_FDRq0.05_logFC_0.25")
  # if (A %in% c("LBD_S","LBD_AS","LBD_ATS") && B == "AD_AT")    return("TYPE_LBDvsAD_gene_DEGs_FDRq0.05_logFC_0.25")
  
  # No bulk analogue for subtype-vs-subtype (per your file list)
  return(NA_character_)
}

# Venn plot with counts + % labels for 2-set case using ggVennDiagram if available.
plot_venn_2set_counts_perc <- function(sets,
                                       title_text,
                                       set_labels = c("bulk RNAseq", "all cell types"),
                                       set_colors = c("bulk RNAseq" = "lightblue",
                                                      "all cell types" = "lightpink")) {
  
  stopifnot(is.list(sets), length(sets) == 2)
  names(sets) <- set_labels
  
  s_bulk <- unique(as.character(sets[[set_labels[1]]]))
  s_all  <- unique(as.character(sets[[set_labels[2]]]))
  
  n_bulk <- length(s_bulk)
  n_all  <- length(s_all)
  n_int  <- length(intersect(s_bulk, s_all))
  n_only_bulk <- n_bulk - n_int
  n_only_all  <- n_all  - n_int
  n_union <- length(union(s_bulk, s_all))
  
  lab_only_bulk <- sprintf("%d\n(%.1f%%)", n_only_bulk, 100 * n_only_bulk / n_union)
  lab_only_all  <- sprintf("%d\n(%.1f%%)", n_only_all,  100 * n_only_all  / n_union)
  lab_int       <- sprintf("%d\n(%.1f%%)", n_int,       100 * n_int       / n_union)
  
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    plot.new()
    title(main = "Install VennDiagram to render venns")
    return(invisible(NULL))
  }
  
  grid::grid.newpage()
  
  vd <- VennDiagram::draw.pairwise.venn(
    area1 = n_bulk,
    area2 = n_all,
    cross.area = n_int,
    category = set_labels,
    fill = c(set_colors[set_labels[1]], set_colors[set_labels[2]]),
    alpha = c(0.5, 0.5),
    cex = 0,           # hide default region counts
    cat.cex = 1.1,     # category label size
    ind = FALSE
  )
  
  grid::grid.draw(vd)
  
  # ---- Find where the category labels were actually drawn (left vs right) ----
  # vd is a grob list; category labels are text grobs with label == set label.
  get_label_x <- function(lbl) {
    idx <- which(vapply(vd, function(g) {
      inherits(g, "text") && !is.null(g$label) && identical(as.character(g$label), lbl)
    }, logical(1)))
    if (length(idx) == 0) return(NA_real_)
    as.numeric(vd[[idx[1]]]$x)
  }
  
  x_bulk <- get_label_x(set_labels[1])
  x_all  <- get_label_x(set_labels[2])
  
  # If we can't detect positions (rare), fall back to conventional left/right
  if (is.na(x_bulk) || is.na(x_all)) {
    x_left <- 0.35
    x_right <- 0.65
    bulk_is_left <- TRUE
  } else {
    x_left <- min(x_bulk, x_all)
    x_right <- max(x_bulk, x_all)
    bulk_is_left <- (x_bulk < x_all)
  }
  
  # Put the correct "only" label on the side where that set is drawn
  left_lab  <- if (bulk_is_left) lab_only_bulk else lab_only_all
  right_lab <- if (bulk_is_left) lab_only_all  else lab_only_bulk
  
  # Slightly inset from circle centers (helps readability)
  x_left_text  <- x_left  + 0.03
  x_right_text <- x_right - 0.03
  x_mid_text   <- (x_left + x_right) / 2
  
  grid::grid.text(left_lab,
                  x = grid::unit(x_left_text, "npc"),
                  y = grid::unit(0.50, "npc"),
                  gp = grid::gpar(fontsize = 11, fontface = "bold"))
  grid::grid.text(right_lab,
                  x = grid::unit(x_right_text, "npc"),
                  y = grid::unit(0.50, "npc"),
                  gp = grid::gpar(fontsize = 11, fontface = "bold"))
  grid::grid.text(lab_int,
                  x = grid::unit(x_mid_text, "npc"),
                  y = grid::unit(0.50, "npc"),
                  gp = grid::gpar(fontsize = 11, fontface = "bold"))
  
  grid::grid.text(title_text,
                  y = grid::unit(0.97, "npc"),
                  gp = grid::gpar(fontsize = 11, fontface = "bold"))
  
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# Main loop over comparisons
# ------------------------------------------------------------------------------

for (comp in comparisons) {
  
  A <- comp[[1]]
  B <- comp[[2]]
  comp_tag <- paste0(A, "_vs_", B)
  
  comp_dir <- file.path(out_dir, comp_tag)
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- get bulk genes for this comparison (if available) ----
  bulk_key <- bulk_key_for_comp(A, B)
  if (is.na(bulk_key) || !bulk_key %in% names(bulk_deg_tables)) {
    message("Skipping comparison (no matching bulk DEG table found): ", comp_tag,
            " ; bulk_key=", bulk_key %||% "NA")
    
    write.table(
      data.frame(comparison = comp_tag, reason = "no_bulk_table", bulk_key = bulk_key),
      file = file.path(comp_dir, "SKIPPED_noBulkTable.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    next
  }
  
  bulk_genes <- bulk_deg_tables[[bulk_key]]$gene %>%
    as.character() %>%
    unique()
  
  # ---- sn union across all cell types ----
  sn_union <- unique(unlist(lapply(cell_types, function(ct) {
    get_sn_deg_genes(deg_tables, ct, A, B, qval, FC)
  }), use.names = FALSE))
  
  # ---- shared gene list (bulk ∩ sn_union) ----
  shared_genes <- sort(intersect(bulk_genes, sn_union))
  
  shared_df <- data.frame(gene = shared_genes, stringsAsFactors = FALSE)
  write.table(shared_df,
              file = file.path(comp_dir, paste0("sharedGenes_bulk_vs_allCellTypes_", comp_tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # ---- Venn (one per comparison) ----
  venn_pdf <- file.path(comp_dir, paste0("Venn_bulkRNAseq_vs_allCellTypes_", comp_tag, ".pdf"))
  pdf(venn_pdf, width = 6.5, height = 6)
  
  # Order the sets: bulk first, union second (matches labels/colors below)
  sets <- list(
    bulk = bulk_genes,
    all_cells = sn_union
  )
  
  plot_venn_2set_counts_perc(
    sets = sets,
    title_text = paste0(comp_tag, ": bulk RNAseq vs all cell types (snRNA)"),
    set_labels = c("bulk RNAseq", "all cell types"),
    set_colors = c("bulk RNAseq" = "lightblue", "all cell types" = "lightpink")
  )
  dev.off()
  
  # ---- per cell type overlap summary ----
  overlap_df <- lapply(cell_types, function(ct) {
    sn_ct <- get_sn_deg_genes(deg_tables, ct, A, B, qval, FC)
    ov <- intersect(sn_ct, bulk_genes)
    
    n_sn   <- length(sn_ct)
    n_bulk <- length(bulk_genes)
    n_ov   <- length(ov)
    
    data.frame(
      cell_type = ct,
      n_sn = n_sn,
      n_bulk = n_bulk,
      n_overlap = n_ov,
      frac_sn_in_bulk = if (n_sn   > 0) n_ov / n_sn   else NA_real_,
      frac_bulk_in_sn = if (n_bulk > 0) n_ov / n_bulk else NA_real_,
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    mutate(
      cell_type = factor(cell_type, levels = cell_types),
      # reverse for plotting so neuron is on top
      cell_type_plot = factor(as.character(cell_type), levels = rev(cell_types))
    )
  
  write.table(overlap_df %>% select(-cell_type_plot),
              file = file.path(comp_dir, paste0("overlap_perCellType_", comp_tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # ---- common plotting theme ----
  base_theme <- theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(size = 9),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  # ---- Plot A: overlap / sn cell-type DEGs ----
  pA <- ggplot(overlap_df, aes(x = frac_sn_in_bulk, y = cell_type_plot, fill = cell_type)) +
    geom_col() +
    scale_fill_manual(values = color_panel) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste0(comp_tag, ": fraction of cell-type DEGs recapitulated in bulk"),
      x = "Overlap / sn cell-type DEGs",
      y = NULL
    ) +
    base_theme +
    theme(legend.position = "none")
  
  # ---- Plot B: overlap / bulk DEGs ----
  pB <- ggplot(overlap_df, aes(x = frac_bulk_in_sn, y = cell_type_plot, fill = cell_type)) +
    geom_col() +
    scale_fill_manual(values = color_panel) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste0(comp_tag, ": fraction of bulk DEGs captured by each cell type"),
      x = "Overlap / bulk DEGs",
      y = NULL
    ) +
    base_theme +
    theme(legend.position = "none")
  
  # ---- Plot C: total overlap counts ----
  pC <- ggplot(overlap_df, aes(x = n_overlap, y = cell_type_plot, fill = cell_type)) +
    geom_col() +
    scale_fill_manual(values = color_panel) +
    labs(
      title = paste0(comp_tag, ": total overlapping DEGs (count)"),
      x = "Number of overlapping DEGs",
      y = NULL
    ) +
    base_theme +
    theme(legend.position = "none")
  
  # Save all 3 plots into one PDF per comparison
  plot_pdf <- file.path(comp_dir, paste0("bulk_vs_sn_perCellType_overlapPlots_", comp_tag, ".pdf"))
  pdf(plot_pdf, width = 5, height = 3.0)
  print(pA)
  print(pB)
  print(pC)
  dev.off()
  
  # Save gene lists used (optional)
  write.table(sort(unique(bulk_genes)),
              file = file.path(comp_dir, paste0("bulk_DEGs_", comp_tag, "_geneList.txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(sort(unique(sn_union)),
              file = file.path(comp_dir, paste0("sn_allCellTypes_union_DEGs_", comp_tag, "_geneList.txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  message("Done: ", comp_tag,
          " | bulk_key=", bulk_key,
          " | shared=", length(shared_genes))
}
