#!/usr/bin/env Rscript
# ==============================================================================
# hdWGCNA GO Enrichr barplots (loop all cell types)
# - Reads:  ../rObjects/hdWGCNA_<celltype>_final.qs
# - Runs:   RunEnrichr() for GO BP/CC/MF 2023
# - Saves:  ../results/hdWGCNA_GO/enrichr_table_<celltype>.tsv
# - Style:  text size 8 (global ggplot base size)
# ==============================================================================

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
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(FSA)
  library(enrichR)
  library(tibble)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------
wgcna_name <- "CWOW"
obj_dir    <- "../rObjects"
out_base   <- "../results/hdWGCNA_GO"
cell_types <- c("neuron", "interneuron", "opc") # "astrocyte", "microglia", "oligodendrocyte", 

dbs <- c(
  "GO_Biological_Process_2023",
  "GO_Cellular_Component_2023",
  "GO_Molecular_Function_2023"
)

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 8))

# ------------------------------------------------------------------------------
# Helper: coalesce across columns that exist (NO tidy-eval)
# ------------------------------------------------------------------------------
coalesce_cols <- function(df, cols) {
  cols <- cols[cols %in% colnames(df)]
  if (length(cols) == 0) return(rep(NA_character_, nrow(df)))

  vecs <- lapply(cols, function(cn) df[[cn]])
  # Reduce coalesce across vectors
  out <- vecs[[1]]
  if (length(vecs) > 1) {
    for (i in 2:length(vecs)) out <- dplyr::coalesce(out, vecs[[i]])
  }
  out
}

# ------------------------------------------------------------------------------
# Helper: add module colors + clean columns (robust module handling)
# ------------------------------------------------------------------------------
add_module_colors_and_clean <- function(enrich_df, dataObject.ct, wgcna_name = "CWOW") {

  # ---- module -> color map from object ----
  modules_df <- dataObject.ct@misc[[wgcna_name]]$wgcna_modules %>%
    as.data.frame() %>%
    as_tibble() %>%
    dplyr::select(module, color) %>%
    distinct() %>%
    mutate(module = gsub("\u2013", "-", module))  # normalize en-dash

  if (!all(c("module", "color") %in% colnames(modules_df))) {
    stop("wgcna_modules is missing required columns: 'module' and/or 'color'.")
  }

  enrich_tbl <- enrich_df %>%
    as.data.frame() %>%
    as_tibble()

  # Drop unwanted columns if present
  enrich_tbl <- enrich_tbl %>%
    select(-any_of(c("Old.P.value", "Old.Adjusted.P.value")))

  # ---- Find module-like columns in Enrichr output ----
  module_candidates <- c(
    "module", "Module", "module_name", "ModuleName", "modules", "Modules",
    "module.x", "module.y"
  )
  present_module_cols <- intersect(module_candidates, colnames(enrich_tbl))

  if (length(present_module_cols) == 0) {
    warning("No recognizable module column found in Enrichr table; returning table without colors.")
    enrich_tbl$module <- NA_character_
    enrich_tbl$color  <- NA_character_
    return(enrich_tbl)
  }

  # ---- Create canonical module vector safely (outside mutate) ----
  module_vec <- coalesce_cols(enrich_tbl, present_module_cols)
  module_vec <- gsub("\u2013", "-", module_vec) # normalize en-dash
  module_vec <- gsub("–", "-", module_vec)      # normalize literal en-dash just in case

  # Add canonical module column (overwrite if already exists)
  enrich_tbl$module <- module_vec

  # Drop any extra module columns (keep canonical "module" only)
  drop_mod_cols <- setdiff(present_module_cols, "module")
  if (length(drop_mod_cols) > 0) {
    enrich_tbl <- enrich_tbl %>% select(-any_of(drop_mod_cols))
  }

  # ---- Join module -> color ----
  enrich_tbl <- enrich_tbl %>% left_join(modules_df, by = "module")

  # If color.x/color.y happen, collapse them safely
  if (any(c("color.x", "color.y") %in% colnames(enrich_tbl))) {
    enrich_tbl$color <- coalesce_cols(enrich_tbl, c("color", "color.y", "color.x"))
    enrich_tbl <- enrich_tbl %>% select(-any_of(c("color.x", "color.y")))
  }

  enrich_tbl
}

# ------------------------------------------------------------------------------
# Loop cell types
# ------------------------------------------------------------------------------
for (ct in cell_types) {

  message("------------------------------------------------------------")
  message("Cell type: ", ct)

  seu_path <- file.path(obj_dir, paste0("hdWGCNA_", ct, "_final.qs"))
  if (!file.exists(seu_path)) {
    warning("Missing object: ", seu_path)
    next
  }

  dataObject.ct <- qread(seu_path)

  dataObject.ct <- RunEnrichr(
    dataObject.ct,
    dbs = dbs,
    max_genes = 100
  )

  enrich_df <- GetEnrichrTable(dataObject.ct)
  enrich_df <- add_module_colors_and_clean(enrich_df, dataObject.ct, wgcna_name = wgcna_name)

  write.table(
    enrich_df,
    file = file.path(out_base, paste0("enrichr_table_", ct, ".tsv")),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

  EnrichrBarPlot(
    dataObject.ct,
    outdir    = out_base,
    n_terms   = 10,
    plot_size = c(6, 7),
    logscale  = TRUE
  )

  message("Saved table + plots to: ", out_base)
}

message("Done.")
