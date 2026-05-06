#!/usr/bin/env Rscript
# ==============================================================================
# hdWGCNA module eigengene vs DISEASE GROUP plots (donor-level), BOTH SEXES TOGETHER
# - X axis = disease group (metadata column "group")
# - Boxplot by group + donor points
# - Kruskal-Wallis p-value shown at top-center
# - Pairwise posthoc within module: Dunn test (BH-adjusted) written to TSV
# - Hub genes bubble underneath (top 9, 3 per line)
# - Skip module == "grey"
#
# Output size/style: width=2.2, height=2.5, font_size=8
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
})


theme_set(cowplot::theme_cowplot(font_size = 8))

# ----------------------------- paths ------------------------------------------
obj_dir      <- "../rObjects"
base_out_dir <- "../results/hdWGCNA"

out_dir <- file.path(base_out_dir, "ME_boxplot_disease")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

summary_kw_file   <- file.path(out_dir, "ME_disease_kw_p_summary.tsv")
summary_pair_file <- file.path(out_dir, "ME_disease_pairwise_dunn_BH.tsv")

# ----------------------------- settings ---------------------------------------
cell_types <- c("astrocyte", "neuron", "interneuron", "opc", "microglia", "oligodendrocyte")

min_donors_total     <- 8
jitter_width_default <- 0.18

# Group order + labels
group_levels_raw <- c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS")
group_labels_map <- c(
  CONTROL = "Control",
  AD_AT   = "AD (AT)",
  LBD_S   = "LBD (S)",
  LBD_AS  = "LBD (AS)",
  LBD_ATS = "LBD (ATS)"
)
group_levels_label <- unname(group_labels_map[group_levels_raw])

# ---- label formatter ----
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 1e-3) return(paste0("p = ", format(p, scientific = TRUE, digits = 2)))
  paste0("p = ", signif(p, 3))
}

wrap_genes_3_per_line <- function(genes) {
  if (length(genes) == 0) return("(hub genes not found)")
  idx <- ceiling(seq_along(genes) / 3)
  paste(tapply(genes, idx, paste, collapse = "  "), collapse = "\n")
}

is_blackish <- function(col) {
  if (is.na(col) || length(col) == 0) return(FALSE)
  if (tolower(col) == "black") return(TRUE)
  rgb <- tryCatch(col2rgb(col), error = function(e) NULL)
  if (is.null(rgb)) return(FALSE)
  (0.2126 * rgb[1] + 0.7152 * rgb[2] + 0.0722 * rgb[3]) / 255 < 0.08
}

# ----------------------------- collectors -------------------------------------
kw_rows       <- list()
pairwise_rows <- list()

# helper: safe dunn posthoc (BH)
run_dunn_bh <- function(df, value_col = "value", group_col = "group_label") {
  # df must have columns: value_col numeric, group_col factor
  df2 <- df %>% filter(!is.na(.data[[value_col]]), !is.na(.data[[group_col]]))
  df2[[group_col]] <- droplevels(df2[[group_col]])
  if (nlevels(df2[[group_col]]) < 2) return(NULL)

  # FSA::dunnTest returns list with $res: Comparison, Z, P.unadj, P.adj
  out <- tryCatch(
    FSA::dunnTest(stats::as.formula(paste0(value_col, " ~ ", group_col)),
                  data = df2, method = "bh"),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)

  res <- out$res
  # Normalize comparison formatting "A - B"
  res <- res %>%
    as.data.frame() %>%
    mutate(
      Comparison = as.character(Comparison),
      group1 = sub("\\s*-\\s*.*$", "", Comparison),
      group2 = sub("^.*\\s*-\\s*", "", Comparison)
    ) %>%
    select(group1, group2, Z, P.unadj, P.adj)

  res
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
  md  <- seu[[]]
  md$cell_id <- rownames(md)

  if (!("Sample_ID" %in% colnames(md))) stop(ct, ": missing Sample_ID")
  if (!("group" %in% colnames(md)))     stop(ct, ": missing group")

  # donor-level group table (stable; take first non-NA per donor)
  donor_group <- md %>%
    group_by(Sample_ID) %>%
    summarise(
      group = { x <- group[!is.na(group)]; if (length(x) == 0) NA_character_ else as.character(x[1]) },
      .groups = "drop"
    )
donor_group
  # module colors (skip grey)
  mods_raw <- GetModules(seu) %>%
    distinct(module, color) %>%
    filter(module != "grey")

  # MEs (harmonized)
  MEs <- as.data.frame(GetMEs(seu, harmonized = TRUE))
  MEs$cell_id <- rownames(MEs)

  # hubs
  hub_file <- file.path(base_out_dir, paste0("hubs_", ct, ".tsv"))
  hubs_tbl <- if (file.exists(hub_file)) read.delim(hub_file, stringsAsFactors = FALSE) else NULL

  # long ME table (skip grey) + join Sample_ID for donor aggregation
  ME_long <- pivot_longer(MEs, -cell_id, names_to = "module", values_to = "value") %>%
    filter(module != "grey") %>%
    inner_join(md[, c("cell_id", "Sample_ID")], by = "cell_id") %>%
    filter(!is.na(value))

  # donor-level ME per module
  donor_me <- ME_long %>%
    group_by(module, Sample_ID) %>%
    summarise(value = mean(value), .groups = "drop")

  donor_df <- donor_me %>%
    left_join(donor_group, by = "Sample_ID") %>%
    left_join(mods_raw, by = "module")

  for (m in unique(donor_df$module)) {

    if (m == "grey") next
    df_m <- donor_df %>% filter(module == m)

    mod_col <- df_m$color[which(!is.na(df_m$color))[1]]
    if (is.na(mod_col) || length(mod_col) == 0) mod_col <- "grey70"

    txt_col <- if (is_blackish(mod_col)) "white" else "black"

    # hub bubble text (top 9, 3 per line)
    hubs_text <- "(hub genes not found)"
    if (!is.null(hubs_tbl) && all(c("module", "gene_name") %in% colnames(hubs_tbl))) {
      hubs <- hubs_tbl %>%
        filter(module == m) %>%
        pull(gene_name) %>%
        unique() %>%
        head(9)
      if (length(hubs) > 0) hubs_text <- wrap_genes_3_per_line(hubs)
    }

    # prep group factor + labels
    df_g <- df_m %>%
      mutate(group_raw = as.character(group)) %>%
      filter(!is.na(group_raw)) %>%
      filter(group_raw %in% group_levels_raw) %>%
      mutate(
        group_label = factor(group_labels_map[group_raw], levels = group_levels_label)
      ) %>%
      filter(!is.na(group_label))

    if (nrow(df_g) < min_donors_total) next
    if (nlevels(droplevels(df_g$group_label)) < 2) next

    # Kruskal-Wallis p
    kw_p <- tryCatch(
      kruskal.test(value ~ group_label, data = df_g)$p.value,
      error = function(e) NA_real_
    )

    kw_rows[[length(kw_rows) + 1]] <- data.frame(
      cell_type = ct,
      module = m,
      module_color = mod_col,
      n = nrow(df_g),
      kw_p = kw_p,
      stringsAsFactors = FALSE
    )

    # Pairwise Dunn BH
    dunn_res <- run_dunn_bh(df_g, value_col = "value", group_col = "group_label")
    if (!is.null(dunn_res) && nrow(dunn_res) > 0) {
      pairwise_rows[[length(pairwise_rows) + 1]] <- dunn_res %>%
        mutate(
          cell_type = ct,
          module = m,
          module_color = mod_col,
          n = nrow(df_g)
        ) %>%
        select(cell_type, module, module_color, n, group1, group2, Z, P.unadj, P.adj)
    }

    # label at top-center (middle category)
    levs <- levels(df_g$group_label)
    x_mid <- levs[ceiling(length(levs) / 2)]
    label_df <- data.frame(x = x_mid, y = Inf, label = fmt_p(kw_p))

    # low-N overlay (if group has <3 donors)
    lowN_df <- df_g %>%
      group_by(group_label) %>%
      summarise(n_bin = n(), med = median(value), .groups = "drop") %>%
      filter(n_bin < 3)

    p_main <- ggplot(df_g, aes(x = group_label, y = value)) +
      geom_boxplot(fill = mod_col, outlier.shape = NA, width = 0.65) +
      geom_crossbar(
        data = lowN_df,
        aes(x = group_label, y = med, ymin = med, ymax = med),
        inherit.aes = FALSE,
        width = 0.50,
        linewidth = 0.8,
        color = "black"
      ) +
      geom_point(
        position = position_jitter(width = jitter_width_default, height = 0),
        size = 1.25,
        alpha = 0.85,
        color = "black"
      ) +
      geom_text(
        data = label_df,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        vjust = 1.15,
        hjust = 0.5,
        size = 2.25
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
      labs(
        title = paste(m),
        x = NULL,
        y = "Module eigengene"
      ) +
      theme(
        plot.title = element_text(face = "plain", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(4, 6, 4, 6)
      ) +
      coord_cartesian(clip = "off")

    p_bubble <- ggplot(
      data.frame(x = 1, y = 1, label = hubs_text),
      aes(x, y, label = label)
    ) +
      geom_label(
        fill = mod_col,
        color = txt_col,
        label.size = 0,
        size = 2.25,
        lineheight = 0.95,
        label.padding = unit(0.45, "lines"),
        label.r = unit(0.8, "lines")
      ) +
      theme_void() +
      theme(plot.margin = margin(0, 10, 1, 10))

    final <- cowplot::plot_grid(
      p_main, p_bubble, ncol = 1, rel_heights = c(1, 0.22)
    )

    ggsave(
      filename = file.path(out_dir, paste0("ME_group_", m, ".pdf")),
      plot = final,
      width = 2.2,
      height = 2.5,
      device = cairo_pdf
    )
  }

  rm(seu)
  invisible(gc())
}

# ==============================================================================
# Write summary tables
# ==============================================================================
if (length(kw_rows) > 0) {
  kw_df <- bind_rows(kw_rows) %>%
    filter(module != "grey") %>%
    arrange(cell_type, module)
  write.table(kw_df, summary_kw_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote KW summary: ", summary_kw_file)
} else {
  warning("No KW rows collected; KW summary file not written.")
}

if (length(pairwise_rows) > 0) {
  pair_df <- bind_rows(pairwise_rows) %>%
    filter(module != "grey") %>%
    arrange(cell_type, module, P.adj)
  write.table(pair_df, summary_pair_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote pairwise summary: ", summary_pair_file)
} else {
  warning("No pairwise rows collected; pairwise summary file not written.")
}

message("DONE.")
