#!/usr/bin/env Rscript
# ==============================================================================
# hdWGCNA module eigengene vs pathology plots (donor-level), split by sex
# - Facet by sex, boxplot by score, donor points, linear fit line (lm)
# - Spearman rho + p per sex at top-center
# - Hub genes bubble underneath (top 9, 3 per line), text white only if module color black
# ==============================================================================
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")
#!/usr/bin/env Rscript
# ==============================================================================
# hdWGCNA module eigengene vs pathology plots (DONOR-level), split by sex
# FIXES IN THIS VERSION:
#   - Cing.LB bins include 20 and axis ALWAYS shows 0,5,10,15,20,25
#   - Cing binning uses your count-based rule:
#       0 -> 0
#       1-5 -> 5
#       6-10 -> 10
#       11-15 -> 15
#       16-20 -> 20
#       21+ -> 25
#   - Linear fit is drawn with geom_abline() (no "hook")
#   - Edge bins not clipped: symmetric x padding for ALL traits
#   - Skip module == "grey"
#   - No "Hub genes" title (bubble only)
# ==============================================================================

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
  library(scales)
})

theme_set(cowplot::theme_cowplot(font_size = 8))

# ----------------------------- paths ------------------------------------------
obj_dir <- "../rObjects"
base_out_dir <- "../results/hdWGCNA"
out_dir <- file.path(base_out_dir, "ME_pathology_boxplot_by_sex")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

summary_file <- file.path(out_dir, "ME_pathology_rho_p_summary.tsv")

# ----------------------------- settings ---------------------------------------
cell_types <- c(
"astrocyte", "neuron", "interneuron", "opc", "microglia", "oligodendrocyte"
)

traits <- c("Braak.NFT", "Thal.amyloid", "Cing.LB")

braak_levels <- 0:6
thal_levels  <- 0:5
cing_bins    <- c(0, 5, 10, 15, 20, 25)   # KEEP 20

min_donors_total <- 8
min_donors_sex   <- 4
jitter_width_default <- 0.18

# x padding so first/last bins aren't clipped
pad_small <- 0.6
pad_cing  <- 2.0

# ---- label formatters ----
fmt_p_vec <- function(p) {
  ifelse(is.na(p), "p = NA",
         ifelse(p < 1e-3,
                paste0("p = ", format(p, scientific = TRUE, digits = 2)),
                paste0("p = ", signif(p, 3))))
}
fmt_rho_vec <- function(r) {
  ifelse(is.na(r), "rho = NA", paste0("rho = ", signif(r, 3)))
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

# ---- Cing.LB binning (YOUR RULE) ----
bin_cing <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  out <- rep(NA_real_, length(x))
  out[x == 0]                  <- 0
  out[x >= 1  & x <= 5]         <- 5
  out[x >= 6  & x <= 10]        <- 10
  out[x >= 11 & x <= 15]        <- 15
  out[x >= 16 & x <= 20]        <- 20
  out[x >= 21]                  <- 25
  out
}

# ---- straight-line fit per sex for geom_abline() ----
make_abline_df <- function(df, min_n) {
  out <- lapply(split(df, df$sex_inferred), function(d) {
    if (nrow(d) < min_n) return(NULL)
    fit <- tryCatch(lm(value ~ score_num, d), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    co <- coef(fit)
    data.frame(
      sex_inferred = unique(d$sex_inferred),
      intercept = unname(co[1]),
      slope = unname(co[2])
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out
}

# ----------------------------- summary collector ------------------------------
summary_rows <- list()

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
  if (!("sex_inferred" %in% colnames(md))) stop(ct, ": missing sex_inferred")

  md$sex_inferred <- factor(tolower(as.character(md$sex_inferred)), levels = c("female", "male"))
  md$Braak.NFT    <- suppressWarnings(as.numeric(as.character(md$Braak.NFT)))
  md$Thal.amyloid <- suppressWarnings(as.numeric(as.character(md$Thal.amyloid)))
  md$Cing.LB_bin  <- bin_cing(md$Cing.LB)

  # donor-level pathology table (stable)
  donor_path <- md %>%
    group_by(Sample_ID, sex_inferred) %>%
    summarise(
      Braak.NFT = { x <- Braak.NFT[!is.na(Braak.NFT)]; if (length(x)==0) NA_real_ else x[1] },
      Thal.amyloid = { x <- Thal.amyloid[!is.na(Thal.amyloid)]; if (length(x)==0) NA_real_ else x[1] },
      Cing.LB = { x <- Cing.LB_bin[!is.na(Cing.LB_bin)]; if (length(x)==0) NA_real_ else x[1] },
      .groups = "drop"
    )

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

  # long ME table (skip grey)
  ME_long <- pivot_longer(MEs, -cell_id, names_to = "module", values_to = "value") %>%
    filter(module != "grey") %>%
    inner_join(md[, c("cell_id", "Sample_ID", "sex_inferred")], by = "cell_id") %>%
    filter(!is.na(value))

  # donor-level ME per module
  donor_me <- ME_long %>%
    group_by(module, Sample_ID, sex_inferred) %>%
    summarise(value = mean(value), .groups = "drop")

  donor_df <- donor_me %>%
    left_join(donor_path, by = c("Sample_ID", "sex_inferred")) %>%
    left_join(mods_raw, by = "module")

  for (m in unique(donor_df$module)) {

    if (m == "grey") next

    df_m <- donor_df %>% filter(module == m)
    mod_col <- df_m$color[which(!is.na(df_m$color))[1]]
    if (is.na(mod_col) || length(mod_col) == 0) mod_col <- "grey70"

    # bubble text color only
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

    for (trait in traits) {

      score_levels <- switch(trait,
        "Braak.NFT" = braak_levels,
        "Thal.amyloid" = thal_levels,
        "Cing.LB" = cing_bins
      )

      df_t <- df_m %>%
        mutate(score_num = .data[[trait]]) %>%
        filter(!is.na(score_num)) %>%
        filter(score_num %in% score_levels)

      if (nrow(df_t) < min_donors_total) next
      if (length(unique(df_t$score_num)) < 2) next

      # spearman stats per sex
      stats <- df_t %>%
        group_by(sex_inferred) %>%
        summarise(
          n = n(),
          rho = if (n() >= min_donors_sex) suppressWarnings(cor(value, score_num, method = "spearman")) else NA_real_,
          p   = if (n() >= min_donors_sex) suppressWarnings(cor.test(value, score_num, method = "spearman", exact = FALSE)$p.value) else NA_real_,
          .groups = "drop"
        )

      summary_rows[[length(summary_rows) + 1]] <-
        stats %>%
        mutate(cell_type = ct, module = m, trait = trait, module_color = mod_col) %>%
        select(cell_type, module, module_color, trait, sex_inferred, n, rho, p)

      # label DF
      x_mid <- score_levels[ceiling(length(score_levels) / 2)]
      label_df <- stats %>%
        mutate(label = paste0(fmt_rho_vec(rho), "  ", fmt_p_vec(p)),
               x = x_mid, y = Inf)

      # straight-line fit for each sex
      abline_df <- make_abline_df(df_t, min_donors_sex)

      # low-N overlay (if bin has <3 donors)
      lowN_df <- df_t %>%
        group_by(sex_inferred, score_num) %>%
        summarise(n_bin = n(), med = median(value), .groups = "drop") %>%
        filter(n_bin < 3)

      # symmetric padding for ALL traits (prevents cutoff)
      pad <- if (trait == "Cing.LB") pad_cing else pad_small
      x_limits <- c(min(score_levels) - pad, max(score_levels) + pad)

      box_w <- if (trait == "Cing.LB") 2.4 else 0.7
      cb_w  <- if (trait == "Cing.LB") 1.8 else 0.5

      p_main <- ggplot(df_t, aes(x = score_num, y = value)) +
        geom_boxplot(
          aes(group = score_num),
          fill = mod_col,
          outlier.shape = NA,
          width = box_w
        ) +
        geom_crossbar(
          data = lowN_df,
          aes(x = score_num, y = med, ymin = med, ymax = med),
          inherit.aes = FALSE,
          width = cb_w,
          linewidth = 0.8,
          color = "black"
        ) +
        geom_point(
          position = position_jitter(width = jitter_width_default, height = 0),
          size = 1.25,
          alpha = 0.85,
          color = "black"
        ) +
        { if (!is.null(abline_df)) geom_abline(
          data = abline_df,
          aes(intercept = intercept, slope = slope),
          inherit.aes = FALSE,
          color = "darkgray",
          linewidth = 0.9
        ) } +
        facet_wrap(~ sex_inferred, nrow = 1) +
        geom_text(
          data = label_df,
          aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          vjust = 1.15,
          hjust = 0.5,
          size = 2.25
        ) +
        # CRITICAL: force breaks to include 20 for Cing.LB
        scale_x_continuous(
          breaks = score_levels,
          limits = x_limits,
          oob = scales::oob_squish
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
        labs(
          title = paste(m, trait),
          x = trait,
          y = "Module eigengene"
        ) +
        theme(
          plot.title = element_text(face = "plain", size = 8),
          plot.margin = margin(4, 6, 4, 6)
        ) +
        coord_cartesian(clip = "off")

      # bubble only (no "Hub genes" title)
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

      final <- plot_grid(
        p_main,
        p_bubble,
        ncol = 1,
        rel_heights = c(1, 0.22)
      )

      ggsave(
        filename = file.path(out_dir, paste0("ME_pathology_", m, "_", trait, "_by_sex.pdf")),
        plot = final,
        width = 3.5,
        height = 2.5,
        device = cairo_pdf
      )
    }
  }

  rm(seu)
  invisible(gc())
}

# ==============================================================================
# Write summary table
# ==============================================================================
if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows) %>%
    filter(module != "grey") %>%
    arrange(cell_type, module, trait, sex_inferred)

  write.table(summary_df, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote summary: ", summary_file)
} else {
  warning("No summary rows collected; summary file not written.")
}

message("DONE.")
