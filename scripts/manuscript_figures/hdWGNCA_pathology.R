#!/usr/bin/env Rscript

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
  library(scales)
})

theme_set(cowplot::theme_cowplot(font_size = 7))

# ----------------------------- paths ------------------------------------------
obj_dir <- "../rObjects"
base_out_dir <- "../results/hdWGCNA"

out_dir <- file.path(base_out_dir, "ME_pathology_no_hub_bubble")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

summary_file <- file.path(out_dir, "ME_pathology_rho_p_summary.tsv")

# ----------------------------- settings ---------------------------------------
cell_types <- c(
  "astrocyte"
)

traits <- c("Braak.NFT", "Thal.amyloid", "Cing.LB")

braak_levels <- 0:6
thal_levels  <- 0:5
cing_bins    <- c(0, 5, 10, 15, 20, 25)

min_donors_total <- 8
jitter_width_default <- 0.18
jitter_width_cing <- 0.55

pad_small <- 0.6
pad_cing  <- 2.0

# ----------------------------- helper functions -------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 1e-3) return(paste0("p = ", format(p, scientific = TRUE, digits = 2)))
  paste0("p = ", signif(p, 3))
}

fmt_rho <- function(r) {
  if (is.na(r)) return("rho = NA")
  paste0("rho = ", signif(r, 3))
}

bin_cing <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  out <- rep(NA_real_, length(x))
  out[x == 0]                  <- 0
  out[x >= 1  & x <= 5]        <- 5
  out[x >= 6  & x <= 10]       <- 10
  out[x >= 11 & x <= 15]       <- 15
  out[x >= 16 & x <= 20]       <- 20
  out[x >= 21]                 <- 25
  out
}

make_abline_df_all <- function(df, min_n) {
  if (nrow(df) < min_n) return(NULL)

  fit <- tryCatch(
    lm(value ~ score_num, data = df),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  co <- coef(fit)

  data.frame(
    intercept = unname(co[1]),
    slope = unname(co[2])
  )
}

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

  if (!("Sample_ID" %in% colnames(md))) {
    stop(ct, ": missing Sample_ID")
  }

  md$Braak.NFT    <- suppressWarnings(as.numeric(as.character(md$Braak.NFT)))
  md$Thal.amyloid <- suppressWarnings(as.numeric(as.character(md$Thal.amyloid)))
  md$Cing.LB_bin  <- bin_cing(md$Cing.LB)

  donor_path <- md %>%
    group_by(Sample_ID) %>%
    summarise(
      Braak.NFT = {
        x <- Braak.NFT[!is.na(Braak.NFT)]
        if (length(x) == 0) NA_real_ else x[1]
      },
      Thal.amyloid = {
        x <- Thal.amyloid[!is.na(Thal.amyloid)]
        if (length(x) == 0) NA_real_ else x[1]
      },
      Cing.LB = {
        x <- Cing.LB_bin[!is.na(Cing.LB_bin)]
        if (length(x) == 0) NA_real_ else x[1]
      },
      .groups = "drop"
    )

  mods_raw <- GetModules(seu) %>%
    distinct(module, color) %>%
    filter(module != "grey")

  MEs <- as.data.frame(GetMEs(seu, harmonized = TRUE))
  MEs$cell_id <- rownames(MEs)

  ME_long <- pivot_longer(
    MEs,
    cols = -cell_id,
    names_to = "module",
    values_to = "value"
  ) %>%
    filter(module != "grey") %>%
    inner_join(md[, c("cell_id", "Sample_ID")], by = "cell_id") %>%
    filter(!is.na(value))

  donor_me <- ME_long %>%
    group_by(module, Sample_ID) %>%
    summarise(value = mean(value), .groups = "drop")

  donor_df <- donor_me %>%
    left_join(donor_path, by = "Sample_ID") %>%
    left_join(mods_raw, by = "module")

  for (m in unique(donor_df$module)) {

    if (m == "grey") next

    df_m <- donor_df %>%
      filter(module == m)

    mod_col <- df_m$color[which(!is.na(df_m$color))[1]]
    if (is.na(mod_col) || length(mod_col) == 0) mod_col <- "grey70"

    for (trait in traits) {

      score_levels <- switch(
        trait,
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

      rho <- suppressWarnings(
        cor(df_t$value, df_t$score_num, method = "spearman")
      )

      p <- suppressWarnings(
        cor.test(
          df_t$value,
          df_t$score_num,
          method = "spearman",
          exact = FALSE
        )$p.value
      )

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        cell_type = ct,
        module = m,
        module_color = mod_col,
        trait = trait,
        n = nrow(df_t),
        rho = rho,
        p = p,
        stringsAsFactors = FALSE
      )

      x_mid <- score_levels[ceiling(length(score_levels) / 2)]

      label_df <- data.frame(
        x = x_mid,
        y = Inf,
        label = paste0(fmt_rho(rho), "  ", fmt_p(p))
      )

      abline_df <- make_abline_df_all(df_t, min_donors_total)

      lowN_df <- df_t %>%
        group_by(score_num) %>%
        summarise(
          n_bin = n(),
          med = median(value),
          .groups = "drop"
        ) %>%
        filter(n_bin < 3)

      pad <- if (trait == "Cing.LB") pad_cing else pad_small

      x_limits <- c(
        min(score_levels) - pad,
        max(score_levels) + pad
      )

      x_axis_label <- if (trait == "Cing.LB") {
        expression("Lewy bodies per 0.25 mm"^2)
      } else {
        trait
      }

      if (trait == "Cing.LB") {

        p_main <- ggplot(df_t, aes(x = score_num, y = value)) +
          geom_point(
            shape = 21,
            fill = mod_col,
            colour = "black",
            stroke = 0.7,
            size = 1.8,
            alpha = 1,
            position = position_jitter(width = jitter_width_cing, height = 0)
          )

      } else {

        p_main <- ggplot(df_t, aes(x = score_num, y = value)) +
          geom_boxplot(
            aes(group = score_num),
            fill = mod_col,
            color = "black",
            outlier.shape = NA,
            width = 0.7
          ) +
          geom_crossbar(
            data = lowN_df,
            aes(x = score_num, y = med, ymin = med, ymax = med),
            inherit.aes = FALSE,
            width = 0.5,
            linewidth = 0.8,
            color = "black"
          ) +
          geom_point(
            shape = 21,
            fill = mod_col,
            colour = "black",
            stroke = 0.7,
            size = 1.6,
            alpha = 1,
            position = position_jitter(width = jitter_width_default, height = 0)
          )
      }

      p_main <- p_main +
        {
          if (!is.null(abline_df)) {
            geom_abline(
              data = abline_df,
              aes(intercept = intercept, slope = slope),
              inherit.aes = FALSE,
              color = "darkgray",
              linewidth = 0.9
            )
          }
        } +
        geom_text(
          data = label_df,
          aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          vjust = 1.15,
          hjust = 0.5,
          size = 2.25
        ) +
        scale_x_continuous(
          breaks = score_levels,
          limits = x_limits,
          oob = scales::oob_squish
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
        labs(
          title = paste(m, trait),
          x = x_axis_label,
          y = "Module eigengene"
        ) +
        theme(
          plot.title = element_text(face = "plain", size = 7),
          plot.margin = margin(4, 6, 4, 6)
        ) +
        coord_cartesian(clip = "off")

      ggsave(
        filename = file.path(out_dir, paste0("ME_pathology_", m, "_", trait, ".pdf")),
        plot = p_main,
        width = 2.25,
        height = 1.85,
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
    arrange(cell_type, module, trait)

  write.table(
    summary_df,
    summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("Wrote summary: ", summary_file)

} else {
  warning("No summary rows collected; summary file not written.")
}

message("DONE.")