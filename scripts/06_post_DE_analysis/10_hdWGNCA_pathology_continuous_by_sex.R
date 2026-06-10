#!/usr/bin/env Rscript
# ==============================================================================
# hdWGCNA module eigengene vs pathology plots (donor-level), split by sex
# - Continuous scatterplots for Braak.NFT, Thal.amyloid, and Cing.LB
# - NO boxplots
# - NO Cing.LB binning
# - NO hub gene bubbles
# - Donor points + straight linear fit line per sex
# - Spearman rho + p per sex
# - Points filled by module color, outlined dark gray
# - Female = circle, male = square
# - Skip module == "grey"
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
  library(scales)
})

theme_set(cowplot::theme_cowplot(font_size = 8))

# ----------------------------- paths ------------------------------------------
obj_dir <- "../rObjects"
base_out_dir <- "../results/hdWGCNA"

out_dir <- file.path(base_out_dir, "ME_pathology_continuous_by_sex")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

summary_file <- file.path(out_dir, "ME_pathology_rho_p_summary.tsv")

# ----------------------------- settings ---------------------------------------
cell_types <- c(
  "astrocyte", "neuron", "interneuron", "opc", "microglia", "oligodendrocyte"
)

traits <- c("Braak.NFT", "Thal.amyloid", "Cing.LB")

min_donors_total <- 8
min_donors_sex   <- 4
jitter_width_default <- 0.05

fmt_p_vec <- function(p) {
  ifelse(
    is.na(p), "p = NA",
    ifelse(
      p < 1e-3,
      paste0("p = ", format(p, scientific = TRUE, digits = 2)),
      paste0("p = ", signif(p, 3))
    )
  )
}

fmt_rho_vec <- function(r) {
  ifelse(is.na(r), "rho = NA", paste0("rho = ", signif(r, 3)))
}

make_abline_df <- function(df, min_n) {
  out <- lapply(split(df, df$sex_inferred), function(d) {
    if (nrow(d) < min_n) return(NULL)
    if (length(unique(d$score_num)) < 2) return(NULL)

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

trait_x_label <- function(trait) {
  if (trait == "Cing.LB") return("Lewy bodies per 0.25 mm²")
  trait
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

  if (!("Sample_ID" %in% colnames(md))) stop(ct, ": missing Sample_ID")
  if (!("sex_inferred" %in% colnames(md))) stop(ct, ": missing sex_inferred")

  md$sex_inferred <- factor(
    tolower(as.character(md$sex_inferred)),
    levels = c("female", "male")
  )

  md$Braak.NFT    <- suppressWarnings(as.numeric(as.character(md$Braak.NFT)))
  md$Thal.amyloid <- suppressWarnings(as.numeric(as.character(md$Thal.amyloid)))
  md$Cing.LB      <- suppressWarnings(as.numeric(as.character(md$Cing.LB)))

  donor_path <- md %>%
    group_by(Sample_ID, sex_inferred) %>%
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
        x <- Cing.LB[!is.na(Cing.LB)]
        if (length(x) == 0) NA_real_ else x[1]
      },
      .groups = "drop"
    )

  mods_raw <- GetModules(seu) %>%
    distinct(module, color) %>%
    filter(module != "grey")

  MEs <- as.data.frame(GetMEs(seu, harmonized = TRUE))
  MEs$cell_id <- rownames(MEs)

  ME_long <- pivot_longer(MEs, -cell_id, names_to = "module", values_to = "value") %>%
    filter(module != "grey") %>%
    inner_join(
      md[, c("cell_id", "Sample_ID", "sex_inferred")],
      by = "cell_id"
    ) %>%
    filter(!is.na(value), !is.na(sex_inferred))

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

    for (trait in traits) {

      df_t <- df_m %>%
        mutate(score_num = .data[[trait]]) %>%
        filter(!is.na(score_num), !is.na(value), !is.na(sex_inferred))

      if (nrow(df_t) < min_donors_total) next
      if (length(unique(df_t$score_num)) < 2) next

      stats <- df_t %>%
        group_by(sex_inferred) %>%
        summarise(
          n = n(),
          rho = if (n() >= min_donors_sex && length(unique(score_num)) >= 2) {
            suppressWarnings(cor(value, score_num, method = "spearman"))
          } else {
            NA_real_
          },
          p = if (n() >= min_donors_sex && length(unique(score_num)) >= 2) {
            suppressWarnings(
              cor.test(value, score_num, method = "spearman", exact = FALSE)$p.value
            )
          } else {
            NA_real_
          },
          .groups = "drop"
        )

      summary_rows[[length(summary_rows) + 1]] <- stats %>%
        mutate(
          cell_type = ct,
          module = m,
          module_color = mod_col,
          trait = trait
        ) %>%
        select(cell_type, module, module_color, trait, sex_inferred, n, rho, p)

      x_mid <- mean(range(df_t$score_num, na.rm = TRUE))

      label_df <- stats %>%
        mutate(
          x = x_mid,
          y = Inf,
          label = paste0(fmt_rho_vec(rho), "\n", fmt_p_vec(p))
        )

      abline_df <- make_abline_df(df_t, min_donors_sex)

      x_range <- range(df_t$score_num, na.rm = TRUE)
      x_pad <- diff(x_range) * 0.08
      if (x_pad == 0) x_pad <- 0.5

      if (trait == "Braak.NFT") {
        x_breaks <- 0:6
      } else if (trait == "Thal.amyloid") {
        x_breaks <- 0:5
      } else {
        x_breaks <- pretty(df_t$score_num)
      }

      p_main <- ggplot(df_t, aes(x = score_num, y = value)) +
        geom_point(
          aes(shape = sex_inferred),
          position = position_jitter(width = jitter_width_default, height = 0),
          size = 1.8,
          stroke = 0.35,
          fill = mod_col,
          color = "grey20",
          alpha = 0.95
        ) +
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
        facet_wrap(~ sex_inferred, nrow = 1) +
        geom_text(
          data = label_df,
          aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          vjust = 1.5,
          hjust = 0.5,
          size = 2,
          lineheight = 0.8
        ) +
        scale_shape_manual(
          values = c(
            female = 21,
            male = 22
          ),
          drop = FALSE
        ) +
        scale_x_continuous(
          breaks = x_breaks,
          limits = c(x_range[1] - x_pad, x_range[2] + x_pad),
          oob = scales::oob_squish
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
        labs(
          title = paste(m, trait),
          x = trait_x_label(trait),
          y = "Module eigengene"
        ) +
        theme(
          legend.position = "none",
          plot.title = element_text(face = "plain", size = 8),
          plot.margin = margin(4, 6, 4, 6)
        ) +
        coord_cartesian(clip = "off")

      ggsave(
        filename = file.path(out_dir, paste0("ME_pathology_", m, "_", trait, "_by_sex.pdf")),
        plot = p_main,
        width = 3.5,
        height = 1.75,
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