#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(readr)
  library(ragg)
  library(scales)
  library(stringr)
  library(svglite)
})

repo_root <- normalizePath(Sys.getenv("NHANES_REPO", "."), mustWork = TRUE)
out_dir <- file.path(repo_root, "results", "exwas", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

description_file <- file.path(repo_root, "data", "data_dictionary", "pollutant_description.csv")
cvd_file <- file.path(repo_root, "results", "exwas", "survey_weighted_exwas_CVD.csv")
le8_file <- file.path(repo_root, "results", "exwas", "survey_weighted_exwas_LE8_primary.csv")
stopifnot(file.exists(description_file), file.exists(cvd_file), file.exists(le8_file))

description <- read_csv(description_file, show_col_types = FALSE) %>%
  select(exposure, family, matrix)

# Reproduce the original manuscript mapping exactly. In the source ExposomeSet,
# Set1 and the added colours were assigned in familyNames(exp) order. Keep the
# mapping explicit so filtering or reordering the description table cannot shift
# a colour onto a different toxicant family.
original_family_order <- c(
  "Aromatic Amines",
  "Heavy Metals",
  "Volatile Organic Compounds (VOC)",
  "Brominated Flame Retardants",
  "Personal Care Products",
  "Pesticides",
  "Gaseous Pollutants",
  "Herbicides",
  "Phthalates",
  "Polyfluoroalkyl Substances (PFAS)",
  "Inorganic Anions",
  "Polychlorinated Biphenyls",
  "Polycyclic Aromatic Hydrocarbons (PAH)"
)
original_family_colors <- c(
  brewer.pal(9, "Set1"),
  "#FB9A99", "#66C2A5", "#FFD92F", "#A6D854", "#E78AC3"
)[seq_along(original_family_order)]
family_colors <- setNames(original_family_colors, original_family_order)
# Use the warm yellow sampled from the author-provided PAH colour reference.
family_colors["Polycyclic Aromatic Hydrocarbons (PAH)"] <- "#F9D957"
family_order <- original_family_order[original_family_order %in% unique(description$family)]

prepare_results <- function(path, outcome) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(method == "survey_weighted", status == "ok") %>%
    left_join(description, by = "exposure") %>%
    mutate(
      outcome = outcome,
      family = factor(family, levels = family_order),
      nominal = p_value < 0.05,
      fdr_significant = p_fdr_bh < 0.05,
      # R's syntactic-name conversion can prepend X to names beginning with a
      # digit. Remove only that artificial leading X for figure labels while
      # preserving the original exposure key used in the analysis.
      display_label = sub("^X([0-9])", "\\1", exposure)
    )
}

cvd <- prepare_results(cvd_file, "CVD") %>%
  mutate(label_point = nominal & estimate > 1)

le8 <- prepare_results(le8_file, "LE8") %>%
  arrange(p_value) %>%
  mutate(
    significance_rank = row_number(),
    label_point = fdr_significant & estimate < 0
  )

# On a raw-P axis, show the largest P value that still passes BH-FDR < 0.05.
le8_bh_p_cutoff <- max(le8$p_value[le8$fdr_significant], na.rm = TRUE)

stopifnot(nrow(cvd) == 149, nrow(le8) == 149)
stopifnot(!anyNA(cvd$family), !anyNA(le8$family))

p_axis <- scale_y_continuous(
  trans = scales::transform_compose("log10", "reverse"),
  breaks = c(1, 0.5, 0.1, 0.05, 0.01, 0.001, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14),
  labels = function(x) vapply(x, function(value) {
    if (is.na(value)) return(NA_character_)
    if (value >= 0.1) return(sub("\\.0$", "", sprintf("%.1f", value)))
    if (value >= 0.01) return(sub("0$", "", sprintf("%.2f", value)))
    if (value >= 0.001) return(sprintf("%.3f", value))
    format(value, scientific = TRUE, digits = 1)
  }, character(1)),
  name = "P value"
)

theme_original_weighted <- function() {
  theme_classic(base_size = 10, base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9, color = "black"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 9.5, color = "#4D4D4D", hjust = 0),
      legend.title = element_text(size = 9.5),
      legend.text = element_text(size = 8.5),
      legend.key.height = grid::unit(4.2, "mm"),
      legend.position = "right",
      plot.margin = margin(8, 10, 8, 8)
    )
}

common_layers <- list(
  scale_color_manual(values = family_colors, drop = FALSE),
  p_axis,
  guides(color = guide_legend(title = "Toxicant family", override.aes = list(size = 3, alpha = 1))),
  theme_original_weighted()
)

figure1 <- ggplot(cvd, aes(x = estimate, y = p_value, color = family)) +
  common_layers +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray45", linewidth = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.5) +
  geom_point(size = 2.6, alpha = 0.78) +
  geom_point(
    data = filter(cvd, label_point),
    shape = 1, size = 3.5, stroke = 0.55, color = "black",
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = filter(cvd, label_point),
    aes(label = display_label),
    size = 2.6,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.2,
    min.segment.length = 0,
    segment.color = "gray55",
    max.overlaps = Inf,
    seed = 20260825,
    show.legend = FALSE
  ) +
  labs(
    title = "Adjusted association between environmental toxicants and CVD",
    subtitle = sprintf(
      "Survey-weighted logistic regression; nominal P < 0.05: %d/149; BH-FDR < 0.05: %d/149",
      sum(cvd$nominal), sum(cvd$fdr_significant)
    ),
    x = "Odds ratio"
  )

figure5a <- ggplot(le8, aes(x = estimate, y = p_value, color = family)) +
  common_layers +
  geom_hline(
    yintercept = le8_bh_p_cutoff,
    linetype = "dashed", color = "gray45", linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.5) +
  geom_point(size = 2.6, alpha = ifelse(le8$fdr_significant, 0.9, 0.45)) +
  geom_point(
    data = filter(le8, label_point),
    shape = 1, size = 3.5, stroke = 0.5, color = "black",
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = filter(le8, label_point),
    aes(label = display_label),
    size = 2.15,
    color = "black",
    box.padding = 0.2,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.color = "gray55",
    max.overlaps = 22,
    force = 1.2,
    seed = 20260825,
    show.legend = FALSE
  ) +
  labs(
    title = "Adjusted association between environmental toxicants and LE8",
    subtitle = sprintf(
      "Survey-weighted linear regression; BH-FDR < 0.05: %d/149; negative associations labelled: %d",
      sum(le8$fdr_significant), sum(le8$label_point)
    ),
    x = expression(beta~"coefficient for LE8 score")
  )

save_figure <- function(plot, stem, width_mm = 183, height_mm = 145, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(file.path(out_dir, paste0(stem, ".svg")), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()

  grDevices::pdf(
    file.path(out_dir, paste0(stem, ".pdf")),
    width = width_in, height = height_in, family = "Helvetica",
    useDingbats = FALSE
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    file.path(out_dir, paste0(stem, ".tiff")),
    width = width_in, height = height_in, units = "in", res = dpi,
    compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    file.path(out_dir, paste0(stem, "_preview.png")),
    width = width_in, height = height_in, units = "in", res = 220
  )
  print(plot)
  grDevices::dev.off()
}

write_csv(cvd, file.path(out_dir, "Figure1_weighted_source_data.csv"), na = "")
write_csv(le8, file.path(out_dir, "Figure5A_weighted_source_data.csv"), na = "")
write_csv(
  tibble::tibble(
    family = family_order,
    color_hex = unname(family_colors[family_order])
  ),
  file.path(out_dir, "toxicant_family_palette_mapping.csv")
)
save_figure(figure1, "Figure1_weighted")
save_figure(figure5a, "Figure5A_weighted")

cat("Figure 1: nominal", sum(cvd$nominal), "FDR", sum(cvd$fdr_significant), "\n")
cat("Figure 5A: nominal", sum(le8$nominal), "FDR", sum(le8$fdr_significant), "\n")
cat("Output:", out_dir, "\n")
