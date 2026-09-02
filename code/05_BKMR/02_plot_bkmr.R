suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: plot_bkmr_cvd_overall_10groups.R INPUT_CSV OUTPUT_DIR")
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- sprintf("C%02d", 1:11)
group_labels <- c(
  C01 = "C01\nFlame retardants",
  C02 = "C02\nSerum metals",
  C03 = "C03\nUrinary metals",
  C04 = "C04\nPAH metabolites",
  C05 = "C05\nPhenols/parabens",
  C06 = "C06\nPerchlorate/nitrate/SCN",
  C07 = "C07\nPFAS",
  C08 = "C08\nPhthalates",
  C09 = "C09\nSerum VOCs",
  C10 = "C10\nVOC metabolites",
  C11 = "C11\nPesticides"
)

dat <- read_csv(input_file, show_col_types = FALSE)
required <- c("category_group_id", "outcome", "quantile", "est", "sd")
missing_cols <- setdiff(required, names(dat))
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cvd <- dat %>%
  filter(outcome == "CVD", category_group_id %in% group_levels) %>%
  mutate(
    category_group_id = factor(category_group_id, levels = group_levels),
    quantile = as.numeric(quantile),
    est = as.numeric(est),
    sd = as.numeric(sd),
    lower = est - 1.96 * sd,
    upper = est + 1.96 * sd,
    interval_status = factor(
      if_else(lower > 0 | upper < 0,
              "95% CrI excludes 0", "95% CrI includes 0"),
      levels = c("95% CrI includes 0", "95% CrI excludes 0")
    )
  ) %>%
  arrange(category_group_id, quantile)

expected_rows <- length(group_levels) * 11L
if (nrow(cvd) != expected_rows) {
  warning("Expected ", expected_rows, " CVD rows but found ", nrow(cvd), ".")
}

source_data <- cvd %>%
  transmute(
    category_group_id = as.character(category_group_id),
    outcome,
    mixture_quantile = quantile,
    posterior_estimate = est,
    posterior_sd = sd,
    credible_interval_95_low = lower,
    credible_interval_95_high = upper,
    interval_status
  )
write_csv(source_data, file.path(output_dir, "BKMR_CVD_overall_effect_10groups_source_data.csv"))

palette <- c(
  "95% CrI includes 0" = "#4E6677",
  "95% CrI excludes 0" = "#B24745"
)

p <- ggplot(cvd, aes(x = quantile, y = est, group = category_group_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35,
             colour = "#6B6B6B") +
  geom_line(linewidth = 0.38, colour = "#7A7A7A") +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, colour = interval_status),
    width = 0.012, linewidth = 0.42
  ) +
  geom_point(aes(colour = interval_status), size = 1.65) +
  facet_wrap(
    ~ category_group_id,
    ncol = 5,
    labeller = as_labeller(group_labels),
    scales = "fixed"
  ) +
  scale_colour_manual(values = palette, drop = FALSE) +
  scale_x_continuous(
    breaks = c(0.25, 0.50, 0.75),
    labels = c("25", "50", "75"),
    limits = c(0.235, 0.765),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = c(-0.10, -0.05, 0, 0.05, 0.10, 0.15),
    labels = scales::label_number(accuracy = 0.01),
    limits = c(-0.12, 0.15),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Mixture percentile",
    y = "Overall mixture effect on CVD (probit scale)",
    colour = NULL
  ) +
  theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.32, colour = "black"),
    axis.ticks = element_line(linewidth = 0.32, colour = "black"),
    axis.ticks.length = grid::unit(1.2, "mm"),
    axis.title.x = element_text(size = 7.2, margin = margin(t = 5)),
    axis.title.y = element_text(size = 7.2, margin = margin(r = 5)),
    axis.text = element_text(size = 6.2, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.3, face = "bold", lineheight = 0.95,
                              margin = margin(b = 3)),
    panel.spacing.x = grid::unit(4.0, "mm"),
    panel.spacing.y = grid::unit(5.0, "mm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 6.2),
    legend.key.width = grid::unit(4.0, "mm"),
    legend.spacing.x = grid::unit(1.5, "mm"),
    plot.margin = margin(4, 5, 3, 4)
  ) +
  guides(colour = guide_legend(
    override.aes = list(linewidth = 0.7, size = 1.8),
    nrow = 1,
    byrow = TRUE
  ))

width_mm <- 183
height_mm <- 108
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
base <- file.path(output_dir, "BKMR_CVD_overall_effect_10groups")

svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in,
                 bg = "white")
print(p)
dev.off()

p_pdf <- p + theme(text = element_text(family = "Helvetica"))
grDevices::pdf(paste0(base, ".pdf"), width = width_in, height = height_in,
               family = "Helvetica", bg = "white", useDingbats = FALSE,
               compress = TRUE)
print(p_pdf)
dev.off()

ragg::agg_tiff(paste0(base, ".tiff"), width = width_in, height = height_in,
               units = "in", res = 600, background = "white",
               scaling = 1)
print(p)
dev.off()

ragg::agg_png(paste0(base, ".png"), width = width_in, height = height_in,
              units = "in", res = 600, background = "white",
              scaling = 1)
print(p)
dev.off()

legend_text <- c(
  "Fig. X | Overall mixture effects on cardiovascular disease across ten pollutant groups.",
  "Points show posterior mean differences in the BKMR latent probit function when all",
  "exposures in a group are jointly fixed at the indicated percentile, relative to the",
  "50th percentile. Error bars denote 95% credible intervals calculated as posterior",
  "mean ± 1.96 × posterior SD. Red points and intervals exclude zero; blue-grey intervals",
  "include zero. Models used 20,000 MCMC iterations, a 10,000-iteration burn-in, 30 knots,",
  "and component-wise variable selection. Binary detection indicators were excluded.",
  "The BKMR models did not incorporate NHANES survey weights, strata, or PSUs and are",
  "therefore interpreted as exploratory analyses."
)
writeLines(legend_text, file.path(output_dir, "BKMR_CVD_overall_effect_10groups_legend.txt"),
           useBytes = TRUE)

qa_notes <- c(
  "BKMR CVD overall-effect figure QA",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Backend: R only (ggplot2, svglite, grDevices::pdf, ragg)",
  "Archetype: quantitative grid, 5 columns x 2 rows",
  "Final size: 183 mm x 108 mm",
  "Reference: all exposures at the 50th percentile; null line at 0",
  "Interval: posterior estimate ± 1.96 x posterior SD",
  "Y-axis: fixed across all ten panels (-0.12 to 0.15) for valid comparison",
  "Exports: editable SVG, PDF, 600-dpi TIFF, 600-dpi PNG",
  "PDF device note: the native R PDF device was used because Cairo/XQuartz libraries",
  "were unavailable in the local R runtime; Helvetica was embedded as a standard font.",
  "Reviewer risk: MCMC diagnostics should be reported separately; significant colour does",
  "not override convergence concerns or the exploratory non-survey-weighted design."
)
writeLines(qa_notes, file.path(output_dir, "BKMR_CVD_overall_effect_10groups_QA.txt"),
           useBytes = TRUE)

message("Saved figure bundle to: ", normalizePath(output_dir))
