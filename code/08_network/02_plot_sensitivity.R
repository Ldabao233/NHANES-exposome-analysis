# 06_plot_supplementary_sensitivity.R
# Expected input:
# tables/sensitivity/Sensitivity_enrichment_long.csv
#
# Required columns:
# Scenario, Category, Description, FDR
#
# Scenario values:
# Full, Minus_CTD, Minus_Swiss, Minus_SEA

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

repo_root <- normalizePath(Sys.getenv("NHANES_REPO", "."), mustWork = TRUE)
setwd(file.path(repo_root, "results", "network"))
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

sens <- read_csv(
  "tables/sensitivity/Sensitivity_enrichment_long.csv",
  show_col_types = FALSE
) %>%
  mutate(
    FDR = as.numeric(FDR),
    Score = -log10(pmax(FDR, 1e-300)),
    Scenario = factor(
      Scenario,
      levels = c("Full", "Minus_CTD", "Minus_Swiss", "Minus_SEA"),
      labels = c("Full", "-CTD", "-Swiss", "-SEA")
    )
  )

# Select pathways that are significant in the Full model, then retain top 20.
terms_keep <- sens %>%
  filter(Scenario == "Full", FDR < 0.05) %>%
  arrange(FDR) %>%
  distinct(Description) %>%
  slice_head(n = 20) %>%
  pull(Description)

plot_dat <- sens %>%
  filter(Description %in% terms_keep) %>%
  mutate(
    Description = factor(Description, levels = rev(terms_keep))
  )

p <- ggplot(plot_dat, aes(x = Scenario, y = Description, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(
    low = "#F0F8FF",
    high = "#1F78C8",
    name = "-log10(FDR)"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 0, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 183, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = width_in, height = height_in)
  print(plot)
  dev.off()
  grDevices::pdf(
    paste0(filename, ".pdf"), width = width_in, height = height_in,
    family = "Helvetica", useDingbats = FALSE
  )
  print(plot)
  dev.off()
  ragg::agg_tiff(
    paste0(filename, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw"
  )
  print(plot)
  dev.off()
  ragg::agg_png(
    paste0(filename, ".png"), width = width_in, height = height_in,
    units = "in", res = 180
  )
  print(plot)
  dev.off()
}

save_pub_r(
  p,
  "figures/FigureS_leave_one_source_out_enrichment_heatmap",
  width_mm = 183,
  height_mm = 183
)
