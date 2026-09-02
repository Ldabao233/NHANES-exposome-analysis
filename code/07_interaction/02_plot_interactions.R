# ==============================================================================
# Figure 6. Toxicant–lifestyle effect-modification analyses
#
# Required input files:
#   1. component_primary_continuous_interaction.csv
#   2. component_RERI_bootstrap.csv
#   3. individual_toxicant_primary_interaction_256.csv
#   4. LE8_interaction_exploratory_signals.csv
#
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
  library(scales)
})

repo_root <- normalizePath(Sys.getenv("NHANES_REPO", "."), mustWork = TRUE)
input_dir <- file.path(repo_root, "results", "interaction")

# ---- CVD interaction files ----
component_primary_file <- file.path(
  input_dir, "component_primary_continuous_interaction.csv"
)

component_reri_file <- file.path(
  input_dir, "component_RERI_bootstrap.csv"
)

individual_primary_file <- file.path(
  input_dir, "individual_toxicant_primary_interaction_256.csv"
)

# ---- LE8 interaction file ----
le8_file <- file.path(
  input_dir, "LE8_interaction_exploratory_signals.csv"
)

out_dir <- file.path(
  repo_root, "results", "interaction", "figures"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c(
  component_primary_file,
  component_reri_file,
  individual_primary_file,
  le8_file
)) {
  if (!file.exists(f)) stop("Cannot find file: ", f)
}

# ------------------------------------------------------------------------------
# 1. READ DATA
# ------------------------------------------------------------------------------

cvd_component <- read_csv(component_primary_file, show_col_types = FALSE)
cvd_reri      <- read_csv(component_reri_file, show_col_types = FALSE)
cvd_individual <- read_csv(individual_primary_file, show_col_types = FALSE)
le8 <- read_csv(le8_file, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 2. COMMON THEME
# ------------------------------------------------------------------------------

theme_pub <- theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9.5, color = "gray30"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 9.5),
    plot.margin = margin(8, 12, 8, 8)
  )

palette <- c(fdr05 = "#D73027", fdr10 = "#2C7FB8")

# ==============================================================================
# PANEL A
# CVD MULTIPLICATIVE INTERACTION FOREST PLOT
# ==============================================================================

# 2-NP comes from the component-level analysis.
a_2np <- cvd_component %>%
  filter(
    Pollutant == "X2_NP_ln_U",
    Modifier == "Smoking"
  ) %>%
  transmute(
    Chemical = "2-NP",
    Modifier = "Smoking",
    Estimate = OR_Interaction_Ratio,
    Low = OR_Interaction_Low,
    High = OR_Interaction_High,
    P = P_Interaction,
    FDR = P_FDR_BH_within_family
  )

# MiBP and thiocyanate come from the individual-toxicant analysis.
a_other <- cvd_individual %>%
  filter(
    (Pollutant == "MiBP_ln_U" & Modifier == "Smoking") |
    (Pollutant == "Thiocyanate_ln_U" & Modifier == "Smoking") |
    (Pollutant == "Thiocyanate_ln_U" & Modifier == "Diet")
  ) %>%
  transmute(
    Chemical = case_when(
      Pollutant == "MiBP_ln_U" ~ "MiBP",
      Pollutant == "Thiocyanate_ln_U" ~ "Thiocyanate",
      TRUE ~ Pollutant_Label
    ),
    Modifier,
    Estimate = OR_Interaction_Ratio,
    Low = OR_Interaction_Low,
    High = OR_Interaction_High,
    P = P_Interaction,
    FDR = P_FDR_BH_within_family
  )

plotA_data <- bind_rows(a_2np, a_other) %>%
  mutate(
    Label = paste0(Chemical, " \u00D7 ", Modifier),
    Significance = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      FDR < 0.10 ~ "FDR < 0.10",
      TRUE ~ "FDR \u2265 0.10"
    ),
    FDR_label = paste0("FDR = ", sprintf("%.3f", FDR)),
    Direction = if_else(Estimate > 1, "Positive", "Negative")
  ) %>%
  mutate(
    Label = factor(
      Label,
      levels = rev(c(
        "2-NP \u00D7 Smoking",
        "MiBP \u00D7 Smoking",
        "Thiocyanate \u00D7 Smoking",
        "Thiocyanate \u00D7 Diet"
      ))
    )
  )

pA <- ggplot(
  plotA_data,
  aes(
    y = Label,
    x = Estimate,
    xmin = Low,
    xmax = High,
    shape = Direction
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.55,
    color = "gray55"
  ) +
  geom_errorbarh(
    aes(color = Significance),
    height = 0.16,
    linewidth = 0.75
  ) +
  geom_point(
    aes(color = Significance),
    size = 3.0,
    stroke = 0.9
  ) +
  geom_text(
    aes(
      y = Label,
      label = FDR_label
    ),
    x = max(plotA_data$High) * 1.08,
    hjust = 0,
    size = 3.1,
    inherit.aes = FALSE,
    data = plotA_data
  ) +
  scale_shape_manual(
    values = c("Positive" = 16, "Negative" = 17),
    name = "Effect direction"
  ) +
  scale_color_manual(
    values = c("FDR < 0.05" = unname(palette["fdr05"]), "FDR < 0.10" = unname(palette["fdr10"]), "FDR ≥ 0.10" = "grey70"),
    name = "FDR"
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
    labels = number_format(accuracy = 0.01),
    expand = expansion(mult = c(0.05, 0.34))
  ) +
  labs(
    title = "A  CVD: multiplicative interactions",
    x = "Interaction OR ratio (95% CI)",
    y = NULL,
    shape = NULL
  ) +
  theme_pub

# ==============================================================================
# PANEL B
# CVD ADDITIVE INTERACTIONS: RERI FOREST PLOT
# ==============================================================================

# Smoking additive interactions surviving within-family BH-FDR < 0.05.
plotB_data <- cvd_reri %>%
  filter(Modifier == "Smoking", RERI_P_FDR_BH_within_family < 0.05,
         Any_Sparse_Cell == FALSE, status == "ok") %>%
  mutate(
    Chemical = case_when(
      Pollutant == "X2_NP_ln_U" ~ "2-NP", Pollutant == "X1_NP_ln_U" ~ "1-NP",
      Pollutant == "X2_FL_ln_U" ~ "2-FL", Pollutant == "MB3_Cys_ln_U" ~ "MB3-Cys",
      Pollutant == "NPMA_ln_U" ~ "NPMA", Pollutant == "X3_HPMA_ln_U" ~ "3-HPMA",
      Pollutant == "X34_MHA_ln_U" ~ "3/4-MHA", Pollutant == "X2_MHA_ln_U" ~ "2-MHA",
      TRUE ~ Pollutant_Label
    ),
    Label = paste0(Chemical, " × Smoking"),
    FDR_label = paste0("FDR = ", sprintf("%.3f", RERI_P_FDR_BH_within_family)),
    Significance = if_else(RERI_P_FDR_BH_within_family < 0.05, "FDR < 0.05", "FDR < 0.10"),
    Direction = if_else(RERI > 0, "Positive", "Negative")
  ) %>% arrange(RERI) %>% mutate(Label = factor(Label, levels = Label))

pB <- ggplot(plotB_data, aes(y = Label, x = RERI, xmin = RERI_Low, xmax = RERI_High)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.55, color = "gray55") +
  geom_errorbarh(aes(color = Significance), height = 0.16, linewidth = 0.75) +
  geom_point(aes(shape = Direction, color = Significance), size = 3, stroke = 0.9) +
  geom_text(aes(label = FDR_label), x = max(plotB_data$RERI_High) + 0.35, hjust = 0, size = 3.0) +
  scale_shape_manual(
    values = c("Positive" = 16, "Negative" = 17),
    guide = "none"
  ) +
  scale_color_manual(values = c("FDR < 0.05" = unname(palette["fdr05"]), "FDR < 0.10" = unname(palette["fdr10"])), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.24))) +
  labs(title = "B  CVD: additive interactions",
       subtitle = "Smoking-related interactions surviving within-family BH-FDR correction",
       x = "RERI (95% CI)", y = NULL, shape = NULL) +
  theme_pub

# ==============================================================================
# PANEL C
# LE8 INTERACTION FOREST PLOT: NDMA × PA and MeP × PA
# ==============================================================================

plotC_data <- le8 %>%
  filter(
    Modifier == "PhysicalActivity",
    Pollutant %in% c("NDMA_binary_U", "MeP_ln_U")
  ) %>%
  transmute(
    Chemical = case_when(
      Pollutant == "NDMA_binary_U" ~ "NDMA",
      Pollutant == "MeP_ln_U" ~ "MeP",
      TRUE ~ Pollutant
    ),
    Estimate = Beta_Interaction,
    Low = Interaction_Low,
    High = Interaction_High,
    FDR = P_FDR_BH_within_modifier
  ) %>%
  mutate(
    Label = paste0(Chemical, " \u00D7 Physical activity"),
    Significance = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      FDR < 0.10 ~ "FDR < 0.10",
      TRUE ~ "FDR \u2265 0.10"
    ),
    FDR_label = paste0("FDR = ", sprintf("%.3f", FDR)),
    Direction = if_else(Estimate > 0, "Positive", "Negative"),
    Label = factor(
      Label,
      levels = rev(c(
        "NDMA \u00D7 Physical activity",
        "MeP \u00D7 Physical activity"
      ))
    )
  )

pC <- ggplot(
  plotC_data,
  aes(
    y = Label,
    x = Estimate,
    xmin = Low,
    xmax = High,
    shape = Direction
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "gray55"
  ) +
  geom_errorbarh(
    aes(color = Significance),
    height = 0.16,
    linewidth = 0.8
  ) +
  geom_point(
    aes(color = Significance),
    size = 3.1,
    stroke = 0.9
  ) +
  geom_text(
    aes(
      label = FDR_label
    ),
    x = max(plotC_data$High) + 1.8,
    hjust = 0,
    size = 3.15
  ) +
  scale_shape_manual(
    values = c("Positive" = 16, "Negative" = 17),
    guide = "none"
  ) +
  scale_color_manual(
    values = c("FDR < 0.05" = unname(palette["fdr05"]), "FDR < 0.10" = unname(palette["fdr10"]), "FDR ≥ 0.10" = "grey70"),
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.06, 0.30))
  ) +
  labs(
    title = "C  LE8: lifestyle effect modification",
    subtitle = "Outcome is LE8 excluding the corresponding lifestyle component",
    x = expression(beta[interaction]~"(95% CI)"),
    y = NULL,
    shape = NULL
  ) +
  theme_pub +
  theme(legend.position = "top")

# ==============================================================================
# PANEL D
# NDMA × PHYSICAL ACTIVITY: STRATIFIED LE8_no_PA ASSOCIATIONS
# ==============================================================================

row_ndma <- le8 %>%
  filter(
    Pollutant == "NDMA_binary_U",
    Modifier == "PhysicalActivity"
  )

if (nrow(row_ndma) != 1) {
  stop("Expected exactly one NDMA_binary_U × PhysicalActivity row.")
}

plotD_data <- tibble(
  Group = c(
    "Ideal physical activity",
    "Non-ideal physical activity"
  ),
  Estimate = c(
    row_ndma$Toxicant_Beta_Healthy,
    row_ndma$Toxicant_Beta_Risk
  ),
  Low = c(
    row_ndma$Toxicant_Beta_Healthy_Low,
    row_ndma$Toxicant_Beta_Risk_Low
  ),
  High = c(
    row_ndma$Toxicant_Beta_Healthy_High,
    row_ndma$Toxicant_Beta_Risk_High
  )
) %>%
  mutate(
    Group = factor(
      Group,
      levels = rev(c(
        "Ideal physical activity",
        "Non-ideal physical activity"
      ))
    ),
    Text = sprintf("%.2f (%.2f-%.2f)", Estimate, Low, High),
    Significance = "FDR < 0.05",
    Direction = if_else(Estimate > 0, "Positive", "Negative")
  )

d_annotation <- paste0(
  "Beta_interaction = ",
  sprintf("%.2f", row_ndma$Beta_Interaction),
  " (95% CI ",
  sprintf("%.2f", row_ndma$Interaction_Low),
  " to ",
  sprintf("%.2f", row_ndma$Interaction_High),
  "); FDR = ",
  sprintf("%.3f", row_ndma$P_FDR_BH_within_modifier)
)

pD <- ggplot(
  plotD_data,
  aes(
    y = Group,
    x = Estimate,
    xmin = Low,
    xmax = High
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "gray55"
  ) +
  geom_errorbarh(
    aes(color = Significance),
    height = 0.16,
    linewidth = 0.8
  ) +
  geom_point(
    aes(shape = Direction, color = Significance),
    size = 3.2
  ) +
  geom_text(
    aes(
      label = Text
    ),
    x = max(plotD_data$High) + 1.0,
    hjust = 0,
    size = 3.15
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.08, 0.34))
  ) +
  scale_color_manual(
    values = c("FDR < 0.05" = unname(palette["fdr05"])),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("Positive" = 16, "Negative" = 17),
    guide = "none"
  ) +
  labs(
    title = "D  NDMA association with LE8 by physical activity",
    subtitle = d_annotation,
    x = expression(beta~"for NDMA exposure on LE8"[no-PA]),
    y = NULL
  ) +
  theme_pub +
  theme(
    legend.position = "none"
  )

# ==============================================================================
# 3. SAVE INDIVIDUAL PANELS
# ==============================================================================

ggsave(
  file.path(out_dir, "Figure6A_CVD_interaction_forest.pdf"),
  pA,
  width = 7.2,
  height = 4.2,
  device = "pdf"
)

ggsave(
  file.path(out_dir, "Figure6B_2NP_smoking_stratified.pdf"),
  pB,
  width = 7.2,
  height = 4.2,
  device = "pdf"
)

ggsave(
  file.path(out_dir, "Figure6C_LE8_interaction_forest.pdf"),
  pC,
  width = 7.2,
  height = 4.2,
  device = "pdf"
)

ggsave(
  file.path(out_dir, "Figure6D_NDMA_PA_stratified.pdf"),
  pD,
  width = 7.2,
  height = 4.2,
  device = "pdf"
)

# ==============================================================================
# 4. COMBINE A–D
# ==============================================================================

Figure6 <- (pA | pB) /
           (pC | pD) +
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "top"
  )

print(Figure6)

ggsave(
  file.path(out_dir, "Figure6_interaction_combined.pdf"),
  Figure6,
  width = 14,
  height = 9,
  device = "pdf"
)

ggsave(
  file.path(out_dir, "Figure6_interaction_combined.tiff"),
  Figure6,
  width = 14,
  height = 9,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ==============================================================================
# 5. EXPORT PLOTTED NUMBERS FOR AUDIT
# ==============================================================================

write_csv(
  plotA_data,
  file.path(out_dir, "Figure6A_plot_data.csv")
)

write_csv(
  plotB_data,
  file.path(out_dir, "Figure6B_plot_data.csv")
)

write_csv(
  plotC_data,
  file.path(out_dir, "Figure6C_plot_data.csv")
)

write_csv(
  plotD_data,
  file.path(out_dir, "Figure6D_plot_data.csv")
)

cat("\nFigure 6 completed.\n")
cat("Output directory:\n", out_dir, "\n")
