#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mediation)
  library(readr)
})

set.seed(20260827)
root <- normalizePath(Sys.getenv("NHANES_REPO", "."))
long <- readRDS(file.path(root, "data", "processed", "exwas_analysis_long.rds"))
exwas <- read_csv(file.path(root, "results", "exwas", "survey_weighted_exwas_LE8_primary.csv"),
                  show_col_types = FALSE)
description <- read_csv(file.path(root, "data", "data_dictionary", "pollutant_description.csv"),
                        show_col_types = FALSE) |>
  select(exposure, family)
out_dir <- file.path(root, "results", "mediation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidates <- exwas |>
  filter(method == "survey_weighted", p_fdr_bh < 0.05, estimate < 0) |>
  pull(exposure)
covariates <- c("age", "sex", "eth", "marital", "PIR_level", "edu", "alcohol.user")
rows <- list()

for (i in seq_along(candidates)) {
  chemical <- candidates[i]
  dat <- long |>
    filter(exposure == chemical, age >= 20) |>
    transmute(
      CVD = as.numeric(CVD > 0), LE8, treatment = as.numeric(exposure_value),
      model_weight = analysis_weight, across(all_of(covariates))
    )
  dat <- dat[complete.cases(dat), ]
  dat$model_weight <- dat$model_weight / mean(dat$model_weight)

  mediator_formula <- reformulate(c("treatment", covariates), "LE8")
  outcome_formula <- reformulate(c("treatment", "LE8", covariates), "CVD")
  mediator_model <- lm(mediator_formula, dat, weights = model_weight)
  outcome_model <- glm(outcome_formula, dat, weights = model_weight,
                       family = binomial("probit"))
  mediator_model$call$formula <- mediator_formula
  outcome_model$call$formula <- outcome_formula

  fit <- mediate(
    mediator_model, outcome_model,
    treat = "treatment", mediator = "LE8",
    control.value = 0, treat.value = 1,
    boot = TRUE, sims = 1000
  )
  result <- summary(fit)

  rows[[i]] <- tibble(
    exposure = chemical,
    n = nrow(dat),
    CVD_events = sum(dat$CVD == 1),
    indirect_effect = unname(result$d.avg),
    indirect_low = unname(result$d.avg.ci[1]),
    indirect_high = unname(result$d.avg.ci[2]),
    indirect_p = unname(result$d.avg.p),
    direct_effect = unname(result$z.avg),
    direct_low = unname(result$z.avg.ci[1]),
    direct_high = unname(result$z.avg.ci[2]),
    direct_p = unname(result$z.avg.p),
    total_effect = unname(result$tau.coef),
    proportion_mediated = unname(result$n.avg)
  )
}

results <- bind_rows(rows) |>
  left_join(description, by = "exposure") |>
  mutate(indirect_fdr = p.adjust(indirect_p, "BH")) |>
  arrange(indirect_fdr)

write_csv(results, file.path(out_dir, "TableS6_estimated_indirect_effects.csv"), na = "")

plot_data <- results |>
  slice_min(indirect_p, n = 20) |>
  mutate(exposure = reorder(exposure, indirect_effect))

p <- ggplot(plot_data, aes(indirect_effect, exposure)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = indirect_low, xmax = indirect_high), height = 0.15) +
  geom_point(aes(colour = family), size = 2.5) +
  labs(x = "Estimated indirect effect through LE8", y = NULL, colour = NULL) +
  theme_classic(base_size = 10)

ggsave(file.path(out_dir, "Figure5B_weighted_mediation.pdf"), p, width = 7, height = 6)
