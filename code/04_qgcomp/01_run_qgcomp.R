#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(qgcomp)
  library(readr)
  library(tibble)
})

root <- normalizePath(Sys.getenv("NHANES_REPO", "."))
groups <- read_csv(file.path(root, "data", "data_dictionary", "mixture_groups.csv"),
                   show_col_types = FALSE)$category_group_id
data_dir <- file.path(root, "data", "processed", "mixtures")
out_dir <- file.path(root, "results", "qgcomp")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()
weight_rows <- list()
k <- 1

for (group_id in groups) {
  object <- readRDS(file.path(data_dir, paste0(group_id, "_mixture_data.rds")))

  for (outcome in c("CVD", "LE8")) {
    exposures <- object$exposures
    covariates <- if (outcome == "CVD") object$cvd_covariates else object$le8_covariates
    needed <- c(outcome, exposures, covariates, "analysis_weight")
    dat <- object$data |>
      filter(age >= 20) |>
      filter(if_all(all_of(needed), ~ !is.na(.x))) |>
      mutate(model_weight = analysis_weight / mean(analysis_weight)) |>
      droplevels()

    form <- reformulate(c(exposures, covariates), response = outcome)
    family <- if (outcome == "CVD") binomial() else gaussian()
    set.seed(20260827 + k)

    point <- qgcomp.glm.noboot(
      form, dat, expnms = exposures, q = 4,
      family = family, weights = dat$model_weight, bayes = TRUE
    )
    fit <- qgcomp.glm.boot(
      form, dat, expnms = exposures, q = 4,
      family = family, weights = dat$model_weight,
      B = 1000, seed = 20260827 + k, bayes = TRUE
    )

    psi <- unname(fit$psi[1])
    ci <- unname(fit$ci[1:2])
    summary_rows[[k]] <- tibble(
      category_group_id = group_id,
      family = first(object$data$family),
      outcome = outcome,
      n_exposures = length(exposures),
      n_complete = nrow(dat),
      CVD_events = if (outcome == "CVD") sum(dat$CVD == 1) else NA_integer_,
      estimate = psi,
      ci_low = ci[1],
      ci_high = ci[2],
      p_value = unname(fit$pval["psi1"]),
      effect = if (outcome == "CVD") exp(psi) else psi,
      effect_ci_low = if (outcome == "CVD") exp(ci[1]) else ci[1],
      effect_ci_high = if (outcome == "CVD") exp(ci[2]) else ci[2]
    )

    weight_rows[[k]] <- bind_rows(
      enframe(point$pos.weights, "exposure", "weight") |> mutate(direction = "positive"),
      enframe(point$neg.weights, "exposure", "weight") |> mutate(direction = "negative")
    ) |>
      mutate(category_group_id = group_id, outcome = outcome, .before = 1)
    k <- k + 1
  }
}

summary_table <- bind_rows(summary_rows) |>
  group_by(outcome) |>
  mutate(BH_FDR = p.adjust(p_value, "BH")) |>
  ungroup()

write_csv(summary_table, file.path(out_dir, "qgcomp_model_summary.csv"), na = "")
write_csv(bind_rows(weight_rows), file.path(out_dir, "qgcomp_directional_weights.csv"), na = "")
