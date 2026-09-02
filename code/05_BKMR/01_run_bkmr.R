#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(bkmr)
  library(dplyr)
  library(fields)
  library(readr)
})

root <- normalizePath(Sys.getenv("NHANES_REPO", "."))
group_table <- read_csv(
  file.path(root, "data", "data_dictionary", "mixture_groups.csv"),
  show_col_types = FALSE
) |>
  filter(n_continuous >= 3)
data_dir <- file.path(root, "data", "processed", "mixtures")
out_dir <- file.path(root, "results", "BKMR")
model_dir <- file.path(out_dir, "models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

iterations <- 5000
burn_in <- 2500
summary_rows <- list()
pip_rows <- list()
overall_rows <- list()
single_rows <- list()
k <- 1

for (group_id in group_table$category_group_id) {
  object <- readRDS(file.path(data_dir, paste0(group_id, "_mixture_data.rds")))

  for (outcome in c("CVD", "LE8")) {
    exposures <- object$continuous_exposures
    covariates <- if (outcome == "CVD") object$cvd_covariates else object$le8_covariates
    needed <- c(outcome, exposures, covariates)
    dat <- object$data |>
      filter(age >= 20) |>
      filter(if_all(all_of(needed), ~ !is.na(.x))) |>
      droplevels()

    Z <- scale(as.matrix(dat[, exposures]))
    X <- model.matrix(reformulate(covariates), dat)[, -1, drop = FALSE]
    X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
    y <- as.numeric(dat[[outcome]])
    knots <- cover.design(unique(Z), nd = min(50, nrow(unique(Z))))$design

    set.seed(20260828 + k)
    fit <- kmbayes(
      y = y, Z = Z, X = X, iter = iterations,
      family = if (outcome == "CVD") "binomial" else "gaussian",
      varsel = TRUE, knots = knots, verbose = TRUE
    )
    keep <- (burn_in + 1):iterations

    summary_rows[[k]] <- tibble(
      category_group_id = group_id,
      family = first(object$data$family),
      outcome = outcome,
      n_exposures = length(exposures),
      n_complete = nrow(dat),
      CVD_events = if (outcome == "CVD") sum(y == 1) else NA_integer_,
      iterations = iterations,
      burn_in = burn_in,
      knots = nrow(knots)
    )
    pip_rows[[k]] <- ExtractPIPs(fit, sel = keep, z.names = exposures) |>
      mutate(category_group_id = group_id, outcome = outcome, .before = 1)
    overall_rows[[k]] <- OverallRiskSummaries(fit, qs = seq(0.25, 0.75, 0.05),
                                              q.fixed = 0.50, sel = keep) |>
      mutate(category_group_id = group_id, outcome = outcome, .before = 1)
    single_rows[[k]] <- SingVarRiskSummaries(fit, qs.diff = c(0.25, 0.75),
                                             q.fixed = 0.50, sel = keep) |>
      mutate(category_group_id = group_id, outcome = outcome, .before = 1)

    saveRDS(fit, file.path(model_dir, paste0(group_id, "_", outcome, "_fit.rds")))
    k <- k + 1
  }
}

write_csv(bind_rows(summary_rows), file.path(out_dir, "BKMR_model_summary.csv"), na = "")
write_csv(bind_rows(pip_rows), file.path(out_dir, "BKMR_PIPs.csv"), na = "")
write_csv(bind_rows(overall_rows), file.path(out_dir, "BKMR_overall_risk.csv"), na = "")
write_csv(bind_rows(single_rows), file.path(out_dir, "BKMR_single_variable_risk.csv"), na = "")
