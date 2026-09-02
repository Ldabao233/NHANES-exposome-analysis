#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(gWQS)
  library(readr)
  library(survey)
})

options(survey.lonely.psu = "adjust")
root <- normalizePath(Sys.getenv("NHANES_REPO", "."))
groups <- read_csv(file.path(root, "data", "data_dictionary", "mixture_groups.csv"),
                   show_col_types = FALSE)$category_group_id
data_dir <- file.path(root, "data", "processed", "mixtures")
out_dir <- file.path(root, "results", "WQS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_rows <- list()
estimate_rows <- list()
weight_rows <- list()
count_rows <- list()
k <- 1

for (group_id in groups) {
  object <- readRDS(file.path(data_dir, paste0(group_id, "_mixture_data.rds")))

  for (outcome in c("CVD", "LE8")) {
    exposures <- object$exposures
    covariates <- if (outcome == "CVD") object$cvd_covariates else object$le8_covariates
    needed <- c(outcome, exposures, covariates, "analysis_weight", "strata_cycle", "psu_cycle")
    dat <- object$data |>
      filter(age >= 20) |>
      filter(if_all(all_of(needed), ~ !is.na(.x))) |>
      mutate(wqs_weight = analysis_weight / mean(analysis_weight)) |>
      droplevels()

    for (direction in c("positive", "negative")) {
      set.seed(20260826 + k)
      fit <- gwqs(
        reformulate(c("wqs", covariates), outcome),
        data = dat,
        mix_name = exposures,
        q = 4,
        validation = 0.40,
        rh = 50,
        b = 1000,
        b1_pos = direction == "positive",
        b_constr = TRUE,
        family = if (outcome == "CVD") "binomial" else "gaussian",
        signal = "t2",
        weights = "wqs_weight",
        seed = 20260826 + k,
        plan_strategy = "sequential",
        plots = FALSE,
        tables = FALSE
      )

      job_estimates <- list()
      job_weights <- list()
      job_counts <- list()

      for (h in seq_along(fit$gwqslist)) {
        one <- fit$gwqslist[[h]]
        validation <- one$validation_rows
        weights <- one$final_weights |>
          transmute(exposure = mix_name, weight = mean_weight)

        scores <- matrix(NA_real_, nrow(dat), length(exposures),
                         dimnames = list(NULL, exposures))
        for (chemical in exposures) {
          scores[, chemical] <- as.numeric(cut(
            dat[[chemical]], one$qi[[chemical]], labels = FALSE,
            include.lowest = TRUE
          )) - 1
        }
        dat$wqs <- as.numeric(scores %*% setNames(weights$weight, weights$exposure)[exposures])

        validation_data <- dat[validation, ]
        design <- svydesign(~psu_cycle, strata = ~strata_cycle,
                            weights = ~analysis_weight, nest = TRUE,
                            data = validation_data)
        survey_fit <- svyglm(
          reformulate(c("wqs", covariates), outcome), design,
          family = if (outcome == "CVD") quasibinomial() else gaussian()
        )
        beta <- coef(survey_fit)["wqs"]
        se <- sqrt(vcov(survey_fit)["wqs", "wqs"])
        p <- 2 * pt(abs(beta / se), degf(design), lower.tail = FALSE)

        job_estimates[[h]] <- tibble(
          category_group_id = group_id, outcome = outcome, direction = direction,
          holdout = h, estimate = beta, std_error = se, p_value = p,
          effect = if (outcome == "CVD") exp(beta) else beta
        )
        job_weights[[h]] <- weights |>
          mutate(category_group_id = group_id, outcome = outcome,
                 direction = direction, holdout = h, .before = 1)
        job_counts[[h]] <- tibble(
          category_group_id = group_id, outcome = outcome, direction = direction,
          holdout = h, n_training = sum(!validation), n_validation = sum(validation),
          training_CVD_events = if (outcome == "CVD") sum(dat$CVD[!validation] == 1) else NA_integer_,
          validation_CVD_events = if (outcome == "CVD") sum(dat$CVD[validation] == 1) else NA_integer_
        )
      }

      estimates <- bind_rows(job_estimates)
      qs <- quantile(estimates$estimate, c(0.025, 0.5, 0.975))
      model_rows[[k]] <- tibble(
        category_group_id = group_id, outcome = outcome, direction = direction,
        n_complete = nrow(dat), repeated_holdouts = 50, bootstraps = 1000,
        mean_estimate = mean(estimates$estimate), median_estimate = qs[2],
        percentile_ci_low = qs[1], percentile_ci_high = qs[3],
        empirical_two_sided_p = min(1, 2 * min(mean(estimates$estimate <= 0),
                                               mean(estimates$estimate >= 0)))
      )
      estimate_rows[[k]] <- estimates
      weight_rows[[k]] <- bind_rows(job_weights)
      count_rows[[k]] <- bind_rows(job_counts)
      k <- k + 1
    }
  }
}

models <- bind_rows(model_rows) |>
  group_by(outcome, direction) |>
  mutate(BH_FDR = p.adjust(empirical_two_sided_p, "BH")) |>
  ungroup()
weights <- bind_rows(weight_rows)
weight_summary <- weights |>
  group_by(category_group_id, outcome, direction, exposure) |>
  summarise(mean_weight = mean(weight), sd_weight = sd(weight), .groups = "drop")

write_csv(models, file.path(out_dir, "formal_rhWQS_model_summary.csv"), na = "")
write_csv(bind_rows(estimate_rows), file.path(out_dir, "formal_rhWQS_holdout_estimates.csv"), na = "")
write_csv(weights, file.path(out_dir, "formal_rhWQS_holdout_weights.csv"), na = "")
write_csv(weight_summary, file.path(out_dir, "formal_rhWQS_weight_summary.csv"), na = "")
write_csv(bind_rows(count_rows), file.path(out_dir, "formal_rhWQS_training_validation_counts.csv"), na = "")
