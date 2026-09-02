#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survey)
})

options(survey.lonely.psu = "adjust")
root <- normalizePath(Sys.getenv("NHANES_REPO", "."))
dat <- readRDS(file.path(root, "data", "processed", "exwas_analysis_long.rds"))
out_dir <- file.path(root, "results", "exwas")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

models <- list(
  CVD = c("age", "sex", "eth", "marital", "PIR_level", "edu", "smoke",
          "alcohol.user", "PA_level", "BMI_level"),
  LE8 = c("age", "sex", "eth", "marital", "PIR_level", "edu", "alcohol.user")
)

results <- list()
k <- 1

for (outcome in names(models)) {
  covariates <- models[[outcome]]

  for (chemical in sort(unique(dat$exposure))) {
    d <- dat |>
      filter(exposure == chemical, age >= 20) |>
      mutate(
        y = if (outcome == "CVD") as.numeric(CVD > 0) else as.numeric(LE8),
        x = if (first(transform) == "binary_above_detection_threshold") factor(exposure_value)
            else as.numeric(exposure_value),
        strata = interaction(cycle, sdmvstra, drop = TRUE),
        psu = interaction(cycle, sdmvpsu, drop = TRUE)
      )

    keep <- c("y", "x", covariates, "analysis_weight", "strata", "psu")
    d <- d[complete.cases(d[, keep]), ] |> droplevels()
    form <- reformulate(c("x", covariates), response = "y")
    design <- svydesign(~psu, strata = ~strata, weights = ~analysis_weight,
                        nest = TRUE, data = d)
    fit <- svyglm(form, design,
                  family = if (outcome == "CVD") quasibinomial() else gaussian())

    co <- summary(fit)$coefficients
    term <- grep("^x", rownames(co))[1]
    beta <- co[term, 1]
    se <- co[term, 2]
    df <- degf(design)
    p <- 2 * pt(abs(beta / se), df, lower.tail = FALSE)
    ci <- beta + c(-1, 1) * qt(0.975, df) * se

    results[[k]] <- tibble(
      exposure = chemical,
      outcome = outcome,
      method = "survey_weighted",
      effect_scale = if (outcome == "CVD") "odds_ratio" else "beta",
      estimate = if (outcome == "CVD") exp(beta) else beta,
      conf_low = if (outcome == "CVD") exp(ci[1]) else ci[1],
      conf_high = if (outcome == "CVD") exp(ci[2]) else ci[2],
      p_value = p,
      n = nrow(d),
      events = if (outcome == "CVD") sum(d$y == 1) else NA_integer_,
      status = "ok"
    )
    k <- k + 1
  }
}

results <- bind_rows(results) |>
  group_by(outcome) |>
  mutate(p_fdr_bh = p.adjust(p_value, "BH")) |>
  ungroup()

write_csv(filter(results, outcome == "CVD"),
          file.path(out_dir, "survey_weighted_exwas_CVD.csv"), na = "")
write_csv(filter(results, outcome == "LE8"),
          file.path(out_dir, "survey_weighted_exwas_LE8_primary.csv"), na = "")
write_csv(results, file.path(out_dir, "survey_weighted_exwas_all_results.csv"), na = "")
