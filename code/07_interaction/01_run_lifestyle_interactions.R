#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survey)
  library(tidyr)
})

options(survey.lonely.psu = "adjust")
root <- normalizePath(Sys.getenv("NHANES_REPO", "."))
long <- readRDS(file.path(root, "data", "processed", "exwas_analysis_long.rds"))
phenotype <- read_csv(file.path(root, "data", "processed", "phenotype.csv"),
                      show_col_types = FALSE) |>
  select(id, score_hei, score_pa, score_smoke, score_sleep) |>
  rename(SEQN = id)
long <- left_join(long, phenotype, by = "SEQN")
out_dir <- file.path(root, "results", "interaction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

modifiers <- tribble(
  ~Modifier, ~score, ~threshold,
  "Diet", "score_hei", 50,
  "Physical activity", "score_pa", 100,
  "Smoking", "score_smoke", 100,
  "Sleep", "score_sleep", 100
)
covariates <- c("age", "sex", "eth", "marital", "PIR_level", "edu", "alcohol.user")
model_rows <- list()
cell_rows <- list()
k <- 1

for (chemical in sort(unique(long$exposure))) {
  source <- filter(long, exposure == chemical, age >= 20)

  for (j in seq_len(nrow(modifiers))) {
    score <- modifiers$score[j]
    dat <- source |>
      transmute(
        CVD = as.numeric(CVD > 0),
        exposure_raw = as.numeric(exposure_value),
        modifier_risk = as.numeric(.data[[score]] < modifiers$threshold[j]),
        analysis_weight, cycle, sdmvpsu, sdmvstra,
        across(all_of(covariates))
      )
    dat <- dat[complete.cases(dat), ] |> droplevels()

    values <- unique(dat$exposure_raw)
    binary <- grepl("_binary_", chemical) || length(values) == 2
    dat$exposure_model <- if (binary) as.numeric(dat$exposure_raw == max(values))
                          else as.numeric(scale(dat$exposure_raw))
    dat$strata <- interaction(dat$cycle, dat$sdmvstra, drop = TRUE)
    dat$psu <- interaction(dat$cycle, dat$sdmvstra, dat$sdmvpsu, drop = TRUE)

    design <- svydesign(~psu, strata = ~strata, weights = ~analysis_weight,
                        nest = TRUE, data = dat)
    fit <- svyglm(
      reformulate(c("exposure_model * modifier_risk", covariates), "CVD"),
      design, family = quasibinomial()
    )
    co <- summary(fit)$coefficients
    term <- "exposure_model:modifier_risk"
    beta <- co[term, 1]
    se <- co[term, 2]

    model_rows[[k]] <- tibble(
      Pollutant = chemical,
      Modifier = modifiers$Modifier[j],
      n = nrow(dat),
      CVD_events = sum(dat$CVD == 1),
      exposure_scale = if (binary) "binary" else "per SD",
      OR_Interaction_Ratio = exp(beta),
      OR_Interaction_Low = exp(beta - 1.96 * se),
      OR_Interaction_High = exp(beta + 1.96 * se),
      P_Interaction = co[term, 4]
    )

    cutpoint <- if (binary) min(values) else median(dat$exposure_raw)
    dat$chemical_high <- if (binary) dat$exposure_model else as.numeric(dat$exposure_raw >= cutpoint)
    cell_rows[[k]] <- dat |>
      count(chemical_high, modifier_risk, name = "N") |>
      left_join(
        dat |> group_by(chemical_high, modifier_risk) |>
          summarise(CVD_events = sum(CVD == 1), .groups = "drop"),
        by = c("chemical_high", "modifier_risk")
      ) |>
      mutate(Pollutant = chemical, Modifier = modifiers$Modifier[j], .before = 1)
    k <- k + 1
  }
}

results <- bind_rows(model_rows) |>
  mutate(P_FDR_BH_all = p.adjust(P_Interaction, "BH")) |>
  group_by(Modifier) |>
  mutate(P_FDR_BH_within_family = p.adjust(P_Interaction, "BH")) |>
  ungroup() |>
  arrange(P_FDR_BH_all)

write_csv(results, file.path(out_dir, "interaction_primary_survey_weighted_all.csv"), na = "")
write_csv(bind_rows(cell_rows), file.path(out_dir, "interaction_cell_specific_event_counts.csv"), na = "")
