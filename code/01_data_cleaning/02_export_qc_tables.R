#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

root <- normalizePath(Sys.getenv("NHANES_REPO", "."))
groups <- read_csv(file.path(root, "data", "data_dictionary", "mixture_groups.csv"),
                   show_col_types = FALSE)$category_group_id
data_dir <- file.path(root, "data", "processed", "mixtures")
rows <- list()

for (i in seq_along(groups)) {
  object <- readRDS(file.path(data_dir, paste0(groups[i], "_mixture_data.rds")))
  dat <- object$data
  cvd_ok <- complete.cases(dat[, c("CVD", object$exposures, object$cvd_covariates)])
  le8_ok <- complete.cases(dat[, c("LE8", object$exposures, object$le8_covariates)])

  rows[[i]] <- tibble(
    category_group_id = groups[i],
    family = first(dat$family),
    n_start = nrow(dat),
    n_cvd_final = sum(cvd_ok),
    cvd_events = sum(dat$CVD[cvd_ok] == 1),
    n_le8_final = sum(le8_ok)
  )
}

write_csv(bind_rows(rows), file.path(root, "results", "mixture_sample_sizes.csv"), na = "")
