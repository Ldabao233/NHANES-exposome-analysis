#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

repo_root <- normalizePath(Sys.getenv("NHANES_REPO", "."), mustWork = TRUE)
data_dir <- file.path(repo_root, "data", "processed")
long <- readRDS(file.path(data_dir, "exwas_analysis_long.rds"))

required <- c(
  "SEQN", "cycle", "exposure", "exposure_value", "analysis_weight",
  "sdmvpsu", "sdmvstra", "age", "CVD", "LE8"
)
stopifnot(!length(setdiff(required, names(long))))
stopifnot(all(long$age >= 20, na.rm = TRUE))

groups <- read_csv(
  file.path(repo_root, "data", "data_dictionary", "mixture_groups.csv"),
  show_col_types = FALSE
)

for (group_id in groups$category_group_id) {
  path <- file.path(data_dir, "mixtures", paste0(group_id, "_mixture_data.rds"))
  stopifnot(file.exists(path))
  object <- readRDS(path)
  stopifnot(all(c("data", "exposures", "cvd_covariates", "le8_covariates") %in% names(object)))
  stopifnot(all(object$exposures %in% names(object$data)))
  stopifnot(all(object$data$age >= 20, na.rm = TRUE))
}

module_files <- list.files(file.path(data_dir, "exposure_modules"), pattern = "[.]rds$")
stopifnot(length(module_files) == 16L)
inventory <- read_csv(
  file.path(repo_root, "data", "data_dictionary", "module_inventory.csv"),
  show_col_types = FALSE
)
module_paths <- file.path(data_dir, "exposure_modules", inventory$file)
stopifnot(all(unname(tools::md5sum(module_paths)) == inventory$md5))

cat("Adult participants:", n_distinct(long$SEQN), "\n")
cat("Exposures:", n_distinct(long$exposure), "\n")
cat("Mixture groups:", nrow(groups), "\n")
cat("Exposure modules:", length(module_files), "\n")
