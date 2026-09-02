# NHANES Exposome Analysis Repository

This repository contains the processed analysis datasets, finalized statistical code, and main results for the adult NHANES 2013–2016 analyses used in the study:

**“Environmental Toxicants Associated with Cardiovascular Disease and a Novel LE8 Indicator Approach for Screening Potential Cardiovascular Toxicants.”**



## Directory Structure

```text
repository/
├── README.md
├── data/
│   ├── processed/
│   │   ├── exposure_modules/
│   │   ├── mixtures/
│   │   └── network/
│   └── data_dictionary/
├── code/
│   ├── 01_data_cleaning/
│   ├── 02_ExWAS/
│   ├── 03_WQS/
│   ├── 04_qgcomp/
│   ├── 05_BKMR/
│   ├── 06_mediation/
│   ├── 07_interaction/
│   └── 08_network/
└── results/
```

## Data

- `data/processed/exposure_modules/`: Intermediate datasets for the 16 environmental exposure modules, including `HM_urine.rds`.
- `data/processed/exwas_analysis_long.rds`: Non-imputed long-format dataset used for single-exposure analyses among adults aged ≥20 years.
- `data/processed/mixtures/`: Analysis datasets for the 11 chemical mixture groups.
- `data/processed/phenotype.csv`: LE8 component scores used in the lifestyle interaction analyses.
- `data/processed/network/`: Downloaded and filtered target information from CTD, SEA, SuperPred, and SwissTargetPrediction, together with the cardiovascular disease (CVD) gene list.
- `data/data_dictionary/mixture_groups.csv`: Chemical mixture groups are numbered consecutively from C01 to C11. The `source_group_id` variable retains the identifiers used in the original analyses for traceability.

The original NHANES XPT files are not redistributed in this repository.

`module_inventory.csv` provides the dimensions, file sizes, and MD5 checksums for each exposure module. `processed_data_manifest.csv` records checksums for all processed analysis datasets.

## Reproducing the Analyses

Run the following commands from the repository root directory. Alternatively, the environment variable `NHANES_REPO` can be set to the repository path.

```bash
XPT_DIR=/path/to/flat_xpt_files RDS_DIR=data/processed/exposure_modules \
  Rscript code/01_data_cleaning/00_build_exposure_modules_from_xpt.R

Rscript code/01_data_cleaning/01_validate_processed_data.R
Rscript code/01_data_cleaning/02_export_qc_tables.R

Rscript code/02_ExWAS/01_run_weighted_exwas.R
Rscript code/02_ExWAS/02_plot_exwas.R

Rscript code/03_WQS/01_run_repeated_holdout_wqs.R

Rscript code/04_qgcomp/01_run_qgcomp.R

Rscript code/05_BKMR/01_run_bkmr.R
Rscript code/05_BKMR/02_plot_bkmr.R results/BKMR/BKMR_overall_risk.csv results/BKMR

Rscript code/06_mediation/01_run_mediation.R

Rscript code/07_interaction/01_run_lifestyle_interactions.R
Rscript code/07_interaction/02_plot_interactions.R

Rscript code/08_network/01_plot_main_figure.R
Rscript code/08_network/02_plot_sensitivity.R
```

The script `00_build_exposure_modules_from_xpt.R` generates 16 exposure-module RDS files from the NHANES 2013–2014 and 2015–2016 XPT files.

Only the required XPT filenames are specified in the script. Before running the script, place all required XPT files in the same directory and specify that directory using `XPT_DIR`.

Urinary creatinine adjustment is performed directly using the `ALB_CR` XPT files from the two NHANES cycles.

# 