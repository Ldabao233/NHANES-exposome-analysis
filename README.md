# NHANES exposome analysis repository

本仓库整理了论文返修所用的成人 NHANES 2013–2016 分析数据、正式统计代码和主要结果。原始目录未被修改；这里保留的是可复现主分析所需的精简版本。

## Directory structure

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

- `data/processed/exposure_modules/`：16 个污染物模块中间文件，包括 `HM_urine.rds`。
- `data/processed/exwas_analysis_long.rds`：未插补、年龄 ≥20 岁的单污染物长数据。
- `data/processed/mixtures/`：11 个混合暴露组的分析数据。
- `data/processed/phenotype.csv`：生活方式交互分析所需的 LE8 分项评分。
- `data/processed/network/`：已经下载和筛选后的 CTD、SEA、SuperPred、SwissTargetPrediction 及 CVD 基因列表。
- `data/data_dictionary/mixture_groups.csv`：连续编号 C01–C11；`source_group_id` 保留旧分析编号，便于追溯。

原始 NHANES XPT 下载文件没有重复打包。`module_inventory.csv` 提供每个模块的维度、文件大小和 MD5 校验值；`processed_data_manifest.csv` 记录所有分析数据文件的校验值。

## Reproduce the analyses

在仓库根目录运行。也可以把环境变量 `NHANES_REPO` 指向本目录。

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

`00_build_exposure_modules_from_xpt.R` 从 2013–2014 和 2015–2016 的 XPT 文件生成 16 个模块 RDS。脚本中只记录 XPT 文件名；运行前将所需 XPT 放在同一目录，并通过 `XPT_DIR` 指定该目录。尿肌酐校正直接使用两个周期的 `ALB_CR` XPT 文件。

正式参数为：WQS 50 次 repeated holdout、每次 1,000 次 bootstrap、40% 验证集并同时检验两个方向；qgcomp 1,000 次 bootstrap；BKMR 5,000 次迭代，前 2,500 次为 burn-in。

BKMR 固定删除二元检出指标，并运行至少保留 3 个连续暴露的组。

核心统计脚本均按从上到下的线性流程编写，不再包含服务器调度、续跑、复杂错误捕获或冗长诊断函数。绘图脚本仅保留少量主题和多格式导出函数，以避免重复代码。已有 BKMR 诊断结果仍保留在 `results/BKMR/`，但精简版主模型脚本不再重复生成这些诊断文件。

主要 R 包包括 `tidyverse`、`survey`、`gWQS`、`qgcomp`、`bkmr`、`coda`、`mediation`、`ggrepel`、`patchwork`、`ggalluvial`、`igraph`、`svglite` 和 `ragg`。

## Network toxicology

网络毒理部分不包含任何靶点下载或网页查询代码。下载后的靶点表位于 `data/processed/network/`；富集、PPI 和 MCODE 结果位于 `results/network/`；`code/08_network/` 仅保留论文图及敏感性图的绘图代码。
