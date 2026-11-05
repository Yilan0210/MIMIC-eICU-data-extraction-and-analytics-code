# Model development and validation (R)

R code for the analyses reported in the paper: development of the 16-predictor model, internal–external cross-validation (IECV) across seven clusters, comparison with existing scores, immune–inflammation–nutrition markers, sensitivity analyses, additional validation, figures, tables and the web calculator.

The scripts read the cleaned data produced by the Python pipeline in this repository (`aecopd_clean_<TAG>.csv`, or the `clean_all` sheet of `aecopd_clean_<TAG>_LOCAL_ONLY.xlsx`). **No data are included.**

- Analysis population: all eligible patients in both databases (`pop == "main"` and `in_s2 == 1`; MIMIC-IV 1,276 + eICU-CRD 3,141 = 4,417 patients, 456 deaths).
- Marker subset: patients with albumin and a differential leukocyte count measured (`in_analysis == 1`; 1,858 patients), used only for the marker analyses.

## Requirements

R 4.5.1 (x86_64-pc-linux-gnu) was used, with mice 3.18.0, rms 8.1.1, Hmisc 5.2.6, metafor 5.2.1, glmnet 5.2, gbm 2.3.1, lme4 2.1.0, ggplot2 4.0.3, patchwork 1.3.2, ggalluvial 0.12.6, ragg 1.4.0, systemfonts 1.2.3, readxl 1.5.0 and writexl 1.5.4.

```r
pkgs <- c("mice", "rms", "Hmisc", "metafor", "glmnet", "gbm", "lme4", "ggplot2", "patchwork", "ggalluvial",
          "ragg", "systemfonts", "readxl", "writexl")
install.packages(pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)])
```

## Running

1. Open `run_all.R` and set `CODE_DIR` (this `analysis` folder), `DATA_DIR` (the folder with the cleaned data) and `CORES`.
2. Test run: `TEST_RUN <- TRUE` (m = 3, B = 20, about 10–20 min; output in `DATA_DIR/v4_out_test`; the numbers are not usable).
3. Final run: `TEST_RUN <- FALSE` (m = 33 imputations, 200 bootstrap samples per imputation, 7-cluster IECV; about 3–5 h with 4 cores; output in `DATA_DIR/v4_out`).
4. Then run `run_A12.R` (sensitivity analysis A12) and `run_extra.R` (additional validation V1–V9).

`run_chain.R` runs all steps in sequence in one job (test runs, final run, A12, additional validation) and stops at the first error. `run_figures.R` redraws the figures from saved results without re-analysis. `run_all.sh` runs the scripts from the command line (`nohup bash run_all.sh > v4_run.log 2>&1 &`; settings through the environment variables `AECOPD_DATA_DIR`, `AECOPD_OUT_DIR`, `AECOPD_CORES`, `AECOPD_M`, `AECOPD_B`).

In RStudio, run the scripts as a Background Job (Background Jobs -> Start Background Job -> R Script = the script -> Working Directory = this folder -> Start). When a script is sourced directly in the console, parallel IECV is switched off.

Sensitivity analysis A10 (90-day death, MIMIC-IV only) needs `aecopd_mimic_death90_<TAG>_LOCAL_ONLY.csv` from `python/export_death90.py` in `DATA_DIR`; without it, A10 is skipped.

## Scripts

| Script | Content | Main output |
|---|---|---|
| `R/00_setup.R` | settings; core model (16 variables, 24 df: 8 spline + 8 linear/binary); imputation, splines, Rubin's rules, performance measures, random-effects pooling, model approximation | – |
| `R/01_describe.R` | checks before modelling: numbers, clusters, missingness, comparison of patients with/without blood gas and markers, Riley sample size | `result/01_describe.xlsx` |
| `R/02_develop.R` | multiple imputation; full model by ML + uniform bootstrap shrinkage; simplified model (model approximation, R² ≥ 0.95); ORs, partial effects, importance, non-linearity tests, stability | `result/02_development.xlsx` |
| `R/03_iecv.R` | 7-cluster IECV with re-imputation, knot placement, estimation and simplification in each fold; paired comparison with scores; risk groups; case mix | `result/03_iecv.xlsx` |
| `R/04_crossdb_update.R` | cross-database validation in both directions (MIMIC-IV -> eICU-CRD with ridge regression); intercept and slope updating within each database | `result/04_crossdb_update.xlsx` |
| `R/05_markers.R` | added value of markers over the core linear predictor (LR test, bootstrap and IECV dC, NRI, net benefit); markers over ICU scores; NLR/PNI constraints and 7 indices; inverse probability weighting; dose–response | `result/05_markers.xlsx` |
| `R/06_sensitivity.R` | sensitivity analyses A1–A14 (A12 via `run_A12.R`); subgroups and fairness | `result/06_sensitivity.xlsx` |
| `R/07_figures_main.R` | main figures and draft legends | `figure/Figure*` |
| `R/08_figures_supp.R` | supplementary figures and draft legends | `figure/eFigure*` |
| `R/09_tables_calculator.R` | tables and the web calculator (full model with splines) | `result/09_tables.xlsx`, `calculator/` |
| `R/10_extra_validation.R` | additional validation V1–V9: temporal validation, hospital benchmarking, hospital type and case definition, case-mix decomposition, value of information, fairness, reclassification, multiverse | `result/10_extra_validation.xlsx` |
| `R/11_figures_extra.R` | figures for V1–V9 | `figure/FigX*` |
| `R/iecv_functions.R` | one validation fold (shared by 03 and 04) | – |
| `R/fig_style.R` | common figure style | – |

## Outputs and data protection

`result/`, `figure/` and `calculator/` contain aggregate results only. `local/` contains patient-level data (imputed datasets, held-out predictions, hospital-level results) and must stay on the analysis server; it must not be uploaded to GitHub or any online service.

## Figure numbers

The scripts use working figure names. They correspond to the figures in the paper as follows.

| Paper | Working figure(s) |
|---|---|
| Figure 1 | `Figure1_study_flow` |
| Figure 2 | `eFigure1_IECV_performance` |
| Figure 3 | `Figure3_transportability` (panel a), `FigX3_casemix_decomposition` (panel b), `FigX1_temporal_validation` (panels c–e) |
| Figure 4 | `Figure4_comparison_scores` (panels a–d), `FigX6_reclassification` (panels e, f) |
| Figure 5 | `Figure6_risk_stratification` (panels a, b), `FigX2_hospital_benchmarking` (panels c–e) |
| Supplementary Figures S1–S15 | `eFigure9_missing`, `Figure2_predictor_effects`, `Figure5_nomogram`, `eFigure2_simplified_vs_full`, `eFigure4_stability`, `FigX9_model_updating`, `FigX8_hospital_type_case_definition`, `eFigure7_ROC`, `FigX4_decision_uncertainty`, `eFigure5_sensitivity`, `FigX7_specification_curve`, `FigX5_fairness_calibration`, `eFigure6_subgroups`, `Figure7_markers`, `eFigure8_markers_IECV` |
