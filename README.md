# Analysis Code — Reproduction Guide

Manuscript: *Hemoglobin glycation index and incident stroke in middle-aged and older adults:
prospective evidence from the CHARLS and ELSA cohorts with internal and external consistency
evaluation* (submitted to Diabetes Research and Clinical Practice).

## 1. Data access (raw data are NOT included)

This study is a secondary analysis of two publicly available cohorts whose user agreements
prohibit redistribution of microdata:

- **CHARLS** (China Health and Retirement Longitudinal Study): register at
  https://charls.pku.edu.cn/ ; download Waves 1–6 (2011–2023) questionnaire data and the
  2011/2015 blood assay data.
- **ELSA** (English Longitudinal Study of Ageing): UK Data Service study number **5050**;
  nurse datasets (Waves 2/4/6/8/9) and main interview data (Waves 1–11).

The de-identified **analytic datasets** produced by the extraction scripts
(`analytic_*.csv`) accompany this code package so that every downstream step can be
verified without access to the raw data. Set the path variables at the top of each script
(`OUT_DIR`, default `D:/CHARLS/HGI_stroke`) to the folder holding these CSV files.

## 2. Environment

- R ≥ 4.0 (analyses run under R 4.6 on Windows; scripts are platform-independent)
- Packages: `survival`, `cmprsk`, `mice`, `tableone`, `dplyr`, `splines`, `metafor`
- Random seed wherever resampling is involved (multiple imputation, bootstrap): **20260903**

## 3. Pipeline overview

### 3.1 CHARLS main analysis (`R_scripts/charls/`)

| Order | Script | Purpose | Key outputs |
|---|---|---|---|
| 1 | `CHARLS_HGI_S01_data_extraction_screening.R` | Variable extraction from raw CHARLS; baseline analytic set; HGI derivation | `analytic_hgi_baseline.csv`, `analytic_delta_hgi.csv` |
| 2 | `CHARLS_HGI_W6_S01_data_extension.R` | Extension to Wave 6 (2023) follow-up; event-time rules fu_v2 / fu_v2b | `analytic_hgi_baseline_w6.csv` |
| 3 | `CHARLS_HGI_W6_S02_main_analysis.R` | Main analysis: Cox M1–M3, prespecified HbA1c key test, RCS, Fine–Gray, subgroups | main-set estimates (Table 2; Figs 2–3) |
| 4 | `CHARLS_HGI_S04_sensitivity_S1_S6.R` | Prespecified sensitivity analyses S1–S6 (incl. multiple imputation, m = 20) | Fig. S2b, Table S8 |
| 5 | `CHARLS_HGI_2015_internal_validation.R` | 2015 landmark internal-consistency set (independent blood batch, post-2015 window) | landmark estimates (Table 2) |
| 6 | `CHARLS_HGI_W6_S06_figure_table_patch_20261119.R` | Table/figure numeric patches | — |
| 7 | `CHARLS_HGI_W6_S07_HGI_SHR_mutual_20261119.R` | Mutual HGI × SHR adjustment (S7) | Table S8 |
| 8 | `CHARLS_HGI_W6_S08_Evalue_interval_censoring_20261119.R` | E-values; interval-censoring refits (S8) | Table S8 |
| 9 | `CHARLS_HGI_W6_S09_revision_supplementary_20261120.R` | Post hoc analyses P-series (kidney function, early events, etc.) | Table S8 |
| 10 | `CHARLS_HGI_W6_S10_2015_covariate_refit_nonDM_20261120.R` | Landmark refit with 2015-measured covariates; non-diabetic layer | Tables S3, S8 |
| 11 | `CHARLS_HGI_W6_S11_2015_definition_switch_20261120.R` | Definition-switch values (S11) | Table S8 |
| 12 | `CHARLS_HGI_W6_S12_nonDM_stratification_landmark_20261120.R` | Non-diabetic reconciliation; landmark extensions (S12) | Table S8 |
| 13 | `CHARLS_HGI_W6_S13_treatment_strata_eGFR_early_20261120.R` | Diabetes-treatment strata; eGFR; early-event exclusion | Table S8 |
| 14 | `CHARLS_HGI_W6_S14_landmark_FG_absolute_risk_Evalue_20261121.R` | Landmark Fine–Gray; absolute risk (Table S5); E-values; interaction power (MDE 1.81) | Table S5, S8 |
| 15 | `CHARLS_HGI_W6_S15_equivalence_TyG_mutual_20261121.R` | Algebraic-equivalence verification (HGI+HbA1c ≡ FPG+HbA1c); TyG mutual adjustment; AIC ladder | Table S6 |

### 3.2 ELSA external evaluation (`R_scripts/elsa/`)

| Order | Script | Purpose |
|---|---|---|
| 1 | `ELSA_S01_build_analytic_set.R` | Build ELSA Wave 2 analytic set (fasting FPG + HbA1c, stroke-free, ≥ 1 follow-up) |
| 2 | `ELSA_S02_main_analysis.R` | Cox ladder, key test, Fine–Gray, non-diabetic layer (Table 2; Fig. 2) |
| 3 | `ELSA_S03_pooled_baseline.R` | Pooled Waves 2/4/6 stratified sensitivity analysis (Table S1) |

### 3.3 Cross-cohort meta-analysis (`R_scripts/meta/`)

- `CHARLS_ELSA_cross_cohort_meta.R`: fixed/random-effects pooling of the two independent sets
  (primary pooled estimate 1.363), the three sets, and the three non-diabetic sets (Fig. 2).

### 3.4 Figures and tables

- `figure_scripts/fig1_strobe.R` … `fig5_robust.R`: Figures 1–3 and Supplementary
  Figures S1–S2 (600-dpi TIFF).
- `table_scripts/build_tables.R`: all manuscript and supplementary tables (CSV).

## 4. Verification anchors

Every number in the manuscript was machine-checked against the frozen result files
(audit report: `DRCP版数字核对说明_20261121.md`). Primary anchors:

- Main set: n = 8,684; 852 events; model 3 HR 1.352 (1.108–1.650); key test HR 1.337 (1.039–1.719).
- Non-diabetic: n = 7,283; 681 events; model 3 HR 1.475 (1.147–1.896).
- Pooled CHARLS 2011 + ELSA: HR 1.363 (1.136–1.635), I² = 0%.

## 5. Notes

- A few legacy pipeline scripts (notably the early CHARLS extraction and patch scripts)
  contain Chinese-language comments; all executable code, variable names, and output
  files are in English, and every script was syntax-checked under R 4.6.2 before release.
- Path defaults at the top of some scripts point to local Windows folders (the ELSA
  scripts use a Chinese folder name, `ELSA验证`); set `DATA_DIR` / `OUT_DIR` (ELSA:
  `D` / `O`) to your own folders before running.

Contact: Jing Tian, jingscorpio@163.com
