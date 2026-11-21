####################################################################################
# Integrated tables for the HGI x stroke manuscript remake (frozen 2026-11-20)
# Table 1: baseline characteristics of the three analytic sets (by HGI group)
# Table 2: main results across sets + cross-cohort pooling
# Table 3: five-marker head-to-head
# Table 4: robustness battery
# Table S1: ELSA pooled-baseline (W2/W4/W6) transparency report
# Table S2: ELSA subgroups
# Table S3: CHARLS 2015 landmark extra analyses (trait-marker evidence)
# Table S4: STROBE participant flow numbers
# All effect sizes are hard-coded from the frozen result files (results/);
# Table 1 is recomputed from the frozen analytic csvs with the exact frozen
# set definitions. Output: freeze_20261120/tables/*.csv
####################################################################################
rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(dplyr); library(haven) })
FRZ <- "D:/CHARLS/HGI_stroke/freeze_20261120"
OUT <- file.path(FRZ, "tables"); dir.create(OUT, showWarnings = FALSE)

f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
fmt_msd  <- function(x, d = 1) sprintf("%s (%s)", f2(mean(x, na.rm = TRUE), d),
                                       f2(sd(x, na.rm = TRUE), d))
fmt_np   <- function(x) sprintf("%d (%.1f)", sum(x == 1, na.rm = TRUE),
                                  mean(x == 1, na.rm = TRUE) * 100)
fmt_mdiqr <- function(x, d = 2) sprintf("%s [%s, %s]", f2(median(x, na.rm = TRUE), d),
                                        f2(quantile(x, .25, na.rm = TRUE), d),
                                        f2(quantile(x, .75, na.rm = TRUE), d))
# DRCP journal rule: report HbA1c in NGSP (%) AND IFCC (mmol/mol, no decimals)
# conversion: mmol/mol = 10.93 * % - 23.5 (mean); SD scales by 10.93 only
fmt_ifcc <- function(x) sprintf("%d (%d)", round(10.93 * mean(x, na.rm = TRUE) - 23.5),
                                round(10.93 * sd(x, na.rm = TRUE)))
smd_cont <- function(x, g) { x0 <- x[g == 0]; x1 <- x[g == 1]
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) /
    sqrt((sd(x1, na.rm = TRUE)^2 + sd(x0, na.rm = TRUE)^2) / 2) }
smd_bin <- function(x, g) { p0 <- mean(x[g == 0] == 1, na.rm = TRUE)
  p1 <- mean(x[g == 1] == 1, na.rm = TRUE)
  (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2) }

## ================================================================================
## 1. Rebuild the three analytic sets exactly as frozen
## ================================================================================
## ---- CHARLS 2011 main set (W6, n = 8,684) ---------------------------------------
d <- read.csv(file.path(FRZ, "data/charls/analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
c11 <- d %>%
  mutate(hgi_hi = as.integer(hgi >= CUT_P90)) %>%
  filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c",
                                        "incident_stroke6", "fu_v2")))))
stopifnot(nrow(c11) == 8684, sum(c11$incident_stroke6) == 852)

## ---- CHARLS 2015 landmark set (2015-measured covariates, cc n = 6,581) ------------
to_w1_id <- function(id) { id <- as.character(id)
  hit <- nchar(id) == 12 & substr(id, 10, 10) == "0"
  id[hit] <- paste0(substr(id[hit], 1, 9), substr(id[hit], 11, 12)); id }
hs <- read_dta("D:/CHARLS/2015/Health_Status_and_Functioning.dta")
mk <- function(z, a) as.integer((!is.na(z) & z == 1) | (!is.na(a) & a == 1))
w3_dis <- data.frame(ID = to_w1_id(hs$ID),
                     htn15 = mk(hs$zda007_1_, hs$da007_1_),
                     dyslip15 = mk(hs$zda007_2_, hs$da007_2_),
                     dm15 = mk(hs$zda007_3_, hs$da007_3_),
                     heart15 = mk(hs$zda007_7_, hs$da007_7_))
bm <- read_dta("D:/CHARLS/2015/Biomarker.dta")
sbp_mat <- sapply(bm[c("qa003","qa007","qa011")], function(x) {
  x <- as.numeric(x); x[x < 60 | x > 260] <- NA; x })
w3_bp <- data.frame(ID = to_w1_id(bm$ID), sbp15 = rowMeans(sbp_mat, na.rm = TRUE))
w3_bp$sbp15[is.nan(w3_bp$sbp15)] <- NA
bl <- read_dta("D:/CHARLS/2015/Blood.dta")
w3_bl <- data.frame(ID = to_w1_id(bl$ID),
                    crp15 = as.numeric(bl$bl_crp), tg15 = as.numeric(bl$bl_tg),
                    hdl15 = as.numeric(bl$bl_hdl), glu3 = as.numeric(bl$bl_glu),
                    hba1c3 = as.numeric(bl$bl_hbalc))
dd <- read.csv(file.path(FRZ, "data/charls/analytic_delta_hgi.csv"),
               colClasses = c(ID = "character"))
dd <- dd %>%
  select(-any_of(intersect(cov_m3, setdiff(names(dd), "ID")))) %>%
  left_join(d %>% select(ID, stroke_w6, last_contact, w6_report_yr, w6_dxyear,
                         death6, death_time, all_of(cov_m3)), by = "ID") %>%
  left_join(w3_bl, by = "ID") %>%
  filter(is.na(first_stroke_wave) | first_stroke_wave >= 4) %>%
  mutate(
    age15 = age + 4,
    ev15 = as.integer((!is.na(stroke_w4) & stroke_w4 == 1) |
                      (!is.na(stroke_w5) & stroke_w5 == 1) |
                      (!is.na(stroke_w6) & stroke_w6 == 1)),
    rep_yr = case_when(!is.na(stroke_w4) & stroke_w4 == 1 ~ 2018,
                       !is.na(stroke_w5) & stroke_w5 == 1 ~ 2020,
                       !is.na(stroke_w6) & stroke_w6 == 1 ~ w6_report_yr),
    prev_yr = case_when(!is.na(stroke_w4) & stroke_w4 == 1 ~ 2015,
                        !is.na(stroke_w5) & stroke_w5 == 1 ~ 2018,
                        !is.na(stroke_w6) & stroke_w6 == 1 ~ 2020),
    t15 = if_else(ev15 == 1,
                  ifelse(!is.na(w6_dxyear) & !is.na(stroke_w6) & stroke_w6 == 1 &
                           w6_dxyear > 2015 & w6_dxyear <= rep_yr,
                         w6_dxyear - 2015, (prev_yr + rep_yr) / 2 - 2015),
                  last_contact - 2015)) %>%
  filter(!is.na(t15) & t15 > 0) %>%
  select(-any_of(c("htn","dm","dyslip","heart","sbp","crp","tg","hdl"))) %>%
  left_join(w3_dis, by = "ID") %>% left_join(w3_bp, by = "ID")
CUT15 <- as.numeric(quantile(dd$hgi_2015, 0.9, na.rm = TRUE))
dd <- dd %>% mutate(hgi_hi = as.integer(hgi_2015 >= CUT15))
cov15 <- c("age15","sex","urban","edu","htn15","dm15","dyslip15","heart15",
           "sbp15","crp15","tg15","hdl15")
c15 <- dd %>% filter(complete.cases(across(all_of(
  c(cov15, "hgi_2015", "ev15", "t15")))))
stopifnot(nrow(dd) == 6730, nrow(c15) == 6581, sum(c15$ev15) == 706)

## ---- ELSA W2 main set (n = 2,577; csv is already complete-case) -------------------
els <- read.csv(file.path(FRZ, "data/elsa/analytic_elsa_hgi_baseline_w2.csv")) %>%
  rename(hgi_hi = hgi_hi_p90)
stopifnot(nrow(els) == 2577, sum(els$incident_stroke) == 176)

## ================================================================================
## 2. Table 1 — baseline characteristics (Overall / Low / High / SMD per set)
## ================================================================================
## spec: label, type (m = mean SD; p = n %; q = median IQR),
##       acc = list of 3 accessors (C2011, C2015, ELSA), digits
rows_spec <- list(
  list("Age, years", "m", list(~ .x$age, ~ .x$age15, ~ .x$age), 1),
  list("Female, n (%)", "p", list(~ as.integer(.x$sex == 2), ~ as.integer(.x$sex == 2),
       ~ as.integer(.x$sex == 1))),
  list("Urban residence, n (%)", "p", list(~ .x$urban, ~ .x$urban, ~ .x$urban)),
  list("Highest education category, n (%)a", "p",
       list(~ as.integer(.x$edu >= 5), ~ as.integer(.x$edu >= 5),
            ~ as.integer(.x$edu == 3))),
  list("Married or partnered, n (%)b", "p",
       list(~ as.integer(.x$marital %in% c(1, 2)),
            ~ as.integer(.x$marital %in% c(1, 2)), NULL)),
  list("Body-mass index, kg/m2", "m", list(~ .x$bmi, ~ .x$bmi, NULL), 1),
  list("Systolic blood pressure, mm Hg", "m", list(~ .x$sbp, ~ .x$sbp15, ~ .x$sbp), 0),
  list("Fasting plasma glucose, mmol/L", "m",
       list(~ .x$glu / 18.0, ~ .x$glu3 / 18.0, ~ .x$fglu), 2),
  list("HbA1c, %", "m", list(~ .x$hba1c, ~ .x$hba1c3, ~ .x$hba1c), 2),
  list("HbA1c, mmol/mol", "ifcc", list(~ .x$hba1c, ~ .x$hba1c3, ~ .x$hba1c)),
  list("Triglycerides, mmol/L", "m",
       list(~ .x$tg / 88.57, ~ .x$tg15 / 88.57, ~ .x$trig), 2),
  list("HDL cholesterol, mmol/L", "m",
       list(~ .x$hdl / 38.67, ~ .x$hdl15 / 38.67, ~ .x$hdl / 38.67), 2),
  list("C-reactive protein, mg/L", "q", list(~ .x$crp, ~ .x$crp15, ~ .x$hscrp), 2),
  list("Hypertension, n (%)", "p", list(~ .x$htn, ~ .x$htn15, ~ .x$htn)),
  list("Diabetes, n (%)", "p", list(~ .x$dm, ~ .x$dm15, ~ .x$dm)),
  list("Dyslipidemia, n (%)", "p", list(~ .x$dyslip, ~ .x$dyslip15, ~ .x$dyslip)),
  list("Heart disease, n (%)", "p", list(~ .x$heart, ~ .x$heart15, ~ .x$heart)),
  list("HGIc", "m", list(~ .x$hgi, ~ .x$hgi_2015, ~ .x$hgi), 2),
  list("Anemia, n (%)d", "p", list(~ .x$anemia, NULL, NULL))
)

summ_one <- function(dat, sp, set_idx) {
  acc <- sp[[3]][[set_idx]]
  if (is.null(acc)) return(c(Overall = "\u2014", Low = "\u2014", High = "\u2014", SMD = "\u2014"))
  g <- dat$hgi_hi
  x <- eval(acc[[2]], list(.x = dat))
  typ <- sp[[2]]; ddig <- if (length(sp) >= 4) sp[[4]] else 1
  if (typ == "m") {
    c(fmt_msd(x, ddig), fmt_msd(x[g == 0], ddig), fmt_msd(x[g == 1], ddig),
      f2(smd_cont(x, g), 2))
  } else if (typ == "p") {
    c(fmt_np(x), fmt_np(x[g == 0]), fmt_np(x[g == 1]), f2(smd_bin(x, g), 2))
  } else if (typ == "ifcc") {
    c(fmt_ifcc(x), fmt_ifcc(x[g == 0]), fmt_ifcc(x[g == 1]), "\u2014")
  } else {
    c(fmt_mdiqr(x, ddig), fmt_mdiqr(x[g == 0], ddig), fmt_mdiqr(x[g == 1], ddig),
      f2(smd_cont(x, g), 2))
  }
}
summ_set <- function(dat, spec, set_idx)
  as.data.frame(t(sapply(spec, summ_one, dat = dat, set_idx = set_idx))) |>
    setNames(c("Overall", "Low", "High", "SMD"))

s11 <- summ_set(c11, rows_spec, 1); s15 <- summ_set(c15, rows_spec, 2)
sel <- summ_set(els, rows_spec, 3)
tab1 <- data.frame(
  Characteristic = sapply(rows_spec, `[[`, 1),
  C2011_Overall = s11$Overall, C2011_Low = s11$Low, C2011_High = s11$High, C2011_SMD = s11$SMD,
  C2015_Overall = s15$Overall, C2015_Low = s15$Low, C2015_High = s15$High, C2015_SMD = s15$SMD,
  ELSA_Overall  = sel$Overall, ELSA_Low  = sel$Low, ELSA_High  = sel$High, ELSA_SMD  = sel$SMD,
  check.names = FALSE)
miqr <- function(x, d = 1) sprintf("%s [%s, %s]", f2(median(x, na.rm = TRUE), d),
                                   f2(quantile(x, .25, na.rm = TRUE), d),
                                   f2(quantile(x, .75, na.rm = TRUE), d))
hdr_vals <- list(
  list("Participants, n",
       nrow(c11), sum(c11$hgi_hi == 0), sum(c11$hgi_hi == 1),
       nrow(c15), sum(c15$hgi_hi == 0), sum(c15$hgi_hi == 1),
       nrow(els), sum(els$hgi_hi == 0), sum(els$hgi_hi == 1)),
  list("Incident stroke, n (%)",
       sprintf("%d (%.1f)", sum(c11$incident_stroke6), mean(c11$incident_stroke6) * 100),
       sprintf("%d (%.1f)", sum(c11$incident_stroke6[c11$hgi_hi == 0]),
               mean(c11$incident_stroke6[c11$hgi_hi == 0]) * 100),
       sprintf("%d (%.1f)", sum(c11$incident_stroke6[c11$hgi_hi == 1]),
               mean(c11$incident_stroke6[c11$hgi_hi == 1]) * 100),
       sprintf("%d (%.1f)", sum(c15$ev15), mean(c15$ev15) * 100),
       sprintf("%d (%.1f)", sum(c15$ev15[c15$hgi_hi == 0]),
               mean(c15$ev15[c15$hgi_hi == 0]) * 100),
       sprintf("%d (%.1f)", sum(c15$ev15[c15$hgi_hi == 1]),
               mean(c15$ev15[c15$hgi_hi == 1]) * 100),
       sprintf("%d (%.1f)", sum(els$incident_stroke), mean(els$incident_stroke) * 100),
       sprintf("%d (%.1f)", sum(els$incident_stroke[els$hgi_hi == 0]),
               mean(els$incident_stroke[els$hgi_hi == 0]) * 100),
       sprintf("%d (%.1f)", sum(els$incident_stroke[els$hgi_hi == 1]),
               mean(els$incident_stroke[els$hgi_hi == 1]) * 100)),
  list("Deaths, n",
       sum(c11$death6 == 1, na.rm = TRUE), "", "",
       sum(c15$death6 == 1 & (c15$death_time - 4) > 0, na.rm = TRUE), "", "",
       sum(els$died == TRUE, na.rm = TRUE), "", ""),
  list("Follow-up, median [IQR], years",
       miqr(c11$fu_v2), "", "", miqr(c15$t15), "", "", miqr(els$fu_years_total), "", ""),
  list("HGI P90 cut-off",
       f2(CUT_P90, 3), "", "", f2(CUT15, 3), "", "",
       "0.455", "", ""))  # ELSA threshold from frozen build log (fitted on the
                           # analysis-eligible set n = 3,062), not on the cc set
hdr <- as.data.frame(t(sapply(hdr_vals, function(r) {
  c(r[[1]], r[[2]], r[[3]], r[[4]], "\u2014", r[[5]], r[[6]], r[[7]], "\u2014",
    r[[8]], r[[9]], r[[10]], "\u2014")
})))
names(hdr) <- names(tab1)
tab1 <- rbind(hdr, tab1)
write.csv(tab1, file.path(OUT, "Table1_baseline_integrated.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## ================================================================================
## 3. Tables 2–4 and supplementary (frozen numbers)
## ================================================================================
hr <- function(HR, lo, hi, p) sprintf("%.3f (%.3f\u2013%.3f), P %s", HR, lo, hi, p)
tab2 <- data.frame(Analysis = c(
  "Model 1 (unadjusted)",
  "Model 2 (age, sex, urban/rural, education)",
  "Model 3 (fully adjusted, 12 covariates)",
  "Model 3 + HbA1c (key test)",
  "Fine\u2013Gray competing risk (death)",
  "Pooled: CHARLS 2011 + CHARLS 2015 (internal, random effects)",
  "Pooled: CHARLS 2011 + ELSA W2 (external; fixed-effect shown, fixed = random at I2 = 0%)",
  "Pooled: all three sets (random effects)",
  "Pooled non-diabetic: CHARLS 2011 + CHARLS 2015 (fixed = random, I2 = 0%)",
  "Pooled non-diabetic: CHARLS 2011 + ELSA W2 (fixed = random, I2 = 0%)",
  "Pooled non-diabetic: all three sets (fixed = random, I2 = 0%)"),
  `CHARLS 2011 (n=8,684)` = c(
    hr(1.582, 1.303, 1.920, "< 0.001"), hr(1.553, 1.279, 1.885, "< 0.001"),
    hr(1.352, 1.108, 1.650, "= 0.003"), hr(1.337, 1.039, 1.719, "= 0.024"),
    hr(1.345, 1.102, 1.640, "= 0.004"), "", "", "", "", "", ""),
  `CHARLS 2011 non-diabetic (n=7,283)` = c(
    hr(1.582, 1.231, 2.032, "< 0.001"), hr(1.552, 1.208, 1.995, "< 0.001"),
    hr(1.475, 1.147, 1.896, "= 0.002"), hr(1.524, 1.127, 2.060, "= 0.006"),
    hr(1.479, 1.158, 1.889, "= 0.002"), "", "", "", "", "", ""),
  `CHARLS 2015 landmark (n=6,581)` = c(
    hr(1.494, 1.207, 1.850, "< 0.001"), hr(1.463, 1.181, 1.812, "< 0.001"),
    hr(1.152, 0.920, 1.442, "= 0.219"), hr(1.179, 0.889, 1.563, "= 0.253"),
    hr(1.139, 0.912, 1.422, "= 0.250"), "", "", "", "", "", ""),
  `CHARLS 2015 landmark non-diabetic (n=5,267)` = c(
    hr(1.476, 0.933, 2.333, "= 0.096"), hr(1.518, 0.960, 2.401, "= 0.074"),
    hr(1.492, 0.942, 2.362, "= 0.088"), hr(1.536, 0.954, 2.474, "= 0.077"),
    hr(1.465, 0.926, 2.317, "= 0.100"), "", "", "", "", "", ""),
  `ELSA W2 (n=2,577)` = c(
    hr(1.813, 1.187, 2.770, "= 0.006"), hr(1.681, 1.098, 2.572, "= 0.017"),
    hr(1.420, 0.903, 2.235, "= 0.129"), hr(1.440, 0.845, 2.455, "= 0.180"),
    hr(1.326, 0.830, 2.119, "= 0.240"), "", "", "", "", "", ""),
  `ELSA W2 non-diabetic (n=2,474)` = c(
    hr(1.742, 1.067, 2.844, "= 0.026"), hr(1.655, 1.011, 2.707, "= 0.045"),
    hr(1.511, 0.918, 2.485, "= 0.104"), hr(1.425, 0.799, 2.539, "= 0.230"),
    "\u2014", "", "", "", "", "", ""),
  `Pooled HR (95% CI); P; I2` = c("", "", "", "", "",
    "1.259 (1.077\u20131.472); P = 0.004; I2 = 8%",
    "1.363 (1.136\u20131.635); P < 0.001; I2 = 0%",
    "1.275 (1.107\u20131.469); P < 0.001; I2 = 0%",
    "1.479 (1.186\u20131.844); P < 0.001; I2 = 0%",
    "1.482 (1.184\u20131.855); P < 0.001; I2 = 0%",
    "1.484 (1.213\u20131.816); P < 0.001; I2 = 0%"),
  check.names = FALSE)
write.csv(tab2, file.path(OUT, "Table2_main_results.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

tab3 <- data.frame(Marker = c("HGI", "HbA1c", "FPG", "SHR", "TyG"),
  `CHARLS C-index` = c("0.6809 (reference)", "0.6814", "0.6804", "0.6801", "0.6829"),
  `CHARLS delta-C vs HGI (95% CI)` = c("\u2014", "+0.0005 (\u22120.0010 to 0.0024)",
    "\u22120.0005 (\u22120.0035 to 0.0024)", "\u22120.0009 (\u22120.0030 to 0.0007)",
    "+0.0022 (\u22120.0022 to 0.0072)"),
  `CHARLS per-SD HR` = c(1.042, 1.065, 1.046, 0.994, 1.272),
  `ELSA C-index` = c(0.7037, 0.7022, 0.7022, "0.7044 (highest)", 0.7023),
  `ELSA per-SD HR` = c("1.118 (only marker > 1)", 1.051, 0.903, 0.850, 0.825),
  check.names = FALSE)
write.csv(tab3, file.path(OUT, "Table3_marker_headtohead.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

tab4 <- data.frame(Analysis = c(
  "Main analysis (Model 3)",
  "Fine\u2013Gray competing risk (death)",
  "S1  HGI equation refit in non-diabetic subset",
  "S2  Published (NHANES) HGI equation",
  "S3  Excluding anaemia (low-HGI arm HR 0.96, P = 0.73)",
  "S4a Alternative follow-up rule (report wave)",
  "S4b Diagnosis-year refined follow-up",
  "S5a Without SBP covariate",
  "S5b Multiple imputation (m = 20)",
  "S7  Mutual adjustment for stress hyperglycaemia ratio",
  "S8  Interval-censoring model (Weibull)",
  "S9  Joint model: high HGI + FPG (post hoc)",
  "P1  Non-diabetic analysis, within-subset P90 threshold (post hoc)",
  "P2  Non-diabetic analysis, narrow self-report definition (post hoc)",
  "P3  Non-diabetic analysis, equation refit within non-diabetic subset (post hoc)",
  "P4  Untreated/undiagnosed diabetes stratum (post hoc; treated n = 54, too small)",
  "P5  Excluding CKD (eGFR < 60, CKD-EPI 2021) (post hoc)",
  "P6  Model 3 + continuous eGFR (post hoc)",
  "P7  Excluding strokes within first 2 years (post hoc)",
  "P8  Normoglycaemia only (excluding prediabetes; within-subset P90) (post hoc)",
  "P9  Mutual adjustment for TyG index, non-diabetic set (post hoc)",
  "E-value (point estimate; CI limit)",
  "E-value, non-diabetic set (point estimate; CI limit)"),
  `CHARLS 2011` = c(
    hr(1.352, 1.108, 1.650, "= 0.003"), hr(1.345, 1.102, 1.640, "= 0.004"),
    hr(1.275, 1.039, 1.563, "= 0.020"), hr(1.362, 1.117, 1.660, "= 0.002"),
    hr(1.280, 1.036, 1.582, "= 0.022"), hr(1.333, 1.092, 1.627, "= 0.005"),
    hr(1.352, 1.108, 1.650, "= 0.003"), hr(1.346, 1.118, 1.620, "= 0.002"),
    hr(1.306, 1.089, 1.567, "= 0.004"), hr(1.397, 1.129, 1.729, "= 0.002"),
    hr(1.352, 1.108, 1.650, "= 0.003"), hr(1.334, 1.091, 1.632, "= 0.005"),
    hr(1.280, 1.023, 1.603, "= 0.031"), hr(1.431, 1.154, 1.775, "= 0.001"),
    hr(1.381, 1.110, 1.719, "= 0.004"), hr(1.305, 0.931, 1.831, "= 0.123"),
    hr(1.367, 1.119, 1.670, "= 0.002"), hr(1.357, 1.112, 1.656, "= 0.003"),
    hr(1.317, 1.074, 1.616, "= 0.008"), hr(1.389, 0.992, 1.944, "= 0.055"),
    hr(1.531, 1.190, 1.969, "< 0.001"),
    "2.04; 1.46", "2.31; 1.56 (three-set pooled non-diabetic: 2.33; 1.72)"),
  `ELSA W2` = c(
    hr(1.420, 0.903, 2.235, "= 0.129"), hr(1.326, 0.830, 2.119, "= 0.240"),
    "\u2014", "\u2014", "\u2014", "\u2014", "\u2014", "\u2014", "\u2014",
    hr(1.164, 0.693, 1.954, "= 0.567"), hr(1.433, 0.910, 2.255, "= 0.120"),
    "\u2014", "\u2014", "\u2014", "\u2014", "\u2014", "\u2014", "\u2014",
    "\u2014", "\u2014", hr(1.432, 0.864, 2.373, "= 0.163"), "2.19; 1.00 (CI crosses null)", "2.39; 1.00 (CI crosses null)"),
  check.names = FALSE)
write.csv(tab4, file.path(OUT, "Table4_robustness.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

## ---- Table S1: ELSA pooled baseline ---------------------------------------------------
tabS1 <- data.frame(Analysis = c(
  "Pooled W2/W4/W6 (n = 5,455; 314 events) \u2014 M1",
  "Pooled \u2014 M2", "Pooled \u2014 M3 (stratified by entry wave)",
  "Pooled \u2014 M3 + HbA1c: HGI high", "Pooled \u2014 M3 + HbA1c: HbA1c",
  "Pooled \u2014 per-SD HGI",
  "W2 subset (n = 2,577; 176 events) \u2014 M2", "W2 \u2014 M3",
  "W4 subset (n = 1,907; 106 events) \u2014 M2", "W4 \u2014 M3",
  "W6 subset (n = 971; 32 events) \u2014 M2", "W6 \u2014 M3",
  "W2+W4 (n = 4,484; 282 events) \u2014 M1", "W2+W4 \u2014 M2", "W2+W4 \u2014 M3",
  "Pooled Fine\u2013Gray (death data incomplete for W4/W6)"),
  `HR (95% CI), P` = c(
    hr(1.426, 1.016, 2.002, "= 0.041"), hr(1.308, 0.930, 1.839, "= 0.123"),
    hr(1.120, 0.782, 1.605, "= 0.536"), hr(1.009, 0.668, 1.523, "= 0.967"),
    hr(1.135, 0.903, 1.427, "= 0.278"),
    "1.086 (0.974\u20131.212), P = 0.140",
    hr(1.681, 1.098, 2.572, "= 0.017"), hr(1.420, 0.903, 2.235, "= 0.129"),
    hr(0.710, 0.344, 1.463, "= 0.353"), hr(0.619, 0.293, 1.309, "= 0.209"),
    hr(1.596, 0.599, 4.250, "= 0.350"), hr(1.181, 0.409, 3.413, "= 0.758"),
    hr(1.372, 0.954, 1.973, "= 0.088"), hr(1.276, 0.886, 1.838, "= 0.191"),
    hr(1.104, 0.751, 1.622, "= 0.615"), hr(1.071, 0.751, 1.526, "= 0.706")),
  check.names = FALSE)
write.csv(tabS1, file.path(OUT, "TableS1_ELSA_pooled.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

## ---- Table S2: ELSA subgroups ---------------------------------------------------------
tabS2 <- data.frame(
  Subgroup = rep(c("Diabetes", "Sex", "Age", "Residence"), each = 2),
  Level = c("No", "Yes", "Male", "Female", "< 65 y", "\u2265 65 y", "Rural", "Urban"),
  n = c(2520, 57, 1203, 1374, 1509, 1068, 706, 1871),
  Events = c(169, 7, 96, 80, 78, 98, 52, 124),
  `HR (95% CI)` = c("1.482 (0.931\u20132.361)", "1.045 (0.120\u20139.065)",
                    "1.114 (0.595\u20132.086)", "1.674 (0.859\u20133.260)",
                    "1.589 (0.826\u20133.058)", "1.299 (0.683\u20132.472)",
                    "1.657 (0.669\u20134.105)", "1.420 (0.835\u20132.414)"),
  P = c(0.097, 0.968, 0.735, 0.130, 0.166, 0.426, 0.275, 0.196),
  `P interaction` = rep(c(0.536, 0.843, 0.819, 0.620), each = 2),
  check.names = FALSE)
write.csv(tabS2, file.path(OUT, "TableS2_ELSA_subgroups.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

## ---- Table S3: 2015 landmark extras ----------------------------------------------------
cc_delta <- dd %>% filter(!is.na(hgi) & !is.na(hgi_2015))
cor_hgi <- cor(cc_delta$hgi, cc_delta$hgi_2015)
ret <- mean(cc_delta$hgi_2015[cc_delta$hgi >= CUT_P90] >= CUT15, na.rm = TRUE)
tabS3 <- data.frame(Analysis = c(
  "Per-SD HGI-2015 (Model 3)",
  "RCS non-linearity (likelihood-ratio test), P",
  "HGI-2015 high vs low, adjusted for continuous HGI-2011",
  "Continuous HGI-2011 in the same model",
  "Pearson correlation, HGI-2011 vs HGI-2015",
  "Retention in the high-HGI group, 2011 \u2192 2015, %",
  "High-HGI threshold (P90), 2011",
  "High-HGI threshold (P90), 2015",
  "Sensitivity: Model 3 with 2011-carried covariates (former primary)",
  "Sensitivity: key test (HbA1c-adjusted) with 2011-carried covariates"),
  Estimate = c(
    hr(0.989, 0.921, 1.062, "= 0.758"), "0.202",
    hr(1.138, 0.902, 1.436, "= 0.276"), hr(1.025, 0.909, 1.157, "= 0.687"),
    f2(cor_hgi, 2), f2(ret * 100, 1), f2(CUT_P90, 4), f2(CUT15, 4),
    hr(1.117, 0.874, 1.428, "= 0.376"), hr(1.236, 0.908, 1.682, "= 0.178")),
  check.names = FALSE)
write.csv(tabS3, file.path(OUT, "TableS3_landmark_extra.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

## ---- Table S4: STROBE flow ---------------------------------------------------------------
fl_c <- read.csv(file.path(FRZ, "results/charls/strobe_flow.csv"))
tabS4 <- data.frame(
  Set = c(rep("CHARLS 2011", nrow(fl_c) + 1), rep("CHARLS 2015 landmark", 3),
          rep("ELSA W2", 7)),
  Step = c(as.character(fl_c$step),
           "Complete-case main analysis set",
           "Stroke-free at 2015 landmark with 2015 blood draw",
           "Complete-case (12 covariates measured in 2015 + HGI-2015)",
           "Landmark main analysis set",
           "Harmonised ELSA total", "W2 fasting blood subsample (FPG + HbA1c)",
           "Age \u2265 50 at W2", "No prevalent stroke at W2",
           "\u2265 1 interview after W2", "Complete-case main analysis set",
           "Incident strokes W3\u2013W11"),
  n = c(fl_c$n, 8684, 6730, 6581, 6581, 21679, 3428, 3428, 3349, 3062, 2577, 176))
write.csv(tabS4, file.path(OUT, "TableS4_strobe_flow.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

## ---- Table S5: absolute risk by HGI group (from results/S14_absolute_risk.csv) ------
ar <- read.csv(file.path(FRZ, "results/S14_absolute_risk.csv"))
tabS5 <- data.frame(
  Set = ar$set,
  `High HGI: n` = ar$n_high, `High HGI: events` = ar$ev_high,
  `High HGI: cumulative incidence, %` = ar$ci_high_pct,
  `High HGI: rate per 1000 person-years` = ar$rate_high,
  `Low HGI: n` = ar$n_low, `Low HGI: events` = ar$ev_low,
  `Low HGI: cumulative incidence, %` = ar$ci_low_pct,
  `Low HGI: rate per 1000 person-years` = ar$rate_low,
  `Risk difference, percentage points (95% CI)` =
    sprintf("%.1f (%.1f to %.1f)", ar$RD_pp, ar$RD_lo, ar$RD_hi),
  `Rate difference per 1000 person-years (95% CI)` =
    sprintf("%.1f (%.1f to %.1f)", ar$rate_diff_1000py, ar$rateD_lo, ar$rateD_hi),
  `Median follow-up, years` = ar$median_fu,
  check.names = FALSE)
write.csv(tabS5, file.path(OUT, "TableS5_absolute_risk.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

## ---- Table S6: algebraic equivalence check + AIC ladder (S15) -----------------------
eqv <- read.csv(file.path(FRZ, "results/S15_equivalence.csv"))
aic <- read.csv(file.path(FRZ, "results/S15_AIC_ladder.csv"))
trm <- read.csv(file.path(FRZ, "results/S15_AIC_ladder_terms.csv"))
trm$str <- sprintf("%s: HR %.3f (%.3f–%.3f), P %s", trm$term, trm$HR,
                   trm$lo, trm$hi, ifelse(trm$P < 0.001, "< 0.001",
                                          paste("=", trm$P)))
tabS6 <- data.frame(
  Section = c(rep("Algebraic equivalence (continuous specification)", nrow(eqv)),
              rep("AIC ladder (binary high-HGI exposure)", nrow(aic)),
              rep("Model term estimates", nrow(trm))),
  Item = c(eqv$item, aic$model, paste(trm$model, "-", trm$term)),
  Value = c(as.character(eqv$value),
            sprintf("AIC %.1f (dAIC vs base %.2f)", aic$AIC, aic$dAIC_vs_base),
            trm$str),
  check.names = FALSE)
write.csv(tabS6, file.path(OUT, "TableS6_equivalence_AIC.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

cat("cor(HGI-2011, HGI-2015) =", round(cor_hgi, 3), "; retention =",
    round(ret * 100, 1), "%\n")
cat("CUT_P90(2011) =", round(CUT_P90, 4), "; CUT15 =", round(CUT15, 4), "\n")
cat("Tables written to", OUT, "\nDONE\n")
