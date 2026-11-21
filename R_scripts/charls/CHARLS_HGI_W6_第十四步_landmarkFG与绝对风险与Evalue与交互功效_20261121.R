####################################################################################
# S14 (2026-11-21): revision round 5 — reviewer-driven additions
#  A. CHARLS 2015 landmark Fine-Gray (fills Table 2 "—"; reviewer R1-M3)
#  B. Absolute risk by HGI group (high vs low) for six sets (reviewer R2-M1):
#     crude cumulative incidence, person-year rates, RD and rate difference + 95% CI
#  C. E-values for the non-diabetic layer (reviewer R1-m2)
#  D. Diabetes-interaction power analysis (user question 2026-11-21):
#     interaction coefficient SE -> minimum detectable ratio-of-HRs (80% power),
#     post-hoc power at the observed ratio, cross-cohort pooled interaction,
#     and events multiplier needed for 80% power
# Canonical sets identical to S11/S12 (anchors re-verified by stopifnot).
# Output: results/S14_landmark_FG.csv, S14_absolute_risk.csv,
#         S14_nonDM_Evalues.csv, S14_interaction_power.csv + console
####################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(haven); library(survival); library(dplyr) })
BASE <- "D:/CHARLS/HGI_stroke/freeze_20261120"
OUT  <- file.path(BASE, "results")
to_w1_id <- function(id) { id <- as.character(id)
  hit <- nchar(id) == 12 & substr(id, 10, 10) == "0"
  id[hit] <- paste0(substr(id[hit], 1, 9), substr(id[hit], 11, 12)); id }
grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(term = term, HR = round(ci[,1],3), lo = round(ci[,3],3),
             hi = round(ci[,4],3), P = signif(p,3))
}
evalue <- function(hr) ifelse(hr >= 1, hr + sqrt(hr * (hr - 1)),
                              1/hr + sqrt(1/hr * (1/hr - 1)))

## ================================================================================
## 0. Rebuild frozen sets (identical to S12)
## ================================================================================
d <- read.csv(file.path(BASE, "data/charls/analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
dmain <- d %>% mutate(hgi_hi_p90 = as.integer(hgi >= CUT_P90)) %>%
  filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c",
                                        "incident_stroke6", "fu_v2")))))
stopifnot(nrow(dmain) == 8684, sum(dmain$incident_stroke6) == 852)

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
dd <- read.csv(file.path(BASE, "data/charls/analytic_delta_hgi.csv"),
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
cov15 <- c("age15","sex","urban","edu","htn15","dm15","dyslip15","heart15",
           "sbp15","crp15","tg15","hdl15")
lm15 <- dd %>% mutate(hgi15_hi_p90 = as.integer(hgi_2015 >= CUT15)) %>%
  filter(complete.cases(across(all_of(c(cov15, "hgi_2015", "hba1c3", "ev15", "t15")))))
stopifnot(nrow(lm15) == 6581, sum(lm15$ev15) == 706)
cat("anchors ok: 2011 cc 8684/852; landmark cc 6581/706\n")

adj2 <- "age + sex + urban + edu"
adj3_ndm <- paste(setdiff(cov_m3, "dm"), collapse = " + ")
adj15_ndm <- paste(setdiff(cov15, "dm15"), collapse = " + ")
FS <- "Surv(fu_v2, incident_stroke6)"
S  <- "Surv(t15, ev15)"
ndm    <- dmain %>% filter(diabetes == 0)
stopifnot(nrow(ndm) == 7283, sum(ndm$incident_stroke6) == 681)
lm_ndm <- lm15 %>% filter(dm15 == 0 & hba1c3 < 6.5 & glu3 < 126)
stopifnot(nrow(lm_ndm) == 5267, sum(lm_ndm$ev15) == 521)

e <- read.csv(file.path(BASE, "data/elsa/analytic_elsa_hgi_baseline_w2.csv"))
cov_e <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
SE <- "Surv(fu_mid, incident_stroke)"
en <- e %>% filter(dm == 0 & hba1c < 6.5 & fglu < 7.0)
stopifnot(nrow(e) == 2577, nrow(en) == 2474)

## ================================================================================
## A. Landmark Fine-Gray (death as competing event, 2015-2023 window)
## ================================================================================
cat("\n==== A. landmark Fine-Gray ====\n")
lmfg <- lm15 %>% mutate(
  death_t15 = death_time - 4,                       # death_time is vs 2011
  has_death = !is.na(death6) & death6 == 1 & !is.na(death_t15) & death_t15 > 0,
  has_stroke = ev15 == 1,
  fg_event = case_when(
    has_stroke & (!has_death | t15 <= death_t15) ~ 1,
    has_death & (!has_stroke | death_t15 < t15)  ~ 2,
    TRUE ~ 0),
  fg_time = case_when(fg_event == 1 ~ t15,
                      fg_event == 2 ~ death_t15,
                      TRUE ~ t15),
  fg_time = pmax(fg_time, 0.25))
cat("landmark FG status (0/1/2):\n"); print(table(lmfg$fg_event))
fg_lm <- NA
if (requireNamespace("cmprsk", quietly = TRUE)) {
  library(cmprsk)
  X <- model.matrix(~ hgi15_hi_p90 + age15 + sex + urban + edu + htn15 + dm15 +
                      dyslip15 + heart15 + sbp15 + crp15 + tg15 + hdl15,
                    data = lmfg)[, -1]
  fm <- crr(lmfg$fg_time, lmfg$fg_event, X, failcode = 1, cencode = 0)
  b <- summary(fm)$coef["hgi15_hi_p90", ]
  fg_lm <- data.frame(set = "CHARLS 2015 landmark", model = "Fine-Gray",
                      term = "hgi15_hi_p90",
                      HR = round(exp(b[1]),3), lo = round(exp(b[1]-1.96*b[3]),3),
                      hi = round(exp(b[1]+1.96*b[3]),3), P = signif(b[5],3),
                      n = nrow(lmfg), events = sum(lmfg$fg_event == 1),
                      competing_deaths = sum(lmfg$fg_event == 2))
  print(fg_lm)
  write.csv(fg_lm, file.path(OUT, "S14_landmark_FG.csv"), row.names = FALSE)
  # non-DM landmark FG (for completeness of the two-layer table)
  lf2 <- lmfg %>% filter(dm15 == 0 & hba1c3 < 6.5 & glu3 < 126)
  X2 <- model.matrix(~ hgi15_hi_p90 + age15 + sex + urban + edu + htn15 +
                       dyslip15 + heart15 + sbp15 + crp15 + tg15 + hdl15,
                     data = lf2)[, -1]
  fm2 <- crr(lf2$fg_time, lf2$fg_event, X2, failcode = 1, cencode = 0)
  b2 <- summary(fm2)$coef["hgi15_hi_p90", ]
  fg_lm2 <- data.frame(set = "CHARLS 2015 landmark non-DM", model = "Fine-Gray",
                       term = "hgi15_hi_p90",
                       HR = round(exp(b2[1]),3), lo = round(exp(b2[1]-1.96*b2[3]),3),
                       hi = round(exp(b2[1]+1.96*b2[3]),3), P = signif(b2[5],3),
                       n = nrow(lf2), events = sum(lf2$fg_event == 1),
                       competing_deaths = sum(lf2$fg_event == 2))
  print(fg_lm2)
  write.csv(rbind(fg_lm, fg_lm2), file.path(OUT, "S14_landmark_FG.csv"),
            row.names = FALSE)
} else cat("cmprsk not available; skipped\n")

## ================================================================================
## B. Absolute risk: high vs low HGI, six sets
## ================================================================================
cat("\n==== B. absolute risk ====\n")
absrisk <- function(set, event, time, grp) {
  z <- qnorm(0.975)
  one <- function(e, t, g) {
    n <- length(e); ev <- sum(e); py <- sum(t)
    ci <- ev / n; rt <- ev / py * 1000
    c(n = n, events = ev, py = round(py, 0), ci_pct = round(ci * 100, 1),
      rate1000 = round(rt, 1))
  }
  hi <- one(event[grp == 1], time[grp == 1], NULL)
  lo <- one(event[grp == 0], time[grp == 0], NULL)
  # cumulative incidence difference (pp), Wald
  p1 <- hi["events"] / hi["n"]; p0 <- lo["events"] / lo["n"]
  rd <- (p1 - p0) * 100
  se_rd <- sqrt(p1 * (1 - p1) / hi["n"] + p0 * (1 - p0) / lo["n"]) * 100
  # rate difference per 1000 py, Wald on rates
  r1 <- hi["events"] / hi["py"] * 1000; r0 <- lo["events"] / lo["py"] * 1000
  rdiff <- r1 - r0
  se_rdiff <- sqrt(hi["events"] / hi["py"]^2 + lo["events"] / lo["py"]^2) * 1000
  data.frame(set = set,
             n_high = hi["n"], ev_high = hi["events"], py_high = hi["py"],
             ci_high_pct = hi["ci_pct"], rate_high = hi["rate1000"],
             n_low = lo["n"], ev_low = lo["events"], py_low = lo["py"],
             ci_low_pct = lo["ci_pct"], rate_low = lo["rate1000"],
             RD_pp = round(rd, 1), RD_lo = round(rd - z * se_rd, 1),
             RD_hi = round(rd + z * se_rd, 1),
             rate_diff_1000py = round(rdiff, 1),
             rateD_lo = round(rdiff - z * se_rdiff, 1),
             rateD_hi = round(rdiff + z * se_rdiff, 1),
             median_fu = round(median(time), 1))
}
ar <- rbind(
  absrisk("CHARLS 2011 main", dmain$incident_stroke6, dmain$fu_v2, dmain$hgi_hi_p90),
  absrisk("CHARLS 2011 non-DM", ndm$incident_stroke6, ndm$fu_v2, ndm$hgi_hi_p90),
  absrisk("CHARLS 2015 landmark", lm15$ev15, lm15$t15, lm15$hgi15_hi_p90),
  absrisk("CHARLS 2015 landmark non-DM", lm_ndm$ev15, lm_ndm$t15, lm_ndm$hgi15_hi_p90),
  absrisk("ELSA W2", e$incident_stroke, e$fu_mid, e$hgi_hi_p90),
  absrisk("ELSA W2 non-DM", en$incident_stroke, en$fu_mid, en$hgi_hi_p90))
print(ar)
write.csv(ar, file.path(OUT, "S14_absolute_risk.csv"), row.names = FALSE)

## ================================================================================
## C. E-values for the non-diabetic layer (frozen estimates)
## ================================================================================
cat("\n==== C. non-DM E-values ====\n")
ev_rows <- data.frame(
  estimate = c("CHARLS 2011 non-DM M3", "CHARLS 2011 non-DM key",
               "CHARLS 2015 landmark non-DM M3", "ELSA W2 non-DM M3",
               "non-DM pooled (3 sets)"),
  HR  = c(1.475, 1.524, 1.492, 1.511, 1.484),
  lo  = c(1.147, 1.127, 0.942, 0.918, 1.213))
ev_rows$E_point <- round(evalue(ev_rows$HR), 2)
ev_rows$E_CI <- ifelse(ev_rows$lo <= 1, 1.00, round(evalue(ev_rows$lo), 2))
print(ev_rows)
write.csv(ev_rows, file.path(OUT, "S14_nonDM_Evalues.csv"), row.names = FALSE)

## ================================================================================
## D. Diabetes-interaction power analysis
## ================================================================================
cat("\n==== D. interaction power ====\n")
# CHARLS: full-set model with interaction (reproduces P_interaction = 0.398)
m_int_c <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 * diabetes +",
                                  paste(setdiff(cov_m3, "dm"), collapse = " + "))),
                 data = dmain)
b_c <- summary(m_int_c)$coefficients["hgi_hi_p90:diabetes", ]
# ELSA: full-set model with interaction (diab_grp == dm==1; reproduces P = 0.536)
e$diab_grp <- as.integer(!is.na(e$dm) & e$dm == 1)
m_int_e <- coxph(as.formula(paste(SE, "~ hgi_hi_p90 * diab_grp +",
                                  paste(setdiff(cov_e, "dm"), collapse = " + "))),
                 data = e)
b_e <- summary(m_int_e)$coefficients["hgi_hi_p90:diab_grp", ]

pw <- function(b, se, label, events) {
  za <- qnorm(0.975); zb <- qnorm(0.8)
  ratio <- exp(b)
  mde <- exp((za + zb) * se)                 # min detectable ratio-of-HRs, 80% power
  posthoc <- pnorm(abs(b) / se - za)         # power at observed ratio (approx)
  mult <- ((za + zb) * se / abs(b))^2        # events multiplier for 80% at observed
  data.frame(analysis = label, ratio_HR = round(ratio, 3),
             P_interaction = signif(2 * pnorm(-abs(b / se)), 3),
             se = round(se, 3),
             MDE_ratio_80pct = round(mde, 2),
             posthoc_power = round(posthoc, 3),
             events_now = events,
             events_needed_80pct = round(events * mult, 0))
}
d1 <- pw(b_c[1], b_c[3], "CHARLS 2011 (diabetes interaction)", sum(dmain$incident_stroke6))
d2 <- pw(b_e[1], b_e[3], "ELSA W2 (diabetes interaction)", sum(e$incident_stroke))
# pooled interaction coefficient (inverse variance)
w <- 1 / c(b_c[3], b_e[3])^2
b_p <- sum(w * c(b_c[1], b_e[1])) / sum(w); se_p <- sqrt(1 / sum(w))
d3 <- pw(b_p, se_p, "Pooled CHARLS+ELSA interaction",
         sum(dmain$incident_stroke6) + sum(e$incident_stroke))
pwres <- rbind(d1, d2, d3)
print(pwres)
write.csv(pwres, file.path(OUT, "S14_interaction_power.csv"), row.names = FALSE)

cat("\nS14 done. saved: S14_landmark_FG.csv, S14_absolute_risk.csv,",
    "S14_nonDM_Evalues.csv, S14_interaction_power.csv\n")
