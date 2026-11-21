####################################################################################
# Revision supplement 2 (2026-11-20), user-question driven
#  S10a: CHARLS 2015 landmark refit with 2015-measured covariates (R1-M4 resolution)
#  S10b: Non-diabetic full analysis ladder (CHARLS main set + ELSA W2)
# Output: results/S10_*.csv and console log
####################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(haven); library(survival); library(dplyr) })
BASE <- "D:/CHARLS/HGI_stroke/freeze_20261120"
OUT  <- file.path(BASE, "results")
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
adj3 <- paste(cov_m3, collapse = " + ")
adj2 <- "age + sex + urban + edu"
to_w1_id <- function(id) { id <- as.character(id)
  hit <- nchar(id) == 12 & substr(id, 10, 10) == "0"
  id[hit] <- paste0(substr(id[hit], 1, 9), substr(id[hit], 11, 12)); id }

grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(term = term, HR = round(ci[,1],3), lo = round(ci[,3],3),
             hi = round(ci[,4],3), P = signif(p,3))
}

## =================================================================================
## S10a: 2015 landmark with 2015-measured covariates
## =================================================================================
cat("========== S10a: 2015 covariate refit ==========\n")

# -- W3 disease status: inspect z vs da semantics -----------------------------------
hs <- read_dta("D:/CHARLS/2015/Health_Status_and_Functioning.dta")
for (x in c("1_","2_","3_","7_")) {
  z <- hs[[paste0("zda007_", x)]]; a <- hs[[paste0("da007_", x)]]
  cat("cond", x, " z==1:", sum(z == 1, na.rm = TRUE), " da==1:", sum(a == 1, na.rm = TRUE),
      " z|da:", sum((!is.na(z) & z == 1) | (!is.na(a) & a == 1)), "\n")
}
mk <- function(z, a) as.integer((!is.na(z) & z == 1) | (!is.na(a) & a == 1))
w3_dis <- data.frame(ID = to_w1_id(hs$ID),
                     htn15    = mk(hs$zda007_1_, hs$da007_1_),
                     dyslip15 = mk(hs$zda007_2_, hs$da007_2_),
                     dm15     = mk(hs$zda007_3_, hs$da007_3_),
                     heart15  = mk(hs$zda007_7_, hs$da007_7_))

# -- W3 SBP (three readings, same range filter as 2011 pipeline) ---------------------
bm <- read_dta("D:/CHARLS/2015/Biomarker.dta")
sbp_mat <- sapply(bm[c("qa003","qa007","qa011")], function(x) {
  x <- as.numeric(x); x[x < 60 | x > 260] <- NA; x
})
w3_bp <- data.frame(ID = to_w1_id(bm$ID), sbp15 = rowMeans(sbp_mat, na.rm = TRUE))
w3_bp$sbp15[is.nan(w3_bp$sbp15)] <- NA

# -- W3 blood covariates --------------------------------------------------------------
bl <- read_dta("D:/CHARLS/2015/Blood.dta")
w3_bl <- data.frame(ID = to_w1_id(bl$ID),
                    crp15 = as.numeric(bl$bl_crp), tg15 = as.numeric(bl$bl_tg),
                    hdl15 = as.numeric(bl$bl_hdl), glu15 = as.numeric(bl$bl_glu),
                    hba1c15 = as.numeric(bl$bl_hbalc))

# -- rebuild the landmark set exactly as the frozen script -----------------------------
dd <- read.csv(file.path(BASE, "data/charls/analytic_delta_hgi.csv"),
               colClasses = c(ID = "character"))
d  <- read.csv(file.path(BASE, "data/charls/analytic_hgi_baseline_w6.csv"),
               colClasses = c(ID = "character"))
adjv <- strsplit(gsub(" ", "", adj3), "\\+")[[1]]
dd <- dd %>%
  select(-any_of(intersect(adjv, setdiff(names(dd), "ID")))) %>%
  left_join(d %>% select(ID, stroke_w6, last_contact, w6_report_yr, w6_dxyear,
                         all_of(adjv)), by = "ID") %>%
  left_join(w3_bl %>% select(ID, glu15, hba1c15), by = "ID") %>%
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
  filter(!is.na(t15) & t15 > 0)

# rename 2015 blood to match frozen names used downstream
dd <- dd %>% rename(glu3 = glu15, hba1c3 = hba1c15)

CUT15 <- as.numeric(quantile(dd$hgi_2015, 0.9, na.rm = TRUE))
dd <- dd %>% mutate(hgi15_hi_p90 = as.integer(hgi_2015 >= CUT15))
cat("landmark set rebuilt:", nrow(dd), "n; events:", sum(dd$ev15),
    " (frozen: 6,730 / 716 pre-cc)\n")

# -- anchor: frozen carried-covariate model --------------------------------------------
dmain_old <- dd %>% filter(complete.cases(across(all_of(c(cov_m3, "hgi_2015", "ev15", "t15")))))
S <- "Surv(t15, ev15)"
adj3_old <- paste(c("age15","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl"), collapse = " + ")
m_old <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 +", adj3_old)), data = dmain_old)
cat("anchor (2011-carried covariates): n =", nrow(dmain_old), "; events =", sum(dmain_old$ev15), "; ")
print(grab(m_old, "hgi15_hi_p90"))

# -- merge W3 covariates ---------------------------------------------------------------
dd2 <- dd %>%
  select(-any_of(c("htn","dm","dyslip","heart","sbp","crp","tg","hdl"))) %>%
  left_join(w3_dis, by = "ID") %>% left_join(w3_bp, by = "ID") %>% left_join(w3_bl, by = "ID")
cov15 <- c("age15","sex","urban","edu","htn15","dm15","dyslip15","heart15","sbp15","crp15","tg15","hdl15")
adj15 <- paste(cov15, collapse = " + ")
dmain15 <- dd2 %>% filter(complete.cases(across(all_of(c(cov15, "hgi_2015", "hba1c3", "ev15", "t15")))))
cat("2015-covariate complete-case set:", nrow(dmain15), "n; events", sum(dmain15$ev15), "\n")

# prevalence comparison
cmp <- data.frame(
  var = c("htn","dm","dyslip","heart"),
  carried_2011 = c(mean(dmain_old$htn), mean(dmain_old$dm), mean(dmain_old$dyslip), mean(dmain_old$heart)),
  measured_2015 = c(mean(dmain15$htn15), mean(dmain15$dm15), mean(dmain15$dyslip15), mean(dmain15$heart15)))
print(cmp)

m15_1 <- coxph(as.formula(paste(S, "~ hgi15_hi_p90")), data = dmain15)
m15_2 <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 + age15 + sex + urban + edu")), data = dmain15)
m15_3 <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 +", adj15)), data = dmain15)
m15_k <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 + hba1c3 +", adj15)), data = dmain15)
res10a <- rbind(
  cbind(model = "M1 unadjusted", grab(m15_1, "hgi15_hi_p90")),
  cbind(model = "M2 demographic", grab(m15_2, "hgi15_hi_p90")),
  cbind(model = "M3 full (2015 covariates)", grab(m15_3, "hgi15_hi_p90")),
  cbind(model = "key test + HbA1c-2015", grab(m15_k, "hgi15_hi_p90")),
  cbind(model = "key test: HbA1c term", grab(m15_k, "hba1c3")))
print(res10a)
write.csv(res10a, file.path(OUT, "S10a_landmark_2015covariates_refit.csv"), row.names = FALSE)
write.csv(cmp, file.path(OUT, "S10a_covariate_vintage_comparison.csv"), row.names = FALSE)

## =================================================================================
## S10b: non-diabetic ladder (CHARLS main + ELSA W2)
## =================================================================================
cat("\n========== S10b: non-diabetic analyses ==========\n")
dd0 <- d %>% mutate(hgi_hi_p90 = as.integer(hgi >= quantile(hgi, 0.9, na.rm = TRUE)))
dmain <- dd0 %>% filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c",
                                                       "incident_stroke6", "fu_v2")))))
ndm <- dmain %>% filter(diabetes == 0)
cat("CHARLS non-diabetic set:", nrow(ndm), "n; events", sum(ndm$incident_stroke6),
    " (of main 8,684 / 852)\n")
adj3_ndm <- paste(setdiff(cov_m3, "dm"), collapse = " + ")
FS <- "Surv(fu_v2, incident_stroke6)"
n1 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90")), data = ndm)
n2 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj2)), data = ndm)
n3 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj3_ndm)), data = ndm)
nk <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 + hba1c +", adj3_ndm)), data = ndm)
res10b1 <- rbind(
  cbind(cohort = "CHARLS non-DM", model = "M1", grab(n1, "hgi_hi_p90")),
  cbind(cohort = "CHARLS non-DM", model = "M2", grab(n2, "hgi_hi_p90")),
  cbind(cohort = "CHARLS non-DM", model = "M3", grab(n3, "hgi_hi_p90")),
  cbind(cohort = "CHARLS non-DM", model = "key test", grab(nk, "hgi_hi_p90")),
  cbind(cohort = "CHARLS non-DM", model = "key: HbA1c", grab(nk, "hba1c")))
print(res10b1)
write.csv(res10b1, file.path(OUT, "S10b_CHARLS_nondiabetic_ladder.csv"), row.names = FALSE)

# ELSA non-diabetic
e <- read.csv(file.path(BASE, "data/elsa/analytic_elsa_hgi_baseline_w2.csv"))
e <- e %>% mutate(ndm = as.integer(dm == 0 & hba1c < 6.5 & fglu < 7.0))
en <- e %>% filter(ndm == 1)
cat("ELSA non-diabetic set:", nrow(en), "n; events", sum(en$incident_stroke), "\n")
cov_e <- c("age","sex","urban","edu","htn","dyslip","heart","sbp","crp","tg","hdl")
adj_e <- paste(cov_e, collapse = " + ")
SE <- "Surv(fu_mid, incident_stroke)"
e1 <- coxph(as.formula(paste(SE, "~ hgi_hi_p90")), data = en)
e2 <- coxph(as.formula(paste(SE, "~ hgi_hi_p90 + age + sex + urban + edu")), data = en)
e3 <- coxph(as.formula(paste(SE, "~ hgi_hi_p90 +", adj_e)), data = en)
ek <- coxph(as.formula(paste(SE, "~ hgi_hi_p90 + hba1c +", adj_e)), data = en)
res10b2 <- rbind(
  cbind(cohort = "ELSA non-DM", model = "M1", grab(e1, "hgi_hi_p90")),
  cbind(cohort = "ELSA non-DM", model = "M2", grab(e2, "hgi_hi_p90")),
  cbind(cohort = "ELSA non-DM", model = "M3", grab(e3, "hgi_hi_p90")),
  cbind(cohort = "ELSA non-DM", model = "key test", grab(ek, "hgi_hi_p90")),
  cbind(cohort = "ELSA non-DM", model = "key: HbA1c", grab(ek, "hba1c")))
print(res10b2)
write.csv(res10b2, file.path(OUT, "S10b_ELSA_nondiabetic_ladder.csv"), row.names = FALSE)
cat("\nSaved S10 outputs to", OUT, "\n")
