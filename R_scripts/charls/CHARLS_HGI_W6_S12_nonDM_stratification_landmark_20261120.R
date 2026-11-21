####################################################################################
# S12 (2026-11-20): non-diabetic layer — reconciliation, landmark extension, pooling
#  1. Anchor-reproduce frozen S10b (CHARLS n=7,283/681, M3 1.475, key 1.524;
#     ELSA n=2,474/164, M3 1.511)
#  2. Reconcile with external AI numbers (n=7,144, M3 1.281; refit A2 1.357/key 1.464):
#     variants on the FROZEN cc set — threshold source, diabetes definition,
#     equation source — plus collinearity diagnostics for within-subset refit
#  3. NEW: CHARLS 2015 landmark non-diabetic ladder (2015-measured covariates)
#  4. Non-diabetic cross-set pooling (two-CHARLS, 2011+ELSA, three-set)
#  5. Normoglycemia sensitivity (exclude prediabetes), CHARLS 2011
# Canonical (frozen) definition: broad diabetes exclusion; full-population HGI
# equation; full-population P90 threshold. Output: results/S12_*.csv
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
pool <- function(est, lo, hi) {
  y <- log(est); se <- (log(hi) - log(lo)) / (2 * 1.959964)
  w <- 1 / se^2
  yF <- sum(w * y) / sum(w); seF <- sqrt(1 / sum(w))
  Q <- sum(w * (y - yF)^2); dfr <- length(y) - 1
  I2 <- max(0, (Q - dfr) / Q) * 100
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - dfr) / C)
  wR <- 1 / (se^2 + tau2)
  yR <- sum(wR * y) / sum(wR); seR <- sqrt(1 / sum(wR))
  z <- qnorm(0.975)
  data.frame(fixed_HR = exp(yF), fixed_lo = exp(yF - z*seF), fixed_hi = exp(yF + z*seF),
             fixed_P = 2*pnorm(-abs(yF/seF)),
             re_HR = exp(yR), re_lo = exp(yR - z*seR), re_hi = exp(yR + z*seR),
             re_P = 2*pnorm(-abs(yR/seR)), I2 = round(I2, 1))
}

## ================================================================================
## 0. Rebuild frozen sets (2011 main cc; 2015 landmark cc with 2015 covariates)
## ================================================================================
d <- read.csv(file.path(BASE, "data/charls/analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
dmain <- d %>% mutate(hgi_hi_p90 = as.integer(hgi >= CUT_P90)) %>%
  filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c",
                                        "incident_stroke6", "fu_v2")))))
stopifnot(nrow(dmain) == 8684, sum(dmain$incident_stroke6) == 852)

## landmark rebuild (identical to S10a/S11)
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
                         all_of(cov_m3)), by = "ID") %>%
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
cat("anchors ok: 2011 cc 8684/852; landmark cc 6581/706; CUTs",
    round(CUT_P90,4), round(CUT15,4), "\n")

out <- list()

## ================================================================================
## 1. Anchor: frozen CHARLS 2011 non-DM (canonical)
## ================================================================================
adj2 <- "age + sex + urban + edu"
adj3_ndm <- paste(setdiff(cov_m3, "dm"), collapse = " + ")
FS <- "Surv(fu_v2, incident_stroke6)"
ndm <- dmain %>% filter(diabetes == 0)
stopifnot(nrow(ndm) == 7283, sum(ndm$incident_stroke6) == 681)
a3 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj3_ndm)), data = ndm)
ak <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 + hba1c +", adj3_ndm)), data = ndm)
cat("anchor CHARLS non-DM M3:", grab(a3, "hgi_hi_p90")$HR,
    "; key:", grab(ak, "hgi_hi_p90")$HR, "\n")

## ================================================================================
## 2. Reconciliation variants on the FROZEN cc set
## ================================================================================
## V2: within-subset P90 (external AI choice)
ndm <- ndm %>% mutate(h_within = as.integer(hgi >= quantile(hgi, 0.9, na.rm = TRUE)))
v2 <- coxph(as.formula(paste(FS, "~ h_within +", adj3_ndm)), data = ndm)
## V3: narrow diabetes (self-report dm only), canonical full-P90
ndm_narrow <- dmain %>% filter(dm == 0)
v3 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj3_ndm)), data = ndm_narrow)
## V4: within-subset refit equation + within-subset P90 (external AI A2)
eq_ndm <- lm(hba1c ~ glu, data = ndm)
ndm <- ndm %>% mutate(hgi_refit = hba1c - predict(eq_ndm, newdata = .),
                      h2 = as.integer(hgi_refit >= quantile(hgi_refit, 0.9, na.rm = TRUE)))
v4  <- coxph(as.formula(paste(FS, "~ h2 +", adj3_ndm)), data = ndm)
v4k <- coxph(as.formula(paste(FS, "~ h2 + hba1c +", adj3_ndm)), data = ndm)
## collinearity diagnostics for the refit exposure
cor_refit_a1c <- cor(ndm$hgi_refit, ndm$hba1c)
cor_orig_a1c  <- cor(ndm$hgi, ndm$hba1c)
cat("non-DM eq: R2 =", round(summary(eq_ndm)$r.squared, 3),
    "; cor(hgi_refit, HbA1c) =", round(cor_refit_a1c, 3),
    "; cor(hgi_orig, HbA1c) =", round(cor_orig_a1c, 3), "\n")
recon <- rbind(
  cbind(variant = "V1 canonical: broad DM excl, full-pop P90 (frozen)",
        n = nrow(ndm), events = sum(ndm$incident_stroke6), grab(a3, "hgi_hi_p90")),
  cbind(variant = "V2 broad DM excl, within-subset P90",
        n = nrow(ndm), events = sum(ndm$incident_stroke6), grab(v2, "h_within")),
  cbind(variant = "V3 narrow DM excl (self-report only), full-pop P90",
        n = nrow(ndm_narrow), events = sum(ndm_narrow$incident_stroke6),
        grab(v3, "hgi_hi_p90")),
  cbind(variant = "V4 within-subset refit eq + within-subset P90",
        n = nrow(ndm), events = sum(ndm$incident_stroke6), grab(v4, "h2")),
  cbind(variant = "V4 key test (+ HbA1c)",
        n = nrow(ndm), events = sum(ndm$incident_stroke6), grab(v4k, "h2")))
print(recon)
out$recon <- recon
out$collinear <- data.frame(
  metric = c("within-nonDM equation R2", "cor(refit HGI, HbA1c)",
             "cor(canonical HGI, HbA1c)"),
  value = c(round(summary(eq_ndm)$r.squared, 3), round(cor_refit_a1c, 3),
            round(cor_orig_a1c, 3)))

## ================================================================================
## 3. NEW: 2015 landmark non-DM ladder (2015 covariates; full-landmark CUT15)
## ================================================================================
lm_ndm <- lm15 %>% filter(dm15 == 0 & hba1c3 < 6.5 & glu3 < 126)
cat("landmark non-DM:", nrow(lm_ndm), "n; events", sum(lm_ndm$ev15), "\n")
S <- "Surv(t15, ev15)"
adj15_ndm <- paste(setdiff(cov15, "dm15"), collapse = " + ")
l1 <- coxph(as.formula(paste(S, "~ hgi15_hi_p90")), data = lm_ndm)
l2 <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 + age15 + sex + urban + edu")),
            data = lm_ndm)
l3 <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 +", adj15_ndm)), data = lm_ndm)
lk <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 + hba1c3 +", adj15_ndm)),
            data = lm_ndm)
out$landmark_ndm <- rbind(
  cbind(set = "CHARLS 2015 landmark non-DM", model = "M1", grab(l1, "hgi15_hi_p90")),
  cbind(set = "CHARLS 2015 landmark non-DM", model = "M2", grab(l2, "hgi15_hi_p90")),
  cbind(set = "CHARLS 2015 landmark non-DM", model = "M3", grab(l3, "hgi15_hi_p90")),
  cbind(set = "CHARLS 2015 landmark non-DM", model = "key", grab(lk, "hgi15_hi_p90")),
  cbind(set = "CHARLS 2015 landmark non-DM", model = "key: HbA1c", grab(lk, "hba1c3")))
print(out$landmark_ndm)
lm3 <- grab(l3, "hgi15_hi_p90")

## ================================================================================
## 4. Non-DM pooling (M3 estimates)
## ================================================================================
c11_ndm <- grab(a3, "hgi_hi_p90")                       # 1.475 1.147 1.896
elsa_ndm <- data.frame(HR = 1.511, lo = 0.918, hi = 2.485)  # frozen S10b
p_2c <- pool(c(c11_ndm$HR, lm3$HR), c(c11_ndm$lo, lm3$lo), c(c11_ndm$hi, lm3$hi))
p_ce <- pool(c(c11_ndm$HR, elsa_ndm$HR), c(c11_ndm$lo, elsa_ndm$lo),
             c(c11_ndm$hi, elsa_ndm$hi))
p_3  <- pool(c(c11_ndm$HR, lm3$HR, elsa_ndm$HR),
             c(c11_ndm$lo, lm3$lo, elsa_ndm$lo),
             c(c11_ndm$hi, lm3$hi, elsa_ndm$hi))
cat("non-DM pooled: two-CHARLS\n"); print(p_2c)
cat("non-DM pooled: CHARLS2011 + ELSA\n"); print(p_ce)
cat("non-DM pooled: three-set\n"); print(p_3)
out$pool2c <- cbind(set = "two-CHARLS non-DM", p_2c)
out$poolce <- cbind(set = "CHARLS2011+ELSA non-DM", p_ce)
out$pool3  <- cbind(set = "three-set non-DM", p_3)

## ================================================================================
## 5. Normoglycemia sensitivity (CHARLS 2011; exclude prediabetes; within-subset P90)
## ================================================================================
ng <- dmain %>% filter(diabetes == 0 & glu < 100 & hba1c < 5.7) %>%
  mutate(h_ng = as.integer(hgi >= quantile(hgi, 0.9, na.rm = TRUE)))
cat("normoglycemia set:", nrow(ng), "n; events", sum(ng$incident_stroke6),
    "; within-subset P90 =", round(quantile(ng$hgi, 0.9, na.rm = TRUE), 4), "\n")
g3 <- coxph(as.formula(paste(FS, "~ h_ng +", adj3_ndm)), data = ng)
gk <- coxph(as.formula(paste(FS, "~ h_ng + hba1c +", adj3_ndm)), data = ng)
out$normogly <- rbind(
  cbind(set = "CHARLS normoglycemia (within-subset P90)", model = "M3", grab(g3, "h_ng")),
  cbind(set = "CHARLS normoglycemia (within-subset P90)", model = "key", grab(gk, "h_ng")))
print(out$normogly)

write.csv(out$recon, file.path(OUT, "S12_nonDM_reconciliation.csv"), row.names = FALSE)
write.csv(out$landmark_ndm, file.path(OUT, "S12_landmark_nonDM.csv"), row.names = FALSE)
write.csv(out$normogly, file.path(OUT, "S12_normoglycemia.csv"), row.names = FALSE)
write.csv(rbind(out$pool2c, out$poolce, out$pool3),
          file.path(OUT, "S12_nonDM_pooled.csv"), row.names = FALSE)
write.csv(out$collinear, file.path(OUT, "S12_refit_collinearity.csv"), row.names = FALSE)
cat("saved.\n")
