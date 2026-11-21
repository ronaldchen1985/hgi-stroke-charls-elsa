####################################################################################
# S13 (2026-11-20): treatment-stratified diabetes analysis + eGFR + early events
#  A. Diabetes treatment strata (CHARLS 2011 main cc set):
#     non-DM / untreated-or-undiagnosed DM / treated DM (da014s1-s3 endorsed)
#     -> tests the mechanism claim "treatment distorts HGI"
#  B. eGFR sensitivity (CKD-EPI 2021 race-free, crea mg/dL):
#     B1 exclude eGFR < 60; B2 M3 + continuous eGFR
#  C. Reverse causation: exclude incident strokes within first 2 years (fu_v2 <= 2)
# Output: results/S13_treatment_eGFR_early.csv + console
####################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(haven); library(survival); library(dplyr) })
BASE <- "D:/CHARLS/HGI_stroke/freeze_20261120"
OUT  <- file.path(BASE, "results")
grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(term = term, HR = round(ci[,1],3), lo = round(ci[,3],3),
             hi = round(ci[,4],3), P = signif(p,3))
}

d <- read.csv(file.path(BASE, "data/charls/analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
dmain <- d %>% mutate(hgi_hi_p90 = as.integer(hgi >= CUT_P90)) %>%
  filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c",
                                        "incident_stroke6", "fu_v2")))))
stopifnot(nrow(dmain) == 8684, sum(dmain$incident_stroke6) == 852)
FS <- "Surv(fu_v2, incident_stroke6)"
adj3_ndm <- paste(setdiff(cov_m3, "dm"), collapse = " + ")
adj3 <- paste(cov_m3, collapse = " + ")
out <- list()

## ================================================================================
## A. Treatment strata (diabetes treatment from W1 da014s1-s4)
## ================================================================================
hs1 <- read_dta("D:/CHARLS/2011/Health_Status_and_Functioning.dta")
trt <- data.frame(ID = as.character(hs1$ID),
                  dm_treated = as.integer(
                    (!is.na(as.numeric(hs1$da014s1)) & as.numeric(hs1$da014s1) == 1) |
                    (!is.na(as.numeric(hs1$da014s2)) & as.numeric(hs1$da014s2) == 1) |
                    (!is.na(as.numeric(hs1$da014s3)) & as.numeric(hs1$da014s3) == 1)))
dd <- dmain %>% left_join(trt, by = "ID") %>%
  mutate(dm_treated = ifelse(is.na(dm_treated), 0L, dm_treated),
         stratum = case_when(diabetes == 0 ~ "non-diabetic",
                             diabetes == 1 & dm_treated == 0 ~ "untreated/undiagnosed DM",
                             TRUE ~ "treated DM"))
print(table(dd$stratum))
print(tapply(dd$incident_stroke6, dd$stratum, sum))
for (s in c("non-diabetic", "untreated/undiagnosed DM", "treated DM")) {
  sub <- dd %>% filter(stratum == s)
  if (nrow(sub) < 300) { cat(s, ": too small\n"); next }
  m <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj3_ndm)), data = sub)
  g <- grab(m, "hgi_hi_p90")
  cat(s, ": n =", nrow(sub), ", events =", sum(sub$incident_stroke6),
      ", M3 HR =", g$HR, "(", g$lo, "-", g$hi, ") P =", g$P, "\n")
  out[[paste0("A_", s)]] <- cbind(analysis = paste("treatment stratum:", s),
                                  n = nrow(sub), events = sum(sub$incident_stroke6), g)
}
## interaction: high HGI x treated-DM vs the rest (descriptive)
dd <- dd %>% mutate(dm_any = as.integer(diabetes == 1))
mi <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 * dm_treated +", adj3_ndm)),
            data = dd %>% filter(diabetes == 1))
cat("within-diabetics interaction HGI x treated P =",
    signif(summary(mi)$coefficients["hgi_hi_p90:dm_treated", "Pr(>|z|)"], 3), "\n")
out$A_interact <- cbind(analysis = "within-DM interaction: high HGI x treated",
                        n = sum(dd$diabetes == 1),
                        events = sum(dd$incident_stroke6[dd$diabetes == 1]),
                        grab(mi, "hgi_hi_p90:dm_treated"))

## ================================================================================
## B. eGFR sensitivity (CKD-EPI 2021 race-free; crea mg/dL; sex 1=male 2=female)
## ================================================================================
egfr <- function(scr, age, female) {
  k  <- ifelse(female == 1, 0.7, 0.9)
  a  <- ifelse(female == 1, -0.241, -0.302)
  e  <- 142 * pmin(scr / k, 1)^a * pmax(scr / k, 1)^(-1.200) * 0.9938^age
  ifelse(female == 1, e * 1.012, e)
}
dd <- dd %>% mutate(female = as.integer(sex == 2),
                    egfr = egfr(crea, age, female))
cat("eGFR available n =", sum(!is.na(dd$egfr)),
    "; eGFR<60 n =", sum(dd$egfr < 60, na.rm = TRUE), "\n")
## B1: exclude CKD (eGFR<60), require measured creatinine
b1 <- dd %>% filter(!is.na(egfr) & egfr >= 60)
mb1 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj3)), data = b1)
g <- grab(mb1, "hgi_hi_p90")
cat("B1 exclude eGFR<60: n =", nrow(b1), ", events =", sum(b1$incident_stroke6),
    ", M3 HR =", g$HR, "(", g$lo, "-", g$hi, ") P =", g$P, "\n")
out$B1 <- cbind(analysis = "exclude CKD (eGFR<60, creatinine-measured)",
                n = nrow(b1), events = sum(b1$incident_stroke6), g)
## B2: M3 + continuous eGFR
b2 <- dd %>% filter(!is.na(egfr))
mb2 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 + egfr +", adj3)), data = b2)
g <- grab(mb2, "hgi_hi_p90")
cat("B2 M3 + eGFR: n =", nrow(b2), ", events =", sum(b2$incident_stroke6),
    ", HR =", g$HR, "(", g$lo, "-", g$hi, ") P =", g$P, "\n")
out$B2 <- cbind(analysis = "M3 + continuous eGFR", n = nrow(b2),
                events = sum(b2$incident_stroke6), g)

## ================================================================================
## C. Exclude incident strokes within first 2 years (reverse causation)
## ================================================================================
c1 <- dd %>% filter(!(incident_stroke6 == 1 & fu_v2 <= 2))
mc1 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj3)), data = c1)
g <- grab(mc1, "hgi_hi_p90")
cat("C exclude 2y events: n =", nrow(c1), ", events =", sum(c1$incident_stroke6),
    ", M3 HR =", g$HR, "(", g$lo, "-", g$hi, ") P =", g$P, "\n")
out$C <- cbind(analysis = "exclude strokes within first 2 years",
               n = nrow(c1), events = sum(c1$incident_stroke6), g)

res <- do.call(rbind, out)
print(res)
write.csv(res, file.path(OUT, "S13_treatment_eGFR_early.csv"), row.names = FALSE)
cat("saved.\n")
