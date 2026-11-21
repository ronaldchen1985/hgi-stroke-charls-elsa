####################################################################################
# Revision supplement (mock-review driven), 2026-11-20
#  S9a: ELSA W2 calibration-equation R-squared (review item R1-m5 / AUD-06)
#  S9b: CHARLS joint FPG + HGI models (review item R1-M1)
#       HbA1c = predicted HbA1c + HGI by construction; the complementary
#       specification enters FPG and HGI jointly so all components of HbA1c
#       appear in one model.
# Inputs: frozen analytic CSVs in freeze_20261120/data
####################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(survival); library(dplyr) })
BASE <- "D:/CHARLS/HGI_stroke/freeze_20261120"
OUT  <- file.path(BASE, "results")

## ---- S9a: ELSA equation R-squared ------------------------------------------------
e <- read.csv(file.path(BASE, "data/elsa/analytic_elsa_hgi_full_w2.csv"))  # n = 3,062 eligible
fit_e <- lm(hba1c ~ fglu, data = e)
cat("ELSA W2 equation (n =", nrow(e), "): intercept =", round(coef(fit_e)[1], 3),
    " slope =", round(coef(fit_e)[2], 3),
    " R2 =", round(summary(fit_e)$r.squared, 3), "\n")
cat("  (frozen manuscript equation: 3.827 + 0.330 x FPG)\n\n")

## ---- S9b: CHARLS joint FPG + HGI models ------------------------------------------
d <- read.csv(file.path(BASE, "data/charls/analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
CUT <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
d$hgi_hi_p90 <- as.integer(d$hgi >= CUT)
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
adj3 <- paste(cov_m3, collapse = " + ")
dmain <- d %>% filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c",
                                                     "incident_stroke6", "fu_v2")))))
cat("CHARLS main complete-case set:", nrow(dmain), "n; events", sum(dmain$incident_stroke6), "\n")
cat("cor(HGI, FPG) =", round(cor(dmain$hgi, dmain$glu), 4),
    " (orthogonality check; OLS residual vs predictor)\n\n")

grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  c(term = term, HR = round(ci[,1],3), lo = round(ci[,3],3), hi = round(ci[,4],3),
    P = signif(p,3))
}
FS <- "Surv(fu_v2, incident_stroke6)"

# (i) dichotomous high HGI + FPG (replaces HbA1c key test with the FPG component)
m_j1 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 + glu +", adj3)), data = dmain)
# (ii) continuous HGI + FPG, both scaled per SD
dmain <- dmain %>% mutate(hgi_z = as.numeric(scale(hgi)), glu_z = as.numeric(scale(glu)))
m_j2 <- coxph(as.formula(paste(FS, "~ hgi_z + glu_z +", adj3)), data = dmain)

res <- rbind(
  data.frame(model = "high HGI (>=P90) + FPG", t(grab(m_j1, "hgi_hi_p90"))),
  data.frame(model = "high HGI (>=P90) + FPG", t(grab(m_j1, "glu"))),
  data.frame(model = "continuous joint (per SD)", t(grab(m_j2, "hgi_z"))),
  data.frame(model = "continuous joint (per SD)", t(grab(m_j2, "glu_z")))
)
print(res, row.names = FALSE)

write.csv(res, file.path(OUT, "S9_revision_supplement_CHARLS_joint_FPG_HGI.csv"),
          row.names = FALSE)
write.csv(data.frame(n = nrow(e), intercept = coef(fit_e)[1], slope = coef(fit_e)[2],
                     r2 = summary(fit_e)$r.squared),
          file.path(OUT, "S9_revision_supplement_ELSA_equation_R2.csv"), row.names = FALSE)
cat("\nSaved to", OUT, "\n")
