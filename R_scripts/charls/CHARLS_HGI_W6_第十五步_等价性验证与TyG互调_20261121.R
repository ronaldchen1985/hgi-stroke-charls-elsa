####################################################################################
# S15 (2026-11-21): pre-freeze decisive checks (user-approved, from strategy doc)
#  A. Algebraic-equivalence package (CHARLS 2011 main set):
#     A1. continuous HGI + HbA1c model vs FPG + HbA1c model:
#         log-likelihood identity (linear reparametrization) and
#         coef identity  beta_HGI = -beta_FPG / b  (b = slope of HbA1c ~ FPG)
#     A2. binary-exposure AIC ladder: base / +HbA1c / +highHGI / +FPG /
#         +HbA1c+HGI(key) / +HbA1c+FPG / +HGI+FPG(S9)
#  B. Non-diabetic HGI x TyG mutual adjustment (stop-loss item 2):
#     CHARLS non-DM (n=7,283) and ELSA W2 non-DM (n=2,474):
#     M3 + high HGI alone / + TyG(per SD) alone / both (mutual);
#     continuous HGI per-SD x TyG per-SD mutual
# Output: results/S15_equivalence.csv, S15_AIC_ladder.csv, S15_TyG_mutual.csv
####################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(survival); library(dplyr) })
BASE <- "D:/CHARLS/HGI_stroke/freeze_20261120"
OUT  <- file.path(BASE, "results")
grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(term = term, HR = round(ci[,1],3), lo = round(ci[,3],3),
             hi = round(ci[,4],3), P = signif(p,3))
}

## ---- frozen sets ---------------------------------------------------------------
d <- read.csv(file.path(BASE, "data/charls/analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
dmain <- d %>% mutate(hgi_hi_p90 = as.integer(hgi >= CUT_P90)) %>%
  filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c",
                                        "incident_stroke6", "fu_v2")))))
stopifnot(nrow(dmain) == 8684, sum(dmain$incident_stroke6) == 852)
adj3 <- paste(cov_m3, collapse = " + ")
adj3_ndm <- paste(setdiff(cov_m3, "dm"), collapse = " + ")
FS <- "Surv(fu_v2, incident_stroke6)"
ndm <- dmain %>% filter(diabetes == 0)
stopifnot(nrow(ndm) == 7283, sum(ndm$incident_stroke6) == 681)

e <- read.csv(file.path(BASE, "data/elsa/analytic_elsa_hgi_baseline_w2.csv"))
cov_e <- c("age","sex","urban","edu","htn","dyslip","heart","sbp","crp","tg","hdl")
adj_e <- paste(cov_e, collapse = " + ")
SE <- "Surv(fu_mid, incident_stroke)"
en <- e %>% filter(dm == 0 & hba1c < 6.5 & fglu < 7.0)
stopifnot(nrow(en) == 2474, sum(en$incident_stroke) == 164)

## ================================================================================
## A1. Equivalence: continuous HGI + HbA1c  ==  FPG + HbA1c (reparametrization)
## ================================================================================
cat("\n==== A1. equivalence check (continuous) ====\n")
eq <- lm(hba1c ~ glu, data = dmain)
b_slope <- coef(eq)["glu"]          # within-set slope (should be ~0.0151)
m_hgi <- coxph(as.formula(paste(FS, "~ hgi + hba1c +", adj3)), data = dmain)
m_fpg <- coxph(as.formula(paste(FS, "~ glu + hba1c +", adj3)), data = dmain)
ll_hgi <- as.numeric(logLik(m_hgi)); ll_fpg <- as.numeric(logLik(m_fpg))
bh <- coef(m_hgi)["hgi"]; bf <- coef(m_fpg)["glu"]
cat("slope b =", round(b_slope, 6), "\n")
cat("logLik HGI model =", round(ll_hgi, 4), "; FPG model =", round(ll_fpg, 4),
    "; diff =", round(ll_hgi - ll_fpg, 6), "\n")
cat("coef HGI =", round(bh, 5), "; -coef FPG / b =", round(-bf / b_slope, 5), "\n")
eq_out <- data.frame(
  item = c("slope b (HbA1c ~ FPG, main set)",
           "logLik: continuous HGI + HbA1c + covs",
           "logLik: FPG + HbA1c + covs",
           "logLik difference",
           "coef HGI (continuous, HGI model)",
           "-coef FPG / b (FPG model)",
           "identity residual"),
  value = c(round(b_slope, 6), round(ll_hgi, 4), round(ll_fpg, 4),
            round(ll_hgi - ll_fpg, 6), round(bh, 5),
            round(-bf / b_slope, 5), round(bh + bf / b_slope, 8)))
print(eq_out)
write.csv(eq_out, file.path(OUT, "S15_equivalence.csv"), row.names = FALSE)

## ================================================================================
## A2. AIC ladder (binary high-HGI exposure, CHARLS 2011 main set)
## ================================================================================
cat("\n==== A2. AIC ladder ====\n")
fits <- list(
  base        = coxph(as.formula(paste(FS, "~", adj3)), data = dmain),
  plus_HbA1c  = coxph(as.formula(paste(FS, "~ hba1c +", adj3)), data = dmain),
  plus_highHGI= coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj3)), data = dmain),
  plus_FPG    = coxph(as.formula(paste(FS, "~ glu +", adj3)), data = dmain),
  key_HbA1c_HGI = coxph(as.formula(paste(FS, "~ hgi_hi_p90 + hba1c +", adj3)),
                        data = dmain),
  HbA1c_FPG   = coxph(as.formula(paste(FS, "~ hba1c + glu +", adj3)), data = dmain),
  S9_HGI_FPG  = coxph(as.formula(paste(FS, "~ hgi_hi_p90 + glu +", adj3)),
                      data = dmain))
aic_tab <- data.frame(
  model = names(fits),
  AIC = sapply(fits, AIC),
  logLik = sapply(fits, function(f) as.numeric(logLik(f))))
aic_tab$dAIC_vs_base <- round(aic_tab$AIC - aic_tab$AIC[1], 2)
print(aic_tab)
write.csv(aic_tab, file.path(OUT, "S15_AIC_ladder.csv"), row.names = FALSE)
term_rows <- rbind(
  cbind(model = "plus_HbA1c", grab(fits$plus_HbA1c, "hba1c")),
  cbind(model = "plus_highHGI", grab(fits$plus_highHGI, "hgi_hi_p90")),
  cbind(model = "plus_FPG", grab(fits$plus_FPG, "glu")),
  cbind(model = "key_HbA1c_HGI", grab(fits$key_HbA1c_HGI, "hgi_hi_p90")),
  cbind(model = "key_HbA1c_HGI", grab(fits$key_HbA1c_HGI, "hba1c")),
  cbind(model = "HbA1c_FPG", grab(fits$HbA1c_FPG, "hba1c")),
  cbind(model = "HbA1c_FPG", grab(fits$HbA1c_FPG, "glu")),
  cbind(model = "S9_HGI_FPG", grab(fits$S9_HGI_FPG, "hgi_hi_p90")),
  cbind(model = "S9_HGI_FPG", grab(fits$S9_HGI_FPG, "glu")))
print(term_rows)
write.csv(term_rows, file.path(OUT, "S15_AIC_ladder_terms.csv"), row.names = FALSE)

## ================================================================================
## B. Non-diabetic HGI x TyG mutual adjustment
## ================================================================================
cat("\n==== B. HGI x TyG mutual adjustment, non-diabetic ====\n")
run_block <- function(dat, S, adj, tag) {
  dat$tyg_z <- as.numeric(scale(log(dat$tg * dat$glu / 2)))
  dat$hgi_sd <- as.numeric(scale(dat$hgi))
  m_h  <- coxph(as.formula(paste(S, "~ hgi_hi_p90 +", adj)), data = dat)
  m_t  <- coxph(as.formula(paste(S, "~ tyg_z +", adj)), data = dat)
  m_ht <- coxph(as.formula(paste(S, "~ hgi_hi_p90 + tyg_z +", adj)), data = dat)
  m_c  <- coxph(as.formula(paste(S, "~ hgi_sd +", adj)), data = dat)
  m_ct <- coxph(as.formula(paste(S, "~ hgi_sd + tyg_z +", adj)), data = dat)
  out <- rbind(
    cbind(set = tag, model = "high HGI alone", grab(m_h, "hgi_hi_p90")),
    cbind(set = tag, model = "TyG per-SD alone", grab(m_t, "tyg_z")),
    cbind(set = tag, model = "mutual: high HGI | TyG", grab(m_ht, "hgi_hi_p90")),
    cbind(set = tag, model = "mutual: TyG | high HGI", grab(m_ht, "tyg_z")),
    cbind(set = tag, model = "continuous HGI per-SD alone", grab(m_c, "hgi_sd")),
    cbind(set = tag, model = "mutual: HGI per-SD | TyG", grab(m_ct, "hgi_sd")),
    cbind(set = tag, model = "mutual: TyG | HGI per-SD", grab(m_ct, "tyg_z")))
  out$n <- nrow(dat); out$events <- sum(dat[[ifelse(S == FS, "incident_stroke6",
                                                    "incident_stroke")]])
  out
}
b1 <- run_block(ndm, FS, adj3_ndm, "CHARLS 2011 non-DM")
## ELSA TyG: file tg is mg/dL; convert fglu (mmol/L) to mg/dL for the TyG formula
b2 <- run_block(en %>% mutate(glu = fglu * 18.016), SE, adj_e, "ELSA W2 non-DM")
tyg_out <- rbind(b1, b2)
print(tyg_out)
write.csv(tyg_out, file.path(OUT, "S15_TyG_mutual.csv"), row.names = FALSE)

cat("\nS15 done. saved: S15_equivalence.csv, S15_AIC_ladder.csv, S15_TyG_mutual.csv\n")
