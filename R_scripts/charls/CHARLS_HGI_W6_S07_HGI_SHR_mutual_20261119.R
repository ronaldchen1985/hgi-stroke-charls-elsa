####################################################################################
# CHARLS HGI x stroke -- Step 7: HGI x SHR mutual-adjustment & interaction (W6)
# Aim: respond to the "is HGI independent of acute glucose ratio (SHR)" critique
#   (Tong 2026-style). SHR = (FPG mmol/L) / (1.59*HbA1c - 2.59), same as M4.
# Sample: identical to main analysis (dmain, 12-covariate complete case, fu_v2).
# Output: HGI_SHR_mutual_W6.txt / .rds / 表S7_HGI_SHR互调_W6版.csv
# Seed 20260903. Date 2026-11-19.
####################################################################################

rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(survival); library(dplyr) })
OUT_DIR <- "D:/CHARLS/HGI_stroke"

d <- read.csv(file.path(OUT_DIR, "analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
d <- d %>% mutate(
  hgi_hi_p90 = as.integer(hgi >= CUT_P90),
  shr = (glu / 18) / (1.59 * hba1c - 2.59))

cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
adj3 <- "age + sex + urban + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"
dmain <- d %>%
  filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c", "shr",
                                        "incident_stroke6", "fu_v2"))))) %>%
  mutate(hgi_z = as.numeric(scale(hgi)),
         shr_z = as.numeric(scale(shr)),
         shr_hi_med = as.integer(shr >= median(shr)))
message("Step7 sample: ", nrow(dmain), " n; events ", sum(dmain$incident_stroke6),
        "; deaths ", sum(dmain$death6 == 1, na.rm = TRUE))

grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(term = term, HR = round(ci[,1],3), lo = round(ci[,3],3),
             hi = round(ci[,4],3), P = signif(p,3))
}
FS <- "Surv(fu_v2, incident_stroke6)"
res <- list()

## ---- 0. correlation & collinearity --------------------------------------------------
cat("==== 0. correlation / collinearity ====\n")
res$cor <- data.frame(
  pearson_hgi_shr  = round(cor(dmain$hgi, dmain$shr), 3),
  spearman_hgi_shr = round(cor(dmain$hgi, dmain$shr, method = "spearman"), 3))
print(res$cor)
# VIF for the two exposures in the fully adjusted mutual model (linear analog)
fm_hgi <- lm(as.formula(paste("hgi_z ~ shr_z +", adj3)), data = dmain)
fm_shr <- lm(as.formula(paste("shr_z ~ hgi_z +", adj3)), data = dmain)
res$vif <- data.frame(exposure = c("hgi_z","shr_z"),
                      VIF = round(c(1/(1-summary(fm_hgi)$r.squared),
                                    1/(1-summary(fm_shr)$r.squared)), 2))
print(res$vif)

## ---- 1. single-exposure anchors ------------------------------------------------------
cat("\n==== 1. single-exposure models ====\n")
a1 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj3)), data = dmain)
res$A1_hgiP90_alone <- grab(a1, "hgi_hi_p90"); print(res$A1_hgiP90_alone)
a2 <- coxph(as.formula(paste(FS, "~ shr_z +", adj3)), data = dmain)
res$A2_shrSD_alone <- grab(a2, "shr_z"); print(res$A2_shrSD_alone)
a3 <- coxph(as.formula(paste(FS, "~ hgi_z +", adj3)), data = dmain)
res$A3_hgiSD_alone <- grab(a3, "hgi_z"); print(res$A3_hgiSD_alone)

## ---- 2. mutual adjustment ------------------------------------------------------------
cat("\n==== 2. mutual adjustment ====\n")
b1 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 + shr_z +", adj3)), data = dmain)
res$B1_mutual_P90 <- bind_rows(grab(b1,"hgi_hi_p90"), grab(b1,"shr_z")); print(res$B1_mutual_P90)
b2 <- coxph(as.formula(paste(FS, "~ hgi_z + shr_z +", adj3)), data = dmain)
res$B2_mutual_SD  <- bind_rows(grab(b2,"hgi_z"), grab(b2,"shr_z")); print(res$B2_mutual_SD)

## ---- 3. interaction tests -------------------------------------------------------------
cat("\n==== 3. interaction ====\n")
c_add <- coxph(as.formula(paste(FS, "~ hgi_z + shr_z +", adj3)), data = dmain)
c_int <- coxph(as.formula(paste(FS, "~ hgi_z * shr_z +", adj3)), data = dmain)
lrt_c <- anova(c_add, c_int)
res$C1_interaction_SD <- data.frame(P_interaction = signif(lrt_c[["Pr(>|Chi|)"]][2], 3))
print(res$C1_interaction_SD)

d_add <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 + shr_z +", adj3)), data = dmain)
d_int <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 * shr_z +", adj3)), data = dmain)
lrt_d <- anova(d_add, d_int)
res$C2_interaction_P90 <- data.frame(P_interaction = signif(lrt_d[["Pr(>|Chi|)"]][2], 3))
print(res$C2_interaction_P90)

## ---- 4. stratified estimates ----------------------------------------------------------
cat("\n==== 4. stratified ====\n")
strat <- lapply(c(0,1), function(L) {
  dd <- filter(dmain, shr_hi_med == L)
  m <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj3)), data = dd)
  g <- grab(m, "hgi_hi_p90"); g$stratum <- ifelse(L==0, "SHR below median", "SHR at/above median")
  g$n <- nrow(dd); g$events <- sum(dd$incident_stroke6); g
}) %>% bind_rows()
res$D1_hgi_by_SHRstrata <- strat; print(strat)
strat2 <- lapply(c(0,1), function(L) {
  dd <- filter(dmain, hgi_hi_p90 == L)
  m <- coxph(as.formula(paste(FS, "~ shr_z +", adj3)), data = dd)
  g <- grab(m, "shr_z"); g$stratum <- ifelse(L==0, "HGI below P90", "HGI at/above P90")
  g$n <- nrow(dd); g$events <- sum(dd$incident_stroke6); g
}) %>% bind_rows()
res$D2_shr_by_HGIstrata <- strat2; print(strat2)

## ---- export ---------------------------------------------------------------------------
sink(file.path(OUT_DIR, "HGI_SHR_mutual_W6.txt"))
cat("CHARLS HGI x SHR mutual adjustment & interaction (", format(Sys.time(), "%F %T"), ")\n",
    "Sample = W6 main set: ", nrow(dmain), " n; events ", sum(dmain$incident_stroke6),
    "; deaths ", sum(dmain$death6 == 1, na.rm = TRUE),
    "\nSHR = (FPG mmol/L)/(1.59*HbA1c - 2.59); per-SD = z-scored; adj3 = 12 covariates; time = fu_v2\n", sep = "")
for (nm in names(res)) { cat("\n##### ", nm, " #####\n"); print(res[[nm]]) }
sink()
saveRDS(res, file.path(OUT_DIR, "HGI_SHR_mutual_W6.rds"))

tab <- bind_rows(
  res$A1_hgiP90_alone %>% mutate(model = "A1 HGI>=P90 alone"),
  res$A2_shrSD_alone  %>% mutate(model = "A2 SHR per SD alone"),
  res$A3_hgiSD_alone  %>% mutate(model = "A3 HGI per SD alone"),
  res$B1_mutual_P90   %>% mutate(model = "B1 mutual (P90 + SHR-SD)"),
  res$B2_mutual_SD    %>% mutate(model = "B2 mutual (HGI-SD + SHR-SD)"),
  res$D1_hgi_by_SHRstrata %>% select(term, HR, lo, hi, P, stratum) %>%
    mutate(model = paste0("D1 HGI-P90 within ", stratum)) %>% select(-stratum),
  res$D2_shr_by_HGIstrata %>% select(term, HR, lo, hi, P, stratum) %>%
    mutate(model = paste0("D2 SHR-SD within ", stratum)) %>% select(-stratum))
write.csv(tab, file.path(OUT_DIR, "表S7_HGI_SHR互调_W6版.csv"),
          row.names = FALSE, fileEncoding = "GB18030")
message("\nDONE: HGI_SHR_mutual_W6.txt / .rds / 表S7_HGI_SHR互调_W6版.csv")
