####################################################################################
# S11 (2026-11-20): values needed for the landmark-2015-covariate switch
#  on the new 2015-covariate complete-case set (6,581 / 706):
#   - per-SD HGI-2015, RCS LRT, mutual adjustment with HGI-2011
#   - Fine-Gray for the CHARLS non-diabetic set (Table 2 completeness)
#   - exact pooled estimates (two-CHARLS; three-set), fixed + DL, I2
#   - median follow-up of the new landmark set
# Output: results/S11_switch_values.csv + console
####################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(haven); library(survival); library(dplyr); library(splines) })
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

## ---- rebuild new landmark set (same as S10a) -------------------------------------
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
d  <- read.csv(file.path(BASE, "data/charls/analytic_hgi_baseline_w6.csv"),
               colClasses = c(ID = "character"))
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
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
adj15 <- paste(cov15, collapse = " + ")
S <- "Surv(t15, ev15)"
dmain15 <- dd %>% mutate(hgi15_hi_p90 = as.integer(hgi_2015 >= CUT15),
                         hgi15_sd = as.numeric(scale(hgi_2015))) %>%
  filter(complete.cases(across(all_of(c(cov15, "hgi_2015", "hba1c3", "hgi", "ev15", "t15")))))
cat("new landmark cc set:", nrow(dmain15), "; events:", sum(dmain15$ev15),
    "; median t15:", median(dmain15$t15), "yrs\n")

out <- list()
## per-SD + RCS
m_sd <- coxph(as.formula(paste(S, "~ hgi15_sd +", adj15)), data = dmain15)
out$persd <- cbind(model = "per-SD HGI-2015 (2015 covariates)", grab(m_sd, "hgi15_sd"))
m_lin <- coxph(as.formula(paste(S, "~ hgi_2015 +", adj15)), data = dmain15)
m_rcs <- coxph(as.formula(paste(S, "~ ns(hgi_2015, df = 4) +", adj15)), data = dmain15)
lrt <- anova(m_lin, m_rcs)
cat("RCS LRT P =", signif(lrt$`Pr(>|Chi|)`[2], 3), "\n")
out$rcs <- data.frame(model = "RCS nonlinearity LRT P (2015 covariates)",
                      term = "ns df=4", HR = NA, lo = NA, hi = NA,
                      P = signif(lrt$`Pr(>|Chi|)`[2], 3))
## mutual adjustment
m_mut <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 + hgi +", adj15)), data = dmain15)
out$mut1 <- cbind(model = "mutual adj: high HGI-2015", grab(m_mut, "hgi15_hi_p90"))
out$mut2 <- cbind(model = "mutual adj: continuous HGI-2011", grab(m_mut, "hgi"))
## non-DM Fine-Gray
dmain <- d %>% mutate(hgi_hi_p90 = as.integer(hgi >= quantile(hgi, 0.9, na.rm = TRUE))) %>%
  filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c", "incident_stroke6", "fu_v2")))))
ndm <- dmain %>% filter(diabetes == 0)
if (requireNamespace("cmprsk", quietly = TRUE)) {
  library(cmprsk)
  fstat <- ifelse(ndm$incident_stroke6 == 1, 1, ifelse(ndm$death6 == 1, 2, 0))
  ftime <- ndm$fu_v2
  x <- model.matrix(~ hgi_hi_p90 + age + sex + urban + edu + htn + dyslip +
                    heart + sbp + crp + tg + hdl, data = ndm)[, -1]
  fg <- crr(ftime, fstat, x, failcode = 1, cencode = 0)
  b <- summary(fg)$coef["hgi_hi_p90", ]
  cat("non-DM Fine-Gray SHR:", round(exp(b[1]),3), " (",
      round(exp(b[1] - 1.96*b[3]),3), "-", round(exp(b[1] + 1.96*b[3]),3),
      ") P =", signif(b[5],3), "\n")
  out$fg <- data.frame(model = "non-DM Fine-Gray SHR", term = "hgi_hi_p90",
                       HR = round(exp(b[1]),3), lo = round(exp(b[1]-1.96*b[3]),3),
                       hi = round(exp(b[1]+1.96*b[3]),3), P = signif(b[5],3))
}
## exact pooling
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
p2 <- pool(c(1.352, 1.152), c(1.108, 0.920), c(1.650, 1.442))
cat("two-CHARLS pooled (2011 + 2015-new):\n"); print(p2)
p3 <- pool(c(1.352, 1.152, 1.420), c(1.108, 0.920, 0.903), c(1.650, 1.442, 2.235))
cat("three-set pooled:\n"); print(p3)
p1e <- pool(c(1.352, 1.420), c(1.108, 0.903), c(1.650, 2.235))
cat("2011 + ELSA pooled (anchor, expect 1.363):\n"); print(p1e)

res <- do.call(rbind, out)
write.csv(res, file.path(OUT, "S11_switch_values.csv"), row.names = FALSE)
write.csv(p2, file.path(OUT, "S11_pooled_two_CHARLS.csv"), row.names = FALSE)
write.csv(p3, file.path(OUT, "S11_pooled_three_sets.csv"), row.names = FALSE)
cat("saved.\n")
