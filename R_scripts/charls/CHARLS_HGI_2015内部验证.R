####################################################################################
# CHARLS HGI x stroke -- internal validation with the 2015 (W3) blood draw
# Landmark analysis: baseline = 2015; exposure = HGI derived within the 2015
# (KingMed) blood batch; outcome = incident stroke 2015-2023 (W4/W5/W6).
# Mirrors the ELSA external-validation logic, inside CHARLS.
# Output: internal_validation_2015.txt / .rds  (in D:/CHARLS/HGI_stroke)
####################################################################################
rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(survival); library(dplyr); library(splines);
                                 library(haven) })
OUT_DIR <- "D:/CHARLS/HGI_stroke"
adj3 <- "age15 + sex + urban + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"
adj2 <- "age15 + sex + urban + edu"
cov_m3 <- c("age15","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")

grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(term = term, HR = round(ci[,1],3), lo = round(ci[,3],3),
             hi = round(ci[,4],3), P = signif(p,3))
}

## ---- 1. 2015 raw blood: refit the within-batch equation ------------------------------
w3_blood <- read_dta("D:/CHARLS/2015/Blood.dta")
v_id3 <- "ID"; v_glu3 <- "bl_glu"; v_a1c3 <- "bl_hbalc"
blood3 <- w3_blood %>%
  select(ID = all_of(v_id3), glu3 = all_of(v_glu3), hba1c3 = all_of(v_a1c3)) %>%
  mutate(across(-ID, as.numeric), ID = as.character(ID))
to_w1_id <- function(id) { id <- as.character(id)
  hit <- nchar(id) == 12 & substr(id, 10, 10) == "0"
  id[hit] <- paste0(substr(id[hit], 1, 9), substr(id[hit], 11, 12)); id }
blood3$ID <- to_w1_id(blood3$ID)
fit_2015 <- lm(hba1c3 ~ glu3, data = blood3)
cat("2015 equation refit: HbA1c =", round(coef(fit_2015)[1],4), "+",
    round(coef(fit_2015)[2],4), "x FPG; R2 =", round(summary(fit_2015)$r.squared,3), "\n")

## ---- 2. landmark analytic set ----------------------------------------------------------
dd <- read.csv(file.path(OUT_DIR, "analytic_delta_hgi.csv"), colClasses = c(ID = "character"))
d  <- read.csv(file.path(OUT_DIR, "analytic_hgi_baseline_w6.csv"), colClasses = c(ID = "character"))
adjv <- strsplit(gsub(" ", "", gsub("age15","age",adj3)), "\\+")[[1]]
dd <- dd %>%
  select(-any_of(intersect(adjv, setdiff(names(dd), "ID")))) %>%
  left_join(d %>% select(ID, stroke_w6, last_contact, w6_report_yr, w6_dxyear,
                         all_of(adjv)), by = "ID") %>%
  left_join(blood3 %>% select(ID, glu3, hba1c3), by = "ID") %>%
  filter(is.na(first_stroke_wave) | first_stroke_wave >= 4) %>%   # stroke-free at 2015
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

CUT15 <- as.numeric(quantile(dd$hgi_2015, 0.9, na.rm = TRUE))
dd <- dd %>% mutate(hgi15_hi_p90 = as.integer(hgi_2015 >= CUT15),
                    hgi15_sd = as.numeric(scale(hgi_2015)))
cat("landmark set (2015 baseline):", nrow(dd), "n; events 2015-2023:", sum(dd$ev15),
    "; high-HGI(2015 P90 =", round(CUT15,4), "):", sum(dd$hgi15_hi_p90, na.rm=TRUE), "\n")
cat("crude event rate: high", round(mean(dd$ev15[dd$hgi15_hi_p90==1], na.rm=TRUE)*100,1),
    "% vs low", round(mean(dd$ev15[dd$hgi15_hi_p90==0], na.rm=TRUE)*100,1), "%\n")

res <- list()
sink(file.path(OUT_DIR, "internal_validation_2015.txt"), split = TRUE)
cat("CHARLS internal validation: 2015 (W3) landmark analysis,", format(Sys.time(), "%F %T"), "\n")
cat("n =", nrow(dd), "; events =", sum(dd$ev15), "; high-HGI =",
    sum(dd$hgi15_hi_p90, na.rm = TRUE), "; P90(2015) =", round(CUT15, 4), "\n\n")

## ---- 3. main models ----------------------------------------------------------------------
dmain <- dd %>% filter(complete.cases(across(all_of(c(cov_m3, "hgi_2015", "ev15", "t15")))))
cat("complete-case set:", nrow(dmain), "n; events", sum(dmain$ev15), "\n\n")
S <- paste("Surv(t15, ev15)")
m1 <- coxph(as.formula(paste(S, "~ hgi15_hi_p90")), data = dmain)
m2 <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 +", adj2)), data = dmain)
m3 <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 +", adj3)), data = dmain)
cat("== M1/M2/M3 (high HGI-2015 vs low) ==\n")
res$M1 <- grab(m1, "hgi15_hi_p90"); res$M2 <- grab(m2, "hgi15_hi_p90"); res$M3 <- grab(m3, "hgi15_hi_p90")
print(do.call(rbind, res)); cat("PH test (M3):\n"); print(cox.zph(m3)$table["hgi15_hi_p90", ])

## ---- 4. key test (+ HbA1c-2015) -----------------------------------------------------------
m3k <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 + hba1c3 +", adj3)), data = dmain)
cat("\n== key test (+ HbA1c-2015) ==\n")
print(summary(m3k)$conf.int[c("hgi15_hi_p90","hba1c3"), ])
print(summary(m3k)$coefficients[c("hgi15_hi_p90","hba1c3"), "Pr(>|z|)"])
res$key <- grab(m3k, "hgi15_hi_p90")

## ---- 5. per-SD + RCS -----------------------------------------------------------------------
m3s <- coxph(as.formula(paste(S, "~ hgi15_sd +", adj3)), data = dmain)
cat("\nper-SD HGI-2015:\n"); print(grab(m3s, "hgi15_sd"))
res$perSD <- grab(m3s, "hgi15_sd")
m_lin <- coxph(as.formula(paste(S, "~ hgi_2015 +", adj3)), data = dmain)
m_rcs <- coxph(as.formula(paste(S, "~ ns(hgi_2015, df = 4) +", adj3)), data = dmain)
lrt <- anova(m_lin, m_rcs)
cat("\nRCS nonlinearity LRT:\n"); print(lrt)

## ---- 6. mutual adjustment with 2011 HGI (trait vs state) ----------------------------------
m3b <- coxph(as.formula(paste(S, "~ hgi15_hi_p90 + hgi +", adj3)), data = dmain)
cat("\n== HGI-2015 adjusted for HGI-2011 ==\n")
print(summary(m3b)$conf.int[c("hgi15_hi_p90","hgi"), ])
res$adj2011 <- grab(m3b, "hgi15_hi_p90")

saveRDS(res, file.path(OUT_DIR, "internal_validation_2015.rds"))
cat("\nDONE\n")
sink()
