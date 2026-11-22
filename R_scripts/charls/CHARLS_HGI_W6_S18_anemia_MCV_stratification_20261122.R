####################################################################################
# R4 revision analysis (2026-11-22), mock-review driven (clinical-chemistry seat W1)
# Anemia etiology stratification by MCV in the CHARLS 2011 main set:
#   microcytic anemia (MCV < 80 fL, iron-deficiency-like; expected to RAISE HbA1c)
#   vs non-microcytic anemia (MCV >= 80 fL; hemolysis/blood-loss-like, LOWER HbA1c)
# Tests: (A) biomarker shift by subtype; (B) HGI-tail enrichment;
#        (C) crude event rates (comorbidity-marker alternative);
#        (D) Cox high-HGI and low-end tests under subtype-specific exclusions;
#        (E) five-group mid-reference model under each exclusion set.
# Input : freeze_20261120/data/charls/analytic_hgi_baseline_w6.csv
# Output: results/R4_anemia_mcv_patterns_20261122.csv
#         results/R4_anemia_mcv_cox_20261122.csv
#         results/R4_anemia_mcv_fivegroup_20261122.csv
#         results/R4_log_20261122.txt
####################################################################################

rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(survival); library(dplyr) })
BASE <- "D:/CHARLS/HGI_stroke/freeze_20261120"
OUT  <- file.path(BASE, "results")

sink(file.path(OUT, "R4_log_20261122.txt"), split = TRUE)
on.exit(sink(NULL), add = TRUE)

cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
adj3   <- paste(cov_m3, collapse = " + ")
FS     <- "Surv(fu_v2, incident_stroke6)"

grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(term = term, HR = round(ci[,1], 4), lo = round(ci[,3], 4),
             hi = round(ci[,4], 4), P = signif(p, 4))
}

d <- read.csv(file.path(BASE, "data/charls/analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
CUT_P10 <- as.numeric(quantile(d$hgi, 0.1, na.rm = TRUE))
dmain <- d %>%
  filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c",
                                        "incident_stroke6", "fu_v2"))))) %>%
  mutate(hgi_hi_p90 = as.integer(hgi >= CUT_P90),
         hgi_lo_p10 = as.integer(hgi <  CUT_P10),
         anem_sub = case_when(is.na(anemia) | is.na(mcv) ~ NA_character_,
                              anemia == 0                ~ "no_anemia",
                              anemia == 1 & mcv <  80    ~ "microcytic",
                              anemia == 1 & mcv >= 80    ~ "non_microcytic"),
         hgi_g5 = cut(hgi, breaks = quantile(hgi, c(0,.1,.25,.75,.9,1), na.rm = TRUE),
                      include.lowest = TRUE, labels = c("<P10","P10-25","P25-75","P75-90",">=P90")))
cat("main set:", nrow(dmain), "n; events", sum(dmain$incident_stroke6), "\n")
cat("P10 =", round(CUT_P10, 3), " P90 =", round(CUT_P90, 3), "\n")
print(table(dmain$anem_sub, useNA = "ifany"))
da <- dmain %>% filter(!is.na(anem_sub)) %>% mutate(anem_sub = factor(anem_sub))
cat("analysis set with MCV-classified anemia:", nrow(da), "n; events",
    sum(da$incident_stroke6), "\n\n")

## ---- A. biomarker pattern by anemia subtype ---------------------------------------
cat("==== A. biomarker pattern (unadjusted means) ====\n")
pat <- da %>% group_by(anem_sub) %>%
  summarise(n = n(), events = sum(incident_stroke6),
            hba1c = round(mean(hba1c), 3), fpg = round(mean(glu), 2),
            hgi_mean = round(mean(hgi), 3), hgi_sd = round(sd(hgi), 3),
            hb = round(mean(hb), 1), mcv = round(mean(mcv), 1),
            pct_low_P10 = round(mean(hgi_lo_p10) * 100, 1),
            pct_high_P90 = round(mean(hgi_hi_p90) * 100, 1), .groups = "drop")
print(as.data.frame(pat))

cat("\n[age/sex-adjusted HGI difference vs no_anemia (lm)]\n")
f_adj <- lm(hgi ~ anem_sub + age + factor(sex), data = da)
print(round(summary(f_adj)$coefficients[2:3, ], 4))

## ---- B. tail composition -----------------------------------------------------------
cat("\n==== B. composition of HGI tails by anemia subtype ====\n")
for (g in c("hgi_lo_p10", "hgi_hi_p90")) {
  sub <- da %>% filter(.data[[g]] == 1)
  cat(g, " group: n =", nrow(sub), "\n")
  print(round(prop.table(table(sub$anem_sub)) * 100, 1))
}

## ---- C. crude event rates: subtype x HGI tertile-of-interest ----------------------
cat("\n==== C. crude stroke rates (per 1000 person-years) ====\n")
rate <- da %>% mutate(hgi_g3 = case_when(hgi_lo_p10 == 1 ~ "low_<P10",
                                         hgi_hi_p90 == 1 ~ "high_>=P90",
                                         TRUE ~ "mid")) %>%
  group_by(anem_sub, hgi_g3) %>%
  summarise(n = n(), events = sum(incident_stroke6),
            py = sum(fu_v2), rate_per1000py = round(events / py * 1000, 1),
            .groups = "drop")
print(as.data.frame(rate))

## ---- D. Cox models under subtype-specific exclusions ------------------------------
cat("\n==== D. Cox: high-HGI and low-end tests by exclusion set ====\n")
run_set <- function(dd, label) {
  cat("\n-- ", label, ": n =", nrow(dd), ", events =", sum(dd$incident_stroke6), "\n")
  m_hi  <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj3)), data = dd)
  m_key <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 + hba1c +", adj3)), data = dd)
  m_lo  <- coxph(as.formula(paste(FS, "~ hgi_lo_p10 +", adj3)), data = dd)
  r <- bind_rows(cbind(set = label, analysis = "M3 high HGI",      grab(m_hi,  "hgi_hi_p90")),
                 cbind(set = label, analysis = "key test high HGI", grab(m_key, "hgi_hi_p90")),
                 cbind(set = label, analysis = "key test HbA1c",    grab(m_key, "hba1c")),
                 cbind(set = label, analysis = "low-end (<P10 vs rest)", grab(m_lo, "hgi_lo_p10")))
  print(r); r
}
res_cox <- bind_rows(
  run_set(dmain,                                     "full main set (anchor)"),
  run_set(da,                                        "MCV-classified set"),
  run_set(da %>% filter(anem_sub != "microcytic"),     "exclude microcytic anemia"),
  run_set(da %>% filter(anem_sub != "non_microcytic"), "exclude non-microcytic anemia"),
  run_set(da %>% filter(anem_sub == "no_anemia"),      "exclude all anemia (anchor S3)"))
# adjustment-instead-of-exclusion
m_adj <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 + anem_sub +", adj3)), data = da)
r_adj <- cbind(set = "MCV-classified set", analysis = "M3 + anemia subtype adjusted",
               grab(m_adj, "hgi_hi_p90"))
m_adj_lo <- coxph(as.formula(paste(FS, "~ hgi_lo_p10 + anem_sub +", adj3)), data = da)
r_adj_lo <- cbind(set = "MCV-classified set", analysis = "low-end + anemia subtype adjusted",
                  grab(m_adj_lo, "hgi_lo_p10"))
print(r_adj); print(r_adj_lo)
res_cox <- bind_rows(res_cox, r_adj, r_adj_lo)

## ---- E. five-group mid-reference model under each exclusion set --------------------
cat("\n==== E. five-group model (P25-75 reference) ====\n")
run_5g <- function(dd, label) {
  dd <- dd %>% mutate(hgi_g5 = relevel(factor(hgi_g5), ref = "P25-75"))
  m <- coxph(as.formula(paste(FS, "~ hgi_g5 +", adj3)), data = dd)
  r <- bind_rows(lapply(c("<P10","P10-25","P75-90",">=P90"), function(g)
    cbind(set = label, grab(m, paste0("hgi_g5", g)))))
  print(r); r
}
res_5g <- bind_rows(
  run_5g(dmain,                                     "full main set"),
  run_5g(da %>% filter(anem_sub != "microcytic"),     "exclude microcytic"),
  run_5g(da %>% filter(anem_sub != "non_microcytic"), "exclude non-microcytic"),
  run_5g(da %>% filter(anem_sub == "no_anemia"),      "exclude all anemia"))

## ---- write out ----------------------------------------------------------------------
write.csv(as.data.frame(pat), file.path(OUT, "R4_anemia_mcv_patterns_20261122.csv"),
          row.names = FALSE)
write.csv(as.data.frame(rate), file.path(OUT, "R4_anemia_mcv_rates_20261122.csv"),
          row.names = FALSE)
write.csv(res_cox, file.path(OUT, "R4_anemia_mcv_cox_20261122.csv"), row.names = FALSE)
write.csv(res_5g,  file.path(OUT, "R4_anemia_mcv_fivegroup_20261122.csv"), row.names = FALSE)
cat("\n[done] outputs written to", OUT, "\n")
