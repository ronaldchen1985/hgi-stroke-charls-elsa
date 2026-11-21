####################################################################################
# CHARLS HGI x stroke -- Step 8: E-value + interval-censoring robustness (W6)
#   Part A: E-values (VanderWeele & Ding) for the main / key-test / Fine-Gray / S7a
#           estimates -- HR treated as RR approximation (rare-ish outcome caveat).
#   Part B: interval-censored PH model (icenReg::ic_sp; events interval-censored
#           between the pre-reporting wave and the reporting wave; non-events
#           right-censored at last contact).
# Output: Evalue_icenReg_W6.txt / .rds
####################################################################################

rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(survival); library(dplyr) })
has_icen <- requireNamespace("icenReg", quietly = TRUE)
OUT_DIR <- "D:/CHARLS/HGI_stroke"

d <- read.csv(file.path(OUT_DIR, "analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
d <- d %>% mutate(hgi_hi_p90 = as.integer(hgi >= CUT_P90),
                  shr = (glu / 18) / (1.59 * hba1c - 2.59))
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
adj3 <- "age + sex + urban + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"
dmain <- d %>% filter(complete.cases(across(all_of(c(cov_m3,"hgi","hba1c",
                                                     "incident_stroke6","fu_v2")))))
message("Step8 sample: ", nrow(dmain), " n; events ", sum(dmain$incident_stroke6))
res <- list()

## ---- Part A: E-values -----------------------------------------------------------------
cat("==== A. E-values ====\n")
evalue <- function(hr) ifelse(hr >= 1, hr + sqrt(hr * (hr - 1)),
                              1/hr + sqrt(1/hr * (1/hr - 1)))
m_main <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ hgi_hi_p90 +", adj3)),
                data = dmain)
m_key  <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ hgi_hi_p90 + hba1c +", adj3)),
                data = dmain)
m_s7   <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ hgi_hi_p90 + shr +", adj3)),
                data = dmain)
get_ci <- function(fit) summary(fit)$conf.int["hgi_hi_p90", c(1,3,4)]
ci_main <- get_ci(m_main); ci_key <- get_ci(m_key); ci_s7 <- get_ci(m_s7)
ci_fg   <- c(1.345, 1.102, 1.640)   # Fine-Gray from M5 (cmprsk machine output)
ev <- data.frame(
  analysis = c("Main Cox (M3 full)", "Key test (+HbA1c)",
               "Fine-Gray competing risk", "S7a mutual adjustment for SHR"),
  HR  = round(c(ci_main[1], ci_key[1], ci_fg[1], ci_s7[1]), 3),
  HR_lo = round(c(ci_main[2], ci_key[2], ci_fg[2], ci_s7[2]), 3),
  E_value_point = round(evalue(c(ci_main[1], ci_key[1], ci_fg[1], ci_s7[1])), 2),
  E_value_CI    = round(evalue(c(ci_main[2], ci_key[2], ci_fg[2], ci_s7[2])), 2))
res$Evalues <- ev; print(ev)

## ---- Part B: interval-censored PH model ----------------------------------------------
cat("\n==== B. interval censoring (icenReg) ====\n")
if (has_icen) {
  library(icenReg)
  dic <- dmain %>% mutate(
    rep_yr = case_when(!is.na(stroke_w2) & stroke_w2 == 1 ~ 2013,
                       !is.na(stroke_w3) & stroke_w3 == 1 ~ 2015,
                       !is.na(stroke_w4) & stroke_w4 == 1 ~ 2018,
                       !is.na(stroke_w5) & stroke_w5 == 1 ~ 2020,
                       !is.na(stroke_w6) & stroke_w6 == 1 ~ w6_report_yr),
    prev_yr = case_when(!is.na(stroke_w2) & stroke_w2 == 1 ~ 2011,
                        !is.na(stroke_w3) & stroke_w3 == 1 ~ 2013,
                        !is.na(stroke_w4) & stroke_w4 == 1 ~ 2015,
                        !is.na(stroke_w5) & stroke_w5 == 1 ~ 2018,
                        !is.na(stroke_w6) & stroke_w6 == 1 ~ 2020),
    L = if_else(incident_stroke6 == 1, prev_yr - 2011, last_contact - 2011),
    R = if_else(incident_stroke6 == 1, rep_yr - 2011, Inf)) %>%
    filter(!is.na(L) & !is.na(R) & R > L)
  message("icenReg set: ", nrow(dic), " n; interval-censored events ",
          sum(dic$incident_stroke6))
  fml <- as.formula(paste("cbind(L, R) ~ hgi_hi_p90 +", adj3))
  # semiparametric PH: point estimate only (ic_sp gives no model-based SE)
  t0 <- Sys.time()
  ic <- ic_sp(fml, data = dic, model = "ph")
  message("ic_sp elapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  s <- summary(ic)
  row <- s$summaryParameters["hgi_hi_p90", ]
  res$icenReg_sp <- data.frame(term = "hgi_hi_p90",
                               HR = round(exp(row[["Estimate"]]), 3))
  print(res$icenReg_sp)
  # parametric PH (Weibull + log-logistic): model-based SEs for the CI
  grab_ic <- function(fit) {
    r <- summary(fit)$summaryParameters["hgi_hi_p90", ]
    data.frame(term = "hgi_hi_p90",
               HR = round(exp(r[["Estimate"]]), 3),
               lo = round(exp(r[["Estimate"]] - 1.96 * r[["Std.Error"]]), 3),
               hi = round(exp(r[["Estimate"]] + 1.96 * r[["Std.Error"]]), 3),
               P  = signif(r[["p"]], 3))
  }
  ic_w <- ic_par(fml, data = dic, model = "ph", dist = "weibull")
  res$icenReg_par_weibull <- grab_ic(ic_w); print(res$icenReg_par_weibull)
  ic_l <- ic_par(fml, data = dic, model = "ph", dist = "loglogistic")
  res$icenReg_par_loglogistic <- grab_ic(ic_l); print(res$icenReg_par_loglogistic)
} else cat("skipped (no icenReg)\n")

## ---- export ----------------------------------------------------------------------------
sink(file.path(OUT_DIR, "Evalue_icenReg_W6.txt"))
cat("CHARLS HGI x stroke: E-value + interval-censoring robustness (",
    format(Sys.time(), "%F %T"), ")\n",
    "Sample = W6 main set: ", nrow(dmain), " n; events ",
    sum(dmain$incident_stroke6), "\n", sep = "")
cat("E-value: unmeasured confounder would need associations of this magnitude\n",
    "with BOTH exposure and outcome, beyond the 12 measured covariates, to fully\n",
    "explain away the estimate (VanderWeele & Ding 2017; HR as RR approximation).\n", sep = "")
for (nm in names(res)) { cat("\n##### ", nm, " #####\n"); print(res[[nm]]) }
sink()
saveRDS(res, file.path(OUT_DIR, "Evalue_icenReg_W6.rds"))
message("\nDONE: Evalue_icenReg_W6.txt / .rds")
