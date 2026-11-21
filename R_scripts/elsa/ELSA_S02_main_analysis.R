####################################################################################
# ELSA external validation -- Step 2: formal analysis (mirrors CHARLS W6 step 2)
# Input: analytic_elsa_hgi_baseline_w2.csv (complete-case main set) +
#        analytic_elsa_hgi_full_w2.csv (for FG death timing)
# Output: ELSA_results.txt / .rds
####################################################################################
rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(survival); library(dplyr); library(splines) })
has_cmprsk <- requireNamespace("cmprsk", quietly = TRUE)
has_icen   <- requireNamespace("icenReg", quietly = TRUE)
O <- "D:/CHARLS/HGI_stroke/ELSA验证/"
d <- read.csv(paste0(O, "analytic_elsa_hgi_baseline_w2.csv"))
cat("main set:", nrow(d), "n; events", sum(d$incident_stroke), "; deaths",
    sum(d$died == 1), "\n")

covs <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
adj2 <- "age + sex + urban + edu"
adj3 <- paste(covs, collapse = " + ")
FS <- "Surv(fu_mid, incident_stroke)"
grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(term = term, HR = round(ci[,1],3), lo = round(ci[,3],3),
             hi = round(ci[,4],3), P = signif(p,3))
}
res <- list()

## ---- M1 sequentially adjusted ----------------------------------------------------------
cat("\n==== M1 ====\n")
m1 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90")), data = d)
m2 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj2)), data = d)
m3 <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj3)), data = d)
res$M1 <- bind_rows(lapply(list(m1, m2, m3), grab, term = "hgi_hi_p90"), .id = "model")
res$M1$model <- c("m1_unadj", "m2_demo", "m3_full"); print(res$M1)
cat("PH test (m3):\n"); print(cox.zph(m3))

## ---- Key test: + HbA1c ------------------------------------------------------------------
cat("\n==== key test ====\n")
mk <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 + hba1c +", adj3)), data = d)
res$KeyTest <- bind_rows(grab(mk, "hgi_hi_p90"), grab(mk, "hba1c")); print(res$KeyTest)
hr_key <- summary(mk)$conf.int["hgi_hi_p90",1]
p_key  <- summary(mk)$coefficients["hgi_hi_p90","Pr(>|z|)"]
cat("[decision] HR =", round(hr_key,3), "P =", signif(p_key,3), "->",
    ifelse(hr_key > 1.2 && p_key < 0.05, "HOLDS\n", "attenuated\n"))

## ---- RCS ---------------------------------------------------------------------------------
cat("\n==== RCS ====\n")
m_lin <- coxph(as.formula(paste(FS, "~ hgi +", adj3)), data = d)
m_rcs <- coxph(as.formula(paste(FS, "~ ns(hgi, df = 4) +", adj3)), data = d)
lrt <- anova(m_lin, m_rcs); print(lrt)
res$RCS <- data.frame(P_nonlinear = signif(lrt[["Pr(>|Chi|)"]][2], 3)); print(res$RCS)
res$M1_cont <- grab(m_lin, "hgi"); print(res$M1_cont)

## ---- Head-to-head C-index -----------------------------------------------------------------
cat("\n==== head-to-head ====\n")
expos <- c(hgi = "hgi", HbA1c = "hba1c", FPG = "fglu", SHR = "shr", TyG = "tyg")
fits <- lapply(expos, function(x)
  coxph(as.formula(paste(FS, "~", x, "+", adj3)), data = d))
cidx <- sapply(fits, function(f) concordance(f)$concordance)
hr_sd <- sapply(seq_along(expos), function(k)
  exp(coef(fits[[k]])[1] * sd(d[[expos[k]]], na.rm = TRUE)))
res$M4 <- data.frame(indicator = names(expos), C_index = round(cidx,4),
                     HR_perSD = round(hr_sd,3)); print(res$M4)

## ---- Fine-Gray ----------------------------------------------------------------------------
cat("\n==== Fine-Gray ====\n")
if (has_cmprsk) {
  library(cmprsk)
  full <- read.csv(paste0(O, "analytic_elsa_hgi_full_w2.csv"))
  fg <- full %>% filter(idauniq %in% d$idauniq)
  # death wave = first iwstat in {5,6}; death year: raxyear else wave-year rule
  iws <- paste0("r", 3:10, "iwstat"); iwys <- paste0("r", 3:10, "iwy")
  dw <- apply(fg[, iws], 1, function(r) { v <- which(r %in% c(5,6)); if (length(v)) min(v) else NA_integer_ })
  fg$death_wave <- ifelse(is.na(dw), NA, dw + 2)
  fg$death_wave[!is.na(fg$w11indout) & fg$w11indout == 95 & is.na(fg$death_wave)] <- 11
  dy <- rep(NA_real_, nrow(fg))
  for (i in seq_len(nrow(fg))) {
    w <- fg$death_wave[i]; if (is.na(w)) next
    if (!is.na(fg$raxyear[i])) { dy[i] <- fg$raxyear[i]; next }
    y_now <- if (w <= 10) fg[[paste0("r", w, "iwy")]][i] else fg$w11_iwy[i]
    st <- fg[[paste0("r", w, "iwstat")]][i]
    if (!is.na(st) && st == 5) { dy[i] <- y_now } else {
      y_prev <- if (w - 1 >= 3) fg[[paste0("r", w - 1, "iwy")]][i] else fg$r2iwy[i]
      dy[i] <- (y_now + y_prev) / 2 }
  }
  fg$death_time <- dy - fg$r2iwy
  fg <- fg %>% mutate(
    fg_event = case_when(
      incident_stroke == 1 & (died == 0 | is.na(death_time) | fu_mid <= death_time) ~ 1,
      died == 1 & (incident_stroke == 0 | (!is.na(death_time) & death_time < fu_mid)) ~ 2,
      TRUE ~ 0),
    fg_time = case_when(fg_event == 1 ~ fu_mid, fg_event == 2 ~ death_time,
                        TRUE ~ fu_mid),
    fg_time = pmax(fg_time, 0.5, na.rm = TRUE))
  cat("FG status (0/1/2):\n"); print(table(fg$fg_event))
  X <- model.matrix(as.formula(paste("~ hgi_hi_p90 +", adj3)), data = fg)[, -1]
  fgm <- crr(fg$fg_time, fg$fg_event, X, failcode = 1, cencode = 0)
  sfg <- summary(fgm)
  res$FG <- data.frame(term = rownames(sfg$coef),
                       SHR = round(exp(sfg$coef[,"coef"]),3),
                       lo = round(exp(sfg$coef[,"coef"] - 1.96*sfg$coef[,"se(coef)"]),3),
                       hi = round(exp(sfg$coef[,"coef"] + 1.96*sfg$coef[,"se(coef)"]),3),
                       P = signif(sfg$coef[,"p-value"],3))
  print(res$FG[res$FG$term == "hgi_hi_p90", ])
}

## ---- S7 mutual SHR -----------------------------------------------------------------------
cat("\n==== S7 mutual SHR ====\n")
d$shr_z <- as.numeric(scale(d$shr))
s7a <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 + shr_z +", adj3)), data = d)
res$S7 <- bind_rows(grab(s7a, "hgi_hi_p90"), grab(s7a, "shr_z")); print(res$S7)
s7i_add <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 + shr_z +", adj3)), data = d)
s7i <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 * shr_z +", adj3)), data = d)
res$S7_int <- data.frame(P_interaction = signif(anova(s7i_add, s7i)[["Pr(>|Chi|)"]][2], 3))
print(res$S7_int)
cat("cor(hgi, shr):", round(cor(d$hgi, d$shr), 3), "\n")

## ---- S8 interval censoring -----------------------------------------------------------------
cat("\n==== S8 interval censoring ====\n")
if (has_icen) {
  library(icenReg)
  d$L <- ifelse(d$incident_stroke == 1, d$event_year_prev - d$r2iwy,
                d$last_contact - d$r2iwy)
  d$R <- ifelse(d$incident_stroke == 1, d$event_year_rep - d$r2iwy, Inf)
  dic <- d %>% filter(!is.na(L) & !is.na(R) & R > L)
  fml <- as.formula(paste("cbind(L, R) ~ hgi_hi_p90 +", adj3))
  ic_sp <- ic_sp(fml, data = dic, model = "ph")
  ssp <- summary(ic_sp)$summaryParameters["hgi_hi_p90", ]
  ic_w <- ic_par(fml, data = dic, model = "ph", dist = "weibull")
  sw_ <- summary(ic_w)$summaryParameters["hgi_hi_p90", ]
  res$S8 <- data.frame(
    model = c("semiparametric (ic_sp)", "Weibull PH (ic_par)"),
    HR = round(exp(c(ssp[["Estimate"]], sw_[["Estimate"]])), 3),
    lo = c(NA, round(exp(sw_[["Estimate"]] - 1.96 * sw_[["Std.Error"]]), 3)),
    hi = c(NA, round(exp(sw_[["Estimate"]] + 1.96 * sw_[["Std.Error"]]), 3)),
    P  = c(NA, signif(sw_[["p"]], 3)))
  print(res$S8)
}

## ---- Subgroups ------------------------------------------------------------------------------
cat("\n==== subgroups ====\n")
d$diab_grp <- ifelse(!is.na(d$dm) & d$dm == 1, 1, 0)
d$age65 <- as.integer(d$age >= 65)
run_sub <- function(var) {
  lev <- sort(unique(d[[var]]))
  adj <- paste(setdiff(covs, var), collapse = " + ")
  out <- lapply(lev, function(L) {
    dd <- d[d[[var]] == L, ]
    m <- try(coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", adj)), data = dd), silent = TRUE)
    if (inherits(m, "try-error")) return(NULL)
    g <- grab(m, "hgi_hi_p90"); g$level <- L; g$n <- nrow(dd)
    g$events <- sum(dd$incident_stroke); g
  }) %>% bind_rows()
  m_add <- coxph(as.formula(paste(FS, "~ hgi_hi_p90 +", var, "+", adj)), data = d)
  m_int <- try(coxph(as.formula(paste(FS, "~ hgi_hi_p90 *", var, "+", adj)), data = d), silent = TRUE)
  out$P_interaction <- signif(anova(m_add, m_int)[["Pr(>|Chi|)"]][2], 3)
  out$subgroup <- var; out
}
res$Subgroups <- bind_rows(lapply(c("diab_grp", "sex", "age65", "urban"), run_sub))
print(res$Subgroups)

## ---- E-value ---------------------------------------------------------------------------------
evalue <- function(hr) ifelse(hr >= 1, hr + sqrt(hr * (hr - 1)),
                              1/hr + sqrt(1/hr * (1/hr - 1)))
ci_m3 <- summary(m3)$conf.int["hgi_hi_p90", c(1,3,4)]
res$Evalue <- data.frame(HR = round(ci_m3[1],3),
                         E_point = round(evalue(ci_m3[1]), 2),
                         E_CI = round(evalue(ci_m3[2]), 2))
print(res$Evalue)

## ---- export ------------------------------------------------------------------------------------
sink(paste0(O, "ELSA_results.txt"))
cat("ELSA external validation results (", format(Sys.time(), "%F %T"), ")\n",
    "Baseline Wave 2 (2004/05); main set ", nrow(d), " n; events ",
    sum(d$incident_stroke), "; deaths ", sum(d$died == 1),
    "; median fu ", round(median(d$fu_years_total), 1), " y\n", sep = "")
for (nm in names(res)) { cat("\n##### ", nm, " #####\n"); print(res[[nm]]) }
sink()
saveRDS(res, paste0(O, "ELSA_results.rds"))
message("DONE")
