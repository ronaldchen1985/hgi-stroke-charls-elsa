####################################################################################
# CHARLS HGI x stroke -- W6 extended formal analysis (M1-M6 + sensitivities)
# Input: analytic_hgi_baseline_w6.csv (follow-up 2011-2023; 1031 events; 815 deaths)
# Primary time rule = fu_v2 (interval midpoint + censoring at last observed wave,
#   consistent with the Methods text of the manuscript).
#   fu_v1  = extended published rule (event = reporting-wave year; all non-events
#            censored at 2023) as comparison anchor.
#   fu_v2b = diagnosis-year refined variant.
# Key test decision rule unchanged: HR>1.2 and P<0.05. Seed 20260903.
####################################################################################

rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(survival); library(dplyr); library(splines) })
has_cmprsk <- requireNamespace("cmprsk", quietly = TRUE)
has_mice   <- requireNamespace("mice", quietly = TRUE)
OUT_DIR <- "D:/CHARLS/HGI_stroke"

## ---- 0. Anchor validation with the ORIGINAL csv and original time rule ----------
d0 <- read.csv(file.path(OUT_DIR, "analytic_hgi_baseline.csv"), colClasses = c(ID = "character"))
CUT0 <- as.numeric(quantile(d0$hgi, 0.9, na.rm = TRUE))
d0$hgi_hi_p90 <- as.integer(d0$hgi >= CUT0)
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
d0m <- d0 %>% filter(complete.cases(across(all_of(c(cov_m3,"hgi","hba1c","incident_stroke","fu_years")))))
adj3 <- "age + sex + urban + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"
m0  <- coxph(as.formula(paste("Surv(fu_years, incident_stroke) ~ hgi_hi_p90 +", adj3)), data = d0m)
m0k <- coxph(as.formula(paste("Surv(fu_years, incident_stroke) ~ hgi_hi_p90 + hba1c +", adj3)), data = d0m)
cat("[Anchor-check] original main HR = ", round(summary(m0)$conf.int["hgi_hi_p90",1],3),
    " (expect 1.392); key test HR = ",
    round(summary(m0k)$conf.int["hgi_hi_p90",1],3), " (expect 1.333)\n\n", sep = "")

## ---- 1. W6 data ------------------------------------------------------------------
d <- read.csv(file.path(OUT_DIR, "analytic_hgi_baseline_w6.csv"), colClasses = c(ID = "character"))
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
d <- d %>% mutate(
  hgi_hi_p90 = as.integer(hgi >= CUT_P90),
  hgi_hi_05  = as.integer(hgi >= 0.5),
  shr = (glu / 18) / (1.59 * hba1c - 2.59),
  tyg = log(tg * glu / 2))
dmain <- d %>% filter(complete.cases(across(all_of(c(cov_m3,"hgi","hba1c",
                                                     "incident_stroke6","fu_v2"))))) %>%
  mutate(age65 = as.integer(age >= 65))
message("W6 main complete-case set: ", nrow(dmain), " n; events ",
        sum(dmain$incident_stroke6), " (", round(mean(dmain$incident_stroke6)*100,1), "%); deaths ",
        sum(dmain$death6 == 1, na.rm = TRUE))
adj2 <- "age + sex + urban + edu"

grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(term = term, HR = round(ci[,1],3), lo = round(ci[,3],3),
             hi = round(ci[,4],3), P = signif(p,3))
}
res <- list()
FS <- function(tv) paste("Surv(", tv, ", incident_stroke6)")

## ---- M1 main analysis --------------------------------------------------------------
cat("\n==== M1 main (fu_v2) ====\n")
m1_1 <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_p90")), data = dmain)
m1_2 <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_p90 +", adj2)), data = dmain)
m1_3 <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_p90 +", adj3)), data = dmain)
res$M1 <- bind_rows(lapply(list(m1_1, m1_2, m1_3), grab, term = "hgi_hi_p90"), .id = "model")
res$M1$model <- c("m1_unadj","m2_demo","m3_full"); print(res$M1)
cat("PH test (m3):\n"); print(cox.zph(m1_3))
m1b_3 <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_05 +", adj3)), data = dmain)
res$M1b <- grab(m1b_3, "hgi_hi_05"); print(res$M1b)

## ---- M3 key test: + HbA1c ----------------------------------------------------------
cat("\n==== M3 key test ====\n")
m3_a1c   <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_p90 + hba1c +", adj3)), data = dmain)
m3_a1c_b <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_05  + hba1c +", adj3)), data = dmain)
res$M3 <- bind_rows(grab(m3_a1c,"hgi_hi_p90"), grab(m3_a1c,"hba1c"),
                    grab(m3_a1c_b,"hgi_hi_05"), grab(m3_a1c_b,"hba1c"))
res$M3$model <- rep(c("P90+A1c","0.5+A1c"), each = 2); print(res$M3)
hr_key <- summary(m3_a1c)$conf.int["hgi_hi_p90",1]
p_key  <- summary(m3_a1c)$coefficients["hgi_hi_p90","Pr(>|z|)"]
cat("[Prespecified decision] HR =", round(hr_key,3), ", P =", signif(p_key,3), "->",
    ifelse(hr_key > 1.2 && p_key < 0.05, " independent association HOLDS\n",
           ifelse(p_key < 0.05, " significant but HR<=1.2\n", " attenuated to non-significant\n")))

## ---- M2 RCS ------------------------------------------------------------------------
cat("\n==== M2 RCS ====\n")
m2_lin <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi +", adj3)), data = dmain)
m2_rcs <- coxph(as.formula(paste(FS("fu_v2"), "~ ns(hgi, df = 4) +", adj3)), data = dmain)
lrt2 <- anova(m2_lin, m2_rcs); print(lrt2)
p_nl <- lrt2[["Pr(>|Chi|)"]][2]
dn <- filter(dmain, diabetes == 0)
m_nl  <- coxph(Surv(fu_v2, incident_stroke6) ~ hgi + age + sex, data = dn)
m_nlr <- coxph(Surv(fu_v2, incident_stroke6) ~ ns(hgi, df = 4) + age + sex, data = dn)
lrt_ndm <- anova(m_nl, m_nlr); print(lrt_ndm)
res$M2 <- data.frame(overall_P_nonlinear = signif(p_nl,3),
                     nondiabetic_P_nonlinear = signif(lrt_ndm[["Pr(>|Chi|)"]][2],3))
print(res$M2)
b <- ns(dmain$hgi, df = 4)
grid <- seq(quantile(dmain$hgi,.01), quantile(dmain$hgi,.99), length.out = 200)
Xc <- sweep(predict(b, grid), 2, as.numeric(predict(b, median(dmain$hgi))), "-")
co <- coef(m2_rcs); V <- vcov(m2_rcs); idx <- grep("ns\\(hgi", names(co))
eta <- as.numeric(Xc %*% co[idx])
se  <- sqrt(pmax(diag(Xc %*% V[idx, idx, drop = FALSE] %*% t(Xc)), 0))
png(file.path(OUT_DIR, "Figure2_RCS_W6.png"), width = 1800, height = 1400, res = 220)
par(mar = c(4.5,4.5,3,1))
plot(grid, exp(eta), type = "n", ylim = c(0.6, 2.0), xlab = "HGI", ylab = "HR (95% CI)",
     main = "HGI and incident stroke dose-response (CHARLS 2011-2023, fully adjusted RCS)")
polygon(c(grid, rev(grid)), c(exp(eta-1.96*se), rev(exp(eta+1.96*se))),
        col = adjustcolor("steelblue", .25), border = NA)
lines(grid, exp(eta), lwd = 2.5, col = "steelblue4")
abline(h = 1, lty = 2, col = "grey40"); abline(v = CUT_P90, lty = 3, col = "firebrick")
text(CUT_P90, 1.9, "P90", pos = 4, col = "firebrick", cex = .9); dev.off()

## ---- M4 head-to-head C-index --------------------------------------------------------
cat("\n==== M4 head-to-head ====\n")
expos <- c(hgi = "hgi", HbA1c = "hba1c", FPG = "glu", SHR = "shr", TyG = "tyg")
dmain4 <- dmain %>% filter(complete.cases(across(all_of(c("shr","tyg")))))
message("M4 sample ", nrow(dmain4), " n; events ", sum(dmain4$incident_stroke6))
fits4 <- lapply(expos, function(x)
  coxph(as.formula(paste(FS("fu_v2"), "~", x, "+", adj3)), data = dmain4))
cidx <- sapply(fits4, function(f) concordance(f)$concordance)
hr_sd <- sapply(seq_along(expos), function(k)
  exp(coef(fits4[[k]])[1] * sd(dmain4[[expos[k]]], na.rm = TRUE)))
res$M4 <- data.frame(indicator = names(expos), C_index = round(cidx,4), HR_perSD = round(hr_sd,3))
print(res$M4)
set.seed(20260903); B <- 500; n <- nrow(dmain4)
boot <- matrix(NA_real_, B, length(expos), dimnames = list(NULL, names(expos)))
for (b_ in seq_len(B)) {
  ii <- sample.int(n, replace = TRUE); db <- dmain4[ii, ]
  for (k in seq_along(expos)) {
    fb <- try(coxph(as.formula(paste(FS("fu_v2"), "~", expos[k], "+", adj3)), data = db), silent = TRUE)
    if (!inherits(fb,"try-error")) boot[b_, k] <- concordance(fb)$concordance
  }
}
dc <- sweep(boot, 1, boot[,"hgi"], "-")
ci_d <- apply(dc, 2, quantile, probs = c(.025,.975), na.rm = TRUE)
res$M4dC <- data.frame(indicator = names(expos),
                       dC_vs_HGI = round(colMeans(dc, na.rm = TRUE),4),
                       lo = round(ci_d[1,],4), hi = round(ci_d[2,],4))
print(res$M4dC)

## ---- M5 Fine-Gray -------------------------------------------------------------------
cat("\n==== M5 Fine-Gray ====\n")
if (has_cmprsk) {
  library(cmprsk)
  dfg <- dmain %>% filter(!is.na(death6) | incident_stroke6 == 1)
  dfg$death_time[is.na(dfg$death_time) & dfg$death6 == 0] <-
    dfg$fu_v2[dfg$death_time[is.na(dfg$death_time) & dfg$death6 == 0]]
  dfg <- dfg %>% mutate(
    has_stroke = incident_stroke6 == 1,
    has_death  = !is.na(death6) & death6 == 1,
    fg_event = case_when(
      has_stroke & (!has_death | is.na(death_time) | fu_v2 <= death_time) ~ 1,
      has_death & (!has_stroke | (!is.na(death_time) & death_time < fu_v2)) ~ 2,
      TRUE ~ 0),
    fg_time = case_when(fg_event == 1 ~ fu_v2, fg_event == 2 ~ death_time, TRUE ~ fu_v2),
    fg_time = pmax(fg_time, 0.5))
  cat("Status counts (0 censored / 1 stroke / 2 death):\n"); print(table(dfg$fg_event))
  X <- model.matrix(~ hgi_hi_p90 + age + sex + urban + edu + htn + dm + dyslip +
                      heart + sbp + crp + tg + hdl, data = dfg)[, -1]
  fgm <- crr(dfg$fg_time, dfg$fg_event, X, failcode = 1, cencode = 0)
  sfg <- summary(fgm)
  res$M5 <- data.frame(term = rownames(sfg$coef),
                       SHR = round(exp(sfg$coef[,"coef"]),3),
                       lo = round(exp(sfg$coef[,"coef"]-1.96*sfg$coef[,"se(coef)"]),3),
                       hi = round(exp(sfg$coef[,"coef"]+1.96*sfg$coef[,"se(coef)"]),3),
                       P = signif(sfg$coef[,"p-value"],3))
  print(res$M5[res$M5$term == "hgi_hi_p90", ])
} else cat("skipped (no cmprsk)\n")

## ---- M6 subgroups and interaction ---------------------------------------------------
cat("\n==== M6 subgroups ====\n")
subs <- list(diabetes = "diabetes", sex = "sex", age65 = "age65", urban = "urban")
run_sub <- function(var) {
  dmain$sub <- dmain[[var]]
  lev <- sort(unique(dmain$sub))
  adj <- paste(setdiff(strsplit(gsub(" ","",adj3), "\\+")[[1]],
                       c(var, if (var == "diabetes") "dm" else "")), collapse = " + ")
  out <- lapply(lev, function(L) {
    dd <- filter(dmain, sub == L)
    m <- try(coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_p90 +", adj)), data = dd), silent = TRUE)
    if (inherits(m,"try-error")) return(NULL)
    gg <- grab(m, "hgi_hi_p90"); gg$level <- L; gg$n <- nrow(dd)
    gg$events <- sum(dd$incident_stroke6); gg
  }) %>% bind_rows()
  m_add <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_p90 +", var, "+", adj)), data = dmain)
  m_int <- try(coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_p90 *", var, "+", adj)), data = dmain), silent = TRUE)
  out$P_interaction <- signif(anova(m_add, m_int)[["Pr(>|Chi|)"]][2], 3)
  out$subgroup <- var; out
}
res$M6 <- bind_rows(lapply(subs, run_sub)); print(res$M6)
rd <- res$M6
lab <- with(rd, paste0(subgroup, "=", level, " (n=", n, ", events ", events, ")"))
y <- nrow(rd):1
png(file.path(OUT_DIR, "Figure3_forest_W6.png"), width = 1800, height = 1400, res = 220)
par(mar = c(4.5, 16, 2, 1))
plot(rd$HR, y, xlim = c(0.5, 2.5), yaxt = "n", ylab = "", pch = 19, log = "x",
     xlab = "HR (95% CI, log)", main = "High HGI (>=P90) and incident stroke: subgroups (to 2023)")
segments(rd$lo, y, rd$hi, y)
axis(2, at = y, labels = lab, las = 1, cex.axis = .8); abline(v = 1, lty = 2, col = "grey40")
dev.off()

## ---- S1-S6 sensitivities ------------------------------------------------------------
cat("\n==== S1-S6 ====\n")
Sres <- list()
cut_ndm <- as.numeric(quantile(dmain$hgi_ndm, .9, na.rm = TRUE))
dmain$hgi_ndm_hi <- as.integer(dmain$hgi_ndm >= cut_ndm)
s1 <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_ndm_hi +", adj3)), data = dmain)
Sres$S1 <- grab(s1, "hgi_ndm_hi"); print(Sres$S1)

cut_pub <- as.numeric(quantile(dmain$hgi_pub, .9, na.rm = TRUE))
dmain$hgi_pub_hi <- as.integer(dmain$hgi_pub >= cut_pub)
s2 <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_pub_hi +", adj3)), data = dmain)
Sres$S2 <- grab(s2, "hgi_pub_hi"); print(Sres$S2)

dnoa <- dmain %>% filter(is.na(anemia) | anemia == 0)
s3 <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_p90 +", adj3)), data = dnoa)
Sres$S3 <- grab(s3, "hgi_hi_p90"); print(Sres$S3)
cut_p10 <- as.numeric(quantile(dmain$hgi, .1, na.rm = TRUE))
dnoa$hgi_lo <- as.integer(dnoa$hgi < cut_p10)
s3b <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_lo +", adj3)), data = dnoa)
Sres$S3_lowarm <- grab(s3b, "hgi_lo"); print(Sres$S3_lowarm)

s4a <- coxph(as.formula(paste(FS("fu_v1"), "~ hgi_hi_p90 +", adj3)), data = dmain)
Sres$S4a_fu_v1 <- grab(s4a, "hgi_hi_p90"); print(Sres$S4a_fu_v1)
s4b <- coxph(as.formula(paste(FS("fu_v2b"), "~ hgi_hi_p90 +", adj3)), data = dmain)
Sres$S4b_fu_v2b <- grab(s4b, "hgi_hi_p90"); print(Sres$S4b_fu_v2b)

adj_nosbp <- "age + sex + urban + edu + htn + dm + dyslip + heart + crp + tg + hdl"
d_nosbp <- d %>% filter(complete.cases(across(all_of(
  c(setdiff(cov_m3,"sbp"),"hgi","hba1c","incident_stroke6","fu_v2")))))
s5a <- coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_p90 +", adj_nosbp)), data = d_nosbp)
Sres$S5a_noSBP <- grab(s5a, "hgi_hi_p90"); print(Sres$S5a_noSBP)

if (has_mice) {
  library(mice)
  imp_vars <- c("hgi_hi_p90","age","sex","urban","edu","htn","dm","dyslip","heart",
                "sbp","crp","tg","hdl","hba1c","glu","hb","bmi",
                "incident_stroke6","fu_v2","death6")
  dimp <- d[, imp_vars] %>%
    mutate(across(c(sex, urban, htn, dm, dyslip, heart, incident_stroke6, death6,
                    hgi_hi_p90), ~as.integer(round(.))))
  ini <- mice(dimp, maxit = 0, printFlag = FALSE)
  meth <- ini$meth
  meth[c("sex","urban","htn","dm","dyslip","heart","incident_stroke6","death6")] <- "logreg"
  meth["hgi_hi_p90"] <- ""
  imp <- mice(dimp, m = 20, maxit = 5, meth = meth, printFlag = FALSE, seed = 20260903)
  fits_imp <- with(imp, coxph(as.formula(paste(FS("fu_v2"), "~ hgi_hi_p90 +", adj3))))
  pooled <- pool(fits_imp); s5b <- summary(pooled, conf.int = TRUE)
  rh <- s5b[s5b$term == "hgi_hi_p90", ]
  Sres$S5b_mice <- data.frame(term = "hgi_hi_p90",
                              HR = round(exp(rh$estimate),3),
                              lo = round(exp(rh$`2.5 %`),3), hi = round(exp(rh$`97.5 %`),3),
                              P = signif(rh$p.value,3))
  print(Sres$S5b_mice)
}

## S6 delta-HGI with outcome window 2015-2023
f_delta <- file.path(OUT_DIR, "analytic_delta_hgi.csv")
if (file.exists(f_delta)) {
  dd <- read.csv(f_delta, colClasses = c(ID = "character"))
  # delta csv already carries stroke_w4/stroke_w5; only W6 columns must be joined
  dw6 <- d %>% select(ID, stroke_w6, last_contact, w6_report_yr, w6_dxyear)
  adjv <- strsplit(gsub(" ","",adj3),"\\+")[[1]]
  dd <- dd %>% select(-any_of(intersect(adjv, setdiff(names(dd), "ID")))) %>%
    left_join(dw6, by = "ID") %>%
    left_join(d %>% select(ID, all_of(adjv)), by = "ID") %>%
    filter(is.na(first_stroke_wave) | first_stroke_wave >= 4) %>%
    mutate(
      ev15 = as.integer((!is.na(stroke_w4) & stroke_w4 == 1) |
                          (!is.na(stroke_w5) & stroke_w5 == 1) |
                          (!is.na(stroke_w6) & stroke_w6 == 1)),
      rep_yr = case_when(!is.na(stroke_w4) & stroke_w4 == 1 ~ 2018,
                         !is.na(stroke_w5) & stroke_w5 == 1 ~ 2020,
                         !is.na(stroke_w6) & stroke_w6 == 1 ~ w6_report_yr),
      prev_yr = case_when(!is.na(stroke_w4) & stroke_w4 == 1 ~ 2015,
                          !is.na(stroke_w5) & stroke_w5 == 1 ~ 2018,
                          !is.na(stroke_w6) & stroke_w6 == 1 ~ 2020),
      last15 = last_contact,
      t15 = if_else(ev15 == 1,
                    ifelse(!is.na(w6_dxyear) & !is.na(stroke_w6) & stroke_w6 == 1 &
                             w6_dxyear > 2015 & w6_dxyear <= rep_yr,
                           w6_dxyear - 2015, (prev_yr + rep_yr) / 2 - 2015),
                    last15 - 2015),
      delta_sd = as.numeric(scale(delta_hgi)),
      delta_q4 = as.integer(delta_hgi >= quantile(delta_hgi, .75, na.rm = TRUE))) %>%
    filter(!is.na(t15) & t15 > 0)
  message("S6 sample (no stroke pre-2015, both bloods): ", nrow(dd),
          " n; events 2015-2023 ", sum(dd$ev15))
  s6a <- coxph(as.formula(paste("Surv(t15, ev15) ~ delta_sd + hgi +", adj3)), data = dd)
  s6b <- coxph(as.formula(paste("Surv(t15, ev15) ~ delta_q4 + hgi +", adj3)), data = dd)
  Sres$S6_deltaSD <- grab(s6a, "delta_sd"); print(Sres$S6_deltaSD)
  Sres$S6_deltaQ4 <- grab(s6b, "delta_q4"); print(Sres$S6_deltaQ4)
}

## ---- export --------------------------------------------------------------------------
sink(file.path(OUT_DIR, "formal_results_W6.txt"))
cat("CHARLS HGI x stroke W6 formal results (", format(Sys.time(), "%F %T"), ")\n",
    "Follow-up 2011-2023; main set ", nrow(dmain), " n; events ",
    sum(dmain$incident_stroke6), "; deaths ", sum(dmain$death6 == 1, na.rm = TRUE), "\n", sep = "")
cat("Time rules: primary fu_v2 (midpoint + last-contact censoring); S4a fu_v1 (published rule); S4b fu_v2b (dx-year refined)\n\n")
for (nm in names(res)) { cat("\n##### ", nm, " #####\n"); print(res[[nm]]) }
cat("\n\n===== Sensitivities =====\n")
for (nm in names(Sres)) { cat("\n##### ", nm, " #####\n"); print(Sres[[nm]]) }
sink()
saveRDS(c(res, Sres), file.path(OUT_DIR, "formal_results_W6.rds"))
message("\nDONE: formal_results_W6.txt / .rds; Figure2 and Figure3 W6 versions updated.")
