####################################################################################
# CHARLS HGI x stroke -- W6 patch: figures + PH/time-varying robustness + W6 tables
# Date: 2026-11-19. Reads analytic_hgi_baseline_w6.csv and formal_results_W6.rds.
# Does NOT modify raw outputs; overwrites Figure2_RCS_W6.png (title fix) and adds
# Figure3_dC_forest_W6.png / Figure4_subgroup_forest_W6.png / 图1_STROBE流程图_W6版.png,
# plus 表1-4 W6版 csv and patch log.
####################################################################################
rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(survival); library(dplyr); library(splines) })
OUT_DIR <- "D:/CHARLS/HGI_stroke"
LOG <- file(file.path(OUT_DIR, "_code_w6", "patch_log_20261119.txt"), open = "wt")
sink(LOG, split = TRUE)

d <- read.csv(file.path(OUT_DIR, "analytic_hgi_baseline_w6.csv"), colClasses = c(ID = "character"))
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
d$hgi_hi_p90 <- as.integer(d$hgi >= CUT_P90)
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
adj3 <- "age + sex + urban + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"
dmain <- d %>% filter(complete.cases(across(all_of(c(cov_m3,"hgi","hba1c",
                                                     "incident_stroke6","fu_v2")))))
cat("main set n =", nrow(dmain), "; events =", sum(dmain$incident_stroke6), "\n")

## ---- 1. PH test + time-varying robustness (W6) --------------------------------------
m3_full <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ hgi_hi_p90 +", adj3)), data = dmain)
cat("\n== PH test (W6, m3_full) ==\n"); print(cox.zph(m3_full))
m_tv <- coxph(as.formula(paste(
  "Surv(fu_v2, incident_stroke6) ~ hgi_hi_p90 +", adj3, "+ tt(htn) + tt(sbp)")),
  data = dmain, tt = function(x, t, ...) x * log(t))
ci_tv <- summary(m_tv)$conf.int["hgi_hi_p90", ]
p_tv  <- summary(m_tv)$coefficients["hgi_hi_p90", "Pr(>|z|)"]
cat(sprintf("\n== Time-varying robustness (W6): hgi_hi_p90 HR = %.3f (%.3f-%.3f), P = %.4f ==\n",
            ci_tv[1], ci_tv[3], ci_tv[4], p_tv))

## ---- 2. Figure 2 RCS (title fix + rug) ----------------------------------------------
m2_lin <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ hgi +", adj3)), data = dmain)
m2_rcs <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ ns(hgi, df = 4) +", adj3)), data = dmain)
b <- ns(dmain$hgi, df = 4)
grid <- seq(quantile(dmain$hgi, .01), quantile(dmain$hgi, .99), length.out = 200)
Xc <- sweep(predict(b, grid), 2, as.numeric(predict(b, median(dmain$hgi))), "-")
co <- coef(m2_rcs); V <- vcov(m2_rcs); idx <- grep("ns\\(hgi", names(co))
eta <- as.numeric(Xc %*% co[idx])
se  <- sqrt(pmax(diag(Xc %*% V[idx, idx, drop = FALSE] %*% t(Xc)), 0))
png(file.path(OUT_DIR, "Figure2_RCS_W6.png"), width = 1800, height = 1400, res = 220)
par(mar = c(4.5, 4.5, 3.2, 1))
plot(grid, exp(eta), type = "n", ylim = c(0.6, 2.0), xlab = "HGI", ylab = "HR (95% CI)",
     main = "HGI and incident stroke dose-response\n(CHARLS 2011-2023, fully adjusted RCS)",
     cex.main = 0.95)
polygon(c(grid, rev(grid)), c(exp(eta - 1.96 * se), rev(exp(eta + 1.96 * se))),
        col = adjustcolor("steelblue", .25), border = NA)
lines(grid, exp(eta), lwd = 2.5, col = "steelblue4")
abline(h = 1, lty = 2, col = "grey40"); abline(v = CUT_P90, lty = 3, col = "firebrick")
text(CUT_P90, 1.9, "P90", pos = 4, col = "firebrick", cex = .9)
rug(dmain$hgi, col = adjustcolor("grey40", 0.25))
dev.off()
cat("\nFigure2_RCS_W6.png redrawn (two-line title, rug added)\n")

## ---- 3. Figures 3/4 from saved rds ----------------------------------------------------
res <- readRDS(file.path(OUT_DIR, "formal_results_W6.rds"))
# 3a. deltaC forest (manuscript Fig 3)
dc <- res$M4dC %>% filter(indicator != "hgi")
dc$indicator <- factor(dc$indicator, levels = c("TyG","SHR","FPG","HbA1c"))
dc <- dc[order(dc$indicator), ]
png(file.path(OUT_DIR, "Figure3_dC_forest_W6.png"), width = 1800, height = 1100, res = 220)
par(mar = c(4.5, 6, 3, 1))
plot(dc$dC_vs_HGI, seq_len(nrow(dc)), xlim = c(-0.008, 0.012), yaxt = "n", ylab = "",
     pch = 19, cex = 1.4, xlab = "Delta C-index vs HGI (95% CI, bootstrap)",
     main = "Discrimination head-to-head (reference: HGI)\nCHARLS 2011-2023", cex.main = 0.95)
segments(dc$lo, seq_len(nrow(dc)), dc$hi, seq_len(nrow(dc)), lwd = 1.5)
axis(2, at = seq_len(nrow(dc)), labels = as.character(dc$indicator), las = 1)
abline(v = 0, lty = 2, col = "grey40")
dev.off()
cat("Figure3_dC_forest_W6.png created (manuscript Fig 3, W6 version)\n")
print(dc)

# 3b. subgroup forest (manuscript Fig 4), fixed title
rd <- res$M6
lab_map <- c("diabetes" = "diabetes", "sex" = "sex", "age65" = "age65", "urban" = "urban")
lab <- with(rd, paste0(subgroup, "=", level, " (n=", n, ", events ", events, ")"))
y <- nrow(rd):1
png(file.path(OUT_DIR, "Figure4_subgroup_forest_W6.png"), width = 1800, height = 1400, res = 220)
par(mar = c(4.5, 16, 3, 1))
plot(rd$HR, y, xlim = c(0.5, 2.5), yaxt = "n", ylab = "", pch = 19, log = "x",
     xlab = "HR (95% CI, log)",
     main = "High HGI (>=P90) and incident stroke:\nsubgroup analyses (CHARLS 2011-2023)", cex.main = 0.95)
segments(rd$lo, y, rd$hi, y)
axis(2, at = y, labels = lab, las = 1, cex.axis = .8); abline(v = 1, lty = 2, col = "grey40")
dev.off()
cat("Figure4_subgroup_forest_W6.png created (manuscript Fig 4, W6 version, title fixed)\n")

## ---- 4. STROBE flow W6 ---------------------------------------------------------------
png(file.path(OUT_DIR, "图1_STROBE流程图_W6版.png"), width = 2000, height = 2500, res = 200)
par(mar = c(1, 1, 1, 1))
plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 12.5))
main_x <- c(0.6, 6.0); excl_x <- c(7.0, 9.9); bh <- 0.85
steps <- list(
  list(t = "2011 baseline sample\nN = 17,705",                              excl = "Excluded (n = 6,173):\nno blood draw / missing biomarkers"),
  list(t = "Blood test with FPG & HbA1c\nN = 11,532",                       excl = "Excluded (n = 227):\nage < 45 years"),
  list(t = "Age >= 45 years\nN = 11,305",                                   excl = "Excluded (n = 267):\nprevalent stroke at baseline"),
  list(t = "Free of stroke at baseline\nN = 11,038",                        excl = "Excluded (n = 387):\nlost to all follow-up waves (2013-2023)"),
  list(t = ">= 1 follow-up outcome record\n(waves 2-6, 2013-2023)\nN = 10,651", excl = "Excluded (n = 8):\nmissing key covariates"),
  list(t = "Final analytic sample\nN = 10,643",                             excl = NULL)
)
ys <- seq(11.6, 1.8, length.out = 6)
for (i in seq_along(steps)) {
  fill <- if (i %in% c(1, 6)) "grey90" else "white"
  if (i == 6) fill <- "#DCE9F5"
  rect(main_x[1], ys[i] - bh/2, main_x[2], ys[i] + bh/2, col = fill, border = "grey30")
  text(mean(main_x), ys[i], steps[[i]]$t, cex = 0.85)
  if (!is.null(steps[[i]]$excl)) {
    rect(excl_x[1], ys[i] - bh/2 + 0.15, excl_x[2], ys[i] + bh/2 + 0.55, col = "grey95", border = "grey50")
    text(mean(excl_x), ys[i] + 0.35, steps[[i]]$excl, cex = 0.72)
    arrows(main_x[2], ys[i], excl_x[1], ys[i] + 0.35, lty = 2, col = "grey40", length = 0.08)
  }
  if (i < 6) arrows(mean(main_x), ys[i] - bh/2, mean(main_x), ys[i + 1] + bh/2, length = 0.1, lwd = 1.4)
}
text(5, 0.7, "Incident stroke 2011-2023: 1,031 (9.7%); all-cause deaths: 815;\ncomplete-case sample for main models: 8,684 (852 events)",
     cex = 0.9, font = 2)
dev.off()
cat("图1_STROBE流程图_W6版.png created\n")

## ---- 5. W6 tables (CSV, GB18030) ------------------------------------------------------
hi <- dmain$hgi_hi_p90 == 1
ev_hi <- sum(dmain$incident_stroke6[hi]); n_hi <- sum(hi)
ev_lo <- sum(dmain$incident_stroke6[!hi]); n_lo <- sum(!hi)
cat(sprintf("\nTable1 event row: overall %d (%.1f%%); low %d (%.1f%%); high %d (%.1f%%)\n",
            nrow(dmain) * 0 + sum(dmain$incident_stroke6), 100 * mean(dmain$incident_stroke6),
            ev_lo, 100 * ev_lo / n_lo, ev_hi, 100 * ev_hi / n_hi))

grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, ]; p <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  sprintf("%.3f (%.3f-%.3f)", ci[1], ci[3], ci[4])
}
m1_1 <- coxph(Surv(fu_v2, incident_stroke6) ~ hgi_hi_p90, data = dmain)
m1_2 <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ hgi_hi_p90 + age + sex + urban + edu")), data = dmain)
d$hgi_hi_05 <- as.integer(d$hgi >= 0.5); dmain$hgi_hi_05 <- as.integer(dmain$hgi >= 0.5)
m1b <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ hgi_hi_05 +", adj3)), data = dmain)
m3k <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ hgi_hi_p90 + hba1c +", adj3)), data = dmain)
m3k5 <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ hgi_hi_05 + hba1c +", adj3)), data = dmain)
tab2 <- data.frame(
  模型 = c("模型1 未调整","模型2 人口学调整","模型3 完全调整","稳健性: 切点≥0.5（完全调整）",
           "关键检验: P90切点+HbA1c","关键检验: 0.5切点+HbA1c"),
  HR_95CI = c(grab(m1_1,"hgi_hi_p90"), grab(m1_2,"hgi_hi_p90"), grab(m3_full,"hgi_hi_p90"),
              grab(m1b,"hgi_hi_05"), grab(m3k,"hgi_hi_p90"), grab(m3k5,"hgi_hi_05")),
  P = c(signif(summary(m1_1)$coefficients["hgi_hi_p90",5],3),
        signif(summary(m1_2)$coefficients["hgi_hi_p90",5],3),
        signif(summary(m3_full)$coefficients["hgi_hi_p90",5],3),
        signif(summary(m1b)$coefficients["hgi_hi_05",5],3),
        signif(summary(m3k)$coefficients["hgi_hi_p90",5],3),
        signif(summary(m3k5)$coefficients["hgi_hi_05",5],3)),
  样本量 = nrow(dmain), 事件数 = sum(dmain$incident_stroke6))
write.csv(tab2, file.path(OUT_DIR, "表2_主模型与关键检验_W6版.csv"), row.names = FALSE, fileEncoding = "GB18030")

tab3 <- res$M6 %>% transmute(
  亚组 = recode(subgroup, diabetes = "糖尿病状态", sex = "性别", age65 = "年龄组", urban = "城乡"),
  分层 = case_when(
    subgroup == "diabetes" & level == 0 ~ "非糖尿病", subgroup == "diabetes" & level == 1 ~ "糖尿病",
    subgroup == "sex" & level == 1 ~ "男", subgroup == "sex" & level == 2 ~ "女",
    subgroup == "age65" & level == 0 ~ "年龄<65岁", subgroup == "age65" & level == 1 ~ "年龄≥65岁",
    subgroup == "urban" & level == 0 ~ "农村", subgroup == "urban" & level == 1 ~ "城市"),
  n = n, events = events,
  `HR (95%CI)` = sprintf("%.3f (%.3f-%.3f)", HR, lo, hi),
  P = signif(P, 3), P_interaction = signif(P_interaction, 3))
write.csv(tab3, file.path(OUT_DIR, "表3_亚组与交互_W6版.csv"), row.names = FALSE, fileEncoding = "GB18030")

tab4 <- data.frame(
  分析 = c("主分析（锚点，M3完全调整）","M5 Fine-Gray竞争风险","S1 HGI换非糖尿病人群方程","S2 HGI换已发表NHANES方程",
           "S3 排除贫血者（高HGI端）","S3附 低HGI左臂（排除贫血后）","S4a 原发表时间规则(fu_v1)",
           "S4b 确诊年细化(fu_v2b)","S5a 剔除SBP扩大样本","S5b mice多重插补(m=20)",
           "S6 ΔHGI每升高1SD（2015-2023窗口）","S6 ΔHGI最高四分位vs其余"),
  HR_95CI = c(grab(m3_full,"hgi_hi_p90"),
              sprintf("%.3f (%.3f-%.3f)", res$M5$SHR[1], res$M5$lo[1], res$M5$hi[1]),
              sprintf("%.3f (%.3f-%.3f)", res$S1$HR, res$S1$lo, res$S1$hi),
              sprintf("%.3f (%.3f-%.3f)", res$S2$HR, res$S2$lo, res$S2$hi),
              sprintf("%.3f (%.3f-%.3f)", res$S3$HR, res$S3$lo, res$S3$hi),
              sprintf("%.3f (%.3f-%.3f)", res$S3_lowarm$HR, res$S3_lowarm$lo, res$S3_lowarm$hi),
              sprintf("%.3f (%.3f-%.3f)", res$S4a_fu_v1$HR, res$S4a_fu_v1$lo, res$S4a_fu_v1$hi),
              sprintf("%.3f (%.3f-%.3f)", res$S4b_fu_v2b$HR, res$S4b_fu_v2b$lo, res$S4b_fu_v2b$hi),
              sprintf("%.3f (%.3f-%.3f)", res$S5a_noSBP$HR, res$S5a_noSBP$lo, res$S5a_noSBP$hi),
              sprintf("%.3f (%.3f-%.3f)", res$S5b_mice$HR, res$S5b_mice$lo, res$S5b_mice$hi),
              sprintf("%.3f (%.3f-%.3f)", res$S6_deltaSD$HR, res$S6_deltaSD$lo, res$S6_deltaSD$hi),
              sprintf("%.3f (%.3f-%.3f)", res$S6_deltaQ4$HR, res$S6_deltaQ4$lo, res$S6_deltaQ4$hi)),
  P = c(signif(summary(m3_full)$coefficients["hgi_hi_p90",5],3), signif(res$M5$P[1],3),
        signif(res$S1$P,3), signif(res$S2$P,3), signif(res$S3$P,3), signif(res$S3_lowarm$P,3),
        signif(res$S4a_fu_v1$P,3), signif(res$S4b_fu_v2b$P,3), signif(res$S5a_noSBP$P,3),
        signif(res$S5b_mice$P,3), signif(res$S6_deltaSD$P,3), signif(res$S6_deltaQ4$P,3)))
write.csv(tab4, file.path(OUT_DIR, "表4_敏感性分析汇总_W6版.csv"), row.names = FALSE, fileEncoding = "GB18030")
cat("\nW6 tables written: 表2/表3/表4 W6版 csv\n")
print(tab2); print(tab3); print(tab4)
sink(); close(LOG)
message("PATCH DONE")
