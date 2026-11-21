# ==============================================================================
# CHARLS × HGI × 新发卒中 —— 第四步：敏感性分析 S1—S6（按SAP v1.0 第5.7节）
#
# 输入：analytic_hgi_baseline.csv + analytic_delta_hgi.csv（第一步产物）
# 前提：第三步 M1—M6 已完成，主分析锚点 HR = 1.392 (1.122—1.728)
# 每条S分析都与主锚点对比；结果汇总写入 sensitivity_S1_S6.txt/.rds
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(survival); library(dplyr); library(splines) })
has_mice <- requireNamespace("mice", quietly = TRUE)
if (!has_mice) message("提示：S5需要mice包，请先 install.packages(\"mice\")")

OUT_DIR <- "D:/CHARLS/HGI_stroke"
d <- read.csv(file.path(OUT_DIR, "analytic_hgi_baseline.csv"),
              colClasses = c(ID = "character"))
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
d$hgi_hi_p90 <- as.integer(d$hgi >= CUT_P90)
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
dmain <- d %>% filter(complete.cases(across(all_of(
  c(cov_m3, "hgi", "hba1c", "incident_stroke", "fu_years")))))
adj3 <- "age + sex + urban + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"
message("主分析集：", nrow(dmain), " 人；锚点 HR 1.392 (1.122—1.728)")

grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(term = term, HR = round(ci[, 1], 3),
             lo = round(ci[, 3], 3), hi = round(ci[, 4], 3), P = signif(p, 3))
}
Sres <- list()
MAIN <- "主分析锚点: HR 1.392 (1.122—1.728), P=0.003"

# ========== S1. HGI换用非糖尿病人群方程（hgi_ndm） ==============================
cat("\n==================== S1 非糖尿病方程HGI ====================\n")
cut_ndm <- as.numeric(quantile(dmain$hgi_ndm, 0.9, na.rm = TRUE))
dmain$hgi_ndm_hi <- as.integer(dmain$hgi_ndm >= cut_ndm)
s1 <- coxph(as.formula(paste("Surv(fu_years, incident_stroke) ~ hgi_ndm_hi +", adj3)),
            data = dmain)
Sres$S1 <- grab(s1, "hgi_ndm_hi"); print(Sres$S1); cat(MAIN, "\n")

# ========== S2. HGI换用已发表CHARLS方程（hgi_pub = HbA1c − 0.013×FPG − 4.804） ==
cat("\n==================== S2 已发表方程HGI ====================\n")
cut_pub <- as.numeric(quantile(dmain$hgi_pub, 0.9, na.rm = TRUE))
dmain$hgi_pub_hi <- as.integer(dmain$hgi_pub >= cut_pub)
s2 <- coxph(as.formula(paste("Surv(fu_years, incident_stroke) ~ hgi_pub_hi +", adj3)),
            data = dmain)
Sres$S2 <- grab(s2, "hgi_pub_hi"); print(Sres$S2); cat(MAIN, "\n")

# ========== S3. 排除贫血者（HbA1c失真人群；检验J型左臂） ========================
cat("\n==================== S3 排除贫血 ====================\n")
dnoa <- dmain %>% filter(is.na(anemia) | anemia == 0)
message("排除贫血后：", nrow(dnoa), " 人，事件 ", sum(dnoa$incident_stroke),
        "（排除 ", nrow(dmain) - nrow(dnoa), " 人）")
s3 <- coxph(as.formula(paste("Surv(fu_years, incident_stroke) ~ hgi_hi_p90 +", adj3)),
            data = dnoa)
Sres$S3 <- grab(s3, "hgi_hi_p90"); print(Sres$S3); cat(MAIN, "\n")
# J型左臂检验：低HGI（<P10）组在排除贫血后是否仍高风险
cut_p10 <- as.numeric(quantile(dmain$hgi, 0.1, na.rm = TRUE))
dnoa$hgi_lo <- as.integer(dnoa$hgi < cut_p10)
s3b <- coxph(as.formula(paste("Surv(fu_years, incident_stroke) ~ hgi_lo +", adj3)),
             data = dnoa)
Sres$S3_lowarm <- grab(s3b, "hgi_lo")
cat("排除贫血后，低HGI(<P10)组的HR（若回落向1→左臂确由贫血驱动）：\n")
print(Sres$S3_lowarm)

# ========== S4. 事件时间区间中点插补（替代波次年份） ============================
cat("\n==================== S4 区间中点插补 ====================\n")
wy <- c("2" = 2013, "3" = 2015, "4" = 2018, "5" = 2020)
prev_wave <- c("2" = 2011, "3" = 2013, "4" = 2015, "5" = 2018)  # 事件发生区间下界
last_obs <- apply(dmain[, paste0("stroke_w", 2:5)], 1, function(r) {
  hit <- which(!is.na(r)); if (length(hit)) wy[max(hit)] else 2011 })
ds4 <- dmain %>% mutate(
  fu_mid = if_else(incident_stroke == 1,
                   (prev_wave[as.character(first_stroke_wave)] +
                      wy[as.character(first_stroke_wave)]) / 2 - 2011,
                   last_obs - 2011))
s4 <- coxph(as.formula(paste("Surv(fu_mid, incident_stroke) ~ hgi_hi_p90 +", adj3)),
            data = ds4)
Sres$S4 <- grab(s4, "hgi_hi_p90"); print(Sres$S4); cat(MAIN, "\n")

# ========== S5. 多重插补（mice，20套） + 无SBP模型（双角度缺失稳健性） ==========
cat("\n==================== S5a 无SBP模型（N回升） ====================\n")
adj_nosbp <- "age + sex + urban + edu + htn + dm + dyslip + heart + crp + tg + hdl"
d_nosbp <- d %>% filter(complete.cases(across(all_of(
  c(setdiff(cov_m3, "sbp"), "hgi", "hba1c", "incident_stroke", "fu_years")))))
message("剔除SBP后完整病例：", nrow(d_nosbp), " 人，事件 ", sum(d_nosbp$incident_stroke))
s5a <- coxph(as.formula(paste("Surv(fu_years, incident_stroke) ~ hgi_hi_p90 +",
                              adj_nosbp)), data = d_nosbp)
Sres$S5a_noSBP <- grab(s5a, "hgi_hi_p90"); print(Sres$S5a_noSBP); cat(MAIN, "\n")

cat("\n==================== S5b mice多重插补（20套） ====================\n")
if (has_mice) {
  library(mice)
  imp_vars <- c("hgi_hi_p90","age","sex","urban","edu","htn","dm","dyslip","heart",
                "sbp","crp","tg","hdl","hba1c","glu","hb","bmi",
                "incident_stroke","fu_years","death")
  dimp <- d[, imp_vars]
  # 结局与随访时间必须参与插补模型（否则插补有偏）
  ini <- mice(dimp, maxit = 0, printFlag = FALSE)
  meth <- ini$meth; meth[c("sex","urban","htn","dm","dyslip","heart",
                           "incident_stroke","death","hgi_hi_p90")] <- "logreg"
  meth["hgi_hi_p90"] <- ""   # 由hgi派生的二分类不插补（hgi本身无缺失）
  imp <- mice(dimp, m = 20, maxit = 5, meth = meth, printFlag = FALSE, seed = 20260903)
  fits_imp <- with(imp, coxph(as.formula(
    paste("Surv(fu_years, incident_stroke) ~ hgi_hi_p90 +", adj3))))
  pooled <- pool(fits_imp)
  s5b <- summary(pooled, conf.int = TRUE)
  row_hgi <- s5b[s5b$term == "hgi_hi_p90", ]
  Sres$S5b_mice <- data.frame(term = "hgi_hi_p90",
                              HR = round(exp(row_hgi$estimate), 3),
                              lo = round(exp(row_hgi$`2.5 %`), 3),
                              hi = round(exp(row_hgi$`97.5 %`), 3),
                              P  = signif(row_hgi$p.value, 3))
  print(Sres$S5b_mice); cat(MAIN, "\n")
} else cat("跳过（未安装mice）\n")

# ========== S6. ΔHGI子分析（2011→2015变化；结局窗口2015—2020） ==================
cat("\n==================== S6 ΔHGI子分析 ====================\n")
f_delta <- file.path(OUT_DIR, "analytic_delta_hgi.csv")
if (file.exists(f_delta)) {
  dd <- read.csv(f_delta, colClasses = c(ID = "character"))
  # 关键：ΔHGI用2011与2015两期血样构建，2013波卒中早于2015血样，必须剔除
  dd <- dd %>% filter(is.na(first_stroke_wave) | first_stroke_wave >= 4)
  dd <- dd %>% mutate(
    ev15 = as.integer((!is.na(stroke_w4) & stroke_w4 == 1) |
                        (!is.na(stroke_w5) & stroke_w5 == 1)),
    # 事件时间自2015年起算；无事件者在2015后最后观测波删失
    last15 = apply(dd[, c("stroke_w4","stroke_w5")], 1, function(r) {
      hit <- which(!is.na(r)); if (length(hit)) c(2018, 2020)[max(hit)] else NA_real_ }),
    t15 = if_else(ev15 == 1,
                  if_else(!is.na(stroke_w4) & stroke_w4 == 1, 2018 - 2015, 2020 - 2015),
                  last15 - 2015),
    delta_sd = as.numeric(scale(delta_hgi)),
    delta_q4 = as.integer(delta_hgi >= quantile(delta_hgi, 0.75, na.rm = TRUE)))
  dd <- dd %>% filter(!is.na(t15) & t15 > 0)
  message("S6样本（2015后无卒中史且两期血样齐全）：", nrow(dd), " 人，2015后事件 ",
          sum(dd$ev15))
  s6a <- coxph(as.formula(paste("Surv(t15, ev15) ~ delta_sd + hgi +", adj3)), data = dd)
  s6b <- coxph(as.formula(paste("Surv(t15, ev15) ~ delta_q4 + hgi +", adj3)), data = dd)
  Sres$S6_deltaSD <- grab(s6a, "delta_sd")
  Sres$S6_deltaQ4 <- grab(s6b, "delta_q4")
  cat("ΔHGI每升高1SD（调整基线HGI）：\n"); print(Sres$S6_deltaSD)
  cat("ΔHGI最高四分位 vs 其余：\n"); print(Sres$S6_deltaQ4)
} else cat("未找到 analytic_delta_hgi.csv，跳过（需先运行第一步脚本的8b节）\n")

# ========== 导出 ================================================================
sink(file.path(OUT_DIR, "sensitivity_S1_S6.txt"))
cat("CHARLS HGI×卒中 敏感性分析（", format(Sys.time(), "%F %T"), "）\n",
    MAIN, "\n", sep = "")
for (nm in names(Sres)) { cat("\n##### ", nm, " #####\n"); print(Sres[[nm]]) }
sink()
saveRDS(Sres, file.path(OUT_DIR, "sensitivity_S1_S6.rds"))
message("\nS1—S6完成，结果已写入 sensitivity_S1_S6.txt/.rds。贴回全部结果做最终核对。")
