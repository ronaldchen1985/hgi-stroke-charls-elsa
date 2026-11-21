################################################################################
# CHARLS × HGI × 新发卒中 —— W6 扩展：数据构建脚本
# 输入：analytic_hgi_baseline.csv（2011-2020 版）+ CHARLS W6 (2021/2022/2023)
# 输出：analytic_hgi_baseline_w6.csv + build_w6_log.txt + strobe_flow_w6.csv
#
# W6 结构（2026-11-18 逐文件核实）：
#   2021/2022/2023 = 同一波次的正卷与补采，个体仅入选一个 res_wave；
#   卒中 = comb01(da011_8_, zdisafter2020_8_, zdisease_8_)（槽位8=卒中，经
#          da030_8_ "Stroke Complication Treatment" 交叉验证）；
#   确诊年份 = da016_1_8_（有效范围 2012—当年）；
#   死亡 = 2023 Sample_Infor.dta 的 died 标志（877 人，res_wave=确认波次），
#          2021-2023 无死亡年份 → 按预设规则区间中点插补。
################################################################################

rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(haven); library(dplyr) })

ROOT <- "D:/CHARLS"; OUT <- "D:/CHARLS/HGI_stroke"
L <- character(); say <- function(...) { x <- paste0(..., collapse = ""); cat(x, "\n"); L <<- c(L, x) }

num <- function(x) suppressWarnings(as.numeric(x))
g   <- function(d, v) if (v %in% names(d)) num(d[[v]]) else rep(NA_real_, nrow(d))
to_w1_id <- function(id) {
  id <- as.character(id); hit <- nchar(id) == 12 & substr(id, 10, 10) == "0"
  id[hit] <- paste0(substr(id[hit], 1, 9), substr(id[hit], 11, 12)); id
}
comb01 <- function(...) {
  z <- list(...)
  positive <- Reduce(`|`, lapply(z, function(x) !is.na(x) & x == 1))
  observed <- Reduce(`|`, lapply(z, function(x) !is.na(x)))
  ifelse(positive, 1L, ifelse(observed, 0L, NA_integer_))
}

## ---- 0. 基线分析集 ----------------------------------------------------------
d <- read.csv(file.path(OUT, "analytic_hgi_baseline.csv"), colClasses = c(ID = "character"))
say("基线分析集 n = ", nrow(d), "；原版事件 = ", sum(d$incident_stroke, na.rm = TRUE),
    "；原版死亡 = ", sum(d$death == 1, na.rm = TRUE),
    "；死亡未知 = ", sum(is.na(d$death)))

## ---- 1. W6 卒中提取（三批合并） ---------------------------------------------
vars6 <- c("ID", "da011_8_", "zdisafter2020_8_", "zdisease_8_", "da016_1_8_")
w6list <- list()
for (y in c(2021, 2022, 2023)) {
  p <- file.path(ROOT, as.character(y), "Health_Status_and_Functioning.dta")
  dd <- as.data.frame(read_dta(p, col_select = all_of(vars6)))
  validyr <- function(x) { x <- num(x); x[x < 2012 | x > y] <- NA; x }
  z <- data.frame(
    ID = to_w1_id(dd$ID), w6_file_year = y,
    st6 = comb01(g(dd, "da011_8_"), g(dd, "zdisafter2020_8_"), g(dd, "zdisease_8_")),
    dx  = validyr(g(dd, "da016_1_8_")))
  say("W6 ", y, " 批次: 记录 ", nrow(z), "，ID ", length(unique(z$ID)),
      "，卒中阳性 ", sum(z$st6 == 1, na.rm = TRUE),
      "，状态可判 ", sum(!is.na(z$st6)))
  w6list[[as.character(y)]] <- z
}
w6 <- bind_rows(w6list)
# 同人多批次：正报告取最早批次年；确诊年取该最早正报告批次内的有效确诊年
w6c <- w6 %>% group_by(ID) %>% arrange(ID, w6_file_year) %>%
  summarise(
    w6_obs      = any(!is.na(st6)),
    w6_pos      = any(st6 == 1, na.rm = TRUE) & any(!is.na(st6)),
    w6_int_yr   = max(w6_file_year),
    w6_rep_yr   = if (any(st6 == 1, na.rm = TRUE)) min(w6_file_year[which(st6 == 1)]) else NA_real_,
    w6_dx_yr    = { k <- which(st6 == 1)[1]
                    if (length(k) && !is.na(w6_file_year[k])) dx[k] else NA_real_ },
    .groups = "drop") %>% as.data.frame()
say("W6 三批合并：唯一受访者 ", nrow(w6c), "；状态可判 ", sum(w6c$w6_obs),
    "；卒中正报告 ", sum(w6c$w6_pos, na.rm = TRUE))

## ---- 2. W6 死亡（Sample_Infor 2023） ----------------------------------------
s <- as.data.frame(read_dta(file.path(ROOT, "2023", "Sample_Infor.dta"),
                            col_select = c("ID", "res_wave", "died", "iyear")))
s$ID <- to_w1_id(s$ID)
s6 <- s %>% group_by(ID) %>%
  summarise(died6 = max(num(died), na.rm = TRUE),
            confirm_yr = suppressWarnings(min(num(res_wave), na.rm = TRUE)),
            .groups = "drop") %>% as.data.frame()
s6$died6[!is.finite(s6$died6)] <- NA
say("Sample_Infor(W6): 唯一 ID ", nrow(s6), "；died=1 ", sum(s6$died6 == 1, na.rm = TRUE))

## ---- 3. 合并到基线分析集并更新结局 -------------------------------------------
d <- d %>% left_join(w6c, by = "ID") %>% left_join(s6, by = "ID")

# 3.1 W6 卒中观测（对基线后无卒中者）
d$stroke_w6 <- NA_integer_
d$stroke_w6[!is.na(d$w6_obs)] <- as.integer(d$w6_pos[!is.na(d$w6_obs)])

# 3.2 分波观测年：w2-w5 的观测年沿用原csv的stroke_wX非缺失；W6观测年=w6_file_year对应批
#     （w6c的w6_obs已汇总；这里只需要"最后观测年"）
wy <- c("2" = 2013, "3" = 2015, "4" = 2015, "4" = 2018, "5" = 2020)
wy <- c("2" = 2013, "3" = 2015, "4" = 2018, "5" = 2020)
last_obs_w15 <- apply(as.matrix(d[, paste0("stroke_w", 2:5)]), 1, function(r) {
  hit <- which(!is.na(r)); if (length(hit)) wy[max(hit)] else NA_real_ })
d$last_contact <- pmax(last_obs_w15, 2011, na.rm = TRUE)
# W6 健康接触（受访）年：该人在任一批次出现（含仅状态未知者）
d$last_contact[!is.na(d$w6_int_yr)] <- pmax(d$last_contact[!is.na(d$w6_int_yr)],
                                            d$w6_int_yr[!is.na(d$w6_int_yr)], na.rm = TRUE)

# 3.3 新发卒中更新
d$n_wave_obs6 <- rowSums(!is.na(d[, c(paste0("stroke_w", 2:5), "stroke_w6")]))
d$first_stroke_wave6 <- d$first_stroke_wave
d$first_stroke_wave6[is.na(d$first_stroke_wave6) & d$stroke_w6 == 1] <- 6
d$incident_stroke6 <- ifelse(d$n_wave_obs6 == 0, NA_real_,
                             ifelse(!is.na(d$first_stroke_wave6), 1, 0))
d$w6_report_yr <- d$w6_rep_yr
d$w6_dxyear <- d$w6_dx_yr

say("结局更新：原事件 ", sum(d$incident_stroke == 1, na.rm = TRUE),
    " → 新事件 ", sum(d$incident_stroke6 == 1, na.rm = TRUE),
    "（W6 新增 ", sum(d$first_stroke_wave6 == 6, na.rm = TRUE), "）")
say("W6 新增事件的报告年分布: ",
    paste(names(table(d$w6_report_yr[d$first_stroke_wave6 == 6])),
          table(d$w6_report_yr[d$first_stroke_wave6 == 6]), sep = "=", collapse = "; "))
say("W6 新增事件确诊年(da016)有效: ", sum(!is.na(d$w6_dxyear[d$first_stroke_wave6 == 6])))
say("观测波次数分布(w2-w6): ", paste(names(table(d$n_wave_obs6)),
                                     table(d$n_wave_obs6), sep = "=", collapse = "; "))

# 3.4 死亡更新（superseded 规则：确认年早于最后健康接触年 → 丢弃该死亡记录）
# 2026-11-18 交叉验证：Sample_Infor 匹配 8635/10643；died=1 共 418 人，全部为
# 原 death==0（2020 波后死亡）；原 death==1（397 例 exit 死亡）与 W6 died 零重叠
# → death6 = 原1 或 (died6==1 且未被推翻)
d$confirm_yr[!is.finite(d$confirm_yr)] <- NA
d$death_superseded <- !is.na(d$confirm_yr) & !is.na(d$last_contact) &
  d$confirm_yr < d$last_contact
say("死亡记录被后续健康接触推翻: ", sum(d$death_superseded, na.rm = TRUE))
w6_died <- !is.na(d$died6) & d$died6 == 1 & !d$death_superseded
d$death6 <- d$death
d$death6[w6_died] <- 1
d$death6[is.na(d$death6) & d$n_wave_obs6 > 0] <- 0
d$death_confirm_yr <- ifelse(w6_died, d$confirm_yr, NA_real_)
say("死亡更新: 原 1 =", sum(d$death == 1, na.rm = TRUE),
    " → 新 1 = ", sum(d$death6 == 1, na.rm = TRUE),
    "；W6 died 标志命中 = ", sum(w6_died, na.rm = TRUE),
    "（其中与原 death=1 重叠 ", sum(w6_died & d$death == 1, na.rm = TRUE),
    "；净新增 ", sum(w6_died & d$death != 1, na.rm = TRUE), "）")
say("死亡更新后：存活 = ", sum(d$death6 == 0, na.rm = TRUE),
    "；未知 = ", sum(is.na(d$death6)))
say("W6 死亡确认年分布: ",
    paste(names(table(d$confirm_yr[w6_died])), table(d$confirm_yr[w6_died]),
          sep = "=", collapse = "; "))

# 3.5 随访时间（两种规则）
prev_wave_y <- c("2" = 2011, "3" = 2013, "4" = 2015, "5" = 2018, "6" = 2020)
rep_yr_all <- ifelse(d$first_stroke_wave6 == 6, d$w6_report_yr, NA_real_)
rep_yr_all[d$first_stroke_wave6 %in% 2:5] <-
  wy[as.character(d$first_stroke_wave6[d$first_stroke_wave6 %in% 2:5])]
d$fu_v1 <- ifelse(d$incident_stroke6 == 1, rep_yr_all - 2011, 2023 - 2011)
# V2 = 中点 + 按最后观测删失（与稿件Methods文字一致的规则）
d$fu_v2 <- ifelse(d$incident_stroke6 == 1,
                  (prev_wave_y[as.character(d$first_stroke_wave6)] + rep_yr_all) / 2 - 2011,
                  d$last_contact - 2011)
# V2b = 确诊年细化（确诊年须晚于W6前最后无卒中观测年且≤报告年，否则回退V2）
last_obs_pre_w6 <- pmax(last_obs_w15, 2011, na.rm = TRUE)
dx_use <- d$w6_dxyear
dx_use[!(dx_use > last_obs_pre_w6 & dx_use <= rep_yr_all)] <- NA
d$fu_v2b <- ifelse(d$incident_stroke6 == 1,
                   ifelse(!is.na(dx_use), dx_use - 2011, d$fu_v2), d$fu_v2)
say("V2b 确诊年细化可用于: ", sum(!is.na(dx_use) & d$incident_stroke6 == 1), " 例事件")
say("时间规则核对: V2 非事件删失年分布: ",
    paste(names(table(d$fu_v2[d$incident_stroke6 == 0])),
          table(d$fu_v2[d$incident_stroke6 == 0]), sep = "=", collapse = "; "))

# 3.6 死亡时间（Fine-Gray 用）：确切年(2013/2020 exit) → 否则中点(最后存活接触, 确认年)
death_yr_exit <- pmax(num(d$deathyr_w2), num(d$deathyr_w5), na.rm = TRUE)
death_yr_exit[!is.finite(death_yr_exit)] <- NA
death_yr_exit[death_yr_exit < 2011 | death_yr_exit > 2023] <- NA
last_alive <- ifelse(!is.na(d$deathyr_w2) | !is.na(d$deathyr_w5),
                     d$last_contact, d$last_contact)
# 确切年死亡者的"最后存活接触"：exit文件口径=出现于exit即死亡；近似用last_contact
d$death_time <- NA_real_
ix <- which(d$death6 == 1)
d$death_time[ix] <- ifelse(!is.na(death_yr_exit[ix]),
                           death_yr_exit[ix] - 2011,
                           (last_alive[ix] + ifelse(!is.na(d$death_confirm_yr[ix]),
                                                    d$death_confirm_yr[ix], 2020)) / 2 - 2011)
d$death_time[ix] <- pmax(d$death_time[ix],
                         pmax(last_alive[ix], 2011, na.rm = TRUE) - 2011 + 0.25, na.rm = TRUE)
say("死亡时间: 确切年 ", sum(d$death6 == 1 & !is.na(death_yr_exit), na.rm = TRUE),
    "；中点插补 ", sum(d$death6 == 1 & is.na(death_yr_exit), na.rm = TRUE))

## ---- 4. STROBE 流程对照 ------------------------------------------------------
d$id_original <- !is.na(d$ID)
flow6 <- data.frame(
  step = c("有2011年FPG与HbA1c", "年龄≥45岁", "排除基线卒中", "至少1次随访(w2-w6)", "关键协变量非缺失"),
  n = c(NA, NA, NA, sum(!is.na(d$incident_stroke6)), NA))
say("W6 版随访可判: ", sum(!is.na(d$incident_stroke6)),
    "（原版 ", sum(!is.na(d$incident_stroke)), "）；分析集不变 n=", nrow(d))

## ---- 5. 导出 -----------------------------------------------------------------
keep <- c("ID","sex","edu","marital","age","urban","glu","hba1c","hb","hct","mcv","plt",
          "wbc","crp","tg","tc","hdl","ldl","crea","bun","ua","htn","dm","dyslip","heart",
          "stroke_w1","dm_med","sbp","dbp","pulse","height","weight","waist","bmi",
          paste0("stroke_w", 2:5), "stroke_w6",
          "deathyr_w2","deathyr_w5","death","death6","n_wave_obs","n_wave_obs6",
          "incident_stroke","incident_stroke6","first_stroke_wave","first_stroke_wave6",
          "w6_report_yr","w6_dxyear","diabetes","hgi_pred","hgi","hgi_ndm","hgi_pub",
          "hgi_tertile","anemia","last_contact","fu_years","fu_v1","fu_v2","fu_v2b","death_time")
keep <- intersect(keep, names(d))
write.csv(d[, keep], file.path(OUT, "analytic_hgi_baseline_w6.csv"), row.names = FALSE, na = "")
writeLines(L, file.path(OUT, "build_w6_log.txt"), useBytes = TRUE)
write.csv(flow6, file.path(OUT, "strobe_flow_w6.csv"), row.names = FALSE)
say("完成：analytic_hgi_baseline_w6.csv 已写出。")
