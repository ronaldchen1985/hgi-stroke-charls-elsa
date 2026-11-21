################################################################################
# CHARLS × HGI × 新发卒中 —— 第一步：数据提取与快速筛查脚本
# 用途：从CHARLS原始数据提取分析变量 → 构建分析样本 → 构建HGI → 快速筛查
# 对应操作手册：《CHARLS_HGI_数据提取操作手册》
# 作者按：CHARLS各版本文件的变量名可能因发布版本略有差异，本脚本用"按变量名正则
#        +变量标签关键词"双重匹配自动识别；每处识别结果都会打印出来，请务必逐条
#        核对打印结果是否与codebook一致，不一致时在【配置区】手工指定。
################################################################################

# ========== 0. 配置区（按你的实际文件路径与变量名修改） =========================
DATA_DIR  <- "D:/CHARLS"          # 存放CHARLS各波次解压后数据的目录
OUT_DIR   <- "D:/CHARLS/HGI_stroke"  # 输出目录
MIN_AGE   <- 45
FOLLOW_WAVES <- c(2, 3, 4, 5)     # 结局随访波次（2013/2015/2018/2020）

# 各数据文件路径（已按你实际下载的文件名配置，2026-09版）
# 注意三个命名陷阱：
#   ① 2011年文件名全小写，且 biomarkers.dta = 体格测量、Blood_20140429.dta = 血液数据（别搞反）
#   ② 2013年起文件名改CamelCase
#   ③ 2020年死亡数据在 Exit_Module.dta（不叫Exit_Interview）
#   ④ 2015/2018年没有独立的Exit_Interview文件（已官网确认）：死亡信息只有
#     2013年Exit_Interview（含确切死亡日期）+ 2020年Exit_Module两路，脚本自动跳过缺失文件
FILES <- list(
  # ---- 2011（W1，全小写文件名）----
  w1_blood   = file.path(DATA_DIR, "2011/Blood_20140429.dta"),            # 血液数据（核心）
  w1_phys    = file.path(DATA_DIR, "2011/biomarkers.dta"),                # 体格测量：血压/身高/体重/腰围
  w1_demo    = file.path(DATA_DIR, "2011/demographic_background.dta"),    # 人口学
  w1_health  = file.path(DATA_DIR, "2011/health_status_and_functioning.dta"), # 慢病+卒中基线+吸烟饮酒
  w1_psu     = file.path(DATA_DIR, "2011/psu.dta"),                       # 城乡属性（urban_rural）备选来源
  # ---- 2013（W2，CamelCase）----
  w2_health  = file.path(DATA_DIR, "2013/Health_Status_and_Functioning.dta"),
  w2_exit    = file.path(DATA_DIR, "2013/Exit_Interview.dta"),
  # ---- 2015（W3）----
  w3_health  = file.path(DATA_DIR, "2015/Health_Status_and_Functioning.dta"),
  w3_blood   = file.path(DATA_DIR, "2015/Blood.dta"),             # 2015血样（已补下）→ ΔHGI子分析用
                                                                # 若实际文件名不同，改这一行即可
  # ---- 2018（W4）----
  w4_health  = file.path(DATA_DIR, "2018/Health_Status_and_Functioning.dta"),
  # ---- 2020（W5）----
  w5_health  = file.path(DATA_DIR, "2020/Health_Status_and_Functioning.dta"),
  w5_exit    = file.path(DATA_DIR, "2020/Exit_Module.dta")       # 2020死亡在此文件
)

# ========== 1. 依赖包 ==========================================================
pkgs <- c("haven", "dplyr", "tidyr", "purrr", "stringr",
          "survival", "tableone")
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

# rms为可选包：装不上（R版本问题）时自动改用 splines::ns() 回退方案，
# RCS的拟合与检验结果与rms等价，仅画图语法不同
has_rms <- requireNamespace("rms", quietly = TRUE)
if (!has_rms) {
  if (!requireNamespace("splines", quietly = TRUE)) stop("splines包缺失（R自带，不应发生）")
  library(splines)
  message("未加载rms：RCS将使用 splines::ns() 回退方案（检验结论等价）")
}
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ========== 2. 工具函数：变量自动识别 =========================================
# 在数据框中按"变量名正则"与"变量标签关键词"双重匹配找列，并打印结果供人工核对
pick_var <- function(df, name_regex = NULL, label_regex = NULL, context = "") {
  nm <- names(df)
  lab <- vapply(df, function(x) { l <- attr(x, "label"); if (is.null(l)) "" else l }, character(1))
  hit <- rep(FALSE, length(nm))
  if (!is.null(name_regex))  hit <- hit | str_detect(tolower(nm),  tolower(name_regex))
  if (!is.null(label_regex)) hit <- hit | str_detect(tolower(lab), tolower(label_regex))
  cand <- nm[hit]
  message(sprintf("[%s] 候选变量：%s", context,
                  if (length(cand)) paste(cand, collapse = ", ") else "（未匹配到，请手工指定）"))
  if (length(cand) >= 1) return(cand[1]) else return(NA_character_)
}

read_dta_safe <- function(path) {
  if (!file.exists(path)) { message("文件不存在：", path); return(NULL) }
  haven::read_dta(path)
}

# 范围过滤：缺失码（999等）与离谱值置NA。血液/体格变量共用
range_filter <- function(x, lo, hi) { x <- as.numeric(x); x[x < lo | x > hi] <- NA; x }

# ========== 3. 基线（2011）血液数据：暴露构建素材 ==============================
w1_blood <- read_dta_safe(FILES$w1_blood)
v_id    <- pick_var(w1_blood, name_regex = "^id$",            context = "ID")
v_glu   <- pick_var(w1_blood, name_regex = "glucose|^glu$|fbg|fasting.*glu", label_regex = "glucose", context = "空腹血糖")
v_hba1c <- pick_var(w1_blood, name_regex = "hba1c|a1c|glyc",  label_regex = "glycated|hba1c", context = "HbA1c")
# 血常规：qc1_vbXXX编码命名，标签已核实（2026-09）→ 直接钉死
# qc1_vb004=血红蛋白, qc1_vb005=红细胞压积, qc1_vb006=MCV, qc1_vb009=血小板, qc1_vb002=白细胞
# 注意：HbA1c标签含"Hemoglobin"会抢血红蛋白的匹配，不可再用标签匹配hb
v_hb    <- pick_var(w1_blood, name_regex = "^qc1_vb004$", context = "血红蛋白")
v_hct   <- pick_var(w1_blood, name_regex = "^qc1_vb005$", context = "红细胞压积")
v_mcv   <- pick_var(w1_blood, name_regex = "^qc1_vb006$", context = "MCV(备用)")
v_plt   <- pick_var(w1_blood, name_regex = "^qc1_vb009$", context = "血小板(备用)")
v_wbc   <- pick_var(w1_blood, name_regex = "^qc1_vb002$", context = "白细胞(备用)")
v_crp   <- pick_var(w1_blood, name_regex = "crp|c-react",     context = "CRP")
v_tg    <- pick_var(w1_blood, name_regex = "triglyceride|^tg$", label_regex = "triglyceride", context = "甘油三酯")
v_tc    <- pick_var(w1_blood, name_regex = "total.*chol|^tc$|tcho", label_regex = "total cholesterol", context = "总胆固醇")
v_hdl   <- pick_var(w1_blood, name_regex = "hdl",             context = "HDL-C")
v_ldl   <- pick_var(w1_blood, name_regex = "ldl",             context = "LDL-C")
v_crea  <- pick_var(w1_blood, name_regex = "creatinine|^crea|^scr", label_regex = "creatinine", context = "肌酐")
v_bun   <- pick_var(w1_blood, name_regex = "bun|urea nitrogen|^urea$", label_regex = "urea nitrogen", context = "BUN")
v_ua    <- pick_var(w1_blood, name_regex = "uric|^ua$",       label_regex = "uric acid", context = "尿酸")

blood_vars <- c(glu = v_glu, hba1c = v_hba1c, hb = v_hb, hct = v_hct,
                mcv = v_mcv, plt = v_plt, wbc = v_wbc, crp = v_crp,
                tg = v_tg, tc = v_tc, hdl = v_hdl, ldl = v_ldl,
                crea = v_crea, bun = v_bun, ua = v_ua)
stopifnot(!is.na(v_glu), !is.na(v_hba1c))   # 核心变量必须识别成功

blood <- w1_blood %>%
  select(ID = all_of(v_id), all_of(na.omit(blood_vars))) %>%
  mutate(across(-ID, as.numeric))

# 范围过滤缺失码（999等）与离谱值；范围参考2011官方用户指南
blood$glu   <- range_filter(blood$glu,   18, 723)
blood$hba1c <- range_filter(blood$hba1c, 3.5, 14.5)
blood$hb    <- range_filter(blood$hb,    5.4, 27.9)
blood$crp   <- range_filter(blood$crp,   0,  180)
blood$tg    <- range_filter(blood$tg,    3,  500)   # 顶编码500
blood$tc    <- range_filter(blood$tc,    50, 400)
blood$hdl   <- range_filter(blood$hdl,   5,  150)
blood$ldl   <- range_filter(blood$ldl,   10, 300)
blood$crea  <- range_filter(blood$crea,  0.1, 15)
blood$bun   <- range_filter(blood$bun,   10, 2000)  # 注意2011年BUN单位为mg/L
blood$ua    <- range_filter(blood$ua,    1,  20)

# --- 血液变量范围核查（2011官方用户指南参考范围；超范围需检查单位）-------------
range_check <- function(x, lo, hi, name) {
  n_out <- sum(x < lo | x > hi, na.rm = TRUE)
  message(sprintf("范围核查 %s: 超出[%g, %g]的观测 %d 个（共%d个非缺失）",
                  name, lo, hi, n_out, sum(!is.na(x))))
}
range_check(blood$glu,   18, 723,  "FPG(mg/dL)")
range_check(blood$hba1c, 3.5, 14.5,"HbA1c(%)")
range_check(blood$hb,    5.4, 27.9,"Hb(g/dL)")
range_check(blood$crp,   0,  180,  "CRP(mg/L)")
range_check(blood$tg,    3,  500,  "TG(mg/dL，注意顶编码500)")
# 注意：2011年BUN单位为 mg/L（不是mg/dL）；如算BUN/Cr比值需先换算

# ========== 4. 基线问卷：人口学与协变量 =======================================
w1_demo <- read_dta_safe(FILES$w1_demo)
v_id_d  <- pick_var(w1_demo, name_regex = "^id$", context = "ID(人口学)")
# 年龄：ba004(标签"Age")2026-09实测仅34人非缺失（仅特殊情况填写）→ 标准做法：2011-出生年(ba002_1)，ba004兜底
v_sex   <- pick_var(w1_demo, name_regex = "sex|gender|ragender", label_regex = "gender|sex", context = "性别")
v_edu   <- pick_var(w1_demo, name_regex = "edu|school", label_regex = "education", context = "教育")
v_marry <- pick_var(w1_demo, name_regex = "marry|marital", label_regex = "marital", context = "婚姻")
v_urban <- pick_var(w1_demo, name_regex = "urban|rural|hukou|resid", label_regex = "urban|rural|residence", context = "城乡(人口学)")

demo <- w1_demo %>%
  select(ID = all_of(v_id_d), sex = all_of(v_sex),
         edu = all_of(v_edu), marital = all_of(v_marry))
demo$age <- 2011 - as.numeric(w1_demo$ba002_1)
demo$age[is.na(demo$age)] <- as.numeric(w1_demo$ba004)[is.na(demo$age)]
if (!is.na(v_urban)) demo$urban <- w1_demo[[v_urban]]

# 城乡变量：demographic文件里的bc002实测仅123个非缺失（2026-09验证，非城乡变量）→
# 触发条件放宽为"缺失率>50%即回退"，从psu.dta（社区层面）经communityID合并
if (is.na(v_urban) || sum(!is.na(demo$urban)) < 0.5 * nrow(demo)) {
  message("人口学文件的城乡变量不可用（缺失率过高），从psu.dta（社区层面）合并……")
  psu <- read_dta_safe(FILES$w1_psu)
  if (!is.null(psu)) {
    v_urban_p <- pick_var(psu, name_regex = "urban|resid", label_regex = "urban|rural", context = "城乡(PSU)")
    v_com_p   <- pick_var(psu, name_regex = "community", context = "社区ID(PSU)")
    v_com_d   <- pick_var(w1_demo, name_regex = "community", context = "社区ID(人口学)")
    if (!any(is.na(c(v_urban_p, v_com_p, v_com_d)))) {
      demo$urban <- NULL   # 清掉错误列，避免合并后变成urban.x/urban.y
      psu_u <- psu %>%
        transmute(.com = as.character(.data[[v_com_p]]),
                  urban = as.numeric(.data[[v_urban_p]])) %>%
        distinct(.com, .keep_all = TRUE)
      demo <- demo %>%
        mutate(.com = as.character(w1_demo[[v_com_d]])) %>%
        left_join(psu_u, by = ".com") %>% select(-.com)
      message("城乡变量已从PSU合并完成")
    } else message("PSU回退失败：请在配置区手工指定城乡变量来源")
  }
}
demo <- demo %>% mutate(across(-ID, as.numeric))

# 人口学变量分布核查（候选有干扰项，必须看分布确认识别正确）
message("年龄分布（预期中位数约59岁，范围45—100+）："); print(summary(demo$age))
message("性别分布（预期1=男、2=女）："); print(table(demo$sex, useNA = "ifany"))
message("教育分布："); print(table(demo$edu, useNA = "ifany"))
message("婚姻分布："); print(table(demo$marital, useNA = "ifany"))
if ("urban" %in% names(demo)) {
  message("城乡分布："); print(table(demo$urban, useNA = "ifany"))
}

w1_health <- read_dta_safe(FILES$w1_health)
v_id_h <- pick_var(w1_health, name_regex = "^id$", context = "ID(健康)")
# 慢病：高血压/糖尿病/血脂异常用da007系列标签匹配（已验证命中）；
# 心脏病/卒中/糖尿病用药按2026-09实测结果钉死变量名（标签句式与本数据版本不符，靠标签匹配会漏）
v_htn    <- pick_var(w1_health, label_regex = "hypertension|high blood pressure", context = "高血压(基线)")
v_dm     <- pick_var(w1_health, label_regex = "diabetes|high blood sugar", context = "糖尿病(基线)")
v_dyslip <- pick_var(w1_health, label_regex = "dyslipidemia|blood lipid", context = "血脂异常(基线)")
v_heart  <- pick_var(w1_health, name_regex = "^da007_7_$", context = "心脏病(基线)")      # 心梗/冠心病/心绞痛/心衰等
v_stroke_w1 <- pick_var(w1_health, name_regex = "^da007_8_$", context = "卒中(基线)")     # 卒中诊断题

# 糖尿病用药：da009_1_3_(中药) 或 da009_2_3_(西药) 任一=1 即视为用药
dm_med_raw1 <- if ("da009_1_3_" %in% names(w1_health)) as.numeric(w1_health$da009_1_3_) else NA
dm_med_raw2 <- if ("da009_2_3_" %in% names(w1_health)) as.numeric(w1_health$da009_2_3_) else NA
v_dm_med <- "dm_med_any"

# da007系列编码：1=是(被诊断)，2=否——先验证编码，再建数据
message("da007_8_(卒中)取值分布："); print(table(w1_health$da007_8_, useNA = "ifany"))
message("da007_7_(心脏病)取值分布："); print(table(w1_health$da007_7_, useNA = "ifany"))

health1 <- w1_health %>%
  select(ID = all_of(v_id_h), htn = all_of(v_htn), dm = all_of(v_dm),
         dyslip = all_of(v_dyslip), heart = all_of(v_heart),
         stroke_w1 = all_of(v_stroke_w1)) %>%
  mutate(across(-ID, as.numeric),
         dm_med = as.integer(dm_med_raw1 == 1 | dm_med_raw2 == 1),
         dm_med = if_else(is.na(dm_med_raw1) & is.na(dm_med_raw2), NA_integer_, dm_med))

w1_phys <- read_dta_safe(FILES$w1_phys)
v_id_p <- pick_var(w1_phys, name_regex = "^id$", context = "ID(体格)")
# 血压/脉搏：自动抓取全部测量列（2011年通常3次：systolic1—3、diastolic1—3），取行均值
lab_p <- vapply(w1_phys, function(x) { l <- attr(x, "label"); if (is.null(l)) "" else l }, character(1))
sbp_cols   <- names(w1_phys)[str_detect(tolower(names(w1_phys)), "systolic|^sbp") |
                               str_detect(tolower(lab_p), "systolic")]
dbp_cols   <- names(w1_phys)[str_detect(tolower(names(w1_phys)), "diastolic|^dbp") |
                               str_detect(tolower(lab_p), "diastolic")]
pulse_cols <- names(w1_phys)[str_detect(tolower(names(w1_phys)), "pulse") |
                               str_detect(tolower(lab_p), "pulse")]
message("SBP测量列: ", paste(sbp_cols, collapse = ", "))
message("DBP测量列: ", paste(dbp_cols, collapse = ", "))
message("脉搏测量列: ", paste(pulse_cols, collapse = ", "), "（提取备用，可用于RHR类选题）")
phys <- data.frame(ID = w1_phys[[v_id_p]])
# 血压/脉搏3次测量先各自过滤缺失码（999等），再取行均值（range_filter定义见第2节）
sbp_mat   <- as.data.frame(lapply(w1_phys[sbp_cols],   range_filter, lo = 60, hi = 260))
dbp_mat   <- as.data.frame(lapply(w1_phys[dbp_cols],   range_filter, lo = 30, hi = 160))
pulse_mat <- as.data.frame(lapply(w1_phys[pulse_cols], range_filter, lo = 30, hi = 200))
phys$sbp   <- rowMeans(sbp_mat, na.rm = TRUE)
phys$dbp   <- rowMeans(dbp_mat, na.rm = TRUE)
phys$pulse <- rowMeans(pulse_mat, na.rm = TRUE)
phys[c("sbp","dbp","pulse")] <- lapply(phys[c("sbp","dbp","pulse")],
                                       function(x) { x[is.nan(x)] <- NA; x })
phys$height <- range_filter(w1_phys$qi002, 120, 210)
phys$weight <- range_filter(w1_phys$ql002, 30, 200)
phys$waist  <- range_filter(w1_phys$qm002, 40, 200)
phys$bmi    <- phys$weight / (phys$height / 100)^2
message("BMI分布核查："); print(summary(phys$bmi))   # 预期中位数约23—24，极端值<15或>40应极少

# ========== 5. 结局提取：各随访波次卒中 + 退出访谈死亡 =========================
# 卒中变量结构（2026-09逐波核实）：
#   w2/w3(2013/2015)：da019_w2_1=追访新诊断("Doctor Tell Stroke") + da007_8_=曾被诊断(新样本)
#   w4(2018)：da007_8_（该波标签规范，曾被诊断题）
#   w5(2020)：慢病结构大改——da002_X_=旧病转归题（非诊断，不可用），da003_X_=新诊断题，
#             卒中=da003_8_（2026-09交叉验证通过：2018无卒中者中363例新诊断，跳答逻辑自洽）
#   阳性=任一题=1；ID出现在该波文件中=该波已访，未报告阳性记0
#   注意：da019_w2_1的NA为跳答（无新诊断才跳过），按"未报告"处理——正式分析前请按codebook复核跳答逻辑
get_wave_stroke <- function(path, wave_label) {
  d <- read_dta_safe(path); if (is.null(d)) return(NULL)
  idv <- pick_var(d, name_regex = "^id$", context = paste0("ID(", wave_label, ")"))
  use_new  <- "da019_w2_1" %in% names(d)   # w2/w3新诊断题
  use_ever <- "da007_8_"  %in% names(d)    # w3/w4曾被诊断题（w4标签规范）
  use_2020 <- "da003_8_"  %in% names(d)    # w5新诊断题（2020结构）
  if (!use_new && !use_ever && !use_2020) {
    message("  ", wave_label, "波未找到卒中变量，该波结局缺失"); return(NULL)
  }
  pos <- rep(FALSE, nrow(d))
  if (use_new)  { v <- as.numeric(d$da019_w2_1); pos <- pos | (!is.na(v) & v == 1) }
  if (use_ever) { v <- as.numeric(d$da007_8_);  pos <- pos | (!is.na(v) & v == 1) }
  if (use_2020) { v <- as.numeric(d$da003_8_);  pos <- pos | (!is.na(v) & v == 1) }
  message("  ", wave_label, "波卒中阳性 ", sum(pos), " 例（该波记录数 ", nrow(d), "）")
  out <- data.frame(ID = d[[idv]])
  out[[paste0("stroke_", wave_label)]] <- as.integer(pos)
  out
}
stroke_waves <- list(
  get_wave_stroke(FILES$w2_health, "w2"),
  get_wave_stroke(FILES$w3_health, "w3"),
  get_wave_stroke(FILES$w4_health, "w4"),
  get_wave_stroke(FILES$w5_health, "w5")
) %>% purrr::compact() %>% purrr::reduce(full_join, by = "ID")

# 死亡：退出访谈仅对去世者进行——出现在Exit文件中即代表死亡（无需判断取值）；
# 顺带提取死亡年份(exb001_1)备用。死亡日期精度：2013波有确切年月，2020波见codebook
get_wave_death <- function(path, wave_label) {
  d <- read_dta_safe(path); if (is.null(d)) return(NULL)
  idv <- pick_var(d, name_regex = "^id$", context = paste0("ID(退出", wave_label, ")"))
  message("  ", wave_label, "退出访谈记录数(≈该时段死亡数): ", nrow(d))
  out <- d %>% transmute(ID = .data[[idv]], !!paste0("death_", wave_label) := 1L)
  if ("exb001_1" %in% names(d)) out[[paste0("deathyr_", wave_label)]] <- as.numeric(d$exb001_1)
  distinct(out, ID, .keep_all = TRUE)
}
death_waves <- list(
  get_wave_death(FILES$w2_exit, "w2"), get_wave_death(FILES$w5_exit, "w5")
) %>% purrr::compact() %>% purrr::reduce(full_join, by = "ID")

# ---- ID格式统一（2026-09验证，匹配率93.5%）------------------------------------
# 2011基线为11位ID；2013+为12位（=前9位+"0"+末2位）。反向转换：12位删去第10位的"0"
to_w1_id <- function(id) {
  id <- as.character(id)
  hit <- nchar(id) == 12 & substr(id, 10, 10) == "0"
  id[hit] <- paste0(substr(id[hit], 1, 9), substr(id[hit], 11, 12))
  id
}
stroke_waves$ID <- to_w1_id(stroke_waves$ID)
death_waves$ID  <- to_w1_id(death_waves$ID)
# 转换后若有重复ID（罕见：多个12位映射同一11位），保阳性去重
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
stroke_waves <- stroke_waves %>% group_by(ID) %>%
  summarise(across(everything(), safe_max), .groups = "drop")
death_waves  <- death_waves  %>% group_by(ID) %>%
  summarise(across(everything(), safe_max), .groups = "drop")
message("ID转换后，基线与随访可匹配人数: ",
        sum(stroke_waves$ID %in% demo$ID), " / ", nrow(stroke_waves))

# ========== 6. 合并并构建分析变量 =============================================
analytic <- demo %>%
  full_join(blood,   by = "ID") %>%
  full_join(health1, by = "ID") %>%
  full_join(phys,    by = "ID") %>%
  full_join(stroke_waves, by = "ID") %>%
  { if (!is.null(death_waves)) left_join(., death_waves, by = "ID") else . }

# 结局1：随访期新发卒中（任一后续波次报告=1；全NA=失访）
stroke_cols <- grep("^stroke_w[2-5]$", names(analytic), value = TRUE)
analytic <- analytic %>%
  mutate(
    n_wave_obs   = rowSums(!is.na(across(all_of(stroke_cols)))),
    incident_stroke = case_when(
      n_wave_obs == 0 ~ NA_real_,                       # 无任何随访→失访
      rowSums(across(all_of(stroke_cols)), na.rm = TRUE) > 0 ~ 1,
      TRUE ~ 0),
    first_stroke_wave = apply(across(all_of(stroke_cols)), 1, function(r) {
      hit <- which(r == 1); if (length(hit)) hit[1] + 1 else NA_real_ })  # 2/3/4/5波
  )
# 结局2：随访期死亡（竞争事件）。在任一Exit文件出现=死亡=1；
# 有随访记录且未出现在Exit文件=存活=0；完全失访=NA
death_cols <- grep("^death_w[2-5]$", names(analytic), value = TRUE)
if (length(death_cols)) {
  analytic$death <- as.integer(rowSums(analytic[, death_cols, drop = FALSE] == 1, na.rm = TRUE) > 0)
  analytic$death[analytic$death == 0 & analytic$n_wave_obs == 0] <- NA
  message("死亡人数: ", sum(analytic$death == 1, na.rm = TRUE),
          "；存活: ", sum(analytic$death == 0, na.rm = TRUE),
          "；状态未知(失访): ", sum(is.na(analytic$death)))
}

# ========== 7. 纳排流程（STROBE数字，逐步打印） ===============================
flow <- tibble(step = character(), n = integer(), excluded = integer())
add_step <- function(df_before, df_after, label) {
  flow <<- bind_rows(flow, tibble(
    step = label, n = nrow(df_after), excluded = nrow(df_before) - nrow(df_after)))
  df_after
}
d0 <- analytic
message("Wave1合并总样本: ", nrow(d0))
d1 <- d0 %>% filter(!is.na(glu) & !is.na(hba1c)) %>% add_step(d0, ., "有2011年FPG与HbA1c")
d2 <- d1 %>% filter(age >= MIN_AGE)              %>% add_step(d1, ., "年龄≥45岁")
d3 <- d2 %>% filter(is.na(stroke_w1) | stroke_w1 != 1) %>% add_step(d2, ., "排除基线卒中")
d4 <- d3 %>% filter(!is.na(incident_stroke))     %>% add_step(d3, ., "至少有1次随访结局")
d5 <- d4 %>% filter(!is.na(sex))                 %>% add_step(d4, ., "关键协变量(性别)非缺失")
print(flow); message("最终分析样本: ", nrow(d5), "；事件数: ", sum(d5$incident_stroke))

# ========== 8. HGI构建 =========================================================
# 慢病变量编码统一：da007系列1=是/2=否 → 1/0，便于模型系数解释（不改变任何判定逻辑）
d5 <- d5 %>% mutate(across(c(htn, dm, dyslip, heart),
                           ~ if_else(is.na(.x), NA_real_, as.numeric(.x == 1))))
d5 <- d5 %>% mutate(
  diabetes = case_when(dm == 1 | (!is.na(dm_med) & dm_med == 1) |
                         hba1c >= 6.5 | glu >= 126 ~ 1, TRUE ~ 0)
)

# 主方程：全分析样本内 HbA1c ~ FPG
fit_main <- lm(hba1c ~ glu, data = d5)
d5$hgi_pred <- predict(fit_main)
d5$hgi      <- d5$hba1c - d5$hgi_pred
message("主方程: HbA1c = ", round(coef(fit_main)[1], 4), " + ",
        round(coef(fit_main)[2], 4), " × FPG；R² = ", round(summary(fit_main)$r.squared, 3))

# 敏感性方程A：仅在非糖尿病人群建方程，外推全人群
fit_ndm <- lm(hba1c ~ glu, data = filter(d5, diabetes == 0))
d5$hgi_ndm <- d5$hba1c - predict(fit_ndm, newdata = d5)

# 敏感性方程B：已发表CHARLS方程（肝病研究）：pred = 0.013 × FPG + 4.804
d5$hgi_pub <- d5$hba1c - (0.013 * d5$glu + 4.804)

# 分层HGI三分位 + 贫血标记（WHO：男<13.0、女<12.0 g/dL）
d5 <- d5 %>% mutate(
  hgi_tertile = ntile(hgi, 3),
  anemia = case_when(sex == 1 & hb < 13.0 ~ 1,   # 按性别编码调整：CHARLS通常1=男
                     sex == 2 & hb < 12.0 ~ 1,
                     !is.na(hb) ~ 0, TRUE ~ NA_real_))

# ========== 8b. ΔHGI子分析：2015血样 → HGI变化（可选，主分析跑通后再做） =========
# 关键：2015年血样由金域(KingMed)检测、2011年由右安门/首医大检测，存在批次效应，
# 因此2015年HGI必须在2015年样本内部重新建 HbA1c~FPG 方程，不能沿用2011年方程系数
w3_blood <- read_dta_safe(FILES$w3_blood)
if (!is.null(w3_blood)) {
  v_id3   <- pick_var(w3_blood, name_regex = "^id$", context = "ID(2015血)")
  v_glu3  <- pick_var(w3_blood, name_regex = "glucose|^glu$|fbg|fasting.*glu", label_regex = "glucose", context = "FPG(2015)")
  v_a1c3  <- pick_var(w3_blood, name_regex = "hba1c|a1c|glyc", label_regex = "glycated|hba1c", context = "HbA1c(2015)")
  stopifnot(!is.na(v_glu3), !is.na(v_a1c3))
  blood3 <- w3_blood %>%
    select(ID = all_of(v_id3), glu3 = all_of(v_glu3), hba1c3 = all_of(v_a1c3)) %>%
    mutate(across(-ID, as.numeric), ID = to_w1_id(ID))   # 2015血样也是12位ID，需转换
  range_check(blood3$glu3, 18, 723, "FPG2015(mg/dL)"); range_check(blood3$hba1c3, 3.5, 14.5, "HbA1c2015(%)")

  fit_2015 <- lm(hba1c3 ~ glu3, data = blood3)     # 2015内部方程（抵消批次效应）
  # 不用predict()：它会丢缺失行导致错位——用系数手工计算，长度天然对齐
  blood3$hgi_2015 <- blood3$hba1c3 - (coef(fit_2015)[1] + coef(fit_2015)[2] * blood3$glu3)
  message("2015方程: HbA1c = ", round(coef(fit_2015)[1], 4), " + ",
          round(coef(fit_2015)[2], 4), " × FPG；R² = ", round(summary(fit_2015)$r.squared, 3))

  ddelta <- d5 %>% inner_join(blood3 %>% select(ID, hgi_2015), by = "ID") %>%
    mutate(delta_hgi = hgi_2015 - hgi)
  message("ΔHGI子样本（两期血样齐全且在分析集内）: ", nrow(ddelta), "人")
  write.csv(ddelta, file.path(OUT_DIR, "analytic_delta_hgi.csv"), row.names = FALSE)
}

# ========== 9. 快速筛查（立项决策用，1—2天工作量） =============================
# 9.1 随访时间：区间删失简化处理=基线至首次报告波次的间隔年数（快速筛查可接受；
#     正式分析用中点插补或区间删失模型，见操作手册）
wave_year <- c(w2 = 2013, w3 = 2015, w4 = 2018, w5 = 2020)
d5 <- d5 %>% mutate(
  fu_years = if_else(incident_stroke == 1,
                     wave_year[paste0("w", first_stroke_wave)] - 2011,
                     2020 - 2011))

# 9.2 Cox三模型
f <- Surv(fu_years, incident_stroke) ~ hgi
m1 <- coxph(f, data = d5)
m2 <- coxph(update(f, . ~ . + age + sex + urban + edu), data = d5)
m3 <- coxph(update(f, . ~ . + age + sex + urban + edu + htn + dm + dyslip + heart +
                     sbp + crp + tg + hdl), data = d5)
for (m in list(m1, m2, m3)) print(summary(m)$conf.int["hgi", , drop = FALSE])

# 9.3 三分位趋势 + RCS形态
d5$hgi_tertile <- factor(d5$hgi_tertile)
print(coxph(update(f, . ~ hgi_tertile + age + sex + urban + edu + htn + dm +
                     dyslip + heart), data = d5))

if (has_rms) {
  library(rms)
  dd <- datadist(d5); options(datadist = "dd")
  m_rcs <- cph(Surv(fu_years, incident_stroke) ~ rcs(hgi, 4) + age + sex,
               data = d5, x = TRUE, y = TRUE)
  print(anova(m_rcs))   # 看非线性P值：判断U/L/J型的直接证据
} else {
  # --- 回退方案：splines::ns() 实现RCS，检验结论与rms等价 ---
  m_lin <- coxph(Surv(fu_years, incident_stroke) ~ hgi + age + sex, data = d5)
  m_rcs <- coxph(Surv(fu_years, incident_stroke) ~ ns(hgi, df = 4) + age + sex, data = d5)
  lrt <- anova(m_lin, m_rcs)          # 似然比检验：非线性是否成立
  print(lrt)
  message("非线性检验 P = ", signif(lrt[["Pr(>|Chi|)"]][2], 3),
          "（<0.05提示存在非线性剂量反应，结合图形判断U/L/J型）")

  # 画剂量反应曲线（参照点=HGI中位数），存为rcs_hgi.png
  b    <- ns(d5$hgi, df = 4)
  grid <- seq(quantile(d5$hgi, .01, na.rm = TRUE),
              quantile(d5$hgi, .99, na.rm = TRUE), length.out = 100)
  Xn   <- predict(b, grid)
  Xr   <- predict(b, median(d5$hgi, na.rm = TRUE))
  Xc   <- sweep(Xn, 2, as.numeric(Xr), "-")
  co   <- coef(m_rcs); V <- vcov(m_rcs)
  idx  <- grep("ns\\(hgi", names(co))
  eta  <- as.numeric(Xc %*% co[idx])
  se   <- sqrt(pmax(diag(Xc %*% V[idx, idx, drop = FALSE] %*% t(Xc)), 0))
  png(file.path(OUT_DIR, "rcs_hgi.png"), width = 1600, height = 1200, res = 200)
  plot(grid, exp(eta), type = "l", lwd = 2,
       ylim = range(exp(pmax(eta - 1.96*se, -5)), exp(pmin(eta + 1.96*se, 5))),
       xlab = "HGI", ylab = "HR (95% CI), ref = median",
       main = "HGI与新发卒中的剂量反应（RCS, ns回退方案）")
  lines(grid, exp(pmax(eta - 1.96*se, -5)), lty = 2)
  lines(grid, exp(pmin(eta + 1.96*se, 5)), lty = 2)
  abline(h = 1, col = "grey50")
  dev.off()
}

# 9.4 描述性基线表（按HGI三分位）
tab1_vars <- c("age","sex","urban","edu","htn","dm","dyslip","heart",
               "glu","hba1c","hgi","hb","crp","tg","hdl","sbp","anemia")
print(CreateTableOne(vars = tab1_vars, strata = "hgi_tertile",
                     data = as.data.frame(d5), factorVars = c("sex","urban","edu","htn","dm","dyslip","heart","anemia")))

# 9.5 关键稳健性筛查（立项决策的另外三条标准 + 审稿人必问）--------------------
# ① 审稿人必问：调整HbA1c本身后，HGI是否仍独立？（对应Meta分析争议点）
m3_a1c <- coxph(update(f, . ~ . + age + sex + urban + edu + htn + dm + dyslip +
                         heart + hba1c), data = d5)
cat("【标准⑤·审稿人必问】调整HbA1c后 HGI的HR：\n")
print(summary(m3_a1c)$conf.int["hgi", , drop = FALSE])

# ② 标准③：非糖尿病亚组（diabetes==0）
d5_ndm <- filter(d5, diabetes == 0)
m_ndm <- coxph(update(f, . ~ . + age + sex + urban + edu + htn), data = d5_ndm)
cat("【标准③】非糖尿病亚组（n =", nrow(d5_ndm), "，事件", sum(d5_ndm$incident_stroke), "）HGI的HR：\n")
print(summary(m_ndm)$conf.int["hgi", , drop = FALSE])

# ③ 标准④：排除贫血（HbA1c失真人群）
d5_noanemia <- filter(d5, is.na(anemia) | anemia == 0)
m_noanemia <- coxph(update(f, . ~ . + age + sex + urban + edu + htn + dm + dyslip + heart),
                    data = d5_noanemia)
cat("【标准④】排除贫血后（n =", nrow(d5_noanemia), "，事件", sum(d5_noanemia$incident_stroke), "）HGI的HR：\n")
print(summary(m_noanemia)$conf.int["hgi", , drop = FALSE])

# 9.6 J型故事的最终检验（线性故事已被HbA1c调整否决，非线性故事需自己的证据链）----
# ① 充分调整后的非线性是否仍在（此前LRT只调整了年龄性别）
m3_rcs <- coxph(Surv(fu_years, incident_stroke) ~ ns(hgi, df = 4) + age + sex + urban +
                  edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl, data = d5)
cat("【检验A】充分调整后非线性LRT：\n"); print(anova(m3, m3_rcs))

# ② J型右臂二分化：高HGI组 vs 其余（两个切点：上10%分位 + 绝对值0.5）
for (cut in c(quantile(d5$hgi, 0.9, na.rm = TRUE), 0.5)) {
  d5$hgi_hi <- as.integer(d5$hgi >= cut)
  m_hi <- coxph(update(f, . ~ hgi_hi + age + sex + urban + edu + htn + dm + dyslip + heart),
                data = d5)
  cat("【检验B】高HGI(≥", round(cut, 3), ")组HR（占", mean(d5$hgi_hi) * 100, "%）：\n", sep = "")
  print(summary(m_hi)$conf.int["hgi_hi", , drop = FALSE])
}

# ③ 非糖尿病亚组的非线性（J型是否在无糖尿病人群仍成立）
m_ndm_lin <- coxph(Surv(fu_years, incident_stroke) ~ hgi + age + sex, data = d5_ndm)
m_ndm_rcs <- coxph(Surv(fu_years, incident_stroke) ~ ns(hgi, df = 4) + age + sex, data = d5_ndm)
cat("【检验C】非糖尿病亚组非线性LRT：\n"); print(anova(m_ndm_lin, m_ndm_rcs))

# ========== 10. 导出 ===========================================================
write.csv(d5, file.path(OUT_DIR, "analytic_hgi_baseline.csv"), row.names = FALSE)
write.csv(flow, file.path(OUT_DIR, "strobe_flow.csv"), row.names = FALSE)
saveRDS(list(fit_main = fit_main, fit_ndm = fit_ndm, m1 = m1, m2 = m2, m3 = m3,
             m_rcs = m_rcs), file.path(OUT_DIR, "screen_models.rds"))
message("完成。请核对：①所有pick_var打印的变量识别结果 ②flow表各步排除人数 ",
        "③主方程R²(预期0.3—0.5) ④Cox与RCS结果 → 决定是否立项。")
