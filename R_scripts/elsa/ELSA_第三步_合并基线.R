####################################################################################
# ELSA external validation -- Step 3: pooled-baseline sensitivity analysis
# New entrants with fasting blood at Wave 2, 4, or 6 each contribute one baseline
# (first eligible wave). HGI equation + P90 fitted per entry wave.
# Output: build_log_pooled.txt, analytic_elsa_hgi_pooled.csv, ELSA_pooled_results.txt
####################################################################################
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(haven); library(dplyr); library(survival) })
D <- "D:/CHARLS/HGI_stroke/ELSA验证/data_raw/stata/stata13_se/"
O <- "D:/CHARLS/HGI_stroke/ELSA验证/"
logf <- file(paste0(O, "build_log_pooled.txt"), open = "wt"); sink(logf, split = TRUE)
na_neg <- function(x) { x[x < 0] <- NA; x }

## ---- 1. harmonized core (all waves we need) --------------------------------------------
gh <- read_dta(paste0(D, "gh_elsa_h.dta"),
               col_select = c("idauniq", "ragender", "raeducl",
                              "r2agey", "r4agey", "r6agey",
                              "r2hibpe", "r4hibpe", "r6hibpe",
                              "r2diabe", "r4diabe", "r6diabe",
                              "r2hchole", "r4hchole", "r6hchole",
                              "r2hearte", "r4hearte", "r6hearte",
                              "r2iwy", "r4iwy", "r6iwy",
                              paste0("r", 1:10, "stroke"),
                              paste0("r", 1:10, "iwstat"),
                              paste0("r", 1:10, "iwy"))) %>%
  mutate(across(-idauniq, as.numeric), across(-idauniq, na_neg))
cat("harmonized:", nrow(gh), "\n")

## ---- 2. nurse biomarkers per wave --------------------------------------------------------
get_nurse <- function(file, w) {
  nu <- read_dta(paste0(D, file),
                 col_select = tidyselect::any_of(c("idauniq", "fglu", "hba1c", "chol", "hdl",
                                                   "trig", "hscrp", "sysval", "SYSVAL"))) %>%
    mutate(across(-idauniq, as.numeric), across(-idauniq, na_neg))
  if ("SYSVAL" %in% names(nu)) nu <- nu %>% rename(sysval = SYSVAL)
  nu$bw <- w
  # W6 reports HbA1c in IFCC mmol/mol; convert to NGSP % (NGSP = 0.09148*IFCC + 2.152)
  if (median(nu$hba1c, na.rm = TRUE) > 20) {
    nu$hba1c <- 0.09148 * nu$hba1c + 2.152
    cat("note: wave", w, "HbA1c converted IFCC -> NGSP %\n")
  }
  nu %>% select(idauniq, bw, fglu, hba1c, tidyselect::any_of(c("chol","hdl","trig","hscrp","sysval")))
}
nu_all <- bind_rows(get_nurse("wave_2_nurse_data_v2.dta", 2),
                    get_nurse("wave_4_nurse_data.dta", 4),
                    get_nurse("wave_6_elsa_nurse_data_v2.dta", 6))
nu_all <- nu_all %>% filter(!is.na(fglu) & !is.na(hba1c)) %>%
  group_by(idauniq) %>% arrange(bw) %>% slice(1) %>% ungroup()   # first eligible wave
cat("persons with fglu+hba1c at any of W2/4/6 (first wave kept):", nrow(nu_all), "\n")
print(table(nu_all$bw))

## ---- 3. W11 stroke, urban, death ---------------------------------------------------------
w11 <- read_dta(paste0(D, "wave_11_elsa_data_eul.dta"),
                col_select = c("idauniq", "hediia7", "w11indout", "iintdaty")) %>%
  mutate(across(-idauniq, as.numeric), across(-idauniq, na_neg))
geo <- read_dta(paste0(D, "elsa_geog_urindewr_2001_eul.dta"),
                col_select = c("idauniq", "w2_urindewr_2001", "w4_urindewr_2001", "w6_urindewr_2001")) %>%
  mutate(across(-idauniq, as.numeric), across(-idauniq, na_neg))
eol <- read_dta(paste0(D, "h_elsa_eol_a2.dta"),
                col_select = c("idauniq", "raxyear")) %>%
  mutate(across(-idauniq, as.numeric), across(-idauniq, na_neg))

df <- nu_all %>% left_join(gh, by = "idauniq") %>% left_join(w11, by = "idauniq") %>%
  left_join(geo, by = "idauniq") %>% left_join(eol, by = "idauniq")
cat("merged:", nrow(df), "\n")

## ---- 4. per-person baseline extraction ----------------------------------------------------
g <- function(col, bw) ifelse(bw == 2, df[[paste0("r2", col)]],
                       ifelse(bw == 4, df[[paste0("r4", col)]], df[[paste0("r6", col)]]))
df <- df %>% mutate(
  age   = g("agey", bw), htn = g("hibpe", bw), dm = g("diabe", bw),
  dyslip = g("hchole", bw), heart = g("hearte", bw), base_iwy = g("iwy", bw),
  sex = ragender - 1, edu = raeducl,
  urban = as.integer(ifelse(bw == 2, w2_urindewr_2001,
                     ifelse(bw == 4, w4_urindewr_2001, w6_urindewr_2001)) == 1),
  sbp = sysval, crp = hscrp, tg = trig * 88.57, hdl = hdl * 38.67,
  death_year = raxyear)

## ---- 5. flow -----------------------------------------------------------------------------
cat("\n== STROBE flow (pooled) ==\n")
df <- df %>% filter(age >= 50); cat("1. age>=50 at entry:", nrow(df), "\n")
# stroke-free at entry: no stroke reported in any wave <= bw
stroke_free <- rep(TRUE, nrow(df))
for (j in 1:6) {
  s <- df[[paste0("r", j, "stroke")]]
  stroke_free <- stroke_free & !(j <= df$bw & !is.na(s) & s == 1)
}
df <- df[stroke_free, ]; cat("2. stroke-free at entry:", nrow(df), "\n")

# last contact after entry (waves bw+1 .. 10 interviewed, or W11)
df$w11_iwy <- df$iintdaty
df$w11_int <- as.integer(!is.na(df$w11indout) & df$w11indout >= 11 & df$w11indout <= 25)
lastc <- rep(NA_real_, nrow(df))
for (i in seq_len(nrow(df))) {
  yrs <- c()
  for (j in 3:10) if (j > df$bw[i] && !is.na(df[[paste0("r", j, "iwstat")]][i]) &&
                      df[[paste0("r", j, "iwstat")]][i] == 1)
    yrs <- c(yrs, df[[paste0("r", j, "iwy")]][i])
  if (df$w11_int[i] == 1 && df$bw[i] < 11) yrs <- c(yrs, df$w11_iwy[i])
  if (length(yrs)) lastc[i] <- max(yrs)
}
df$last_contact <- lastc
df <- df %>% filter(!is.na(last_contact) & last_contact > base_iwy)
cat("3. at least one interview after entry:", nrow(df), "\n")

## ---- 6. incident stroke -------------------------------------------------------------------
ev_wave <- rep(NA_integer_, nrow(df))
for (i in seq_len(nrow(df))) {
  for (j in max(3, df$bw[i] + 1):10) {
    s <- df[[paste0("r", j, "stroke")]][i]
    if (!is.na(s) && s == 1) { ev_wave[i] <- j; break }
  }
  if (is.na(ev_wave[i]) && !is.na(df$hediia7[i]) && df$hediia7[i] == 1) ev_wave[i] <- 11L
}
df$ev_wave <- ev_wave
df$incident_stroke <- as.integer(!is.na(ev_wave))

rep_y <- prev_y <- rep(NA_real_, nrow(df))
for (i in seq_len(nrow(df))) {
  w <- df$ev_wave[i]; if (is.na(w)) next
  rep_y[i] <- if (w <= 10) df[[paste0("r", w, "iwy")]][i] else df$w11_iwy[i]
  pw <- w - 1
  while (pw > df$bw[i]) {
    s <- df[[paste0("r", pw, "stroke")]][i]
    if (!is.na(s) && s == 0) { prev_y[i] <- df[[paste0("r", pw, "iwy")]][i]; break }
    pw <- pw - 1
  }
  if (is.na(prev_y[i])) prev_y[i] <- df$base_iwy[i]
}
df$event_year_rep <- rep_y; df$event_year_prev <- prev_y
df$fu_mid <- ifelse(df$incident_stroke == 1, (rep_y + prev_y) / 2 - df$base_iwy,
                    df$last_contact - df$base_iwy)
df <- df %>% filter(is.na(fu_mid) | fu_mid > 0)
cat("4. analysis-eligible:", nrow(df), "; events:", sum(df$incident_stroke), "\n")
print(table(df$bw, df$incident_stroke))

## ---- 7. HGI per entry wave ------------------------------------------------------------------
df$hgi <- NA_real_; df$hgi_hi_p90 <- NA_integer_
for (w in c(2, 4, 6)) {
  idx <- which(df$bw == w)
  eq <- lm(hba1c ~ fglu, data = df[idx, ])
  cat("\nW", w, " HGI eq: intercept", round(coef(eq)[1], 3), " slope", round(coef(eq)[2], 3),
      " n =", length(idx), "\n", sep = "")
  df$hgi[idx] <- df$hba1c[idx] - predict(eq, newdata = df[idx, ])
  cut <- as.numeric(quantile(df$hgi[idx], 0.9, na.rm = TRUE))
  df$hgi_hi_p90[idx] <- as.integer(df$hgi[idx] >= cut)
  cat("W", w, " P90 =", round(cut, 4), " high n =", sum(df$hgi_hi_p90[idx]), "\n", sep = "")
}
df$shr <- df$fglu / (1.59 * df$hba1c - 2.59)
df$tyg <- log((df$trig * 88.57) * (df$fglu * 18) / 2)
df$died <- as.integer(apply(df[, paste0("r", 1:10, "iwstat")], 1, function(r) any(r %in% c(5, 6)))) |
  as.integer(!is.na(df$w11indout) & df$w11indout == 95)

covs <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
main <- df %>% filter(complete.cases(across(all_of(c(covs, "hgi", "hba1c",
                                                     "incident_stroke", "fu_mid", "bw")))))
cat("\n== pooled main complete-case set ==\n")
cat("n =", nrow(main), "; events =", sum(main$incident_stroke),
    "; deaths =", sum(main$died, na.rm = TRUE),
    "; high HGI =", sum(main$hgi_hi_p90), "\n")
print(table(main$bw, main$incident_stroke))
write.csv(df, paste0(O, "analytic_elsa_hgi_pooled_full.csv"), row.names = FALSE)
write.csv(main, paste0(O, "analytic_elsa_hgi_pooled.csv"), row.names = FALSE)

## ---- 8. analysis ----------------------------------------------------------------------------
resf <- file(paste0(O, "ELSA_pooled_results.txt"), open = "wt"); sink(resf, split = TRUE)
cat("ELSA pooled-baseline (W2/W4/W6) sensitivity analysis\n")
cat("n =", nrow(main), " events =", sum(main$incident_stroke),
    " deaths =", sum(main$died), "\n\n")

m1 <- coxph(Surv(fu_mid, incident_stroke) ~ hgi_hi_p90 + strata(bw), data = main)
m2 <- coxph(Surv(fu_mid, incident_stroke) ~ hgi_hi_p90 + age + sex + urban + edu + strata(bw), data = main)
f3 <- as.formula(paste("Surv(fu_mid, incident_stroke) ~ hgi_hi_p90 +",
                       paste(covs, collapse = "+"), "+ strata(bw)"))
m3 <- coxph(f3, data = main)
pr <- function(m, tag) {
  s <- summary(m)$conf.int["hgi_hi_p90", ]; p <- summary(m)$coefficients["hgi_hi_p90", "Pr(>|z|)"]
  cat(sprintf("%-10s HR %.3f (%.3f-%.3f) P = %.4f\n", tag, s[1], s[3], s[4], p))
}
cat("== Cox (stratified by entry wave) ==\n")
pr(m1, "M1_unadj"); pr(m2, "M2_demo"); pr(m3, "M3_full")
cat("PH test (m3):\n"); print(cox.zph(m3)$table["hgi_hi_p90", ])

f3b <- as.formula(paste("Surv(fu_mid, incident_stroke) ~ hgi_hi_p90 + hba1c +",
                        paste(covs, collapse = "+"), "+ strata(bw)"))
m3b <- coxph(f3b, data = main)
cat("\n== key test (+ HbA1c) ==\n"); print(summary(m3b)$coefficients[c("hgi_hi_p90","hba1c"), ])
print(summary(m3b)$conf.int[c("hgi_hi_p90","hba1c"), ])

# per-SD continuous
m3c <- coxph(update(f3, ~ . - hgi_hi_p90 + I(scale(hgi))), data = main)
cat("\nper-SD HGI:\n"); print(summary(m3c)$conf.int["I(scale(hgi))", , drop = FALSE])

# wave-specific M2/M3 + descriptives
cat("\n== per-wave descriptives & estimates ==\n")
for (w in c(2, 4, 6)) {
  sub <- main %>% filter(bw == w)
  cat(sprintf("\nW%d: n=%d, events=%d, deaths=%d, median age=%.0f, median fu=%.1f y\n",
              w, nrow(sub), sum(sub$incident_stroke), sum(sub$died, na.rm=TRUE),
              median(sub$age), median(sub$fu_mid)))
  m2w <- tryCatch(coxph(Surv(fu_mid, incident_stroke) ~ hgi_hi_p90 + age + sex + urban + edu,
                        data = sub), error = function(e) NULL)
  if (!is.null(m2w)) pr(m2w, paste0("W", w, "_M2"))
  m3w <- tryCatch(coxph(update(f3, ~ . - strata(bw)), data = sub), error = function(e) NULL)
  if (!is.null(m3w)) pr(m3w, paste0("W", w, "_M3"))
}

# W2+W4 only variant (both NGSP %-era cohorts)
sub24 <- main %>% filter(bw %in% c(2, 4))
f24 <- as.formula(paste("Surv(fu_mid, incident_stroke) ~ hgi_hi_p90 +",
                        paste(covs, collapse = "+"), "+ strata(bw)"))
cat("\n== W2+W4 only (n=", nrow(sub24), ", events=", sum(sub24$incident_stroke), ") ==\n", sep = "")
pr(coxph(Surv(fu_mid, incident_stroke) ~ hgi_hi_p90 + strata(bw), data = sub24), "M1")
pr(coxph(Surv(fu_mid, incident_stroke) ~ hgi_hi_p90 + age + sex + urban + edu + strata(bw),
         data = sub24), "M2")
pr(coxph(f24, data = sub24), "M3")

# Fine-Gray: multi-state status 0=censored, 1=stroke, 2=death (death time = raxyear,
# else wave year of first iwstat==5, else midpoint around iwstat==6)
fgd <- main %>%
  mutate(fg_time = ifelse(died == 1 & incident_stroke == 0,
                          ifelse(!is.na(death_year) & death_year > base_iwy,
                                 death_year - base_iwy, pmax(fu_mid, 0.5)),
                          fu_mid),
         fg_stat = ifelse(incident_stroke == 1, 1, ifelse(died == 1, 2, 0)))
fgd$fg_stat[is.na(fgd$fg_time)] <- 0; fgd$fg_time[is.na(fgd$fg_time)] <- median(fgd$fg_time, na.rm = TRUE)
cat("\nFG status (0/1/2):\n"); print(table(fgd$fg_stat))
fgd$fg_stat <- factor(fgd$fg_stat)
fg <- finegray(Surv(fg_time, fg_stat) ~ ., data = fgd %>%
                 select(fg_time, fg_stat, hgi_hi_p90, all_of(covs), bw), etype = "1")
mfg <- coxph(Surv(fgstart, fgstop, fgstatus) ~ hgi_hi_p90 + age + sex + urban + edu +
               htn + dm + dyslip + heart + sbp + crp + tg + hdl + strata(bw),
             data = fg, weight = fgwt)
cat("\n== Fine-Gray (death data incomplete for W4/W6; interpret with caution) ==\n"); pr(mfg, "FG_SHR")

sink(); close(resf)
cat("\nDONE\n")
sink(); close(logf)
