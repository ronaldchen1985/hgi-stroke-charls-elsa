####################################################################################
# S5 revision analysis (2026-11-22), mock-review driven
#  Unified ABSOLUTE cut-off sensitivity: HGI >= 0.5 (same threshold in all sets),
#  M3-adjusted Cox in CHARLS 2011 main / CHARLS 2015 landmark / ELSA W2,
#  plus fixed-effect pooled estimate; non-diabetic counterparts included.
#  Purpose: address "cohort-specific P90 cut-offs hamper comparability" (S5).
# Inputs : freeze_20261120/data/{charls,elsa} (+ raw CHARLS 2015 dta for landmark)
# Outputs: results/S5_absolute_cutoff_20261122.csv
#          results/S5_log_20261122.txt
####################################################################################

rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(haven); library(survival); library(dplyr) })
BASE <- "D:/CHARLS/HGI_stroke/freeze_20261120"
OUT  <- file.path(BASE, "results")
set.seed(20260903)

logf <- file.path(OUT, "S5_log_20261122.txt")
sink(logf, split = TRUE)
on.exit({ sink(NULL) }, add = TRUE)

to_w1_id <- function(id) { id <- as.character(id)
  hit <- nchar(id) == 12 & substr(id, 10, 10) == "0"
  id[hit] <- paste0(substr(id[hit], 1, 9), substr(id[hit], 11, 12)); id }

grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(HR = round(ci[,1], 4), lo = round(ci[,3], 4),
             hi = round(ci[,4], 4), P = signif(p, 4))
}

pool_fe <- function(betas, ses) {
  w <- 1 / ses^2
  b <- sum(w * betas) / sum(w); se <- sqrt(1 / sum(w))
  z <- b / se
  q <- sum(w * (betas - b)^2); dfq <- length(betas) - 1
  i2 <- max(0, (q - dfq) / q) * 100
  data.frame(HR = round(exp(b), 4), lo = round(exp(b - 1.96 * se), 4),
             hi = round(exp(b + 1.96 * se), 4), P = signif(2 * pnorm(-abs(z)), 4),
             I2 = round(i2, 1))
}

cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")

run_one <- function(tag, data, surv, hgi_var, covs) {
  adj <- paste(covs, collapse = " + ")
  f <- as.formula(paste(surv, "~ hgi_ge05 +", adj))
  m <- coxph(f, data = data)
  g <- grab(m, "hgi_ge05")
  n_hi <- sum(data$hgi_ge05 == 1); ev_hi <- sum(data$hgi_ge05 == 1 & data[[trimws(sub(".*,", "", gsub("[()]", "", surv)))]] == 1)
  cat(sprintf("%-32s n=%6d  exposed(HGI>=0.5)=%4d (%.1f%%)  exposed events=%3d  HR %.3f (%.3f-%.3f) P=%.4f\n",
              tag, nrow(data), n_hi, 100 * n_hi / nrow(data), ev_hi, g$HR, g$lo, g$hi, g$P))
  data.frame(set = tag, n = nrow(data), n_exposed = n_hi,
             pct_exposed = round(100 * n_hi / nrow(data), 1),
             events_exposed = ev_hi, g, beta = log(g$HR),
             se = (log(g$hi) - log(g$lo)) / (2 * 1.96))
}

res <- list()

## ---- CHARLS 2011 main ---------------------------------------------------------
d <- read.csv(file.path(BASE, "data/charls/analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
dmain <- d %>% mutate(hgi_ge05 = as.integer(hgi >= 0.5)) %>%
  filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c", "incident_stroke6", "fu_v2")))))
res[["CHARLS 2011 main (>=0.5)"]] <- run_one("CHARLS 2011 main (>=0.5)", dmain,
  "Surv(fu_v2, incident_stroke6)", "hgi", cov_m3)
ndm <- dmain %>% filter(diabetes == 0)
res[["CHARLS 2011 non-DM (>=0.5)"]] <- run_one("CHARLS 2011 non-DM (>=0.5)", ndm,
  "Surv(fu_v2, incident_stroke6)", "hgi", setdiff(cov_m3, "dm"))

## ---- CHARLS 2015 landmark (rebuilt per frozen step-10 code) -------------------
hs <- read_dta("D:/CHARLS/2015/Health_Status_and_Functioning.dta")
mk <- function(z, a) as.integer((!is.na(z) & z == 1) | (!is.na(a) & a == 1))
w3_dis <- data.frame(ID = to_w1_id(hs$ID),
                     htn15    = mk(hs$zda007_1_, hs$da007_1_),
                     dyslip15 = mk(hs$zda007_2_, hs$da007_2_),
                     dm15     = mk(hs$zda007_3_, hs$da007_3_),
                     heart15  = mk(hs$zda007_7_, hs$da007_7_))
bm <- read_dta("D:/CHARLS/2015/Biomarker.dta")
sbp_mat <- sapply(bm[c("qa003","qa007","qa011")], function(x) {
  x <- as.numeric(x); x[x < 60 | x > 260] <- NA; x })
w3_bp <- data.frame(ID = to_w1_id(bm$ID), sbp15 = rowMeans(sbp_mat, na.rm = TRUE))
w3_bp$sbp15[is.nan(w3_bp$sbp15)] <- NA
bl <- read_dta("D:/CHARLS/2015/Blood.dta")
w3_bl <- data.frame(ID = to_w1_id(bl$ID),
                    crp15 = as.numeric(bl$bl_crp), tg15 = as.numeric(bl$bl_tg),
                    hdl15 = as.numeric(bl$bl_hdl), glu15 = as.numeric(bl$bl_glu),
                    hba1c15 = as.numeric(bl$bl_hbalc))

dd <- read.csv(file.path(BASE, "data/charls/analytic_delta_hgi.csv"),
               colClasses = c(ID = "character"))
dd <- dd %>%
  select(-any_of(intersect(cov_m3, setdiff(names(dd), "ID")))) %>%
  left_join(d %>% select(ID, stroke_w6, last_contact, w6_report_yr, w6_dxyear,
                         all_of(cov_m3)), by = "ID") %>%
  left_join(w3_bl %>% select(ID, glu15, hba1c15), by = "ID") %>%
  filter(is.na(first_stroke_wave) | first_stroke_wave >= 4) %>%
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
  filter(!is.na(t15) & t15 > 0) %>%
  rename(glu3 = glu15, hba1c3 = hba1c15)
dd2 <- dd %>%
  select(-any_of(c("htn","dm","dyslip","heart","sbp","crp","tg","hdl"))) %>%
  left_join(w3_dis, by = "ID") %>% left_join(w3_bp, by = "ID") %>%
  left_join(w3_bl, by = "ID")
cov15 <- c("age15","sex","urban","edu","htn15","dm15","dyslip15","heart15",
           "sbp15","crp15","tg15","hdl15")
dmain15 <- dd2 %>% mutate(hgi_ge05 = as.integer(hgi_2015 >= 0.5)) %>%
  filter(complete.cases(across(all_of(c(cov15, "hgi_2015", "hba1c3", "ev15", "t15")))))
res[["CHARLS 2015 landmark (>=0.5)"]] <- run_one("CHARLS 2015 landmark (>=0.5)", dmain15,
  "Surv(t15, ev15)", "hgi_2015", cov15)
ndm15 <- dmain15 %>% filter(dm15 == 0 & hba1c3 < 6.5 & glu3 < 126)
res[["CHARLS 2015 landmark non-DM (>=0.5)"]] <- run_one("CHARLS 2015 landmark non-DM (>=0.5)", ndm15,
  "Surv(t15, ev15)", "hgi_2015", setdiff(cov15, "dm15"))

## ---- ELSA W2 -------------------------------------------------------------------
e <- read.csv(file.path(BASE, "data/elsa/analytic_elsa_hgi_baseline_w2.csv")) %>%
  mutate(hgi_ge05 = as.integer(hgi >= 0.5))
res[["ELSA W2 (>=0.5)"]] <- run_one("ELSA W2 (>=0.5)", e,
  "Surv(fu_mid, incident_stroke)", "hgi", cov_m3)
en <- e %>% filter(dm == 0 & hba1c < 6.5 & fglu < 7.0)
res[["ELSA W2 non-DM (>=0.5)"]] <- run_one("ELSA W2 non-DM (>=0.5)", en,
  "Surv(fu_mid, incident_stroke)", "hgi", setdiff(cov_m3, "dm"))

## ---- pooled (full sets / non-DM sets) ------------------------------------------
rr <- bind_rows(res)
full <- rr[!grepl("non-DM", rr$set), ]; nondm <- rr[grepl("non-DM", rr$set), ]
p_full  <- cbind(set = "Pooled (3 sets, fixed effect)", pool_fe(full$beta, full$se))
p_nondm <- cbind(set = "Pooled non-DM (3 sets, fixed effect)", pool_fe(nondm$beta, nondm$se))
cat("\n--- pooled ---\n"); print(p_full); print(p_nondm)

out <- bind_rows(rr %>% select(-beta, -se),
                 p_full %>% mutate(n = NA, n_exposed = NA, pct_exposed = NA, events_exposed = NA),
                 p_nondm %>% mutate(n = NA, n_exposed = NA, pct_exposed = NA, events_exposed = NA))
write.csv(out, file.path(OUT, "S5_absolute_cutoff_20261122.csv"), row.names = FALSE)
cat("\n[done] S5_absolute_cutoff_20261122.csv written\n")
