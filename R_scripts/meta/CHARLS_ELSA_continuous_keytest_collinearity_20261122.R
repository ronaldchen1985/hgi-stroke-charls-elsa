####################################################################################
# R1+R2 revision analyses (2026-11-22), mock-review driven
#  R1: continuous-form key test (exact FPG + HbA1c algebraic equivalent) in all
#      three analytic sets + non-diabetic counterparts; numeric equivalence check
#  R2: collinearity diagnostics (VIF, condition index) for binary and continuous
#      key-test models; fragility index for the significant binary key tests
# Inputs : freeze_20261120/data/{charls,elsa}  (+ raw CHARLS 2015 dta for landmark)
# Outputs: results/R1R2_continuous_keytest_20261122.csv
#          results/R1R2_collinearity_20261122.csv
#          results/R1R2_fragility_20261122.csv
#          results/R1R2_log_20261122.txt  (console mirror)
####################################################################################

rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(haven); library(survival); library(dplyr) })
BASE <- "D:/CHARLS/HGI_stroke/freeze_20261120"
OUT  <- file.path(BASE, "results")
set.seed(20260903)

logf <- file.path(OUT, "R1R2_log_20261122.txt")
sink(logf, split = TRUE)
on.exit({ sink(NULL) }, add = TRUE)

to_w1_id <- function(id) { id <- as.character(id)
  hit <- nchar(id) == 12 & substr(id, 10, 10) == "0"
  id[hit] <- paste0(substr(id[hit], 1, 9), substr(id[hit], 11, 12)); id }

grab <- function(fit, term) {
  ci <- summary(fit)$conf.int[term, , drop = FALSE]
  p  <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(term = term, HR = round(ci[,1], 4), lo = round(ci[,3], 4),
             hi = round(ci[,4], 4), P = signif(p, 4))
}

## generic VIF on the model's RHS design (linear-model analogue, standard practice)
vif_rhs <- function(data, rhs_terms) {
  X <- model.matrix(as.formula(paste("~", paste(rhs_terms, collapse = "+"))),
                    data = data)[, -1, drop = FALSE]
  vifs <- sapply(seq_len(ncol(X)), function(j) {
    r2 <- summary(lm(X[, j] ~ X[, -j]))$r.squared
    1 / (1 - r2)
  })
  names(vifs) <- colnames(X); vifs
}

## condition index of the scaled RHS design
condition_index <- function(data, rhs_terms) {
  X <- model.matrix(as.formula(paste("~", paste(rhs_terms, collapse = "+"))),
                    data = data)[, -1, drop = FALSE]
  Xs <- scale(X); sv <- svd(Xs)$d
  ci <- max(sv) / sv
  c(max = max(ci), second = sort(ci, decreasing = TRUE)[2])
}

cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
adj3   <- paste(cov_m3, collapse = " + ")

res_key  <- list()   # binary + continuous key tests, all sets
res_coll <- list()   # VIF / condition index
res_frag <- list()   # fragility index

run_block <- function(tag, data, surv, hi_var, hgi_var, a1c_var, fpg_var, covs,
                      fragility = FALSE) {
  adj <- paste(covs, collapse = " + ")
  cat("\n==================== ", tag, " ====================\n")
  cat("n =", nrow(data), "\n")

  # -- binary key test (anchor) -------------------------------------------------------
  f_bin <- as.formula(paste(surv, "~", hi_var, "+", a1c_var, "+", adj))
  m_bin <- coxph(f_bin, data = data)
  print(bind_rows(grab(m_bin, hi_var), grab(m_bin, a1c_var)))
  out_bin <- cbind(set = tag, model = "binary key test",
                   grab(m_bin, hi_var), note = "")

  # -- continuous key test: per-SD HGI + HbA1c ---------------------------------------
  sd_hgi <- sd(data[[hgi_var]], na.rm = TRUE)
  data$..hgi_z <- data[[hgi_var]] / sd_hgi
  f_con <- as.formula(paste(surv, "~ ..hgi_z +", a1c_var, "+", adj))
  m_con <- coxph(f_con, data = data)
  cat("continuous key test (HGI per SD =", round(sd_hgi, 4), "):\n")
  print(bind_rows(grab(m_con, "..hgi_z"), grab(m_con, a1c_var)))
  out_con <- cbind(set = tag, model = "continuous key test (per SD)",
                   grab(m_con, "..hgi_z"), note = paste0("SD=", round(sd_hgi, 4)))
  out_con$term <- "hgi_perSD"
  out_a1c  <- cbind(set = tag, model = "continuous key test: HbA1c term",
                    grab(m_con, a1c_var), note = "")

  # -- algebraic-equivalent form: FPG + HbA1c -----------------------------------------
  f_fpg <- as.formula(paste(surv, "~", fpg_var, "+", a1c_var, "+", adj))
  m_fpg <- coxph(f_fpg, data = data)
  bg_z  <- coef(m_con)["..hgi_z"]; gf <- coef(m_fpg)[fpg_var]
  bg    <- bg_z / sd_hgi                      # unscaled per-unit HGI coefficient
  b0    <- -gf / bg                           # implied construction slope of the HGI equation
  lp_d  <- max(abs(predict(m_con) - predict(m_fpg)))
  cat(sprintf("[equivalence check] logLik HGI-form = %.4f ; FPG-form = %.4f\n",
              logLik(m_con), logLik(m_fpg)))
  cat(sprintf("[equivalence check] max |linear predictor difference| = %.3e\n", lp_d))
  cat(sprintf("[equivalence check] implied construction slope b0 = -gamma_FPG/beta_HGI = %.6f (documented equation slope should match)\n",
              b0))
  out_fpg <- cbind(set = tag, model = "equivalent FPG form: FPG term",
                   grab(m_fpg, fpg_var), note = "FPG HR per 1 unit")
  out_fpgA <- cbind(set = tag, model = "equivalent FPG form: HbA1c term",
                    grab(m_fpg, a1c_var), note = "")

  res_key[[tag]] <- bind_rows(out_bin, out_con, out_a1c, out_fpg, out_fpgA)

  # -- collinearity diagnostics --------------------------------------------------------
  for (lbl in c("binary", "continuous")) {
    rhs <- if (lbl == "binary") c(hi_var, a1c_var, covs) else c(hgi_var, a1c_var, covs)
    dd  <- data[, c(rhs)]
    dd  <- dd[complete.cases(dd), , drop = FALSE]
    v   <- vif_rhs(dd, rhs)
    ci  <- condition_index(dd, rhs)
    key_vars <- c(hi_var, hgi_var, a1c_var)
    cat("[", lbl, "] VIF (key terms):\n"); print(round(v[intersect(names(v), key_vars)], 2))
    cat("[", lbl, "] max VIF =", round(max(v), 2), " (", names(which.max(v)), ") ;",
        " condition index max =", round(ci["max"], 1), "\n")
    res_coll[[paste(tag, lbl)]] <- data.frame(
      set = tag, model_form = lbl,
      VIF_HGI_term = round(v[if (lbl == "binary") hi_var else hgi_var], 3),
      VIF_HbA1c    = round(v[a1c_var], 3),
      max_VIF      = round(max(v), 3), max_VIF_term = names(which.max(v)),
      condition_index_max = round(ci["max"], 2),
      condition_index_2nd = round(ci["second"], 2))
  }

  # -- fragility index (only where the binary key test is significant) ----------------
  if (fragility) {
    p0 <- summary(m_bin)$coefficients[hi_var, "Pr(>|z|)"]
    tv <- sub("Surv\\((.*),.*", "\\1", surv); tv <- trimws(tv)
    ev <- trimws(sub(".*,", "", gsub("[()]", "", surv)))
    is_event <- data[[ev]] == 1 & data[[hi_var]] == 1
    idx_events <- which(is_event)
    ord_early <- idx_events[order(data[[tv]][idx_events])]
    ord_late  <- rev(ord_early)
    count_fi <- function(order_idx) {
      d2 <- data
      for (k in seq_along(order_idx)) {
        d2[[ev]][order_idx[k]] <- 0   # event -> censored at same time
        m2 <- try(coxph(f_bin, data = d2), silent = TRUE)
        if (inherits(m2, "try-error")) return(NA_integer_)
        p2 <- summary(m2)$coefficients[hi_var, "Pr(>|z|)"]
        hr2 <- summary(m2)$conf.int[hi_var, 1]
        if (p2 >= 0.05) return(c(k = k, HR = hr2, P = p2))
      }
      c(k = length(order_idx), HR = NA_real_, P = NA_real_)
    }
    fe <- count_fi(ord_early); fl <- count_fi(ord_late)
    n_exp_ev <- length(idx_events)
    cat(sprintf("[fragility] exposed events = %d; base P = %.4f\n", n_exp_ev, p0))
    cat(sprintf("[fragility] earliest-first: %d event(s) -> P = %.4f (HR %.3f)\n",
                fe["k"], fe["P"], fe["HR"]))
    cat(sprintf("[fragility] latest-first:   %d event(s) -> P = %.4f (HR %.3f)\n",
                fl["k"], fl["P"], fl["HR"]))
    res_frag[[tag]] <- data.frame(
      set = tag, exposed_events = n_exp_ev, base_P = round(p0, 4),
      FI_earliest_first = fe["k"], P_at_FI_early = round(fe["P"], 4),
      HR_at_FI_early = round(fe["HR"], 3),
      FI_latest_first = fl["k"], P_at_FI_late = round(fl["P"], 4),
      HR_at_FI_late = round(fl["HR"], 3),
      fragility_quotient_early = round(fe["k"] / n_exp_ev, 4))
  }
  invisible(TRUE)
}

## =================================================================================
## Set 1: CHARLS 2011 main set (full + non-diabetic)
## =================================================================================
d <- read.csv(file.path(BASE, "data/charls/analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
CUT_P90 <- as.numeric(quantile(d$hgi, 0.9, na.rm = TRUE))
d <- d %>% mutate(hgi_hi_p90 = as.integer(hgi >= CUT_P90))
dmain <- d %>% filter(complete.cases(across(all_of(
  c(cov_m3, "hgi", "hba1c", "incident_stroke6", "fu_v2")))))
cat("CHARLS main set:", nrow(dmain), "n; events", sum(dmain$incident_stroke6),
    "(anchor expect 8,684 / 852; binary key test HR 1.337, P = 0.024)\n")

run_block("CHARLS 2011 main", dmain, "Surv(fu_v2, incident_stroke6)",
          "hgi_hi_p90", "hgi", "hba1c", "glu", cov_m3, fragility = TRUE)

ndm_c <- dmain %>% filter(diabetes == 0)
cat("\nCHARLS non-DM set:", nrow(ndm_c), "n; events", sum(ndm_c$incident_stroke6),
    "(anchor expect 7,283 / 681; binary key test HR 1.524, P = 0.006)\n")
run_block("CHARLS 2011 non-DM", ndm_c, "Surv(fu_v2, incident_stroke6)",
          "hgi_hi_p90", "hgi", "hba1c", "glu", setdiff(cov_m3, "dm"),
          fragility = TRUE)

## =================================================================================
## Set 2: CHARLS 2015 landmark set (rebuilt per frozen step-10 code)
## =================================================================================
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
adjv <- cov_m3
dd <- dd %>%
  select(-any_of(intersect(adjv, setdiff(names(dd), "ID")))) %>%
  left_join(d %>% select(ID, stroke_w6, last_contact, w6_report_yr, w6_dxyear,
                         all_of(adjv)), by = "ID") %>%
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
CUT15 <- as.numeric(quantile(dd$hgi_2015, 0.9, na.rm = TRUE))
dd <- dd %>% mutate(hgi15_hi_p90 = as.integer(hgi_2015 >= CUT15))
dd2 <- dd %>%
  select(-any_of(c("htn","dm","dyslip","heart","sbp","crp","tg","hdl"))) %>%
  left_join(w3_dis, by = "ID") %>% left_join(w3_bp, by = "ID") %>%
  left_join(w3_bl, by = "ID")
cov15 <- c("age15","sex","urban","edu","htn15","dm15","dyslip15","heart15",
           "sbp15","crp15","tg15","hdl15")
dmain15 <- dd2 %>% filter(complete.cases(across(all_of(
  c(cov15, "hgi_2015", "hba1c3", "ev15", "t15")))))
cat("\nCHARLS 2015 landmark set:", nrow(dmain15), "n; events", sum(dmain15$ev15),
    "(anchor expect 6,581 / 706; binary key test HR 1.179, P = 0.253)\n")

run_block("CHARLS 2015 landmark", dmain15, "Surv(t15, ev15)",
          "hgi15_hi_p90", "hgi_2015", "hba1c3", "glu3", cov15)

ndm15 <- dmain15 %>% filter(dm15 == 0 & hba1c3 < 6.5 & glu3 < 126)
cat("\nCHARLS 2015 landmark non-DM set:", nrow(ndm15), "n; events", sum(ndm15$ev15), "\n")
run_block("CHARLS 2015 landmark non-DM", ndm15, "Surv(t15, ev15)",
          "hgi15_hi_p90", "hgi_2015", "hba1c3", "glu3", setdiff(cov15, "dm15"))

## =================================================================================
## Set 3: ELSA W2 (full + non-diabetic)
## =================================================================================
e <- read.csv(file.path(BASE, "data/elsa/analytic_elsa_hgi_baseline_w2.csv"))
cat("\nELSA W2 set:", nrow(e), "n; events", sum(e$incident_stroke),
    "(anchor expect 2,577 / 176; binary key test HR 1.440, P = 0.180)\n")
run_block("ELSA W2", e, "Surv(fu_mid, incident_stroke)",
          "hgi_hi_p90", "hgi", "hba1c", "fglu", cov_m3)

en <- e %>% mutate(ndm = as.integer(dm == 0 & hba1c < 6.5 & fglu < 7.0)) %>%
  filter(ndm == 1)
cat("\nELSA W2 non-DM set:", nrow(en), "n; events", sum(en$incident_stroke),
    "(anchor expect 2,474 / 164; binary key test HR 1.425)\n")
run_block("ELSA W2 non-DM", en, "Surv(fu_mid, incident_stroke)",
          "hgi_hi_p90", "hgi", "hba1c", "fglu", setdiff(cov_m3, "dm"))

## =================================================================================
## write out
## =================================================================================
rk <- bind_rows(res_key);  rownames(rk) <- NULL
rc <- bind_rows(res_coll); rownames(rc) <- NULL
write.csv(rk, file.path(OUT, "R1R2_continuous_keytest_20261122.csv"), row.names = FALSE)
write.csv(rc, file.path(OUT, "R1R2_collinearity_20261122.csv"),  row.names = FALSE)
if (length(res_frag)) {
  rf <- bind_rows(res_frag); rownames(rf) <- NULL
  write.csv(rf, file.path(OUT, "R1R2_fragility_20261122.csv"), row.names = FALSE)
  print(rf)
}
cat("\n[done] outputs written to", OUT, "\n")
