####################################################################################
# ELSA external validation -- Step 1: build analytic dataset (baseline = Wave 2)
# Mirrors CHARLS HGI x stroke pipeline. Output:
#   analytic_elsa_hgi_baseline_w2.csv + build_log_ELSA.txt
# Sources:
#   gh_elsa_h.dta (Harmonized ELSA, waves 1-10): outcome/covariates/interview info
#   wave_2_nurse_data_v2.dta: biomarkers (fglu, hba1c, lipids, hscrp, BP)
#   wave_11_elsa_data_eul.dta: W11 incident stroke (hediiam7), interview date
#   elsa_geog_urindewr_2011_eul.dta: urban/rural (w6 proxy for baseline)
#   h_elsa_eol_a2.dta: death year (raxyear)
####################################################################################
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(haven); library(dplyr) })
D <- "D:/CHARLS/HGI_stroke/ELSA验证/data_raw/stata/stata13_se/"
O <- "D:/CHARLS/HGI_stroke/ELSA验证/"
logf <- file(paste0(O, "build_log_ELSA.txt"), open = "wt"); sink(logf, split = TRUE)

na_neg <- function(x) { x[x < 0] <- NA; x }

## ---- 1. Harmonized core -------------------------------------------------------------
gh <- read_dta(paste0(D, "gh_elsa_h.dta"),
               col_select = c("idauniq", "ragender", "r2agey", "raeducl",
                              paste0("r", 1:10, "stroke"),
                              "r2hibpe", "r2diabe", "r2hearte", "r2hchole",
                              paste0("r", 1:10, "iwstat"), paste0("r", 1:10, "iwy"))) %>%
  mutate(across(-idauniq, as.numeric),
         across(-idauniq, na_neg))
cat("harmonized:", nrow(gh), "rows\n")

## ---- 2. W2 nurse biomarkers -----------------------------------------------------------
nu <- read_dta(paste0(D, "wave_2_nurse_data_v2.dta"),
               col_select = tidyselect::any_of(c("idauniq", "fglu", "hba1c", "chol", "hdl", "ldl",
                              "trig", "hscrp", "sysval", "diaval", "fasteli",
                              "wbc", "cfib"))) %>%
  mutate(across(-idauniq, as.numeric), across(-idauniq, na_neg))
cat("W2 nurse:", nrow(nu), "rows\n")
cat("fglu: n =", sum(!is.na(nu$fglu)), " median =", median(nu$fglu, na.rm = TRUE),
    " (expect mmol/L ~5-6)\n")
cat("fglu missing by fasteli:\n"); print(table(nu$fasteli, is.na(nu$fglu), useNA = "always"))
cat("hba1c: n =", sum(!is.na(nu$hba1c)), " median =", median(nu$hba1c, na.rm = TRUE),
    " range =", paste(range(nu$hba1c, na.rm = TRUE), collapse = "-"), " (expect % ~5-7)\n")
cat("trig median =", median(nu$trig, na.rm = TRUE), "; hdl median =",
    median(nu$hdl, na.rm = TRUE), " (expect mmol/L)\n")
cat("hscrp median =", median(nu$hscrp, na.rm = TRUE), " (expect mg/L)\n")
cat("sysval median =", median(nu$sysval, na.rm = TRUE), "\n")

## ---- 3. W11 incident stroke + interview date -----------------------------------------
w11 <- read_dta(paste0(D, "wave_11_elsa_data_eul.dta"),
                col_select = c("idauniq", "hediia7", "hediacst", "w11indout",
                               "iintdaty", "iintdatm")) %>%
  mutate(across(-idauniq, as.numeric), across(-idauniq, na_neg))
cat("W11: incident stroke (hediia7==1) n =", sum(w11$hediia7 == 1, na.rm = TRUE), "\n")

## ---- 4. urban/rural + death year -------------------------------------------------------
geo <- read_dta(paste0(D, "elsa_geog_urindewr_2001_eul.dta"),
                col_select = c("idauniq", "w2_urindewr_2001")) %>%
  mutate(across(-idauniq, as.numeric), across(-idauniq, na_neg))
eol <- read_dta(paste0(D, "h_elsa_eol_a2.dta"),
                col_select = c("idauniq", "raxyear")) %>%
  mutate(across(-idauniq, as.numeric), across(-idauniq, na_neg))
cat("EOL death year non-missing:", sum(!is.na(eol$raxyear)), "\n")

## ---- 5. merge -------------------------------------------------------------------------
df <- gh %>%
  left_join(nu, by = "idauniq") %>%
  left_join(w11, by = "idauniq") %>%
  left_join(geo, by = "idauniq") %>%
  left_join(eol, by = "idauniq")
cat("merged:", nrow(df), "\n")

## ---- 6. flow --------------------------------------------------------------------------
cat("\n== STROBE flow ==\n")
n0 <- nrow(df); cat("0. harmonized total:", n0, "\n")
df <- df %>% filter(!is.na(fglu) & !is.na(hba1c))
cat("1. with W2 fglu + hba1c:", nrow(df),
    " (fasting-eligible blood subsample by design)\n")
df <- df %>% filter(r2agey >= 50)
cat("2. age >= 50 at W2:", nrow(df), "\n")
df <- df %>% filter(is.na(r2stroke) | r2stroke == 0)
cat("3. no prevalent stroke at W2:", nrow(df), "\n")

# last interview year after W2 (waves 3..11)
iwy_cols <- paste0("r", 3:10, "iwy")
iwstat_cols <- paste0("r", 3:10, "iwstat")
df$w11_iwy <- df$iintdaty
get_last <- function(yrs) {
  apply(yrs, 1, function(r) { v <- which(!is.na(r)); if (length(v)) max(r[v]) else NA_real_ })
}
# last interviewed wave year (iwstat==1 counts as interviewed; W11: w11indout 11-25 interviewed)
iy <- df[, iwy_cols]
iy[df[, iwstat_cols] != 1] <- NA
last10 <- get_last(iy)
df$w11_int <- as.integer(!is.na(df$w11indout) & df$w11indout >= 11 & df$w11indout <= 25)
df$last_contact <- pmax(last10,
                        ifelse(df$w11_int == 1, df$w11_iwy, NA_real_), na.rm = TRUE)
df$last_contact[is.infinite(df$last_contact)] <- NA
cat("4. at least one interview after W2:", sum(!is.na(df$last_contact)), "\n")
df <- df %>% filter(!is.na(last_contact))

## ---- 7. incident stroke ----------------------------------------------------------------
# first post-baseline wave with stroke==1 (ever-measure); W11 from hediiam7
sw <- df[, paste0("r", 3:10, "stroke")]
first_ev10 <- apply(sw, 1, function(r) { v <- which(r == 1); if (length(v)) min(v) else NA_integer_ })
df$ev_wave10 <- ifelse(is.na(first_ev10), NA, first_ev10 + 2)  # wave number 3..10
df$ev_wave <- df$ev_wave10
df$ev_wave[is.na(df$ev_wave) & !is.na(df$hediia7) & df$hediia7 == 1] <- 11
df$incident_stroke <- as.integer(!is.na(df$ev_wave))
cat("5. incident strokes W3-W11:", sum(df$incident_stroke), "\n")

# reporting year and previous stroke-free year
wave_yr <- function(w, i) {
  if (is.na(w)) return(NA_real_)
  if (w <= 10) return(as.numeric(df[[paste0("r", w, "iwy")]][i]))
  return(as.numeric(df$w11_iwy[i]))
}
prev_free_yr <- function(w, i) {
  # last wave before event wave with stroke==0 (interviewed, no stroke)
  prevw <- w - 1
  while (prevw >= 3) {
    s <- df[[paste0("r", prevw, "stroke")]][i]
    if (!is.na(s) && s == 0) return(as.numeric(df[[paste0("r", prevw, "iwy")]][i]))
    prevw <- prevw - 1
  }
  return(as.numeric(df$r2iwy[i]))  # W2 interview year as pre-baseline reference
}
rep_y <- mapply(wave_yr, df$ev_wave, seq_len(nrow(df)))
prev_y <- mapply(function(w, i) if (is.na(w)) NA_real_ else prev_free_yr(w, i),
                 df$ev_wave, seq_len(nrow(df)))
df$event_year_rep <- rep_y
df$event_year_prev <- prev_y
df$fu_mid <- ifelse(df$incident_stroke == 1,
                    (rep_y + prev_y) / 2 - df$r2iwy,
                    df$last_contact - df$r2iwy)
df <- df %>% filter(is.na(fu_mid) | fu_mid > 0)
cat("6. analysis-eligible:", nrow(df), "events:", sum(df$incident_stroke), "\n")

## ---- 8. HGI ----------------------------------------------------------------------------
f_eq <- lm(hba1c ~ fglu, data = df)
cat("\n== HGI equation (ELSA, W2) ==\n"); print(coef(f_eq))
df$hgi_pred <- predict(f_eq, newdata = df)
df$hgi <- df$hba1c - df$hgi_pred
CUT <- as.numeric(quantile(df$hgi, 0.9, na.rm = TRUE))
df$hgi_hi_p90 <- as.integer(df$hgi >= CUT)
cat("P90 cutoff:", round(CUT, 4), "; high-HGI n =", sum(df$hgi_hi_p90), "\n")
df$shr <- df$fglu / (1.59 * df$hba1c - 2.59)
df$tyg <- log(df$trig * df$fglu * 88.57 / 2)  # trig mmol/L->mg/dL (x88.57), glu mmol/L->mg/dL (x18)
df$tyg <- log((df$trig * 88.57) * (df$fglu * 18) / 2)

## ---- 9. covariates + death -------------------------------------------------------------
died_flag <- as.integer(apply(df[, iwstat_cols], 1, function(r) any(r %in% c(5, 6)))) |
  as.integer(!is.na(df$w11indout) & df$w11indout == 95)
df <- df %>% mutate(
  age = r2agey, sex = ragender - 1,  # 0 man, 1 woman (CHARLS-style numeric)
  edu = raeducl, urban = as.integer(w2_urindewr_2001 == 1),
  htn = r2hibpe, dm = r2diabe, dyslip = r2hchole, heart = r2hearte,
  sbp = sysval, crp = hscrp, tg = trig * 88.57, hdl = hdl * 38.67,  # to mg/dL
  died = died_flag,
  death_year = raxyear,
  fu_years_total = last_contact - r2iwy)
cat("deaths:", sum(df$died, na.rm = TRUE), "\n")

covs <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
main <- df %>% filter(complete.cases(across(all_of(c(covs, "hgi", "hba1c",
                                                     "incident_stroke", "fu_mid")))))
cat("\n== main complete-case set ==\n")
cat("n =", nrow(main), "; events =", sum(main$incident_stroke),
    " (", round(mean(main$incident_stroke) * 100, 1), "%); deaths =",
    sum(main$died, na.rm = TRUE),
    "; median fu =", round(median(main$fu_years_total), 1), "y\n")
cat("high HGI n =", sum(main$hgi_hi_p90), "\n")
cat("missingness by covariate (analysis-eligible):\n")
print(sapply(df[covs], function(x) sum(is.na(x))))

write.csv(df, paste0(O, "analytic_elsa_hgi_full_w2.csv"), row.names = FALSE)
write.csv(main, paste0(O, "analytic_elsa_hgi_baseline_w2.csv"), row.names = FALSE)
cat("\nDONE\n")
sink(); close(logf)
