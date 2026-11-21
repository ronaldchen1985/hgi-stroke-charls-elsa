####################################################################################
# Cross-cohort pooling (fixed/random effects) of fully-adjusted HGI-stroke HRs
# Reads frozen result objects; writes cross_cohort_meta.txt + cross_cohort_meta.rds
# Part of freeze_20261120. Seed-free (deterministic).
####################################################################################
rm(list = ls()); gc()
suppressPackageStartupMessages(library(dplyr))
F <- "D:/CHARLS/HGI_stroke/freeze_20261120"

## ---- extract log-HR and SE from frozen results -----------------------------------------
# CHARLS 2011 main (M3 full adjustment)
charls <- readRDS(file.path(F, "results/charls/formal_results_W6.rds"))$M1 %>%
  filter(model == "m3_full")
# CHARLS 2015 landmark (M3)
iv <- readRDS(file.path(F, "results/charls/internal_validation_2015.rds"))$M3
# ELSA W2 external (M3 full adjustment)
elsa <- readRDS(file.path(F, "results/elsa/ELSA_results.rds"))$M1 %>%
  filter(model == "m3_full")

est <- bind_rows(
  charls %>% transmute(cohort = "CHARLS 2011 main",      HR, lo, hi),
  iv     %>% transmute(cohort = "CHARLS 2015 landmark",  HR, lo, hi),
  elsa   %>% transmute(cohort = "ELSA W2 external",      HR, lo, hi)) %>%
  mutate(b = log(HR), se = (log(hi) - log(lo)) / (2 * 1.96))
print(est)

pool <- function(x, tag) {
  w  <- 1 / x$se^2
  b_fe <- sum(w * x$b) / sum(w); se_fe <- sqrt(1 / sum(w))
  Q  <- sum(w * (x$b - b_fe)^2); dfq <- nrow(x) - 1
  pQ <- pchisq(Q, dfq, lower.tail = FALSE)
  I2 <- max(0, (Q - dfq) / Q) * 100
  tausq <- max(0, (Q - dfq) / (sum(w) - sum(w^2) / sum(w)))
  wr <- 1 / (x$se^2 + tausq)
  b_re <- sum(wr * x$b) / sum(wr); se_re <- sqrt(1 / sum(wr))
  z  <- b_re / se_re; p <- 2 * pnorm(-abs(z))
  data.frame(set = tag,
             FE_HR = round(exp(b_fe), 3), FE_lo = round(exp(b_fe - 1.96 * se_fe), 3),
             FE_hi = round(exp(b_fe + 1.96 * se_fe), 3),
             RE_HR = round(exp(b_re), 3), RE_lo = round(exp(b_re - 1.96 * se_re), 3),
             RE_hi = round(exp(b_re + 1.96 * se_re), 3),
             P = signif(p, 3), Q = round(Q, 2), P_het = signif(pQ, 3), I2 = round(I2))
}

res <- bind_rows(pool(est[1:2, ], "CHARLS 2011+2015 (internal)"),
                 pool(est[c(1, 3), ], "CHARLS 2011 + ELSA W2 (external)"),
                 pool(est, "All three sets"))
print(res)

sink(file.path(F, "results/meta/cross_cohort_meta.txt"))
cat("Cross-cohort pooling of fully-adjusted HGI(high P90)-stroke HRs,", format(Sys.time(), "%F %T"), "\n\n")
cat("Input estimates:\n"); print(est)
cat("\nPooled:\n"); print(res)
sink()
saveRDS(list(estimates = est, pooled = res), file.path(F, "results/meta/cross_cohort_meta.rds"))
cat("DONE\n")
