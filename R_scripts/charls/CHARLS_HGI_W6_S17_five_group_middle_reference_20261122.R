####################################################################################
# S17 (2026-11-22): five-group exposure with MIDDLE reference (P25-P75)
# Rationale: the binary high-vs-rest contrast lets the artifact-contaminated
# left tail (anemia-driven spuriously low HbA1c) inflate the reference risk,
# diluting the etiologic contrast of interest (right tail vs HGI ~ 0).
# Groups: <P10 | P10-25 | P25-75 (ref) | P75-90 | >=P90
# Sets: CHARLS main, CHARLS non-diabetic, ELSA W2, ELSA non-diabetic.
####################################################################################
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(dplyr); library(survival) })
F <- "D:/CHARLS/HGI_stroke/freeze_20261120"
OUT <- file.path(F, "results")

mk5 <- function(dd) {
  q <- quantile(dd$hgi, c(.10, .25, .75, .90), na.rm = TRUE)
  dd %>% mutate(g5 = cut(hgi, breaks = c(-Inf, q[1], q[2], q[3], q[4], Inf),
                         labels = c("<P10", "P10-25", "P25-75", "P75-90", ">=P90")) %>%
                  relevel(factor(.), ref = "P25-75"))
}
adj_c <- "age + sex + urban + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"
adj_c_nd <- "age + sex + urban + edu + htn + dyslip + heart + sbp + crp + tg + hdl"
adj_e <- "age + sex + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"
adj_e_nd <- "age + sex + edu + htn + dyslip + heart + sbp + crp + tg + hdl"

run_set <- function(dd, adj, label, fml_surv = "Surv(fu, ev)") {
  dd <- mk5(dd)
  tab <- dd %>% group_by(g5) %>% summarize(n = n(), ev = sum(ev), .groups = "drop")
  cat(sprintf("\n==== %s ====\n", label)); print(as.data.frame(tab))
  fit <- coxph(as.formula(paste(fml_surv, "~ g5 +", adj)), data = dd)
  s <- summary(fit)$conf.int
  rn <- rownames(s)
  out <- data.frame(set = label, group = sub("^g5", "", rn),
                    HR = s[, "exp(coef)"], lo = s[, "lower .95"], hi = s[, "upper .95"],
                    P = summary(fit)$coefficients[, "Pr(>|z|)"], row.names = NULL)
  print(out, digits = 3)
  list(tab = tab, est = out)
}
fmt <- function(e) with(e, sprintf("%s: HR %.3f (%.3f-%.3f), P=%.4g", group, HR, lo, hi, P))

res <- list()

## ---- CHARLS main set ----
d <- read.csv(file.path(F, "data/charls/analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
d <- d %>% filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c",
                                                 "incident_stroke6", "fu_v2"))))) %>%
  rename(ev = incident_stroke6, fu = fu_v2)
stopifnot(nrow(d) == 8684, sum(d$ev) == 852)
res$charls <- run_set(d, adj_c, "CHARLS 2011 main (M3)")

## ---- CHARLS non-diabetic ----
dn <- d %>% filter(diabetes == 0)
stopifnot(nrow(dn) == 7283, sum(dn$ev) == 681)
res$charls_nd <- run_set(dn, adj_c_nd, "CHARLS non-diabetic (M3)")

## ---- ELSA W2 ----
e <- read.csv(file.path(F, "data/elsa/analytic_elsa_hgi_baseline_w2.csv")) %>%
  rename(ev = incident_stroke, fu = fu_mid)
adj_e_full <- "age + sex + urban + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"
adj_e_nd   <- "age + sex + urban + edu + htn + dyslip + heart + sbp + crp + tg + hdl"
res$elsa <- run_set(e, adj_e_full, "ELSA W2 (M3)")
en <- e %>% filter(dm == 0 & hba1c < 6.5 & fglu < 7.0)
stopifnot(nrow(en) == 2474, sum(en$ev) == 164)
res$elsa_nd <- run_set(en, adj_e_nd, "ELSA W2 non-diabetic (M3)")

## ---- save ----
est_all <- bind_rows(lapply(res, `[[`, "est"))
write.csv(est_all, file.path(OUT, "S17_five_group_middle_reference.csv"), row.names = FALSE)
cat("\nSaved:", file.path(OUT, "S17_five_group_middle_reference.csv"), "\n")
cat("\n==== KEY CONTRASTS (vs P25-75 reference) ====\n")
for (nm in names(res)) { cat("\n", nm, "\n"); cat(fmt(res[[nm]]$est), sep = "\n") }
