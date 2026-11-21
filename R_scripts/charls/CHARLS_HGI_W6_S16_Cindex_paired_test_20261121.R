####################################################################################
# Step 16 | Paired-bootstrap P values for the head-to-head C-index comparison (M4)
# Reviewer request (R2-M4): the delta-C bootstrap CIs existed (B=500) but no P values.
# Same construction as Step 2 (M4): base model = 12 model-3 covariates; each marker
# added separately; Harrell C via survival::concordance; paired bootstrap on subjects.
# Two-sided P = 2 * min(P(dC <= 0), P(dC >= 0)) across bootstrap replicates.
# Output: results/charls/S16_cindex_paired_test.txt (+ csv)
####################################################################################
rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(survival); library(dplyr) })
OUT_DIR  <- "D:/CHARLS/HGI_stroke"
RES_DIR  <- "D:/CHARLS/HGI_stroke/freeze_20261120/results/charls"

d <- read.csv(file.path(OUT_DIR, "analytic_hgi_baseline_w6.csv"), colClasses = c(ID = "character"))
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
d <- d %>% mutate(shr = (glu / 18) / (1.59 * hba1c - 2.59),
                  tyg = log(tg * glu / 2))
dmain <- d %>% filter(complete.cases(across(all_of(c(cov_m3,"hgi","hba1c",
                                                     "incident_stroke6","fu_v2")))))
dmain4 <- dmain %>% filter(complete.cases(across(all_of(c("shr","tyg")))))
adj3 <- "age + sex + urban + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"
FS <- function(tv) paste("Surv(", tv, ", incident_stroke6)")
expos <- c(hgi = "hgi", HbA1c = "hba1c", FPG = "glu", SHR = "shr", TyG = "tyg")

cat("M4 sample:", nrow(dmain4), "n; events", sum(dmain4$incident_stroke6), "\n")

fits4 <- lapply(expos, function(x)
  coxph(as.formula(paste(FS("fu_v2"), "~", x, "+", adj3)), data = dmain4))
cidx <- sapply(fits4, function(f) concordance(f)$concordance)
cat("apparent C:", paste(names(cidx), round(cidx, 4), collapse = " | "), "\n")

suppressPackageStartupMessages(library(parallel))
set.seed(20260903); B <- 2000
seeds <- sample.int(2^30, B)
ncore <- max(1, detectCores() - 1)
cat("parallel bootstrap on", ncore, "cores\n")
one_rep <- function(sd_, dmain4, expos, adj3) {
  suppressPackageStartupMessages(library(survival))
  set.seed(sd_); ii <- sample.int(nrow(dmain4), replace = TRUE); db <- dmain4[ii, ]
  out <- rep(NA_real_, length(expos))
  for (k in seq_along(expos)) {
    fb <- try(coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~", expos[k], "+", adj3)),
                  data = db), silent = TRUE)
    if (!inherits(fb, "try-error")) out[k] <- concordance(fb)$concordance
  }
  out
}
cl <- makeCluster(ncore)
boot_list <- parLapply(cl, seeds, one_rep, dmain4 = dmain4, expos = expos, adj3 = adj3)
stopCluster(cl)
boot <- do.call(rbind, boot_list)
colnames(boot) <- names(expos)
dc <- sweep(boot, 1, boot[, "hgi"], "-")
out <- lapply(colnames(dc), function(mk) {
  x <- dc[, mk]; x <- x[!is.na(x)]
  p2 <- 2 * min(mean(x <= 0), mean(x >= 0))
  data.frame(marker = mk, dC = mean(x), lo = quantile(x, .025), hi = quantile(x, .975),
             P_two_sided = min(1, p2))
})
res <- do.call(rbind, out)
res <- res[res$marker != "hgi", ]
res$dC <- round(res$dC, 4); res$lo <- round(res$lo, 4); res$hi <- round(res$hi, 4)
res$P_two_sided <- signif(res$P_two_sided, 3)
print(res)

sink(file.path(RES_DIR, "S16_cindex_paired_test.txt"))
cat("S16 | Paired-bootstrap test of delta C-index vs HGI (reference)\n")
cat("Base model: 12 model-3 covariates; each marker added separately; same participants.\n")
cat("Sample:", nrow(dmain4), "participants,", sum(dmain4$incident_stroke6), "events. B =", B,
    "paired bootstrap resamples (subject-level), seed 20260903.\n")
cat("Apparent C-indices:", paste(names(cidx), round(cidx, 4), collapse = " | "), "\n\n")
cat("Two-sided P = 2*min(P(dC<=0), P(dC>=0)) across replicates.\n\n")
print(res)
sink()
write.csv(res, file.path(RES_DIR, "S16_cindex_paired_test.csv"), row.names = FALSE)
cat("saved S16_cindex_paired_test.txt/.csv\n")
