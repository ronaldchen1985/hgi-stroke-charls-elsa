####################################################################################
# Fig 3 v3 | Dose-response of the HGI-stroke association (CHARLS main set), 2x2
# Columns = population layer: left all participants, right without diabetes.
# Rows = evidence layer: top fully adjusted RCS, bottom observed decile rates.
# a full RCS; b non-diabetic RCS; c full deciles; d non-diabetic deciles.
####################################################################################
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(survival);
                                 library(splines); library(patchwork) })
SKILL <- "C:/Users/ronaldchen/AppData/Roaming/kimi-desktop/daimon-share/daimon/runtime/kimi-code/home/plugins/managed/nature-skills@nature-skills/skills/nature-figure"
source("D:/CHARLS/HGI_stroke/freeze_20261120/figures/_align_helper.R")
OUT <- "D:/CHARLS/HGI_stroke/freeze_20261120/figures"
F <- "D:/CHARLS/HGI_stroke/freeze_20261120"

d <- read.csv(file.path(F, "data/charls/analytic_hgi_baseline_w6.csv"),
              colClasses = c(ID = "character"))
cov_m3 <- c("age","sex","urban","edu","htn","dm","dyslip","heart","sbp","crp","tg","hdl")
d <- d %>% filter(complete.cases(across(all_of(c(cov_m3, "hgi", "hba1c",
                                                 "incident_stroke6", "fu_v2")))))
stopifnot(nrow(d) == 8684, sum(d$incident_stroke6) == 852)
adj3 <- "age + sex + urban + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"
adj3_nd <- "age + sex + urban + edu + htn + dyslip + heart + sbp + crp + tg + hdl"
dn <- d %>% filter(diabetes == 0)
stopifnot(nrow(dn) == 7283, sum(dn$incident_stroke6) == 681)

make_rcs <- function(dd, adj) {
  B <- ns(dd$hgi, df = 4)
  fit <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ B +", adj)), data = dd)
  ref <- median(dd$hgi)
  grid_x <- seq(quantile(dd$hgi, 0.01), quantile(dd$hgi, 0.99), length.out = 200)
  B_ref <- predict(B, ref); B_new <- predict(B, grid_x)
  bb <- coef(fit)[paste0("B", 1:4)]
  V  <- vcov(fit)[paste0("B", 1:4), paste0("B", 1:4)]
  XB <- sweep(B_new, 2, 1, "*") - sweep(matrix(B_ref, nrow(B_new), 4, byrow = TRUE), 2, 1, "*")
  lp <- as.vector(XB %*% bb); se <- sqrt(rowSums((XB %*% V) * XB))
  m_lin <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ hgi +", adj)), data = dd)
  p_nl <- anova(m_lin, fit)[2, "Pr(>|Chi|)"]
  list(dat = data.frame(hgi = grid_x, HR = exp(lp), lo = exp(lp - 1.96 * se),
                        hi = exp(lp + 1.96 * se)),
       p_nl = p_nl, cut90 = as.numeric(quantile(dd$hgi, 0.9)), dd = dd)
}

plot_rcs <- function(obj, hdr, ylo = 0.7) {
  rcs <- obj$dat
  ggplot(rcs, aes(hgi, HR)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    geom_vline(xintercept = obj$cut90, linetype = "dotted", colour = "#B2182B", linewidth = 0.35) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#B2182B", alpha = 0.14) +
    geom_line(colour = "#B2182B", linewidth = 0.7) +
    geom_rug(data = { set.seed(7); obj$dd %>% slice_sample(n = 1500) }, aes(x = hgi),
             inherit.aes = FALSE, sides = "b", alpha = 0.15, length = unit(1.5, "mm")) +
    annotate("text", x = obj$cut90, y = ylo + 0.02, label = "P90",
             colour = "#B2182B", size = 2.0, hjust = -0.15) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.6, size = 2.0,
             label = hdr, fontface = "bold") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 3.4, size = 2.0,
             label = sprintf("P for nonlinearity = %.3f", obj$p_nl)) +
    labs(x = "Hemoglobin glycation index (continuous)", y = "Hazard ratio (95% CI)") +
    coord_cartesian(ylim = c(ylo, max(rcs$hi) * 1.02), clip = "off") +
    theme_classic(base_size = 6.5, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
          plot.margin = margin(3, 3, 2, 2))
}

make_dec <- function(dd) {
  dd <- dd %>% mutate(dec = cut(hgi, breaks = quantile(hgi, probs = 0:10 / 10, na.rm = TRUE),
                                include.lowest = TRUE, labels = 1:10))
  dd %>% group_by(dec) %>%
    summarize(n = n(), ev = sum(incident_stroke6), rate = ev / n * 100,
              lo = binom.test(ev, n)$conf.int[1] * 100,
              hi = binom.test(ev, n)$conf.int[2] * 100, .groups = "drop") %>%
    mutate(dec = as.integer(dec), top = dec == 10)
}

plot_dec <- function(dec, hdr, ymax) {
  ggplot(dec, aes(dec, rate)) +
    geom_col(aes(fill = top), width = 0.72, colour = NA) +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.35, colour = "grey25") +
    scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey72"), guide = "none") +
    scale_x_continuous(breaks = 1:10) +
    coord_cartesian(ylim = c(0, ymax), clip = "off") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.6, size = 2.0,
             label = hdr, fontface = "bold") +
    labs(x = "HGI decile (1 = lowest, 10 = highest)", y = "Incident stroke, %") +
    theme_classic(base_size = 6.5, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
          plot.margin = margin(3, 2, 2, 3))
}

## ---- panels -----------------------------------------------------------------------------
ra <- make_rcs(d, adj3)
write.csv(ra$dat, file.path(OUT, "Fig3a_RCS_source.csv"), row.names = FALSE)
cat(sprintf("panel a: n=%d ev=%d P_nl=%.4f\n", nrow(d), sum(d$incident_stroke6), ra$p_nl))
pa <- plot_rcs(ra, sprintf("All participants (n = %s)", format(nrow(d), big.mark = ",")))

rb <- make_rcs(dn, adj3_nd)
write.csv(rb$dat, file.path(OUT, "Fig3b_RCS_nonDM_source.csv"), row.names = FALSE)
cat(sprintf("panel b: n=%d ev=%d P_nl=%.4f\n", nrow(dn), sum(dn$incident_stroke6), rb$p_nl))
pb <- plot_rcs(rb, sprintf("Without diabetes (n = %s)", format(nrow(dn), big.mark = ",")))

dec_a <- make_dec(d);  write.csv(dec_a, file.path(OUT, "Fig3c_deciles_source.csv"), row.names = FALSE)
dec_b <- make_dec(dn); write.csv(dec_b, file.path(OUT, "Fig3d_deciles_nonDM_source.csv"), row.names = FALSE)
cat(sprintf("panel c top-decile rate=%.1f%%; panel d top-decile rate=%.1f%%\n",
            dec_a$rate[dec_a$dec == 10], dec_b$rate[dec_b$dec == 10]))
ymax <- max(dec_a$hi, dec_b$hi) * 1.08
pc <- plot_dec(dec_a, "All participants", ymax)
pd <- plot_dec(dec_b, "Without diabetes", ymax)

fig3 <- pa + pb + pc + pd + plot_layout(nrow = 2, ncol = 2, byrow = TRUE) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "sans"))

w_mm <- 183; h_mm <- 132
manifest <- file.path(OUT, "Fig3_dose_v3.alignment-layout.json")
write_panel_layout_measured(fig3, manifest_path = manifest,
                            width_in = w_mm / 25.4, height_in = h_mm / 25.4)
qa <- system2("python",
  shQuote(c(file.path(SKILL, "scripts/audit_panel_alignment.py"), manifest,
            "--json-out", file.path(OUT, "Fig3_dose_v3.alignment.json"),
            "--overlay-svg", file.path(OUT, "Fig3_dose_v3.alignment.svg"),
            "--tolerance-pt", "1.5", "--gutter-tolerance-pt", "1.5", "--strict")),
  stdout = TRUE, stderr = TRUE)
cat(paste(qa, collapse = "\n"), "\n")
status <- attr(qa, "status"); if (is.null(status)) status <- 0
if (status != 0) stop("alignment gate failed: ", status)

svglite::svglite(file.path(OUT, "Fig3_dose_v3.svg"), width = w_mm / 25.4, height = h_mm / 25.4)
print(fig3); dev.off()
grDevices::pdf(file.path(OUT, "Fig3_dose_v3.pdf"), width = w_mm / 25.4,
               height = h_mm / 25.4, family = "sans", useDingbats = FALSE)
print(fig3); dev.off()
ragg::agg_tiff(file.path(OUT, "Fig3_dose_v3.tiff"), width = w_mm / 25.4,
               height = h_mm / 25.4, units = "in", res = 600)
print(fig3); dev.off()
cat("DONE Fig3 v3\n")
