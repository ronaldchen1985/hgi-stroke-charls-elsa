####################################################################################
# Fig 3 v2 | Dose-response shape of the HGI-stroke association (CHARLS main set)
# a = fully adjusted RCS, full main set; b = fully adjusted RCS, participants
# without diabetes (co-primary layer); c = observed stroke incidence by HGI decile.
# R backend (ggplot2 + patchwork). Adds panel b to v1 (2026-11-22).
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

## helper: fit RCS (ns df=4) and return plotting data + nonlinearity P
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

## ---- panel a: full main set -------------------------------------------------------------
ra <- make_rcs(d, adj3)
write.csv(ra$dat, file.path(OUT, "Fig3a_RCS_source.csv"), row.names = FALSE)
cat(sprintf("panel a: n=%d ev=%d P_nl=%.4f cut90=%.3f\n",
            nrow(d), sum(d$incident_stroke6), ra$p_nl, ra$cut90))
pa <- plot_rcs(ra, sprintf("All participants (n = %s)", format(nrow(d), big.mark = ",")))

## ---- panel b: participants without diabetes (broad definition) -------------------------
dn <- d %>% filter(diabetes == 0)
stopifnot(nrow(dn) == 7283, sum(dn$incident_stroke6) == 681)
adj3_nd <- "age + sex + urban + edu + htn + dyslip + heart + sbp + crp + tg + hdl"
rb <- make_rcs(dn, adj3_nd)
write.csv(rb$dat, file.path(OUT, "Fig3b_RCS_nonDM_source.csv"), row.names = FALSE)
cat(sprintf("panel b: n=%d ev=%d P_nl=%.4f cut90=%.3f\n",
            nrow(dn), sum(dn$incident_stroke6), rb$p_nl, rb$cut90))
pb <- plot_rcs(rb, sprintf("Without diabetes (n = %s)", format(nrow(dn), big.mark = ",")))

## ---- panel c: event rate by HGI decile (full main set) ----------------------------------
d <- d %>% mutate(dec = cut(hgi, breaks = quantile(hgi, probs = 0:10 / 10, na.rm = TRUE),
                            include.lowest = TRUE, labels = 1:10))
dec <- d %>% group_by(dec) %>%
  summarize(n = n(), ev = sum(incident_stroke6), rate = ev / n * 100,
            lo = binom.test(ev, n)$conf.int[1] * 100,
            hi = binom.test(ev, n)$conf.int[2] * 100, .groups = "drop") %>%
  mutate(dec = as.integer(dec), top = dec == 10)
write.csv(dec, file.path(OUT, "Fig3c_deciles_source.csv"), row.names = FALSE)

pc <- ggplot(dec, aes(dec, rate)) +
  geom_col(aes(fill = top), width = 0.72, colour = NA) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.35, colour = "grey25") +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey72"), guide = "none") +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "HGI decile (1 = lowest, 10 = highest)", y = "Incident stroke, %") +
  theme_classic(base_size = 6.5, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        plot.margin = margin(3, 2, 2, 3))

fig3 <- (pa | pb | pc) + plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "sans"))

w_mm <- 183; h_mm <- 72
manifest <- file.path(OUT, "Fig3_dose_v2.alignment-layout.json")
write_panel_layout_measured(fig3, manifest_path = manifest,
                            width_in = w_mm / 25.4, height_in = h_mm / 25.4)
qa <- system2("python",
  shQuote(c(file.path(SKILL, "scripts/audit_panel_alignment.py"), manifest,
            "--json-out", file.path(OUT, "Fig3_dose_v2.alignment.json"),
            "--overlay-svg", file.path(OUT, "Fig3_dose_v2.alignment.svg"),
            "--tolerance-pt", "1.5", "--gutter-tolerance-pt", "1.5", "--strict")),
  stdout = TRUE, stderr = TRUE)
cat(paste(qa, collapse = "\n"), "\n")
status <- attr(qa, "status"); if (is.null(status)) status <- 0
if (status != 0) stop("alignment gate failed: ", status)

svglite::svglite(file.path(OUT, "Fig3_dose_v2.svg"), width = w_mm / 25.4, height = h_mm / 25.4)
print(fig3); dev.off()
grDevices::pdf(file.path(OUT, "Fig3_dose_v2.pdf"), width = w_mm / 25.4,
               height = h_mm / 25.4, family = "sans", useDingbats = FALSE)
print(fig3); dev.off()
ragg::agg_tiff(file.path(OUT, "Fig3_dose_v2.tiff"), width = w_mm / 25.4,
               height = h_mm / 25.4, units = "in", res = 600)
print(fig3); dev.off()
cat("DONE Fig3 v2\n")
