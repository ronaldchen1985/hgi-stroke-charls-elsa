####################################################################################
# Fig 3 | Dose-response shape of the HGI-stroke association (CHARLS main set)
# Claim: the association is concentrated in the top tail of HGI (threshold/J
# pattern), not a linear gradient. a = fully adjusted RCS (hero);
# b = observed stroke incidence by HGI decile. R backend (ggplot2 + patchwork).
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
adj3 <- "age + sex + urban + edu + htn + dm + dyslip + heart + sbp + crp + tg + hdl"

## ---- panel a: RCS spline (ns df = 4, consistent with the frozen LRT) -------------------
B <- ns(d$hgi, df = 4)
fit <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ B +", adj3)), data = d)
ref <- median(d$hgi)
grid_x <- seq(quantile(d$hgi, 0.01), quantile(d$hgi, 0.99), length.out = 200)
B_ref <- predict(B, ref)
B_new <- predict(B, grid_x)
bb <- coef(fit)[paste0("B", 1:4)]
V  <- vcov(fit)[paste0("B", 1:4), paste0("B", 1:4)]
XB <- sweep(B_new, 2, rep(1, 4), "*") - sweep(matrix(B_ref, nrow(B_new), 4, byrow = TRUE),
                                              2, rep(1, 4), "*")
lp <- as.vector(XB %*% bb)
se <- sqrt(rowSums((XB %*% V) * XB))
rcs <- data.frame(hgi = grid_x, HR = exp(lp), lo = exp(lp - 1.96 * se),
                  hi = exp(lp + 1.96 * se))
write.csv(rcs, file.path(OUT, "Fig3a_RCS_source.csv"), row.names = FALSE)
# LRT for nonlinearity (for annotation)
m_lin <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ hgi +", adj3)), data = d)
m_rcs <- coxph(as.formula(paste("Surv(fu_v2, incident_stroke6) ~ B +", adj3)), data = d)
p_nl <- anova(m_lin, m_rcs)[2, "Pr(>|Chi|)"]
cut90 <- as.numeric(quantile(d$hgi, 0.9))

pa <- ggplot(rcs, aes(hgi, HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = cut90, linetype = "dotted", colour = "#B2182B", linewidth = 0.35) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#B2182B", alpha = 0.14) +
  geom_line(colour = "#B2182B", linewidth = 0.7) +
  geom_rug(data = { set.seed(7); d %>% slice_sample(n = 1500) }, aes(x = hgi),
           inherit.aes = FALSE, sides = "b", alpha = 0.15, length = unit(1.5, "mm")) +
  annotate("text", x = cut90, y = 0.72, label = "P90",
           colour = "#B2182B", size = 2.0, hjust = -0.15) +
  annotate("text", x = min(d$hgi) + 0.05, y = max(rcs$hi) * 0.97,
           label = sprintf("P for nonlinearity = %.3f", p_nl), hjust = 0, size = 2.0) +
  labs(x = "Hemoglobin glycation index (continuous)", y = "Hazard ratio (95% CI)") +
  coord_cartesian(ylim = c(0.7, max(rcs$hi) * 1.02)) +
  theme_classic(base_size = 6.5, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        plot.margin = margin(3, 3, 2, 2))

## ---- panel b: event rate by HGI decile --------------------------------------------------
d <- d %>% mutate(dec = cut(hgi, breaks = quantile(hgi, probs = 0:10 / 10, na.rm = TRUE),
                            include.lowest = TRUE, labels = 1:10))
dec <- d %>% group_by(dec) %>%
  summarize(n = n(), ev = sum(incident_stroke6), rate = ev / n * 100,
            lo = binom.test(ev, n)$conf.int[1] * 100,
            hi = binom.test(ev, n)$conf.int[2] * 100, .groups = "drop") %>%
  mutate(dec = as.integer(dec), top = dec == 10)
write.csv(dec, file.path(OUT, "Fig3b_deciles_source.csv"), row.names = FALSE)

pb <- ggplot(dec, aes(dec, rate)) +
  geom_col(aes(fill = top), width = 0.72, colour = NA) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.35, colour = "grey25") +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey72"), guide = "none") +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "HGI decile (1 = lowest, 10 = highest)", y = "Incident stroke, %") +
  theme_classic(base_size = 6.5, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        plot.margin = margin(3, 2, 2, 3))

fig3 <- (pa | pb) + plot_layout(widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "sans"))

w_mm <- 183; h_mm <- 78
manifest <- file.path(OUT, "Fig3_dose.alignment-layout.json")
write_panel_layout_measured(fig3, manifest_path = manifest,
                            width_in = w_mm / 25.4, height_in = h_mm / 25.4)
qa <- system2("python",
  shQuote(c(file.path(SKILL, "scripts/audit_panel_alignment.py"), manifest,
            "--json-out", file.path(OUT, "Fig3_dose.alignment.json"),
            "--overlay-svg", file.path(OUT, "Fig3_dose.alignment.svg"),
            "--tolerance-pt", "1.5", "--gutter-tolerance-pt", "1.5", "--strict")),
  stdout = TRUE, stderr = TRUE)
cat(paste(qa, collapse = "\n"), "\n")
status <- attr(qa, "status"); if (is.null(status)) status <- 0
if (status != 0) stop("alignment gate failed: ", status)

svglite::svglite(file.path(OUT, "Fig3_dose.svg"), width = w_mm / 25.4, height = h_mm / 25.4)
print(fig3); dev.off()
grDevices::pdf(file.path(OUT, "Fig3_dose.pdf"), width = w_mm / 25.4,
               height = h_mm / 25.4, family = "sans", useDingbats = FALSE)
print(fig3); dev.off()
ragg::agg_tiff(file.path(OUT, "Fig3_dose.tiff"), width = w_mm / 25.4,
               height = h_mm / 25.4, units = "in", res = 600)
print(fig3); dev.off()
cat("DONE Fig3\n")
