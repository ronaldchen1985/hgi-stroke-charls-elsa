####################################################################################
# Fig 5 | Subgroup consistency and robustness battery (CHARLS main set)
# Claim: the association is consistent across clinically meaningful subgroups
# (no interaction) and survives a battery of robustness checks.
# a = subgroup forest with interaction P; b = sensitivity-analysis forest
# against the main M3 estimate. R backend (ggplot2 + patchwork).
####################################################################################
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
SKILL <- "C:/Users/ronaldchen/AppData/Roaming/kimi-desktop/daimon-share/daimon/runtime/kimi-code/home/plugins/managed/nature-skills@nature-skills/skills/nature-figure"
source("D:/CHARLS/HGI_stroke/freeze_20261120/figures/_align_helper.R")
OUT <- "D:/CHARLS/HGI_stroke/freeze_20261120/figures"

## ---- panel a: subgroups (frozen Table 3 values) ------------------------------------------
sg <- tribble(
  ~factor, ~level, ~nev, ~HR, ~lo, ~hi, ~p_int,
  "Diabetes",      "No",        "7,283 / 681", 1.475, 1.147, 1.896, "0.40",
  "Diabetes",      "Yes",       "1,401 / 171", 1.223, 0.876, 1.708, "0.40",
  "Sex",           "Male",      "4,017 / 401", 1.290, 0.952, 1.749, "0.77",
  "Sex",           "Female",    "4,667 / 451", 1.368, 1.050, 1.783, "0.77",
  "Age",           "< 65 y",    "6,326 / 582", 1.295, 1.016, 1.652, "0.99",
  "Age",           ">= 65 y",   "2,358 / 270", 1.466, 1.036, 2.073, "0.99",
  "Residence",     "Rural",     "5,684 / 562", 1.502, 1.176, 1.917, "0.21",
  "Residence",     "Urban",     "3,000 / 290", 1.151, 0.816, 1.624, "0.21")
write.csv(sg, file.path(OUT, "Fig5a_subgroups_source.csv"), row.names = FALSE)
sg$y <- c(9, 8, 6.8, 5.8, 4.6, 3.6, 2.4, 1.4)
sg_hdr <- data.frame(factor = unique(sg$factor), y = c(9.75, 7.55, 5.35, 3.15))
sg_int <- data.frame(p_int = c("0.40", "0.77", "0.99", "0.21"), y = c(8.5, 6.3, 4.1, 1.9))
ylims <- c(0.9, 10.6)

## single-panel forest: all text columns live in the x-range beyond the
## widest CI (max 2.07), so no stroke can cross them
x_lab <- 0.53; x_nev <- 3.0; x_hr <- 5.05; x_int <- 10.2
pa <- ggplot(sg) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_linerange(aes(xmin = lo, xmax = hi, y = y), linewidth = 0.55, colour = "grey30") +
  geom_point(aes(x = HR, y = y), size = 1.6, colour = "#B2182B") +
  geom_text(aes(x = x_lab, y = y, label = level), hjust = 0, size = 1.9, family = "sans") +
  geom_text(aes(x = x_nev, y = y, label = nev), hjust = 0, size = 1.7,
            colour = "grey40", family = "sans") +
  geom_text(aes(x = x_hr, y = y, label = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi)),
            hjust = 0, size = 1.7, family = "sans") +
  geom_text(data = sg_hdr, aes(x = x_lab, y = y, label = factor), hjust = 0,
            size = 2.0, fontface = "bold", family = "sans") +
  geom_text(data = sg_int, aes(x = x_int, y = y, label = p_int), hjust = 0,
            size = 1.7, colour = "grey40", family = "sans") +
  annotate("text", x = x_nev, y = 10.4, label = "n / events", hjust = 0,
           size = 1.9, fontface = "bold") +
  annotate("text", x = x_hr, y = 10.4, label = "HR (95% CI)", hjust = 0,
           size = 1.9, fontface = "bold") +
  annotate("text", x = x_int, y = 10.52, label = "P for", hjust = 0,
           size = 1.9, fontface = "bold") +
  annotate("text", x = x_int, y = 10.22, label = "interaction", hjust = 0,
           size = 1.9, fontface = "bold") +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 2), limits = c(0.5, 16)) +
  coord_cartesian(ylim = ylims, clip = "off") +
  labs(x = "Hazard ratio (log scale)", y = NULL) +
  theme_classic(base_size = 6.5, base_family = "sans") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.line.x = element_line(linewidth = 0.35),
        axis.ticks.x = element_line(linewidth = 0.35),
        plot.margin = margin(3, 2, 2, 2))

## ---- panel b: sensitivity battery (frozen Table 4 / S7 / S8 / E-value) -------------------
sa <- tribble(
  ~lab,                                                        ~HR,   ~lo,   ~hi,   ~ref,
  "Main analysis (M3, fully adjusted)",                        1.352, 1.108, 1.650, TRUE,
  "Fine-Gray competing risk (death)",                          1.345, 1.102, 1.640, FALSE,
  "S1  HGI derived in non-diabetic subset",                    1.275, 1.039, 1.563, FALSE,
  "S2  Published (NHANES) HGI equation",                       1.362, 1.117, 1.660, FALSE,
  "S3  Excluding anaemia (low-HGI arm 0.96, P = 0.73)",        1.280, 1.036, 1.582, FALSE,
  "S4a Alternative follow-up rule (report wave)",              1.333, 1.092, 1.627, FALSE,
  "S4b Diagnosis-year refined follow-up",                      1.352, 1.108, 1.650, FALSE,
  "S5a Without SBP covariate",                                 1.346, 1.118, 1.620, FALSE,
  "S5b Multiple imputation (m = 20)",                          1.306, 1.089, 1.567, FALSE,
  "S7  Mutual adjustment for stress hyperglycaemia ratio",     1.397, 1.129, 1.729, FALSE,
  "S8  Interval-censoring model (Weibull)",                    1.352, 1.108, 1.650, FALSE)
write.csv(sa, file.path(OUT, "Fig5b_sensitivity_source.csv"), row.names = FALSE)
sa$y <- nrow(sa):1

pb <- ggplot(sa) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = 1.352, linetype = "dotted", colour = "#B2182B", linewidth = 0.3) +
  geom_linerange(aes(xmin = lo, xmax = hi, y = y), linewidth = 0.55, colour = "grey30") +
  geom_point(aes(x = HR, y = y, colour = ref), size = 1.6) +
  geom_text(aes(x = 2.1, y = y, label = sprintf("%.3f (%.3f-%.3f)", HR, lo, hi)),
            hjust = 0, size = 1.8, family = "sans") +
  scale_y_continuous(breaks = sa$y, labels = sa$lab, expand = c(0.04, 0)) +
  scale_colour_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey35"), guide = "none") +
  scale_x_log10(breaks = c(1, 1.25, 1.5, 2), limits = c(0.95, 3.6)) +
  labs(x = "Hazard ratio (log scale)", y = NULL) +
  theme_classic(base_size = 6.5, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        axis.text.y = element_text(size = 5.6),
        plot.margin = margin(3, 2, 2, 3))

fig5 <- (pa | pb) + plot_layout(widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "sans"))

w_mm <- 183; h_mm <- 92
manifest <- file.path(OUT, "Fig5_robust.alignment-layout.json")
write_panel_layout_measured(fig5, manifest_path = manifest,
                            width_in = w_mm / 25.4, height_in = h_mm / 25.4)
qa <- system2("python",
  shQuote(c(file.path(SKILL, "scripts/audit_panel_alignment.py"), manifest,
            "--json-out", file.path(OUT, "Fig5_robust.alignment.json"),
            "--overlay-svg", file.path(OUT, "Fig5_robust.alignment.svg"),
            "--tolerance-pt", "1.5", "--gutter-tolerance-pt", "1.5", "--strict")),
  stdout = TRUE, stderr = TRUE)
cat(paste(qa, collapse = "\n"), "\n")
status <- attr(qa, "status"); if (is.null(status)) status <- 0
if (status != 0) stop("alignment gate failed: ", status)

svglite::svglite(file.path(OUT, "Fig5_robust.svg"), width = w_mm / 25.4, height = h_mm / 25.4)
print(fig5); dev.off()
grDevices::pdf(file.path(OUT, "Fig5_robust.pdf"), width = w_mm / 25.4,
               height = h_mm / 25.4, family = "sans", useDingbats = FALSE)
print(fig5); dev.off()
ragg::agg_tiff(file.path(OUT, "Fig5_robust.tiff"), width = w_mm / 25.4,
               height = h_mm / 25.4, units = "in", res = 600)
print(fig5); dev.off()
cat("DONE Fig5\n")
