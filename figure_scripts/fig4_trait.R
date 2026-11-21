####################################################################################
# Fig 4 | HGI as a trait-like marker (CHARLS two-wave blood evidence)
# Claim: HGI is a within-person stable, trait-like glycation phenotype: its
# level (not its 4-year change) carries stroke risk, and the predictive signal
# wanes with time since measurement (regression dilution).
# a = within-person stability 2011 vs 2015 (representative example);
# b = level vs change vs time window forest (primary evidence + boundary).
# R backend (ggplot2 + patchwork).
####################################################################################
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
SKILL <- "C:/Users/ronaldchen/AppData/Roaming/kimi-desktop/daimon-share/daimon/runtime/kimi-code/home/plugins/managed/nature-skills@nature-skills/skills/nature-figure"
source("D:/CHARLS/HGI_stroke/freeze_20261120/figures/_align_helper.R")
OUT <- "D:/CHARLS/HGI_stroke/freeze_20261120/figures"
F <- "D:/CHARLS/HGI_stroke/freeze_20261120"

dd <- read.csv(file.path(F, "data/charls/analytic_delta_hgi.csv"),
               colClasses = c(ID = "character"))
ok <- complete.cases(dd$hgi, dd$hgi_2015)
r_ab <- cor(dd$hgi[ok], dd$hgi_2015[ok])
c11 <- dd$hgi >= quantile(dd$hgi, 0.9, na.rm = TRUE)
c15 <- dd$hgi_2015 >= quantile(dd$hgi_2015, 0.9, na.rm = TRUE)
stay <- mean(c15[ok & c11]) * 100
q11 <- as.numeric(quantile(dd$hgi, 0.9, na.rm = TRUE))
q15 <- as.numeric(quantile(dd$hgi_2015, 0.9, na.rm = TRUE))

set.seed(11)
ds <- dd[ok, ] %>% slice_sample(n = 4000)
write.csv(ds %>% select(hgi, hgi_2015), file.path(OUT, "Fig4a_stability_source.csv"),
          row.names = FALSE)

pa <- ggplot(ds, aes(hgi, hgi_2015)) +
  geom_vline(xintercept = q11, linetype = "dotted", colour = "grey55", linewidth = 0.35) +
  geom_hline(yintercept = q15, linetype = "dotted", colour = "grey55", linewidth = 0.35) +
  geom_point(size = 0.35, alpha = 0.28, colour = "grey30", stroke = 0) +
  geom_smooth(method = "lm", se = FALSE, colour = "#B2182B", linewidth = 0.6) +
  annotate("text", x = max(ds$hgi) - 0.15, y = min(ds$hgi_2015) + 0.15,
           hjust = 1, vjust = 0,
           size = 2.0, lineheight = 1.15,
           label = sprintf("r = %.2f (4 years apart,\ndifferent assay batches)\n%.0f%% of 2011-high remain\nhigh in 2015 (10%% by chance)",
                           r_ab, stay)) +
  labs(x = "HGI at 2011 (batch A)", y = "HGI at 2015 (batch B)") +
  theme_classic(base_size = 6.5, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        plot.margin = margin(3, 3, 2, 2))

## ---- panel b: level vs change vs time window ---------------------------------------------
fb <- tribble(
  ~lab,                                              ~HR,   ~lo,   ~hi,   ~kind,
  "HGI 2011, outcome 2011-2023 (main)",              1.352, 1.108, 1.650, "level",
  "HGI 2015, outcome 2015-2023 (landmark)",          1.152, 0.920, 1.442, "level",
  "HGI 2011, outcome 2015-2023 (late window)",       1.025, 0.909, 1.157, "level",
  "Delta HGI 2011-2015, per SD",                     0.958, 0.878, 1.046, "change",
  "Delta HGI, top quartile vs rest",                 0.949, 0.781, 1.155, "change")
write.csv(fb, file.path(OUT, "Fig4b_windows_source.csv"), row.names = FALSE)
fb$kind_lab <- ifelse(fb$kind == "level", "HGI level (high vs low)",
                      "4-year change in HGI")
fb$y <- nrow(fb):1
pb <- ggplot(fb, aes(y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_linerange(aes(xmin = lo, xmax = hi, colour = kind_lab), linewidth = 0.6) +
  geom_point(aes(x = HR, colour = kind_lab), size = 1.7) +
  geom_text(aes(x = 2.05, label = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi)),
            hjust = 0, size = 1.9, family = "sans") +
  scale_y_continuous(breaks = fb$y, labels = fb$lab, expand = c(0.06, 0)) +
  scale_colour_manual(values = c("HGI level (high vs low)" = "#B2182B",
                                 "4-year change in HGI" = "grey45"), guide = "none") +
  scale_x_log10(limits = c(0.7, 3.6), breaks = c(0.75, 1, 1.5, 2), expand = c(0, 0)) +
  labs(x = "Hazard ratio (log scale)", y = NULL) +
  theme_classic(base_size = 6.5, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        axis.text.y = element_text(size = 5.8),
        plot.margin = margin(3, 2, 2, 3))

fig4 <- (pa | pb) + plot_layout(widths = c(1, 1.25)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "sans"))

w_mm <- 183; h_mm <- 72
manifest <- file.path(OUT, "Fig4_trait.alignment-layout.json")
write_panel_layout_measured(fig4, manifest_path = manifest,
                            width_in = w_mm / 25.4, height_in = h_mm / 25.4)
qa <- system2("python",
  shQuote(c(file.path(SKILL, "scripts/audit_panel_alignment.py"), manifest,
            "--json-out", file.path(OUT, "Fig4_trait.alignment.json"),
            "--overlay-svg", file.path(OUT, "Fig4_trait.alignment.svg"),
            "--tolerance-pt", "1.5", "--gutter-tolerance-pt", "1.5", "--strict")),
  stdout = TRUE, stderr = TRUE)
cat(paste(qa, collapse = "\n"), "\n")
status <- attr(qa, "status"); if (is.null(status)) status <- 0
if (status != 0) stop("alignment gate failed: ", status)

svglite::svglite(file.path(OUT, "Fig4_trait.svg"), width = w_mm / 25.4, height = h_mm / 25.4)
print(fig4); dev.off()
grDevices::pdf(file.path(OUT, "Fig4_trait.pdf"), width = w_mm / 25.4,
               height = h_mm / 25.4, family = "sans", useDingbats = FALSE)
print(fig4); dev.off()
ragg::agg_tiff(file.path(OUT, "Fig4_trait.tiff"), width = w_mm / 25.4,
               height = h_mm / 25.4, units = "in", res = 600)
print(fig4); dev.off()
cat("DONE Fig4\n")
