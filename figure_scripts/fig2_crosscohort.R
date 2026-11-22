####################################################################################
# Fig 2 | Cross-cohort forest plot: HGI (high vs low) and incident stroke
# Claim: the fully adjusted association is direction- and magnitude-consistent
# across three analytic sets; pooled across independent cohorts it is
# statistically significant with no between-cohort heterogeneity.
# Hero = M3 (full adjustment) + pooled diamonds; M1/M2 = adjustment gradient.
# Single-panel quantitative figure; R backend (ggplot2).
# Output: Fig2_crosscohort.{pdf,svg,tiff} + source data CSV + QA artifacts
####################################################################################
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
SKILL <- "C:/Users/ronaldchen/AppData/Roaming/kimi-desktop/daimon-share/daimon/runtime/kimi-code/home/plugins/managed/nature-skills@nature-skills/skills/nature-figure"
OUT <- "D:/CHARLS/HGI_stroke/freeze_20261120/figures"

## ---- source data (frozen results; see freeze_20261120/results) ---------------------------
rows <- tribble(
  ~grp,          ~model, ~HR,   ~lo,   ~hi,   ~p_lab,
  "CHARLS 2011 main analysis (n = 8,684; 852 events)",     "M1 unadjusted",        1.582, 1.303, 1.920, "P < 0.001",
  "CHARLS 2011 main analysis (n = 8,684; 852 events)",     "M2 demographic",       1.553, 1.279, 1.885, "P < 0.001",
  "CHARLS 2011 main analysis (n = 8,684; 852 events)",     "M3 fully adjusted",    1.352, 1.108, 1.650, "P = 0.003",
  "CHARLS 2015 landmark (n = 6,581; 706 events)",          "M1 unadjusted",        1.494, 1.207, 1.850, "P < 0.001",
  "CHARLS 2015 landmark (n = 6,581; 706 events)",          "M2 demographic",       1.463, 1.181, 1.812, "P < 0.001",
  "CHARLS 2015 landmark (n = 6,581; 706 events)",          "M3 fully adjusted",    1.152, 0.920, 1.442, "P = 0.219",
  "ELSA Wave 2 external (n = 2,577; 176 events)",          "M1 unadjusted",        1.813, 1.187, 2.770, "P = 0.006",
  "ELSA Wave 2 external (n = 2,577; 176 events)",          "M2 demographic",       1.681, 1.098, 2.572, "P = 0.017",
  "ELSA Wave 2 external (n = 2,577; 176 events)",          "M3 fully adjusted",    1.420, 0.903, 2.235, "P = 0.129")
pooled <- tribble(
  ~grp,                                              ~HR,   ~lo,   ~hi,   ~p_lab,
  "Pooled: CHARLS 2011 + ELSA (fixed = random, I2 = 0%)", 1.363, 1.136, 1.635, "P < 0.001",
  "Pooled: all three sets (fixed = random, I2 = 0%)",      1.275, 1.107, 1.469, "P < 0.001",
  "Pooled non-diabetic: all three sets (I2 = 0%)",         1.484, 1.213, 1.816, "P < 0.001")
write.csv(bind_rows(rows %>% mutate(type = "cohort"),
                    pooled %>% mutate(model = "pooled", type = "pooled")),
          file.path(OUT, "Fig2_crosscohort_source.csv"), row.names = FALSE)

## ---- geometry -----------------------------------------------------------------------------
adjv <- c("M1 unadjusted", "M2 demographic", "M3 fully adjusted")
rows <- rows %>% mutate(model = factor(model, levels = adjv),
                        g = rep(1:3, each = 3), within = rep(0:2, 3))
# group tops at 14 / 10.8 / 7.6 ; rows spaced 1.0 within, 1.8 between groups
rows$y <- c(14, 13, 12, 10.3, 9.3, 8.3, 6.6, 5.6, 4.6)
grp_pos <- data.frame(grp = unique(rows$grp), yt = c(14.7, 11.0, 7.3))
pooled$y <- c(3.4, 2.4, 1.4)
pool_grp_y <- 4.1

col_m <- c("M1 unadjusted" = "grey72", "M2 demographic" = "grey42",
           "M3 fully adjusted" = "#B2182B")

fig2 <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_linerange(data = rows, aes(xmin = lo, xmax = hi, y = y, colour = model),
                 linewidth = 0.6) +
  geom_point(data = rows, aes(x = HR, y = y, colour = model, size = model)) +
  geom_polygon(data = pooled %>% mutate(id = row_number()) %>% rowwise() %>%
                 do(data.frame(id = .$id,
                               x = c(.$lo, .$HR, .$hi, .$HR),
                               y = .$y + c(0, 0.42, 0, -0.42))) %>% ungroup(),
               aes(x = x, y = y, group = id), fill = "#2166AC", colour = NA) +
  geom_text(data = rows, aes(x = 3.35, y = y, label = sprintf("%.2f (%.2f-%.2f)  %s", HR, lo, hi, p_lab)),
            hjust = 0, size = 1.9, family = "sans") +
  geom_text(data = pooled, aes(x = 3.35, y = y, label = sprintf("%.2f (%.2f-%.2f)  %s", HR, lo, hi, p_lab)),
            hjust = 0, size = 1.9, family = "sans", fontface = "bold", colour = "#2166AC") +
  geom_text(data = grp_pos, aes(x = 0.42, y = yt, label = grp), hjust = 0,
            size = 2.1, fontface = "bold", family = "sans") +
  geom_text(data = pooled %>%
              mutate(short = c("CHARLS 2011 + ELSA", "All three sets",
                               "Non-diabetic, all three sets")),
            aes(x = 0.42, y = y, label = short), hjust = 0, size = 1.8,
            colour = "#2166AC", family = "sans") +
  annotate("text", x = 0.42, y = pool_grp_y, hjust = 0, size = 2.1,
           fontface = "bold", label = "Cross-cohort pooling") +
  annotate("text", x = 3.35, y = 15.7, hjust = 0, size = 2.0,
           fontface = "bold", label = "HR (95% CI)") +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3), limits = c(0.42, 6.4),
                expand = c(0, 0)) +
  scale_colour_manual(values = col_m, guide = "none") +
  scale_size_manual(values = c("M1 unadjusted" = 1.1, "M2 demographic" = 1.3,
                               "M3 fully adjusted" = 1.8), guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.5, 16.6), clip = "off") +
  labs(x = "Hazard ratio, high HGI vs low HGI (log scale)", y = NULL) +
  theme_classic(base_size = 6.5, base_family = "sans") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 5.8),
        axis.title.x = element_text(size = 6.2),
        plot.margin = margin(3, 4, 3, 3))

# model key (direct labels; own row above the HR column header to avoid overlap)
fig2 <- fig2 +
  annotate("point", x = c(1.6, 2.05, 2.55), y = 16.0,
           colour = unname(col_m), size = c(1.1, 1.3, 1.8)) +
  annotate("text", x = c(1.7, 2.17, 2.7), y = 16.0, hjust = 0, size = 1.9,
           label = c("M1", "M2", "M3 fully adjusted")) +
  annotate("text", x = 1.08, y = 16.0, hjust = 0, size = 1.9, label = "Adjustment:")

w_mm <- 183; h_mm <- 102
svglite::svglite(file.path(OUT, "Fig2_crosscohort.svg"), width = w_mm / 25.4, height = h_mm / 25.4)
print(fig2); dev.off()
grDevices::pdf(file.path(OUT, "Fig2_crosscohort.pdf"), width = w_mm / 25.4,
               height = h_mm / 25.4, family = "sans", useDingbats = FALSE)
print(fig2); dev.off()
ragg::agg_tiff(file.path(OUT, "Fig2_crosscohort.tiff"), width = w_mm / 25.4,
               height = h_mm / 25.4, units = "in", res = 600)
print(fig2); dev.off()
cat("DONE Fig2\n")
