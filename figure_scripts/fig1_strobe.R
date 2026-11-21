####################################################################################
# Fig 1 | Multi-panel STROBE flow: three analytic cohorts
# Claim: discovery (CHARLS 2011), internal landmark (CHARLS 2015) and external
# validation (ELSA W2) cohorts were assembled with harmonized eligibility logic.
# Panels a/b/c = setup role; quantitative grid; R backend (ggplot2 + patchwork).
# Output: freeze_20261120/figures/Fig1_STROBE.{pdf,svg,tiff} + QA artifacts
####################################################################################
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(grid) })
SKILL <- "C:/Users/ronaldchen/AppData/Roaming/kimi-desktop/daimon-share/daimon/runtime/kimi-code/home/plugins/managed/nature-skills@nature-skills/skills/nature-figure"
source("D:/CHARLS/HGI_stroke/freeze_20261120/figures/_align_helper.R")
OUT <- "D:/CHARLS/HGI_stroke/freeze_20261120/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- flow-column builder: boxes + one text run per line (tight bboxes) ------------------
flow_col <- function(tag, subtitle, stages) {
  # stages: data.frame(l1, l2, n, excl1, excl2)
  k <- nrow(stages)
  pitch <- 1.18; hh <- 0.44                       # slot pitch / box half-height
  y <- (k:1) * pitch
  d <- stages; d$y <- y
  # per-line label runs inside each box
  lines <- do.call(rbind, lapply(seq_len(k), function(i) {
    data.frame(x = 5, y = y[i] + c(0.27, 0, -0.27),
               txt = c(d$l1[i], d$l2[i], paste0("(n = ", d$n[i], ")")))
  }))
  lines <- lines[nzchar(lines$txt), ]
  # per-line exclusion runs in gaps
  ex <- do.call(rbind, lapply(2:k, function(i) {
    gy <- y[i - 1] - pitch / 2
    t1 <- d$excl1[i]; t2 <- d$excl2[i]
    if (!nzchar(t1)) return(NULL)
    if (nzchar(t2)) data.frame(x = 10.45, y = gy + c(0.13, -0.13), txt = c(t1, t2))
    else            data.frame(x = 10.45, y = gy, txt = t1)
  }))
  ggplot(d) +
    geom_rect(aes(xmin = 0, xmax = 10, ymin = y - hh, ymax = y + hh),
              fill = "grey96", colour = "grey30", linewidth = 0.35) +
    geom_text(data = lines, aes(x = x, y = y, label = txt),
              size = 2.0, family = "sans") +
    geom_segment(data = data.frame(y0 = y[-k] - hh, y1 = y[-k] - pitch + hh),
                 aes(x = 5, xend = 5, y = y0, yend = y1),
                 arrow = arrow(length = unit(1.6, "mm"), type = "closed"),
                 linewidth = 0.4, colour = "grey20") +
    geom_text(data = ex, aes(x = x, y = y, label = txt), hjust = 0,
              size = 1.8, colour = "grey45", family = "sans") +
    annotate("text", x = 0, y = k * pitch + hh + 0.55, label = tag, hjust = 0,
             size = 2.9, fontface = "bold", family = "sans") +
    annotate("text", x = 0.35, y = k * pitch + hh + 0.22, label = subtitle, hjust = 0,
             size = 2.0, colour = "grey25", family = "sans") +
    coord_cartesian(xlim = c(0, 15.4), ylim = c(0.3, k * pitch + hh + 0.85), clip = "off") +
    theme_void(base_family = "sans") +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(2, 1, 2, 1))
}

a <- flow_col("a", "CHARLS 2011 main analysis (discovery)",
  data.frame(
    l1 = c("FPG and HbA1c measured", "Age >= 45 years", "Free of stroke",
           ">= 1 follow-up interview", "Complete covariates", "Incident stroke events"),
    l2 = c("at Wave 1 (2011)", "", "at baseline", "during 2013-2023",
           "(main analysis set)", "during 2011-2023"),
    n = c("11,532", "11,305", "11,038", "10,643", "8,684", "852"),
    excl1 = c("", "-227 age < 45", "-267 baseline stroke", "-395 lost",
              "-1,959 missing", ""),
    excl2 = c("", "", "", "to follow-up", "covariates", "")))

b <- flow_col("b", "CHARLS 2015 landmark (internal validation)",
  data.frame(
    l1 = c("Blood drawn at both", "Free of stroke", ">= 1 follow-up interview",
           "Incident stroke events"),
    l2 = c("2011 and 2015", "before 2015", "2018-2023, complete covariates",
           "during 2015-2023"),
    n = c("7,166", "6,730", "6,581", "706"),
    excl1 = c("", "-436 stroke", "-149 missing", ""),
    excl2 = c("", "before 2015", "covariates (2015)", "")))

c <- flow_col("c", "ELSA Wave 2 (external validation)",
  data.frame(
    l1 = c("Harmonized ELSA", "Fasting FPG and HbA1c", "Free of stroke",
           ">= 1 follow-up interview", "Complete covariates", "Incident stroke events"),
    l2 = c("participants", "at Wave 2 (2004/05), age >= 50", "at baseline",
           "during Waves 3-11", "(main analysis set)", "during 2006-2023"),
    n = c("21,679", "3,428", "3,349", "3,062", "2,577", "176"),
    excl1 = c("", "-18,251 no fasting", "-79 baseline stroke", "-287 lost",
              "-485 missing", ""),
    excl2 = c("", "blood / age < 50", "", "to follow-up", "covariates", "")))

fig1 <- (a | b | c) + plot_layout(widths = c(1, 1, 1))

w_mm <- 183; h_mm <- 118
manifest <- file.path(OUT, "Fig1_STROBE.alignment-layout.json")
write_panel_layout_measured(fig1, manifest_path = manifest,
                            width_in = w_mm / 25.4, height_in = h_mm / 25.4)
qa <- system2("python",
  shQuote(c(file.path(SKILL, "scripts/audit_panel_alignment.py"), manifest,
            "--json-out", file.path(OUT, "Fig1_STROBE.alignment.json"),
            "--overlay-svg", file.path(OUT, "Fig1_STROBE.alignment.svg"),
            "--tolerance-pt", "1.5", "--gutter-tolerance-pt", "1.5", "--strict")),
  stdout = TRUE, stderr = TRUE)
cat(paste(qa, collapse = "\n"), "\n")
status <- attr(qa, "status"); if (is.null(status)) status <- 0
if (status != 0) stop("Panel alignment gate failed with exit code ", status)

svglite::svglite(file.path(OUT, "Fig1_STROBE.svg"), width = w_mm / 25.4, height = h_mm / 25.4)
print(fig1); dev.off()
grDevices::pdf(file.path(OUT, "Fig1_STROBE.pdf"), width = w_mm / 25.4,
               height = h_mm / 25.4, family = "sans", useDingbats = FALSE)
print(fig1); dev.off()
ragg::agg_tiff(file.path(OUT, "Fig1_STROBE.tiff"), width = w_mm / 25.4,
               height = h_mm / 25.4, units = "in", res = 600)
print(fig1); dev.off()
cat("DONE Fig1\n")
