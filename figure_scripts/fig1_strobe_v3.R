####################################################################################
# Fig 1 (v3 redesign) | Three-panel STROBE participant flow
# CHARLS 2011 discovery / CHARLS 2015 landmark / ELSA W2 external evaluation
# v3 changes vs v2: larger exclusion labels (>=6pt), larger box text,
# tighter vertical rhythm, exclusion labels closer to arrows, darker grey
# v3.1: exclusion notes are now boxed (light fill + grey border), font 2.35mm,
# arrow connectors point to box edge; long notes re-wrapped to fit the box
# v3.2: panel tag chips removed -> plain bold letters, consistent with Figs 2-4
# Output: Fig1_STROBE_v3.tiff (600 dpi, 183 mm) + pdf + alignment JSON
####################################################################################
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(grid) })
OUT <- "D:/CHARLS/HGI_stroke/freeze_20261120/figures"
SKILL_SCRIPTS <- "C:/Users/ronaldchen/AppData/Roaming/kimi-desktop/daimon-share/daimon/runtime/kimi-code/home/plugins/managed/nature-skills@nature-skills/skills/nature-figure/scripts"

flow_col <- function(tag, subtitle, stages, accent, fill, n_max) {
  k <- nrow(stages)
  pitch <- 1.12; hh <- 0.46
  y <- (k:1) * pitch + (n_max - k) * pitch      # bottom-align columns
  d <- stages; d$y <- y
  is_event <- seq_len(k) == k

  # exclusion notes: boxed labels in the gaps, arrow connector from stage box
  ex_list <- lapply(2:k, function(i) {
    if (!nzchar(d$excl[i])) return(NULL)
    ls <- strsplit(d$excl[i], "\n", fixed = TRUE)[[1]]
    m <- length(ls)
    yc <- y[i - 1] - pitch / 2
    bh <- 0.20 + 0.12 * (m - 1)              # half-height: line block + padding
    list(box = data.frame(xmin = 10.05, xmax = 17.15, ymin = yc - bh, ymax = yc + bh),
         txt = data.frame(x = 13.6,
                          y = yc + ((m - 1) / 2 - seq(0, m - 1)) * 0.24, txt = ls),
         arr = data.frame(x0 = 9.45, x1 = 9.98, y = yc))
  })
  ex_list <- ex_list[!vapply(ex_list, is.null, logical(1))]
  exb <- if (length(ex_list)) do.call(rbind, lapply(ex_list, `[[`, "box")) else NULL
  ext <- if (length(ex_list)) do.call(rbind, lapply(ex_list, `[[`, "txt")) else NULL
  exa <- if (length(ex_list)) do.call(rbind, lapply(ex_list, `[[`, "arr")) else NULL

  box_fill <- ifelse(is_event, accent, fill)
  box_col  <- ifelse(is_event, accent, "grey45")
  txt_col  <- ifelse(is_event, "white", "grey10")

  label_df <- do.call(rbind, lapply(seq_len(k), function(i) {
    l1_lines <- strsplit(d$l1[i], "\n", fixed = TRUE)[[1]]
    all_lines <- c(l1_lines, sprintf("(n = %s)", formatC(d$n[i], format = "f", big.mark = ",", digits = 0)))
    m <- length(all_lines)
    yy <- y[i] + ((m - 1) / 2 - seq(0, m - 1)) * 0.22
    is_n <- seq_len(m) == m
    data.frame(x = 4.7, y = yy, txt = all_lines,
               font = ifelse(is_n, 2, 1), sz = ifelse(is_n, 2.35, 2.2),
               col = txt_col[i], stringsAsFactors = FALSE)
  }))

  ggplot(d) +
    geom_rect(aes(xmin = 0, xmax = 9.5, ymin = y - hh, ymax = y + hh),
              fill = box_fill, colour = box_col, linewidth = 0.5) +
    # 每行一个独立 layer（独立 textGrob），避免 PDF 文本行框互相吞并
    lapply(seq_len(nrow(label_df)), function(i) {
      annotate("text", x = label_df$x[i], y = label_df$y[i], label = label_df$txt[i],
               fontface = label_df$font[i], colour = label_df$col[i],
               size = label_df$sz[i], family = "sans", lineheight = 1.0)
    }) +
    geom_segment(data = data.frame(y0 = y[-k] - hh, y1 = y[-k] - pitch + hh),
                 aes(x = 4.7, xend = 4.7, y = y0, yend = y1),
                 arrow = arrow(length = unit(2.0, "mm"), type = "closed"),
                 linewidth = 0.55, colour = "grey30") +
    { if (!is.null(exb)) geom_rect(data = exb,
                                   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                   fill = "#F6F4F1", colour = "grey55", linewidth = 0.4) } +
    { if (!is.null(exa)) geom_segment(data = exa, aes(x = x0, xend = x1, y = y, yend = y),
                                      arrow = arrow(length = unit(1.7, "mm"), type = "closed"),
                                      linewidth = 0.4, colour = "grey45") } +
    { if (!is.null(ext)) lapply(seq_len(nrow(ext)), function(i) {
        annotate("text", x = ext$x[i], y = ext$y[i], label = ext$txt[i],
                 size = 2.35, colour = "grey20", family = "sans")
      }) } +
    # panel tag (plain bold, consistent with other figures) + subtitle
    annotate("text", x = 0, y = n_max * pitch + hh + 0.73, label = tag, hjust = 0,
             colour = "grey10", fontface = "bold", size = 3.2, family = "sans") +
    annotate("text", x = 1.05, y = n_max * pitch + hh + 0.73, label = subtitle, hjust = 0,
             size = 2.5, colour = "grey15", fontface = "bold", family = "sans") +
    coord_cartesian(xlim = c(0, 17.4), ylim = c(0.3, n_max * pitch + hh + 1.2), clip = "off") +
    theme_void(base_family = "sans") +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(2, 1, 2, 1))
}

a <- flow_col("a", "CHARLS 2011 main analysis (discovery)",
  data.frame(
    l1 = c("FPG and HbA1c measured\nat Wave 1 (2011)", "Age \u2265 45 years",
           "Free of stroke at baseline", "\u2265 1 follow-up interview\nduring 2013\u20132023",
           "Complete covariates\n(main analysis set)", "Incident stroke\nduring 2011\u20132023"),
    n = c(11532, 11305, 11038, 10651, 8684, 852),
    excl = c("", "\u2212227 age < 45", "\u2212267 baseline\nstroke", "\u2212387 lost to\nfollow-up",
             "\u22128 missing sex\n\u22121,959 missing\ncovariates", ""),
    stringsAsFactors = FALSE),
  accent = "#1F6FA8", fill = "#E9F1F8", n_max = 6)

b <- flow_col("b", "CHARLS 2015 landmark (internal validation)",
  data.frame(
    l1 = c("Blood drawn at both\n2011 and 2015", "Free of stroke before 2015",
           "\u2265 1 interview 2018\u20132023,\ncomplete covariates (2015)",
           "Incident stroke\nduring 2015\u20132023"),
    n = c(7166, 6730, 6581, 706),
    excl = c("", "\u2212436 stroke\nbefore 2015", "\u2212149 missing\ncovariates", ""),
    stringsAsFactors = FALSE),
  accent = "#2E8B74", fill = "#EAF5F0", n_max = 6)

c <- flow_col("c", "ELSA Wave 2 (external validation)",
  data.frame(
    l1 = c("Harmonized ELSA\nparticipants",
           "Fasting FPG and HbA1c\nat Wave 2 (2004/05),\nage \u2265 50 years",
           "Free of stroke at baseline", "\u2265 1 follow-up interview\nduring Waves 3\u201311",
           "Complete covariates\n(main analysis set)", "Incident stroke\nduring 2006\u20132023"),
    n = c(21679, 3428, 3349, 3062, 2577, 176),
    excl = c("", "\u221218,251 no fasting\nblood or age < 50", "\u221279 baseline\nstroke",
             "\u2212287 lost to\nfollow-up", "\u2212485 missing\ncovariates", ""),
    stringsAsFactors = FALSE),
  accent = "#9B5B73", fill = "#F7EEF2", n_max = 6)

fig <- a + b + c + plot_layout(nrow = 1)

w_in <- 183 / 25.4
h_in <- 118 / 25.4

# render-time panel alignment gate (mandatory for multi-panel figures)
# 本机 grid null 单位换算失效，改用视口法测量（measure_panels_vp.R），
# 同一后端中立 JSON 审计器复核。
source(file.path(OUT, "measure_panels_vp.R"))
manifest_path <- file.path(OUT, "Fig1_STROBE_v3.alignment-layout.json")
measure_patchwork_panels(fig, w_in, h_in, manifest_path)
py <- "C:/Users/ronaldchen/AppData/Roaming/kimi-desktop/daimon-share/daimon/runtime/python/python.exe"
if (!file.exists(py)) py <- "python"
status <- system2(py, c(file.path(SKILL_SCRIPTS, "audit_panel_alignment.py"),
                        shQuote(manifest_path),
                        "--json-out", shQuote(file.path(OUT, "Fig1_STROBE_v3.alignment.json")),
                        "--overlay-svg", shQuote(file.path(OUT, "Fig1_STROBE_v3.alignment.svg")),
                        "--tolerance-pt", "1.5", "--gutter-tolerance-pt", "1.5", "--strict"))
if (status != 0) stop("Panel alignment gate failed with exit code ", status)

svglite::svglite(file.path(OUT, "Fig1_STROBE_v3.svg"), width = w_in, height = h_in)
print(fig); dev.off()
ggsave(file.path(OUT, "Fig1_STROBE_v3.tiff"), fig, width = w_in, height = h_in,
       dpi = 600, device = "tiff", compression = "lzw", bg = "white")
ggsave(file.path(OUT, "Fig1_STROBE_v3.pdf"), fig, width = w_in, height = h_in,
       device = "pdf", bg = "white")
cat("done\n")
