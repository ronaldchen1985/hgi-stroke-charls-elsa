####################################################################################
# Fig 1 (v2 redesign) | Three-panel STROBE participant flow
# CHARLS 2011 discovery / CHARLS 2015 landmark / ELSA W2 external evaluation
# Design: panel accent colours, highlighted event boxes, elbow exclusion labels
# Output: Fig1_STROBE_v2.tiff (600 dpi, 4322 px wide) + pdf
####################################################################################
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(grid) })
OUT <- "D:/CHARLS/HGI_stroke/freeze_20261120/figures"

panel_theme <- function(accent, fill) {
  list(accent = accent, fill = fill)
}

flow_col <- function(tag, subtitle, stages, accent, fill, n_max) {
  # stages: data.frame(l1, n, excl)  excl = character vector ("" or "−..." possibly with \n)
  k <- nrow(stages)
  pitch <- 1.18; hh <- 0.42
  y <- (k:1) * pitch + (n_max - k) * pitch      # bottom-align columns
  d <- stages; d$y <- y
  is_event <- seq_len(k) == k

  # exclusion labels in the gaps, with elbow connector
  ex <- do.call(rbind, lapply(2:k, function(i) {
    if (!nzchar(d$excl[i])) return(NULL)
    data.frame(x = 10.55, y = y[i - 1] - pitch / 2, txt = d$excl[i])
  }))
  elbows <- if (!is.null(ex)) data.frame(x0 = 10, x1 = 10.45,
             y = sapply(2:k, function(i) y[i - 1] - pitch / 2)[nzchar(d$excl[2:k]) == 1 | nzchar(d$excl[2:k])]) else NULL

  box_fill <- ifelse(is_event, accent, fill)
  box_col  <- ifelse(is_event, accent, "grey55")
  txt_col  <- ifelse(is_event, "white", "grey15")

  label_df <- do.call(rbind, lapply(seq_len(k), function(i) {
    data.frame(x = 5, y = y[i] + c(0.15, -0.15),
               txt = c(d$l1[i], sprintf("(n = %s)", formatC(d$n[i], format = "f", big.mark = ",", digits = 0))),
               font = c(1, 2), sz = c(1.9, 2.0), col = txt_col[i], stringsAsFactors = FALSE)
  }))

  ggplot(d) +
    geom_rect(aes(xmin = 0, xmax = 10, ymin = y - hh, ymax = y + hh),
              fill = box_fill, colour = box_col, linewidth = 0.45) +
    geom_text(data = label_df,
              aes(x = x, y = y, label = txt, fontface = font, colour = I(col)),
              size = label_df$sz, family = "sans", lineheight = 0.95) +
    geom_segment(data = data.frame(y0 = y[-k] - hh, y1 = y[-k] - pitch + hh),
                 aes(x = 5, xend = 5, y = y0, yend = y1),
                 arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
                 linewidth = 0.5, colour = "grey35") +
    { if (!is.null(ex)) geom_segment(data = elbows, aes(x = x0, xend = x1, y = y, yend = y),
                                     linewidth = 0.3, colour = "grey60") } +
    { if (!is.null(ex)) geom_text(data = ex, aes(x = x, y = y, label = txt), hjust = 0,
                                  size = 1.75, colour = "grey40", family = "sans", lineheight = 0.95) } +
    # panel tag chip + subtitle
    annotate("rect", xmin = 0, xmax = 0.62, ymin = n_max * pitch + hh + 0.42,
             ymax = n_max * pitch + hh + 1.02, fill = accent, colour = NA) +
    annotate("text", x = 0.31, y = n_max * pitch + hh + 0.72, label = tag,
             colour = "white", fontface = "bold", size = 3.2, family = "sans") +
    annotate("text", x = 0.85, y = n_max * pitch + hh + 0.72, label = subtitle, hjust = 0,
             size = 2.3, colour = "grey20", fontface = "bold", family = "sans") +
    coord_cartesian(xlim = c(0, 16.8), ylim = c(0.3, n_max * pitch + hh + 1.15), clip = "off") +
    theme_void(base_family = "sans") +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(2, 1, 2, 1))
}

a <- flow_col("a", "CHARLS 2011 main analysis (discovery)",
  data.frame(
    l1 = c("FPG and HbA1c measured
at Wave 1 (2011)", "Age ≥ 45 years",
           "Free of stroke at baseline", "≥ 1 follow-up interview
during 2013–2023",
           "Complete covariates
(main analysis set)", "Incident stroke during 2011–2023"),
    n = c(11532, 11305, 11038, 10651, 8684, 852),
    excl = c("", "−227 age < 45", "−267 baseline stroke", "−387 lost to follow-up",
             "−8 missing sex\n−1,959 missing covariates", ""),
    stringsAsFactors = FALSE),
  accent = "#1F6FA8", fill = "#E9F1F8", n_max = 6)

b <- flow_col("b", "CHARLS 2015 landmark (internal validation)",
  data.frame(
    l1 = c("Blood drawn at both 2011 and 2015", "Free of stroke before 2015",
           "≥ 1 follow-up interview 2018–2023,\ncomplete covariates (2015)",
           "Incident stroke during 2015–2023"),
    n = c(7166, 6730, 6581, 706),
    excl = c("", "−436 stroke before 2015", "−149 missing covariates", ""),
    stringsAsFactors = FALSE),
  accent = "#2E8B74", fill = "#EAF5F0", n_max = 6)

c <- flow_col("c", "ELSA Wave 2 (external validation)",
  data.frame(
    l1 = c("Harmonized ELSA participants",
           "Fasting FPG and HbA1c at Wave 2,\nage ≥ 50 years (2004/05)",
           "Free of stroke at baseline", "≥ 1 follow-up interview
during Waves 3–11",
           "Complete covariates
(main analysis set)", "Incident stroke during 2006–2023"),
    n = c(21679, 3428, 3349, 3062, 2577, 176),
    excl = c("", "−18,251 no fasting blood
or age < 50", "−79 baseline stroke",
             "−287 lost to follow-up", "−485 missing covariates", ""),
    stringsAsFactors = FALSE),
  accent = "#9B5B73", fill = "#F7EEF2", n_max = 6)

fig <- a + b + c + plot_layout(nrow = 1)

w_in <- 7.2047   # 4322 px @600dpi
h_in <- 4.65
ggsave(file.path(OUT, "Fig1_STROBE_v2.tiff"), fig, width = w_in, height = h_in,
       dpi = 600, device = "tiff", compression = "lzw", bg = "white")
ggsave(file.path(OUT, "Fig1_STROBE_v2.pdf"), fig, width = w_in, height = h_in,
       device = "pdf", bg = "white")
cat("done\n")
