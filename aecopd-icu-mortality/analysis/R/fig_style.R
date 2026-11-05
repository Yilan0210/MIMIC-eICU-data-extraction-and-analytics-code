# =====================================================================
# fig_style.R: common style for the figures (following the figure requirements of journals such as Lancet, JAMA and BMJ)
#   * sans-serif font (Arial if available), 7–8 pt; no in-plot titles or grid lines; thin lines
#   * colour-blind-friendly Okabe-Ito palette; panel letters for multi-panel figures
#   * sizes: single column 89 mm, 1.5 columns 120 mm, double column 183 mm
#   * each figure saved as vector PDF (embedded fonts) and 600 dpi PNG
# =====================================================================
suppressPackageStartupMessages({
  for (p in c("ggplot2", "patchwork")) if (!requireNamespace(p, quietly = TRUE))
    stop("Missing R package ", p, ". Install first: install.packages(c('ggplot2','patchwork'))")
  library(ggplot2); library(patchwork)
})

HAS_CAIRO <- isTRUE(capabilities("cairo"))
# Font: PDF (cairo) and PNG (ragg) must use the same font with both regular and bold faces,
# otherwise regular and bold text may come from different fonts. The first installed font in the list is used.
FIG_FONT <- local({
  f <- Sys.getenv("AECOPD_FONT", "")
  if (nzchar(f)) return(f)
  if (!requireNamespace("systemfonts", quietly = TRUE)) return("sans")
  sf <- systemfonts::system_fonts()
  is_bold <- if ("bold" %in% names(sf)) as.logical(sf$bold) else
    as.character(sf$weight) %in% c("semibold", "bold", "ultrabold", "heavy") | grepl("Bold", sf$style)
  is_reg <- !is_bold & (as.character(sf$weight) %in% c("normal", "regular", "book") | grepl("^(Regular|Book|Roman)$", sf$style))
  upright <- !as.logical(sf$italic) %in% TRUE
  ok <- function(fam) { k <- sf$family == fam & upright; any(k & is_reg) && any(k & is_bold) }
  cand <- c("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "Nimbus Sans L", "TeX Gyre Heros", "Noto Sans", "DejaVu Sans")
  hit <- cand[vapply(cand, ok, logical(1))]
  if (length(hit)) hit[1] else "sans"
})
message("Figure font: ", FIG_FONT)
# Unit "x10^9/L": many fonts lack the superscript 9 (shown as a box), so axis titles use a plotmath superscript
u109 <- function(prefix) bquote(.(prefix) * ", ×" * 10^9 * "/L")
W1 <- 89; W15 <- 120; W2 <- 183            # journal column widths (mm)
BASE <- 7.5                                 # base font size (pt, at final size)

PAL <- c(M0 = "#000000", M1 = "#0072B2", M2 = "#D55E00", BAP65 = "#999999", DECAF_mod = "#CC79A7",
         APACHE_IVa = "#009E73", SAPSII = "#E69F00", SOFA = "#56B4E9", treat_all = "#8C8C8C")
LAB <- c(M0 = "M0: clinical severity", M1 = "M1: M0 + ln neutrophils, ln lymphocytes, albumin",
         M2 = "M2: M1 + eosinopenia", BAP65 = "BAP-65", DECAF_mod = "Modified DECAF",
         APACHE_IVa = "APACHE IVa", SAPSII = "SAPS II", SOFA = "SOFA", treat_all = "Treat all", treat_none = "Treat none")

theme_journal <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FIG_FONT) +
    theme(axis.line = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks.length = unit(1.2, "mm"),
          axis.text = element_text(colour = "black", size = base - 0.5),
          axis.title = element_text(colour = "black", size = base),
          legend.title = element_blank(), legend.text = element_text(size = base - 0.5),
          legend.key.size = unit(3.2, "mm"), legend.key.width = unit(6, "mm"),
          legend.background = element_blank(), legend.spacing.y = unit(0.5, "mm"),
          plot.tag = element_text(face = "bold", size = base + 1.5),
          plot.margin = margin(2, 3, 2, 2, "mm"),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = base))
}

save_fig <- function(p, name, w_mm, h_mm) {
  # panel letters of all composite figures in bold (wrap_elements panels do not inherit theme_journal)
  if (inherits(p, "patchwork")) p <- p & theme(plot.tag = element_text(face = "bold", size = BASE + 1.5, family = FIG_FONT))
  pdf_dev <- if (HAS_CAIRO) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(file.path(DIR_FIG, paste0(name, ".pdf")), p, width = w_mm, height = h_mm, units = "mm", device = pdf_dev)
  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(file.path(DIR_FIG, paste0(name, ".png")), p, width = w_mm, height = h_mm, units = "mm",
                  dpi = 600, device = png_dev, bg = "white")
  msg("Figure: ", name, " (PDF + PNG, ", w_mm, "×", h_mm, " mm)")
}

# true minus sign in axis tick labels
lab_minus <- function(x) {      # same number of decimals on all ticks; true minus sign
  ok <- x[is.finite(x)]
  d <- if (length(ok)) max(nchar(sub("^[^.]*\\.?", "", format(ok, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)))) else 0
  out <- formatC(x, format = "f", digits = min(d, 4)); out[!is.finite(x)] <- NA; sub("^-", "\u2212", out) }
# "deaths/patients" with thousands separators: 230/1,858
dn_txt <- function(d, n) paste0(prettyNum(d, big.mark = ","), "/", prettyNum(n, big.mark = ","))
# number format: "0.71 (0.68–0.74)"; "to" when an interval contains negative values; true minus sign
num <- function(x, d) { x <- round(x, d); x[!is.na(x) & x == 0] <- 0   # avoid "−0.000"
  s <- formatC(x, format = "f", digits = d); sub("^-", "−", s) }
# Nimbus Sans has no subscript 2 (U+2082) and would fall back to another font; written as PaCO2 in figures (tables unaffected)
if (exists("VAR_LAB")) VAR_LAB <- gsub("\u2082", "2", VAR_LAB, fixed = TRUE)
ci_txt <- function(e, l, h, d = 2) ifelse(is.na(l), num(e, d),
  ifelse(l < 0 | h < 0, paste0(num(e, d), " (", num(l, d), " to ", num(h, d), ")"),
         paste0(num(e, d), " (", num(l, d), "–", num(h, d), ")")))

# ---------------------------------------------------------------------
# Forest plot: label column on the left + point estimates with 95% CI in the middle + value column on the right
#   d: label, sub (second column, e.g. deaths/patients), est, lo, hi, type ("row" / "pooled" / "header"),
#      pooled rows may have pi_lo, pi_hi (95% prediction interval, drawn in the row below)
# ---------------------------------------------------------------------
forest_plot <- function(d, xlab, ref = NULL, threshold = NULL, log = FALSE, digits = 2,
                        head = c("", "", "Estimate (95% CI)"), tag = NULL, xlim = NULL,
                        widths = c(2.9, 3.3, 1.8), sub_x = 0.86) {
  if (!"pi_lo" %in% names(d)) { d$pi_lo <- NA; d$pi_hi <- NA }
  if (!"sub" %in% names(d)) d$sub <- ""
  rows <- list(); y <- 0
  for (i in rev(seq_len(nrow(d)))) {
    r <- d[i, ]
    if (r$type == "pooled" && !is.na(r$pi_lo)) {
      y <- y + 1; rows[[length(rows) + 1]] <- data.frame(y = y, kind = "pi", label = "   95% prediction interval", sub = "",
        est = NA, lo = r$pi_lo, hi = r$pi_hi,
        txt = paste0("(", num(r$pi_lo, digits), if (r$pi_lo < 0) " to " else "–", num(r$pi_hi, digits), ")"))
    }
    y <- y + 1
    rows[[length(rows) + 1]] <- data.frame(y = y, kind = r$type, label = r$label, sub = r$sub, est = r$est, lo = r$lo, hi = r$hi,
                                           txt = if (r$type == "header") "" else ci_txt(r$est, r$lo, r$hi, digits))
  }
  g <- do.call(rbind, rows); top <- max(g$y) + 1.2
  ylims <- c(0.4, top + 0.4)
  bold <- ifelse(g$kind %in% c("pooled", "header"), "bold", "plain")
  void <- theme_void(base_size = BASE, base_family = FIG_FONT) + theme(plot.margin = margin(2, 0, 2, 1, "mm"),
                                                                       plot.tag = element_text(face = "bold", size = BASE + 1.5))
  pl <- ggplot(g) + geom_text(aes(x = 0, y = y, label = label), hjust = 0, size = BASE / .pt, family = FIG_FONT, fontface = bold) +
    geom_text(aes(x = sub_x, y = y, label = sub), hjust = 0.5, size = BASE / .pt, family = FIG_FONT) +
    annotate("text", x = 0, y = top, label = head[1], hjust = 0, fontface = "bold", size = BASE / .pt, family = FIG_FONT) +
    annotate("text", x = sub_x, y = top, label = head[2], hjust = 0.5, fontface = "bold", size = BASE / .pt, family = FIG_FONT) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = ylims, expand = c(0, 0)) + void
  if (!is.null(tag)) pl <- pl + labs(tag = tag)
  pr <- ggplot(g) + geom_text(aes(x = 0, y = y, label = txt), hjust = 0, size = BASE / .pt, family = FIG_FONT, fontface = bold) +
    annotate("text", x = 0, y = top, label = head[3], hjust = 0, fontface = "bold", size = BASE / .pt, family = FIG_FONT) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = ylims, expand = c(0, 0)) + void
  rr <- g[g$kind == "row" & !is.na(g$est), ]; pp <- g[g$kind == "pooled", ]; pi <- g[g$kind == "pi", ]
  dia <- if (nrow(pp)) do.call(rbind, lapply(seq_len(nrow(pp)), function(i) with(pp[i, ],
    data.frame(id = i, x = c(lo, est, hi, est), y = c(y, y + 0.32, y, y - 0.32))))) else NULL
  pm <- ggplot() + theme_journal() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          axis.title.y = element_blank(), plot.margin = margin(2, 1, 2, 1, "mm"))
  if (!is.null(ref)) pm <- pm + geom_vline(xintercept = ref, linewidth = 0.3, colour = "grey40")
  if (!is.null(threshold)) pm <- pm + geom_vline(xintercept = threshold, linewidth = 0.3, linetype = "dashed", colour = "grey40")
  pm <- pm + geom_errorbar(orientation = "y", data = rr[!is.na(rr$lo), ], aes(xmin = lo, xmax = hi, y = y), width = 0, linewidth = 0.35) +
    geom_point(data = rr, aes(x = est, y = y), shape = 15, size = 1.6)
  if (!is.null(dia)) pm <- pm + geom_polygon(data = dia, aes(x = x, y = y, group = id), fill = "black")
  if (nrow(pi)) pm <- pm + geom_errorbar(orientation = "y", data = pi, aes(xmin = lo, xmax = hi, y = y), width = 0.25, linewidth = 0.35)
  pm <- pm + scale_y_continuous(limits = ylims, expand = c(0, 0)) + labs(x = xlab)
  if (log) { rg <- range(c(g$lo, g$hi, g$est, xlim), na.rm = TRUE); lb <- c(0.1, 0.2, 0.25, 0.33, 0.5, 0.67, 0.8, 1, 1.25, 1.5, 2, 3, 4, 5, 10)
    lb <- lb[lb >= rg[1] / 1.05 & lb <= rg[2] * 1.05]; if (length(lb) > 6) lb <- lb[lb %in% c(0.1, 0.25, 0.5, 0.8, 1, 1.25, 2, 4, 10)] }
  pm <- pm + if (log) scale_x_log10(limits = xlim, labels = lab_minus, breaks = lb) else scale_x_continuous(limits = xlim, labels = lab_minus, breaks = function(l) pretty(l, n = 4))
  pl + pm + pr + plot_layout(widths = widths)
}
forest_height <- function(n_rows) max(42, 13 + 4.2 * n_rows)   # mm
