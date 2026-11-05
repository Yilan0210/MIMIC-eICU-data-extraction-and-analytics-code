# =====================================================================
# 11_figures_extra.R: figures for additional validation V1–V9 (reads local/10_extra.rds); journal style, PDF + 600 dpi PNG
#   FigX1 temporal validation        FigX2 hospital level (funnel plots)     FigX3 decomposition of case mix and coefficient differences
#   FigX4 uncertainty in clinical utility  FigX5 fairness calibration curves   FigX6 reclassification (alluvial plot)
#   FigX7 multiverse analysis        FigX8 performance heat map by hospital type and case definition
#   working figure numbers; see README for the numbering in the paper
# =====================================================================
source("00_setup.R"); source("fig_style.R")
ex  <- readRDS(file.path(DIR_LOCAL, "10_extra.rds"))
r03 <- readRDS(file.path(DIR_LOCAL, "03_iecv_results.rds"))
legends <- character()
pct_lab <- function(x) paste0(round(100 * x), "%")
COL_MIMIC <- "#0072B2"; COL_EICU <- "#E69F00"; COL_NEW <- "#0072B2"
DB_COL <- c(MIMIC = COL_MIMIC, eICU = COL_EICU); DB_LAB <- c(MIMIC = "MIMIC-IV", eICU = "eICU-CRD")
pooledC <- r03$pooled$est[r03$pooled$measure == "logit_C" & r03$pooled$model == "full"][1]
f2 <- function(x) num(x, 2); f3 <- function(x) num(x, 3)
txt_ci <- function(e, l, h, d = 2) ci_txt(e, l, h, d)

# =====================================================================
# FigX1 temporal validation
# =====================================================================
v1 <- ex$v1
if (!is.null(v1$drift) && nrow(v1$drift)) {
  msg("FigX1 temporal validation")
  dr <- v1$drift; dr$dbl <- DB_LAB[dr$db]; dr$x <- dr$year_mid + ifelse(dr$db == "eICU", 0.35, 0)   # offset from the MIMIC points of the same year
  pa <- ggplot(dr, aes(x, C, colour = db)) + geom_hline(yintercept = pooledC, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
    geom_line(linewidth = 0.4) + geom_errorbar(aes(ymin = C_lo, ymax = C_hi), width = 0.25, linewidth = 0.35) + geom_point(size = 1.6) +
    scale_colour_manual(values = DB_COL, labels = DB_LAB) + scale_x_continuous(breaks = sort(unique(round(dr$year_mid)))) +
    labs(x = "Year (midpoint of period)", y = "C statistic (held-out predictions)") + theme_journal() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
  pb <- ggplot(dr, aes(x, OE, colour = db)) + geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey50") +
    geom_line(linewidth = 0.4) + geom_errorbar(aes(ymin = OE_lo, ymax = OE_hi), width = 0.25, linewidth = 0.35) + geom_point(size = 1.6) +
    scale_colour_manual(values = DB_COL, labels = DB_LAB) + scale_y_log10(labels = lab_minus) +
    scale_x_continuous(breaks = sort(unique(round(dr$year_mid)))) +
    labs(x = "Year (midpoint of period)", y = "Observed/expected ratio (log scale)") + theme_journal() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
  pl <- list(pa, pb)
  if (!is.null(v1$fwd_curves) && nrow(v1$fwd_curves)) {
    cv <- v1$fwd_curves; fw <- v1$fwd
    labs_d <- vapply(unique(cv$design), function(d) { s <- fw[fw$design == d & fw$model == "full", ]
      g <- function(m) s[s$measure == m, ]
      sprintf("%s: C %s, O/E %s, slope %s", sub(" \\(developed.*$", "", d), f2(g("logit_C")$est), f2(g("log_OE")$est), f2(g("cal_slope")$est)) }, character(1))
    cv$dl <- factor(labs_d[cv$design], levels = labs_d)
    mx <- min(1, max(cv$pred, cv$obs) * 1.1)
    pc <- ggplot(cv, aes(pred, obs, colour = dl, fill = dl)) + geom_abline(linetype = "dashed", linewidth = 0.3, colour = "grey50") +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.5) +
      scale_colour_manual(values = c("#D55E00", "#009E73")) + scale_fill_manual(values = c("#D55E00", "#009E73")) +
      scale_x_continuous(labels = pct_lab) + scale_y_continuous(labels = pct_lab) + coord_cartesian(xlim = c(0, mx), ylim = c(0, mx)) +
      labs(x = "Predicted risk", y = "Observed proportion") + theme_journal() + theme(legend.position = "bottom", legend.direction = "vertical")
    pl <- c(pl, list(pc))
  }
  fx1 <- wrap_plots(pl, nrow = 1) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom", legend.direction = "vertical")
  save_fig(fx1, "FigX1_temporal_validation", W2, 82)
  legends["FigX1_temporal_validation"] <- paste("Temporal validation. (A, B) C statistic and observed/expected ratio of held-out predictions from",
    "internal-external cross-validation by calendar period (MIMIC-IV anchor-year groups; eICU-CRD discharge year); the dashed line marks the pooled C.",
    if (length(pl) == 3) "(C) Forward temporal validation: the full modelling pipeline was repeated in all other patients and the model was applied to the latest periods; flexible calibration curves with 95% CIs." else "")
}

# =====================================================================
# FigX2 hospital level: funnel plots (new model vs APACHE IVa) and agreement of hospital O/E ratios
# =====================================================================
v2 <- ex$v2
if (!is.null(v2)) {
  msg("FigX2 hospital level")
  funnel_lim <- function(E, p) { r <- stats::qpois(p, E); al <- (stats::ppois(r, E) - p) / stats::dpois(r, E); pmax((r - al) / E, 0) }
  YMAX <- 4
  fun_plot <- function(h, title, xmax, note = "") {
    g <- exp(seq(log(max(0.05, min(h$E) * 0.9)), log(xmax), length.out = 300))
    lim <- data.frame(E = g, l95 = funnel_lim(g, 0.025), u95 = funnel_lim(g, 0.975), l998 = funnel_lim(g, 0.001), u998 = funnel_lim(g, 0.999))
    h$y <- pmin(h$SMR, YMAX); h$clip <- h$SMR > YMAX
    h$flag95 <- factor(h$flag95, c("Lower than expected", "As expected", "Higher than expected"))
    ggplot(h, aes(E, y)) + geom_hline(yintercept = 1, linewidth = 0.3) +
      geom_line(data = lim, aes(E, u95), linetype = "dashed", linewidth = 0.3, inherit.aes = FALSE) +
      geom_line(data = lim, aes(E, l95), linetype = "dashed", linewidth = 0.3, inherit.aes = FALSE) +
      geom_line(data = lim, aes(E, u998), linetype = "dotted", linewidth = 0.35, inherit.aes = FALSE) +
      geom_line(data = lim, aes(E, l998), linetype = "dotted", linewidth = 0.35, inherit.aes = FALSE) +
      geom_point(aes(fill = flag95, shape = clip), size = 1.5, colour = "grey20", stroke = 0.2, alpha = 0.9) +
      scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none") +
      scale_fill_manual(values = c("Lower than expected" = "#56B4E9", "As expected" = "grey80", "Higher than expected" = "#D55E00"), drop = FALSE,
                        name = "Funnel classification (95% limits)") +
      guides(fill = guide_legend(override.aes = list(shape = 21, size = 2, colour = "grey20", fill = c("#56B4E9", "grey80", "#D55E00")))) +
      scale_x_log10(limits = c(NA, xmax)) + coord_cartesian(ylim = c(0, YMAX)) +
      labs(x = "Expected deaths in hospital (log scale)", y = "Observed/expected deaths", subtitle = if (nzchar(note)) paste0(title, "\n", note) else title) + theme_journal() +
      theme(legend.position = "bottom", plot.subtitle = element_text(size = BASE, face = "bold"))
  }
  xmax <- max(c(v2$hm$E, v2$ha_ap$E)) * 1.1
  s <- v2$sum; od <- function(m, p) s$overdispersion[s$model == m & s$patients == p]
  pa <- fun_plot(v2$hm, "New model", xmax, sprintf("%d hospitals; dispersion %.2f", nrow(v2$hm), od("New model", "All eICU-CRD patients")))
  pb <- fun_plot(v2$ha_ap, "APACHE IVa", xmax, sprintf("%d hospitals; dispersion %.2f", nrow(v2$ha_ap), od("APACHE IVa", "Patients with APACHE IVa")))
  ag <- v2$agree; ag <- ag[ag$E_model >= 2 & ag$E_apache >= 2, ]
  rho <- if (nrow(ag) >= 5) stats::cor(ag$SMR_model, ag$SMR_apache, method = "spearman") else NA
  lim2 <- c(0, max(c(ag$SMR_model, ag$SMR_apache, 1.5), na.rm = TRUE) * 1.05)
  pc <- ggplot(ag, aes(SMR_apache, SMR_model)) + geom_abline(linetype = "dashed", linewidth = 0.3, colour = "grey50") +
    geom_hline(yintercept = 1, linewidth = 0.2, colour = "grey70") + geom_vline(xintercept = 1, linewidth = 0.2, colour = "grey70") +
    geom_point(aes(size = n), shape = 21, fill = COL_NEW, colour = "white", alpha = 0.8, stroke = 0.2) +
    scale_size_area(max_size = 3.5, guide = "none") +
    annotate("text", x = 0, y = lim2[2], hjust = 0, vjust = 1, size = (BASE - 1) / .pt, family = FIG_FONT, lineheight = 0.9,
             label = sprintf("Spearman \u03c1 = %s\n%d hospitals with \u22652\nexpected deaths", if (is.na(rho)) "NA" else f2(rho), nrow(ag))) +
    labs(x = "O/E with APACHE IVa", y = "O/E with the new model", subtitle = "Agreement between models") +
    coord_equal(xlim = lim2, ylim = lim2) + theme_journal() + theme(plot.subtitle = element_text(size = BASE, face = "bold"))
  pa <- pa + guides(fill = "none")   # legend taken from panel B only (panel A may have no hospital "below expected", so the legends cannot be merged)
  fx2 <- (pa + pb + pc) + plot_layout(guides = "collect", widths = c(1, 1, 0.9)) + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
  save_fig(fx2, "FigX2_hospital_benchmarking", W2, 85)
  legends["FigX2_hospital_benchmarking"] <- paste("Hospital-level performance and benchmarking in eICU-CRD. (A, B) Funnel plots of observed/expected deaths per hospital",
    "with the new model (held-out predictions from internal-external cross-validation) and with APACHE IVa; dashed and dotted lines are 95% and 99.8%",
    "control limits (Poisson, Spiegelhalter 2005); triangles mark ratios above 4 drawn at 4; dispersion is the mean squared winsorised z score",
    "(values above 1 indicate between-hospital variation beyond chance). (C) Agreement of hospital observed/expected ratios between the two models.")
}

# =====================================================================
# FigX3 decomposition of case mix and coefficient differences (van Klaveren 2016)
# =====================================================================
v4 <- ex$v4
if (!is.null(v4)) {
  msg("FigX3 case-mix decomposition")
  v4$lab <- paste0(v4$cluster, " ", v4$cluster_name)
  v4 <- v4[order(as.integer(as.character(v4$cluster))), ]
  v4$xl <- factor(paste0(v4$cluster, "\n", gsub(" ", "\n", sub("^eICU ", "", sub(" era$", "", v4$cluster_name)))), levels = unique(
    paste0(v4$cluster, "\n", gsub(" ", "\n", sub("^eICU ", "", sub(" era$", "", v4$cluster_name))))))
  pa <- ggplot(v4, aes(xl)) +
    geom_segment(aes(xend = xl, y = C_model_based, yend = C_observed), linewidth = 0.4, colour = "grey60") +
    geom_errorbar(aes(ymin = C_obs_lo, ymax = C_obs_hi, colour = db), width = 0, linewidth = 0.35, alpha = 0.6, show.legend = FALSE) +
    geom_point(aes(y = C_model_based, shape = "Model-based C (expected from case mix)"), size = 2.2, colour = "black", stroke = 0.6) +
    geom_point(aes(y = C_observed, colour = db, shape = "Observed C (95% CI)"), size = 2.4) +
    scale_shape_manual(values = c("Model-based C (expected from case mix)" = 5, "Observed C (95% CI)" = 16), name = NULL) +
    scale_colour_manual(values = DB_COL, labels = DB_LAB) +
    guides(shape = guide_legend(override.aes = list(colour = "black"), ncol = 1)) +
    labs(x = "Held-out cluster", y = "C statistic") + theme_journal() + theme(legend.position = "bottom", axis.text.x = element_text(lineheight = 0.85))
  pb <- ggplot(v4, aes(cal_slope, difference, colour = db)) + geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = diff_lo, ymax = diff_hi), width = 0, linewidth = 0.35, alpha = 0.7, show.legend = FALSE) + geom_point(aes(size = n)) +
    geom_text(aes(label = cluster), colour = "white", size = (BASE - 2) / .pt, family = FIG_FONT, fontface = "bold") +
    scale_colour_manual(values = DB_COL, labels = DB_LAB, guide = "none") + scale_size_area(max_size = 5, guide = "none") + scale_y_continuous(labels = lab_minus) +
    labs(x = "Calibration slope in held-out cluster", y = "Observed minus model-based C") + theme_journal() + theme(legend.position = "bottom")
  key <- paste(paste0(v4$cluster, " ", v4$cluster_name), collapse = "; ")
  fx3 <- (pa + pb) + plot_layout(guides = "collect", widths = c(1.2, 1)) + plot_annotation(tag_levels = "A", caption = key) &
    theme(legend.position = "bottom", plot.caption = element_text(size = BASE - 1.5, hjust = 0, family = FIG_FONT))
  save_fig(fx3, "FigX3_casemix_decomposition", W2, 92)
  legends["FigX3_casemix_decomposition"] <- paste("Why performance differs between clusters. (A) Observed C statistic in each held-out cluster (95% CI) and the model-based",
    "C (van Klaveren et al., 2016), the C expected if the model's predicted risks were correct given that cluster's case mix. Similar model-based C",
    "across clusters indicates similar case-mix heterogeneity; an observed C well below it indicates that predictor effects do not transport.",
    "(B) Observed minus model-based C (bootstrap 95% CI) against the calibration slope. Numbers identify clusters; point size is proportional to cluster size.")
}

# =====================================================================
# FigX4 uncertainty in clinical utility (Bayesian bootstrap) and EVPI
# =====================================================================
v6 <- ex$v6
if (!is.null(v6)) {
  msg("FigX4 value of information")
  SNAME <- c(full = "New model", BAP65 = "BAP-65", DECAF_mod = "Modified DECAF", APACHE_IVa = "APACHE IVa", SAPSII = "SAPS II", SOFA = "SOFA",
             treat_all = "Treat all", treat_none = "Treat none")
  SCOL <- c(full = COL_NEW, BAP65 = "#7F7F7F", DECAF_mod = "#CC79A7", APACHE_IVa = "#009E73", SAPSII = "#E69F00", SOFA = "#D55E00",
            treat_all = "#B0B0B0", treat_none = "#000000")
  shade <- annotate("rect", xmin = 0.10, xmax = 0.30, ymin = -Inf, ymax = Inf, fill = "grey93")
  d1 <- v6[v6$subset == "All patients", ]
  brk <- intersect(names(SNAME), unique(d1$strategy))
  pa <- ggplot(d1, aes(threshold, p_best, colour = strategy)) + shade + geom_line(linewidth = 0.5) +
    scale_colour_manual(values = SCOL[brk], labels = SNAME[brk], breaks = brk) +
    scale_x_continuous(labels = pct_lab) + scale_y_continuous(labels = pct_lab, limits = c(0, 1)) +
    labs(x = "Threshold probability", y = "Probability of highest net benefit", subtitle = "All patients") + theme_journal() +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = BASE, face = "bold")) + guides(colour = guide_legend(ncol = 2))
  ku <- unique(v6[, c("subset", "threshold", "p_model_useful", "EVPI_model_vs_defaults", "EVPI_all")])
  SUBCOL <- c("All patients" = COL_NEW, "eICU-CRD, APACHE IVa available" = "#009E73", "MIMIC-IV" = "#E69F00")
  pb <- ggplot(ku, aes(threshold, p_model_useful, colour = subset)) + shade + geom_line(linewidth = 0.5) +
    scale_colour_manual(values = SUBCOL) + scale_x_continuous(labels = pct_lab) + scale_y_continuous(labels = pct_lab, limits = c(0, 1)) +
    labs(x = "Threshold probability", y = "Probability the model is useful") + theme_journal() + theme(legend.position = "bottom", legend.direction = "vertical")
  pc <- ggplot(ku, aes(threshold, 1000 * EVPI_all, colour = subset)) + shade + geom_line(linewidth = 0.5) +
    scale_colour_manual(values = SUBCOL) + scale_x_continuous(labels = pct_lab) +
    labs(x = "Threshold probability", y = "EVPI (net true positives per 1,000 patients)") + theme_journal() +
    theme(legend.position = "bottom", legend.direction = "vertical")
  inner <- (pb + pc) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.direction = "vertical")
  fx4 <- (pa | inner) + plot_layout(widths = c(1, 2)) + plot_annotation(tag_levels = "A")
  save_fig(fx4, "FigX4_decision_uncertainty", W2, 85)
  legends["FigX4_decision_uncertainty"] <- paste("Uncertainty in clinical utility (Bayesian bootstrap of held-out predictions). (A) Probability that each strategy has the",
    "highest net benefit at each threshold. (B) Probability that the new model has higher net benefit than treating all and treating none.",
    "(C) Expected value of perfect information for validation (EVPI; Sadatsafavi et al.): the expected gain in net benefit from eliminating uncertainty",
    "about which strategy is best; values near zero indicate that further validation data would be unlikely to change the decision. Shading marks 10%-30%.")
}

# =====================================================================
# FigX5 fairness: calibration curves by sex, race and age
# =====================================================================
v7 <- ex$v7
if (!is.null(v7$curves) && nrow(v7$curves)) {
  msg("FigX5 fairness")
  cv <- v7$curves; tb <- v7$tab
  mx <- min(0.6, max(cv$pred, cv$obs) * 1.05)
  pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
  one <- function(g) { z <- cv[cv$group == g, ]; t <- tb[tb$group == g, ]
    lab <- vapply(unique(z$level), function(l) { r <- t[t$level == l, ]
      if (!nrow(r)) return(l)
      sprintf("%s (n = %s)\nO/E %s, slope %s", sub(">=", "\u2265", l), prettyNum(r$n, big.mark = ","), f2(r$OE), f2(r$slope)) }, character(1))
    z$ll <- factor(lab[z$level], levels = lab)
    ggplot(z, aes(pred, obs, colour = ll, fill = ll)) + geom_abline(linetype = "dashed", linewidth = 0.3, colour = "grey50") +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 0.5) +
      scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
      scale_x_continuous(labels = pct_lab) + scale_y_continuous(labels = pct_lab) + coord_cartesian(xlim = c(0, mx), ylim = c(0, mx)) +
      labs(x = "Predicted risk", y = "Observed proportion", subtitle = g) + theme_journal() +
      theme(legend.position = "bottom", legend.direction = "vertical", plot.subtitle = element_text(size = BASE, face = "bold"),
            legend.text = element_text(size = BASE - 1.5, lineheight = 0.9), legend.key.height = unit(5, "mm")) +
      guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) }
  fx5 <- wrap_plots(lapply(unique(cv$group), one), nrow = 1) + plot_annotation(tag_levels = "A")
  save_fig(fx5, "FigX5_fairness_calibration", W2, 105)
  legends["FigX5_fairness_calibration"] <- paste("Calibration by sex, race and ethnicity, and age (held-out predictions averaged across imputations). Flexible calibration",
    "curves (natural splines of the logit of predicted risk) with 95% CIs, restricted to the 1st-99th percentiles of predicted risk in each group;",
    "observed/expected ratios and calibration slopes are pooled across imputations.")
}

# =====================================================================
# FigX6 reclassification: alluvial plot
# =====================================================================
v8 <- ex$v8
if (!is.null(v8) && requireNamespace("ggalluvial", quietly = TRUE)) {
  msg("FigX6 reclassification")
  MAXM <- max(0.4, min(1, stats::quantile(rep(v8$flow$mortality, v8$flow$n), 0.99)))
  alluv <- function(cmp, from_cmp, axis_from, levs_from) {
    d <- v8$flow[v8$flow$comparison == cmp, ]
    # box labels: number of patients and observed mortality, computed for the patients in this panel
    relab <- function(col, lev) { k <- as.character(d[[col]]) == lev; if (!any(k)) return(lev)
      sprintf("%s\n%s (%s%%)", sub("-", "\u2013", sub(">=", "\u2265", lev)), prettyNum(sum(d$n[k]), big.mark = ","), round(100 * sum(d$deaths[k]) / sum(d$n[k]))) }
    lf <- vapply(levs_from, function(l) relab("from", l), character(1)); lt <- vapply(c("<10%", "10-30%", ">=30%"), function(l) relab("to", l), character(1))
    d$from <- factor(lf[as.character(d$from)], levels = lf); d$to <- factor(lt[as.character(d$to)], levels = lt)
    ggplot(d, aes(y = n, axis1 = from, axis2 = to)) +
      ggalluvial::geom_alluvium(aes(fill = mortality), width = 0.32, alpha = 0.95, knot.pos = 0.35, colour = "white", linewidth = 0.1) +
      ggalluvial::geom_stratum(width = 0.32, fill = "grey96", colour = "grey35", linewidth = 0.25) +
      geom_text(stat = ggalluvial::StatStratum, aes(label = after_stat(stratum)), min.y = 0.06 * sum(d$n),
                size = (BASE - 1.8) / .pt, family = FIG_FONT, lineheight = 0.85) +
      geom_text(stat = ggalluvial::StatStratum, aes(label = after_stat(ifelse(x == 1, gsub("\n", " ", as.character(stratum)), NA))), max.y = 0.06 * sum(d$n),
                position = position_nudge(x = -0.19), hjust = 1, size = (BASE - 1.8) / .pt, family = FIG_FONT, lineheight = 0.85, na.rm = TRUE) +
      geom_text(stat = ggalluvial::StatStratum, aes(label = after_stat(ifelse(x == 2, gsub("\n", " ", as.character(stratum)), NA))), max.y = 0.06 * sum(d$n),
                position = position_nudge(x = 0.19), hjust = 0, size = (BASE - 1.8) / .pt, family = FIG_FONT, lineheight = 0.85, na.rm = TRUE) +
      scale_x_discrete(limits = c(axis_from, "New model"), expand = c(0.32, 0.05)) + coord_cartesian(clip = "off") +
      scale_fill_gradient(low = "#FDE0C5", high = "#B2182B", labels = pct_lab, limits = c(0, MAXM), oob = scales::squish, name = "Observed mortality") +
      guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
      labs(x = NULL, y = "Patients") + theme_journal() +
      theme(legend.position = "bottom", legend.title = element_text(size = BASE - 1), legend.key.width = unit(10, "mm"), axis.line = element_blank(),
            axis.ticks = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),   # box labels already give the numbers, so no y-axis is needed
            plot.margin = margin(2, 16, 2, 2, "mm"))                                                      # room on the right for labels of small boxes
  }
  pa <- alluv("BAP-65 class -> new model risk group", "BAP-65 class", "BAP-65 class", c("I", "II", "III", "IV", "V"))
  pl <- list(pa)
  if (any(v8$flow$comparison == "APACHE IVa risk group -> new model risk group (eICU-CRD)"))
    pl <- c(pl, list(alluv("APACHE IVa risk group -> new model risk group (eICU-CRD)", "APACHE IVa risk group", "APACHE IVa", c("<10%", "10-30%", ">=30%"))))
  fx6 <- wrap_plots(pl, nrow = 1) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
  save_fig(fx6, "FigX6_reclassification", W2, 100)
  legends["FigX6_reclassification"] <- paste("Reclassification of patients from (A) BAP-65 classes (all patients) and (B) APACHE IVa risk groups (eICU-CRD) to the",
    "risk groups of the new model (held-out predictions). Box labels give the number of patients and observed in-hospital mortality; ribbon colour",
    "shows the observed mortality of patients following that path. BAP-65 class is the most frequent class across imputed datasets.")
}

# =====================================================================
# FigX7 multiverse (specification curve)
# =====================================================================
v9 <- ex$v9
if (!is.null(v9) && sum(is.finite(v9$C)) > 5) {
  msg("FigX7 multiverse")
  d <- v9[is.finite(v9$C), ]; d <- d[order(d$C), ]; d$rank <- seq_len(nrow(d)); c0 <- d$C[d$primary]
  d$type <- ifelse(d$primary, "Reference specification", "Other specifications")
  pa <- ggplot(d, aes(rank, C)) + geom_hline(yintercept = c0, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
    geom_linerange(aes(ymin = C_lo, ymax = C_hi, colour = type), linewidth = 0.2, alpha = 0.3) + geom_point(aes(colour = type, size = type)) +
    scale_colour_manual(values = c("Reference specification" = "#D55E00", "Other specifications" = COL_NEW)) +
    scale_size_manual(values = c("Reference specification" = 1.6, "Other specifications" = 0.6)) +
    annotate("text", x = d$rank[d$primary], y = d$C_hi[d$primary], vjust = -0.5, label = "Reference", colour = "#D55E00",
             size = (BASE - 1) / .pt, family = FIG_FONT) +
    scale_x_continuous(expand = c(0.01, 0.01)) + labs(x = NULL, y = "Pooled C (IECV)") + theme_journal() +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
  pb <- ggplot(d, aes(rank, slope)) + geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey50") +
    geom_point(aes(colour = type, size = type)) +
    scale_colour_manual(values = c("Reference specification" = "#D55E00", "Other specifications" = "#56B4E9"), guide = "none") +
    scale_size_manual(values = c("Reference specification" = 1.6, "Other specifications" = 0.6), guide = "none") +
    scale_x_continuous(expand = c(0.01, 0.01)) + labs(x = NULL, y = "Calibration slope") + theme_journal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  FAC <- c(missing = "Missing data", form = "Functional form", estimation = "Estimation", added = "Added predictor", outcome = "Outcome", population = "Population")
  long <- do.call(rbind, lapply(names(FAC), function(f) data.frame(rank = d$rank, factor = FAC[[f]], option = d[[f]], primary = d$primary)))
  long$factor <- factor(long$factor, levels = FAC)
  lev <- unlist(lapply(names(FAC), function(f) rev(unique(v9[[f]][order(v9$spec)]))))
  long$option <- factor(long$option, levels = unique(lev))
  pc <- ggplot(long, aes(rank, option)) + geom_point(aes(colour = primary), size = 0.45, shape = 15) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey20"), guide = "none") +
    facet_grid(factor ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(expand = c(0.01, 0.01)) + labs(x = "Specifications ranked by pooled C", y = NULL) + theme_journal() +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = BASE - 1),
          axis.text.y = element_text(size = BASE - 1.5), panel.spacing.y = unit(0.6, "mm"), strip.background = element_blank())
  fx7 <- (pa / pb / pc) + plot_layout(heights = c(1, 0.6, 2.1)) + plot_annotation(tag_levels = "A")
  save_fig(fx7, "FigX7_specification_curve", W2, 165)
  legends["FigX7_specification_curve"] <- paste(sprintf("Specification curve (multiverse) analysis of %d analytic specifications", nrow(d)),
    "(missing data x functional form x estimation x added predictor x outcome x population). Each specification was evaluated by internal-external",
    "cross-validation across clusters with at least 10 deaths; (A) pooled C statistic with 95% CI and (B) pooled calibration slope, ordered by C;",
    "(C) the analytic choices of each specification. The reference specification (orange) uses heuristic shrinkage and imputations from the full data,",
    "so its C differs slightly from the primary analysis.")
}

# =====================================================================
# FigX8 hospital type and case definition (heat map)
# =====================================================================
hv <- rbind(ex$v3, ex$v5)
if (!is.null(hv)) hv <- hv[stats::ave(seq_len(nrow(hv)), hv$group, FUN = length) > 1, ]   # groupings with only one level are not informative
if (!is.null(hv) && nrow(hv)) {
  msg("FigX8 hospital type and case definition")
  hv$row <- sprintf("%s (%s/%s)", sub(">=", "\u2265", hv$level), hv$events, prettyNum(hv$n, big.mark = ","))
  long <- rbind(data.frame(group = hv$group, row = hv$row, metric = "C statistic", txt = txt_ci(hv$C, hv$C_lo, hv$C_hi), score = hv$C),
                data.frame(group = hv$group, row = hv$row, metric = "O/E ratio", txt = txt_ci(hv$OE, hv$OE_lo, hv$OE_hi), score = -abs(log(hv$OE))),
                data.frame(group = hv$group, row = hv$row, metric = "Calibration slope", txt = txt_ci(hv$slope, hv$slope_lo, hv$slope_hi), score = -abs(hv$slope - 1)))
  long$group <- factor(long$group, levels = unique(hv$group))
  long$rel <- stats::ave(long$score, long$metric, FUN = function(z) if (diff(range(z)) > 0) (z - min(z)) / diff(range(z)) else 0.5)
  long$row <- factor(long$row, levels = rev(unique(hv$row))); long$metric <- factor(long$metric, c("C statistic", "O/E ratio", "Calibration slope"))
  fx8 <- ggplot(long, aes(metric, row, fill = rel)) + geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = txt, colour = rel > 0.6), size = (BASE - 1.5) / .pt, family = FIG_FONT) +
    scale_fill_gradient(low = "#F7F4E9", high = "#1F4E79", breaks = c(0, 1), labels = c("Worse", "Better"), name = "Relative performance\n(within metric)") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") + scale_x_discrete(position = "top") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") + labs(x = NULL, y = NULL) + theme_journal() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", legend.title = element_text(size = BASE - 1),
          legend.key.width = unit(8, "mm"), axis.text.x.top = element_text(face = "bold"), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"), strip.background = element_blank(), panel.spacing.y = unit(1.5, "mm"))
  save_fig(fx8, "FigX8_hospital_type_case_definition", W2, 20 + 7 * length(unique(hv$row)))
  legends["FigX8_hospital_type_case_definition"] <- paste("Performance of the full model by hospital characteristics (eICU-CRD) and by case definition (held-out predictions;",
    "Rubin's rules across imputations). Labels give deaths/patients; shading ranks performance within each metric (closeness to 1 for O/E ratio and slope).")
}

# =====================================================================
# FigX9 model updating within each database (intercept-only and intercept-and-slope updating)
# =====================================================================
stk_path <- file.path(DIR_LOCAL, "03_iecv_stacked.rds")
if (file.exists(stk_path)) {
  msg("FigX9 model updating")
  o9 <- lapply(readRDS(stk_path), function(s) s[order(s$row_id), ])
  av <- o9[[1]][, c("row_id", "db", "y")]; av$lp <- rowMeans(sapply(o9, `[[`, "lp"))
  fcal <- function(y, p, n_grid = 80) { p <- pmin(pmax(p, 1e-6), 1 - 1e-6); lpp <- qlogis(p)
    f <- stats::glm(y ~ splines::ns(lpp, df = 3), family = stats::binomial())
    rg <- stats::quantile(p, c(0.01, 0.99)); g <- seq(rg[[1]], rg[[2]], length.out = n_grid)
    data.frame(pred = g, obs = plogis(stats::predict(f, newdata = data.frame(lpp = qlogis(g))))) }
  VER <- c("Original", "Intercept updated", "Intercept and slope updated")
  ann <- list()
  cur <- do.call(rbind, lapply(split(av, av$db), function(z) {
    g1 <- stats::glm(y ~ 1 + offset(lp), family = stats::binomial(), data = z); g2 <- stats::glm(y ~ lp, family = stats::binomial(), data = z)
    ci <- stats::confint.default(g2)["lp", ]
    ann[[z$db[1]]] <<- data.frame(db = z$db[1], txt = sprintf("Recalibration slope %s (%s to %s)\nIntercept update %s", f2(stats::coef(g2)[["lp"]]),
                                                             f2(ci[[1]]), f2(ci[[2]]), f2(stats::coef(g1)[[1]])))
    rbind(cbind(db = z$db[1], version = VER[1], fcal(z$y, plogis(z$lp))),
          cbind(db = z$db[1], version = VER[2], fcal(z$y, stats::fitted(g1))),
          cbind(db = z$db[1], version = VER[3], fcal(z$y, stats::fitted(g2)))) }))
  ann <- do.call(rbind, ann); cur$version <- factor(cur$version, VER); cur$dbl <- DB_LAB[cur$db]; ann$dbl <- DB_LAB[ann$db]
  mx <- min(0.6, max(cur$pred, cur$obs) * 1.05)
  fx9 <- ggplot(cur, aes(pred, obs, colour = version, linetype = version)) + geom_abline(linetype = "dashed", linewidth = 0.3, colour = "grey50") +
    geom_line(linewidth = 0.55) + facet_wrap(~ dbl) +
    geom_text(data = ann, aes(x = 0, y = mx, label = txt), inherit.aes = FALSE, hjust = 0, vjust = 1, size = (BASE - 1.5) / .pt,
              family = FIG_FONT, lineheight = 0.9) +
    scale_colour_manual(values = c("grey55", COL_NEW, "#D55E00"), name = NULL) + scale_linetype_manual(values = c("solid", "solid", "longdash"), name = NULL) +
    scale_x_continuous(labels = pct_lab) + scale_y_continuous(labels = pct_lab) + coord_cartesian(xlim = c(0, mx), ylim = c(0, mx)) +
    labs(x = "Predicted risk (held-out)", y = "Observed proportion") + theme_journal() +
    theme(legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(face = "bold", size = BASE))
  save_fig(fx9, "FigX9_model_updating", W15, 80)
  legends["FigX9_model_updating"] <- paste("Updating the model within each database (held-out predictions from internal-external cross-validation, averaged",
    "across imputations). Flexible calibration curves for the original predictions, after updating the intercept only, and after logistic",
    "recalibration of intercept and slope. A recalibration slope below 1 indicates predictions that are too extreme and above 1 too moderate.")
}

writeLines(paste0(names(legends), ": ", legends, "\n"), file.path(DIR_FIG, "FigX_legends.txt"))
msg("11 done")
