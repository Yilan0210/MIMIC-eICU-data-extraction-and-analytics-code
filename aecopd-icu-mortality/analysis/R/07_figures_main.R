# =====================================================================
# 07_figures_main.R: main figures; journal style, PDF + 600 dpi PNG (figure numbers here are working numbers; see README)
#   Figure 1  study flow                     Figure 2  predictor effects and importance
#   Figure 3  transportability (heat map, clustered calibration, case mix)  Figure 4  comparison with existing scores (C, decision curves)
#   Figure 5  nomogram (simplified model)    Figure 6  risk stratification and clinical utility
#   Figure 7  markers (secondary questions Q3/Q4)
# Requires the results of 01–05; run this script alone to redraw the figures
# =====================================================================
source("00_setup.R"); source("fig_style.R")
SCORES <- c("BAP65", "DECAF_mod", "APACHE_IVa", "SAPSII", "SOFA")
set.seed(SEED)
d02 <- readRDS(file.path(DIR_LOCAL, "02_dev.rds")); core <- d02$core; simp <- d02$simp; knots <- d02$knots
r03 <- readRDS(file.path(DIR_LOCAL, "03_iecv_results.rds")); CL <- r03$cl_label
stk <- readRDS(file.path(DIR_LOCAL, "03_iecv_stacked.rds"))
r05 <- readRDS(file.path(DIR_LOCAL, "05_markers.rds"))
a <- readRDS(file.path(DIR_LOCAL, "01_population.rds"))$a
legends <- character()
pct_lab <- function(x) paste0(round(100 * x), "%")
p_txt <- function(p) if (p < 0.001) "<0.001" else formatC(p, format = "f", digits = 3)
COLS <- c(full = "#0072B2", simplified = "#56B4E9", BAP65 = "#7F7F7F", DECAF_mod = "#CC79A7", APACHE_IVa = "#009E73",
          SAPSII = "#E69F00", SOFA = "#D55E00", treat_all = "#B0B0B0", treat_none = "#000000")
NAMES <- c(full = "New model (full)", simplified = "New model (simplified)", BAP65 = "BAP-65", DECAF_mod = "Modified DECAF",
           APACHE_IVa = "APACHE IVa", SAPSII = "SAPS II", SOFA = "SOFA", treat_all = "Treat all", treat_none = "Treat none")
COL_MIMIC <- "#0072B2"; COL_EICU <- "#E69F00"; COL_DIED <- "#D55E00"; COL_SURV <- "#56B4E9"
ord <- lapply(stk, function(s) s[order(s$row_id), ])
avg <- ord[[1]][, c("row_id", "label", "db", "y")]; names(avg)[2] <- "cluster"
for (m in c("full", "simplified", SCORES)) avg[[m]] <- rowMeans(sapply(ord, `[[`, m))
nbf <- function(y, p, pt) mean(p >= pt & y == 1) - mean(p >= pt & y == 0) * pt / (1 - pt)

# =====================================================================
# Figure 1: study flow
# =====================================================================
msg("Figure 1")
raw <- load_clean(); raw <- raw[raw$pop == "main", ]
nz <- function(z) { z <- suppressWarnings(as.numeric(z)); !is.na(z) & z == 1 }
flow_db <- function(db) {
  r <- raw[raw$db == db, ]; dth <- function(k) sum(nz(r$death[k]))
  k0 <- rep(TRUE, nrow(r)); steps <- list(
    c("Age <18 years", "incl_age", 0), c("ICU stay <24 h", "incl_los24", 0), c("Outcome unknown", "incl_outcome", 0),
    c("COVID-19", "excl_covid", 1), c("Haematological cancer", "excl_heme", 1), c("HIV/AIDS", "excl_hiv", 1),
    c("Transplant recipient", "excl_tx", 1), c("Cirrhosis", "excl_cirr", 1))
  ex1 <- character(); k <- k0
  for (st in steps) { bad <- if (st[3] == "0") !nz(r[[st[2]]]) else nz(r[[st[2]]]); n <- sum(k & bad); k <- k & !bad
    ex1 <- c(ex1, sprintf("%s: %s", st[1], format(n, big.mark = ","))) }
  kt <- k & nz(r$has_alb_qc) & nz(r$has_diff_qc)
  list(n0 = nrow(r), d0 = dth(k0), ex1 = ex1, nex1 = nrow(r) - sum(k), n1 = sum(k), d1 = dth(k), nt = sum(kt), dt = dth(kt),
       hosp = length(unique(r$hospital_id[k])))
}
fm_ <- flow_db("MIMIC"); fe_ <- flow_db("eICU")
cl_n <- tapply(a$death, a$cluster, length); cl_d <- tapply(a$death, a$cluster, sum)
cl_txt <- function(ids) paste(sprintf("%s %s/%s", sub("^(MIMIC|eICU) ", "", CL[ids]), cl_d[ids], prettyNum(cl_n[ids], big.mark = ",")), collapse = "; ")
boxes <- list(); arrows <- list()
bx <- function(x0, x1, y0, y1, lab, fill = "white", hjust = 0.5) boxes[[length(boxes) + 1]] <<- data.frame(x0 = x0, x1 = x1, y0 = y0, y1 = y1, fill = fill, lab = lab, hj = hjust)
ar <- function(x0, y0, x1, y1, lt = "solid") arrows[[length(arrows) + 1]] <<- data.frame(x0 = x0, y0 = y0, x1 = x1, y1 = y1, lt = lt)
fmt <- function(z) format(z, big.mark = ",")
col_flow <- function(f, name, xm0, xm1, xe0, xe1) {
  xc <- (xm0 + xm1) / 2
  bx(xm0, xm1, 132, 150, sprintf("%s\nFirst ICU admission for AECOPD\nn = %s (%s deaths)", name, fmt(f$n0), fmt(f$d0)), fill = "#EAF2F8")
  bx(xe0, xe1, 97, 131, paste(c(sprintf("Excluded (n = %s)", fmt(f$nex1)), paste0("  ", f$ex1)), collapse = "\n"), hjust = 0)
  bx(xm0, xm1, 70, 86, sprintf("Eligible and analysed\nn = %s (%s deaths)%s", fmt(f$n1), fmt(f$d1), if (f$hosp > 1) sprintf("\n%d hospitals", f$hosp) else ""), fill = "#EAF2F8")
  bx(xe0, xe1, 70, 86, sprintf("Marker subset*\nn = %s (%s deaths)", fmt(f$nt), fmt(f$dt)), fill = "#F4F6F6")
  ar(xc, 132, xc, 86); ar(xc, 114, xe0, 114); ar(xm1, 78, xe0, 78, "dashed"); ar(xc, 70, xc, 58)
}
col_flow(fm_, "MIMIC-IV v3.1\n(single centre, 2008–2022)", 2, 44, 48, 90)
col_flow(fe_, "eICU-CRD v2.0\n(208 hospitals, 2014–2015)", 93, 135, 139, 181)
bx(16, 121, 44, 58, sprintf("Development cohort\nn = %s (%s deaths, %.1f%%)", fmt(fm_$n1 + fe_$n1), fmt(fm_$d1 + fe_$d1),
   100 * (fm_$d1 + fe_$d1) / (fm_$n1 + fe_$n1)), fill = "#FDEBD0")
ar(68.5, 44, 68.5, 32)
bx(4, 179, 8, 32, paste0("Internal-external cross-validation: 7 clusters, each held out in turn (deaths/n)\n",
   "MIMIC-IV coding eras: ", cl_txt(c("1", "2")), "\n", "eICU-CRD regions: ", cl_txt(c("3", "4", "5")), "\n", cl_txt(c("6", "7"))), fill = "#F4F6F6")
bx(4, 179, 0.5, 6, "*Albumin and differential count measured; used for the secondary analyses of immune-inflammation-nutrition markers", fill = NA)
B <- do.call(rbind, boxes); A <- do.call(rbind, arrows); B$tx <- ifelse(B$hj == 0, B$x0 + 2, (B$x0 + B$x1) / 2)
f1 <- ggplot() + geom_rect(data = B[!is.na(B$fill), ], aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = I(fill)), colour = "grey25", linewidth = 0.3) +
  geom_text(data = B, aes(x = tx, y = (y0 + y1) / 2, label = lab, hjust = hj, size = ifelse(hj == 0 | is.na(fill), BASE - 1, BASE - 0.5) / .pt),
            family = FIG_FONT, lineheight = 1.05) + scale_size_identity() +
  geom_segment(data = A, aes(x = x0, y = y0, xend = x1, yend = y1, linetype = I(lt)), linewidth = 0.3, arrow = grid::arrow(length = unit(1.4, "mm"), type = "closed")) +
  scale_x_continuous(limits = c(0, 183), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 152), expand = c(0, 0)) + coord_fixed() + theme_void()
save_fig(f1, "Figure1_study_flow", W2, 152)
legends["Figure1_study_flow"] <- paste("Study flow. All eligible patients formed the development cohort; the marker subset (albumin and",
  "differential count measured) was used only for the secondary analyses. Internal-external cross-validation held out each of seven clusters in turn.")

# =====================================================================
# Figure 2: predictor effects and importance
# =====================================================================
msg("Figure 2")
imp <- d02$imp_tab; imp$lab <- VAR_LAB[imp$var]
imp$kind <- ifelse(imp$var %in% SPLINE_VARS, "Spline (2 df)", "Linear or binary (1 df)")
imp$lab <- factor(imp$lab, levels = imp$lab[order(imp$chi2_minus_df)])
pA <- ggplot(imp, aes(chi2_minus_df, lab)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_segment(aes(x = 0, xend = chi2_minus_df, yend = lab), linewidth = 0.4, colour = "grey55") +
  geom_point(aes(shape = kind), size = 1.8, colour = COLS[["full"]], fill = "white") +
  scale_shape_manual(values = c("Spline (2 df)" = 21, "Linear or binary (1 df)" = 16)) +
  scale_x_continuous(labels = lab_minus) + labs(x = "Wald χ² minus df", y = NULL) + theme_journal() +
  theme(legend.position = "bottom", legend.justification = "left", legend.direction = "vertical", legend.title = element_blank())
x1 <- d02$xs[[1]]
or_lab <- function(z) sub("\\.$", "", sub("0+$", "", formatC(z, format = "f", digits = 2)))
eff_vars <- c("age", "ln_bun", "gcs", "hr", "rr", "map", "temp", "ph", "paco2")
cv <- d02$curves
unit_lab <- function(v) { u <- VAR_UNIT[v]; if (is.na(u) || u == "") VAR_LAB[v] else paste0(VAR_LAB[v], ", ", u) }
yr <- range(exp(c(cv$lo[cv$var %in% eff_vars], cv$hi[cv$var %in% eff_vars])))
yr <- c(max(yr[1], 0.1), min(yr[2], 10))
eff_panel <- function(v) {
  d <- cv[cv$var == v, ]; xr <- to_raw(v, x1[[v]]); xr <- xr[xr >= min(d$x_raw) & xr <= max(d$x_raw)]
  dz <- stats::density(xr, n = 256); dd <- data.frame(x = dz$x, y = exp(log(yr[1]) + dz$y / max(dz$y) * 0.2 * diff(log(yr))))
  dd <- dd[dd$x >= min(d$x_raw) & dd$x <= max(d$x_raw), ]
  ggplot() + geom_ribbon(data = dd, aes(x, ymin = yr[1], ymax = y), fill = "grey88") +
    geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey35") +
    geom_ribbon(data = d, aes(x_raw, ymin = exp(lo), ymax = exp(hi)), fill = COLS[["full"]], alpha = 0.18) +
    geom_line(data = d, aes(x_raw, exp(logor)), linewidth = 0.55, colour = COLS[["full"]]) +
    scale_y_log10(labels = or_lab, breaks = c(0.25, 0.5, 1, 2, 4, 8)) + coord_cartesian(ylim = yr, expand = FALSE) +
    scale_x_continuous(labels = lab_minus, breaks = function(l) pretty(l, n = 4)) + labs(x = unit_lab(v), y = NULL) + theme_journal() + theme(axis.title.x = element_text(size = BASE - 0.5))
}
pB <- wrap_plots(lapply(eff_vars, eff_panel), ncol = 3) & theme(plot.margin = margin(1, 2, 1, 1, "mm"))
pB <- wrap_elements(full = pB + plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))) +
  labs(tag = NULL)
ylab_g <- wrap_elements(full = grid::textGrob("Adjusted odds ratio (95% CI) vs median", rot = 90, gp = grid::gpar(fontsize = BASE, fontfamily = FIG_FONT)))
f2 <- wrap_elements(full = pA) + ylab_g + pB + plot_layout(widths = c(0.9, 0.05, 1.9)) + plot_annotation(tag_levels = list(c("A", "B", "")))
save_fig(f2, "Figure2_predictor_effects", W2, 112)
legends["Figure2_predictor_effects"] <- paste("Predictors of in-hospital death in the core model. (A) Importance: pooled Wald χ² minus degrees of",
  "freedom (both spline terms tested jointly). (B) Partial effects of continuous predictors (odds ratio relative to the median, other predictors",
  "held constant) from restricted cubic splines with three knots; shaded bands are 95% CIs combined across imputations; grey areas show the",
  "distribution of each predictor. Curves are restricted to the 1st–99th percentiles. Odds ratios for all predictors are in Table 2.")

# =====================================================================
# Figure 3: transportability
# =====================================================================
msg("Figure 3")
bc <- r03$by_cluster; po <- r03$pooled[r03$pooled$clusters == "all eligible", ]
MET <- c(logit_C = "C statistic", log_OE = "O/E ratio", cal_slope = "Calibration slope", CITL = "Calibration-in-the-large")
hm <- rbind(
  do.call(rbind, lapply(names(MET), function(me) { b <- bc[bc$measure == me & bc$model == "full", ]
    data.frame(row = paste0(b$label, "  ", CL[b$label], " (", dn_txt(b$deaths, b$n), ")"), ord = as.integer(b$label), metric = me,
               value = b$est, txt = num(b$est, 2)) })),
  do.call(rbind, lapply(names(MET), function(me) { p <- po[po$measure == me & po$model == "full", ]
    rbind(data.frame(row = "Pooled (random effects)", ord = 99, metric = me, value = p$est, txt = num(p$est, 2)),
          data.frame(row = "95% prediction interval", ord = 100, metric = me, value = NA, txt = paste0(num(p$pi_lo, 2), " to ", num(p$pi_hi, 2)))) })))
hm$score <- ifelse(hm$metric == "logit_C", hm$value, ifelse(hm$metric %in% c("log_OE", "cal_slope"), -abs(log(pmax(hm$value, 0.01))), -abs(hm$value)))
hm$rel <- ave(hm$score, hm$metric, FUN = function(s) { r <- range(s, na.rm = TRUE); (s - r[1]) / max(diff(r), 1e-9) })
hm$metric <- factor(MET[hm$metric], levels = MET); hm$row <- factor(hm$row, levels = unique(hm$row[order(-hm$ord)]))
pA <- ggplot(hm, aes(metric, row)) + geom_tile(aes(fill = rel), colour = "white", linewidth = 0.6) +
  geom_text(aes(label = txt, colour = !is.na(rel) & rel > 0.6, fontface = ifelse(ord >= 99, "bold", "plain")), size = (BASE - 1) / .pt, family = FIG_FONT) +
  geom_hline(yintercept = 2.5, linewidth = 0.4) +
  scale_fill_gradient(low = "#F7F4EA", high = "#1F4E79", limits = c(0, 1), breaks = c(0, 1), labels = c("Worse", "Better"), na.value = "white",
                      name = "Relative performance\n(within metric)") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_x_discrete(position = "top", labels = c("C statistic" = "C\nstatistic", "O/E ratio" = "O/E\nratio", "Calibration slope" = "Calibration\nslope",
                                                 "Calibration-in-the-large" = "CITL")) + labs(x = NULL, y = NULL) + theme_journal() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x.top = element_text(size = BASE - 0.5, face = "bold"),
        axis.text.y = element_text(size = BASE - 1), legend.position = "bottom", legend.title = element_text(size = BASE - 1),
        legend.key.width = unit(8, "mm"), legend.key.height = unit(2, "mm"), legend.margin = margin(0, 6, 0, 0, "mm")) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))   # legend title above the legend so that "Better" is not clipped
grid_p <- seq(0.02, 0.80, by = 0.01)
cl_curves <- do.call(rbind, lapply(split(avg, avg$cluster), function(s) {
  lp <- qlogis(pmin(pmax(s$full, 1e-4), 1 - 1e-4)); kn <- stats::quantile(lp, c(.1, .5, .9))
  Bm <- Hmisc::rcspline.eval(lp, knots = kn, inclx = TRUE); f <- stats::glm(s$y ~ Bm, family = stats::binomial())
  rg <- stats::quantile(s$full, c(.05, .95)); g <- grid_p[grid_p >= rg[1] & grid_p <= rg[2]]
  X <- cbind(1, Hmisc::rcspline.eval(qlogis(g), knots = kn, inclx = TRUE))
  data.frame(cluster = s$cluster[1], db = s$db[1], pred = g, eta = drop(X %*% stats::coef(f)), se = sqrt(rowSums((X %*% stats::vcov(f)) * X))) }))
K_CL <- length(unique(cl_curves$cluster))
pool_cal <- do.call(rbind, lapply(split(cl_curves, cl_curves$pred), function(s) {
  if (nrow(s) < max(3, ceiling(0.8 * K_CL))) return(NULL)
  r <- tryCatch(metafor::rma(yi = s$eta, sei = s$se, method = "REML", test = "knha"), error = function(e) metafor::rma(yi = s$eta, sei = s$se, method = "PM", test = "knha"))
  q <- stats::qt(0.975, nrow(s) - 2); sp <- sqrt(r$tau2 + r$se^2)
  data.frame(pred = s$pred[1], est = plogis(r$beta[1]), lo = plogis(r$ci.lb), hi = plogis(r$ci.ub), pi_lo = plogis(r$beta[1] - q * sp), pi_hi = plogis(r$beta[1] + q * sp)) }))
xm <- min(0.8, ceiling((max(pool_cal$pred) + 0.03) * 20) / 20)
pB <- ggplot() + geom_ribbon(data = pool_cal, aes(pred, ymin = pi_lo, ymax = pi_hi), fill = "grey87") +
  geom_ribbon(data = pool_cal, aes(pred, ymin = lo, ymax = hi), fill = "grey62") +
  geom_line(data = cl_curves, aes(pred, plogis(eta), group = cluster, colour = db), linewidth = 0.35, alpha = 0.9) +
  geom_line(data = pool_cal, aes(pred, est), linewidth = 0.8, colour = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey20") +
  scale_colour_manual(values = c(MIMIC = COL_MIMIC, eICU = COL_EICU), labels = c(MIMIC = "MIMIC-IV clusters", eICU = "eICU-CRD clusters")) +
  scale_x_continuous(labels = pct_lab, expand = c(0, 0)) + scale_y_continuous(labels = pct_lab, expand = c(0, 0)) +
  coord_fixed(xlim = c(0, xm), ylim = c(0, xm), clip = "on") + labs(x = "Predicted risk", y = "Observed proportion") + theme_journal() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
cm <- r03$case_mix
pC <- ggplot(cm, aes(sd_lp, C)) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.3, colour = "grey60", linetype = "dashed") +
  geom_point(aes(size = n, colour = ifelse(as.integer(as.character(cluster)) <= 2, "MIMIC", "eICU")), alpha = 0.85) +
  geom_text(aes(label = cluster), size = (BASE - 1.5) / .pt, family = FIG_FONT, colour = "white") +
  scale_colour_manual(values = c(MIMIC = COL_MIMIC, eICU = COL_EICU), guide = "none") + scale_size_area(max_size = 5, guide = "none") +
  scale_x_continuous(labels = lab_minus) + scale_y_continuous(labels = lab_minus) +
  labs(x = "SD of linear predictor in held-out cluster\n(case-mix heterogeneity)", y = "C statistic in held-out cluster") + theme_journal()
f3 <- wrap_elements(full = pA) + (wrap_elements(full = pB) / wrap_elements(full = pC) + plot_layout(heights = c(1.35, 1))) +
  plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "A")
save_fig(f3, "Figure3_transportability", W2, 150)
legends["Figure3_transportability"] <- paste("Transportability across hospitals, regions and coding eras (internal-external cross-validation,",
  "full model). (A) Performance in each held-out cluster, the random-effects pooled estimate and its 95% prediction interval; shading ranks",
  "performance within each metric (closeness to 1 for O/E ratio and calibration slope, to 0 for calibration-in-the-large [CITL]). (B) Clustered",
  "flexible calibration: cluster-specific curves, the pooled curve (black) with its 95% CI (dark grey) and 95% prediction interval (light grey)",
  "from a two-stage random-effects meta-analysis (Barreñada et al., 2025), shown where at least 80% of clusters contribute. (C) Cluster C",
  "statistic against case-mix heterogeneity (SD of the linear predictor); point size is proportional to cluster size.")

# =====================================================================
# Figure 4: comparison with existing scores
# =====================================================================
msg("Figure 4")
P <- r03$pooled
db_rows <- do.call(rbind, lapply(SCORES, function(s) {
  cs <- P[P$measure == "logit_C" & P$model == s, ]; cf <- P[P$measure == "logit_C_same_patients" & P$model == paste0("full|", s), ]
  dd <- P[P$measure == "dC_vs_score" & P$model == paste0("full-", s), ]
  if (!nrow(cs) || !nrow(cf)) return(NULL)
  data.frame(score = s, C_score = cs$est, C_model = cf$est, dC = dd$est, lo = dd$lo, hi = dd$hi, k = cs$k) }))
db_rows$lab <- factor(paste0(NAMES[db_rows$score], ifelse(db_rows$k < 3, " †", "")), levels = rev(paste0(NAMES[db_rows$score], ifelse(db_rows$k < 3, " †", ""))))
db_rows$txt <- paste0("ΔC ", ci_txt(db_rows$dC, db_rows$lo, db_rows$hi, 3))
xlo <- floor(min(db_rows$C_score) * 20) / 20; xhi <- ceiling(max(db_rows$C_model) * 20) / 20
pA <- ggplot(db_rows, aes(y = lab)) + geom_segment(aes(x = C_score, xend = C_model, yend = lab), linewidth = 0.6, colour = "grey70") +
  geom_point(aes(x = C_score, colour = score), size = 2.2, shape = 21, fill = "white", stroke = 0.8) +
  geom_point(aes(x = C_model), size = 2.2, colour = COLS[["full"]]) +
  geom_text(aes(x = xhi + 0.005, label = txt), hjust = 0, size = (BASE - 1) / .pt, family = FIG_FONT) +
  scale_colour_manual(values = COLS, guide = "none") + scale_x_continuous(labels = lab_minus, breaks = seq(0.5, 0.9, 0.05)) +
  coord_cartesian(xlim = c(xlo, xhi + 0.12), clip = "off") + labs(x = "C statistic (pooled across held-out clusters)", y = NULL) +
  annotate("point", x = xlo + 0.005, y = nrow(db_rows) + 0.75, size = 2, shape = 21, fill = "white") +
  annotate("text", x = xlo + 0.011, y = nrow(db_rows) + 0.75, label = "Score", hjust = 0, size = (BASE - 1) / .pt, family = FIG_FONT) +
  annotate("point", x = xlo + 0.05, y = nrow(db_rows) + 0.75, size = 2, colour = COLS[["full"]]) +
  annotate("text", x = xlo + 0.056, y = nrow(db_rows) + 0.75, label = "New model in the same patients", hjust = 0, size = (BASE - 1) / .pt, family = FIG_FONT) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.1))) +
  theme_journal() + theme(plot.margin = margin(2, 30, 2, 2, "mm"))
dca_sub <- function(models, rows_fun) {
  do.call(rbind, lapply(DCA_GRID, function(pt) { nb <- sapply(stk, function(s) { k <- rows_fun(s)
    c(vapply(models, function(m) nbf(s$y[k], s[[m]][k], pt), numeric(1)), treat_all = nbf(s$y[k], rep(1, sum(k)), pt)) })
    data.frame(threshold = pt, t(rowMeans(nb))) }))
}
dca_panel <- function(d, models, title) {
  l <- do.call(rbind, lapply(c(models, "treat_all"), function(m) data.frame(threshold = d$threshold, nb = d[[m]], curve = m)))
  l$curve <- factor(l$curve, levels = c(models, "treat_all")); ymax <- max(d[, models], na.rm = TRUE) * 1.1
  ggplot(l, aes(threshold, nb, colour = curve, linetype = curve)) + annotate("rect", xmin = 0.10, xmax = 0.30, ymin = -Inf, ymax = Inf, fill = "grey94") +
    geom_hline(yintercept = 0, linewidth = 0.3) + geom_line(linewidth = 0.55) +
    scale_colour_manual(values = COLS, labels = NAMES) +
    scale_linetype_manual(values = c(full = "solid", simplified = "solid", BAP65 = "longdash", DECAF_mod = "longdash", APACHE_IVa = "longdash",
                                     SAPSII = "longdash", SOFA = "dotdash", treat_all = "solid"), labels = NAMES) +
    scale_x_continuous(labels = pct_lab, breaks = seq(0.1, 0.5, 0.1)) + scale_y_continuous(labels = lab_minus) +
    coord_cartesian(ylim = c(-0.01, ymax)) + labs(x = "Threshold probability", y = "Net benefit", subtitle = title) + theme_journal() +
    theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1), legend.text = element_text(size = BASE - 1.5),
          legend.key.width = unit(4.5, "mm"), plot.subtitle = element_text(size = BASE, face = "bold"))
}
pB <- dca_panel(dca_sub(c("full", "simplified", "BAP65", "DECAF_mod"), function(s) rep(TRUE, nrow(s))), c("full", "simplified", "BAP65", "DECAF_mod"), "All patients")
pC <- dca_panel(dca_sub(c("full", "APACHE_IVa"), function(s) !is.na(s$APACHE_IVa)), c("full", "APACHE_IVa"), "eICU-CRD (APACHE IVa available)")
pD <- dca_panel(dca_sub(c("full", "SAPSII", "SOFA"), function(s) !is.na(s$SAPSII)), c("full", "SAPSII", "SOFA"), "MIMIC-IV")
f4 <- wrap_elements(full = pA) / (pB | pC | pD) + plot_layout(heights = c(0.8, 1)) + plot_annotation(tag_levels = "A")
save_fig(f4, "Figure4_comparison_scores", W2, 130)
legends["Figure4_comparison_scores"] <- paste("Comparison with existing scores in the same held-out patients. (A) Pooled C statistic of each score",
  "(open circle) and of the new model in the same patients (filled circle), with the pooled paired difference (ΔC, DeLong). † two",
  "clusters only (MIMIC-IV): pooled value for reference. Scores other than APACHE IVa were converted to risks by logistic recalibration in the",
  "development folds; modified DECAF lacks the dyspnoea grade (eMRCD). (B–D) Decision curves of held-out predictions in all patients and",
  "in the subsets with APACHE IVa (eICU-CRD) or SAPS II and SOFA (MIMIC-IV); the shaded area marks the prespecified thresholds (10%–30%).")

# =====================================================================
# Figure 5: nomogram (simplified model)
# =====================================================================
msg("Figure 5")
nd <- data.frame(Age = x1$age, BUN = pmax(x1$bun, 1), HR = x1$hr, RR = x1$rr, MAP = x1$map, Temp = x1$temp, pH = x1$ph, PaCO2 = x1$paco2,
                 Female = x1$female, GCS_motor = x1$gcs, Hb = x1$hgb, Platelets = pmax(x1$plt, 1), WBC = pmax(x1$wbc, 0.1),
                 Glucose = pmax(x1$glucose, 10), AF = x1$af, IMV = x1$imv, death = x1$death)
NM <- c(age = "Age", ln_bun = "BUN", hr = "HR", rr = "RR", map = "MAP", temp = "Temp", ph = "pH", paco2 = "PaCO2", female = "Female",
        gcs = "GCS_motor", hgb = "Hb", ln_plt = "Platelets", ln_wbc = "WBC", ln_glu = "Glucose", af = "AF", imv = "IMV")
LOGV <- c("ln_bun", "ln_plt", "ln_wbc", "ln_glu")
term_txt <- function(v) { nm <- NM[[v]]; z <- if (v %in% LOGV) sprintf("log(%s)", nm) else nm
  if (v %in% SPLINE_VARS) sprintf("rcs(%s, c(%s))", z, paste(signif(knots[[v]], 10), collapse = ", ")) else z }
NLAB <- c(Age = "Age, years", BUN = "Blood urea nitrogen, mg/dL", HR = "Heart rate, beats/min", RR = "Respiratory rate, breaths/min",
          MAP = "Mean arterial pressure, mmHg", Temp = "Temperature, °C", pH = "pH", PaCO2 = "PaCO2, mmHg", Female = "Female sex",
          GCS_motor = "GCS motor score", Hb = "Haemoglobin, g/dL", Platelets = "Platelets, ×10\u2079/L", WBC = "White cells, ×10\u2079/L",
          Glucose = "Glucose, mg/dL", AF = "Atrial fibrillation within 24 h", IMV = "Invasive ventilation, day 1")
for (v in names(NLAB)) Hmisc::label(nd[[v]]) <- NLAB[[v]]
dd <- rms::datadist(nd); assign("dd", dd, envir = .GlobalEnv); options(datadist = "dd")
fo <- stats::as.formula(paste("death ~", paste(vapply(simp$vars, term_txt, character(1)), collapse = " + ")))
fn <- rms::lrm(fo, data = nd)
if (length(stats::coef(fn)) != length(simp$coef)) stop("Nomogram model terms differ from the simplified model")
fn$coefficients[] <- unname(simp$coef)
chk <- max(abs(stats::predict(fn, newdata = nd[1:50, ]) - drop(cbind(1, as.matrix(x1[1:50, simp$terms])) %*% simp$coef)))
if (chk > 1e-6) stop("Linear predictors of the nomogram and the simplified model differ: ", chk)
tk <- list(BUN = c(5, 10, 20, 40, 80, 160), Platelets = c(50, 100, 200, 400), WBC = c(2, 5, 10, 20, 40), Glucose = c(60, 100, 150, 250, 400))
tk <- tk[names(tk) %in% all.vars(fo)]
# axes of continuous variables drawn only to the 0.5–99.5th percentiles (extreme values stretch the point scale and compress other variables); ticks outside the range removed for log-scale variables
for (v in setdiff(intersect(all.vars(fo), names(nd)), c("death", "Female", "AF", "IMV"))) {
  q <- stats::quantile(nd[[v]], c(0.005, 0.995), na.rm = TRUE)
  if (v %in% names(tk)) { z <- tk[[v]]; tk[[v]] <- z[z >= q[1] * 0.8 & z <= q[2] * 1.25] }
  else if (v == "GCS_motor") tk[[v]] <- 1:6
  else { z <- pretty(q, n = 6); tk[[v]] <- z[z >= q[1] - diff(range(z)) / 12 & z <= q[2] + diff(range(z)) / 12] }
}
# Non-monotonic spline effects (e.g. respiratory rate reaching a plateau above about 25 breaths/min) make nomogram axes fold back with crowded ticks;
#       ticks are kept only from the lower end to the point of largest effect, and the legend states that higher values are scored at that point
nom_note <- character()
for (v in intersect(names(tk), NM[intersect(simp$vars, SPLINE_VARS)])) {
  z <- sort(tk[[v]]); if (length(z) < 3) next
  nd2 <- as.data.frame(lapply(nd, function(col) stats::median(col, na.rm = TRUE)))[rep(1, length(z)), , drop = FALSE]; nd2[[v]] <- z
  lz <- stats::predict(fn, newdata = nd2)
  k <- which.max(abs(lz - lz[1])); keep <- z[seq_len(k)]; lk <- lz[seq_len(k)]          # non-monotonic: keep ticks up to the point of largest effect
  thr <- 0.04 * diff(range(lz))                                                         # plateau: drop upper ticks whose effect differs by < 4% of the effect range
  while (length(keep) > 3 && abs(lk[length(lk)] - lk[length(lk) - 1]) < thr) { keep <- keep[-length(keep)]; lk <- lk[-length(lk)] }
  if (length(keep) >= 3 && length(keep) < length(z)) { tk[[v]] <- keep
    nom_note <- c(nom_note, sprintf("%s values above %s are scored as %s", sub(",.*$", "", NLAB[[v]]), max(keep), max(keep))) }
}
nom <- do.call(rms::nomogram, c(list(fn, fun = plogis, fun.at = c(0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9), funlabel = "Risk of in-hospital death", lp = FALSE), tk))
draw_nom <- function() { op <- graphics::par(family = if (FIG_FONT == "sans") "sans" else FIG_FONT, cex = 0.6, mar = c(1, 1, 1, 1))
  plot(nom, xfrac = 0.3, cex.axis = 0.85, cex.var = 1, col.grid = grDevices::gray(c(0.85, 0.95))); graphics::par(op) }
pdf_dev <- if (HAS_CAIRO) grDevices::cairo_pdf else grDevices::pdf
hh <- 30 + 9 * (length(simp$vars) + 3)
pdf_dev(file.path(DIR_FIG, "Figure5_nomogram.pdf"), width = W2 / 25.4, height = hh / 25.4); draw_nom(); invisible(grDevices::dev.off())
if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(file.path(DIR_FIG, "Figure5_nomogram.png"), width = W2, height = hh, units = "mm", res = 600, background = "white") else
  grDevices::png(file.path(DIR_FIG, "Figure5_nomogram.png"), width = W2, height = hh, units = "mm", res = 600, bg = "white")
draw_nom(); invisible(grDevices::dev.off()); msg("Figure: Figure5_nomogram (PDF + PNG)")
legends["Figure5_nomogram"] <- sprintf(paste("Nomogram of the simplified model (%d predictors; R² = %.3f against the full model's linear predictor).",
  "Points for each predictor are summed and the total is read on the risk axis. Continuous predictors enter as restricted cubic splines, so",
  "point scales may be non-monotonic. The full model is available as an online calculator."), length(simp$vars), simp$r2)
if (length(nom_note)) legends["Figure5_nomogram"] <- paste0(legends["Figure5_nomogram"], " Axes stop where the effect levels off: ",
  paste(nom_note, collapse = "; "), ".")

# =====================================================================
# Figure 6: risk stratification and clinical utility
# =====================================================================
msg("Figure 6")
ra <- r03$risk_all; rg <- r03$risk_groups
ra$group <- factor(ra$group, levels = c("<10%", "10–30%", "≥30%")); rg$group <- factor(rg$group, levels = levels(ra$group))
rg$db <- ifelse(as.integer(as.character(rg$cluster)) <= 2, "MIMIC", "eICU")
# diamonds (mean predicted risk) in black with a legend so they are not confused with the blue MIMIC cluster points; patient numbers below the x-axis labels
x_lab6 <- setNames(sprintf("%s\nn = %s (%.0f%%)", as.character(ra$group), prettyNum(ra$n, big.mark = ","), ra$pct_of_patients),
                   as.character(ra$group))
pA <- ggplot() + geom_col(data = ra, aes(group, observed, fill = "Observed (95% CI)"), width = 0.6) +
  geom_errorbar(data = ra, aes(group, ymin = obs_lo, ymax = obs_hi), width = 0.15, linewidth = 0.35) +
  geom_point(data = rg[rg$n >= 20, ], aes(group, observed, colour = db), position = position_jitter(width = 0.18, seed = 1), size = 1.1, alpha = 0.8) +
  geom_point(data = ra, aes(group, predicted, shape = "Mean predicted risk"), size = 2.6, fill = "black", colour = "white", stroke = 0.4) +
  scale_fill_manual(values = c("Observed (95% CI)" = "grey85"), name = NULL) +
  scale_shape_manual(values = c("Mean predicted risk" = 23), name = NULL) +
  scale_colour_manual(values = c(MIMIC = COL_MIMIC, eICU = COL_EICU), labels = c(MIMIC = "MIMIC-IV cluster", eICU = "eICU-CRD cluster"),
                      name = NULL) +
  guides(fill = guide_legend(order = 1, nrow = 1), shape = guide_legend(order = 2, nrow = 1, override.aes = list(fill = "black", colour = "white", size = 2.6)),
         colour = guide_legend(order = 3, nrow = 1, override.aes = list(size = 1.4))) +
  scale_x_discrete(labels = x_lab6) +
  scale_y_continuous(labels = pct_lab, limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Predicted risk group (held-out predictions)", y = "In-hospital mortality") + theme_journal() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(0.3, "mm"),
        legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size = BASE - 1))   # three rows so that the legend of panel A is not clipped
nf <- do.call(rbind, lapply(NB_THRESHOLDS, function(pt) { pos <- avg$full >= pt; y <- avg$y; n <- length(y) / 100
  data.frame(threshold = paste0(round(100 * pt), "%"), cat = c("True positive", "False positive", "False negative", "True negative"),
             per100 = c(sum(pos & y == 1), sum(pos & y == 0), sum(!pos & y == 1), sum(!pos & y == 0)) / n,
             sens = mean(pos[y == 1]), spec = mean(!pos[y == 0])) }))
nf$cat <- factor(nf$cat, levels = rev(c("True positive", "False positive", "False negative", "True negative")))
lab_nf <- unique(nf[, c("threshold", "sens", "spec")])
pB <- ggplot(nf, aes(per100, threshold, fill = cat)) + geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(per100 >= 4, sprintf("%.0f", per100), "")), position = position_stack(vjust = 0.5), size = (BASE - 1) / .pt, family = FIG_FONT) +
  geom_text(data = lab_nf, aes(x = 101, y = threshold, label = sprintf("Sens %.0f%%\nSpec %.0f%%", 100 * sens, 100 * spec)), inherit.aes = FALSE,
            hjust = 0, size = (BASE - 1.5) / .pt, family = FIG_FONT, lineheight = 0.9) +
  scale_fill_manual(values = c("True positive" = "#D55E00", "False positive" = "#F4B183", "False negative" = "#7F7F7F", "True negative" = "#D9D9D9"),
                    breaks = c("True positive", "False positive", "False negative", "True negative")) +
  coord_cartesian(xlim = c(0, 100), clip = "off") + labs(x = "Patients per 100 (held-out predictions)", y = "Risk threshold") + theme_journal() +
  theme(legend.position = "bottom", plot.margin = margin(2, 16, 2, 2, "mm")) + guides(fill = guide_legend(nrow = 2))
f6 <- pA + pB + plot_layout(widths = c(1, 1.2)) + plot_annotation(tag_levels = "A")
save_fig(f6, "Figure6_risk_stratification", W2, 85)
legends["Figure6_risk_stratification"] <- paste("Risk stratification with the full model (held-out predictions from internal-external cross-validation,",
  "averaged across imputations). (A) Observed in-hospital mortality (grey bars, 95% exact CI) and mean predicted risk (black diamonds) in the prespecified",
  "risk groups; dots show observed mortality in each cluster (groups with ≥20 patients). (B) Classification of 100 patients at the prespecified",
  "thresholds, with sensitivity and specificity.")

# =====================================================================
# Figure 7: markers (Q3/Q4)
# =====================================================================
msg("Figure 7")
dr <- r05$dr
MKLAB <- list(neut = u109("Neutrophils"), lymph = u109("Lymphocytes"), alb = "Albumin, g/L", eos = u109("Eosinophils"), nlr = "NLR", pni = "PNI")
dr_panel <- function(z) { d <- z$curve; yl <- range(c(d$lo, d$hi, 1)); yl <- exp(log(yl) + c(-0.05, 0.05) * diff(log(yl)))
  dd <- z$dens; dd <- dd[dd$x >= min(d$x) & dd$x <= max(d$x), ]; dd$yy <- exp(log(yl[1]) + dd$y * 0.22 * diff(log(yl)))
  p <- ggplot() + geom_ribbon(data = dd, aes(x, ymin = yl[1], ymax = yy), fill = "grey88") +
    geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey35") +
    geom_ribbon(data = d, aes(x, ymin = lo, ymax = hi), fill = COLS[["full"]], alpha = 0.18) + geom_line(data = d, aes(x, or), linewidth = 0.5, colour = COLS[["full"]]) +
    annotate("text", x = if (z$log) exp(mean(log(range(d$x)))) else mean(range(d$x)), y = yl[2], vjust = 1, size = (BASE - 1.5) / .pt, family = FIG_FONT,
             label = sprintf("P %s; non-linear %s", p_txt(z$p_overall), p_txt(z$p_nonlinear))) +
    scale_y_log10(labels = lab_minus) + coord_cartesian(ylim = yl, expand = FALSE) + labs(x = MKLAB[[z$var]], y = NULL) + theme_journal() +
    theme(axis.title.x = element_text(size = BASE - 0.5))
  if (z$log) p + scale_x_log10(labels = lab_minus) else p + scale_x_continuous(labels = lab_minus) }
pA <- wrap_plots(Map(function(z, i) dr_panel(z) + labs(y = if (i %in% c(1, 4)) "Adjusted OR" else NULL), dr, seq_along(dr)), nrow = 2)
ip <- r05$iecv_pool; ip$model <- sub(" vs base", "", ip$comparison); dl <- r05$delta; dl$model <- sub(" vs base", "", dl$comparison)
inc <- rbind(data.frame(model = dl$model, type = "Bootstrap-corrected", est = dl$dC_corrected, lo = dl$dC_corr_lo, hi = dl$dC_corr_hi),
             data.frame(model = ip$model, type = "IECV pooled", est = ip$est, lo = ip$lo, hi = ip$hi))
inc$lab <- factor(c(E1 = "+ N, L, albumin", E2 = "+ N, L, albumin,\n  eosinopenia", E3 = "+ splines of\n  N, L, albumin")[inc$model],
                  levels = rev(c("+ N, L, albumin", "+ N, L, albumin,\n  eosinopenia", "+ splines of\n  N, L, albumin")))
pB <- ggplot(inc, aes(est, lab, colour = type)) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = DELTA_C_THRESHOLD, linewidth = 0.3, linetype = "dashed", colour = "grey40") +
  geom_errorbar(orientation = "y", aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.45, position = position_dodge(width = 0.5)) +
  geom_point(size = 1.6, position = position_dodge(width = 0.5)) + scale_colour_manual(values = c("Bootstrap-corrected" = "#7F7F7F", "IECV pooled" = COLS[["full"]]), name = NULL) +
  scale_x_continuous(labels = lab_minus, breaks = function(l) seq(ceiling(l[1] / 0.02) * 0.02, floor(l[2] / 0.02) * 0.02, by = 0.02)) +
  expand_limits(x = c(-0.02, 0.04)) + labs(x = "ΔC over the core model", y = NULL) + theme_journal() +
  theme(legend.position = "bottom", legend.direction = "vertical", plot.subtitle = element_blank())
x1t <- d02$xs[[1]][d02$xs[[1]]$tested == 1, ]; gq <- r05$grid_ref$q; lpm <- r05$grid_ref$lp_median
cN <- r05$coef_N; cP <- r05$coef_P
gN <- exp(seq(log(gq$neut[1]), log(gq$neut[3]), length.out = 150)); gL <- exp(seq(log(max(gq$lymph[1], 0.05)), log(gq$lymph[3]), length.out = 150))
sA <- expand.grid(N = gN, L = gL); sA$risk <- plogis(cN[1] + cN["lp"] * lpm + cN["ln_n"] * log(sA$N) + cN["ln_l"] * log(sA$L) + cN["alb"] * gq$alb[2])
gA <- seq(gq$alb[1], gq$alb[3], length.out = 150); gL2 <- seq(gq$lymph[1], gq$lymph[3], length.out = 150)
sB <- expand.grid(A = gA, L = gL2); sB$risk <- plogis(cP[1] + cP["lp"] * lpm + cP["alb"] * sB$A + cP["lymph"] * sB$L)
rr <- range(c(sA$risk, sB$risk)); ibr <- pretty(rr, n = 7)
fillsc <- scale_fill_gradientn(colours = c("#FFF7EC", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"), labels = pct_lab, limits = rr, name = "Predicted risk",
                               guide = guide_colourbar(title.position = "top", barwidth = unit(30, "mm"), barheight = unit(2, "mm")))
lab_at <- function(d, f) if (nrow(d) >= 15) d[max(1, round(f * nrow(d))), ] else NULL
iso <- do.call(rbind, lapply(c(2, 5, 10, 20, 40), function(r) { d <- data.frame(N = gN, L = gN / r, r = r); d[d$L >= min(gL) & d$L <= max(gL), ] }))
isoB <- do.call(rbind, lapply(seq(25, 55, by = 5), function(v) { d <- data.frame(A = gA, L = (v - gA) / 5, v = v); d[d$L >= min(gL2) & d$L <= max(gL2), ] }))
lbl <- function(d, x, y, t) geom_label(data = d, aes(.data[[x]], .data[[y]], label = t), size = (BASE - 1.5) / .pt, family = FIG_FONT, fill = "white",
                                       label.size = 0, label.padding = unit(0.25, "mm"), label.r = unit(0, "mm"))
il <- do.call(rbind, lapply(split(iso, iso$r), lab_at, f = 0.7)); ilB <- do.call(rbind, lapply(split(isoB, isoB$v), lab_at, f = 0.7))
h2 <- r05$h2; ann <- function(i) { pt <- p_txt(h2$p[i]); sprintf("%s constraint %s\n95%% CI %s to %s; P %s", c("NLR", "PNI")[i], num(h2$est[i], 2), num(h2$lo[i], 2), num(h2$hi[i], 2), if (startsWith(pt, "<")) pt else paste("=", pt)) }
pC <- ggplot(sA, aes(N, L)) + geom_raster(aes(fill = risk)) + fillsc + geom_contour(aes(z = risk), breaks = ibr, colour = "white", linewidth = 0.3) +
  geom_line(data = iso, aes(N, L, group = r), linetype = "dashed", linewidth = 0.35, colour = "grey15") + lbl(il, "N", "L", paste0("NLR ", il$r)) +
  scale_x_log10(expand = c(0, 0), labels = lab_minus) + scale_y_log10(expand = c(0, 0), labels = lab_minus) +
  labs(x = u109("Neutrophils"), y = u109("Lymphocytes"), subtitle = ann(1)) + theme_journal() +
  theme(panel.border = element_rect(fill = NA, linewidth = 0.3), axis.line = element_blank(), plot.subtitle = element_text(size = BASE - 1.5))
pD <- ggplot(sB, aes(A, L)) + geom_raster(aes(fill = risk)) + fillsc + geom_contour(aes(z = risk), breaks = ibr, colour = "white", linewidth = 0.3) +
  geom_line(data = isoB, aes(A, L, group = v), linetype = "dashed", linewidth = 0.35, colour = "grey15") + lbl(ilB, "A", "L", paste0("PNI ", ilB$v)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + labs(x = "Albumin, g/L", y = u109("Lymphocytes"), subtitle = ann(2)) + theme_journal() +
  theme(panel.border = element_rect(fill = NA, linewidth = 0.3), axis.line = element_blank(), plot.subtitle = element_text(size = BASE - 1.5))
f7 <- wrap_elements(full = pA) / (pB + pC + pD + plot_layout(widths = c(0.8, 1, 1), guides = "collect") & theme(legend.position = "bottom", legend.title = element_text(size = BASE - 1))) +
  plot_layout(heights = c(1, 1.05)) + plot_annotation(tag_levels = "A")
save_fig(f7, "Figure7_markers", W2, 165)
legends["Figure7_markers"] <- paste(sprintf("Immune-inflammation-nutrition markers in the marker subset (n = %s; %s deaths).", format(r05$n, big.mark = ","), r05$deaths),
  "(A) Associations with in-hospital death after adjustment for the core model's linear predictor (restricted cubic splines, four knots; reference",
  "at the median; P values are multivariate Wald tests). (B) Gain in C statistic over the core model (bootstrap-corrected and pooled across held-out",
  "clusters); the dashed line marks the prespecified threshold of 0.02. (C, D) Predicted risk across neutrophil and lymphocyte counts, and across",
  "albumin and lymphocyte counts, with other predictors at their medians; white lines are iso-risk contours and dashed lines iso-lines of NLR and",
  "PNI (= albumin + 5 × lymphocytes). If the composite form held, iso-risk contours would run parallel to the iso-lines.")

writeLines(paste0(names(legends), ": ", legends, "\n"), file.path(DIR_FIG, "Figure_legends.txt"), useBytes = TRUE)
msg("07 done: main figures and draft legends in", DIR_FIG)
