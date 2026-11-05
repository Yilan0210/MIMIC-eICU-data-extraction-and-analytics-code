# =====================================================================
# 08_figures_supp.R: supplementary figures; journal style, PDF + PNG (working numbers; see README)
#   eFigure 1  C, O/E and calibration slope per IECV cluster (full model)     eFigure 2  dC of the simplified vs the full model
#   eFigure 3  calibration curves before and after updating in each database   eFigure 4  prediction stability
#   eFigure 5  sensitivity analyses A1–A13                                     eFigure 6  subgroups and fairness
#   eFigure 7  ROC curves                                                      eFigure 8  forest plots of marker dC per cluster
#   eFigure 9  missingness
# =====================================================================
source("00_setup.R"); source("fig_style.R")
r03 <- readRDS(file.path(DIR_LOCAL, "03_iecv_results.rds")); CL <- r03$cl_label
bc <- r03$by_cluster; po <- r03$pooled[r03$pooled$clusters == "all eligible", ]
r04 <- readRDS(file.path(DIR_LOCAL, "04_crossdb_update.rds")); r05 <- readRDS(file.path(DIR_LOCAL, "05_markers.rds"))
r06 <- readRDS(file.path(DIR_LOCAL, "06_sensitivity.rds")); st <- readRDS(file.path(DIR_LOCAL, "02_stability.rds"))
stk <- readRDS(file.path(DIR_LOCAL, "03_iecv_stacked.rds")); a <- readRDS(file.path(DIR_LOCAL, "01_population.rds"))$a
legends <- character(); pct_lab <- function(x) paste0(round(100 * x), "%")
hdr <- function(t) data.frame(label = t, sub = "", est = NA, lo = NA, hi = NA, type = "header")

# ---- eFigure 1: IECV clusters ----
iecv_rows <- function(meas, mm) {
  b <- bc[bc$measure == meas & bc$model == mm, ]; b <- b[order(as.integer(b$label)), ]; p <- po[po$measure == meas & po$model == mm, ]
  rbind(data.frame(label = paste0(b$label, "  ", CL[b$label]), sub = dn_txt(b$deaths, b$n), est = b$est, lo = b$lo, hi = b$hi, type = "row", pi_lo = NA, pi_hi = NA),
        data.frame(label = "Pooled (random effects)", sub = "", est = p$est, lo = p$lo, hi = p$hi, type = "pooled", pi_lo = p$pi_lo, pi_hi = p$pi_hi))
}
fa <- forest_plot(iecv_rows("logit_C", "full"), "C statistic", ref = 0.5, head = c("Held-out cluster", "Deaths/N", "C statistic (95% CI)"), tag = "A")
fb <- forest_plot(iecv_rows("log_OE", "full"), "Observed/expected ratio (log scale)", ref = 1, log = TRUE, head = c("Held-out cluster", "Deaths/N", "O/E ratio (95% CI)"), tag = "B")
fc <- forest_plot(iecv_rows("cal_slope", "full"), "Calibration slope", ref = 1, head = c("Held-out cluster", "Deaths/N", "Slope (95% CI)"), tag = "C")
save_fig(wrap_plots(list(wrap_elements(full = fa), wrap_elements(full = fb), wrap_elements(full = fc)), ncol = 1), "eFigure1_IECV_performance", W2, 3 * forest_height(10))
legends["eFigure1_IECV_performance"] <- paste("Internal-external cross-validation of the full model. Each cluster was held out in turn; the model was",
  "developed in the remaining clusters (multiple imputation, spline knots, maximum likelihood with bootstrap uniform shrinkage) and evaluated in the",
  "held-out cluster, imputed separately without the outcome. (A) C statistic; (B) observed/expected ratio; (C) calibration slope. Diamonds: random-",
  "effects pooled estimates (REML with Hartung-Knapp-Sidik-Jonkman adjustment); lines below: 95% prediction intervals (t with k−2 df).")

# ---- eFigure 2: simplified model ----
p <- forest_plot(iecv_rows("dC", "simplified-full"), "Difference in C statistic, simplified vs full model", ref = 0, digits = 3,
                 head = c("Held-out cluster", "Deaths/N", "ΔC (95% CI)"))
save_fig(p, "eFigure2_simplified_vs_full", W2, forest_height(10))
legends["eFigure2_simplified_vs_full"] <- "Loss of discrimination when the simplified model (used for the nomogram) replaces the full model, in each held-out cluster."

# ---- eFigure 3: model updating ----
cu <- r04$curves_upd
p <- ggplot(cu, aes(pred, obs, colour = version)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey30") +
  geom_line(linewidth = 0.6) + facet_wrap(~ db, labeller = labeller(db = c(MIMIC = "MIMIC-IV", eICU = "eICU-CRD"))) +
  scale_colour_manual(values = c("Original" = "#7F7F7F", "Intercept updated" = "#0072B2")) +
  scale_x_continuous(labels = pct_lab, breaks = seq(0, 0.8, 0.1)) + scale_y_continuous(labels = pct_lab, breaks = seq(0, 0.8, 0.1)) + coord_fixed() +
  labs(x = "Predicted risk (held-out)", y = "Observed proportion (loess)", colour = NULL) + theme_journal() +
  theme(legend.position = "bottom", panel.spacing = unit(7, "mm"))
save_fig(p, "eFigure3_model_updating", W15, 72)
up <- r04$upd
legends["eFigure3_model_updating"] <- paste("Calibration before and after updating the intercept within each database (held-out predictions from",
  sprintf("internal-external cross-validation). Intercept updates: %s.", paste(sprintf("%s %s (95%% CI %s to %s)", up$db, num(up$intercept_update, 2),
  num(up$intercept_update_lo, 2), num(up$intercept_update_hi, 2)), collapse = "; ")))

# ---- eFigure 4: stability ----
set.seed(SEED); B <- ncol(st$boot_p); n <- length(st$p_orig)
idx <- sample.int(n * B, min(n * B, 60000))
sc <- data.frame(x = rep(st$p_orig, B)[idx], y = as.vector(st$boot_p)[idx])
qb <- cut(st$p_orig, unique(stats::quantile(st$p_orig, seq(0, 1, length.out = 41))), include.lowest = TRUE)
qq <- do.call(rbind, lapply(split(seq_len(n), qb), function(i) data.frame(x = mean(st$p_orig[i]), lo = stats::quantile(st$boot_p[i, ], .025), hi = stats::quantile(st$boot_p[i, ], .975))))
ma <- min(1, ceiling(stats::quantile(c(sc$x, sc$y), 0.998) * 10) / 10)
pa <- ggplot(sc, aes(x, y)) + geom_point(size = 0.05, alpha = 0.08, colour = "grey20", shape = 16) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "#D55E00") +
  geom_line(data = qq, aes(x, lo), linetype = "dashed", linewidth = 0.35, colour = "#0072B2") + geom_line(data = qq, aes(x, hi), linetype = "dashed", linewidth = 0.35, colour = "#0072B2") +
  scale_x_continuous(labels = pct_lab) + scale_y_continuous(labels = pct_lab) + coord_fixed(xlim = c(0, ma), ylim = c(0, ma)) +
  labs(x = "Predicted risk, developed model", y = "Predicted risk, bootstrap models") + theme_journal()
nb_cal <- min(B, 100); grid <- seq(0.01, 0.8, 0.01)
cl <- do.call(rbind, lapply(seq_len(nb_cal), function(b) data.frame(b = b, calib_curve(st$y, st$boot_p[, b], grid)))); cl <- cl[!is.na(cl$obs), ]
c0 <- calib_curve(st$y, st$p_orig, grid); c0 <- c0[!is.na(c0$obs), ]; mb <- min(0.8, ceiling(max(c(cl$pred, cl$obs, c0$obs), na.rm = TRUE) * 10) / 10)
pb <- ggplot() + geom_line(data = cl, aes(pred, obs, group = b), linewidth = 0.2, colour = "grey60", alpha = 0.5) + geom_line(data = c0, aes(pred, obs), linewidth = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey30") + scale_x_continuous(labels = pct_lab) + scale_y_continuous(labels = pct_lab) +
  coord_fixed(xlim = c(0, mb), ylim = c(0, mb)) + labs(x = "Predicted risk", y = "Observed proportion (loess)") + theme_journal()
save_fig(pa + pb + plot_annotation(tag_levels = "A"), "eFigure4_stability", W2, 92)
legends["eFigure4_stability"] <- sprintf(paste("Stability of individual predictions (Riley & Collins, 2023), %d bootstrap models. (A) Prediction instability",
  "plot with 2.5th and 97.5th percentiles (dashed). (B) Calibration instability: grey lines are bootstrap models evaluated in the original data."), B)

# ---- eFigure 5: sensitivity analyses ----
tb <- r06$tab; tb <- tb[!is.na(tb$C_corrected), ]
# A10 refitted in MIMIC-IV only and its reference shown as a separate group (not directly comparable with the main analysis in all patients)
mo <- tb$id %in% c("A10", "A10ref")
mk5 <- function(t) data.frame(label = paste0(ifelse(t$id %in% c("—", "A10ref"), "", paste0(t$id, ": ")), t$analysis,
                                             "  (slope ", num(t$slope_corrected, 2), ")"),
                              sub = dn_txt(t$events, t$n), est = t$C_corrected, lo = NA, hi = NA, type = "row")
fd <- mk5(tb[!mo, ])
if (any(mo)) { tm <- tb[mo, ]; tm <- tm[order(tm$id != "A10ref"), ]
  fd <- rbind(fd, data.frame(label = "Refitted in MIMIC-IV only", sub = "", est = NA, lo = NA, hi = NA, type = "header"), mk5(tm)) }
p <- forest_plot(fd, "Optimism-corrected C statistic", threshold = tb$C_corrected[1], digits = 3,
                 head = c("Analysis (bootstrap-corrected calibration slope)", "Deaths/N", "C statistic"), widths = c(6.4, 1.9, 1.2), sub_x = 0.9,
                 xlim = round(tb$C_corrected[1], 2) + c(-1, 1) * max(0.03, ceiling(max(abs(tb$C_corrected - tb$C_corrected[1])) * 100 + 0.5) / 100))
save_fig(p, "eFigure5_sensitivity", W2, forest_height(nrow(fd)))
legends["eFigure5_sensitivity"] <- paste("Sensitivity analyses: optimism-corrected C statistic and calibration slope of the core model under each",
  "modification (bootstrap within each imputed dataset). The dashed line marks the primary analysis. Analyses adding predictors are also",
  "summarised by the apparent ΔC in the supplementary table. Analyses under 'Refitted in MIMIC-IV only' use the 1,276 MIMIC-IV patients;",
  "90-day mortality (A10) should be compared with in-hospital death refitted in the same patients, not with the primary analysis.")

# ---- eFigure 6: subgroups ----
sb <- r06$sub; cfin <- po$est[po$measure == "logit_C" & po$model == "full"]
sd <- do.call(rbind, lapply(unique(sb$group), function(g) { s <- sb[sb$group == g, ]
  rbind(hdr(g), data.frame(label = paste0("   ", s$level), sub = dn_txt(s$events, s$n), est = s$C, lo = s$C_lo, hi = s$C_hi, type = "row")) }))
p <- forest_plot(sd, "C statistic (held-out predictions)", threshold = cfin, head = c("Subgroup", "Deaths/N", "C statistic (95% CI)"))
save_fig(p, "eFigure6_subgroups", W2, forest_height(nrow(sd)))
legends["eFigure6_subgroups"] <- paste("C statistic of the full model by subgroup (held-out predictions; Rubin's rules across imputations). The dashed line marks the pooled C.",
  "Pneumonia within 24 h was taken from emergency department diagnoses in MIMIC-IV (missing, and imputed, for patients without a linked ED record)",
  "and from ICU diagnoses within 24 h or the APACHE admission diagnosis in eICU-CRD; the second pneumonia grouping is restricted to patients with recorded status.")

# ---- eFigure 7: ROC ----
COLS <- c(full = "#0072B2", BAP65 = "#7F7F7F", DECAF_mod = "#CC79A7", APACHE_IVa = "#009E73", SAPSII = "#E69F00", SOFA = "#D55E00")
NAMES <- c(full = "New model", BAP65 = "BAP-65", DECAF_mod = "Modified DECAF", APACHE_IVa = "APACHE IVa", SAPSII = "SAPS II", SOFA = "SOFA")
fpr <- seq(0, 1, by = 0.005)
roc_avg <- function(rows, pr) {
  tprs <- sapply(stk, function(s) { k <- rows(s) & !is.na(s[[pr]]); y <- s$y[k]; p <- s[[pr]][k]; o <- order(-p); y <- y[o]
    stats::approx(c(0, cumsum(1 - y) / sum(1 - y)), c(0, cumsum(y) / sum(y)), xout = fpr, ties = max, rule = 2)$y })
  cc <- mean(sapply(stk, function(s) { k <- rows(s) & !is.na(s[[pr]]); delong(s$y[k], s[[pr]][k])$auc }))
  data.frame(fpr = fpr, tpr = rowMeans(tprs), model = pr, C = cc)
}
roc_panel <- function(rows, preds, title) {
  d <- do.call(rbind, lapply(preds, function(pr) roc_avg(rows, pr))); d$model <- factor(d$model, levels = names(COLS))
  ct <- data.frame(model = factor(preds, levels = names(COLS)), k = rev(seq_along(preds)), lab = sprintf("%s  %.2f", NAMES[preds], tapply(d$C, d$model, mean)[preds]))
  ggplot(d, aes(fpr, tpr, colour = model, linetype = model)) + geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.3) + geom_line(linewidth = 0.55) +
    annotate("text", x = 0.98, y = 0.03 + length(preds) * 0.075, label = "C statistic", hjust = 1, vjust = 0, fontface = "bold", size = (BASE - 1) / .pt, family = FIG_FONT) +
    geom_text(data = ct, aes(x = 0.98, y = 0.03 + (k - 1) * 0.075, label = lab, colour = model), inherit.aes = FALSE, hjust = 1, vjust = 0, size = (BASE - 1) / .pt, family = FIG_FONT, show.legend = FALSE) +
    scale_colour_manual(values = COLS, labels = NAMES, limits = names(COLS), drop = FALSE) +
    scale_linetype_manual(values = c(full = "solid", BAP65 = "longdash", DECAF_mod = "longdash", APACHE_IVa = "longdash", SAPSII = "longdash", SOFA = "dotdash"),
                          labels = NAMES, limits = names(COLS), drop = FALSE) +
    scale_x_continuous(labels = lab_minus) + scale_y_continuous(labels = lab_minus) + coord_fixed() +
    labs(x = "1 − specificity", y = "Sensitivity", subtitle = title) + theme_journal() + guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
    theme(legend.key.width = unit(5, "mm"), plot.subtitle = element_text(size = BASE, face = "bold"))
}
p <- roc_panel(function(s) rep(TRUE, nrow(s)), c("full", "BAP65", "DECAF_mod"), "All patients") +
  roc_panel(function(s) !is.na(s$APACHE_IVa), c("full", "APACHE_IVa"), "eICU-CRD") +
  roc_panel(function(s) !is.na(s$SAPSII), c("full", "SAPSII", "SOFA"), "MIMIC-IV") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "none")
save_fig(p, "eFigure7_ROC", W2, 66)
legends["eFigure7_ROC"] <- "Receiver operating characteristic curves of held-out predictions (vertical averaging across imputations) in all patients and in the subsets with each ICU score."

# ---- eFigure 8: marker dC per cluster ----
bk <- r05$iecv_by_cluster; ip <- r05$iecv_pool
mrows <- function(cmp) { b <- bk[bk$comparison == cmp, ]; b <- b[order(as.integer(b$cluster)), ]; p <- ip[ip$comparison == cmp, ]
  rbind(data.frame(label = paste0(b$cluster, "  ", CL[b$cluster]), sub = as.character(b$deaths), est = b$est, lo = b$lo, hi = b$hi, type = "row", pi_lo = NA, pi_hi = NA),
        data.frame(label = "Pooled (random effects)", sub = "", est = p$est, lo = p$lo, hi = p$hi, type = "pooled", pi_lo = p$pi_lo, pi_hi = p$pi_hi)) }
fa <- forest_plot(mrows("E1 vs base"), "ΔC: + neutrophils, lymphocytes, albumin", ref = 0, threshold = DELTA_C_THRESHOLD, digits = 3, head = c("Held-out cluster", "Deaths", "ΔC (95% CI)"), tag = "A")
fb <- forest_plot(mrows("E3 vs base"), "ΔC: + splines of neutrophils, lymphocytes, albumin", ref = 0, threshold = DELTA_C_THRESHOLD, digits = 3, head = c("Held-out cluster", "Deaths", "ΔC (95% CI)"), tag = "B")
save_fig(wrap_plots(list(wrap_elements(full = fa), wrap_elements(full = fb)), ncol = 1), "eFigure8_markers_IECV", W2, 2 * forest_height(10))
legends["eFigure8_markers_IECV"] <- "Incremental discrimination of the markers over the core model's linear predictor in each held-out cluster (marker subset). Dashed line: prespecified threshold ΔC = 0.02."

# ---- eFigure 9: missingness ----
vv <- c(age = "Age", bun = "BUN", gcs = "GCS motor", hr = "Heart rate", rr = "Respiratory rate", map = "MAP", temp = "Temperature", ph = "pH", paco2 = "PaCO2",
        hgb = "Haemoglobin", plt = "Platelets", wbc = "White cells", glucose = "Glucose", af = "Atrial fibrillation", imv = "Invasive ventilation")
ms <- do.call(rbind, lapply(split(a, a$db), function(z) data.frame(db = z$db[1], var = vv, pct = 100 * vapply(names(vv), function(k) mean(is.na(z[[k]])), numeric(1)))))
ms$var <- factor(ms$var, levels = rev(vv))
p <- ggplot(ms, aes(pct, var, fill = db)) + geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.75), hjust = -0.15, size = (BASE - 1.5) / .pt, family = FIG_FONT) +
  scale_fill_manual(values = c(MIMIC = "#0072B2", eICU = "#E69F00"), labels = c(MIMIC = "MIMIC-IV", eICU = "eICU-CRD")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) + labs(x = "Missing, %", y = NULL) + theme_journal() + theme(legend.position = "bottom")
save_fig(p, "eFigure9_missing", W15, 100)
legends["eFigure9_missing"] <- "Proportion of missing values of the core predictors in all eligible patients, by database (handled by multiple imputation)."

writeLines(paste0(names(legends), ": ", legends, "\n"), file.path(DIR_FIG, "eFigure_legends.txt"), useBytes = TRUE)
msg("08 done")
