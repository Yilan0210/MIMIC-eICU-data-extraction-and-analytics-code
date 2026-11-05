# =====================================================================
# 02_develop.R: development of the core model
#   * multiple imputation in all eligible patients (with the outcome); spline knots at the 10th/50th/90th percentiles of observed values
#   * full model: 16 variables, 24 df; ML + uniform bootstrap shrinkage, intercept re-estimated; bootstrap optimism correction
#   * simplified model: Harrell's model approximation (R2 >= 0.95), used for the nomogram
#   * effects: interquartile-range ORs, partial effect curves, variable importance; prediction stability (MAPE)
# Outputs: result/02_development.xlsx; local/02_dev.rds (imputed data and models)
# =====================================================================
source("00_setup.R")
t0 <- Sys.time()
a <- readRDS(file.path(DIR_LOCAL, "01_population.rds"))$a
knots <- make_knots(a)
msg("Imputation: m =", M_IMP, ", maxit =", MAXIT, ", ", nrow(a), "patients")
imp <- run_mice(a)
xs <- prep_xs(imp, knots)
if (sum(vapply(xs[[1]][, c("death", CORE_TERMS)], function(z) sum(is.na(z)), numeric(1))) > 0) stop("Core model variables still missing after imputation")

msg("Full model:", length(CORE_TERMS), "df, B =", B_BOOT)
core <- develop(xs, CORE_TERMS)
msg("  shrinkage factor", round(core$shrink, 3))
msg("Simplified model (model approximation)")
simp <- approx_model(xs, core)
msg("  kept", length(simp$vars), "variables:", paste(simp$vars, collapse = ", "), "; R2 =", round(simp$r2, 3))

# ---- Performance (apparent and optimism-corrected) ----
val <- core$val
c_simp_app <- mean(vapply(xs, function(x) c_stats(x$death, predict_dev(simp, x))$c, numeric(1)))
perf <- data.frame(model = c("Full model", "Simplified model"), parameters = c(length(CORE_TERMS), length(simp$terms)),
  C_apparent = c(mean(val[, "c_app"]), c_simp_app), C_optimism_corrected = c(mean(val[, "c_corr"]), NA),
  calibration_slope_corrected = c(core$shrink, NA),
  scaled_Brier_apparent = c(mean(1 - val[, "brier_app"] / (mean(a$death) * (1 - mean(a$death)))), NA),
  scaled_Brier_corrected = c(mean(1 - val[, "brier_corr"] / (mean(a$death) * (1 - mean(a$death)))), NA),
  R2_vs_full_LP = c(1, simp$r2))

# ---- Equation ----
eq <- data.frame(term = names(core$coef), coef_ML = unname(core$coef_ml[names(core$coef)]),
                 SE_ML = sqrt(diag(core$vcov_ml))[names(core$coef)], coef_shrunken = unname(core$coef))
kn_tab <- data.frame(variable = names(knots), knot1 = sapply(knots, `[`, 1), knot2 = sapply(knots, `[`, 2), knot3 = sapply(knots, `[`, 3),
                     scale = ifelse(names(knots) == "ln_bun", "ln(BUN, mg/dL)", "original"))
eq_simp <- data.frame(term = names(simp$coef), coef = unname(simp$coef))
simp_freq <- data.frame(variable = CORE_VARS, retained_in_share_of_imputations = simp$freq, retained = CORE_VARS %in% simp$vars)

# ---- Effects ----
x1 <- xs[[1]]
orc <- do.call(rbind, lapply(CORE_VARS, function(v) or_contrast(core, v, x1, knots)))
orc$from_raw <- mapply(to_raw, orc$var, orc$from); orc$to_raw <- mapply(to_raw, orc$var, orc$to)
imp_tab <- importance(core); imp_tab <- imp_tab[order(-imp_tab$chi2_minus_df), ]
curves <- do.call(rbind, lapply(setdiff(CORE_VARS, c("female", "af", "imv")), function(v) {
  z <- x1[[v]]; g <- seq(stats::quantile(z, .01), stats::quantile(z, .99), length.out = 100)
  pe <- partial_effect(core, v, g, stats::median(z), knots); pe$var <- v; pe$x_raw <- to_raw(v, pe$x); pe }))
nl_test <- do.call(rbind, lapply(SPLINE_VARS, function(v) { k <- paste0(v, "_nl"); b <- core$coef_ml[k]; se <- sqrt(core$vcov_ml[k, k])
  data.frame(var = v, chi2_nonlinear = (b / se)^2, p_nonlinear = 2 * stats::pnorm(-abs(b / se))) }))

# ---- Prediction stability (Riley & Collins 2023): imputed dataset 1, B bootstrap refits (ML, same shrinkage factor, intercept re-estimated) ----
msg("Stability: B =", B_BOOT)
set.seed(SEED); p_orig <- predict_dev(core, x1)
boot_p <- sapply(seq_len(B_BOOT), function(b) {
  d <- x1[sample.int(nrow(x1), replace = TRUE), ]; f <- fit_glm(CORE_TERMS, d)
  bb <- stats::coef(f)[-1] * core$shrink; off <- drop(as.matrix(d[, CORE_TERMS]) %*% bb)
  int <- stats::coef(stats::glm(d$death ~ 1 + offset(off), family = stats::binomial()))[1]
  plogis(int + drop(as.matrix(x1[, CORE_TERMS]) %*% bb)) })
ape <- abs(boot_p - p_orig)
stab <- data.frame(B = B_BOOT, MAPE = mean(ape), median_individual_MAPE = stats::median(rowMeans(ape)),
                   p90_individual_MAPE = unname(stats::quantile(rowMeans(ape), .9)), max_individual_MAPE = max(rowMeans(ape)))

saveRDS(list(xs = xs, knots = knots, core = core, simp = simp, curves = curves, orc = orc, imp_tab = imp_tab),
        file.path(DIR_LOCAL, "02_dev.rds"))
saveRDS(list(y = x1$death, p_orig = p_orig, boot_p = boot_p), file.path(DIR_LOCAL, "02_stability.rds"))
write_share(list(
  README = data.frame(note = c("Core model (full model): ML coefficients pooled by Rubin's rules; shrunken coefficients (intercept re-estimated) are used for prediction; spline non-linear terms computed by rcs_nl() (identical to rms::rcs)",
    "equation + knots reproduce the full model; the web calculator in calculator/ uses the same coefficients",
    "OR_IQR: OR for the 75th vs 25th percentile of continuous variables (log-transformed variables reported on the original scale), 1 vs 0 for binary variables; based on pooled ML coefficients",
    "importance: pooled Wald chi-square (both spline terms jointly) minus degrees of freedom", "simplified: variables and coefficients kept by model approximation (predicting the linear predictor of the full model)",
    "This file contains summary statistics only")),
  performance = perf, equation = eq, knots = kn_tab, OR_IQR = orc, importance = imp_tab, nonlinearity = nl_test,
  simplified_equation = eq_simp, simplified_selection = simp_freq, stability = stab,
  effect_curves = curves[, c("var", "x_raw", "logor", "lo", "hi")]), "02_development.xlsx")
msg("02 done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min")
