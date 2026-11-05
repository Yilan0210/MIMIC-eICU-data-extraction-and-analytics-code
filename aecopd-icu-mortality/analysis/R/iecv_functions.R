# =====================================================================
# iecv_functions.R: one external validation fold: development and validation parts imputed separately; spline knots from the development part;
# complete modelling in the development part (ML + uniform bootstrap shrinkage + simplified model); core model, simplified model and comparator scores evaluated in the validation part.
# Shared by 03_iecv.R (7-cluster IECV) and 04_crossdb_update.R (cross-database validation in both directions)
# =====================================================================

# comparator scores converted to risks by logistic regression in the development data, so the validation data are not used
score_model <- function(xs_d, var) {
  cf <- lapply(xs_d, function(x) {
    x <- x[!is.na(x[[var]]), ]
    if (nrow(x) < 30 || sum(x$death) < 5) return(NULL)
    stats::coef(stats::glm(stats::as.formula(paste("death ~", var)), data = x, family = stats::binomial()))
  })
  cf <- Filter(Negate(is.null), cf)
  if (!length(cf)) return(NULL)
  rowMeans(do.call(cbind, cf))
}
score_risk <- function(cf, x, var) if (is.null(cf)) rep(NA_real_, nrow(x)) else plogis(cf[1] + cf[2] * x[[var]])
# columns needed for Q3 (saved only for patients with markers, to limit file size)
KEEP_Q3 <- c("row_id", "cluster", "db", "death", "lp", "neut", "lymph", "alb", "eos", "mono", "plt", "glucose", "lactate",
             "ln_n", "ln_l", "eos_low", "eos_zero", "ln_plt", "ln_mono", "ln_glu_mmol", "ln_lac", "ln_albd", "ln_nlr", "ln_plr", "ln_sii",
             "ln_siri", "pni", "ln_glr", "ln_lar", "apache", "apache_logit", "sapsii", "sofa", "steroid_obs", "sex", "age")
SCORES <- c("BAP65", "DECAF_mod", "APACHE_IVa", "SAPSII", "SOFA")

scaled_brier <- function(y, p) { b <- mean((y - p)^2); pb <- mean(y); 1 - b / (pb * (1 - pb)) }

# Case mix (Debray framework): SD of the linear predictor in the held-out part; C of a membership model (predictors only, development vs held-out patients)
case_mix <- function(xd, xh, lp_h) {
  z <- rbind(data.frame(xd[, CORE_TERMS], g = 0), data.frame(xh[, CORE_TERMS], g = 1))
  f <- stats::glm(stats::reformulate(CORE_TERMS, "g"), data = z, family = stats::binomial())
  c(sd_lp = stats::sd(lp_h), membership_c = c_stats(z$g, stats::fitted(f))$c[1])
}

fold_fit <- function(dev, hold, label, B = B_BOOT, method = "ml_shrink", hold_with_outcome = FALSE, seed = SEED,
                     keep_tested = TRUE) {
  msg("  [", label, "] development", nrow(dev), "patients (", sum(dev$death), "deaths) / validation", nrow(hold), "patients: imputation")
  knots <- make_knots(dev)
  imp_d <- run_mice(dev, seed = seed); xs_d <- prep_xs(imp_d, knots)
  imp_h <- run_mice(hold, with_outcome = hold_with_outcome, seed = seed + 1); xs_h <- prep_xs(imp_h, knots)
  msg("  [", label, "] modelling")
  core <- develop(xs_d, CORE_TERMS, B, method, seed)
  simp <- if (method == "ml_shrink") approx_model(xs_d, core) else NULL
  sc <- list(BAP65 = score_model(xs_d, "bap65"), DECAF_mod = score_model(xs_d, "decaf"),
             SAPSII = score_model(xs_d, "sapsii"), SOFA = score_model(xs_d, "sofa"))
  rows <- list(); stacked <- vector("list", M_IMP); icis <- list(); cms <- list()
  add <- function(j, measure, model, est, var) rows[[length(rows) + 1]] <<-
    data.frame(j = j, measure = measure, model = model, est = unname(est), var = unname(var))
  for (j in seq_len(M_IMP)) {
    x <- xs_h[[j]]; y <- x$death
    P <- cbind(full = predict_dev(core, x)); if (!is.null(simp)) P <- cbind(P, simplified = predict_dev(simp, x))
    cs <- c_stats(y, P)
    for (k in seq_len(ncol(P))) {
      nm <- colnames(P)[k]
      add(j, "logit_C", nm, cs$logit_c[k], cs$var_logit_c[k])
      cl <- calib(y, P[, k])
      add(j, "log_OE", nm, cl["log_oe"], cl["var_log_oe"]); add(j, "cal_slope", nm, cl["slope"], cl["var_slope"])
      add(j, "CITL", nm, cl["citl"], cl["var_citl"])
      add(j, "scaled_Brier", nm, scaled_brier(y, P[, k]), 0)
      for (pt in NB_THRESHOLDS) { nb <- net_benefit(y, P[, k], pt); add(j, sprintf("NB_%02d", pt * 100), nm, nb["nb"], nb["var"]) }
      icis[[length(icis) + 1]] <- data.frame(j = j, model = nm, t(ici(y, P[, k])))
    }
    for (pt in NB_THRESHOLDS) { nb <- net_benefit(y, rep(1, length(y)), pt); add(j, sprintf("NB_%02d", pt * 100), "treat_all", nb["nb"], nb["var"]) }
    if (ncol(P) > 1) add(j, "dC", "simplified-full", cs$dc[1], cs$var_dc[1])
    risks <- list(BAP65 = score_risk(sc$BAP65, x, "bap65"), DECAF_mod = score_risk(sc$DECAF_mod, x, "decaf"),
                  APACHE_IVa = x$apache, SAPSII = score_risk(sc$SAPSII, x, "sapsii"), SOFA = score_risk(sc$SOFA, x, "sofa"))
    for (s in names(risks)) {
      r <- risks[[s]]; idx <- which(!is.na(r))
      if (length(idx) < 30 || sum(y[idx]) < 2 || sum(1 - y[idx]) < 2) next
      cc <- c_stats(y[idx], cbind(r[idx], P[idx, "full"]))
      add(j, "logit_C", s, cc$logit_c[1], cc$var_logit_c[1])
      add(j, "logit_C_same_patients", paste0("full|", s), cc$logit_c[2], cc$var_logit_c[2])
      add(j, "dC_vs_score", paste0("full-", s), cc$dc[1], cc$var_dc[1])
      clr <- calib(y[idx], r[idx]); add(j, "log_OE", s, clr["log_oe"], clr["var_log_oe"]); add(j, "cal_slope", s, clr["slope"], clr["var_slope"])
      for (pt in NB_THRESHOLDS) {
        nb <- net_benefit(y[idx], r[idx], pt); add(j, sprintf("NB_%02d", pt * 100), s, nb["nb"], nb["var"])
        nf <- net_benefit(y[idx], P[idx, "full"], pt); add(j, sprintf("NB_%02d", pt * 100), paste0("full|", s), nf["nb"], nf["var"])
      }
    }
    lp <- qlogis(P[, "full"])
    if (j <= 3) cms[[j]] <- case_mix(xs_d[[j]], x, lp)
    stacked[[j]] <- data.frame(row_id = x$row_id, label = label, db = x$db, y = y, P, lp = lp,
                               BAP65 = risks$BAP65, DECAF_mod = risks$DECAF_mod, APACHE_IVa = risks$APACHE_IVa,
                               SAPSII = risks$SAPSII, SOFA = risks$SOFA, sex = as.character(x$sex), race4 = x$race4,
                               age = x$age, imv = x$imv, pneumonia24 = x$pneumonia24, tested = x$tested,
                               death_s3 = x$death_s3)
  }
  long <- do.call(rbind, rows)
  pooled <- do.call(rbind, lapply(split(long, list(long$measure, long$model), drop = TRUE), function(s) {
    r <- rubin(s$est, s$var); data.frame(measure = s$measure[1], model = s$model[1], est = r$est, se = r$se)
  }))
  ici_tab <- aggregate(cbind(ICI, E50, E90) ~ model, data = do.call(rbind, icis), FUN = mean)
  cm <- colMeans(do.call(rbind, cms))
  # Marker subset (for the Q3 IECV): patients with markers in the development and held-out parts, with the linear predictor of this fold's core model
  sub_t <- function(x) { x <- x[x$tested == 1, ]; x$lp <- qlogis(predict_dev(core, x)); x[, intersect(KEEP_Q3, names(x))] }
  tested <- if (keep_tested) list(dev = lapply(xs_d, sub_t), hold = lapply(xs_h, sub_t)) else NULL
  list(label = label, n = nrow(hold), deaths = sum(hold$death), n_dev = nrow(dev), deaths_dev = sum(dev$death),
       pooled = pooled, ici = ici_tab, case_mix = cm, stacked = stacked, knots = knots,
       core = core[c("terms", "coef", "shrink", "method")], simp = simp[c("terms", "vars", "coef", "r2")], tested = tested)
}

# Format the Rubin-pooled results of one fold as a readable table (C, O/E back-transformed; normal-approximation 95% CI)
fold_table <- function(f) {
  p <- f$pooled
  tr <- function(meas, back) {
    s <- p[p$measure == meas, ]; if (!nrow(s)) return(NULL)
    data.frame(label = f$label, n = f$n, deaths = f$deaths, measure = meas, model = s$model,
               est = back(s$est), lo = back(s$est - 1.96 * s$se), hi = back(s$est + 1.96 * s$se))
  }
  do.call(rbind, list(tr("logit_C", plogis), tr("logit_C_same_patients", plogis), tr("log_OE", exp), tr("cal_slope", identity),
        tr("CITL", identity), tr("scaled_Brier", identity), tr("dC", identity), tr("dC_vs_score", identity),
        do.call(rbind, lapply(sprintf("NB_%02d", NB_THRESHOLDS * 100), tr, back = identity))))
}
