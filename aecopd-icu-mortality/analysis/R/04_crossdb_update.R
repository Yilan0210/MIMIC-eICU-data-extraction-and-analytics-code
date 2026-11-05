# =====================================================================
# 04_crossdb_update.R: cross-database validation in both directions and model updating
#   eICU -> MIMIC: ML + bootstrap shrinkage; MIMIC -> eICU: MIMIC is too small for 24 df, so ridge regression is used
#   Model updating: with IECV held-out predictions, intercept-only and intercept-and-slope recalibration within each database;
#            updating parameters (bootstrap 95% CI) and calibration before and after updating (after updating = apparent) are reported
# =====================================================================
source("00_setup.R"); source("iecv_functions.R")
t0 <- Sys.time()
a <- readRDS(file.path(DIR_LOCAL, "01_population.rds"))$a
stk <- readRDS(file.path(DIR_LOCAL, "03_iecv_stacked.rds"))

msg("Cross-database: eICU -> MIMIC")
f1 <- fold_fit(a[a$db == "eICU", ], a[a$db == "MIMIC", ], "eICU → MIMIC", seed = SEED + 701, keep_tested = FALSE)
msg("Cross-database: MIMIC -> eICU (ridge)")
f2 <- fold_fit(a[a$db == "MIMIC", ], a[a$db == "eICU", ], "MIMIC → eICU (ridge)", method = "ridge", seed = SEED + 702, keep_tested = FALSE)
cross <- rbind(fold_table(f1), fold_table(f2))

# ---- Model updating (held-out predictions averaged across imputations) ----
ord <- lapply(stk, function(s) s[order(s$row_id), ]); av <- ord[[1]][, c("row_id", "db", "y")]
av$lp <- rowMeans(sapply(ord, `[[`, "lp"))
set.seed(SEED)
upd <- do.call(rbind, lapply(split(av, av$db), function(s) {
  fi <- stats::glm(y ~ 1 + offset(lp), data = s, family = stats::binomial()); fs <- stats::glm(y ~ lp, data = s, family = stats::binomial())
  bt <- replicate(500, { d <- s[sample.int(nrow(s), replace = TRUE), ]
    c(stats::coef(stats::glm(y ~ 1 + offset(lp), data = d, family = stats::binomial())), stats::coef(stats::glm(y ~ lp, data = d, family = stats::binomial()))) })
  p0 <- plogis(s$lp); p1 <- plogis(stats::coef(fi) + s$lp); p2 <- stats::fitted(fs)
  cal <- function(p) { c0 <- calib(s$y, p); c(OE = exp(c0[["log_oe"]]), CITL = c0[["citl"]], slope = c0[["slope"]], ICI = ici(s$y, p)[["ICI"]]) }
  data.frame(db = s$db[1], n = nrow(s), deaths = sum(s$y),
             intercept_update = stats::coef(fi), intercept_update_lo = stats::quantile(bt[1, ], .025), intercept_update_hi = stats::quantile(bt[1, ], .975),
             recal_intercept = stats::coef(fs)[1], recal_slope = stats::coef(fs)[2],
             recal_slope_lo = stats::quantile(bt[3, ], .025), recal_slope_hi = stats::quantile(bt[3, ], .975),
             t(setNames(cal(p0), paste0(c("OE", "CITL", "slope", "ICI"), "_original"))),
             t(setNames(cal(p1), paste0(c("OE", "CITL", "slope", "ICI"), "_intercept_updated"))),
             t(setNames(cal(p2), paste0(c("OE", "CITL", "slope", "ICI"), "_logistic_recalibrated"))))
}))
curves_upd <- do.call(rbind, lapply(split(av, av$db), function(s) {
  fi <- stats::glm(y ~ 1 + offset(lp), data = s, family = stats::binomial())
  rbind(data.frame(db = s$db[1], version = "Original", calib_curve(s$y, plogis(s$lp))),
        data.frame(db = s$db[1], version = "Intercept updated", calib_curve(s$y, plogis(stats::coef(fi) + s$lp)))) }))
saveRDS(list(cross = cross, upd = upd, curves_upd = curves_upd), file.path(DIR_LOCAL, "04_crossdb_update.rds"))
write_share(list(README = data.frame(note = c("cross_database: cross-database validation in both directions (Rubin's rules across imputations)",
    "model_updating: recalibration within each database using IECV held-out predictions; bootstrap 95% CI for the updating parameters (500 samples); calibration after updating is apparent and shows how much local recalibration can correct",
    "Local calibration before use in a new hospital: logit(p_new) = intercept update + model linear predictor (intercept-only update)", "This file contains summary statistics only")),
  cross_database = cross, model_updating = upd, calibration_curves = curves_upd), "04_crossdb_update.xlsx")
msg("04 done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min")
