# =====================================================================
# 06_sensitivity.R: sensitivity analyses A1–A14, subgroups and fairness
#   Each: same modelling approach; reports apparent C, bootstrap-corrected C and calibration slope (= shrinkage factor), and apparent dC vs the core model in the same patients
#   Computation: first AECOPD_SENS_M imputations (default 10), AECOPD_SENS_B bootstrap samples per imputation (default 100)
#   A12 = rerun 03_iecv.R with AECOPD_IECV_HOLDOUT_OUTCOME=1; A10 needs the 90-day death file (see README)
# =====================================================================
source("00_setup.R")
t0 <- Sys.time()
SENS_M <- as.integer(Sys.getenv("AECOPD_SENS_M", "10")); SENS_B <- as.integer(Sys.getenv("AECOPD_SENS_B", "100"))
d02 <- readRDS(file.path(DIR_LOCAL, "02_dev.rds")); knots <- d02$knots; core <- d02$core
xs <- d02$xs[seq_len(min(SENS_M, length(d02$xs)))]
a <- readRDS(file.path(DIR_LOCAL, "01_population.rds"))$a
raw <- load_clean(); num <- function(z) suppressWarnings(as.numeric(z))

variant <- function(id, name, xs_v, terms = CORE_TERMS, y = "death", note = "") {
  xs_v <- lapply(xs_v, function(x) x[!is.na(x[[y]]), ])
  per <- t(sapply(seq_along(xs_v), function(i) { x <- xs_v[[i]]
    v <- validate_one(terms, x, B = SENS_B, seed = SEED + 900 + i, y = y)
    f0 <- fit_glm(CORE_TERMS, x, y); f1 <- fit_glm(terms, x, y)
    cs <- c_stats(x[[y]], cbind(stats::fitted(f0), stats::fitted(f1)))
    c(v[c("c_app", "c_corr", "shrink")], dc = cs$dc, var = cs$var_dc) }))
  rd <- if (identical(terms, CORE_TERMS)) data.frame(est = NA, lo = NA, hi = NA) else rubin(per[, "dc"], per[, "var"])
  data.frame(id = id, analysis = name, n = nrow(xs_v[[1]]), events = sum(xs_v[[1]][[y]]), parameters = length(terms),
             C_apparent = mean(per[, "c_app"]), C_corrected = mean(per[, "c_corr"]), slope_corrected = mean(per[, "shrink"]),
             dC_vs_core_apparent = rd$est, dC_lo = rd$lo, dC_hi = rd$hi, note = note)
}
res <- list(variant("—", "Primary analysis (core model)", xs))

msg("A1 complete cases")
need <- c("death", "age", "sex", "bun", "gcs", "hr", "rr", "map", "temp", "ph", "paco2", "hgb", "plt", "wbc", "glucose", "af", "imv")
acc <- a[stats::complete.cases(a[, need]), ]
res[[length(res) + 1]] <- variant("A1", "Complete-case analysis", list(add_design(derive(acc), knots)),
                                  note = sprintf("%d of %d patients with complete core predictors", nrow(acc), nrow(a)))
msg("A2-A4, A7: added variables")
res[[length(res) + 1]] <- variant("A2", "Adding lactate (log)", xs, c(CORE_TERMS, "ln_lac"))
res[[length(res) + 1]] <- variant("A3", "Adding day-1 vasopressors", xs, c(CORE_TERMS, "vaso"))
res[[length(res) + 1]] <- variant("A4", "Adding pneumonia within 24 h", xs, c(CORE_TERMS, "pneumonia24"))
res[[length(res) + 1]] <- variant("A7", "Adding a database indicator", xs, c(CORE_TERMS, "db_mimic"))
msg("A5 treatment limitation")
res[[length(res) + 1]] <- variant("A5", "Excluding treatment limitation within 24 h",
  lapply(xs, function(x) x[!(x$limit_s18 %in% 1), ]), note = "patients without code-status records kept")
msg("A6 missing blood gas indicator")
miss_abg <- is.na(a$ph) | is.na(a$paco2)
xs6 <- lapply(xs, function(x) { x$abg_missing <- as.numeric(miss_abg); x$ph[miss_abg] <- 7.40; x$paco2[miss_abg] <- 40; add_design(x, knots) })
res[[length(res) + 1]] <- variant("A6", "No-blood-gas indicator (pH 7.40, PaCO2 40 imputed)", xs6, c(CORE_TERMS, "abg_missing"),
                                  note = sprintf("%.1f%% without blood gas", 100 * mean(miss_abg)))
msg("A9 outcome including hospice")
res[[length(res) + 1]] <- variant("A9", "Outcome: in-hospital death or discharge to hospice", xs, y = "death_s3",
                                  note = "eICU has no hospice category (equals in-hospital death)")
msg("A11 age truncation")
gt89 <- if ("age_gt89" %in% names(raw)) num(raw$age_gt89) else rep(NA_real_, nrow(raw))
res[[length(res) + 1]] <- variant("A11", "Excluding age-capped patients (>89 years)",
  lapply(xs, function(x) x[if (all(is.na(gt89))) x$age < 90 else !(gt89[x$row_id] %in% 1), ]))
msg("A13 fitted in the marker subset only")
res[[length(res) + 1]] <- variant("A13", "Core model fitted in the marker-tested subset only", lapply(xs, function(x) x[x$tested == 1, ]),
                                  note = "24 df with ~230 deaths: below Riley sample size; for comparison with the marker analyses only")
msg("A10 90-day death")
f90 <- file.path(DATA_DIR, sprintf("aecopd_mimic_death90_%s_LOCAL_ONLY.csv", TAG))
if (file.exists(f90)) {
  d90 <- utils::read.csv(f90, stringsAsFactors = FALSE)
  xs10 <- lapply(xs, function(x) { x <- x[x$db == "MIMIC", ]; x$death_90d <- d90$death_90d[match(as.character(raw$stay_id[x$row_id]), as.character(d90$stay_id))]; x })
  # Reference for A10: in-hospital death refitted in the same MIMIC patients (A10 cannot be compared directly with the main analysis in all patients)
  res[[length(res) + 1]] <- variant("A10ref", "In-hospital death, same MIMIC-IV patients (reference for A10)", xs10, y = "death")
  res[[length(res) + 1]] <- variant("A10", "90-day all-cause mortality (MIMIC-IV only)", xs10, y = "death_90d")
} else res[[length(res) + 1]] <- data.frame(id = "A10", analysis = "90-day mortality", n = NA, events = NA, parameters = NA, C_apparent = NA,
  C_corrected = NA, slope_corrected = NA, dC_vs_core_apparent = NA, dC_lo = NA, dC_hi = NA, note = paste("not run: missing", basename(f90)))
tab <- do.call(rbind, res)

# ---- A8: ridge and LASSO ----
msg("A8 ridge / LASSO")
a8 <- do.call(rbind, lapply(c(ridge = 0, lasso = 1), function(al) {
  per <- t(sapply(seq_along(xs), function(i) { x <- xs[[i]]; X <- as.matrix(x[, CORE_TERMS]); set.seed(SEED + i)
    cv <- glmnet::cv.glmnet(X, x$death, family = "binomial", alpha = al, nfolds = 10)
    p <- drop(stats::predict(cv, X, s = "lambda.min", type = "response")); b <- as.numeric(stats::coef(cv, s = "lambda.min"))[-1]
    c(C = c_stats(x$death, p)$c, cor_lp = stats::cor(stats::qlogis(p), qlogis(predict_dev(core, x))), nonzero = sum(b != 0)) }))
  data.frame(method = if (al == 0) "Ridge" else "LASSO", C_apparent = mean(per[, "C"]), correlation_with_core_LP = mean(per[, "cor_lp"]),
             nonzero_terms = stats::median(per[, "nonzero"]))
}))
a8 <- rbind(data.frame(method = "Core (ML + uniform shrinkage)", C_apparent = mean(core$val[, "c_app"]), correlation_with_core_LP = 1,
                       nonzero_terms = length(CORE_TERMS)), a8)

# ---- A14: gradient boosting (cluster by cluster held out; first 5 imputations; same predictors as the logistic model) ----
a14 <- NULL
if (requireNamespace("gbm", quietly = TRUE)) {
  msg("A14 gradient boosting")
  ks <- levels(a$cluster); M5 <- min(5L, length(d02$xs))
  one <- function(j, k) { x <- d02$xs[[j]]; dev <- x[x$cluster != k, ]; ho <- x[x$cluster == k, ]
    lg <- fit_glm(CORE_TERMS, dev); set.seed(SEED + 10 * j + as.integer(k)); gb <- NULL
    utils::capture.output(gb <- gbm::gbm(stats::reformulate(CORE_VARS, "death"), data = dev[, c("death", CORE_VARS)], distribution = "bernoulli",
      n.trees = 1500, interaction.depth = 2, shrinkage = 0.01, n.minobsinnode = 20, bag.fraction = 0.7, cv.folds = 5, n.cores = 1, verbose = FALSE))
    nt <- suppressMessages(gbm::gbm.perf(gb, method = "cv", plot.it = FALSE))
    P <- cbind(stats::predict(lg, newdata = ho, type = "response"), stats::predict(gb, newdata = ho, n.trees = nt, type = "response"))
    cs <- c_stats(ho$death, P); cl <- calib(ho$death, P[, 2])
    data.frame(j = j, cluster = k, lr = cs$logit_c[1], v_lr = cs$var_logit_c[1], gbm = cs$logit_c[2], v_gbm = cs$var_logit_c[2],
               dc = cs$dc, v_dc = cs$var_dc, oe = cl[["log_oe"]], v_oe = cl[["var_log_oe"]], sl = cl[["slope"]], v_sl = cl[["var_slope"]]) }
  grid <- expand.grid(j = seq_len(M5), k = ks, stringsAsFactors = FALSE)
  r <- do.call(rbind, if (N_CORES > 1) parallel::mclapply(seq_len(nrow(grid)), function(i) one(grid$j[i], grid$k[i]), mc.cores = N_CORES)
               else lapply(seq_len(nrow(grid)), function(i) one(grid$j[i], grid$k[i])))
  byk <- do.call(rbind, lapply(split(r, r$cluster), function(s) { g <- function(e, v) rubin(s[[e]], s[[v]])[, c("est", "se")]
    data.frame(cluster = s$cluster[1], lr = g("lr", "v_lr"), gbm = g("gbm", "v_gbm"), dc = g("dc", "v_dc"), oe = g("oe", "v_oe"), sl = g("sl", "v_sl")) }))
  a14 <- rbind(re_pool(byk$lr.est, byk$lr.se, plogis, "C, logistic regression (ML, same predictors)"),
               re_pool(byk$gbm.est, byk$gbm.se, plogis, "C, gradient boosting"),
               re_pool(byk$dc.est, byk$dc.se, identity, "ΔC, gradient boosting − logistic"),
               re_pool(byk$oe.est, byk$oe.se, exp, "O/E, gradient boosting"), re_pool(byk$sl.est, byk$sl.se, identity, "Calibration slope, gradient boosting"))
  a14$note <- sprintf("leave-one-cluster-out on the first %d imputed datasets; depth 2, learning rate 0.01, trees by 5-fold CV", M5)
}

# ---- Subgroups and fairness (IECV held-out predictions) ----
msg("Subgroups")
stk <- readRDS(file.path(DIR_LOCAL, "03_iecv_stacked.rds"))
grp_def <- list(Sex = function(s) ifelse(s$sex == "F", "Female", "Male"), Race = function(s) s$race4,
                Age = function(s) ifelse(s$age >= 80, "≥80 years", "<80 years"), `Invasive ventilation` = function(s) ifelse(s$imv == 1, "Yes", "No"),
                `Pneumonia within 24 h` = function(s) ifelse(s$pneumonia24 == 1, "Yes", "No"),
                # In MIMIC, pneumonia status is imputed for patients without an emergency department record; a second grouping uses patients with a recorded status only
                `Pneumonia within 24 h (recorded only)` = function(s) { z <- a$pneumonia24[match(s$row_id, a$row_id)]
                  ifelse(z %in% 1, "Yes", ifelse(z %in% 0, "No", NA_character_)) },
                Database = function(s) s$db)
sub <- do.call(rbind, lapply(names(grp_def), function(g) {
  per <- do.call(rbind, lapply(seq_along(stk), function(j) { s <- stk[[j]]; s$lev <- grp_def[[g]](s)
    do.call(rbind, lapply(split(s, s$lev), function(z) { if (sum(z$y) < 2 || sum(1 - z$y) < 2) return(NULL); cs <- c_stats(z$y, z$full)
      data.frame(j = j, level = z$lev[1], n = nrow(z), events = sum(z$y), lc = cs$logit_c, v = cs$var_logit_c, oe = sum(z$y) / sum(z$full)) })) }))
  do.call(rbind, lapply(split(per, per$level), function(z) { r <- rubin(z$lc, z$v)
    data.frame(group = g, level = z$level[1], n = round(mean(z$n)), events = round(mean(z$events)), C = plogis(r$est), C_lo = plogis(r$lo), C_hi = plogis(r$hi),
               O_E = mean(z$oe), note = ifelse(mean(z$events) < 20, "events < 20: unstable", "")) })) }))

saveRDS(list(tab = tab, a8 = a8, a14 = a14, sub = sub), file.path(DIR_LOCAL, "06_sensitivity.rds"))
write_share(list(README = data.frame(note = c(sprintf("Each analysis uses the first %d imputations and %d bootstrap samples per imputation; no multiplicity correction; the main conclusions are not revised on this basis", SENS_M, SENS_B),
    "dC_vs_core_apparent: apparent dC of the variant vs the core model in the same patients (DeLong + Rubin); for added-variable analyses this is the main column",
    "A12: see 03_iecv_A12.xlsx; A8 compares only the apparent C and the agreement of linear predictors", "subgroups: C and O/E of IECV held-out predictions by subgroup", "This file contains summary statistics only")),
  sensitivity = tab, A8_penalised = a8, A14_gradient_boosting = if (is.null(a14)) data.frame(note = "gbm not installed") else a14, subgroups = sub),
  "06_sensitivity.xlsx")
msg("06 done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min")
