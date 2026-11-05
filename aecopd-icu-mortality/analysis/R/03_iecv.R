# =====================================================================
# 03_iecv.R: internal-external cross-validation over 7 clusters (Q1, Q2)
#   Each cluster held out in turn: remaining clusters imputed (with outcome) -> knots, full model (ML + shrinkage), simplified model, comparator score calibration
#                held-out cluster imputed separately (without outcome) -> evaluation
#   Pooling: logit C, log O/E, calibration slope, CITL, scaled Brier, net benefit (10/20/30%); REML + HKSJ; 95% PI with t(k-2)
#   Comparator scores: paired comparison with the full model in the same patients (dC, calibration, net benefit)
#   Case mix: SD of the linear predictor per cluster, membership model C
# A12 (held-out clusters imputed with the outcome): set AECOPD_IECV_HOLDOUT_OUTCOME=1 before running; outputs get the suffix _A12
# =====================================================================
source("00_setup.R"); source("iecv_functions.R")
t0 <- Sys.time()
A12 <- Sys.getenv("AECOPD_IECV_HOLDOUT_OUTCOME", "0") == "1"; sfx <- if (A12) "_A12" else ""
a <- readRDS(file.path(DIR_LOCAL, "01_population.rds"))$a
cls <- levels(a$cluster)
CL_LABEL <- c("1" = "MIMIC ICD-9 era", "2" = "MIMIC ICD-10 era", "3" = "eICU Midwest", "4" = "eICU Northeast",
              "5" = "eICU South", "6" = "eICU West", "7" = "eICU region missing")
msg("IECV:", length(cls), "clusters", if (A12) " (A12: held-out clusters imputed with the outcome)" else "")
run_k <- function(k) {
  dev <- a[a$cluster != k, ]; dev$cluster <- droplevels(dev$cluster)
  fold_fit(dev, a[a$cluster == k, ], label = k, hold_with_outcome = A12, seed = SEED + 100 * as.integer(k), keep_tested = !A12)
}
folds <- if (N_CORES > 1) parallel::mclapply(cls, run_k, mc.cores = N_CORES, mc.preschedule = FALSE) else lapply(cls, run_k)
bad <- vapply(folds, function(f) inherits(f, "try-error") || is.null(f$pooled), logical(1))
if (any(bad)) stop("IECV failed for clusters: ", paste(cls[bad], collapse = ", "), ". Set AECOPD_CORES=1 and rerun to see the error")
names(folds) <- cls
saveRDS(folds, file.path(DIR_LOCAL, paste0("03_iecv_folds", sfx, ".rds")))

by_cluster <- do.call(rbind, lapply(folds, fold_table)); by_cluster$cluster_name <- CL_LABEL[by_cluster$label]
ev <- vapply(folds, function(f) f$deaths, numeric(1)); eligible <- names(ev)[ev >= MIN_EVENTS_POOL]
msg("Clusters pooled: ", paste(eligible, collapse = ", "))
pool_measure <- function(meas, model, back = identity, use = eligible, tag = "all eligible") {
  s <- do.call(rbind, lapply(folds[use], function(f) { p <- f$pooled; p[p$measure == meas & p$model == model, c("est", "se")] }))
  if (is.null(s) || nrow(s) < 2) return(NULL)
  out <- re_pool(s$est, s$se, back = back, label = meas); out$model <- model; out$clusters <- tag; out
}
pooled <- list()
for (mm in c("full", "simplified")) {
  for (z in list(c("logit_C", "plogis"), c("log_OE", "exp"), c("cal_slope", "identity"), c("CITL", "identity")))
    pooled[[length(pooled) + 1]] <- pool_measure(z[1], mm, match.fun(z[2]))
  for (pt in NB_THRESHOLDS) pooled[[length(pooled) + 1]] <- pool_measure(sprintf("NB_%02d", pt * 100), mm)
}
for (pt in NB_THRESHOLDS) pooled[[length(pooled) + 1]] <- pool_measure(sprintf("NB_%02d", pt * 100), "treat_all")
pooled[[length(pooled) + 1]] <- pool_measure("dC", "simplified-full")
for (s in SCORES) {
  has <- eligible[vapply(folds[eligible], function(f) any(f$pooled$model == s), logical(1))]
  if (length(has) < 2) next
  tag <- paste("clusters", paste(has, collapse = ","))
  pooled[[length(pooled) + 1]] <- pool_measure("logit_C", s, plogis, has, tag)
  pooled[[length(pooled) + 1]] <- pool_measure("logit_C_same_patients", paste0("full|", s), plogis, has, tag)
  pooled[[length(pooled) + 1]] <- pool_measure("dC_vs_score", paste0("full-", s), identity, has, tag)
  pooled[[length(pooled) + 1]] <- pool_measure("log_OE", s, exp, has, tag)
  pooled[[length(pooled) + 1]] <- pool_measure("cal_slope", s, identity, has, tag)
  for (pt in NB_THRESHOLDS) for (mm in c(s, paste0("full|", s)))
    pooled[[length(pooled) + 1]] <- pool_measure(sprintf("NB_%02d", pt * 100), mm, identity, has, tag)
}
pooled <- do.call(rbind, Filter(Negate(is.null), pooled))
pooled$note <- ifelse(pooled$k == 2, "Only 2 clusters: fixed-effect pooling, no prediction interval", ifelse(pooled$k < 2, "Fewer than 2 clusters can be pooled (e.g. no patient classified positive at this threshold); not pooled", ""))

# ---- Stacked held-out predictions: decision curves, calibration curves, risk groups ----
stk <- lapply(seq_len(M_IMP), function(j) do.call(rbind, lapply(folds, function(f) f$stacked[[j]])))
curves <- c("full", "simplified", SCORES)
dca <- do.call(rbind, lapply(DCA_GRID, function(pt) {
  nb <- sapply(stk, function(s) c(vapply(curves, function(cv) { ok <- !is.na(s[[cv]])
    if (sum(ok) < 30) NA_real_ else unname(net_benefit(s$y[ok], s[[cv]][ok], pt)["nb"]) }, numeric(1)), treat_all = nb_treat_all(s$y, pt)))
  data.frame(threshold = pt, t(rowMeans(nb)))
}))
cal_curve <- do.call(rbind, lapply(c("full", "simplified"), function(mm) {
  cc <- lapply(stk, function(s) calib_curve(s$y, s[[mm]])); g <- Reduce(function(x, y) merge(x, y, by = "pred"), cc)
  data.frame(model = mm, pred = g$pred, obs = rowMeans(g[, -1, drop = FALSE])) }))
ici_stacked <- do.call(rbind, lapply(c("full", "simplified"), function(mm) data.frame(model = mm, t(rowMeans(sapply(stk, function(s) ici(s$y, s[[mm]])))))))
ord <- lapply(stk, function(s) s[order(s$row_id), ]); avg <- ord[[1]][, c("row_id", "label", "db", "y")]
avg$full <- rowMeans(sapply(ord, `[[`, "full")); avg$simplified <- rowMeans(sapply(ord, `[[`, "simplified"))
avg$grp <- cut(avg$full, c(-Inf, RISK_CUTS, Inf), right = FALSE, labels = c("<10%", "10–30%", "≥30%"))
risk_groups <- do.call(rbind, lapply(split(avg, list(avg$label, avg$grp), drop = TRUE), function(s)
  data.frame(cluster = s$label[1], group = s$grp[1], n = nrow(s), deaths = sum(s$y), predicted = mean(s$full), observed = mean(s$y))))
risk_all <- do.call(rbind, lapply(split(avg, avg$grp), function(s) data.frame(group = s$grp[1], n = nrow(s), pct_of_patients = 100 * nrow(s) / nrow(avg),
  deaths = sum(s$y), predicted = mean(s$full), observed = mean(s$y), obs_lo = stats::binom.test(sum(s$y), nrow(s))$conf.int[1],
  obs_hi = stats::binom.test(sum(s$y), nrow(s))$conf.int[2])))
cm <- data.frame(cluster = names(folds), cluster_name = CL_LABEL[names(folds)], n = sapply(folds, `[[`, "n"), deaths = sapply(folds, `[[`, "deaths"),
                 t(sapply(folds, `[[`, "case_mix")),
                 C = plogis(sapply(folds, function(f) f$pooled$est[f$pooled$measure == "logit_C" & f$pooled$model == "full"])),
                 cal_slope = sapply(folds, function(f) f$pooled$est[f$pooled$measure == "cal_slope" & f$pooled$model == "full"]))
fold_models <- do.call(rbind, lapply(folds, function(f) data.frame(cluster = f$label, n_dev = f$n_dev, deaths_dev = f$deaths_dev,
  shrinkage = f$core$shrink, simplified_vars = paste(f$simp$vars, collapse = ", "), simplified_R2 = f$simp$r2)))

saveRDS(list(pooled = pooled, by_cluster = by_cluster, dca = dca, cal_curve = cal_curve, ici_stacked = ici_stacked, eligible = eligible,
             risk_groups = risk_groups, risk_all = risk_all, case_mix = cm, cl_label = CL_LABEL, fold_models = fold_models),
        file.path(DIR_LOCAL, paste0("03_iecv_results", sfx, ".rds")))
saveRDS(stk, file.path(DIR_LOCAL, paste0("03_iecv_stacked", sfx, ".rds")))
write_share(list(
  README = data.frame(note = c(sprintf("IECV: %d clusters; pooled (deaths >= %d): %s", length(cls), MIN_EVENTS_POOL, paste(eligible, collapse = ",")),
    "by_cluster: performance in each held-out cluster (Rubin's rules across imputations, normal-approximation 95% CI); pooled: random-effects pooling (REML + HKSJ), 95% PI with t(k-2)",
    "full = full model, simplified = simplified model; 'full|score' = full model in the same patients that have the score; dC_vs_score = full model minus score",
    "decision_curve, calibration_curve, ICI_stacked: stacked predictions of all held-out clusters, averaged across imputations; risk_groups: predicted and observed mortality in the three risk groups",
    "case_mix: case mix of each cluster (SD of the linear predictor, membership model C) and performance, used to explain heterogeneity", "This file contains summary statistics only")),
  pooled = pooled, by_cluster = by_cluster, risk_groups_all = risk_all, risk_groups_by_cluster = risk_groups, case_mix = cm,
  fold_models = fold_models, decision_curve = dca, calibration_curve = cal_curve, ICI_stacked = ici_stacked),
  paste0("03_iecv", sfx, ".xlsx"))
msg("03 done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min")
