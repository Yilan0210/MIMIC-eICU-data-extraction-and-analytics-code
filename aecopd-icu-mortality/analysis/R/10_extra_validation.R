# =====================================================================
# 10_extra_validation.R: additional validation V1–V9 (run after the final run of run_all.R; reads its saved results only and does not change the main analysis)
#   V1 temporal validation: performance of IECV held-out predictions by year; forward temporal validation (earlier years -> later years)
#   V2 hospital level: O/E funnel plots of eICU hospitals (Spiegelhalter 2005), random-intercept recalibration, agreement with APACHE IVa
#   V3 hospital type: teaching/non-teaching, bed size
#   V4 decomposition of case mix and coefficient differences: model-based C (van Klaveren 2016) vs observed C per cluster
#   V5 case definition: AECOPD as primary diagnosis or not; eICU identification method (ICD code / diagnosis string)
#   V6 uncertainty in clinical utility: Bayesian bootstrap probability of the best strategy and expected value of perfect information (EVPI)
#   V7 fairness: flexible calibration curves by sex, race and age
#   V8 reclassification: BAP-65 classes, APACHE IVa risk groups -> risk groups of the new model
#   V9 multiverse (specification curve): IECV pooled C and calibration slope for 216 specifications of 6 analytic choices
# Outputs: result/10_extra_validation.xlsx (summary); local/10_extra.rds (for figures; contains hospital-level data; local only)
# =====================================================================
source("00_setup.R"); source("iecv_functions.R")
t0 <- Sys.time()
EXTRA_B <- as.integer(Sys.getenv("AECOPD_EXTRA_B", "1000"))   # V6 Bayesian bootstrap samples
V4_B    <- as.integer(Sys.getenv("AECOPD_V4_B", "200"))       # V4 bootstrap samples
SPEC_M  <- as.integer(Sys.getenv("AECOPD_SPEC_M", "5"))       # V9 imputations per specification
RUN_V1_FORWARD <- Sys.getenv("AECOPD_V1_FORWARD", "1") == "1" # V1 forward temporal validation (slow: each design takes about as long as one IECV fold)

a   <- readRDS(file.path(DIR_LOCAL, "01_population.rds"))$a
d02 <- readRDS(file.path(DIR_LOCAL, "02_dev.rds"))
stk <- readRDS(file.path(DIR_LOCAL, "03_iecv_stacked.rds"))
raw <- load_clean()
EXTRA_COLS <- c("hospital_id", "teaching", "bed_cat", "year_grp", "coding_era", "aecopd_primary", "dx_by_icd", "dx_by_string")
miss_cols <- setdiff(EXTRA_COLS, names(raw))
if (length(miss_cols)) warning("Columns not in the snapshot; related analyses skipped: ", paste(miss_cols, collapse = ", "))
ex <- data.frame(row_id = a$row_id)
for (v in EXTRA_COLS) ex[[v]] <- if (v %in% names(raw)) raw[[v]][a$row_id] else NA
ex$hospital_id <- as.character(ex$hospital_id); ex$year_grp <- trimws(as.character(ex$year_grp))
ex$bed_cat <- trimws(as.character(ex$bed_cat))
num <- function(z) suppressWarnings(as.numeric(z))
for (v in c("teaching", "aecopd_primary", "dx_by_icd", "dx_by_string")) ex[[v]] <- num(ex[[v]])
attach_ex <- function(s) cbind(s, ex[match(s$row_id, ex$row_id), setdiff(names(ex), "row_id"), drop = FALSE])
stk <- lapply(stk, attach_ex)
a$year_grp <- ex$year_grp[match(a$row_id, ex$row_id)]

# IECV held-out predictions averaged across imputations (as in Figure 6 of 07)
ord <- lapply(stk, function(s) s[order(s$row_id), ])
avg <- ord[[1]][, c("row_id", "label", "db", "y", "sex", "race4", "age", EXTRA_COLS)]
for (m in c("full", SCORES)) avg[[m]] <- rowMeans(sapply(ord, `[[`, m))

# ---- General: C, O/E and calibration slope by group (computed per imputation, pooled with Rubin's rules) ----
perf_by <- function(fun, name, min_events = 5) {
  per <- do.call(rbind, lapply(seq_along(stk), function(j) { s <- stk[[j]]; s$lev <- fun(s); s <- s[!is.na(s$lev), ]
    do.call(rbind, lapply(split(s, s$lev), function(z) {
      if (sum(z$y) < min_events || sum(1 - z$y) < min_events) return(NULL)
      cs <- c_stats(z$y, z$full); cl <- calib(z$y, z$full)
      data.frame(j = j, level = z$lev[1], n = nrow(z), events = sum(z$y), lc = cs$logit_c, v = cs$var_logit_c,
                 oe = cl[["log_oe"]], v_oe = cl[["var_log_oe"]], sl = cl[["slope"]], v_sl = cl[["var_slope"]]) })) }))
  if (is.null(per) || !nrow(per)) return(NULL)
  do.call(rbind, lapply(split(per, per$level), function(z) {
    rc <- rubin(z$lc, z$v); ro <- rubin(z$oe, z$v_oe); rs <- rubin(z$sl, z$v_sl)
    data.frame(group = name, level = z$level[1], n = round(mean(z$n)), events = round(mean(z$events)),
               C = plogis(rc$est), C_lo = plogis(rc$lo), C_hi = plogis(rc$hi),
               OE = exp(ro$est), OE_lo = exp(ro$lo), OE_hi = exp(ro$hi), slope = rs$est, slope_lo = rs$lo, slope_hi = rs$hi) }))
}
auc_fast <- function(y, p) { r <- rank(p); n1 <- sum(y); n0 <- length(y) - n1; (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
# flexible calibration curve (logistic regression on a natural spline of the logit prediction; logistic recalibration when events are few) with 95% CI
flex_cal <- function(y, p, n_grid = 80) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6); lp <- qlogis(p); df <- if (sum(y) >= 50) 3 else 1
  f <- if (df > 1) stats::glm(y ~ splines::ns(lp, df = df), family = stats::binomial()) else stats::glm(y ~ lp, family = stats::binomial())
  rg <- stats::quantile(p, c(0.01, 0.99)); g <- seq(rg[[1]], rg[[2]], length.out = n_grid)
  pr <- stats::predict(f, newdata = data.frame(lp = qlogis(g)), se.fit = TRUE)
  data.frame(pred = g, obs = plogis(pr$fit), lo = plogis(pr$fit - 1.96 * pr$se.fit), hi = plogis(pr$fit + 1.96 * pr$se.fit))
}
out <- list(); share <- list()

# =====================================================================
# V1 temporal validation
# =====================================================================
msg("V1 temporal validation")
mid_year <- function(g) vapply(g, function(s) { y <- num(regmatches(s, gregexpr("[0-9]{4}", s))[[1]]); if (length(y)) mean(y) else NA_real_ }, numeric(1))
v1_drift <- perf_by(function(s) ifelse(is.na(s$year_grp) | s$year_grp == "", NA, paste(s$db, s$year_grp, sep = " | ")), "Period")
if (!is.null(v1_drift)) {
  v1_drift$db <- sub(" \\|.*$", "", v1_drift$level); v1_drift$period <- sub("^.*\\| ", "", v1_drift$level)
  v1_drift$year_mid <- mid_year(v1_drift$period); v1_drift <- v1_drift[order(v1_drift$db, v1_drift$year_mid), ]
}
v1_fwd <- NULL; v1_fwd_period <- NULL; v1_fwd_curves <- NULL
if (RUN_V1_FORWARD && any(!is.na(a$year_grp))) {
  mg <- sort(unique(a$year_grp[a$db == "MIMIC" & !is.na(a$year_grp)])); eg <- sort(unique(a$year_grp[a$db == "eICU" & !is.na(a$year_grp)]))
  designs <- list()
  if (length(mg) >= 3) designs$T1 <- list(label = sprintf("MIMIC-IV %s to %s (developed on all other patients)", sub(" - .*$", "", mg[length(mg) - 1]),
                                                          sub("^.* - ", "", mg[length(mg)])),
                                          hold = a$db == "MIMIC" & a$year_grp %in% tail(mg, 2))
  if (length(eg) >= 2) designs$T2 <- list(label = sprintf("eICU-CRD %s (developed on all other patients)", tail(eg, 1)),
                                          hold = a$db == "eICU" & a$year_grp %in% tail(eg, 1))
  designs <- Filter(function(d) sum(a$death[d$hold]) >= MIN_EVENTS_POOL, designs)
  run_t <- function(nm) { d <- designs[[nm]]; dev <- a[!d$hold, ]; dev$cluster <- droplevels(dev$cluster)
    tryCatch(fold_fit(dev, a[d$hold, ], label = nm, hold_with_outcome = FALSE, seed = SEED + 7000 + which(names(designs) == nm), keep_tested = FALSE),
             error = function(e) { message("V1 ", nm, " failed: ", conditionMessage(e)); NULL }) }
  fw <- if (N_CORES > 1 && length(designs) > 1) parallel::mclapply(names(designs), run_t, mc.cores = min(N_CORES, length(designs)), mc.preschedule = FALSE)
        else lapply(names(designs), run_t)
  names(fw) <- names(designs); fw <- Filter(function(f) !is.null(f) && !inherits(f, "try-error"), fw)
  if (length(fw)) {
    v1_fwd <- do.call(rbind, lapply(names(fw), function(nm) { tb <- fold_table(fw[[nm]]); tb$design <- designs[[nm]]$label
      tb[tb$model %in% c("full", "simplified", SCORES, paste0("full-", SCORES)) & tb$measure %in% c("logit_C", "log_OE", "cal_slope", "CITL", "dC_vs_score"), ] }))
    v1_fwd_period <- do.call(rbind, lapply(names(fw), function(nm) { st <- lapply(fw[[nm]]$stacked, attach_ex)
      per <- do.call(rbind, lapply(seq_along(st), function(j) { s <- st[[j]]
        do.call(rbind, lapply(split(s, s$year_grp), function(z) { if (sum(z$y) < 5) return(NULL)
          cs <- c_stats(z$y, z$full); cl <- calib(z$y, z$full)
          data.frame(j = j, period = z$year_grp[1], n = nrow(z), events = sum(z$y), lc = cs$logit_c, v = cs$var_logit_c, oe = cl[["log_oe"]], v_oe = cl[["var_log_oe"]]) })) }))
      if (is.null(per)) return(NULL)
      do.call(rbind, lapply(split(per, per$period), function(z) { rc <- rubin(z$lc, z$v); ro <- rubin(z$oe, z$v_oe)
        data.frame(design = designs[[nm]]$label, period = z$period[1], n = round(mean(z$n)), events = round(mean(z$events)),
                   C = plogis(rc$est), C_lo = plogis(rc$lo), C_hi = plogis(rc$hi), OE = exp(ro$est), OE_lo = exp(ro$lo), OE_hi = exp(ro$hi)) })) }))
    v1_fwd_curves <- do.call(rbind, lapply(names(fw), function(nm) { st <- fw[[nm]]$stacked
      o <- lapply(st, function(s) s[order(s$row_id), ]); y <- o[[1]]$y; p <- rowMeans(sapply(o, `[[`, "full"))
      cbind(design = designs[[nm]]$label, flex_cal(y, p)) }))
  }
}
share$V1_drift <- v1_drift; share$V1_forward <- v1_fwd; share$V1_forward_by_period <- v1_fwd_period
out$v1 <- list(drift = v1_drift, fwd = v1_fwd, fwd_period = v1_fwd_period, fwd_curves = v1_fwd_curves)

# =====================================================================
# V2 hospital level (eICU): funnel plots, random-intercept recalibration, agreement with APACHE IVa
# =====================================================================
msg("V2 hospital level")
funnel_lim <- function(E, p) {   # Spiegelhalter 2005: Poisson distribution, interpolated control limits (O/E scale)
  r <- stats::qpois(p, E); al <- (stats::ppois(r, E) - p) / stats::dpois(r, E); pmax((r - al) / E, 0) }
flag3 <- function(O, E) ifelse(O / E > funnel_lim(E, 0.975), "Higher than expected", ifelse(O / E < funnel_lim(E, 0.025), "Lower than expected", "As expected"))
hosp_tab <- function(d, pcol) {
  d <- d[!is.na(d[[pcol]]), ]
  h <- do.call(rbind, lapply(split(d, d$hospital_id), function(z) data.frame(hospital_id = z$hospital_id[1], n = nrow(z), O = sum(z$y), E = sum(z[[pcol]]))))
  h$SMR <- h$O / h$E; h$flag95 <- flag3(h$O, h$E)
  h$out998 <- h$SMR > funnel_lim(h$E, 0.999) | h$SMR < funnel_lim(h$E, 0.001)
  h
}
overdisp <- function(h) { z <- (h$O - h$E) / sqrt(h$E); q <- stats::quantile(z, c(0.1, 0.9)); zw <- pmin(pmax(z, q[[1]]), q[[2]]); mean(zw^2) }
mixed_cal <- function(d, pcol, label) {
  d <- d[!is.na(d[[pcol]]), ]; d$off <- qlogis(pmin(pmax(d[[pcol]], 1e-6), 1 - 1e-6))
  base <- data.frame(model = label, n = nrow(d), hospitals = length(unique(d$hospital_id)), events = sum(d$y))
  if (!requireNamespace("lme4", quietly = TRUE)) return(cbind(base, note = "lme4 not installed"))
  f0 <- stats::glm(y ~ 1 + offset(off), data = d, family = stats::binomial())
  f1 <- tryCatch(suppressMessages(lme4::glmer(y ~ 1 + offset(off) + (1 | hospital_id), data = d, family = stats::binomial())), error = function(e) NULL)
  if (is.null(f1)) return(cbind(base, note = "glmer did not converge"))
  sdh <- as.numeric(attr(lme4::VarCorr(f1)$hospital_id, "stddev"))
  ci <- tryCatch(suppressWarnings(suppressMessages(lme4::confint.merMod(f1, parm = "theta_", method = "profile", oldNames = FALSE))), error = function(e) NULL)
  lr <- max(0, 2 * (as.numeric(stats::logLik(f1)) - as.numeric(stats::logLik(f0))))
  cbind(base, intercept = unname(lme4::fixef(f1)[1]), intercept_se = sqrt(diag(as.matrix(stats::vcov(f1))))[1],
        sd_hospital = sdh, sd_lo = if (is.null(ci)) NA else ci[1, 1], sd_hi = if (is.null(ci)) NA else ci[1, 2],
        MOR = exp(sqrt(2) * sdh * stats::qnorm(0.75)), LR_p = 0.5 * stats::pchisq(lr, 1, lower.tail = FALSE), note = "")
}
v2 <- NULL
de <- avg[avg$db == "eICU" & !is.na(avg$hospital_id) & avg$hospital_id != "", ]
if (nrow(de) > 0 && length(unique(de$hospital_id)) >= 10) {
  hm <- hosp_tab(de, "full")
  dap <- de[!is.na(de$APACHE_IVa), ]
  hm_ap <- hosp_tab(dap, "full"); ha_ap <- hosp_tab(dap, "APACHE_IVa")
  agree <- merge(hm_ap[, c("hospital_id", "n", "O", "E", "SMR", "flag95")], ha_ap[, c("hospital_id", "E", "SMR", "flag95")], by = "hospital_id", suffixes = c("_model", "_apache"))
  big <- agree$E_model >= 2 & agree$E_apache >= 2
  v2_sum <- rbind(
    data.frame(model = "New model", patients = "All eICU-CRD patients", hospitals = nrow(hm), events = sum(hm$O), expected = sum(hm$E),
               pct_outside_95 = 100 * mean(hm$flag95 != "As expected"), pct_outside_998 = 100 * mean(hm$out998), overdispersion = overdisp(hm)),
    data.frame(model = "New model", patients = "Patients with APACHE IVa", hospitals = nrow(hm_ap), events = sum(hm_ap$O), expected = sum(hm_ap$E),
               pct_outside_95 = 100 * mean(hm_ap$flag95 != "As expected"), pct_outside_998 = 100 * mean(hm_ap$out998), overdispersion = overdisp(hm_ap)),
    data.frame(model = "APACHE IVa", patients = "Patients with APACHE IVa", hospitals = nrow(ha_ap), events = sum(ha_ap$O), expected = sum(ha_ap$E),
               pct_outside_95 = 100 * mean(ha_ap$flag95 != "As expected"), pct_outside_998 = 100 * mean(ha_ap$out998), overdispersion = overdisp(ha_ap)))
  v2_mixed <- rbind(mixed_cal(de, "full", "New model, all eICU-CRD patients"), mixed_cal(dap, "full", "New model, patients with APACHE IVa"),
                    mixed_cal(dap, "APACHE_IVa", "APACHE IVa, patients with APACHE IVa"))
  tabf <- as.data.frame.matrix(table(model = factor(agree$flag95_model, c("Lower than expected", "As expected", "Higher than expected")),
                                     apache = factor(agree$flag95_apache, c("Lower than expected", "As expected", "Higher than expected"))))
  tabf <- cbind(new_model_classification = rownames(tabf), tabf); rownames(tabf) <- NULL
  v2_agree <- data.frame(item = c("Hospitals compared", "Hospitals with >=2 expected deaths under both", "Spearman correlation of O/E (>=2 expected deaths)",
                                  "Same funnel classification (%)"),
                         value = c(nrow(agree), sum(big), if (sum(big) >= 5) round(stats::cor(agree$SMR_model[big], agree$SMR_apache[big], method = "spearman"), 3) else NA,
                                   round(100 * mean(agree$flag95_model == agree$flag95_apache), 1)))
  share$V2_summary <- v2_sum; share$V2_random_intercept <- v2_mixed; share$V2_agreement <- v2_agree; share$V2_flag_crosstab <- tabf
  v2 <- list(hm = hm, hm_ap = hm_ap, ha_ap = ha_ap, agree = agree, sum = v2_sum, mixed = v2_mixed)
} else share$V2_summary <- data.frame(note = "No eICU hospital IDs; V2 not run")
out$v2 <- v2

# =====================================================================
# V3 hospital type, V5 case definition (IECV held-out predictions; C, O/E and slope by group)
# =====================================================================
msg("V3, V5 hospital type and case definition")
v3 <- rbind(
  perf_by(function(s) ifelse(s$db == "eICU" & s$teaching %in% 1, "Teaching", ifelse(s$db == "eICU" & s$teaching %in% 0, "Non-teaching", NA)), "Teaching status (eICU-CRD)"),
  perf_by(function(s) { b <- ifelse(s$db == "eICU" & !is.na(s$bed_cat) & s$bed_cat != "", s$bed_cat, NA)
    lo <- ifelse(grepl("^\\s*<", b), 0, num(sub("^[^0-9]*([0-9]+).*$", "\\1", b)))   # lower bound of the bed category: "<100" coded 0
    ifelse(is.na(lo), NA, ifelse(lo < 250, "<250 beds", ">=250 beds")) }, "Hospital size (eICU-CRD)"))
v5 <- rbind(
  perf_by(function(s) ifelse(s$aecopd_primary %in% 1, "Primary diagnosis", ifelse(s$aecopd_primary %in% 0, "Secondary diagnosis", NA)), "AECOPD as diagnosis"),
  perf_by(function(s) ifelse(s$db != "eICU", NA, ifelse(s$dx_by_icd %in% 1 & s$dx_by_string %in% 1, "ICD code and diagnosis string",
                             ifelse(s$dx_by_icd %in% 1, "ICD code only", ifelse(s$dx_by_string %in% 1, "Diagnosis string only", NA)))), "Case identification (eICU-CRD)"))
share$V3_hospital_type <- v3; share$V5_case_definition <- v5; out$v3 <- v3; out$v5 <- v5

# =====================================================================
# V4 model-based C (van Klaveren 2016): expected C if the model is correct in the cluster, determined by case mix only
# =====================================================================
msg("V4 case-mix decomposition")
mbc <- function(p) { o <- outer(p, p, ">") + 0.5 * outer(p, p, "=="); diag(o) <- 0; w <- outer(p, 1 - p); diag(w) <- 0; sum(o * w) / sum(w) }
set.seed(SEED + 400)
v4 <- do.call(rbind, lapply(split(avg, avg$label), function(z) {
  ob <- c_stats(z$y, z$full); m0 <- mbc(z$full)
  bt <- t(replicate(V4_B, { i <- sample.int(nrow(z), replace = TRUE); c(auc_fast(z$y[i], z$full[i]), mbc(z$full[i])) }))
  dd <- bt[, 1] - bt[, 2]
  data.frame(cluster = z$label[1], db = z$db[1], n = nrow(z), events = sum(z$y), C_observed = ob$c, C_obs_lo = plogis(ob$logit_c - 1.96 * sqrt(ob$var_logit_c)),
             C_obs_hi = plogis(ob$logit_c + 1.96 * sqrt(ob$var_logit_c)), C_model_based = m0,
             C_mb_lo = unname(stats::quantile(bt[, 2], 0.025)), C_mb_hi = unname(stats::quantile(bt[, 2], 0.975)),
             difference = ob$c - m0, diff_lo = unname(stats::quantile(dd, 0.025)), diff_hi = unname(stats::quantile(dd, 0.975)),
             sd_lp = stats::sd(qlogis(z$full))) }))
r03 <- readRDS(file.path(DIR_LOCAL, "03_iecv_results.rds"))
sl <- r03$by_cluster[r03$by_cluster$measure == "cal_slope" & r03$by_cluster$model == "full", c("label", "est")]
v4$cal_slope <- sl$est[match(as.character(v4$cluster), as.character(sl$label))]
v4$cluster_name <- r03$cl_label[as.character(v4$cluster)]
share$V4_casemix_decomposition <- v4; out$v4 <- v4

# =====================================================================
# V6 Bayesian bootstrap: probability that each strategy is best at each threshold, probability that the model is useful, EVPI
# =====================================================================
msg("V6 value of information, B =", EXTRA_B)
bb_nb <- function(y, P, grid = DCA_GRID, B = EXTRA_B, seed = SEED + 600) {
  set.seed(seed); n <- length(y); W <- matrix(stats::rexp(B * n), B, n); W <- W / rowSums(W)
  odds <- matrix(grid / (1 - grid), B, length(grid), byrow = TRUE); res <- list()
  for (s in colnames(P)) { pos <- outer(P[, s], grid, ">="); res[[s]] <- (W %*% (pos * y)) - (W %*% (pos * (1 - y))) * odds }
  res$treat_all <- matrix(drop(W %*% y), B, length(grid)) - matrix(drop(W %*% (1 - y)), B, length(grid)) * odds
  res$treat_none <- matrix(0, B, length(grid)); res
}
voi <- function(y, P, label) {
  ok <- stats::complete.cases(P); y <- y[ok]; P <- P[ok, , drop = FALSE]
  nb <- bb_nb(y, P); S <- c("treat_none", "treat_all", colnames(P))          # ties go to the default strategy
  arr <- simplify2array(nb[S]); dimnames(arr)[[3]] <- S                     # B × T × S
  do.call(rbind, lapply(seq_along(DCA_GRID), function(t) {
    m <- arr[, t, , drop = TRUE]; best <- S[apply(m, 1, which.max)]; en <- colMeans(m)
    dflt <- m[, c("full", "treat_all", "treat_none"), drop = FALSE]
    data.frame(subset = label, threshold = DCA_GRID[t], strategy = S, expected_NB = en,
               p_best = vapply(S, function(s) mean(best == s), numeric(1)),
               p_model_useful = mean(m[, "full"] > pmax(m[, "treat_all"], 0)),
               EVPI_all = mean(apply(m, 1, max)) - max(en),
               EVPI_model_vs_defaults = mean(apply(dflt, 1, max)) - max(colMeans(dflt)), row.names = NULL) }))
}
v6 <- rbind(voi(avg$y, as.matrix(avg[, c("full", "BAP65", "DECAF_mod")]), "All patients"),
            { z <- avg[!is.na(avg$APACHE_IVa), ]; voi(z$y, as.matrix(z[, c("full", "APACHE_IVa")]), "eICU-CRD, APACHE IVa available") },
            { z <- avg[avg$db == "MIMIC" & !is.na(avg$SAPSII) & !is.na(avg$SOFA), ]; if (nrow(z) > 100) voi(z$y, as.matrix(z[, c("full", "SAPSII", "SOFA")]), "MIMIC-IV") })
v6_key <- unique(v6[v6$threshold %in% c(0.10, 0.20, 0.30), c("subset", "threshold", "p_model_useful", "EVPI_all", "EVPI_model_vs_defaults")])
v6_key <- merge(v6_key, v6[v6$strategy == "full" & v6$threshold %in% c(0.10, 0.20, 0.30), c("subset", "threshold", "p_best")], by = c("subset", "threshold"))
names(v6_key)[names(v6_key) == "p_best"] <- "p_model_best_of_all"
share$V6_key_thresholds <- v6_key; share$V6_full_grid <- v6; out$v6 <- v6

# =====================================================================
# V7 fairness: flexible calibration curves and measures by sex, race and age
# =====================================================================
msg("V7 fairness")
grp7 <- list(Sex = function(s) ifelse(s$sex == "F", "Female", "Male"), `Race and ethnicity` = function(s) s$race4,
             Age = function(s) ifelse(s$age >= 80, ">=80 years", "<80 years"))
v7_curves <- do.call(rbind, lapply(names(grp7), function(g) { lev <- grp7[[g]](avg)
  do.call(rbind, lapply(sort(unique(stats::na.omit(lev))), function(l) { z <- avg[lev %in% l, ]
    if (sum(z$y) < 10 || nrow(z) < 50) return(NULL)
    cbind(group = g, level = l, n = nrow(z), events = sum(z$y), flex_cal(z$y, z$full)) })) }))
v7_tab <- do.call(rbind, lapply(names(grp7), function(g) { pb <- perf_by(grp7[[g]], g); lev <- grp7[[g]](avg)
  if (is.null(pb)) return(NULL)
  pb$ICI <- vapply(pb$level, function(l) { z <- avg[lev %in% l, ]; ici(z$y, z$full)[["ICI"]] }, numeric(1)); pb }))
share$V7_fairness <- v7_tab; out$v7 <- list(curves = v7_curves, tab = v7_tab)

# =====================================================================
# V8 reclassification: BAP-65 class (mode across imputations), APACHE IVa risk groups -> risk groups of the new model
# =====================================================================
msg("V8 reclassification")
cut3 <- function(p) cut(p, c(-Inf, RISK_CUTS, Inf), labels = c("<10%", "10-30%", ">=30%"), right = FALSE)
B65 <- sapply(d02$xs, function(x) x$bap65[match(avg$row_id, x$row_id)])
avg$bap_class <- factor(c("I", "II", "III", "IV", "V")[apply(as.matrix(B65), 1, function(r) as.integer(names(which.max(table(r)))))], levels = c("I", "II", "III", "IV", "V"))
avg$grp_model <- cut3(avg$full); avg$grp_apache <- cut3(avg$APACHE_IVa)
flow <- function(d, from, label) { d <- d[!is.na(d[[from]]), ]
  r <- stats::aggregate(cbind(n = 1, deaths = y) ~ f + t, data = data.frame(f = d[[from]], t = d$grp_model, y = d$y), FUN = sum)
  data.frame(comparison = label, from = r$f, to = r$t, n = r$n, deaths = r$deaths, mortality = r$deaths / r$n) }
v8 <- rbind(flow(avg, "bap_class", "BAP-65 class -> new model risk group"),
            flow(avg[avg$db == "eICU", ], "grp_apache", "APACHE IVa risk group -> new model risk group (eICU-CRD)"))
v8_marg <- rbind(
  do.call(rbind, lapply(split(avg, avg$bap_class), function(z) if (nrow(z)) data.frame(comparison = "BAP-65 class", level = as.character(z$bap_class[1]), n = nrow(z), deaths = sum(z$y)))),
  do.call(rbind, lapply(split(avg[!is.na(avg$grp_apache), ], avg$grp_apache[!is.na(avg$grp_apache)]), function(z) if (nrow(z)) data.frame(comparison = "APACHE IVa risk group", level = as.character(z$grp_apache[1]), n = nrow(z), deaths = sum(z$y)))),
  do.call(rbind, lapply(split(avg, avg$grp_model), function(z) if (nrow(z)) data.frame(comparison = "New model risk group", level = as.character(z$grp_model[1]), n = nrow(z), deaths = sum(z$y)))))
v8_marg$mortality <- v8_marg$deaths / v8_marg$n
share$V8_reclassification <- v8; share$V8_margins <- v8_marg; out$v8 <- list(flow = v8, marg = v8_marg)

# =====================================================================
# V9 multiverse: imputation x functional form x estimation x added predictor x outcome x population (2 x 3 x 3 x 3 x 2 x 2 = 216)
#   each specification evaluated by 7-cluster IECV (clusters with >= 10 deaths); Rubin's rules within clusters, then REML pooling of C and calibration slope
#   to limit computation: imputations from 02 (first SPEC_M); knots from all data; ML with heuristic shrinkage (van Houwelingen)
# =====================================================================
msg("V9 multiverse, imputations", SPEC_M)
K4 <- local({ d <- derive(a); setNames(lapply(SPLINE_VARS, function(v) unname(stats::quantile(d[[v]], c(0.05, 0.35, 0.65, 0.95), na.rm = TRUE))), SPLINE_VARS) })
add4 <- function(x) { for (v in SPLINE_VARS) { b <- Hmisc::rcspline.eval(x[[v]], knots = K4[[v]], inclx = FALSE)
  x[[paste0(v, "_s1")]] <- b[, 1]; x[[paste0(v, "_s2")]] <- b[, 2] }; x }
TERMS <- list(`3-knot splines` = CORE_TERMS, Linear = c(SPLINE_VARS, LINEAR_VARS),
              `4-knot splines` = c(as.vector(rbind(SPLINE_VARS, paste0(SPLINE_VARS, "_s1"), paste0(SPLINE_VARS, "_s2"))), LINEAR_VARS))
ADD <- list(None = character(), `+ Lactate` = "ln_lac", `+ Vasopressors` = "vaso")
grid9 <- expand.grid(missing = c("Multiple imputation", "Complete cases"), form = names(TERMS), estimation = c("ML + shrinkage", "Ridge", "LASSO"),
                     added = names(ADD), outcome = c("In-hospital death", "Death or hospice discharge"),
                     population = c("All patients", "Excluding treatment limitation"), stringsAsFactors = FALSE)
xs_mi <- lapply(d02$xs[seq_len(min(SPEC_M, length(d02$xs)))], add4)
NEED_CC <- c("age", "sex", "bun", "gcs", "hr", "rr", "map", "temp", "ph", "paco2", "hgb", "plt", "wbc", "glucose", "af", "imv")
cc_data <- function(added) { need <- c(NEED_CC, if ("ln_lac" %in% ADD[[added]]) "lactate", if ("vaso" %in% ADD[[added]]) "vaso")
  acc <- a[stats::complete.cases(a[, need]), ]; acc$cluster <- droplevels(acc$cluster); add4(add_design(derive(acc), d02$knots)) }
fit_predict <- function(est, terms, dev, ho, y, seed) {
  if (est == "ML + shrinkage") {
    f <- fit_glm(terms, dev, y); chi <- f$null.deviance - f$deviance; s <- max(0, min(1, (chi - length(terms)) / chi))
    b <- stats::coef(f)[-1] * s; b[is.na(b)] <- 0
    off <- drop(as.matrix(dev[, terms]) %*% b)
    int <- stats::coef(stats::glm(dev[[y]] ~ 1 + offset(off), family = stats::binomial()))[1]
    return(plogis(int + drop(as.matrix(ho[, terms]) %*% b)))
  }
  set.seed(seed)
  cv <- glmnet::cv.glmnet(as.matrix(dev[, terms]), dev[[y]], family = "binomial", alpha = if (est == "LASSO") 1 else 0, nfolds = 5)
  drop(stats::predict(cv, as.matrix(ho[, terms]), s = "lambda.min", type = "response"))
}
eval_spec <- function(i) {
  sp <- grid9[i, ]; terms <- c(TERMS[[sp$form]], ADD[[sp$added]]); y <- if (sp$outcome == "In-hospital death") "death" else "death_s3"
  dl <- if (sp$missing == "Multiple imputation") xs_mi else list(cc_data(sp$added))
  dl <- lapply(dl, function(x) { x <- x[!is.na(x[[y]]), ]; if (sp$population != "All patients") x <- x[!(x$limit_s18 %in% 1), ]; x })
  per <- list()
  for (j in seq_along(dl)) { x <- dl[[j]]
    for (k in levels(droplevels(x$cluster))) {
      dev <- x[x$cluster != k, ]; ho <- x[x$cluster == k, ]
      if (sum(ho[[y]]) < MIN_EVENTS_POOL || sum(1 - ho[[y]]) < MIN_EVENTS_POOL) next
      p <- fit_predict(sp$estimation, terms, dev, ho, y, SEED + 9000 + 100 * i + 10 * j + as.integer(k))
      cs <- c_stats(ho[[y]], p); cl <- calib(ho[[y]], p)
      per[[length(per) + 1]] <- data.frame(j = j, k = k, lc = cs$logit_c, v = cs$var_logit_c, sl = cl[["slope"]], vs = cl[["var_slope"]])
    } }
  per <- do.call(rbind, per)
  base <- cbind(spec = i, sp, n = nrow(dl[[1]]), events = sum(dl[[1]][[y]]))
  if (is.null(per) || length(unique(per$k)) < 3) return(cbind(base, k = 0, C = NA, C_lo = NA, C_hi = NA, slope = NA, slope_lo = NA, slope_hi = NA))
  byk <- do.call(rbind, lapply(split(per, per$k), function(z) { r1 <- rubin(z$lc, z$v); r2 <- rubin(z$sl, z$vs)
    data.frame(lc = r1$est, se = r1$se, sl = r2$est, sse = r2$se) }))
  pc <- re_pool(byk$lc, byk$se, plogis, "C"); ps <- re_pool(byk$sl, byk$sse, identity, "slope")
  cbind(base, k = nrow(byk), C = pc$est, C_lo = pc$lo, C_hi = pc$hi, slope = ps$est, slope_lo = ps$lo, slope_hi = ps$hi)
}
safe_spec <- function(i) tryCatch(eval_spec(i), error = function(e) { message("V9 specification ", i, " failed: ", conditionMessage(e))
  cbind(spec = i, grid9[i, ], n = NA, events = NA, k = 0, C = NA, C_lo = NA, C_hi = NA, slope = NA, slope_lo = NA, slope_hi = NA) })
SPEC_N <- as.integer(Sys.getenv("AECOPD_SPEC_N", "0"))       # test runs only: if > 0, only the reference specification and a random sample of specifications are run
idx9 <- if (SPEC_N > 0 && SPEC_N < nrow(grid9)) { set.seed(SEED); sort(unique(c(1, sample(2:nrow(grid9), SPEC_N - 1)))) } else seq_len(nrow(grid9))
v9 <- do.call(rbind, if (N_CORES > 1) parallel::mclapply(idx9, safe_spec, mc.cores = N_CORES, mc.preschedule = TRUE)
                     else lapply(idx9, safe_spec))
v9$primary <- v9$spec == 1
okc <- is.finite(v9$C); c0 <- v9$C[v9$primary]
v9_sum <- data.frame(item = c("Specifications", "Specifications evaluated", "Median pooled C", "Range of pooled C", "Within ±0.01 of the reference specification (%)",
                              "Median calibration slope", "Range of calibration slope", "Reference specification C (heuristic shrinkage, imputations from full data)"),
                     value = c(nrow(v9), sum(okc), round(stats::median(v9$C[okc]), 3), paste(round(range(v9$C[okc]), 3), collapse = " to "),
                               round(100 * mean(abs(v9$C[okc] - c0) <= 0.01), 1), round(stats::median(v9$slope[okc]), 2),
                               paste(round(range(v9$slope[okc]), 2), collapse = " to "), round(c0, 3)))
share$V9_specification_summary <- v9_sum; share$V9_specifications <- v9; out$v9 <- v9

# =====================================================================
saveRDS(out, file.path(DIR_LOCAL, "10_extra.rds"))
readme <- c(
  "Additional validation V1-V9; all based on IECV held-out predictions or re-fitted validation models; the main analysis is not changed; same patients as the main analysis (4,417)",
  "V1: drift = C, O/E and slope of IECV held-out predictions by database and year group; forward = model built in all other patients (complete process: imputation, knots, ML + bootstrap shrinkage) and validated in later years",
  "V2: eICU hospitals; O/E funnel plot control limits as in Spiegelhalter 2005 (Poisson, interpolated); overdispersion = mean of z^2 after 10% winsorisation (>1 indicates between-hospital variation beyond chance); random_intercept = glmer(y ~ offset(logit p) + (1|hospital))",
  "V3, V5: C, O/E and slope of IECV held-out predictions by group (Rubin's rules across imputations)",
  sprintf("V4: model-based C (van Klaveren 2016) = expected C if the model is correct in the cluster, determined by case mix only; difference = observed C - model-based C (bootstrap %d samples, 95%% interval); clearly negative values suggest the coefficients do not hold in that cluster", V4_B),
  sprintf("V6: Bayesian bootstrap %d samples; p_best = probability that the strategy has the highest net benefit; p_model_useful = probability that the model has higher net benefit than 'treat all' and 'treat none'; EVPI = E[max NB] - max E[NB] (net benefit units, x1000 = per 1000 patients)", EXTRA_B),
  "V7: flexible calibration curve = natural spline of logit(prediction) (3 df; logistic recalibration when events < 50)",
  "V8: BAP-65 class = mode across imputed datasets; risk groups <10%, 10-30%, >=30%",
  sprintf("V9: 216 specifications; each evaluated by 7-cluster IECV (first %d imputations of the full data, knots from all data, ML with heuristic shrinkage), so the C of the reference specification differs slightly from the main analysis", SPEC_M),
  "This file contains summary statistics only; hospital-level data are saved only in local/10_extra.rds")
write_share(c(list(README = data.frame(note = readme)), Filter(Negate(is.null), share)), "10_extra_validation.xlsx")
msg("10 done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min")
