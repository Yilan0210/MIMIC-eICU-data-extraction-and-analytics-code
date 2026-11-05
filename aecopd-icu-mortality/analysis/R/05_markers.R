# =====================================================================
# 05_markers.R: secondary questions Q3 (added value of markers) and Q4 (ratio forms), only in the subset with albumin and differential measured
#   Base model: intercept + linear predictor LP of the core model (slope estimated freely)
#   Extension 1: + ln neutrophils + ln lymphocytes + albumin (primary; assessed with the decision rule); extension 2: + eosinophils < 0.05; extension 3: 3-knot splines of the three
#   Measures: Meng-Rubin likelihood ratio test, fraction of new information, bootstrap-corrected dC, IECV pooled dC (LP from each development fold's own core model), change in scaled Brier, categorical NRI
#   H1b: added value of the markers over APACHE IVa (eICU) and SAPS II and SOFA (MIMIC)
#   Q4: NLR and PNI linear constraints; 7 indices compared with their components; inverse probability weighting
# =====================================================================
source("00_setup.R")
t0 <- Sys.time()
d02 <- readRDS(file.path(DIR_LOCAL, "02_dev.rds")); core <- d02$core
folds <- readRDS(file.path(DIR_LOCAL, "03_iecv_folds.rds"))
D1_test <- function(fits, idx) {   # Multiparameter Wald test (Li et al. 1991)
  m <- length(fits); k <- length(idx)
  Q <- matrix(sapply(fits, function(f) stats::coef(f)[idx]), nrow = k)
  qb <- rowMeans(Q); Ub <- Reduce(`+`, lapply(fits, function(f) stats::vcov(f)[idx, idx, drop = FALSE])) / m
  Bm <- if (m > 1) stats::cov(t(Q)) else Ub * 0; if (k == 1) Bm <- matrix(stats::var(Q[1, ]), 1, 1)
  r <- (1 + 1 / m) * sum(diag(Bm %*% solve(Ub))) / k
  D <- drop(t(qb) %*% solve(Ub) %*% qb) / (k * (1 + r)); t <- k * (m - 1)
  v <- if (r <= 0) Inf else if (t > 4) 4 + (t - 4) * (1 + (1 - 2 / t) / r)^2 else t * (1 + 1 / k) * (1 + 1 / r)^2 / 2
  if (is.finite(v)) stats::pf(D, k, v, lower.tail = FALSE) else stats::pchisq(D * k, k, lower.tail = FALSE)
}
xt <- lapply(d02$xs, function(x) { x <- x[x$tested == 1, ]; x$lp <- qlogis(predict_dev(core, x)); x })
MK_SPL <- c("ln_n", "ln_l", "alb")
kmk <- setNames(lapply(MK_SPL, function(v) unname(stats::quantile(xt[[1]][[v]], KNOT_P, na.rm = TRUE))), MK_SPL)
add_mk <- function(x, k = kmk) { for (v in names(k)) x[[paste0(v, "_nl")]] <- rcs_nl(x[[v]], k[[v]]); x }
xt <- lapply(xt, add_mk)
MODS <- list(base = "lp", E1 = c("lp", MARK_E1), E2 = c("lp", MARK_E2), E3 = c("lp", MARK_SPL))
msg("Marker subset: ", nrow(xt[[1]]), "patients, deaths", sum(xt[[1]]$death))

# ---- 1. Likelihood ratio tests, bootstrap-corrected dC ----
msg("1. Likelihood ratio tests and bootstrap correction (B =", B_BOOT, ")")
vals <- lapply(MODS, function(tm) t(sapply(seq_along(xt), function(i) validate_one(tm, xt[[i]], seed = SEED + 30 + i))))
pbar <- mean(xt[[1]]$death)
delta <- do.call(rbind, lapply(c("E1", "E2", "E3"), function(e) {
  lr <- D3_test(xt, MODS[[e]], MODS$base)
  per <- t(sapply(seq_along(xt), function(i) { x <- xt[[i]]
    P <- cbind(stats::fitted(fit_glm(MODS$base, x)), stats::fitted(fit_glm(MODS[[e]], x))); cs <- c_stats(x$death, P)
    opt <- (vals[[e]][i, "c_app"] - vals[[e]][i, "c_corr"]) - (vals$base[i, "c_app"] - vals$base[i, "c_corr"])
    c(dc_app = cs$dc, dc_corr = cs$dc - unname(opt), var = cs$var_dc) }))
  ra <- rubin(per[, "dc_app"], per[, "var"]); rc <- rubin(per[, "dc_corr"], per[, "var"])
  sb <- function(m) mean(1 - vals[[m]][, "brier_corr"] / (pbar * (1 - pbar)))
  data.frame(comparison = paste(e, "vs base"), df = lr$df1, LR_F = lr$F, LR_p = lr$p, median_LR_chi2 = lr$median_LR_chi2,
             new_info_fraction = lr$new_info_fraction, dC_apparent = ra$est, dC_app_lo = ra$lo, dC_app_hi = ra$hi,
             dC_corrected = rc$est, dC_corr_lo = rc$lo, dC_corr_hi = rc$hi, d_scaled_Brier_corrected = sb(e) - sb("base"))
}))

# ---- 2. IECV: fit in the development part of each fold (ML, <= 8 parameters), evaluate in the held-out part ----
msg("2. IECV")
iecv_rows <- list(); stk <- vector("list", M_IMP)
for (k in names(folds)) {
  td <- folds[[k]]$tested; if (is.null(td)) next
  kk <- setNames(lapply(MK_SPL, function(v) unname(stats::quantile(td$dev[[1]][[v]], KNOT_P, na.rm = TRUE))), MK_SPL)
  for (j in seq_len(M_IMP)) {
    dv <- add_mk(td$dev[[j]], kk); ho <- add_mk(td$hold[[j]], kk)
    P <- sapply(MODS, function(tm) stats::predict(fit_glm(tm, dv), newdata = ho, type = "response"))
    cs <- c_stats(ho$death, P)
    for (e in 2:4) iecv_rows[[length(iecv_rows) + 1]] <- data.frame(cluster = k, j = j, comparison = paste(names(MODS)[e], "vs base"),
                                                                       est = cs$dc[e - 1], var = cs$var_dc[e - 1])
    stk[[j]] <- rbind(stk[[j]], data.frame(row_id = ho$row_id, cluster = k, y = ho$death, P))
  }
}
ir <- do.call(rbind, iecv_rows)
byk <- do.call(rbind, lapply(split(ir, list(ir$cluster, ir$comparison), drop = TRUE), function(s) { r <- rubin(s$est, s$var)
  data.frame(cluster = s$cluster[1], comparison = s$comparison[1], deaths = sum(folds[[s$cluster[1]]]$tested$hold[[1]]$death), est = r$est, se = r$se) }))
iecv_pool <- do.call(rbind, lapply(split(byk, byk$comparison), function(s) { s <- s[s$deaths >= MIN_EVENTS_POOL, ]
  o <- re_pool(s$est, s$se, label = "dC"); o$comparison <- s$comparison[1]; o }))
byk$lo <- byk$est - 1.96 * byk$se; byk$hi <- byk$est + 1.96 * byk$se

# ---- 3. Decision (rule applied to extension 1) ----
e1 <- delta[delta$comparison == "E1 vs base", ]; ie1 <- iecv_pool[iecv_pool$comparison == "E1 vs base", ]
dc_judge <- min(e1$dC_corrected, ie1$est)
verdict <- if (e1$LR_p < ALPHA && dc_judge >= DELTA_C_THRESHOLD) "Clinically relevant added value" else if (e1$LR_p < ALPHA) "Statistically significant, limited clinical relevance" else "No added value found"
decision <- data.frame(item = c("Likelihood ratio test P (extension 1 vs base)", "Bootstrap-corrected dC (95% CI)", "IECV pooled dC (95% CI)", "dC used for the decision (smaller of the two)", "Conclusion (decision rule)"),
  value = c(formatC(e1$LR_p, format = "g", digits = 3), fmt_ci(e1$dC_corrected, e1$dC_corr_lo, e1$dC_corr_hi), fmt_ci(ie1$est, ie1$lo, ie1$hi),
            sprintf("%.4f", dc_judge), verdict))
msg("  Q3 conclusion: ", verdict)

# ---- 4. Categorical NRI and difference in net benefit (stacked IECV predictions, averaged across imputations) ----
ord <- lapply(stk, function(s) s[order(s$row_id), ]); avg <- ord[[1]][, c("row_id", "cluster", "y")]
for (m in names(MODS)) avg[[m]] <- rowMeans(sapply(ord, `[[`, m))
cat3 <- function(p) cut(p, c(-Inf, RISK_CUTS, Inf), right = FALSE, labels = FALSE)
nri_fun <- function(y, p0, p1) { c0 <- cat3(p0); c1 <- cat3(p1); up <- c1 > c0; dn <- c1 < c0
  c(event_NRI = mean(up[y == 1]) - mean(dn[y == 1]), nonevent_NRI = mean(dn[y == 0]) - mean(up[y == 0])) }
set.seed(SEED)
nri <- do.call(rbind, lapply(c("E1", "E2", "E3"), function(e) { est <- nri_fun(avg$y, avg$base, avg[[e]])
  bt <- replicate(1000, { i <- sample.int(nrow(avg), replace = TRUE); nri_fun(avg$y[i], avg$base[i], avg[[e]][i]) })
  data.frame(comparison = paste(e, "vs base"), component = names(est), est = unname(est),
             lo = apply(bt, 1, stats::quantile, .025), hi = apply(bt, 1, stats::quantile, .975)) }))
nbf <- function(y, p, pt) mean(p >= pt & y == 1) - mean(p >= pt & y == 0) * pt / (1 - pt)
dnb <- do.call(rbind, lapply(c("E1", "E2", "E3"), function(e) { bt <- replicate(500, { i <- sample.int(nrow(avg), replace = TRUE)
  sapply(DCA_GRID, function(pt) nbf(avg$y[i], avg[[e]][i], pt) - nbf(avg$y[i], avg$base[i], pt)) })
  data.frame(comparison = paste(e, "vs base"), threshold = DCA_GRID,
             est = 100 * sapply(DCA_GRID, function(pt) nbf(avg$y, avg[[e]], pt) - nbf(avg$y, avg$base, pt)),
             lo = 100 * apply(bt, 1, stats::quantile, .025), hi = 100 * apply(bt, 1, stats::quantile, .975)) }))

# ---- 5. H1b: markers over ICU scores ----
msg("5. H1b")
h1b_one <- function(db, term, label) {
  xx <- lapply(xt, function(x) x[x$db == db & !is.na(x[[term]]), ])
  if (!nrow(xx[[1]]) || sum(xx[[1]]$death) < 20) return(NULL)
  t0_ <- term; t1_ <- c(term, MARK_E2); lr <- D3_test(xx, t1_, t0_)
  per <- t(sapply(seq_along(xx), function(i) { x <- xx[[i]]
    v0 <- validate_one(t0_, x, seed = SEED + 60 + i); v1 <- validate_one(t1_, x, seed = SEED + 60 + i)
    cs <- c_stats(x$death, cbind(stats::fitted(fit_glm(t0_, x)), stats::fitted(fit_glm(t1_, x))))
    c(c0 = unname(v0["c_corr"]), dc = unname(cs$dc - ((v1["c_app"] - v1["c_corr"]) - (v0["c_app"] - v0["c_corr"]))), var = cs$var_dc) }))
  r <- rubin(per[, "dc"], per[, "var"])
  data.frame(comparison = label, db = db, n = nrow(xx[[1]]), deaths = sum(xx[[1]]$death), df = lr$df1, LR_p = lr$p,
             C_score_corrected = mean(per[, "c0"]), dC_corrected = r$est, dC_lo = r$lo, dC_hi = r$hi,
             verdict = if (lr$p < ALPHA && r$est >= DELTA_C_THRESHOLD) "Clinically relevant added value" else if (lr$p < ALPHA) "Statistically significant, limited clinical relevance" else "No added value found")
}
h1b <- rbind(h1b_one("eICU", "apache_logit", "APACHE IVa + markers vs APACHE IVa"),
             h1b_one("MIMIC", "sapsii", "SAPS II + markers vs SAPS II"), h1b_one("MIMIC", "sofa", "SOFA + markers vs SOFA"))

# ---- 6. Q4: ratio constraints and index comparison ----
msg("6. Q4")
fit_all <- function(tm) lapply(xt, function(x) fit_glm(tm, x))
fN <- fit_all(c("lp", "ln_n", "ln_l", "alb")); fP <- fit_all(c("lp", "alb", "lymph"))
h2 <- rbind(cbind(constraint = "NLR: β[ln N] + β[ln L] = 0", lincomb(fN, c(ln_n = 1, ln_l = 1))),
            cbind(constraint = "PNI: β[L] − 5×β[albumin] = 0 (albumin g/L, L ×10⁹/L)", lincomb(fP, c(lymph = 1, alb = -5))))
h2$verdict <- ifelse(h2$p < ALPHA, "Constraint rejected (ratio/weighted form loses information)", "Constraint not rejected")
coef_N <- pool_coef(fN)$coef; coef_P <- pool_coef(fP)$coef
IDX <- list(NLR = list("ln_nlr", c("ln_n", "ln_l")), PLR = list("ln_plr", c("ln_plt", "ln_l")),
            SII = list("ln_sii", c("ln_plt", "ln_n", "ln_l")), SIRI = list("ln_siri", c("ln_n", "ln_mono", "ln_l")),
            PNI = list("pni", c("alb", "lymph")), GLR = list("ln_glr", c("ln_glu_mmol", "ln_l")), LAR = list("ln_lar", c("ln_lac", "ln_albd")))
s16 <- do.call(rbind, lapply(names(IDX), function(nm) { ix <- IDX[[nm]]; ti <- c("lp", ix[[1]]); tc <- c("lp", ix[[2]])
  per <- t(sapply(xt, function(x) { f0 <- fit_glm("lp", x); fi <- fit_glm(ti, x); fc <- fit_glm(tc, x)
    cs <- c_stats(x$death, cbind(stats::fitted(f0), stats::fitted(fi), stats::fitted(fc)))
    c(LR_index = f0$deviance - fi$deviance, LR_comps = f0$deviance - fc$deviance, AIC_index = stats::AIC(fi), AIC_comps = stats::AIC(fc),
      dC_index = cs$dc[1], dC_comps = cs$dc[2]) }))
  t_ic <- D3_test(xt, tc, ti)
  data.frame(index = nm, components = paste(ix[[2]], collapse = " + "), LR_chi2_index = stats::median(per[, "LR_index"]),
             LR_chi2_components = stats::median(per[, "LR_comps"]), AIC_index = mean(per[, "AIC_index"]), AIC_components = mean(per[, "AIC_comps"]),
             dC_index = mean(per[, "dC_index"]), dC_components = mean(per[, "dC_comps"]), df_constraint = t_ic$df1, p_constraint = t_ic$p,
             note = if (nm == "LAR") "Lactate often missing; relies on imputation" else "") }))

# ---- 7. Inverse probability weighting (whether measured; the testing model does not include the outcome) ----
msg("7. Inverse probability weighting")
wauc <- function(y, p, w) { pos <- which(y == 1); neg <- which(y == 0)
  sum(outer(w[pos], w[neg]) * (outer(p[pos], p[neg], ">") + 0.5 * outer(p[pos], p[neg], "=="))) / sum(outer(w[pos], w[neg])) }
ipw <- t(sapply(seq_along(d02$xs), function(j) { x <- d02$xs[[j]]
  det <- stats::glm(stats::reformulate(c(CORE_TERMS, "cluster"), "tested"), data = x, family = stats::binomial())
  ps <- stats::fitted(det)[x$tested == 1]; w <- mean(x$tested) / ps; q <- stats::quantile(w, c(.01, .99)); w <- pmin(pmax(w, q[1]), q[2])
  xx <- xt[[j]]; f0 <- fit_glm("lp", xx, w = w); f1 <- fit_glm(MODS$E1, xx, w = w)
  c(dC_weighted = wauc(xx$death, stats::fitted(f1), w) - wauc(xx$death, stats::fitted(f0), w), stats::coef(f1)[MARK_E1], w_min = min(w), w_max = max(w)) }))
ipw_tab <- data.frame(dC_weighted_apparent = mean(ipw[, "dC_weighted"]),
                      OR_ln_n = exp(mean(ipw[, "ln_n"])), OR_ln_l = exp(mean(ipw[, "ln_l"])), OR_alb = exp(mean(ipw[, "alb"])),
                      weight_range = sprintf("%.2f–%.2f", min(ipw[, "w_min"]), max(ipw[, "w_max"])))

# ---- 8. Dose-response (adjusted for LP) ----
msg("8. Dose-response curves")
MKC <- list(neut = TRUE, lymph = TRUE, alb = FALSE, eos = FALSE, nlr = TRUE, pni = FALSE)
xt2 <- lapply(xt, function(x) { x$nlr <- exp(x$ln_nlr); x })
dr <- lapply(names(MKC), function(v) { lg <- MKC[[v]]; tr <- if (lg) function(z) log(pmax(z, ZERO_FLOOR)) else identity
  t1 <- tr(xt2[[1]][[v]]); kn <- unique(stats::quantile(t1, c(.05, .35, .65, .95), na.rm = TRUE))
  if (length(kn) < 3) kn <- stats::quantile(t1[t1 > min(t1, na.rm = TRUE)], c(.1, .5, .9), na.rm = TRUE)
  rng <- stats::quantile(t1, c(.01, .99), na.rm = TRUE); grid <- seq(rng[1], rng[2], length.out = 120); ref <- stats::median(t1, na.rm = TRUE)
  basis <- function(z) { b <- Hmisc::rcspline.eval(z, knots = kn, inclx = TRUE); colnames(b) <- paste0("s", seq_len(ncol(b))); b }
  fits <- lapply(xt2, function(x) stats::glm(stats::reformulate(c("lp", colnames(basis(0))), "death"), data = cbind(x[, c("death", "lp")], basis(tr(x[[v]]))), family = stats::binomial()))
  sc <- colnames(basis(0)); Cg <- sweep(basis(grid), 2, basis(ref))
  est <- sapply(fits, function(f) drop(Cg %*% stats::coef(f)[sc])); vv <- sapply(fits, function(f) rowSums((Cg %*% stats::vcov(f)[sc, sc]) * Cg))
  pr <- do.call(rbind, lapply(seq_along(grid), function(g) rubin(est[g, ], vv[g, ])))
  xv <- xt2[[1]][[v]]; xv <- xv[!is.na(xv)]; xg <- if (lg) exp(grid) else grid; xv <- xv[xv >= min(xg) & xv <= max(xg)]
  dz <- stats::density(if (lg) log(pmax(xv, ZERO_FLOOR)) else xv, n = 256)
  list(var = v, log = lg, curve = data.frame(x = xg, or = exp(pr$est), lo = exp(pr$lo), hi = exp(pr$hi)),
       p_overall = D1_test(fits, sc), p_nonlinear = D1_test(fits, sc[-1]),
       dens = data.frame(x = if (lg) exp(dz$x) else dz$x, y = dz$y / max(dz$y))) })
names(dr) <- names(MKC)
grid_ref <- list(lp_median = stats::median(xt[[1]]$lp), q = lapply(c(neut = "neut", lymph = "lymph", alb = "alb"), function(v) stats::quantile(xt[[1]][[v]], c(.01, .5, .99))))

saveRDS(list(delta = delta, iecv_by_cluster = byk, iecv_pool = iecv_pool, decision = decision, verdict = verdict, nri = nri, dnb = dnb,
             h1b = h1b, h2 = h2, s16 = s16, ipw = ipw_tab, dr = dr, coef_N = coef_N, coef_P = coef_P, grid_ref = grid_ref, kmk = kmk,
             n = nrow(xt[[1]]), deaths = sum(xt[[1]]$death)), file.path(DIR_LOCAL, "05_markers.rds"))
write_share(list(
  README = data.frame(note = c("Q3/Q4 only in the subset with albumin and differential measured; base model = intercept + core model linear predictor LP",
    "E1 = + ln neutrophils, ln lymphocytes, albumin (primary); E2 = E1 + eosinophils < 0.05; E3 = splines with 3 knots for each of the three",
    "decision: decision rule (LR P < 0.05 and dC >= 0.02; dC = smaller of the bootstrap-corrected and IECV pooled estimates)",
    "In IECV, LP comes from the core model of each development fold; extended models are fitted by ML in the subset of the development fold (<= 8 parameters, not shrunk)",
    "H2: tests of the NLR and PNI linear constraints on top of the base model; S16: 7 indices compared with their components", "This file contains summary statistics only")),
  decision = decision, increment = delta, IECV_pooled = iecv_pool, IECV_by_cluster = byk, NRI = nri, H1b = h1b,
  H2_constraints = h2, index_comparison = s16, IPW = ipw_tab,
  dose_response_p = data.frame(marker = names(dr), p_overall = sapply(dr, `[[`, "p_overall"), p_nonlinear = sapply(dr, `[[`, "p_nonlinear")),
  delta_net_benefit = dnb), "05_markers.xlsx")
msg("05 done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min")
