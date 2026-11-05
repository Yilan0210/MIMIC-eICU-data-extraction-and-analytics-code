# =====================================================================
# 00_setup.R: shared settings and functions for modelling
# Analysis population: all eligible patients in both databases (pop == "main" and in_s2 == 1; 4,417 patients); marker subset: in_analysis == 1 (1,858 patients)
# Reads only: <DATA_DIR>/aecopd_clean_<TAG>.csv (or the clean_all sheet of _LOCAL_ONLY.xlsx)
# Outputs: <OUT_DIR>/local/ (patient-level, local only), result/ (summary tables), figure/ (figures), calculator/
# =====================================================================
# ======================== Settings ========================
DATA_DIR <- Sys.getenv("AECOPD_DATA_DIR", "/path/to/data")
TAG      <- Sys.getenv("AECOPD_TAG", "20261102")          # tag of the final data snapshot
OUT_DIR  <- Sys.getenv("AECOPD_OUT_DIR", file.path(DATA_DIR, "model_out"))

M_IMP     <- as.integer(Sys.getenv("AECOPD_M", "33"))      # number of imputations: m = 33 in the final run (>= percentage of incomplete cases, 32.1%)
MAXIT     <- as.integer(Sys.getenv("AECOPD_MAXIT", "20"))  # mice iterations
B_BOOT    <- as.integer(Sys.getenv("AECOPD_B", "200"))     # bootstrap samples (shrinkage, optimism correction, stability)
SEED      <- 20261103
N_CORES   <- as.integer(Sys.getenv("AECOPD_CORES", "1"))   # > 1: IECV clusters run in parallel (Linux)

# ---- Core model (16 variables, 24 degrees of freedom; all entered, no selection) ----
SPLINE_VARS <- c("age", "ln_bun", "hr", "rr", "map", "temp", "ph", "paco2")    # restricted cubic splines with 3 knots (10th/50th/90th percentiles)
LINEAR_VARS <- c("female", "gcs", "hgb", "ln_plt", "ln_wbc", "ln_glu", "af", "imv")
CORE_TERMS  <- c(as.vector(rbind(SPLINE_VARS, paste0(SPLINE_VARS, "_nl"))), LINEAR_VARS)
CORE_VARS   <- c(SPLINE_VARS, LINEAR_VARS)                                       # 16 variables
KNOT_P      <- c(0.10, 0.50, 0.90)
APPROX_R2   <- 0.95       # simplified model: lower bound of R2 for approximating the full model's linear predictor
VAR_LAB <- c(age = "Age", ln_bun = "Blood urea nitrogen", hr = "Heart rate", rr = "Respiratory rate",
             map = "Mean arterial pressure", temp = "Temperature", ph = "pH", paco2 = "PaCO\u2082",
             female = "Female sex", gcs = "GCS motor score", hgb = "Haemoglobin", ln_plt = "Platelets",
             ln_wbc = "White cell count", ln_glu = "Glucose", af = "Atrial fibrillation", imv = "Invasive ventilation")
VAR_UNIT <- c(age = "years", ln_bun = "mg/dL", hr = "beats/min", rr = "breaths/min", map = "mmHg", temp = "\u00b0C",
              ph = "", paco2 = "mmHg", gcs = "points", hgb = "g/dL", ln_plt = "\u00d710\u2079/L", ln_wbc = "\u00d710\u2079/L",
              ln_glu = "mg/dL")
# variables on the original scale (log terms back-transformed in figures and tables)
RAW_OF <- c(age = "age", ln_bun = "bun", hr = "hr", rr = "rr", map = "map", temp = "temp", ph = "ph", paco2 = "paco2",
            gcs = "gcs", hgb = "hgb", ln_plt = "plt", ln_wbc = "wbc", ln_glu = "glucose")

# ---- Markers (Q3/Q4; only in the subset with in_analysis == 1) ----
MARK_E1   <- c("ln_n", "ln_l", "alb")          # extension 1: 3 df (primary; assessed with the decision rule)
MARK_E2   <- c(MARK_E1, "eos_low")              # extension 2
MARK_SPL  <- c("ln_n", "ln_n_nl", "ln_l", "ln_l_nl", "alb", "alb_nl")   # extension 3: splines
ZERO_FLOOR <- 0.01         # value substituted for N, L or monocytes = 0 before log transformation
EOS_CUT    <- 0.05         # eosinopenia
DELTA_C_THRESHOLD <- 0.02  # Q3 decision threshold
ALPHA      <- 0.05
NB_THRESHOLDS  <- c(0.10, 0.20, 0.30)       # thresholds for random-effects pooling of net benefit
DCA_GRID       <- seq(0.05, 0.50, by = 0.01) # decision curves
RISK_CUTS      <- c(0.10, 0.30)              # risk groups <10%, 10–30%, >=30%
MIN_EVENTS_POOL <- 10                         # clusters with < 10 deaths are not pooled
# ===============================================================

# Console encoding: switch to a UTF-8 locale if possible (does not affect the result files)
for (lc in c("C.UTF-8", "en_US.UTF-8", "zh_CN.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", lc)))) break

suppressPackageStartupMessages({
  for (p in c("mice", "rms", "metafor", "glmnet", "readxl", "writexl")) {
    if (!requireNamespace(p, quietly = TRUE))
      stop("Missing R package ", p, ". Install first: install.packages(c('mice','rms','metafor','glmnet','readxl','writexl'))")
  }
  library(mice); library(rms)
})
options(warn = 1, datadist = NULL)

DIR_LOCAL <- file.path(OUT_DIR, "local")
DIR_RESULT <- file.path(OUT_DIR, "result")
DIR_FIG    <- file.path(OUT_DIR, "figure")
dir.create(DIR_LOCAL, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "", ..., "\n")

# ---------------------------------------------------------------------
# 1. Read the data
# ---------------------------------------------------------------------
load_clean <- function() {
  csv  <- file.path(DATA_DIR, sprintf("aecopd_clean_%s.csv", TAG))
  xlsx <- file.path(DATA_DIR, sprintf("aecopd_clean_%s_LOCAL_ONLY.xlsx", TAG))
  if (file.exists(csv)) {
    d <- utils::read.csv(csv, stringsAsFactors = FALSE, na.strings = c("", "NA"))
    src <- csv
  } else if (file.exists(xlsx)) {
    d <- as.data.frame(readxl::read_excel(xlsx, sheet = "clean_all", guess_max = 100000))
    src <- xlsx
  } else stop("In ", DATA_DIR, ": cannot find aecopd_clean_", TAG, ".csv or _LOCAL_ONLY.xlsx")
  msg("Read", src, ": ", nrow(d), "rows")
  d
}

# Study population and analysis variables (before imputation, original scale)
IMP_VARS <- c(
  # outcome and predictors
  "death", "age", "bun", "gcs", "hr", "ph", "af", "imv", "alb", "neut", "lymph", "eos",
  # auxiliary variables (imputation only, and sensitivity analyses / comparator scores)
  "sex", "cluster", "gcs_valid", "paco2", "lactate", "creatinine", "hgb", "plt", "wbc",
  "spo2", "rr", "map", "temp", "steroid", "pneumonia24", "glucose", "mono", "vaso")
CARRY_VARS <- c("row_id", "db", "pop", "apache", "sapsii", "sofa", "af_any", "death_s3",
                "limit_s18", "race4", "steroid_obs", "in_s2", "in_analysis")

prep_vars <- function(d) {
  num <- function(x) suppressWarnings(as.numeric(x))
  out <- data.frame(
    row_id = seq_len(nrow(d)), db = as.character(d$db), pop = as.character(d$pop),
    death = num(d$death), age = num(d$age), bun = num(d$bun), gcs = num(d$gcs_motor_first_d1),
    hr = num(d$hr_first_d1), ph = num(d$ph), af = num(d$af_24h), imv = num(d$invasive_vent_d1),
    alb = num(d$alb_gl), neut = num(d$neut_abs), lymph = num(d$lymph_abs), eos = num(d$eos_abs),
    sex = factor(d$sex, levels = c("M", "F")), cluster = factor(num(d$cluster)),
    gcs_valid = num(d$gcs_valid_first_d1), paco2 = num(d$paco2), lactate = num(d$lactate),
    creatinine = num(d$creatinine), hgb = num(d$hgb), plt = num(d$plt), wbc = num(d$wbc),
    spo2 = num(d$spo2_first_d1), rr = num(d$rr_first_d1), map = num(d$map_first_d1),
    temp = num(d$temp_first_d1), steroid = num(d$steroid_before_diff),
    pneumonia24 = num(d$pneumonia_24h), glucose = num(d$glucose), mono = num(d$mono_abs),
    vaso = num(d$vasopressor_d1),
    apache = num(d$apache_iva_prob), sapsii = num(d$sapsii), sofa = num(d$sofa_d1),
    af_any = num(d$af_any), death_s3 = num(d$death_s3), limit_s18 = num(d$limit_tx_s18),
    race4 = as.character(d$race4), steroid_obs = num(d$steroid_before_diff),
    in_s2 = num(d$in_s2), in_analysis = num(d$in_analysis),
    stringsAsFactors = FALSE)
  # eICU has no hospice category: when death_s3 is empty it equals in-hospital death
  out$death_s3 <- ifelse(is.na(out$death_s3), out$death, out$death_s3)
  # MIMIC pneumonia_24h comes from MIMIC-IV-ED diagnoses (sql/mimic/03_ed_steroid.sql);
  # missing = no linked emergency department record (direct or transfer admission, or a year not covered by the ED module), i.e. "unknown" rather than "no pneumonia", so it stays missing and is imputed
  out
}

# Analysis population: all eligible patients; tested = albumin and differential measured (Q3/Q4 subset)
analysis_pop <- function(v, pop = "main") {
  a <- v[v$pop %in% pop & v$in_s2 %in% 1, ]
  a$tested <- as.numeric(a$in_analysis %in% 1)
  a$cluster <- droplevels(a$cluster)
  a
}

# ---------------------------------------------------------------------
# 2. Imputation and derived variables
# ---------------------------------------------------------------------
# with_outcome = FALSE: imputation of an IECV held-out cluster (without the outcome)
run_mice <- function(a, m = M_IMP, with_outcome = TRUE, seed = SEED) {
  dat <- a[, intersect(IMP_VARS, names(a)), drop = FALSE]
  if (!with_outcome) dat$death <- NULL
  # drop constant columns (e.g. cluster within a single cluster)
  keep <- vapply(dat, function(x) length(unique(x[!is.na(x)])) > 1, logical(1))
  dropped <- dat[, !keep, drop = FALSE]
  dat <- dat[, keep, drop = FALSE]
  meth <- mice::make.method(dat)
  meth[meth != ""] <- "pmm"                    # predictive mean matching for numeric and binary variables (imputes observed values)
  for (f in names(dat)[vapply(dat, is.factor, logical(1))]) if (meth[f] != "") meth[f] <- "polyreg"
  imp <- mice::mice(dat, m = m, maxit = MAXIT, method = meth, seed = seed, printFlag = FALSE)
  attr(imp, "dropped") <- dropped
  attr(imp, "carry") <- a[, setdiff(names(a), names(dat)), drop = FALSE]
  imp
}

# Take imputed dataset i, add back the non-imputed columns and create the derived model variables
complete_i <- function(imp, i) {
  x <- mice::complete(imp, i)
  dr <- attr(imp, "dropped"); if (ncol(dr)) x <- cbind(x, dr)
  ca <- attr(imp, "carry");   x <- cbind(x, ca[, setdiff(names(ca), names(x)), drop = FALSE])
  derive(x)
}

floor0 <- function(x) ifelse(!is.na(x) & x <= 0, ZERO_FLOOR, x)
derive <- function(x) {
  x$ln_bun  <- log(pmax(x$bun, 1))
  x$ln_n    <- log(floor0(x$neut))
  x$ln_l    <- log(floor0(x$lymph))
  x$eos_low <- as.numeric(x$eos < EOS_CUT)
  x$eos_zero <- as.numeric(x$eos == 0)
  # S16 indices and components
  x$ln_plt  <- log(pmax(x$plt, 1)); x$ln_mono <- log(floor0(x$mono))
  x$ln_glu  <- log(pmax(x$glucose / 18.016, 0.1)); x$ln_lac <- log(pmax(x$lactate, 0.1))
  x$ln_albd <- log(x$alb / 10)
  x$ln_nlr  <- x$ln_n - x$ln_l
  x$ln_plr  <- x$ln_plt - x$ln_l
  x$ln_sii  <- x$ln_plt + x$ln_n - x$ln_l
  x$ln_siri <- x$ln_n + x$ln_mono - x$ln_l
  x$pni     <- x$alb + 5 * x$lymph
  x$ln_glr  <- x$ln_glu - x$ln_l
  x$ln_lar  <- x$ln_lac - x$ln_albd
  # comparator scores (computed on the imputed data)
  if ("gcs_valid" %in% names(x)) {
    nb <- (x$bun >= 25) + (x$gcs_valid < 14) + (x$hr >= 109)
    x$bap65 <- ifelse(nb == 0, ifelse(x$age < 65, 1, 2), nb + 2)
  }
  if ("pneumonia24" %in% names(x)) x$decaf <- x$eos_low + x$pneumonia24 + (x$ph < 7.30) + x$af
  x$db_mimic <- as.numeric(x$db == "MIMIC")
  x$female  <- as.numeric(x$sex == "F")
  x$ln_plt  <- log(pmax(x$plt, 1)); x$ln_wbc <- log(pmax(x$wbc, 0.1))
  x$ln_glu  <- log(pmax(x$glucose, 10))       # core model: log of glucose in mg/dL (the ln_glu used for GLR was converted above; see next line)
  x$ln_glu_mmol <- log(pmax(x$glucose / 18.016, 0.1))
  x$ln_glr  <- x$ln_glu_mmol - x$ln_l
  if ("apache" %in% names(x)) x$apache_logit <- qlogis(pmin(pmax(x$apache, 1e-4), 1 - 1e-4))
  x
}

# ---------------------------------------------------------------------
# 3. Models, pooling and performance measures
# ---------------------------------------------------------------------
fml <- function(terms, y = "death") stats::as.formula(paste(y, "~", paste(terms, collapse = " + ")))
# ---- Splines: non-linear term of a restricted cubic spline (3 knots); identical to Hmisc::rcspline.eval(norm = 2) and rms::rcs ----
rcs_nl <- function(x, k) { p <- function(z) pmax(z, 0)^3
  (p(x - k[1]) - p(x - k[2]) * (k[3] - k[1]) / (k[3] - k[2]) + p(x - k[3]) * (k[2] - k[1]) / (k[3] - k[2])) / (k[3] - k[1])^2 }
# Knots: 10th/50th/90th percentiles of the observed (pre-imputation) values in the development data; determined separately in each IECV development fold
make_knots <- function(a, vars = SPLINE_VARS) {
  d <- derive(a); setNames(lapply(vars, function(v) unname(stats::quantile(d[[v]], KNOT_P, na.rm = TRUE))), vars)
}
add_design <- function(x, knots) { for (v in names(knots)) x[[paste0(v, "_nl")]] <- rcs_nl(x[[v]], knots[[v]]); x }

fit_glm <- function(terms, x, y = "death", w = NULL) {
  if (is.null(w)) stats::glm(fml(terms, y), data = x, family = stats::binomial())
  else { x$.w <- w; suppressWarnings(stats::glm(fml(terms, y), data = x, family = stats::quasibinomial(), weights = .w)) }
}
loglik_at <- function(beta, terms, x, y = "death") {
  X <- cbind(1, as.matrix(x[, terms, drop = FALSE])); p <- plogis(drop(X %*% beta))
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12); yy <- x[[y]]
  sum(yy * log(p) + (1 - yy) * log(1 - p))
}

# Rubin's rules: scalar
rubin <- function(est, var) {
  m <- length(est); qb <- mean(est); ub <- mean(var); b <- stats::var(est)
  if (m < 2 || !is.finite(b)) b <- 0
  tv <- ub + (1 + 1 / m) * b
  df <- if (b > 0) (m - 1) * (1 + ub / ((1 + 1 / m) * b))^2 else Inf
  q <- if (is.finite(df)) stats::qt(0.975, df) else stats::qnorm(0.975)
  data.frame(est = qb, se = sqrt(tv), lo = qb - q * sqrt(tv), hi = qb + q * sqrt(tv),
             p = 2 * (if (is.finite(df)) stats::pt(-abs(qb / sqrt(tv)), df) else stats::pnorm(-abs(qb / sqrt(tv)))),
             df = df)
}
# Rubin's rules: coefficients and covariance of a set of glm fits
pool_coef <- function(fits) {
  B <- sapply(fits, stats::coef); m <- ncol(B)
  qb <- rowMeans(B); Ub <- Reduce(`+`, lapply(fits, stats::vcov)) / m
  Bm <- if (m > 1) stats::cov(t(B)) else Ub * 0
  list(coef = qb, vcov = Ub + (1 + 1 / m) * Bm, Ubar = Ub, B = Bm, m = m)
}
lincomb <- function(fits, L) {
  est <- vapply(fits, function(f) sum(L * stats::coef(f)[names(L)]), numeric(1))
  var <- vapply(fits, function(f) drop(t(L) %*% stats::vcov(f)[names(L), names(L)] %*% L), numeric(1))
  rubin(est, var)
}

# Meng-Rubin likelihood ratio test (D3) for nested logistic models
D3_test <- function(xs, terms1, terms0, y = "death") {
  m <- length(xs); k <- length(terms1) - length(terms0)
  f1 <- lapply(xs, function(x) fit_glm(terms1, x, y)); f0 <- lapply(xs, function(x) fit_glm(terms0, x, y))
  b1 <- rowMeans(sapply(f1, stats::coef)); b0 <- rowMeans(sapply(f0, stats::coef))
  dL <- mapply(function(a, b) 2 * (as.numeric(stats::logLik(a)) - as.numeric(stats::logLik(b))), f1, f0)
  dT <- vapply(xs, function(x) 2 * (loglik_at(b1, terms1, x, y) - loglik_at(b0, terms0, x, y)), numeric(1))
  rL <- max(0, (m + 1) / (k * (m - 1)) * (mean(dL) - mean(dT)))
  DL <- mean(dT) / (k * (1 + rL)); t <- k * (m - 1)
  v <- if (rL == 0) Inf else if (t > 4) 4 + (t - 4) * (1 + (1 - 2 / t) / rL)^2 else t * (1 + 1 / k) * (1 + 1 / rL)^2 / 2
  p <- if (is.finite(v)) stats::pf(DL, k, v, lower.tail = FALSE) else stats::pchisq(DL * k, k, lower.tail = FALSE)
  newinfo <- stats::median(vapply(seq_len(m), function(i) {
    lr1 <- f1[[i]]$null.deviance - f1[[i]]$deviance; lr0 <- f0[[i]]$null.deviance - f0[[i]]$deviance
    1 - lr0 / lr1 }, numeric(1)))
  data.frame(df1 = k, df2 = v, F = DL, p = p, rL = rL, median_LR_chi2 = stats::median(dL),
             new_info_fraction = newinfo)
}

# DeLong: C statistics of several predictions in the same patients and their covariance
delong <- function(y, P) {
  P <- as.matrix(P); pos <- which(y == 1); neg <- which(y == 0)
  n1 <- length(pos); n0 <- length(neg); K <- ncol(P)
  V10 <- matrix(NA, n1, K); V01 <- matrix(NA, n0, K)
  for (k in seq_len(K)) {
    a <- P[pos, k]; b <- P[neg, k]
    psi <- outer(a, b, ">") + 0.5 * outer(a, b, "==")
    V10[, k] <- rowMeans(psi); V01[, k] <- colMeans(psi)
  }
  auc <- colMeans(V10)
  S <- stats::cov(V10) / n1 + stats::cov(V01) / n0
  list(auc = auc, S = S)
}
# C (variance on the logit scale) and paired dC
c_stats <- function(y, P) {
  dl <- delong(y, P); auc <- dl$auc; S <- as.matrix(dl$S)
  out <- list(c = auc, var_c = diag(S), logit_c = qlogis(auc), var_logit_c = diag(S) / (auc * (1 - auc))^2)
  if (ncol(as.matrix(P)) >= 2) {
    out$dc <- auc[-1] - auc[1]
    out$var_dc <- vapply(2:length(auc), function(k) S[k, k] + S[1, 1] - 2 * S[1, k], numeric(1))
  }
  out
}

# Calibration: O/E, calibration-in-the-large (CITL) and slope
calib <- function(y, p) {
  lp <- qlogis(pmin(pmax(p, 1e-8), 1 - 1e-8)); O <- sum(y); E <- sum(p)
  f_s <- stats::glm(y ~ lp, family = stats::binomial())
  f_i <- stats::glm(y ~ 1 + offset(lp), family = stats::binomial())
  c(log_oe = log(O / E), var_log_oe = (1 - O / length(y)) / O,
    slope = unname(stats::coef(f_s)[2]), var_slope = stats::vcov(f_s)[2, 2],
    citl = unname(stats::coef(f_i)[1]), var_citl = stats::vcov(f_i)[1, 1])
}
# Flexible calibration curve measures (Austin & Steyerberg 2019)
ici <- function(y, p) {
  fit <- stats::loess(y ~ p, span = 0.75, degree = 2)
  pc <- stats::predict(fit, newdata = data.frame(p = p)); d <- abs(pc - p)
  c(ICI = mean(d, na.rm = TRUE), E50 = stats::median(d, na.rm = TRUE), E90 = unname(stats::quantile(d, 0.9, na.rm = TRUE)))
}
calib_curve <- function(y, p, grid = seq(0.01, 0.8, by = 0.01)) {
  fit <- stats::loess(y ~ p, span = 0.75, degree = 2)
  rg <- stats::quantile(p, c(0.01, 0.99))           # plot only within the 1st–99th percentiles of predictions (too few data at the extremes)
  g <- grid[grid >= rg[1] & grid <= rg[2]]
  data.frame(pred = g, obs = stats::predict(fit, newdata = data.frame(p = g)))
}
# Net benefit and its variance (multinomial approximation)
net_benefit <- function(y, p, pt) {
  n <- length(y); w <- pt / (1 - pt); tp <- sum(p >= pt & y == 1) / n; fp <- sum(p >= pt & y == 0) / n
  nb <- tp - w * fp
  v <- (tp * (1 - tp) + w^2 * fp * (1 - fp) + 2 * w * tp * fp) / n
  c(nb = nb, var = v)
}
nb_treat_all <- function(y, pt) { pr <- mean(y); pr - (1 - pr) * pt / (1 - pt) }

# Random-effects pooling: REML + Hartung-Knapp-Sidik-Jonkman; 95% prediction interval with t(k-2)
re_pool <- function(yi, sei, back = identity, label = "") {
  ok <- is.finite(yi) & is.finite(sei) & sei > 0; yi <- yi[ok]; sei <- sei[ok]; k <- length(yi)
  if (k < 2) return(data.frame(measure = label, k = k, est = NA, lo = NA, hi = NA, pi_lo = NA, pi_hi = NA, tau2 = NA, method = NA))
  # If REML does not converge (common with few clusters and tau^2 near 0), REML with a larger step and then Paule-Mandel are tried; the method used is stored in column method
  fit <- function(m, ctrl = list()) tryCatch(metafor::rma(yi = yi, sei = sei, method = m, test = "knha", control = ctrl), error = function(e) NULL)
  # With only 2 clusters, random-effects (HKSJ) intervals are not meaningful (t with 1 df), so fixed-effect pooling is used
  if (k < 3) { r <- metafor::rma(yi = yi, sei = sei, method = "FE"); meth <- "fixed effect (k = 2)" } else { r <- fit("REML"); meth <- "REML" }
  if (is.null(r)) { r <- fit("REML", list(stepadj = 0.5, maxiter = 1000)); meth <- "REML (step 0.5)" }
  if (is.null(r)) { r <- fit("PM"); meth <- "Paule-Mandel" }
  mu <- as.numeric(r$beta); se <- r$se
  if (k >= 3) { q <- stats::qt(0.975, k - 2); pil <- mu - q * sqrt(r$tau2 + se^2); pih <- mu + q * sqrt(r$tau2 + se^2) }
  else { pil <- NA; pih <- NA }
  data.frame(measure = label, k = k, est = back(mu), lo = back(r$ci.lb), hi = back(r$ci.ub),
             pi_lo = back(pil), pi_hi = back(pih), tau2 = r$tau2, method = meth)
}

# ---------------------------------------------------------------------
# 4. Modelling: maximum likelihood + uniform bootstrap shrinkage, with optimism-corrected performance
# ---------------------------------------------------------------------
validate_one <- function(terms, x, B = B_BOOT, seed = SEED, y = "death") {
  # The outcome is given by argument y (so that alternative outcomes, e.g. A9/A10, are evaluated correctly)
  f <- rms::lrm(fml(terms, y), data = x, x = TRUE, y = TRUE)
  set.seed(seed)
  v <- rms::validate(f, B = B, method = "boot")
  c(c_app = (v["Dxy", "index.orig"] + 1) / 2, c_corr = (v["Dxy", "index.corrected"] + 1) / 2,
    shrink = v["Slope", "index.corrected"], brier_app = v["B", "index.orig"], brier_corr = v["B", "index.corrected"])
}

# Fit one model on a set of imputed datasets: returns pooled shrunken coefficients, re-estimated intercept and validation results per imputation
develop <- function(xs, terms, B = B_BOOT, method = c("ml_shrink", "ridge"), seed = SEED) {
  method <- match.arg(method)
  if (method == "ridge") {
    cf <- sapply(seq_along(xs), function(i) {
      X <- as.matrix(xs[[i]][, terms]); set.seed(seed + i)
      cv <- glmnet::cv.glmnet(X, xs[[i]]$death, family = "binomial", alpha = 0, nfolds = 10)
      as.numeric(stats::coef(cv, s = "lambda.min")) })
    rownames(cf) <- c("(Intercept)", terms)
    return(list(terms = terms, coef = rowMeans(cf), shrink = NA, val = NULL, method = "ridge"))
  }
  fits <- lapply(xs, function(x) fit_glm(terms, x))
  pc <- pool_coef(fits)
  val <- t(sapply(seq_along(xs), function(i) validate_one(terms, xs[[i]], B = B, seed = seed + i)))
  s <- mean(val[, "shrink"])
  b <- pc$coef[-1] * s
  int <- mean(vapply(xs, function(x) {
    off <- drop(as.matrix(x[, terms]) %*% b)
    stats::coef(stats::glm(x$death ~ 1 + offset(off), family = stats::binomial()))[1] }, numeric(1)))
  list(terms = terms, coef = c("(Intercept)" = unname(int), b), coef_ml = pc$coef, vcov_ml = pc$vcov,
       shrink = s, val = val, fits = fits, method = "ml_shrink")
}
predict_dev <- function(model, x) plogis(drop(cbind(1, as.matrix(x[, model$terms])) %*% model$coef))

# Excel output (summary tables only)
write_share <- function(sheets, file) {
  names(sheets) <- substr(names(sheets), 1, 31)
  sheets <- lapply(sheets, function(x) { x <- as.data.frame(x); rownames(x) <- NULL; x })
  writexl::write_xlsx(sheets, file.path(DIR_RESULT, file))
  msg("Written", file.path(DIR_RESULT, file))
}
fmt_ci <- function(e, l, h, d = 3) sprintf(paste0("%.", d, "f (%.", d, "f, %.", d, "f)"), e, l, h)

# ---------------------------------------------------------------------
# 5. Design matrix of the core model, simplified model, effects and importance
# ---------------------------------------------------------------------
cols_of <- function(v) if (v %in% SPLINE_VARS) c(v, paste0(v, "_nl")) else v
# Prepare a set of imputed datasets: derived variables + spline columns
prep_xs <- function(imp, knots) lapply(seq_len(imp$m), function(i) add_design(complete_i(imp, i), knots))

# Simplified model (Harrell's model approximation): least squares with the (shrunken) linear predictor of the full model as the outcome;
# variables (both spline columns together) are removed one at a time, each time the one whose removal lowers R2 least, until R2 would fall below APPROX_R2; variables kept in more than half of the imputations are retained
approx_model <- function(xs, model, vars = CORE_VARS) {
  lp_of <- function(x) qlogis(predict_dev(model, x))
  sel <- lapply(xs, function(x) {
    lp <- lp_of(x); keep <- vars
    r2 <- function(vv) { X <- as.matrix(x[, unlist(lapply(vv, cols_of)), drop = FALSE]); summary(stats::lm(lp ~ X))$r.squared }
    repeat {
      if (length(keep) <= 1) break
      cand <- vapply(keep, function(v) r2(setdiff(keep, v)), numeric(1))
      best <- names(which.max(cand)); if (cand[best] < APPROX_R2) break
      keep <- setdiff(keep, best)
    }
    keep })
  freq <- table(factor(unlist(sel), levels = vars))
  keep <- vars[as.numeric(freq) > length(xs) / 2]
  terms <- unlist(lapply(keep, cols_of))
  cf <- rowMeans(sapply(xs, function(x) { X <- as.matrix(x[, terms, drop = FALSE]); b <- stats::coef(stats::lm(lp_of(x) ~ X)); b }))
  names(cf) <- c("(Intercept)", terms)
  r2 <- mean(vapply(xs, function(x) stats::cor(lp_of(x), drop(cbind(1, as.matrix(x[, terms])) %*% cf))^2, numeric(1)))
  list(terms = terms, vars = keep, coef = cf, r2 = r2, freq = as.numeric(freq) / length(xs), method = "approx")
}

# Variable importance: pooled Wald chi-square (both spline terms jointly) minus degrees of freedom
importance <- function(model, vars = CORE_VARS) {
  b <- model$coef_ml; V <- model$vcov_ml
  do.call(rbind, lapply(vars, function(v) { k <- cols_of(v); bb <- b[k]; VV <- V[k, k, drop = FALSE]
    chi <- drop(t(bb) %*% solve(VV) %*% bb)
    data.frame(var = v, df = length(k), chi2 = chi, p = stats::pchisq(chi, length(k), lower.tail = FALSE), chi2_minus_df = chi - length(k)) }))
}

# Partial effect of one variable (log OR vs the reference value with 95% CI): pooled ML coefficients and Rubin covariance
partial_effect <- function(model, v, grid, ref, knots) {
  k <- cols_of(v); b <- model$coef_ml[k]; V <- model$vcov_ml[k, k, drop = FALSE]
  B <- function(z) if (v %in% SPLINE_VARS) cbind(z, rcs_nl(z, knots[[v]])) else cbind(z)
  L <- sweep(B(grid), 2, B(ref)); est <- drop(L %*% b); se <- sqrt(rowSums((L %*% V) * L))
  data.frame(x = grid, logor = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
}
# OR for the interquartile-range contrast (continuous) or 1 vs 0 (binary)
or_contrast <- function(model, v, x1, knots) {
  z <- x1[[v]]
  if (v %in% c("female", "af", "imv")) { lo <- 0; hi <- 1 } else { q <- stats::quantile(z, c(.25, .75), na.rm = TRUE); lo <- q[1]; hi <- q[2]; if (hi == lo) { lo <- hi - 1 } }   # GCS motor score quartiles are often equal: use a 1-point decrease instead
  pe <- partial_effect(model, v, hi, lo, knots)
  data.frame(var = v, from = unname(lo), to = unname(hi), OR = exp(pe$logor), lo = exp(pe$lo), hi = exp(pe$hi))
}
# Back-transform log-scale variables to the original scale (figures, tables)
to_raw <- function(v, z) if (v %in% c("ln_bun", "ln_plt", "ln_wbc", "ln_glu")) exp(z) else z
