# =====================================================================
# 01_describe.R: description and checks before modelling (no predictor-outcome model is fitted)
#   * patients and deaths in all eligible patients and in the marker subset; patients and deaths per cluster
#   * missingness of core candidate and auxiliary variables (by database)
#   * baseline comparison of patients with/without blood gas and with/without markers (outcome not examined)
#   * Riley sample size criteria (core model 24 df; two approaches for Q3)
# Outputs: result/01_describe.xlsx; local/01_population.rds
# =====================================================================
source("00_setup.R")
t0 <- Sys.time()
raw <- load_clean(); v <- prep_vars(raw); a <- analysis_pop(v)
saveRDS(list(a = a), file.path(DIR_LOCAL, "01_population.rds"))

coh <- do.call(rbind, lapply(list("All eligible patients" = a, "Marker subset" = a[a$tested == 1, ]), function(s)
  do.call(rbind, lapply(split(s, s$db), function(z) data.frame(db = z$db[1], n = nrow(z), deaths = sum(z$death), mortality_pct = round(100 * mean(z$death), 1))))))
coh$population <- sub("\\..*$", "", rownames(coh)); rownames(coh) <- NULL
msg("All eligible patients: ", nrow(a), "patients, deaths", sum(a$death), "; marker subset: ", sum(a$tested), "patients")
expect <- c(MIMIC = 1276, eICU = 3141)
got <- tapply(a$db, a$db, length)
if (any(got[names(expect)] != expect)) warning("Number of eligible patients differs from the QC table (MIMIC 1276, eICU 3141); check the data first!")

CL_LABEL <- c("1" = "MIMIC ICD-9 era", "2" = "MIMIC ICD-10 era", "3" = "eICU Midwest", "4" = "eICU Northeast",
              "5" = "eICU South", "6" = "eICU West", "7" = "eICU region missing")
cl <- do.call(rbind, lapply(split(a, a$cluster), function(z) data.frame(cluster = z$cluster[1], label = CL_LABEL[as.character(z$cluster[1])],
  n = nrow(z), deaths = sum(z$death), n_tested = sum(z$tested), deaths_tested = sum(z$death[z$tested == 1]),
  pooled_in_meta = ifelse(sum(z$death) >= MIN_EVENTS_POOL, "yes", "no (<10 deaths)"))))

RAWV <- c("age", "sex", "bun", "gcs", "hr", "rr", "map", "temp", "ph", "paco2", "hgb", "plt", "wbc", "glucose", "af", "imv",
          "lactate", "creatinine", "spo2", "gcs_valid", "vaso", "pneumonia24", "steroid", "alb", "neut", "lymph", "eos",
          "apache", "sapsii", "sofa")
role <- function(z) ifelse(z %in% c("age", "sex", "bun", "gcs", "hr", "rr", "map", "temp", "ph", "paco2", "hgb", "plt", "wbc", "glucose", "af", "imv"),
                           "Core predictor", ifelse(z %in% c("apache", "sapsii", "sofa"), "Comparator score", "Auxiliary variable (imputation/sensitivity)"))
miss <- do.call(rbind, lapply(split(a, a$db), function(z) data.frame(db = z$db[1], var = RAWV, role = role(RAWV),
  n = nrow(z), n_missing = vapply(RAWV, function(k) sum(is.na(z[[k]])), numeric(1)),
  pct_missing = round(100 * vapply(RAWV, function(k) mean(is.na(z[[k]])), numeric(1)), 1))))
inc <- a[, c("age", "sex", "bun", "gcs", "hr", "rr", "map", "temp", "ph", "paco2", "hgb", "plt", "wbc", "glucose", "af", "imv")]
pct_incomplete <- round(100 * mean(!stats::complete.cases(inc)), 1)
msg("Proportion of incomplete cases for core variables: ", pct_incomplete, "% (m = ", M_IMP, ")")
if (pct_incomplete > M_IMP) warning("The proportion of incomplete cases exceeds the number of imputations m; increase AECOPD_M to ", ceiling(pct_incomplete), " and document the change")

cmp_groups <- function(g, lab) {      # baseline characteristics only, no outcome
  vv <- c("age", "bun", "gcs", "hr", "rr", "map", "temp", "spo2", "hgb", "wbc", "glucose", "imv", "af", "vaso", "pneumonia24")
  do.call(rbind, lapply(vv, function(k) { z <- a[[k]]; f <- function(s) { s <- s[!is.na(s)]
    if (all(s %in% c(0, 1))) sprintf("%.1f%%", 100 * mean(s)) else sprintf("%.1f (%.1f–%.1f)", stats::median(s), stats::quantile(s, .25), stats::quantile(s, .75)) }
    data.frame(comparison = lab, variable = k, group_yes = f(z[g]), group_no = f(z[!g])) }))
}
abg <- !is.na(a$ph); cmp <- rbind(cmp_groups(abg, "Blood gas vs no blood gas"), cmp_groups(a$tested == 1, "Markers measured vs not measured"))
abg_n <- data.frame(db = c("MIMIC", "eICU"), pct_with_abg = round(100 * tapply(abg, a$db, mean)[c("MIMIC", "eICU")], 1),
                    pct_tested = round(100 * tapply(a$tested, a$db, mean)[c("MIMIC", "eICU")], 1))

# Riley criteria (C -> R2 conversion as in pmsampsize: 1 million simulated patients)
riley <- function(p, auc, prev, S = 0.9, delta = 0.05) {
  set.seed(1); n <- 1e6; mu <- sqrt(2) * stats::qnorm(auc); Y <- stats::rbinom(n, 1, prev)
  LP <- ifelse(Y == 1, stats::rnorm(n, mu, 1), stats::rnorm(n, 0, 1))
  f <- stats::glm(Y ~ LP, family = stats::binomial()); R <- 1 - exp(-(f$null.deviance - f$deviance) / n)
  maxR <- 1 - (prev^prev * (1 - prev)^(1 - prev))^2
  n1 <- p / ((S - 1) * log(1 - R / S)); Sd <- R / (R + delta * maxR); n2 <- p / ((Sd - 1) * log(1 - R / Sd))
  n3 <- (1.96 / 0.05)^2 * prev * (1 - prev); nn <- ceiling(max(n1, n2, n3))
  data.frame(parameters = p, assumed_C = auc, prevalence = round(prev, 4), R2_CS = round(R, 4), n_crit1 = ceiling(n1),
             n_crit2 = ceiling(n2), n_crit3 = ceiling(n3), required_n = nn, required_events = ceiling(nn * prev))
}
prev_all <- mean(a$death); prev_t <- mean(a$death[a$tested == 1])
ss <- rbind(cbind(analysis = "Core model (all eligible patients)", riley(length(CORE_TERMS), 0.75, prev_all), available_n = nrow(a), available_events = sum(a$death)),
            cbind(analysis = "Q3 extension 1: LP + 3 markers", riley(5, 0.75, prev_t), available_n = sum(a$tested), available_events = sum(a$death[a$tested == 1])),
            cbind(analysis = "Q3 extension 3: LP + marker splines", riley(8, 0.75, prev_t), available_n = sum(a$tested), available_events = sum(a$death[a$tested == 1])))
ss$satisfied <- ifelse(ss$available_n >= ss$required_n, "yes", "no")
dev_min <- min(vapply(levels(a$cluster), function(k) sum(a$cluster != k), numeric(1)))
ss_fold <- data.frame(note = sprintf("Smallest IECV development fold: %d patients (%d deaths); core model requires %d patients", dev_min,
  min(vapply(levels(a$cluster), function(k) sum(a$death[a$cluster != k]), numeric(1))), ss$required_n[1]))

write_share(list(
  README = data.frame(note = c("Pre-modelling checks: contains no model of predictor-outcome relationships", "cohort: patients and deaths; clusters: IECV clusters",
    "missing: missingness; group_compare: baseline comparison of patients with/without blood gas and with/without markers (no outcome)",
    sprintf("Incomplete cases for core variables %.1f%%; number of imputations m = %d", pct_incomplete, M_IMP), "sample_size: Riley criteria (conservative assumption C = 0.75)")),
  cohort = coh, clusters = cl, missing = miss, abg_tested = abg_n, group_compare = cmp, sample_size = ss, sample_size_fold = ss_fold),
  "01_describe.xlsx")
msg("01 done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min")
