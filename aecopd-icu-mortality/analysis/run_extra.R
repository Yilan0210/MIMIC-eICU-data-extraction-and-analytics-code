# =====================================================================
# run_extra.R: additional validation V1–V9 and its figures (10_extra_validation.R, 11_figures_extra.R)
# Requires a completed run_all.R (reads its saved results in local/; does not change the main analysis)
# Run: Background Jobs -> Start Background Job -> R Script = this file -> Working Directory = the repository folder -> Start
# Packages needed (once): install.packages(c("lme4", "ggalluvial"))
# =====================================================================
CODE_DIR <- "/path/to/repository/analysis"
DATA_DIR <- "/path/to/data"
CORES    <- 4
TEST_RUN <- FALSE   # TRUE = use v4_out_test (test-run results) to check that the code runs, about 10 min; numbers not usable
                    # FALSE = use v4_out (final results), about 1–1.5 h with 4 cores
if (CORES > 1 && interactive() && identical(Sys.getenv("RSTUDIO"), "1")) CORES <- 1
for (p in c("lme4", "ggalluvial")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing R package ", p, "; install first: install.packages(c('lme4', 'ggalluvial'))")
for (v in c("AECOPD_M", "AECOPD_MAXIT", "AECOPD_B", "AECOPD_SENS_M", "AECOPD_SENS_B", "AECOPD_IECV_HOLDOUT_OUTCOME",
            "AECOPD_EXTRA_B", "AECOPD_V4_B", "AECOPD_SPEC_M", "AECOPD_SPEC_N", "AECOPD_V1_FORWARD")) Sys.unsetenv(v)
Sys.setenv(AECOPD_DATA_DIR = DATA_DIR, AECOPD_CORES = as.character(CORES), AECOPD_TAG = "20261102")
if (TEST_RUN) {
  Sys.setenv(AECOPD_OUT_DIR = file.path(DATA_DIR, "v4_out_test"), AECOPD_M = "3", AECOPD_MAXIT = "3", AECOPD_B = "20",
             AECOPD_EXTRA_B = "200", AECOPD_V4_B = "50", AECOPD_SPEC_M = "1", AECOPD_SPEC_N = "12")
} else Sys.setenv(AECOPD_OUT_DIR = file.path(DATA_DIR, "v4_out"))
need <- file.path(Sys.getenv("AECOPD_OUT_DIR"), "local", c("01_population.rds", "02_dev.rds", "03_iecv_stacked.rds", "03_iecv_results.rds"))
if (any(!file.exists(need))) stop("Main analysis results not found; complete run_all.R first:\n", paste(need[!file.exists(need)], collapse = "\n"))
old_wd <- getwd(); setwd(file.path(CODE_DIR, "R")); t0 <- Sys.time()
cat(if (TEST_RUN) "[TEST RUN]" else "[FINAL]", "Output folder: ", Sys.getenv("AECOPD_OUT_DIR"), "\n")
for (f in c("10_extra_validation.R", "11_figures_extra.R")) {
  cat("\n=====", format(Sys.time(), "%F %T"), "start", f, "=====\n"); source(f, local = new.env(parent = globalenv())) }
setwd(old_wd)
cat("\n===== All done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min =====\n")
cat("Aggregate results: result/10_extra_validation.xlsx, figure/FigX*; local/10_extra.rds contains hospital-level data and must stay on the server\n")
