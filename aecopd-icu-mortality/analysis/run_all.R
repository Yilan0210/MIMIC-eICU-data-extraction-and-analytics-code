# =====================================================================
# run_all.R: run all scripts in RStudio (including RStudio Server)
# Recommended: Background Jobs tab -> Start Background Job -> R Script = this file -> Working Directory = the folder of this file -> Start
# =====================================================================

# ======================== Settings ========================
CODE_DIR <- "/path/to/repository/analysis"         # folder of this file (must contain the R/ subfolder)
DATA_DIR <- "/path/to/data"         # folder containing aecopd_clean_<TAG>.csv or _LOCAL_ONLY.xlsx
CORES    <- 4                                      # cores for parallel IECV (used in a Background Job; set to 1 automatically when sourced in the console)
TEST_RUN <- FALSE   # TRUE = test run (m = 3, B = 20, about 10 min, results in v4_out_test; numbers not usable)
                    # FALSE = final run (m = 33, B = 200; about 3–5 h with 4 cores, 3–4 times longer with 1 core; results in v4_out)
# ===============================================================
# forking processes (mclapply) in the RStudio console can hang, so the console run uses one core; Background Jobs are not affected
if (CORES > 1 && interactive() && identical(Sys.getenv("RSTUDIO"), "1")) {
  message("Running directly in the RStudio console: switching to one core (running this file as a Background Job is faster)"); CORES <- 1
}

Sys.setenv(AECOPD_DATA_DIR = DATA_DIR, AECOPD_CORES = as.character(CORES), AECOPD_TAG = "20261102")
Sys.unsetenv("AECOPD_IECV_HOLDOUT_OUTCOME")
if (TEST_RUN) {
  Sys.setenv(AECOPD_M = "3", AECOPD_MAXIT = "3", AECOPD_B = "20", AECOPD_SENS_M = "2", AECOPD_SENS_B = "10",
             AECOPD_OUT_DIR = file.path(DATA_DIR, "v4_out_test"))
} else {
  for (v in c("AECOPD_M", "AECOPD_MAXIT", "AECOPD_B", "AECOPD_SENS_M", "AECOPD_SENS_B")) Sys.unsetenv(v)
  Sys.setenv(AECOPD_OUT_DIR = file.path(DATA_DIR, "v4_out"))
}
if (!dir.exists(file.path(CODE_DIR, "R"))) stop("Cannot find ", file.path(CODE_DIR, "R"), ": set CODE_DIR above to the folder of this repository")
old_wd <- getwd(); setwd(file.path(CODE_DIR, "R"))
scripts <- c("01_describe.R", "02_develop.R", "03_iecv.R", "04_crossdb_update.R", "05_markers.R", "06_sensitivity.R",
             "07_figures_main.R", "08_figures_supp.R", "09_tables_calculator.R")
t_all <- Sys.time()
cat(if (TEST_RUN) "[TEST RUN]" else "[FINAL RUN]", "Output folder: ", Sys.getenv("AECOPD_OUT_DIR"), "\n")
for (f in scripts) { cat("\n=====", format(Sys.time(), "%F %T"), "start", f, "=====\n"); source(f, local = new.env(parent = globalenv())) }
must <- file.path(Sys.getenv("AECOPD_OUT_DIR"), c(paste0("result/", c("01_describe", "02_development", "03_iecv", "04_crossdb_update", "05_markers",
                  "06_sensitivity", "09_tables"), ".xlsx"), "figure/Figure1_study_flow.pdf", "figure/Figure7_markers.pdf",
                  "figure/eFigure9_missing.pdf", "calculator/aecopd_risk_calculator.html"))
if (any(!file.exists(must))) stop("The following result files were not created; see the errors above:\n", paste(must[!file.exists(must)], collapse = "\n"))
cat("\n=====", format(Sys.time(), "%F %T"), "All done in", round(difftime(Sys.time(), t_all, units = "mins"), 1), "min =====\n")
cat("result/, figure/ and calculator/ contain aggregate results only; local/ contains patient-level data and must stay on the server\n")
setwd(old_wd)
