# =====================================================================
# run_A12.R: sensitivity analysis A12 (held-out clusters imputed with the outcome); run after the final run of run_all.R
# Run: Background Jobs -> Start Background Job -> R Script = this file -> Working Directory = the repository folder -> Start
# Output: v4_out/result/03_iecv_A12.xlsx (does not overwrite the main analysis); takes about as long as the 03_iecv step of the main run
# =====================================================================
CODE_DIR <- "/path/to/repository/analysis"
DATA_DIR <- "/path/to/data"
CORES    <- 4
if (CORES > 1 && interactive() && identical(Sys.getenv("RSTUDIO"), "1")) CORES <- 1
for (v in c("AECOPD_M", "AECOPD_MAXIT", "AECOPD_B", "AECOPD_SENS_M", "AECOPD_SENS_B")) Sys.unsetenv(v)
Sys.setenv(AECOPD_DATA_DIR = DATA_DIR, AECOPD_CORES = as.character(CORES), AECOPD_TAG = "20261102",
           AECOPD_OUT_DIR = file.path(DATA_DIR, "v4_out"), AECOPD_IECV_HOLDOUT_OUTCOME = "1")
if (!file.exists(file.path(DATA_DIR, "v4_out", "local", "01_population.rds"))) stop("Complete the final run of run_all.R first")
old_wd <- getwd(); setwd(file.path(CODE_DIR, "R"))
t0 <- Sys.time()
source("03_iecv.R", local = new.env(parent = globalenv()))
Sys.unsetenv("AECOPD_IECV_HOLDOUT_OUTCOME"); setwd(old_wd)
cat("A12 done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
