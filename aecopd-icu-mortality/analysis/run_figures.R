# =====================================================================
# run_figures.R: redraw the figures only (07 main figures, 11 additional validation figures) from the results in v4_out, without re-analysis; about 3–5 min
# Run: Background Jobs -> Start Background Job -> R Script = this file -> Working Directory = the repository folder -> Start
# =====================================================================
CODE_DIR <- "/path/to/repository/analysis"
DATA_DIR <- "/path/to/data"
for (v in c("AECOPD_M", "AECOPD_MAXIT", "AECOPD_B", "AECOPD_SENS_M", "AECOPD_SENS_B", "AECOPD_IECV_HOLDOUT_OUTCOME",
            "AECOPD_EXTRA_B", "AECOPD_V4_B", "AECOPD_SPEC_M", "AECOPD_SPEC_N", "AECOPD_V1_FORWARD")) Sys.unsetenv(v)
Sys.setenv(AECOPD_DATA_DIR = DATA_DIR, AECOPD_CORES = "1", AECOPD_TAG = "20261102", AECOPD_OUT_DIR = file.path(DATA_DIR, "v4_out"))
old_wd <- getwd(); setwd(file.path(CODE_DIR, "R"))
for (f in c("07_figures_main.R", "11_figures_extra.R")) { cat("\n=====", format(Sys.time(), "%T"), f, "=====\n"); source(f, local = new.env(parent = globalenv())) }
setwd(old_wd); cat("\n===== Figures redrawn: see v4_out/figure =====\n")
