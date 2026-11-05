# =====================================================================
# run_chain.R: one Background Job runs all steps in sequence; it stops at the first error and later steps are not run
#   1 test run of the main analysis (v4_out_test)  2 test run of additional validation  3 final run of the main analysis (v4_out)  4 A12  5 final run of additional validation
# Run: Background Jobs -> Start Background Job -> R Script = this file -> Working Directory = the repository folder -> Start
# Expected time: about 6–9 h with 4 cores
# =====================================================================
CODE_DIR  <- "/path/to/repository/analysis"
SKIP_TEST <- FALSE      # TRUE = skip steps 1 and 2 (test runs) and start with the final run
STEPS <- list(c("run_all.R", "TRUE"), c("run_extra.R", "TRUE"), c("run_all.R", "FALSE"), c("run_A12.R", NA), c("run_extra.R", "FALSE"))
if (SKIP_TEST) STEPS <- STEPS[-(1:2)]

pkgs <- c("mice", "rms", "metafor", "glmnet", "readxl", "writexl", "ggplot2", "patchwork", "ragg", "systemfonts", "gbm", "lme4", "ggalluvial")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing R packages: ", paste(miss, collapse = ", "), ". Install first in the Console: install.packages(c(\"", paste(miss, collapse = "\", \""), "\"))")

run_with <- function(f, test) {
  x <- readLines(file.path(CODE_DIR, f), encoding = "UTF-8")
  if (!is.na(test)) x <- sub("^TEST_RUN <- (TRUE|FALSE)", paste("TEST_RUN <-", test), x)
  cat("\n#####", format(Sys.time(), "%F %T"), "Step: ", f, if (!is.na(test)) paste("TEST_RUN =", test) else "", "#####\n")
  eval(parse(text = x, encoding = "UTF-8"), envir = new.env(parent = globalenv()))
}
t0 <- Sys.time()
for (s in STEPS) run_with(s[1], s[2])
cat("\n##### All steps done in", round(difftime(Sys.time(), t0, units = "hours"), 1), "h #####\n")
cat("result/, figure/ and calculator/ in the output folder contain aggregate results only; local/ must stay on the server\n")
