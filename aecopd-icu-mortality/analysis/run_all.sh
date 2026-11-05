#!/usr/bin/env bash
# Command-line run (in the folder of this file): nohup bash run_all.sh > v4_run.log 2>&1 &
# Environment variables: AECOPD_DATA_DIR, AECOPD_OUT_DIR, AECOPD_CORES; for a test run add AECOPD_M=3 AECOPD_B=20
set -euo pipefail
cd "$(dirname "$0")/R"
for f in 01_describe.R 02_develop.R 03_iecv.R 04_crossdb_update.R 05_markers.R 06_sensitivity.R 07_figures_main.R 08_figures_supp.R 09_tables_calculator.R; do
  echo "===== $(date '+%F %T') start $f ====="; Rscript "$f"
done
echo "===== $(date '+%F %T') all done ====="
