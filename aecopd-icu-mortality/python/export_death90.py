# -*- coding: utf-8 -*-
r"""
python/export_death90.py: export MIMIC-IV 90-day death (stay_id, death_90d) for sensitivity analysis A10

Usage on Windows (same as run_all.py; dependencies already installed):
  in a Command Prompt, run:
    py path\to\aecopd-icu-mortality\python\export_death90.py
  and enter the database password when prompted.

Output: OUTPUT_DIR\aecopd_mimic_death90_<TAG>_LOCAL_ONLY.csv
  This file contains patient-level data: keep it on your own analysis machine and do not upload it to any online service.
"""
# ======================== Settings to check (same as run_all.py) ========================
DB_HOST = "localhost"
DB_PORT = 5432                         # port of the PostgreSQL server
DB_USER = "postgres"
MIMIC_DB = ""                          # empty = find the database that contains aecopd.m_cohort automatically
OUTPUT_DIR = r"C:\aecopd_data"         # same folder as OUTPUT_DIR in run_all.py
TAG = "20261102"                       # snapshot tag used in the output file name
# ====================================================================================
import getpass
import sys
import warnings
from pathlib import Path

import pandas as pd
import psycopg2

warnings.filterwarnings("ignore", message=".*SQLAlchemy.*")   # pandas warning about psycopg2 connections; does not affect the results
SQL_FILE = Path(__file__).resolve().parent.parent / "sql" / "mimic" / "04_death90.sql"
OUT = Path(OUTPUT_DIR) / f"aecopd_mimic_death90_{TAG}_LOCAL_ONLY.csv"

pw = getpass.getpass(f"PostgreSQL password ({DB_USER}@{DB_HOST}:{DB_PORT}, input is hidden): ")


def connect(db):
    return psycopg2.connect(host=DB_HOST, port=DB_PORT, user=DB_USER, password=pw, dbname=db, connect_timeout=10)


if not MIMIC_DB:
    dbs = []
    for initial in ("postgres", "template1"):
        try:
            con = connect(initial)
            with con.cursor() as cur:
                cur.execute("SELECT datname FROM pg_database WHERE NOT datistemplate AND datallowconn ORDER BY 1")
                dbs = [r[0] for r in cur.fetchall()]
            con.close()
            break
        except psycopg2.OperationalError as e:
            err = e
    if not dbs:
        sys.exit(f"Cannot connect to PostgreSQL at {DB_HOST}:{DB_PORT}: {err}")
    for db in dbs:
        try:
            con = connect(db)
            with con.cursor() as cur:
                cur.execute("SELECT to_regclass('aecopd.m_cohort') IS NOT NULL AND to_regclass('mimiciv_hosp.patients') IS NOT NULL")
                ok = cur.fetchone()[0]
            con.close()
            if ok:
                MIMIC_DB = db
                break
        except psycopg2.Error:
            pass
    if not MIMIC_DB:
        sys.exit("No database contains both aecopd.m_cohort and mimiciv_hosp.patients: run the pipeline first, or set MIMIC_DB")

con = connect(MIMIC_DB)
d = pd.read_sql(SQL_FILE.read_text(encoding="utf-8"), con)
chk = pd.read_sql("SELECT c.stay_id, a.hospital_expire_flag FROM aecopd.m_cohort c "
                  "JOIN mimiciv_hosp.admissions a ON a.hadm_id = c.hadm_id WHERE c.pop = 'main'", con)
con.close()

d["death_90d"] = d["death_90d"].astype(int)
m = d.merge(chk, on="stay_id", how="left")
print(f"\nDatabase {MIMIC_DB}: {len(d)} ICU stays in the main population; {d['death_90d'].sum()} deaths within 90 days ({100 * d['death_90d'].mean():.1f}%)")
print("Check (rows = in-hospital death, columns = 90-day death): almost all in-hospital deaths should be 1; in-hospital deaths with 90-day death = 0 should only be stays longer than 90 days")
print(pd.crosstab(m["hospital_expire_flag"], m["death_90d"]))
if d["stay_id"].duplicated().any():
    sys.exit("Duplicate stay_id values; no file was written. Please check the output above")
OUT.parent.mkdir(parents=True, exist_ok=True)
d.to_csv(OUT, index=False)
print(f"\nWritten: {OUT}\nNext: copy it to your analysis machine (patient-level data; keep it on your own server only)")
input("\nPress Enter to close this window")
