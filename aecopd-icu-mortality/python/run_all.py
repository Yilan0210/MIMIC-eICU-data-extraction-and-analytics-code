# -*- coding: utf-8 -*-
r"""
python/run_all.py: one-click run: connect to the databases -> name checks -> cohort and analysis tables -> CSV export -> cleaning

Usage on Windows 11:
  1. Install the dependencies once, in a Command Prompt:
       py -m pip install -r path\to\aecopd-icu-mortality\python\requirements.txt
  2. Open this file in a text editor and change only the "Settings" section below (usually only the port)
  3. Double-click this file, or run in a Command Prompt:
       py path\to\aecopd-icu-mortality\python\run_all.py
  4. Enter the database password when prompted (the input is not shown on screen) and press Enter

Database names and the eICU schema are found automatically.
All outputs (snapshots, cleaned data, Excel files) are written to the folder given by OUTPUT_DIR below.
"""

# ======================== Settings ========================
DB_HOST = "localhost"
DB_PORT = 5432          # port of the PostgreSQL server (as in your database client)
DB_USER = "postgres"    # user name (as in your database client)
DB_PASSWORD = ""        # leave empty and type it at run time; never write a password into a file that is uploaded to GitHub

OUTPUT_DIR = r"C:\aecopd_data"   # outputs: raw\ snapshots, clean\ cleaned data, output\ Excel files

MIMIC_DB = ""           # empty = find the database that contains mimiciv_hosp automatically
EICU_DB = ""            # empty = find the database that contains the eICU patient table automatically
EICU_SCHEMA = ""        # empty = detect automatically (usually public or eicu_crd)
EICU_PORT = None        # set only if eICU is on another PostgreSQL server (another port), e.g. 5433

FINAL_SNAPSHOT = False  # False = trial run: file names contain "trial" and can be overwritten
                        # True  = final snapshot: file names contain today's date; existing files are never overwritten
RUN_STEPS = dict(check=True,    # name checks
                 mimic=True,    # MIMIC cohort and analysis table
                 eicu=True,     # eICU cohort and analysis table (skipped automatically if eICU is not found)
                 export=True,   # export CSV snapshots
                 clean=True)    # clean and write Excel files
# ===============================================================

import getpass
import hashlib
import re
import sys
import time
import traceback
from datetime import date, datetime
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))


def pause_and_exit(code):
    if sys.stdin and sys.stdin.isatty():
        input("\nPress Enter to close this window...")
    sys.exit(code)


_missing = []
for _mod, _pip in [("psycopg2", "psycopg2-binary"), ("pandas", "pandas"),
                   ("openpyxl", "openpyxl"), ("numpy", "numpy")]:
    try:
        __import__(_mod)
    except ImportError:
        _missing.append(_pip)
if _missing:
    print("Missing Python packages: " + ", ".join(_missing) + "\nRun in a Command Prompt:\n  py -m pip install "
          + " ".join(_missing) + "\nand run this script again.")
    pause_and_exit(1)

import pandas as pd                       # noqa: E402
import psycopg2                           # noqa: E402
import psycopg2.extensions as pgext       # noqa: E402
from clean import run_clean, write_book   # noqa: E402

# numeric columns are converted to float (the default is Decimal, which is inconvenient for later calculations)
pgext.register_type(pgext.new_type(pgext.DECIMAL.values, "DEC2FLOAT",
                                   lambda v, cur: float(v) if v is not None else None))
INT_OIDS = {20, 21, 23}   # int8 / int2 / int4

REQ_MIMIC = ([f"mimiciv_hosp.{t}" for t in ("patients", "admissions", "diagnoses_icd", "emar", "prescriptions")]
             + [f"mimiciv_icu.{t}" for t in ("icustays", "chartevents", "d_items")]
             + [f"mimiciv_derived.{t}" for t in ("chemistry", "blood_differential", "complete_blood_count",
                                                  "bg", "gcs", "vitalsign", "ventilation", "sapsii",
                                                  "first_day_sofa",
                                                  # vasopressors on D1
                                                  "norepinephrine", "epinephrine", "dopamine",
                                                  "phenylephrine", "vasopressin")])
REQ_EICU = ("patient", "hospital", "diagnosis", "pasthistory", "lab", "nursecharting", "vitalperiodic",
            "vitalaperiodic", "apacheapsvar", "apachepredvar", "apachepatientresult", "medication",
            # admission diagnosis (pneumonia within 24 h), vasopressor infusions, treatment limitation
            "admissiondx", "infusiondrug", "careplangeneral")


# ---------------------------------------------------------------------
# Connection and automatic database discovery
# ---------------------------------------------------------------------
PASSWORD = None


def connect(dbname, port=None):
    try:
        con = psycopg2.connect(host=DB_HOST, port=port or DB_PORT, user=DB_USER,
                               password=PASSWORD, dbname=dbname, connect_timeout=10)
    except psycopg2.OperationalError as e:
        raise SystemExit(
            f"Cannot connect to database {dbname} ({DB_HOST}:{port or DB_PORT}):\n{e}\n"
            "Check that the port and password are correct and that the PostgreSQL service is running"
            " (if your database client can open it, the service is running).")
    con.autocommit = True
    con.set_client_encoding("UTF8")
    with con.cursor() as cur:
        cur.execute("SET client_min_messages TO WARNING")
    return con


def query_df(con, sql):
    with con.cursor() as cur:
        cur.execute(sql)
        cols = [c.name for c in cur.description]
        types = [c.type_code for c in cur.description]
        df = pd.DataFrame(cur.fetchall(), columns=cols)
    for c, t in zip(cols, types):
        if t in INT_OIDS:
            df[c] = df[c].astype("Int64")
    return df


def list_databases(port):
    for initial in ("postgres", "template1"):
        try:
            con = psycopg2.connect(host=DB_HOST, port=port, user=DB_USER, password=PASSWORD,
                                   dbname=initial, connect_timeout=10)
            break
        except psycopg2.OperationalError as e:
            err = e
    else:
        raise SystemExit(f"Cannot connect to PostgreSQL at {DB_HOST}:{port}:\n{err}\n"
                         "Check the port, user name and password.")
    with con.cursor() as cur:
        cur.execute("SELECT datname FROM pg_database WHERE NOT datistemplate AND datallowconn ORDER BY 1")
        names = [r[0] for r in cur.fetchall()]
    con.close()
    return names


def probe(dbname, port):
    """Return (has MIMIC, has MIMIC derived tables, list of eICU schemas)"""
    try:
        con = connect(dbname, port)
    except SystemExit:
        return False, False, []
    with con.cursor() as cur:
        cur.execute("SELECT to_regclass('mimiciv_hosp.patients') IS NOT NULL, "
                    "to_regclass('mimiciv_derived.blood_differential') IS NOT NULL")
        is_mimic, has_derived = cur.fetchone()
        cur.execute("""SELECT DISTINCT table_schema FROM information_schema.columns
                       WHERE table_name = 'patient' AND column_name = 'patientunitstayid'
                       ORDER BY 1""")
        eicu_schemas = [r[0] for r in cur.fetchall()]
    con.close()
    return is_mimic, has_derived, eicu_schemas


def locate():
    """Determine the MIMIC database, the eICU database and the eICU schema"""
    eport = EICU_PORT or DB_PORT
    ports = sorted({DB_PORT, eport})
    found = {}
    for port in ports:
        for db in list_databases(port):
            found[(port, db)] = probe(db, port)
    print("  Databases on the server:")
    for (port, db), (m, dv, es) in found.items():
        tags = (["MIMIC-IV" + (" (with derived tables)" if dv else " (no derived tables)")] if m else []) + \
               [f"eICU (schema {s})" for s in es]
        print(f"    port {port} / {db}: {', '.join(tags) if tags else '-'}")

    mimic_db = MIMIC_DB or next((db for (p, db), v in found.items() if p == DB_PORT and v[0]), "")
    eicu_db, eicu_schema = EICU_DB, EICU_SCHEMA
    if not eicu_db:
        hit = next(((db, v[2]) for (p, db), v in found.items() if p == eport and v[2]), None)
        if hit:
            eicu_db = hit[0]
            eicu_schema = eicu_schema or hit[1][0]
    elif not eicu_schema:
        eicu_schema = (found.get((eport, eicu_db), (0, 0, ["public"]))[2] or ["public"])[0]
    return mimic_db, eicu_db, eicu_schema or "public"


def preflight_mimic(con):
    with con.cursor() as cur:
        lacking = []
        for t in REQ_MIMIC:
            cur.execute("SELECT to_regclass(%s) IS NULL", (t,))
            if cur.fetchone()[0]:
                lacking.append(t)
    if lacking:
        msg = "Tables missing from the MIMIC database:\n  " + "\n  ".join(lacking)
        if any("mimiciv_derived" in t for t in lacking):
            msg += ("\nThe MIMIC derived tables (mimiciv_derived) must first be built with the mimic-code PostgreSQL concepts scripts"
                    " (mimic-iv/concepts_postgres in github.com/MIT-LCP/mimic-code).")
        raise SystemExit(msg)
    print(f"  MIMIC: all {len(REQ_MIMIC)} required tables found")


def preflight_eicu(con, schema):
    with con.cursor() as cur:
        lacking = []
        for t in REQ_EICU:
            cur.execute("SELECT to_regclass(%s) IS NULL", (f"{schema}.{t}",))
            if cur.fetchone()[0]:
                lacking.append(f"{schema}.{t}")
    if lacking:
        raise SystemExit("Tables missing from the eICU database:\n  " + "\n  ".join(lacking))
    print(f"  eICU: all {len(REQ_EICU)} required tables found (schema {schema})")


# ---------------------------------------------------------------------
# Run SQL files
# ---------------------------------------------------------------------
def read_sql_statements(path, schema_from=None, schema_to=None):
    stmts, labels, buf, label = [], [], [], ""
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        m = re.match(r"^--\s*(C[0-9]+)\.", ln)
        if m:
            label = m.group(1)
        code = re.sub(r"--.*$", "", ln)            # the SQL strings in this project contain no "--"
        if not code.strip():
            continue
        buf.append(code)
        if re.search(r";\s*$", code):
            s = re.sub(r";\s*$", "", "\n".join(buf)).strip()
            if schema_from and schema_to and schema_from != schema_to:
                s = s.replace(f"{schema_from}.", f"{schema_to}.")
            stmts.append(s)
            labels.append(label)
            buf = []
    if buf:
        raise ValueError(f"{path} ends with a statement that is not terminated by a semicolon")
    return stmts, labels


def describe(s):
    m = re.match(r"^(CREATE (UNIQUE )?INDEX ON|CREATE SCHEMA( IF NOT EXISTS)?|CREATE TABLE|"
                 r"DROP TABLE( IF EXISTS)?)\s+[A-Za-z0-9_.]+", s)
    return m.group(0) if m else "summary query"


def run_sql_file(con, relpath, schema_from=None, schema_to=None):
    stmts, labels = read_sql_statements(PROJECT_DIR / relpath, schema_from, schema_to)
    print(f"\n>> {relpath} ({len(stmts)} statements)")
    out = {}
    for i, (s, lab) in enumerate(zip(stmts, labels), start=1):
        kw = re.match(r"^[A-Za-z]+", s).group(0).upper()
        t1 = time.time()
        try:
            if kw in ("SELECT", "WITH"):
                res = query_df(con, s)
            else:
                with con.cursor() as cur:
                    cur.execute(s)
                res = None
        except psycopg2.Error as e:
            raise SystemExit(f"\nFailed: {relpath}, statement {i}\n{e}\nStart of statement:\n{s[:400]}")
        print(f"   {i:2d}/{len(stmts)}  {time.time() - t1:7.1f} s  {describe(s)}", flush=True)
        if res is not None:
            name = lab or f"q{i}"
            while name in out:
                name += "b"
            out[name] = res
    return out


def export_table(con, table, path, final):
    if final and path.exists():
        raise SystemExit(f"Final snapshot {path} already exists and will not be overwritten. To export again, rename or delete the old file first.")
    df = query_df(con, f"SELECT * FROM {table} ORDER BY stay_id")
    df.to_csv(path, index=False, na_rep="", encoding="utf-8")
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    print(f"   Exported {table}: {len(df)} rows x {df.shape[1]} columns -> {path}")
    return dict(file=str(path), rows=len(df), cols=df.shape[1], md5=md5,
                exported_at=f"{datetime.now():%Y-%m-%d %H:%M:%S}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    global PASSWORD
    t0 = time.time()
    tag = date.today().strftime("%Y%m%d") if FINAL_SNAPSHOT else "trial"
    data_root = Path(OUTPUT_DIR)
    raw_dir = data_root / "raw"
    out_dir = data_root / "output"
    try:
        raw_dir.mkdir(parents=True, exist_ok=True)
        out_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        raise SystemExit(f"Cannot create the output folder {data_root}: {e}\nCheck OUTPUT_DIR in the settings (for example, that the drive exists).")
    print(f"Code folder: {PROJECT_DIR}\nOutput folder: {data_root}\nFile tag: {tag}")
    if not (PROJECT_DIR / "sql" / "eicu" / "01_cohort.sql").exists():
        raise SystemExit("The sql folder was not found. Keep the repository structure (the python and sql folders at the same level).")

    snaps = [raw_dir / f"mimic_analytic_{tag}.csv", raw_dir / f"eicu_analytic_{tag}.csv"]
    if FINAL_SNAPSHOT and RUN_STEPS["export"] and any(p.exists() for p in snaps):
        raise SystemExit("Today's final snapshot already exists and will not be overwritten. To export again, rename or delete the old files first.")

    PASSWORD = DB_PASSWORD or getpass.getpass("Database password (input is hidden): ")

    print("\n>> Finding databases")
    mimic_db, eicu_db, eicu_schema = locate()
    eport = EICU_PORT or DB_PORT
    if not mimic_db:
        raise SystemExit(f"MIMIC-IV (mimiciv_hosp.patients) was not found on port {DB_PORT}. Check the port, "
                         "or set MIMIC_DB in the settings.")
    print(f"  Using MIMIC: {mimic_db}")
    do_eicu = bool(eicu_db)
    if do_eicu:
        print(f"  Using eICU: {eicu_db} (port {eport}, schema {eicu_schema})")
    else:
        print("  eICU not found: only MIMIC is processed. Run again after eICU has been installed.")

    con_m = connect(mimic_db)
    preflight_mimic(con_m)
    con_e = None
    if do_eicu:
        con_e = connect(eicu_db, eport)
        preflight_eicu(con_e, eicu_schema)

    results = {}
    try:
        if RUN_STEPS["check"]:
            r = run_sql_file(con_m, "sql/mimic/00_check_names.sql")
            results.update({f"MIMIC_check_{k}": v for k, v in r.items()})
            if do_eicu:
                r = run_sql_file(con_e, "sql/eicu/00_check_names.sql", "public", eicu_schema)
                results.update({f"eICU_check_{k}": v for k, v in r.items()})
        if RUN_STEPS["mimic"]:
            r = run_sql_file(con_m, "sql/mimic/01_cohort.sql")
            results.update({f"MIMIC_cohort_{i}": v for i, v in enumerate(r.values(), 1)})
            r = run_sql_file(con_m, "sql/mimic/02_analytic.sql")
            results.update({"MIMIC_flowchart_SQL": v for v in r.values()})
            with con_m.cursor() as cur:
                cur.execute("SELECT to_regclass('mimiciv_ed.pyxis') IS NOT NULL "
                            "AND to_regclass('mimiciv_ed.edstays') IS NOT NULL "
                            "AND to_regclass('mimiciv_ed.diagnosis') IS NOT NULL")
                has_ed = cur.fetchone()[0]
            if has_ed:
                r = run_sql_file(con_m, "sql/mimic/03_ed_steroid.sql")
                # 03 contains several summary queries (corticosteroids, pneumonia); each is saved separately
                names = ["MIMIC_ED_steroid", "MIMIC_ED_pneumonia", "MIMIC_ED_pna_by_year", "MIMIC_ED_pna_codes"]
                results.update({(names[i] if i < len(names) else f"MIMIC_ED_{k}"): v
                                for i, (k, v) in enumerate(r.items())})
            else:
                print("\n  MIMIC-IV-ED module (mimiciv_ed) not found: corticosteroids given in the emergency department cannot be added, "
                      "corticosteroid exposure comes from inpatient medication records only, and pneumonia_24h is empty for all MIMIC stays")
        if RUN_STEPS["eicu"] and do_eicu:
            r = run_sql_file(con_e, "sql/eicu/01_cohort.sql", "public", eicu_schema)
            results.update({f"eICU_cohort_{i}": v for i, v in enumerate(r.values(), 1)})
            r = run_sql_file(con_e, "sql/eicu/02_analytic.sql", "public", eicu_schema)
            results.update({"eICU_flowchart_SQL": v for v in r.values()})
        if results:
            f = out_dir / f"sql_results_{tag}_SHAREABLE.xlsx"
            write_book(f, results)
            print(f"\nName checks and SQL summaries saved: {f}")

        if RUN_STEPS["export"]:
            print("\n>> Exporting CSV snapshots")
            log = [export_table(con_m, "aecopd.m_analytic", snaps[0], FINAL_SNAPSHOT)]
            if do_eicu:
                log.append(export_table(con_e, "aecopd.e_analytic", snaps[1], FINAL_SNAPSHOT))
            pd.DataFrame(log).to_csv(raw_dir / f"export_log_{tag}.csv", index=False, encoding="utf-8")
            print(f"   Export log (with MD5 checksums): {raw_dir / f'export_log_{tag}.csv'}")
    finally:
        con_m.close()
        if con_e is not None:
            con_e.close()

    if RUN_STEPS["clean"]:
        print("\n>> Cleaning")
        run_clean(data_root, tag)

    print(f"\nDone in {(time.time() - t0) / 60:.1f} min.")
    print(f"Outputs are in {data_root}: raw\\ snapshots, clean\\ cleaned data, output\\ Excel files.\n"
          "Only the Excel files in output whose names contain SHAREABLE (aggregate results) may be shared;\n"
          "raw, clean and LOCAL_ONLY files contain patient-level data and must stay on this computer (not in cloud-synced folders).")


if __name__ == "__main__":
    code = 0
    try:
        main()
    except SystemExit as e:
        if e.code not in (None, 0):
            print("\nError:\n" + str(e.code))
            code = 1
    except Exception:
        print("\nError (unexpected exception). Copy the text below when asking for help:")
        traceback.print_exc()
        code = 1
    pause_and_exit(code)
