# -*- coding: utf-8 -*-
"""
python/clean.py: merging, quality control and cleaning of the analysis tables from both databases

Usually called by run_all.py; it can also be run on its own (the output folder defaults to OUTPUT_DIR in run_all.py):
    py clean.py trial          (trial snapshot)
    py clean.py 20261101       (final snapshot)

Run locally only. Reads the snapshots in <output folder>/raw and writes:
    <output folder>/clean/aecopd_clean_<tag>.csv               patient-level, local only
    <output folder>/output/aecopd_clean_<tag>_LOCAL_ONLY.xlsx  patient-level, local only (includes the manual audit sheets)
    <output folder>/output/qc_summary_<tag>_SHAREABLE.xlsx     aggregate numbers only, may be shared
"""
import math
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------
# Parameters (fixed before any model fitting)
# ---------------------------------------------------------------------
GCS_VAR = "gcs_motor_first_d1"      # core model: first GCS motor score (1–6) in D1; sensitivity analysis uses gcs_valid_first_d1
BAP_GCS_VAR = "gcs_valid_first_d1"  # BAP-65 "altered mental status" always uses a valid total score < 14, independent of GCS_VAR
AF_VAR = "af_24h"          # definition available within 24 h of ICU admission (af_any = any diagnosis during the stay)
ZERO_FLOOR = 0.01          # value substituted for N or L = 0 before log transformation (x10^9/L)
EOS_LOW_CUT = 0.05         # eosinopenia: < 0.05 x10^9/L (DECAF definition), binary
REQUIRE_EOS = False        # require an eosinophil count for the analysis population? False = missing eosinophils are imputed (no loss of patients)
SEED = 20261101            # random seed for the manual audit sample

# Clusters for internal-external cross-validation (IECV)
CLUSTER_LABELS = {1: "MIMIC ICD-9 era", 2: "MIMIC ICD-10 era", 3: "eICU Midwest", 4: "eICU Northeast",
                  5: "eICU South", 6: "eICU West", 7: "eICU region missing"}
EICU_REGION_CLUSTER = {"Midwest": 3, "Northeast": 4, "South": 5, "West": 6}

EXPECTED_COLS = [
    "db", "pop", "stay_id", "hadm_id", "subject_id", "hospital_id", "region", "teaching",
    "bed_cat", "coding_era", "year_grp", "incl_age", "incl_los24", "incl_outcome", "excl_covid",
    "excl_heme", "excl_hiv", "excl_tx", "excl_cirr", "has_alb", "has_diff", "aecopd_primary",
    "dx_by_string", "dx_by_icd", "death", "discharge_loc", "age_raw", "age_gt89", "sex",
    "race_raw", "bun", "bun_offset_h", "gcs_min_d1", "gcs_first_d1", "gcs_valid_first_d1",
    "gcs_motor_first_d1", "gcs_unable_d1",
    "gcs_source", "hr_first_d1", "ph", "ph_offset_h", "ph_specimen", "af_any", "af_24h",
    "alb_gdl", "alb_offset_h", "neut_abs", "lymph_abs", "eos_abs", "neut_pct", "lymph_pct",
    "eos_pct", "wbc_diff", "diff_offset_h", "wbc_lag_min", "paco2", "paco2_offset_h", "lactate",
    "lactate_offset_h", "creatinine", "creatinine_offset_h", "hgb", "hgb_offset_h", "plt",
    "plt_offset_h", "wbc", "wbc_offset_h", "spo2_first_d1", "rr_first_d1", "map_first_d1",
    "temp_first_d1", "invasive_vent_d1", "pneumonia_any", "steroid_before_diff",
    "steroid_source", "apache_iva_prob", "sapsii", "sofa_d1",
    # additional variables
    "pneumonia_24h", "vasopressor_d1", "limit_tx_24h", "glucose", "glucose_offset_h", "mono_abs",
    "med_data_avail", "infusion_data_avail", "careplan_data_avail"]

CHR_COLS = ["db", "pop", "stay_id", "hadm_id", "subject_id", "hospital_id", "region", "bed_cat",
            "coding_era", "year_grp", "discharge_loc", "sex", "race_raw", "gcs_source",
            "ph_specimen", "steroid_source"]
INT_COLS = ["teaching", "incl_age", "incl_los24", "incl_outcome", "excl_covid", "excl_heme",
            "excl_hiv", "excl_tx", "excl_cirr", "has_alb", "has_diff", "aecopd_primary",
            "dx_by_string", "dx_by_icd", "death", "age_gt89", "gcs_unable_d1", "af_any", "af_24h",
            "invasive_vent_d1", "pneumonia_any", "steroid_before_diff",
            "pneumonia_24h", "vasopressor_d1", "limit_tx_24h",
            "med_data_avail", "infusion_data_avail", "careplan_data_avail"]

# Quality-control ranges (values outside are set to missing)
QC_RULES = [
    ("age", 18, 90, "years", "core"),
    ("bun", 1, 300, "mg/dL", "core"),
    ("gcs_min_d1", 3, 15, "points", "core"),
    ("gcs_first_d1", 3, 15, "points", "core"),
    ("gcs_valid_first_d1", 3, 15, "points", "core"),
    ("gcs_motor_first_d1", 1, 6, "points", "core model"),
    ("hr_first_d1", 20, 250, "/min", "core"),
    ("ph", 6.8, 7.8, "", "core"),
    ("neut_abs", 0, 100, "x10^9/L", "core"),
    ("lymph_abs", 0, 50, "x10^9/L", "core"),
    ("eos_abs", 0, 10, "x10^9/L", "core"),
    ("alb_gl", 10, 60, "g/L", "core"),
    ("plt", 1, 2000, "x10^9/L", "core"),
    ("neut_pct", 0, 100, "%", "auxiliary"),
    ("lymph_pct", 0, 100, "%", "auxiliary"),
    ("eos_pct", 0, 100, "%", "auxiliary"),
    ("paco2", 10, 200, "mmHg", "auxiliary"),
    ("lactate", 0.1, 30, "mmol/L", "auxiliary"),
    ("creatinine", 0.1, 25, "mg/dL", "auxiliary"),
    ("hgb", 2, 25, "g/dL", "auxiliary"),
    ("wbc", 0.1, 500, "x10^9/L", "auxiliary"),
    ("wbc_diff", 0.1, 500, "x10^9/L", "auxiliary"),
    ("spo2_first_d1", 20, 100, "%", "auxiliary"),
    ("rr_first_d1", 2, 80, "/min", "auxiliary"),
    ("map_first_d1", 20, 250, "mmHg", "auxiliary"),
    ("temp_first_d1", 25, 45, "°C", "auxiliary"),
    ("glucose", 10, 1500, "mg/dL", "additional"),
    ("mono_abs", 0, 50, "x10^9/L", "additional"),
]


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------
def is_val(s, v):
    """Equivalent to R's x %in% v: missing values count as False."""
    return (s == v).fillna(False).astype(bool)


def read_snapshot(path):
    path = Path(path)
    hdr = list(pd.read_csv(path, nrows=0, encoding="utf-8").columns)
    miss = [c for c in EXPECTED_COLS if c not in hdr]
    extra = [c for c in hdr if c not in EXPECTED_COLS]
    if miss or extra:
        raise ValueError(f"{path}: unexpected columns. Missing: {miss}; extra: {extra}")
    dtypes = {c: "string" for c in CHR_COLS}
    dtypes.update({c: "Int64" for c in INT_COLS})
    dtypes.update({c: "float64" for c in EXPECTED_COLS if c not in CHR_COLS + INT_COLS})
    df = pd.read_csv(path, dtype=dtypes, keep_default_na=False, na_values=["", "NA"],
                     encoding="utf-8")
    return df[EXPECTED_COLS]


def map_race(db, raw):
    x0 = raw.fillna("").astype(str).str.strip().str.upper()
    eicu = db == "eICU"
    conds = [
        eicu & (x0 == "CAUCASIAN"),
        eicu & (x0 == "AFRICAN AMERICAN"),
        eicu & x0.isin(["", "OTHER/UNKNOWN"]),
        eicu,
        x0.str.startswith("WHITE"),
        x0.str.startswith("BLACK"),
        x0.isin(["", "UNKNOWN", "UNABLE TO OBTAIN", "PATIENT DECLINED TO ANSWER"]),
    ]
    choices = ["White", "Black", "Unknown", "Other", "White", "Black", "Unknown"]
    return pd.Series(np.select(conds, choices, default="Other"), index=raw.index)


def masked_bool(cond, *inputs):
    """Comparison result as 0/1 float; NaN if any input is missing (like NA propagation in R)."""
    out = cond.astype(float)
    for s in inputs:
        out[s.isna()] = np.nan
    return out


def write_book(path, sheets):
    from openpyxl.styles import Font, PatternFill, Border, Side
    from openpyxl.utils import get_column_letter
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name[:31], index=False)
        head_font = Font(name="Arial", size=10, bold=True)
        body_font = Font(name="Arial", size=10)
        fill = PatternFill("solid", fgColor="D9E1F2")
        border = Border(bottom=Side(style="thin"))
        for name, df in sheets.items():
            ws = xw.sheets[name[:31]]
            for cell in ws[1]:
                cell.font, cell.fill, cell.border = head_font, fill, border
            if ws.max_row * ws.max_column <= 100_000:      # skip cell-by-cell fonts for large sheets to save time
                for row in ws.iter_rows(min_row=2):
                    for cell in row:
                        cell.font = body_font
            ws.freeze_panes = "A2"
            for j, col in enumerate(df.columns, start=1):
                sample = [len(str(v)) for v in df[col].head(200).tolist()]
                width = max([len(str(col))] + sample) + 2
                ws.column_dimensions[get_column_letter(j)].width = min(max(width, 8), 60)


# ---------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------
def run_clean(data_root, tag):
    data_root = Path(data_root)
    raw_dir = data_root / "raw"
    clean_dir = data_root / "clean"
    out_dir = data_root / "output"
    clean_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Read (whichever databases are available, at least one)
    parts = []
    for fn in (f"mimic_analytic_{tag}.csv", f"eicu_analytic_{tag}.csv"):
        p = raw_dir / fn
        if p.exists():
            parts.append(read_snapshot(p))
        else:
            print(f"  Note: {p.name} not found; this database is not included in this run")
    if not parts:
        raise FileNotFoundError(f"No snapshot files with tag {tag} in {raw_dir}")
    d = pd.concat(parts, ignore_index=True)
    dbs = [k for k in ("MIMIC", "eICU") if (d["db"] == k).any()]
    print(f"  Read: " + ", ".join(f"{k} {int((d['db'] == k).sum())} rows" for k in dbs))

    # 2. Structure checks
    if not d["db"].isin(["MIMIC", "eICU"]).all():
        raise ValueError("Unexpected values in column db")
    if d.duplicated(["db", "stay_id"]).any():
        raise ValueError(f"{int(d.duplicated(['db', 'stay_id']).sum())} duplicated ICU stay IDs; stopping")
    if d.duplicated(["db", "pop", "subject_id"]).any():
        raise ValueError("Duplicated patient IDs within a population; stopping")

    # 3. Harmonise coding and units
    d["age"] = d["age_raw"].clip(upper=90)
    d["age_capped"] = (d["age_raw"] > 90).fillna(False).astype(int)
    d["alb_gl"] = d["alb_gdl"] * 10
    d["sex"] = d["sex"].where(d["sex"].isin(["M", "F"]))
    d["race4"] = map_race(d["db"], d["race_raw"])
    d["hospice"] = pd.array([pd.NA] * len(d), dtype="Int64")
    d["death_s3"] = pd.array([pd.NA] * len(d), dtype="Int64")
    for k in dbs:
        idx = d["db"] == k
        has_cat = d.loc[idx, "discharge_loc"].str.contains("hospice", case=False, na=False).any()
        if has_cat:
            h = d.loc[idx, "discharge_loc"].str.contains("hospice", case=False, na=False).astype(int)
            d.loc[idx, "hospice"] = h.values
            dth = d.loc[idx, "death"]
            s3 = pd.Series(np.maximum(dth.astype("float").values, h.values), index=dth.index)
            s3[dth.isna()] = np.nan
            d.loc[idx, "death_s3"] = s3.astype("Int64").values

    # 4. Quality control
    qc_bad = {}
    for var, lo, hi, _, _ in QC_RULES:
        x = d[var]
        bad = x.notna() & ((x < lo) | (x > hi))
        qc_bad[var] = bad.fillna(False).astype(bool)
        d.loc[qc_bad[var], var] = np.nan

    diff_cols = ["neut_abs", "lymph_abs", "eos_abs", "neut_pct", "lymph_pct", "eos_pct", "wbc_diff",
                 "mono_abs"]
    pct_sum = d["neut_pct"] + d["lymph_pct"]
    d["diff_pct_sum_gt100"] = (pct_sum > 100).fillna(False).astype(int)
    pct_bad = (d["db"] == "eICU") & (d["diff_pct_sum_gt100"] == 1)
    d.loc[pct_bad, diff_cols] = np.nan

    # 5. Population flags
    d["has_alb_qc"] = (is_val(d["has_alb"], 1) & d["alb_gl"].notna()).astype(int)
    has_diff_qc = is_val(d["has_diff"], 1) & d["neut_abs"].notna() & d["lymph_abs"].notna()
    if REQUIRE_EOS:
        has_diff_qc = has_diff_qc & d["eos_abs"].notna()
    d["has_diff_qc"] = has_diff_qc.astype(int)
    pass_base = (is_val(d["incl_age"], 1) & is_val(d["incl_los24"], 1) & is_val(d["incl_outcome"], 1)
                 & is_val(d["excl_covid"], 0) & is_val(d["excl_heme"], 0) & is_val(d["excl_hiv"], 0)
                 & is_val(d["excl_tx"], 0) & is_val(d["excl_cirr"], 0))
    d["in_s2"] = pass_base.astype(int)
    d["in_analysis"] = (pass_base & (d["has_alb_qc"] == 1) & (d["has_diff_qc"] == 1)).astype(int)

    # 6. Derived variables
    d["gcs"] = d[GCS_VAR]
    d["af"] = d[AF_VAR]
    n_log = d["neut_abs"].where(d["neut_abs"] != 0, ZERO_FLOOR)
    l_log = d["lymph_abs"].where(d["lymph_abs"] != 0, ZERO_FLOOR)
    with np.errstate(divide="ignore", invalid="ignore"):
        d["ln_bun"] = np.log(d["bun"])
        d["ln_n"] = np.log(n_log)
        d["ln_l"] = np.log(l_log)
        d["ln_e"] = np.log(d["eos_abs"] + 0.01)
        d["nlr"] = n_log / l_log
        d["ln_nlr"] = np.log(d["nlr"])
        d["pni"] = d["alb_gl"] + 5 * d["lymph_abs"]
        d["sii"] = d["plt"] * n_log / l_log
        d["ln_sii"] = np.log(d["sii"])
        d["npar"] = d["neut_pct"] / (d["alb_gl"] / 10)
        # indices for comparison (PLR, SIRI, GLR, LAR)
        d["plr"] = d["plt"] / l_log
        d["siri"] = n_log * d["mono_abs"] / l_log
        d["glr"] = (d["glucose"] / 18.016) / l_log                 # glucose (mmol/L) / lymphocytes
        d["lar"] = d["lactate"] / (d["alb_gl"] / 10)               # lactate (mmol/L) / albumin (g/dL)
    # Treatment-limitation sensitivity analysis: only patients with a documented DNR/DNI/comfort-care record are excluded; no code-status record within 24 h counts as "no limitation"
    d["limit_tx_s18"] = is_val(d["limit_tx_24h"], 1).astype(int)
    # Eosinopenia (binary)
    d["eos_low"] = masked_bool(d["eos_abs"] < EOS_LOW_CUT, d["eos_abs"])
    # BAP-65: altered mental status always uses a valid GCS total score < 14 (BAP_GCS_VAR), independent of the GCS_VAR used in the core model (motor score)
    n_bap = (masked_bool(d["bun"] >= 25, d["bun"])
             + masked_bool(d[BAP_GCS_VAR] < 14, d[BAP_GCS_VAR])
             + masked_bool(d["hr_first_d1"] >= 109, d["hr_first_d1"]))
    bap = pd.Series(np.nan, index=d.index)
    ok = n_bap.notna() & d["age"].notna()
    bap[ok & (n_bap == 0) & (d["age"] < 65)] = 1
    bap[ok & (n_bap == 0) & (d["age"] >= 65)] = 2
    bap[ok & (n_bap > 0)] = n_bap[ok & (n_bap > 0)] + 2
    d["bap65"] = bap.astype("Int64")
    # Modified DECAF: dyspnoea grade omitted; eosinopenia + pneumonia within 24 h + pH < 7.30 + atrial fibrillation within 24 h, 1 point each
    # Complete-case convenience only; recomputed after imputation in the analysis
    d["decaf_mod"] = (d["eos_low"] + masked_bool(is_val(d["pneumonia_24h"], 1), d["pneumonia_24h"])
                      + masked_bool(d["ph"] < 7.30, d["ph"])
                      + masked_bool(is_val(d["af_24h"], 1), d["af_24h"]))
    # IECV clusters: MIMIC by coding era, eICU by hospital region (missing region = cluster 7)
    mim = d["db"] == "MIMIC"
    cl = pd.Series(np.nan, index=d.index)
    cl[mim & is_val(d["coding_era"], "ICD-9")] = 1
    cl[mim & is_val(d["coding_era"], "ICD-10")] = 2
    ei_mask = d["db"] == "eICU"
    cl[ei_mask] = d.loc[ei_mask, "region"].map(EICU_REGION_CLUSTER).fillna(7).astype(float)
    d["cluster"] = cl.astype("Int64")
    for v, off in [("bun", "bun_offset_h"), ("ph", "ph_offset_h"), ("alb_gl", "alb_offset_h"),
                   ("neut_abs", "diff_offset_h"), ("lymph_abs", "diff_offset_h"),
                   ("eos_abs", "diff_offset_h")]:
        d[f"{v}_s4"] = d[v].where(d[off].abs() <= 24)

    # 7. Summary tables (counts and proportions only)
    flow_steps = [
        ("0 First ICU stay with AECOPD (one per patient)", lambda x: pd.Series(True, index=x.index)),
        ("1 Age >= 18 years", lambda x: is_val(x["incl_age"], 1)),
        ("2 ICU stay >= 24 h", lambda x: is_val(x["incl_los24"], 1)),
        ("3 Outcome known", lambda x: is_val(x["incl_outcome"], 1)),
        ("4 Exclude COVID-19", lambda x: is_val(x["excl_covid"], 0)),
        ("5 Exclude haematological malignancy", lambda x: is_val(x["excl_heme"], 0)),
        ("6 Exclude HIV/AIDS", lambda x: is_val(x["excl_hiv"], 0)),
        ("7 Exclude solid organ or bone marrow transplantation", lambda x: is_val(x["excl_tx"], 0)),
        ("8 Exclude cirrhosis", lambda x: is_val(x["excl_cirr"], 0)),
        ("9 Albumin in W4 (after QC)", lambda x: x["has_alb_qc"] == 1),
        ("10 Differential in W4 (after QC)", lambda x: x["has_diff_qc"] == 1),
    ]
    flow_rows = []
    for (k, pop), sub in d.groupby(["db", "pop"], sort=True):
        keep = pd.Series(True, index=sub.index)
        prev = len(sub)
        for name, f in flow_steps:
            keep = keep & f(sub)
            n = int(keep.sum())
            deaths = int(is_val(sub.loc[keep, "death"], 1).sum())
            flow_rows.append(dict(db=k, pop=pop, step=name, n=n, excluded=prev - n, deaths=deaths,
                                  mortality_pct=round(100 * deaths / n, 1) if n else np.nan))
            prev = n
    flowchart = pd.DataFrame(flow_rows)

    rp = d[(d["db"] == "MIMIC") & (d["pop"] == "main") & is_val(d["incl_age"], 1)
           & is_val(d["incl_los24"], 1) & is_val(d["excl_covid"], 0)
           & is_val(d["has_alb"], 1) & is_val(d["has_diff"], 1)]
    repro = pd.DataFrame([dict(n=len(rp), deaths=int(is_val(rp["death"], 1).sum()),
                               previous_n=690, previous_deaths=120,
                               note="If the difference is large, check the code selection, the W4 window and emergency department laboratory results without hadm_id")])

    ana = d[(d["pop"] == "main") & (d["in_analysis"] == 1)].copy()
    events = (ana.groupby("db").agg(n=("death", "size"), deaths=("death", "sum"),
                                    mortality_pct=("death", "mean")).reset_index())
    events["deaths"] = events["deaths"].astype(int)
    events["mortality_pct"] = (100 * events["mortality_pct"].astype(float)).round(1)
    e_by = dict(zip(events["db"], events["deaths"]))
    # Scenario C (MIMIC alone has about 100 events, below the about 125 required by the sample size calculation)
    scenario = "C: development in both databases combined, internal-external cross-validation (IECV, 7 clusters)"
    events_scn = pd.DataFrame(
        [dict(item=f"E({k})", value=str(v)) for k, v in e_by.items()]
        + [dict(item="Development and validation strategy", value=scenario),
           dict(item="Note", value="Strategy chosen before any model fitting")])

    # The core model includes invasive ventilation on D1; GCS is the motor score (GCS_VAR)
    pred_vars = ["age", "bun", GCS_VAR, "hr_first_d1", "ph", AF_VAR, "invasive_vent_d1",
                 "alb_gl", "neut_abs", "lymph_abs", "eos_abs"]
    aux_vars = ["sex", BAP_GCS_VAR, "paco2", "lactate", "creatinine", "hgb", "plt", "wbc",
                "spo2_first_d1", "rr_first_d1", "map_first_d1", "temp_first_d1",
                "steroid_before_diff", "apache_iva_prob", "sapsii", "sofa_d1"]
    new_vars = ["pneumonia_24h", "vasopressor_d1", "limit_tx_24h", "glucose", "mono_abs"]
    structural = {"MIMIC": ["apache_iva_prob"], "eICU": ["sapsii", "sofa_d1"]}
    mrows = []
    for k in sorted(ana["db"].unique()):
        sub = ana[ana["db"] == k]
        for v in pred_vars + aux_vars + new_vars:
            nm = int(sub[v].isna().sum())
            st = v in structural.get(k, [])
            mrows.append(dict(db=k, var=v, n=len(sub), n_missing=nm,
                              pct_missing=np.nan if st else round(100 * nm / max(len(sub), 1), 1),
                              note="not available in this database" if st else "",
                              role=("predictor" if v in pred_vars else
                                    "additional (sensitivity/secondary)" if v in new_vars else "auxiliary/comparator")))
    missing_tab = pd.DataFrame(mrows)
    pm = missing_tab.loc[missing_tab["var"].isin(pred_vars), "pct_missing"].dropna()
    m_imp = max(20, math.ceil(pm.max())) if len(pm) else 20
    missing_tab = pd.concat([missing_tab, pd.DataFrame([dict(
        db="both", var="suggested number of imputations m", n=np.nan, n_missing=np.nan, pct_missing=m_imp,
        note="max(20, highest percentage of missing values among predictors)", role="")])], ignore_index=True)

    qrows = []
    for var, lo, hi, unit, basis in QC_RULES:
        for k in ("MIMIC", "eICU"):
            idx_all = d["db"] == k
            idx_s2 = idx_all & (d["pop"] == "main") & (d["in_s2"] == 1)
            qrows.append(dict(db=k, var=var, range=f"{lo}–{hi}", unit=unit, basis=basis,
                              set_missing_all_rows=int(qc_bad[var][idx_all].sum()),
                              set_missing_main_s2=int(qc_bad[var][idx_s2].sum())))
    qc_summary = pd.DataFrame(qrows)

    s2 = d[(d["pop"] == "main") & (d["in_s2"] == 1)]
    erows = []
    for k, sub in s2.groupby("db", sort=True):
        erows.append(dict(
            db=k, n_main_s2=len(sub),
            age_capped_at_90=int(sub["age_capped"].sum()),
            neut_abs_zero=int((sub["neut_abs"] == 0).sum()),
            lymph_abs_zero=int((sub["lymph_abs"] == 0).sum()),
            eos_abs_zero=int((sub["eos_abs"] == 0).sum()),
            diff_pct_sum_gt100=int(sub["diff_pct_sum_gt100"].sum()),
            wbc_matched_not_same_time=int(((sub["wbc_lag_min"] != 0) & sub["wbc_lag_min"].notna()).sum()),
            alb_lost_by_qc=int((is_val(sub["has_alb"], 1) & (sub["has_alb_qc"] == 0)).sum()),
            diff_lost_by_qc=int((is_val(sub["has_diff"], 1) & (sub["has_diff_qc"] == 0)).sum())))
    edge_counts = pd.DataFrame(erows)

    trows = []
    for v in ["bun_offset_h", "ph_offset_h", "alb_offset_h", "diff_offset_h", "paco2_offset_h",
              "lactate_offset_h", "creatinine_offset_h", "hgb_offset_h", "plt_offset_h", "wbc_offset_h"]:
        for k in sorted(ana["db"].unique()):
            x = ana.loc[ana["db"] == k, v].dropna()
            q = x.quantile([.25, .5, .75]) if len(x) else pd.Series([np.nan] * 3, index=[.25, .5, .75])
            trows.append(dict(db=k, var=v.replace("_offset_h", ""), n=len(x),
                              median_h=round(q[.5], 1), q1_h=round(q[.25], 1), q3_h=round(q[.75], 1),
                              pct_before_icu=round(100 * (x < 0).mean(), 1) if len(x) else np.nan,
                              pct_within_24h=round(100 * (x.abs() <= 24).mean(), 1) if len(x) else np.nan))
    lab_timing = pd.DataFrame(trows)

    def pct(mask):
        return round(100 * float(mask.mean()), 1) if len(mask) else np.nan

    drows = []
    for k, sub in ana.groupby("db", sort=True):
        drows.append(dict(
            db=k, n=len(sub),
            gcs_unable_pct=pct(is_val(sub["gcs_unable_d1"], 1)),
            gcs_from_apache_pct=pct(is_val(sub["gcs_source"], "apacheapsvar")),
            gcs_min_median=sub["gcs_min_d1"].median(),
            gcs_first_median=sub["gcs_first_d1"].median(),
            gcs_valid_median=sub["gcs_valid_first_d1"].median(),
            gcs_valid_missing_pct=pct(sub["gcs_valid_first_d1"].isna()),
            gcs_motor_median=sub["gcs_motor_first_d1"].median(),
            gcs_motor_missing_pct=pct(sub["gcs_motor_first_d1"].isna()),
            eos_low_pct=pct(is_val(sub["eos_low"], 1)),
            pneumonia_24h_pct=pct(is_val(sub["pneumonia_24h"], 1)),
            af_any_pct=pct(is_val(sub["af_any"], 1)),
            af_24h_pct=pct(is_val(sub["af_24h"], 1)),
            af_any_not_24h_n=int((is_val(sub["af_any"], 1) & is_val(sub["af_24h"], 0)).sum()),
            steroid_before_pct=pct(is_val(sub["steroid_before_diff"], 1)),
            steroid_src_emar_pct=pct(sub["steroid_source"].fillna("").str.startswith("emar")),
            steroid_ed_included_pct=pct(sub["steroid_source"].fillna("").str.contains("+ED", regex=False)),
            vent_d1_pct=pct(is_val(sub["invasive_vent_d1"], 1)),
            pneumonia_pct=pct(is_val(sub["pneumonia_any"], 1)),
            ph_arterial_pct=pct(is_val(sub["ph_specimen"], "ART.")) if k == "MIMIC" else np.nan,
            hospice_pct=np.nan if sub["hospice"].isna().all() else pct(is_val(sub["hospice"], 1))))
    diag_tab = pd.DataFrame(drows)

    ei = d[d["db"] == "eICU"].copy()
    if len(ei):
        src = np.select([is_val(ei["dx_by_string"], 1) & is_val(ei["dx_by_icd"], 1),
                         is_val(ei["dx_by_string"], 1), is_val(ei["dx_by_icd"], 1)],
                        ["string and ICD", "string only", "ICD only"], default="")
        ei["source"] = src
        ei["is_ana"] = ((ei["in_analysis"] == 1) & (ei["pop"] == "main")).astype(int)
        dx_source = ei.groupby("source").agg(n_all=("source", "size"),
                                            n_analysis=("is_ana", "sum")).reset_index()
    else:
        dx_source = pd.DataFrame([dict(note="no eICU data in this run")])

    # Additional summary tables (counts, proportions or quantiles only; not by outcome, except cluster_counts)
    # (a) Patients and deaths per IECV cluster (reported after the final snapshot, before modelling)
    crow = []
    for c_id, sub in ana.groupby("cluster", sort=True):
        crow.append(dict(cluster=int(c_id), label=CLUSTER_LABELS.get(int(c_id), ""),
                         db=sub["db"].iloc[0], n=len(sub), deaths=int(is_val(sub["death"], 1).sum()),
                         pooled_in_meta=("no (deaths < 10)" if int(is_val(sub["death"], 1).sum()) < 10 else "yes")))
    cluster_counts = pd.DataFrame(crow)

    # (b) Missing eosinophils: patients with neutrophils and lymphocytes but no eosinophils
    erow2 = []
    for k, sub in ana.groupby("db", sort=True):
        eos_na = sub["eos_abs"].isna()
        erow2.append(dict(db=k, n_analysis=len(sub), eos_abs_missing=int(eos_na.sum()),
                          eos_pct_also_missing=int((eos_na & sub["eos_pct"].isna()).sum()),
                          eos_missing_but_wbc_present=int((eos_na & sub["wbc_diff"].notna()).sum()),
                          require_eos_param=str(REQUIRE_EOS)))
    eos_check = pd.DataFrame(erow2)

    # (c) Overall distribution of the additional variables (distribution and missingness only, not by outcome)
    nrows = []
    for k, sub in ana.groupby("db", sort=True):
        for v in ["pneumonia_24h", "vasopressor_d1", "limit_tx_24h", "limit_tx_s18", "eos_low"]:
            x = sub[v]
            nrows.append(dict(db=k, var=v, n=len(sub), n_missing=int(x.isna().sum()),
                              pct_missing=pct(x.isna()), pct_yes_among_known=pct(is_val(x[x.notna()], 1)),
                              median=np.nan, q1=np.nan, q3=np.nan))
        for v in ["glucose", "mono_abs", "decaf_mod"]:
            x = sub[v].astype(float)
            q = x.quantile([.25, .5, .75]) if x.notna().any() else pd.Series([np.nan] * 3, index=[.25, .5, .75])
            nrows.append(dict(db=k, var=v, n=len(sub), n_missing=int(x.isna().sum()),
                              pct_missing=pct(x.isna()), pct_yes_among_known=np.nan,
                              median=round(q[.5], 2), q1=round(q[.25], 2), q3=round(q[.75], 2)))
    new_var_dist = pd.DataFrame(nrows)

    # (d) eICU hospital-level data availability (for hospitals without data, corticosteroids/vasopressors/treatment limitation are missing, not 0)
    arow = []
    for k, sub in ana.groupby("db", sort=True):
        arow.append(dict(db=k, n=len(sub),
                         pct_in_hosp_with_medication=pct(is_val(sub["med_data_avail"], 1)),
                         pct_in_hosp_with_infusiondrug=pct(is_val(sub["infusion_data_avail"], 1)),
                         pct_in_hosp_with_careplan=pct(is_val(sub["careplan_data_avail"], 1)),
                         steroid_before_diff_missing_pct=pct(sub["steroid_before_diff"].isna()),
                         steroid_before_pct_among_known=pct(is_val(sub["steroid_before_diff"][sub["steroid_before_diff"].notna()], 1))))
    data_avail = pd.DataFrame(arow)

    # 8. Local only: manual audit and unit checks
    audit_vars = ["db", "stay_id", "hadm_id", "subject_id", "age_raw", "bun", GCS_VAR,
                  "hr_first_d1", "ph", AF_VAR, "invasive_vent_d1", "alb_gdl", "neut_abs", "lymph_abs",
                  "eos_abs", "death"]
    audit = pd.concat([g.sample(n=min(20, len(g)), random_state=SEED)
                       for _, g in ana.groupby("db", sort=True)] or [ana.head(0)])[audit_vars].copy()
    for c in ["audit result (match/mismatch)", "mismatched variables", "auditor"]:
        audit[c] = ""
    alb_check = pd.concat([g.sample(n=min(10, len(g)), random_state=SEED)
                           for _, g in ana.groupby("db", sort=True)] or [ana.head(0)])[
        ["db", "stay_id", "alb_gdl", "alb_gl"]].copy()
    for c in ["value in source table (g/dL)", "match"]:
        alb_check[c] = ""

    dictionary = pd.DataFrame([
        ("db / pop", "database; population (main = 491.21/J44.1; related = 491.22/J44.0 only, MIMIC only)"),
        ("stay_id / hadm_id / subject_id", "MIMIC: stay_id/hadm_id/subject_id; eICU: patientunitstayid/patienthealthsystemstayid/uniquepid"),
        ("incl_* / excl_*", "inclusion and exclusion flags (0/1); see SQL 01_cohort"),
        ("has_alb / has_diff", "albumin / differential present in W4 in SQL (before QC)"),
        ("has_alb_qc / has_diff_qc", "albumin / differential still present after QC"),
        ("in_analysis", "marker subset (all inclusion criteria, no exclusions, albumin and differential after QC)"),
        ("in_s2", "analysis population (albumin and differential not required)"),
        ("death", "in-hospital death (primary outcome)"),
        ("hospice / death_s3", "discharge to hospice; outcome including hospice (empty if the database has no such category)"),
        ("age_raw / age", "raw age / age capped at 90 years"),
        ("*_offset_h", "sampling time minus ICU admission time (hours); negative = before ICU admission"),
        ("gcs_min_d1 / gcs_first_d1", "lowest / first GCS in D1 (including intubated records; for comparison only)"),
        ("gcs_valid_first_d1", "first valid GCS total score in D1 (excluding records unable to be scored because of intubation/sedation); used for BAP-65 and a sensitivity analysis"),
        ("gcs_motor_first_d1", "first GCS motor score (1–6) in D1, core model; column gcs is the variable chosen by GCS_VAR"),
        ("gcs_unable_d1", "any record unable to be scored (sedation/intubation) in D1"),
        ("af_any / af_24h", "atrial fibrillation: any diagnosis during the stay / available within 24 h of ICU admission; column af is the variable chosen by AF_VAR"),
        ("alb_gdl / alb_gl", "albumin in original unit g/dL / converted to g/L"),
        ("neut_abs / lymph_abs / eos_abs", "neutrophil / lymphocyte / eosinophil counts from the same draw (x10^9/L)"),
        ("wbc_lag_min", "time between the eICU differential and the matched WBC (minutes); 0 in MIMIC"),
        ("ln_n / ln_l", "ln(N), ln(L); N or L = 0 replaced by ZERO_FLOOR first"),
        ("ln_e", "ln(E + 0.01)"),
        ("nlr / ln_nlr / pni", "NLR, ln(NLR), PNI = Alb (g/L) + 5 x L"),
        ("sii / ln_sii / npar", "SII, ln(SII), NPAR = neutrophil % / albumin (g/dL)"),
        ("bap65", "BAP-65 class (1–5 = I–V); altered mental status = gcs_valid_first_d1 < 14; complete-case convenience only, recomputed after imputation in the analysis"),
        ("*_s4", "laboratory values restricted to within 24 h before or after ICU admission, otherwise empty (sensitivity analysis)"),
        ("race4", "race in four groups: White/Black/Other/Unknown"),
        ("steroid_before_diff", "systemic corticosteroid given before the differential was sampled; steroid_source = data source; "
                                "empty when the eICU hospital has no medication data"),
        ("eos_low", f"eosinopenia: eos_abs < {EOS_LOW_CUT} x10^9/L (binary)"),
        ("decaf_mod", "modified DECAF (0–4): eos_low + pneumonia_24h + pH < 7.30 + af_24h; complete-case convenience only"),
        ("cluster", "IECV cluster: 1 MIMIC ICD-9 era, 2 MIMIC ICD-10 era, 3–6 eICU Midwest/Northeast/South/West, 7 eICU region missing"),
        ("pneumonia_24h", "pneumonia available within 24 h of ICU admission: MIMIC from emergency department diagnoses (empty without an emergency department record), eICU from diagnoses within 24 h or the APACHE admission diagnosis"),
        ("vasopressor_d1", "any vasopressor in D1 (norepinephrine, epinephrine, dopamine, phenylephrine, vasopressin); empty when the eICU hospital has no infusiondrug data"),
        ("limit_tx_24h", "DNR/DNI/comfort care within 24 h of ICU admission: 1 yes, 0 records without limitation, empty = no record"),
        ("limit_tx_s18", "exclusion flag for the treatment-limitation sensitivity analysis: 1 if limit_tx_24h == 1, otherwise 0 (including no record)"),
        ("glucose / mono_abs", "glucose (mg/dL, closest to ICU admission in W4); monocyte count (same differential)"),
        ("plr / siri / glr / lar", "indices: PLT/L, N x M/L, glucose (mmol/L)/L, lactate/albumin (g/dL)"),
        ("*_data_avail", "whether the eICU hospital contributed medication / infusiondrug / careplangeneral data; always 1 in MIMIC"),
    ], columns=["variable", "description"])

    # 9. Write outputs
    params = (f"Parameters: GCS_VAR = {GCS_VAR}; BAP_GCS_VAR = {BAP_GCS_VAR}; AF_VAR = {AF_VAR}; "
              f"ZERO_FLOOR = {ZERO_FLOOR}; EOS_LOW_CUT = {EOS_LOW_CUT}; REQUIRE_EOS = {REQUIRE_EOS}")
    readme_local = pd.DataFrame({"note": [
        f"Snapshot tag: {tag}; created: {datetime.now():%Y-%m-%d %H:%M}",
        "This file contains patient-level data: keep it locally; do not upload it to any online service (including online AI tools) or share it.",
        "clean_all: all extracted rows (with QC and derived variables); filter with pop == 'main' and in_s2 == 1 (analysis population) or in_analysis == 1 (marker subset).",
        "audit_sample: 20 random patients per database, checked against the source tables for 10 predictors and the outcome before modelling.",
        "alb_unit_check: 10 random albumin values per database, checked against the source tables for unit conversion.",
        params]})
    readme_share = pd.DataFrame({"note": [
        f"Snapshot tag: {tag}; aggregate counts and proportions only, no patient-level records.",
        "flowchart: patients at each inclusion/exclusion step by database and population; events_scenario: number of deaths and the development/validation strategy.",
        "repro_check: comparison of MIMIC counts with an earlier extraction (690 patients, 120 deaths).",
        "missing_rates: missingness of each variable and the suggested number of imputations m; qc_set_missing: values set to missing by QC.",
        "lab_timing: median and IQR of laboratory sampling times relative to ICU admission (hours).",
        "diagnostics: distribution of GCS, atrial fibrillation, corticosteroid and ventilation definitions in both databases; edge_counts: counts of edge cases.",
        "cluster_counts: patients and deaths per IECV cluster; eos_check: missing eosinophil check.",
        "new_var_dist: overall distribution of the additional variables (not by outcome); data_avail: eICU hospital-level data availability.",
        params]})

    clean_csv = clean_dir / f"aecopd_clean_{tag}.csv"
    d.to_csv(clean_csv, index=False, na_rep="", encoding="utf-8")
    write_book(out_dir / f"aecopd_clean_{tag}_LOCAL_ONLY.xlsx",
               {"README": readme_local, "clean_all": d, "dictionary": dictionary,
                "audit_sample": audit, "alb_unit_check": alb_check})
    write_book(out_dir / f"qc_summary_{tag}_SHAREABLE.xlsx",
               {"README": readme_share, "flowchart": flowchart, "events_scenario": events_scn,
                "repro_check": repro, "missing_rates": missing_tab, "qc_set_missing": qc_summary,
                "lab_timing": lab_timing, "diagnostics": diag_tab, "edge_counts": edge_counts,
                "eicu_dx_source": dx_source, "cluster_counts": cluster_counts, "eos_check": eos_check,
                "new_var_dist": new_var_dist, "data_avail": data_avail})

    # 10. Console summary
    print("\n==== Flow chart (main population) ====")
    print(flowchart[flowchart["pop"] == "main"][["db", "step", "n", "excluded", "deaths"]]
          .to_string(index=False))
    print("\n==== Deaths and development/validation strategy ====")
    print(events.to_string(index=False))
    print("Development and validation strategy:", scenario)
    print("\n==== IECV clusters ====")
    print(cluster_counts.to_string(index=False) if len(cluster_counts) else "(none)")
    print("Suggested number of imputations m =", m_imp)
    print(f"\nOutputs written to {clean_dir} and {out_dir}")
    return d


if __name__ == "__main__":
    tag_arg = sys.argv[1] if len(sys.argv) > 1 else "trial"
    if len(sys.argv) > 2:
        root_arg = sys.argv[2]
    else:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from run_all import OUTPUT_DIR as root_arg
    run_clean(root_arg, tag_arg)
