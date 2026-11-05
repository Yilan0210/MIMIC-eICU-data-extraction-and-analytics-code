# AECOPD ICU mortality model: data extraction, cleaning and analysis code

Code for cohort assembly, variable extraction, data cleaning and statistical analysis for the study

> **Transportability of a first-day mortality model for patients with acute exacerbation of COPD in intensive care: development and internal–external validation in 187 US hospitals**
> Qian Zhang, Feng Zheng, Jian He, Yan Gu, Yilan Wang. Department of Intensive Care Unit, First People's Hospital of Changshu City, Changshu Hospital Affiliated to Soochow University, Changshu, China.

The SQL scripts build the study cohorts and analysis tables in two public critical care databases, **MIMIC-IV v3.1** and **eICU-CRD v2.0**. The Python scripts run the SQL, export CSV snapshots, merge the two databases and apply quality control. The R scripts in `analysis/` develop and validate the model (see [`analysis/README.md`](analysis/README.md)).

**This repository contains code only. It contains no data.** Access to MIMIC-IV and eICU-CRD requires credentialing on [PhysioNet](https://physionet.org/) and signing the data use agreements.

## Requirements

- PostgreSQL with MIMIC-IV v3.1 and eICU-CRD v2.0 loaded (for example with the build scripts in [mimic-code](https://github.com/MIT-LCP/mimic-code) and [eicu-code](https://github.com/MIT-LCP/eicu-code))
- The MIMIC-IV derived tables (`mimiciv_derived`), built with the PostgreSQL concepts scripts in mimic-code (`mimic-iv/concepts_postgres`)
- Optional: the MIMIC-IV-ED module (`mimiciv_ed`), used for emergency department corticosteroids and pneumonia diagnoses
- Python 3.13 (tested with 3.13.16) and the packages in `python/requirements.txt`
- For the analysis: R 4.5.1 and the packages listed in `analysis/README.md`

```
pip install -r python/requirements.txt
```

## Repository structure

```
python/run_all.py              one-click run: find databases -> name checks -> build tables -> export -> clean
python/clean.py                merging, quality control and cleaning (called by run_all.py)
python/export_death90.py       export of 90-day death in MIMIC-IV (sensitivity analysis)
python/requirements.txt        Python dependencies
sql/mimic/00_check_names.sql   checks of codes, charted values and drug names before extraction
sql/mimic/01_cohort.sql        cohort: AECOPD identification, one ICU stay per patient, inclusion/exclusion flags -> aecopd.m_cohort
sql/mimic/02_analytic.sql      variable extraction and analysis table -> aecopd.m_analytic
sql/mimic/03_ed_steroid.sql    (optional, needs MIMIC-IV-ED) emergency department corticosteroids and pneumonia
sql/mimic/04_death90.sql       90-day death (used by export_death90.py)
sql/eicu/00_check_names.sql    checks of diagnosis strings, laboratory names, nurse charting labels and drug names
sql/eicu/01_cohort.sql         cohort -> aecopd.e_cohort
sql/eicu/02_analytic.sql       variable extraction and analysis table -> aecopd.e_analytic
```

All new tables are created in the schema `aecopd`; the source tables are not modified. The two analysis tables have identical column names and order; variables that are not available in a database are NULL.

## Usage

1. Open `python/run_all.py` and edit the **Settings** section: `DB_HOST`, `DB_PORT`, `DB_USER` and `OUTPUT_DIR` (the folder for all outputs). Database names and the eICU schema are found automatically. Leave `DB_PASSWORD` empty; the password is requested at run time.
2. Run
   ```
   python python/run_all.py
   ```
   On Windows: `py python\run_all.py`, or double-click the file.
3. Start with a trial run (`FINAL_SNAPSHOT = False`): file names contain `trial` and can be overwritten. After checking the flow chart counts and the name checks, set `FINAL_SNAPSHOT = True` to export a final snapshot; final snapshots carry the date and are never overwritten.
4. For the 90-day outcome (MIMIC-IV only), set the same connection settings, `OUTPUT_DIR` and the snapshot `TAG` at the top of `python/export_death90.py`, and run it after `run_all.py`.

If eICU-CRD is not found, only MIMIC-IV is processed. The large eICU tables (`vitalperiodic`, `nursecharting`, `lab`) need an index on `patientunitstayid`; a full run can take tens of minutes.

`clean.py` can also be run on its own: `python python/clean.py <tag> [<output folder>]`.

## Outputs

All outputs are written to `OUTPUT_DIR`. `<tag>` is `trial` or the date of the final snapshot.

| File | Content | Patient-level |
|---|---|---|
| `raw/*_analytic_<tag>.csv` | snapshots exported from the databases | yes |
| `raw/export_log_<tag>.csv` | row counts and MD5 checksums of the snapshots | no |
| `clean/aecopd_clean_<tag>.csv` | cleaned data used for modelling | yes |
| `output/aecopd_clean_<tag>_LOCAL_ONLY.xlsx` | cleaned data, data dictionary, manual audit sheets | yes |
| `output/sql_results_<tag>_SHAREABLE.xlsx` | name checks and SQL flow chart counts | no |
| `output/qc_summary_<tag>_SHAREABLE.xlsx` | flow chart, deaths, missingness, QC counts, sampling times, cluster counts | no |
| `aecopd_mimic_death90_<tag>_LOCAL_ONLY.csv` | 90-day death (from `export_death90.py`) | yes |

**Patient-level files (`raw/`, `clean/` and every file named `LOCAL_ONLY`) must not be uploaded to any online service, including GitHub and online AI tools, as required by the PhysioNet data use agreements.** Only files named `SHAREABLE` contain aggregate results. The `.gitignore` in this repository excludes data files.

## Main definitions

- **Cohort.** MIMIC-IV: hospital admissions with ICD-9-CM 491.21 or ICD-10-CM J44.1; eICU-CRD: diagnosis string containing "acute COPD exacerbation" or ICD code 491.21/J44.1. One ICU stay per patient (the first qualifying stay). Inclusion: age 18 years or older, ICU stay of at least 24 h, known hospital outcome. Exclusions: COVID-19, haematological malignancy, HIV infection, solid organ or bone marrow transplantation, cirrhosis.
- **Time windows.** W4 (laboratory values) runs from hospital arrival (including the emergency department) to 24 h after ICU admission; the value closest to ICU admission is used. D1 runs from ICU admission to 24 h; vital signs and the GCS motor score are the first plausible values in D1.
- **GCS.** The model uses the first motor score in D1 (`gcs_motor_first_d1`), which can be scored in intubated patients. BAP-65 uses the first valid total score below 14 (`gcs_valid_first_d1`), excluding records that cannot be scored because of intubation or sedation.
- **Atrial fibrillation** (`af_24h`): information available within 24 h of ICU admission (MIMIC-IV: rhythm charted as AF in D1 or an AF code in an earlier admission; eICU-CRD: AF diagnosis entered within 24 h or in the past history). ICD-10 codes I48.0, I48.1x, I48.2x and I48.91 match ICD-9 427.31; atrial flutter is not included.
- **Differential leukocyte count.** MIMIC-IV: absolute counts of one specimen; eICU-CRD: percentages multiplied by the white cell count of the same draw or within 120 min.
- **Pneumonia within 24 h** (modified DECAF): MIMIC-IV emergency department diagnoses (missing without a linked emergency department visit); eICU-CRD diagnoses entered within 24 h or the APACHE admission diagnosis.
- **Hospital-level data availability (eICU-CRD).** When a hospital contributed no medication, infusion or care plan data, the related variables are missing, not zero.
- **Quality control.** Values outside plausible ranges (table `QC_RULES` in `clean.py`) are set to missing; age above 89 years is set to 90.
- **Clusters for internal–external cross-validation.** 1, MIMIC-IV ICD-9 era; 2, MIMIC-IV ICD-10 era; 3–6, eICU-CRD Midwest, Northeast, South and West; 7, eICU-CRD region not recorded.

Model development and validation: see [`analysis/README.md`](analysis/README.md).

## Notes

- Never write a database password into a file that is committed to the repository.
- APACHE IVa predicted mortality values below 0.1 are often stored in scientific notation in eICU-CRD (for example `1.16E-2`); the scripts parse them.
- eICU-CRD has no specimen type for pH, so arterial samples cannot be enforced there.

## Citation

If you use this code, please cite the paper above (citation details will be added after publication).
