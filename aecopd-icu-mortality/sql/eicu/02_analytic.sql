-- =====================================================================
-- eICU-CRD 2.0 | 02_analytic.sql
-- Requires: 01_cohort.sql
-- Creates: aecopd.e_lab_w4, e_diff_w4, e_gcs_d1, e_vital_d1, e_misc, e_analytic
-- Variables, time windows and same-draw rule for the differential
-- Time unit: eICU offsets are minutes from ICU admission; the output *_offset_h columns are in hours
-- Windows: W4 = [hospitaladmitoffset, +1440]; D1 = [0, +1440]
-- With repeated measurements, the value with the smallest |offset| is used; ties go to the earlier value
-- If 00_check_names.sql shows different spellings, only the name constants in this file need to be changed
-- =====================================================================

-- ---------------------------------------------------------------------
-- 1. W4 laboratory values (except the differential): for each variable, the value closest to ICU admission
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.e_lab_w4;
CREATE TABLE aecopd.e_lab_w4 AS
WITH raw AS (
  SELECT l.patientunitstayid, l.labname, l.labresultoffset, l.labresult,
         ROW_NUMBER() OVER (PARTITION BY l.patientunitstayid, l.labname, l.labresultoffset
                            ORDER BY l.labresultrevisedoffset DESC NULLS LAST, l.labid DESC) AS rn_dup
  FROM public.lab l
  JOIN aecopd.e_cohort c USING (patientunitstayid)
  WHERE l.labname IN ('albumin','BUN','pH','paCO2','lactate','creatinine','Hgb',
                      'platelets x 1000','WBC x 1000',
                      'glucose')                                  -- glucose (used for the glucose-to-lymphocyte ratio)
    AND l.labresult IS NOT NULL
    AND l.labresultoffset BETWEEN LEAST(COALESCE(c.hospitaladmitoffset, 0), 0) AND 1440
)
SELECT DISTINCT ON (patientunitstayid, labname)
       patientunitstayid, labname, labresult, labresultoffset
FROM raw
WHERE rn_dup = 1
ORDER BY patientunitstayid, labname, ABS(labresultoffset), labresultoffset;

-- ---------------------------------------------------------------------
-- 2. Differential leukocyte count (same draw): -polys/-lymphs/-eos with the same labresultoffset;
--    WBC from the same time, otherwise the closest value within 120 min;
--    the draw closest to ICU admission with neutrophil and lymphocyte percentages and a matched WBC is used
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.e_diff_w4;
CREATE TABLE aecopd.e_diff_w4 AS
WITH raw AS (
  SELECT l.patientunitstayid, l.labname, l.labresultoffset, l.labresult,
         ROW_NUMBER() OVER (PARTITION BY l.patientunitstayid, l.labname, l.labresultoffset
                            ORDER BY l.labresultrevisedoffset DESC NULLS LAST, l.labid DESC) AS rn_dup
  FROM public.lab l
  JOIN aecopd.e_cohort c USING (patientunitstayid)
  WHERE l.labname IN ('-polys','-lymphs','-eos','-monos','WBC x 1000')   -- -monos added for SIRI
    AND l.labresult IS NOT NULL
    AND l.labresultoffset BETWEEN LEAST(COALESCE(c.hospitaladmitoffset, 0), 0) - 120 AND 1440 + 120
),
dedup AS (SELECT * FROM raw WHERE rn_dup = 1),
draws AS (
  SELECT d.patientunitstayid, d.labresultoffset,
         MAX(d.labresult) FILTER (WHERE d.labname = '-polys')  AS neut_pct,
         MAX(d.labresult) FILTER (WHERE d.labname = '-lymphs') AS lymph_pct,
         MAX(d.labresult) FILTER (WHERE d.labname = '-eos')    AS eos_pct,
         MAX(d.labresult) FILTER (WHERE d.labname = '-monos')  AS mono_pct
  FROM dedup d
  JOIN aecopd.e_cohort c USING (patientunitstayid)
  WHERE d.labname IN ('-polys','-lymphs','-eos','-monos')
    AND d.labresultoffset BETWEEN LEAST(COALESCE(c.hospitaladmitoffset, 0), 0) AND 1440
  GROUP BY d.patientunitstayid, d.labresultoffset
),
draws_wbc AS (
  SELECT dr.*, w.labresult AS wbc_diff, w.labresultoffset AS wbc_offset
  FROM draws dr
  LEFT JOIN LATERAL (
    SELECT x.labresult, x.labresultoffset
    FROM dedup x
    WHERE x.patientunitstayid = dr.patientunitstayid
      AND x.labname = 'WBC x 1000'
      AND ABS(x.labresultoffset - dr.labresultoffset) <= 120
    ORDER BY ABS(x.labresultoffset - dr.labresultoffset), x.labresultoffset
    LIMIT 1
  ) w ON TRUE
  WHERE dr.neut_pct IS NOT NULL AND dr.lymph_pct IS NOT NULL
)
SELECT DISTINCT ON (patientunitstayid)
       patientunitstayid,
       labresultoffset                 AS diff_offset,
       wbc_offset - labresultoffset    AS wbc_lag_min,
       neut_pct, lymph_pct, eos_pct, wbc_diff,
       neut_pct  * wbc_diff / 100.0    AS neut_abs,
       lymph_pct * wbc_diff / 100.0    AS lymph_abs,
       eos_pct   * wbc_diff / 100.0    AS eos_abs,
       mono_pct  * wbc_diff / 100.0    AS mono_abs
FROM draws_wbc
WHERE wbc_diff IS NOT NULL
ORDER BY patientunitstayid, ABS(labresultoffset), labresultoffset;

-- ---------------------------------------------------------------------
-- 3. GCS (D1)
--    gcs_min_d1 / gcs_first_d1: lowest / first GCS Total (kept for comparison)
--    gcs_valid_first_d1: first valid total score (3–15) in D1, excluding "unable to score (sedated)" records (used for BAP-65);
--        from GCS Total or Score (Glasgow Coma Scale); otherwise apacheapsvar (not used when meds = 1)
--    gcs_motor_first_d1 (core model): first motor score (1–6) in D1; can be scored in intubated patients
--    Limitation: eICU does not allow GCS records taken during intubation to be identified one by one (MIMIC does); see README
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.e_gcs_d1;
CREATE TABLE aecopd.e_gcs_d1 AS
WITH nc AS (
  SELECT n.patientunitstayid, n.nursingchartoffset AS t,
         n.nursingchartcelltypevalname                       AS nm,
         CASE WHEN TRIM(n.nursingchartvalue) ~ '^[0-9]+(\.[0-9]+)?$'
              THEN TRIM(n.nursingchartvalue)::numeric END    AS val,
         (n.nursingchartvalue ~* 'unable')                   AS unable
  FROM public.nursecharting n
  JOIN aecopd.e_cohort c USING (patientunitstayid)
  WHERE n.nursingchartoffset BETWEEN 0 AND 1440
    AND ((n.nursingchartcelltypevallabel = 'Glasgow coma score'
          AND n.nursingchartcelltypevalname IN ('GCS Total', 'Motor'))
      OR (n.nursingchartcelltypevallabel = 'Score (Glasgow Coma Scale)'
          AND n.nursingchartcelltypevalname = 'Value'))
),
nc_agg AS (
  SELECT patientunitstayid,
    MIN(val) FILTER (WHERE nm = 'GCS Total')                                                  AS gcs_min_nc,
    (ARRAY_AGG(val ORDER BY t) FILTER (WHERE nm = 'GCS Total' AND val IS NOT NULL))[1]        AS gcs_first_nc,
    (ARRAY_AGG(val ORDER BY t) FILTER (WHERE nm IN ('GCS Total', 'Value')
                                         AND val BETWEEN 3 AND 15))[1]                         AS gcs_valid_nc,
    (ARRAY_AGG(val ORDER BY t) FILTER (WHERE nm = 'Motor' AND val BETWEEN 1 AND 6))[1]        AS motor_nc,
    BOOL_OR(unable AND nm = 'GCS Total')                                                      AS unable_nc
  FROM nc GROUP BY patientunitstayid
),
aps AS (
  SELECT patientunitstayid,
         CASE WHEN eyes BETWEEN 1 AND 4 AND motor BETWEEN 1 AND 6 AND verbal BETWEEN 1 AND 5
              THEN eyes + motor + verbal END              AS gcs_aps,
         CASE WHEN motor BETWEEN 1 AND 6 THEN motor END   AS motor_aps,
         (meds = 1)                                       AS meds_aps
  FROM public.apacheapsvar
)
SELECT c.patientunitstayid,
       COALESCE(a.gcs_min_nc,   ap.gcs_aps)                                       AS gcs_min_d1,
       COALESCE(a.gcs_first_nc, ap.gcs_aps)                                       AS gcs_first_d1,
       COALESCE(a.gcs_valid_nc,
                CASE WHEN NOT COALESCE(ap.meds_aps, false) THEN ap.gcs_aps END)   AS gcs_valid_first_d1,
       COALESCE(a.motor_nc, ap.motor_aps)                                         AS gcs_motor_first_d1,
       CASE WHEN a.gcs_valid_nc IS NOT NULL THEN 'nursecharting'
            WHEN ap.gcs_aps IS NOT NULL AND NOT COALESCE(ap.meds_aps, false)
                 THEN 'apacheapsvar' END                                          AS gcs_source,
       (COALESCE(a.unable_nc, false) OR COALESCE(ap.meds_aps, false))             AS gcs_unable_d1
FROM aecopd.e_cohort c
LEFT JOIN nc_agg a  USING (patientunitstayid)
LEFT JOIN aps    ap USING (patientunitstayid);

-- ---------------------------------------------------------------------
-- 4. Vital signs (first plausible value in D1): vitalperiodic first, nursecharting when missing;
--    out-of-range values (mostly monitor artefacts) are skipped and the next value is used;
--    MAP: the earlier of invasive (systemicmean) and non-invasive (noninvasivemean) measurements
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.e_vital_d1;
CREATE TABLE aecopd.e_vital_d1 AS
WITH rng(var, lo, hi) AS (
  VALUES ('hr', 20, 250), ('rr', 2, 80), ('spo2', 20, 100), ('temp', 25, 45), ('map', 20, 250)
),
long AS (
  SELECT v.patientunitstayid, v.observationoffset AS t, 1 AS pri, x.var, x.val
  FROM public.vitalperiodic v
  JOIN aecopd.e_cohort c USING (patientunitstayid)
  CROSS JOIN LATERAL (VALUES ('hr',   v.heartrate::numeric),
                             ('rr',   v.respiration::numeric),
                             ('spo2', v.sao2::numeric),
                             ('temp', v.temperature::numeric),
                             ('map',  v.systemicmean::numeric)) x(var, val)
  WHERE v.observationoffset BETWEEN 0 AND 1440 AND x.val IS NOT NULL
  UNION ALL
  SELECT va.patientunitstayid, va.observationoffset, 1, 'map', va.noninvasivemean::numeric
  FROM public.vitalaperiodic va
  JOIN aecopd.e_cohort c USING (patientunitstayid)
  WHERE va.observationoffset BETWEEN 0 AND 1440 AND va.noninvasivemean IS NOT NULL
  UNION ALL
  SELECT n.patientunitstayid, n.nursingchartoffset, 2,
         CASE n.nursingchartcelltypevalname
              WHEN 'Heart Rate'       THEN 'hr'
              WHEN 'Respiratory Rate' THEN 'rr'
              WHEN 'O2 Saturation'    THEN 'spo2'
              WHEN 'Temperature (C)'  THEN 'temp' END,
         TRIM(n.nursingchartvalue)::numeric
  FROM public.nursecharting n
  JOIN aecopd.e_cohort c USING (patientunitstayid)
  WHERE n.nursingchartcelltypevalname IN ('Heart Rate','Respiratory Rate','O2 Saturation','Temperature (C)')
    AND n.nursingchartoffset BETWEEN 0 AND 1440
    AND TRIM(n.nursingchartvalue) ~ '^[0-9]+(\.[0-9]+)?$'
),
first_val AS (
  SELECT DISTINCT ON (l.patientunitstayid, l.var) l.patientunitstayid, l.var, l.val
  FROM long l
  JOIN rng r ON r.var = l.var AND l.val BETWEEN r.lo AND r.hi
  ORDER BY l.patientunitstayid, l.var, l.pri, l.t
)
SELECT patientunitstayid,
       MAX(val) FILTER (WHERE var = 'hr')   AS hr_first_d1,
       MAX(val) FILTER (WHERE var = 'rr')   AS rr_first_d1,
       MAX(val) FILTER (WHERE var = 'spo2') AS spo2_first_d1,
       MAX(val) FILTER (WHERE var = 'temp') AS temp_first_d1,
       MAX(val) FILTER (WHERE var = 'map')  AS map_first_d1
FROM first_val
GROUP BY patientunitstayid;

-- ---------------------------------------------------------------------
-- 5. Other: atrial fibrillation, pneumonia, invasive ventilation on D1, systemic corticosteroids before sampling, APACHE IVa
--    Two definitions of atrial fibrillation: af_any = any diagnosis during the hospitalisation;
--                  af_24h = diagnosis entered within 24 h of this ICU admission or in the past history (used in the model; avoids information from later in the stay)
--    Codes: 427.31 corresponds to I48.0/I48.1x/I48.2x/I48.91 (atrial flutter I48.3/I48.4/I48.92 excluded, consistent with ICD-9)
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.e_misc;
CREATE TABLE aecopd.e_misc AS
WITH hosp_dx AS (
  SELECT c.patientunitstayid AS index_stay, d.patientunitstayid, d.diagnosisoffset,
         d.diagnosisstring, COALESCE(d.icd9code, '') AS icd
  FROM aecopd.e_cohort c
  JOIN public.patient p ON p.patienthealthsystemstayid = c.patienthealthsystemstayid
  JOIN public.diagnosis d ON d.patientunitstayid = p.patientunitstayid
),
dx AS (
  SELECT index_stay,
    BOOL_OR(icd ~ '(^|,\s*)(427\.31|I48\.(0|1|2|91))' OR diagnosisstring ~* 'atrial fibrillation') AS af_any,
    BOOL_OR(patientunitstayid = index_stay AND diagnosisoffset <= 1440
            AND (icd ~ '(^|,\s*)(427\.31|I48\.(0|1|2|91))' OR diagnosisstring ~* 'atrial fibrillation')) AS af_dx24,
    BOOL_OR(icd ~ '(^|,\s*)(48[0-6]|J1[2-8])' OR diagnosisstring ~* 'pneumonia')                    AS pneumonia_any,
    BOOL_OR(patientunitstayid = index_stay AND diagnosisoffset <= 1440
            AND (icd ~ '(^|,\s*)(48[0-6]|J1[2-8])' OR diagnosisstring ~* 'pneumonia'))              AS pna_dx24
  FROM hosp_dx GROUP BY index_stay
),
-- Additional variables ------------------------------------------------
-- Hospital-level data availability: not all eICU hospitals contributed medication / infusiondrug / careplangeneral data.
-- When a hospital has no records in a table, the corresponding variable is set to NULL (unknown) for its patients, not to "no medication / no limitation"
med_h  AS (SELECT DISTINCT p.hospitalid FROM public.medication m      JOIN public.patient p USING (patientunitstayid)),
inf_h  AS (SELECT DISTINCT p.hospitalid FROM public.infusiondrug i    JOIN public.patient p USING (patientunitstayid)),
cpl_h  AS (SELECT DISTINCT p.hospitalid FROM public.careplangeneral g JOIN public.patient p USING (patientunitstayid)),
-- Pneumonia in the APACHE admission diagnosis (entered at ICU admission, available within 24 h)
adx AS (
  SELECT a.patientunitstayid AS stay_adx,      -- renamed to avoid a clash with USING (patientunitstayid) below
         BOOL_OR(a.admitdxpath ~* 'pneumonia' OR a.admitdxname ~* 'pneumonia') AS pna_adx
  FROM public.admissiondx a JOIN aecopd.e_cohort c USING (patientunitstayid)
  GROUP BY a.patientunitstayid
),
-- Vasopressor on D1: infusiondrug records of vasopressors with a non-zero rate within 0–1440 min
vaso AS (
  SELECT DISTINCT i.patientunitstayid AS stay_vaso
  FROM public.infusiondrug i JOIN aecopd.e_cohort c USING (patientunitstayid)
  WHERE i.infusionoffset BETWEEN 0 AND 1440
    AND i.drugname ~* '(norepinephrine|levophed|epinephrine|adrenalin|dopamine|phenylephrine|neo-?synephrine|vasopressin|pitressin)'
    AND COALESCE(i.drugrate, '') ~ '[1-9]'
),
-- Treatment limitation within 24 h of ICU admission: Care Limitation group of careplangeneral
--   1 = DNR/DNI/comfort care; 0 = records in this group but none of these limitations; NULL = no record in this group within 24 h
--   check the spelling of the values first with C11 in 00_check_names.sql
cpl AS (
  SELECT g.patientunitstayid AS stay_cpl,
         BOOL_OR(g.cplitemvalue ~* '(do not resuscitate|no cpr|\mdnr\M|no intubation|\mdni\M|comfort)') AS limit_tx
  FROM public.careplangeneral g JOIN aecopd.e_cohort c USING (patientunitstayid)
  WHERE g.cplgroup = 'Care Limitation' AND g.cplitemoffset <= 1440
  GROUP BY g.patientunitstayid
),
-- --------------------------------------------------------------------
ph AS (
  SELECT ph.patientunitstayid, BOOL_OR(ph.pasthistorypath ~* 'atrial fibrillation') AS af_hx
  FROM public.pasthistory ph JOIN aecopd.e_cohort c USING (patientunitstayid)
  GROUP BY ph.patientunitstayid
),
ster AS (
  SELECT m.patientunitstayid, MIN(m.drugstartoffset) AS first_steroid_offset
  FROM public.medication m JOIN aecopd.e_cohort c USING (patientunitstayid)
  WHERE m.drugordercancelled = 'No'
    AND m.drugname ~* '(hydrocortisone|solu-?cortef|methylpred|solu-?medrol|\mmedrol|prednison|prednisolon|dexameth|decadron)'
    AND m.drugname !~* '(cream|crea\M|oint|lotion|ophth|otic|enema|foam|supp|nasal|inhal|neb|\mdrops?\M|\mex\M)'
    AND COALESCE(m.routeadmin, '') !~* '(inhal|neb|topical|ophth|otic|\meye|\mear\M|nasal|rectal|\mpr\M|enema)'
  GROUP BY m.patientunitstayid
),
apache AS (
  -- Note: in eICU, predicted values below 0.1 are often stored in scientific notation (e.g. 1.16E-2) and must also be parsed
  SELECT patientunitstayid,
         NULLIF(CASE WHEN TRIM(predictedhospitalmortality::text) ~ '^-?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'
                     THEN TRIM(predictedhospitalmortality::text)::numeric END, -1) AS apache_iva_prob
  FROM public.apachepatientresult
  WHERE apacheversion = 'IVa'
),
vent AS (
  SELECT patientunitstayid, BOOL_OR(oobintubday1 = 1) AS invasive_vent_d1
  FROM public.apachepredvar GROUP BY patientunitstayid
)
SELECT c.patientunitstayid,
       COALESCE(dx.af_any, false)                                     AS af_any,
       COALESCE(dx.af_dx24, false) OR COALESCE(ph.af_hx, false)       AS af_24h,
       COALESCE(dx.pneumonia_any, false)                              AS pneumonia_any,
       v.invasive_vent_d1,
       -- set to NULL when the hospital has no medication data (coding these as 0 would underestimate corticosteroid exposure)
       CASE WHEN df.diff_offset IS NULL OR mh.hospitalid IS NULL THEN NULL
            ELSE COALESCE(s.first_steroid_offset < df.diff_offset, false) END AS steroid_before_diff,
       ap.apache_iva_prob,
       (COALESCE(dx.pna_dx24, false) OR COALESCE(ax.pna_adx, false))   AS pneumonia_24h,
       CASE WHEN ih.hospitalid IS NULL THEN NULL ELSE (vs.stay_vaso IS NOT NULL) END AS vasopressor_d1,
       cp.limit_tx                                                     AS limit_tx_24h,
       (mh.hospitalid IS NOT NULL)                                     AS med_data_avail,
       (ih.hospitalid IS NOT NULL)                                     AS infusion_data_avail,
       (ch.hospitalid IS NOT NULL)                                     AS careplan_data_avail
FROM aecopd.e_cohort c
LEFT JOIN med_h mh ON mh.hospitalid = c.hospitalid
LEFT JOIN inf_h ih ON ih.hospitalid = c.hospitalid
LEFT JOIN cpl_h ch ON ch.hospitalid = c.hospitalid
LEFT JOIN adx  ax  ON ax.stay_adx  = c.patientunitstayid
LEFT JOIN vaso vs  ON vs.stay_vaso = c.patientunitstayid
LEFT JOIN cpl  cp  ON cp.stay_cpl  = c.patientunitstayid
LEFT JOIN dx   ON dx.index_stay = c.patientunitstayid
LEFT JOIN ph   USING (patientunitstayid)
LEFT JOIN ster s  USING (patientunitstayid)
LEFT JOIN aecopd.e_diff_w4 df USING (patientunitstayid)
LEFT JOIN apache ap USING (patientunitstayid)
LEFT JOIN vent v USING (patientunitstayid);

-- ---------------------------------------------------------------------
-- 6. Analysis table: column names and order identical to aecopd.m_analytic (MIMIC)
--    (variables not available in this database are NULL; flag columns are 0/1)
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.e_analytic;
CREATE TABLE aecopd.e_analytic AS
WITH lab AS (
  SELECT patientunitstayid,
    MAX(labresult)       FILTER (WHERE labname = 'BUN')              AS bun,
    MAX(labresultoffset) FILTER (WHERE labname = 'BUN')              AS bun_off,
    MAX(labresult)       FILTER (WHERE labname = 'pH')               AS ph,
    MAX(labresultoffset) FILTER (WHERE labname = 'pH')               AS ph_off,
    MAX(labresult)       FILTER (WHERE labname = 'albumin')          AS alb,
    MAX(labresultoffset) FILTER (WHERE labname = 'albumin')          AS alb_off,
    MAX(labresult)       FILTER (WHERE labname = 'paCO2')            AS paco2,
    MAX(labresultoffset) FILTER (WHERE labname = 'paCO2')            AS paco2_off,
    MAX(labresult)       FILTER (WHERE labname = 'lactate')          AS lactate,
    MAX(labresultoffset) FILTER (WHERE labname = 'lactate')          AS lactate_off,
    MAX(labresult)       FILTER (WHERE labname = 'creatinine')       AS creat,
    MAX(labresultoffset) FILTER (WHERE labname = 'creatinine')       AS creat_off,
    MAX(labresult)       FILTER (WHERE labname = 'Hgb')              AS hgb,
    MAX(labresultoffset) FILTER (WHERE labname = 'Hgb')              AS hgb_off,
    MAX(labresult)       FILTER (WHERE labname = 'platelets x 1000') AS plt,
    MAX(labresultoffset) FILTER (WHERE labname = 'platelets x 1000') AS plt_off,
    MAX(labresult)       FILTER (WHERE labname = 'WBC x 1000')       AS wbc,
    MAX(labresultoffset) FILTER (WHERE labname = 'WBC x 1000')       AS wbc_off,
    MAX(labresult)       FILTER (WHERE labname = 'glucose')          AS glucose,
    MAX(labresultoffset) FILTER (WHERE labname = 'glucose')          AS glucose_off
  FROM aecopd.e_lab_w4 GROUP BY patientunitstayid
)
SELECT
  'eICU'::text                                    AS db,
  'main'::text                                    AS pop,
  c.patientunitstayid::text                       AS stay_id,
  c.patienthealthsystemstayid::text               AS hadm_id,
  c.uniquepid::text                               AS subject_id,
  c.hospitalid::text                              AS hospital_id,
  h.region::text                                  AS region,
  h.teachingstatus::int                           AS teaching,
  h.numbedscategory::text                         AS bed_cat,
  NULL::text                                      AS coding_era,
  c.hospitaldischargeyear::text                   AS year_grp,
  -- inclusion/exclusion
  c.incl_age::int AS incl_age, c.incl_los24::int AS incl_los24, c.incl_outcome::int AS incl_outcome,
  c.excl_covid::int AS excl_covid, c.excl_heme::int AS excl_heme, c.excl_hiv::int AS excl_hiv,
  c.excl_tx::int AS excl_tx, c.excl_cirr::int AS excl_cirr,
  (lab.alb IS NOT NULL)::int                      AS has_alb,
  (df.neut_abs IS NOT NULL AND df.lymph_abs IS NOT NULL)::int AS has_diff,
  c.aecopd_primary::int                           AS aecopd_primary,
  c.dx_by_string::int                             AS dx_by_string,
  c.dx_by_icd::int                                AS dx_by_icd,
  -- outcome
  CASE c.hospitaldischargestatus WHEN 'Expired' THEN 1 WHEN 'Alive' THEN 0 END AS death,
  c.hospitaldischargelocation::text               AS discharge_loc,
  -- demographics
  c.age_raw::numeric                              AS age_raw,
  c.age_gt89::int                                 AS age_gt89,
  CASE c.gender WHEN 'Male' THEN 'M' WHEN 'Female' THEN 'F' END AS sex,
  NULLIF(TRIM(c.ethnicity), '')::text             AS race_raw,
  -- core model variables
  lab.bun::numeric AS bun,               ROUND(lab.bun_off / 60.0, 3)  AS bun_offset_h,
  g.gcs_min_d1::numeric AS gcs_min_d1,   g.gcs_first_d1::numeric AS gcs_first_d1,
  g.gcs_valid_first_d1::numeric AS gcs_valid_first_d1, g.gcs_motor_first_d1::numeric AS gcs_motor_first_d1,
  g.gcs_unable_d1::int AS gcs_unable_d1, g.gcs_source::text AS gcs_source,
  vt.hr_first_d1::numeric AS hr_first_d1,
  lab.ph::numeric AS ph,                 ROUND(lab.ph_off / 60.0, 3)   AS ph_offset_h,
  NULL::text                             AS ph_specimen,
  ms.af_any::int AS af_any,              ms.af_24h::int AS af_24h,
  -- immune-inflammation-nutrition markers (albumin in g/dL; converted to g/L during cleaning)
  lab.alb::numeric AS alb_gdl,           ROUND(lab.alb_off / 60.0, 3)  AS alb_offset_h,
  ROUND(df.neut_abs, 4) AS neut_abs, ROUND(df.lymph_abs, 4) AS lymph_abs, ROUND(df.eos_abs, 4) AS eos_abs,
  df.neut_pct::numeric AS neut_pct, df.lymph_pct::numeric AS lymph_pct, df.eos_pct::numeric AS eos_pct,
  df.wbc_diff::numeric AS wbc_diff,      ROUND(df.diff_offset / 60.0, 3) AS diff_offset_h,
  df.wbc_lag_min::numeric                AS wbc_lag_min,
  -- auxiliary variables
  lab.paco2::numeric AS paco2,           ROUND(lab.paco2_off / 60.0, 3)   AS paco2_offset_h,
  lab.lactate::numeric AS lactate,       ROUND(lab.lactate_off / 60.0, 3) AS lactate_offset_h,
  lab.creat::numeric AS creatinine,      ROUND(lab.creat_off / 60.0, 3)   AS creatinine_offset_h,
  lab.hgb::numeric AS hgb,               ROUND(lab.hgb_off / 60.0, 3)     AS hgb_offset_h,
  lab.plt::numeric AS plt,               ROUND(lab.plt_off / 60.0, 3)     AS plt_offset_h,
  lab.wbc::numeric AS wbc,               ROUND(lab.wbc_off / 60.0, 3)     AS wbc_offset_h,
  vt.spo2_first_d1::numeric AS spo2_first_d1, vt.rr_first_d1::numeric AS rr_first_d1,
  vt.map_first_d1::numeric AS map_first_d1,   vt.temp_first_d1::numeric AS temp_first_d1,
  ms.invasive_vent_d1::int               AS invasive_vent_d1,
  ms.pneumonia_any::int                  AS pneumonia_any,
  ms.steroid_before_diff::int            AS steroid_before_diff,
  'medication'::text                     AS steroid_source,
  -- comparator scores
  ms.apache_iva_prob::numeric            AS apache_iva_prob,
  NULL::numeric                          AS sapsii,
  NULL::numeric                          AS sofa_d1,
  -- variables for sensitivity and secondary analyses (same names and order in both databases)
  ms.pneumonia_24h::int                  AS pneumonia_24h,
  ms.vasopressor_d1::int                 AS vasopressor_d1,
  ms.limit_tx_24h::int                   AS limit_tx_24h,
  lab.glucose::numeric                   AS glucose,
  ROUND(lab.glucose_off / 60.0, 3)       AS glucose_offset_h,
  ROUND(df.mono_abs, 4)                  AS mono_abs,
  ms.med_data_avail::int                 AS med_data_avail,
  ms.infusion_data_avail::int            AS infusion_data_avail,
  ms.careplan_data_avail::int            AS careplan_data_avail
FROM aecopd.e_cohort c
LEFT JOIN public.hospital h      ON h.hospitalid = c.hospitalid
LEFT JOIN lab                    USING (patientunitstayid)
LEFT JOIN aecopd.e_diff_w4 df    USING (patientunitstayid)
LEFT JOIN aecopd.e_gcs_d1  g     USING (patientunitstayid)
LEFT JOIN aecopd.e_vital_d1 vt   USING (patientunitstayid)
LEFT JOIN aecopd.e_misc    ms    USING (patientunitstayid);

CREATE UNIQUE INDEX ON aecopd.e_analytic (stay_id);

-- ---------------------------------------------------------------------
-- 7. Flow chart counts (aggregate numbers only; the cleaning script recomputes them after quality control)
-- ---------------------------------------------------------------------
WITH f AS (SELECT * FROM aecopd.e_analytic)
SELECT step, n, deaths FROM (
  SELECT 0 AS ord, '0 First ICU stay with AECOPD' AS step, COUNT(*) AS n, SUM(death) AS deaths FROM f
  UNION ALL SELECT 1, '1 Age >= 18', COUNT(*), SUM(death) FROM f WHERE incl_age = 1
  UNION ALL SELECT 2, '2 ICU>=24h', COUNT(*), SUM(death) FROM f WHERE incl_age = 1 AND incl_los24 = 1
  UNION ALL SELECT 3, '3 Outcome known', COUNT(*), SUM(death) FROM f WHERE incl_age = 1 AND incl_los24 = 1 AND incl_outcome = 1
  UNION ALL SELECT 4, '4 Exclude four conditions', COUNT(*), SUM(death) FROM f
            WHERE incl_age = 1 AND incl_los24 = 1 AND incl_outcome = 1
              AND excl_heme = 0 AND excl_hiv = 0 AND excl_tx = 0 AND excl_cirr = 0
  UNION ALL SELECT 5, '5 Albumin measured', COUNT(*), SUM(death) FROM f
            WHERE incl_age = 1 AND incl_los24 = 1 AND incl_outcome = 1
              AND excl_heme = 0 AND excl_hiv = 0 AND excl_tx = 0 AND excl_cirr = 0 AND has_alb = 1
  UNION ALL SELECT 6, '6 Differential measured', COUNT(*), SUM(death) FROM f
            WHERE incl_age = 1 AND incl_los24 = 1 AND incl_outcome = 1
              AND excl_heme = 0 AND excl_hiv = 0 AND excl_tx = 0 AND excl_cirr = 0
              AND has_alb = 1 AND has_diff = 1
) s ORDER BY ord;
