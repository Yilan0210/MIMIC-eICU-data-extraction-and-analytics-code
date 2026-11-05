-- =====================================================================
-- MIMIC-IV 3.1 | 02_analytic.sql
-- Requires: 01_cohort.sql; mimiciv_derived built with mimic-code
-- Creates: aecopd.m_lab_w4, m_diff_w4, m_gcs_d1, m_vital_d1, m_misc, m_analytic
-- Windows: W4 = [LEAST(edregtime, admittime), intime + 24h]; D1 = [intime, intime + 24h]
-- Note: emergency department laboratory results often have no hadm_id in the derived tables, so laboratory values are matched by subject_id and time window
-- With repeated measurements, the value closest to intime (absolute difference) is used; ties go to the earlier value
-- =====================================================================

-- ---------------------------------------------------------------------
-- 1. W4 laboratory values (except the differential); arterial specimens preferred for pH and PaCO2 (specimen = 'ART.')
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.m_lab_w4;
CREATE TABLE aecopd.m_lab_w4 AS
WITH long AS (
  SELECT c.stay_id, x.var, x.val, ch.charttime, 1 AS pri, NULL::text AS spec
  FROM aecopd.m_cohort c
  JOIN mimiciv_derived.chemistry ch
    ON ch.subject_id = c.subject_id AND ch.charttime BETWEEN c.arrival AND c.intime + INTERVAL '24 hours'
  CROSS JOIN LATERAL (VALUES ('bun', ch.bun::numeric), ('albumin', ch.albumin::numeric),
                             ('creatinine', ch.creatinine::numeric),
                             ('glucose', ch.glucose::numeric)) x(var, val)   -- glucose (used for the glucose-to-lymphocyte ratio)
  WHERE x.val IS NOT NULL
  UNION ALL
  SELECT c.stay_id, x.var, x.val, cb.charttime, 1, NULL
  FROM aecopd.m_cohort c
  JOIN mimiciv_derived.complete_blood_count cb
    ON cb.subject_id = c.subject_id AND cb.charttime BETWEEN c.arrival AND c.intime + INTERVAL '24 hours'
  CROSS JOIN LATERAL (VALUES ('hgb', cb.hemoglobin::numeric), ('plt', cb.platelet::numeric),
                             ('wbc', cb.wbc::numeric)) x(var, val)
  WHERE x.val IS NOT NULL
  UNION ALL
  SELECT c.stay_id, x.var, x.val, b.charttime,
         CASE WHEN x.var IN ('ph','paco2') AND b.specimen = 'ART.' THEN 0 ELSE 1 END, b.specimen::text
  FROM aecopd.m_cohort c
  JOIN mimiciv_derived.bg b
    ON b.subject_id = c.subject_id AND b.charttime BETWEEN c.arrival AND c.intime + INTERVAL '24 hours'
  CROSS JOIN LATERAL (VALUES ('ph', b.ph::numeric), ('paco2', b.pco2::numeric),
                             ('lactate', b.lactate::numeric)) x(var, val)
  WHERE x.val IS NOT NULL
)
SELECT DISTINCT ON (l.stay_id, l.var)
       l.stay_id, l.var, l.val, l.spec,
       ROUND((EXTRACT(EPOCH FROM (l.charttime - c.intime)) / 3600.0)::numeric, 3) AS offset_h
FROM long l
JOIN aecopd.m_cohort c USING (stay_id)
ORDER BY l.stay_id, l.var, l.pri, ABS(EXTRACT(EPOCH FROM (l.charttime - c.intime))), l.charttime;

-- ---------------------------------------------------------------------
-- 2. Differential leukocyte count: one row (one specimen) of blood_differential
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.m_diff_w4;
CREATE TABLE aecopd.m_diff_w4 AS
SELECT DISTINCT ON (c.stay_id)
       c.stay_id, bd.charttime AS diff_time,
       ROUND((EXTRACT(EPOCH FROM (bd.charttime - c.intime)) / 3600.0)::numeric, 3) AS diff_offset_h,
       bd.neutrophils_abs AS neut_abs, bd.lymphocytes_abs AS lymph_abs, bd.eosinophils_abs AS eos_abs,
       bd.neutrophils     AS neut_pct, bd.lymphocytes     AS lymph_pct, bd.eosinophils     AS eos_pct,
       bd.wbc             AS wbc_diff,
       bd.monocytes_abs   AS mono_abs          -- monocytes (used for SIRI)
FROM aecopd.m_cohort c
JOIN mimiciv_derived.blood_differential bd
  ON bd.subject_id = c.subject_id AND bd.charttime BETWEEN c.arrival AND c.intime + INTERVAL '24 hours'
WHERE bd.neutrophils_abs IS NOT NULL AND bd.lymphocytes_abs IS NOT NULL
ORDER BY c.stay_id, ABS(EXTRACT(EPOCH FROM (bd.charttime - c.intime))), bd.charttime;

-- ---------------------------------------------------------------------
-- 3. GCS (D1): computed from mimiciv_derived.gcs over [intime, intime + 24h]
--    (first_day_gcs starts 6 h before intime, which does not match D1, so it is not used)
--    The derived table records the total score of intubated patients (verbal score 'ETT') as 15 and sets gcs_unable = 1
--    gcs_min_d1 / gcs_first_d1: kept for comparison (include these intubated records)
--    gcs_valid_first_d1: first valid total score (3–15) in D1, excluding intubated records with gcs_unable = 1 (used for BAP-65)
--    gcs_motor_first_d1 (core model): first motor score (1–6) in D1; can be scored in intubated patients
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.m_gcs_d1;
CREATE TABLE aecopd.m_gcs_d1 AS
SELECT c.stay_id,
       MIN(g.gcs)                                                                AS gcs_min_d1,
       (ARRAY_AGG(g.gcs ORDER BY g.charttime) FILTER (WHERE g.gcs IS NOT NULL))[1] AS gcs_first_d1,
       (ARRAY_AGG(g.gcs ORDER BY g.charttime)
          FILTER (WHERE g.gcs BETWEEN 3 AND 15 AND COALESCE(g.gcs_unable, 0) = 0))[1] AS gcs_valid_first_d1,
       (ARRAY_AGG(g.gcs_motor ORDER BY g.charttime)
          FILTER (WHERE g.gcs_motor BETWEEN 1 AND 6))[1]                         AS gcs_motor_first_d1,
       BOOL_OR(g.gcs_unable = 1)                                                 AS gcs_unable_d1
FROM aecopd.m_cohort c
JOIN mimiciv_derived.gcs g
  ON g.stay_id = c.stay_id AND g.charttime BETWEEN c.intime AND c.intime + INTERVAL '24 hours'
GROUP BY c.stay_id;

-- ---------------------------------------------------------------------
-- 4. Vital signs (first plausible value in D1; out-of-range values are skipped and the next value is used)
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.m_vital_d1;
CREATE TABLE aecopd.m_vital_d1 AS
WITH rng(var, lo, hi) AS (
  VALUES ('hr', 20, 250), ('rr', 2, 80), ('spo2', 20, 100), ('temp', 25, 45), ('map', 20, 250)
),
long AS (
  SELECT c.stay_id, v.charttime, x.var, x.val
  FROM aecopd.m_cohort c
  JOIN mimiciv_derived.vitalsign v
    ON v.stay_id = c.stay_id AND v.charttime BETWEEN c.intime AND c.intime + INTERVAL '24 hours'
  CROSS JOIN LATERAL (VALUES ('hr', v.heart_rate::numeric), ('rr', v.resp_rate::numeric),
                             ('spo2', v.spo2::numeric), ('temp', v.temperature::numeric),
                             ('map', v.mbp::numeric)) x(var, val)
  WHERE x.val IS NOT NULL
),
first_val AS (
  SELECT DISTINCT ON (l.stay_id, l.var) l.stay_id, l.var, l.val
  FROM long l
  JOIN rng r ON r.var = l.var AND l.val BETWEEN r.lo AND r.hi
  ORDER BY l.stay_id, l.var, l.charttime
)
SELECT stay_id,
       MAX(val) FILTER (WHERE var = 'hr')   AS hr_first_d1,
       MAX(val) FILTER (WHERE var = 'rr')   AS rr_first_d1,
       MAX(val) FILTER (WHERE var = 'spo2') AS spo2_first_d1,
       MAX(val) FILTER (WHERE var = 'temp') AS temp_first_d1,
       MAX(val) FILTER (WHERE var = 'map')  AS map_first_d1
FROM first_val GROUP BY stay_id;

-- ---------------------------------------------------------------------
-- 5. Other: atrial fibrillation, pneumonia, invasive ventilation on D1, systemic corticosteroids before sampling, SAPS II, SOFA
--    Atrial fibrillation: af_any = any diagnosis during the admission;
--          af_24h = rhythm charted as AF in D1 (itemid 220048) or an AF code in a previous admission (used in the model)
--    Codes: 42731 corresponds to I480/I481x/I482x/I4891 (atrial flutter excluded, consistent with ICD-9)
--    Corticosteroids: administration time from emar when the admission has emar records, otherwise the start time in prescriptions
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.m_misc;
CREATE TABLE aecopd.m_misc AS
WITH dx AS (
  SELECT c.stay_id,
    BOOL_OR((d.icd_version = 9  AND TRIM(d.icd_code) = '42731')
         OR (d.icd_version = 10 AND TRIM(d.icd_code) ~ '^I48(0|1|2|91)'))     AS af_any,
    BOOL_OR((d.icd_version = 9  AND TRIM(d.icd_code) ~ '^48[0-6]')
         OR (d.icd_version = 10 AND TRIM(d.icd_code) ~ '^J1[2-8]'))           AS pneumonia_any
  FROM aecopd.m_cohort c JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id = c.hadm_id
  GROUP BY c.stay_id
),
af_prior AS (
  SELECT DISTINCT c.stay_id
  FROM aecopd.m_cohort c
  JOIN mimiciv_hosp.admissions a ON a.subject_id = c.subject_id AND a.admittime < c.admittime
  JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id = a.hadm_id
  WHERE (d.icd_version = 9 AND TRIM(d.icd_code) = '42731')
     OR (d.icd_version = 10 AND TRIM(d.icd_code) ~ '^I48(0|1|2|91)')
),
af_rhythm AS (
  SELECT DISTINCT c.stay_id
  FROM aecopd.m_cohort c
  JOIN mimiciv_icu.chartevents ce
    ON ce.stay_id = c.stay_id AND ce.itemid = 220048
   AND ce.charttime BETWEEN c.intime AND c.intime + INTERVAL '24 hours'
  WHERE ce.value ~* '^AF'
),
vent AS (
  SELECT c.stay_id,
         BOOL_OR(v.ventilation_status = 'InvasiveVent'
                 AND v.starttime <= c.intime + INTERVAL '24 hours'
                 AND v.endtime   >= c.intime) AS invasive_vent_d1
  FROM aecopd.m_cohort c LEFT JOIN mimiciv_derived.ventilation v ON v.stay_id = c.stay_id
  GROUP BY c.stay_id
),
emar_any AS (
  SELECT DISTINCT e.hadm_id FROM mimiciv_hosp.emar e JOIN aecopd.m_cohort c USING (hadm_id)
),
emar_ster AS (
  SELECT e.hadm_id, MIN(e.charttime) AS t
  FROM mimiciv_hosp.emar e JOIN aecopd.m_cohort c USING (hadm_id)
  WHERE e.medication ~* '(hydrocortisone|solu-?cortef|methylpred|solu-?medrol|\mmedrol|prednison|prednisolon|dexameth|decadron)'
    AND e.medication !~* '(cream|oint|lotion|ophth|otic|enema|foam|supp|nasal|inhal|neb|\mdrops?\M|rectal|\meye\M|\mear\M)'
    AND (e.event_txt ~* 'administered' OR e.event_txt IN ('Started','Restarted','Confirmed'))
    AND e.event_txt !~* '\mnot\M'
  GROUP BY e.hadm_id
),
rx_ster AS (
  SELECT p.hadm_id, MIN(p.starttime) AS t
  FROM mimiciv_hosp.prescriptions p JOIN aecopd.m_cohort c USING (hadm_id)
  WHERE p.drug ~* '(hydrocortisone|solu-?cortef|methylpred|solu-?medrol|\mmedrol|prednison|prednisolon|dexameth|decadron)'
    AND p.drug !~* '(cream|oint|lotion|ophth|otic|enema|foam|supp|nasal|inhal|neb|\mdrops?\M|rectal|\meye\M|\mear\M)'
    AND COALESCE(p.route, '') !~* '^(IH|NU|NEB|INHALATION|TP|TOPICAL|OU|OD|OS|AU|AD|AS|NAS|NASAL|PR|RECTAL|BUCCAL)$'
  GROUP BY p.hadm_id
),
-- Vasopressor on D1: any vasopressor infusion overlapping [intime, intime + 24h]
vaso AS (
  SELECT DISTINCT c.stay_id
  FROM aecopd.m_cohort c
  JOIN (
    SELECT stay_id, starttime, endtime FROM mimiciv_derived.norepinephrine
    UNION ALL SELECT stay_id, starttime, endtime FROM mimiciv_derived.epinephrine
    UNION ALL SELECT stay_id, starttime, endtime FROM mimiciv_derived.dopamine
    UNION ALL SELECT stay_id, starttime, endtime FROM mimiciv_derived.phenylephrine
    UNION ALL SELECT stay_id, starttime, endtime FROM mimiciv_derived.vasopressin
  ) v ON v.stay_id = c.stay_id
     AND v.starttime <= c.intime + INTERVAL '24 hours'
     AND v.endtime   >= c.intime
),
-- Treatment limitation within 24 h of ICU admission. Code status records:
--   1 = DNR/DNI/comfort care recorded; 0 = records present but all full code; NULL = no record within 24 h
--   itemids are found by d_items label; check the spelling first with C9 in 00_check_names.sql
code_status AS (
  SELECT c.stay_id,
         BOOL_OR(ce.value ~* '(dnr|dni|do not|comfort|cmo)') AS limit_tx
  FROM aecopd.m_cohort c
  JOIN mimiciv_icu.chartevents ce
    ON ce.stay_id = c.stay_id
   AND ce.charttime BETWEEN c.intime AND c.intime + INTERVAL '24 hours'
   AND ce.itemid IN (SELECT itemid FROM mimiciv_icu.d_items WHERE label ILIKE 'code status%')  -- includes 223758 'Code Status' and 228687 'Code Status.'
  WHERE ce.value IS NOT NULL
  GROUP BY c.stay_id
)
SELECT c.stay_id,
       COALESCE(dx.af_any, false)                                               AS af_any,
       (ar.stay_id IS NOT NULL OR ap.stay_id IS NOT NULL)                       AS af_24h,
       COALESCE(dx.pneumonia_any, false)                                        AS pneumonia_any,
       COALESCE(vt.invasive_vent_d1, false)                                     AS invasive_vent_d1,
       CASE WHEN ea.hadm_id IS NOT NULL THEN 'emar' ELSE 'prescriptions' END    AS steroid_source,
       CASE WHEN df.diff_time IS NULL THEN NULL
            WHEN ea.hadm_id IS NOT NULL THEN COALESCE(es.t < df.diff_time, false)
            ELSE COALESCE(rx.t < df.diff_time, false) END                       AS steroid_before_diff,
       s.sapsii, f.sofa,
       (vs.stay_id IS NOT NULL)                                                 AS vasopressor_d1,
       cs.limit_tx                                                              AS limit_tx_24h
FROM aecopd.m_cohort c
LEFT JOIN vaso      vs     ON vs.stay_id = c.stay_id
LEFT JOIN code_status cs   ON cs.stay_id = c.stay_id
LEFT JOIN dx               ON dx.stay_id = c.stay_id
LEFT JOIN af_prior  ap     ON ap.stay_id = c.stay_id
LEFT JOIN af_rhythm ar     ON ar.stay_id = c.stay_id
LEFT JOIN vent      vt     ON vt.stay_id = c.stay_id
LEFT JOIN emar_any  ea     ON ea.hadm_id = c.hadm_id
LEFT JOIN emar_ster es     ON es.hadm_id = c.hadm_id
LEFT JOIN rx_ster   rx     ON rx.hadm_id = c.hadm_id
LEFT JOIN aecopd.m_diff_w4 df ON df.stay_id = c.stay_id
LEFT JOIN mimiciv_derived.sapsii         s ON s.stay_id = c.stay_id
LEFT JOIN mimiciv_derived.first_day_sofa f ON f.stay_id = c.stay_id;

-- ---------------------------------------------------------------------
-- 6. Analysis table: column names and order identical to aecopd.e_analytic (eICU)
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.m_analytic;
CREATE TABLE aecopd.m_analytic AS
WITH lab AS (
  SELECT stay_id,
    MAX(val)      FILTER (WHERE var = 'bun')        AS bun,        MAX(offset_h) FILTER (WHERE var = 'bun')        AS bun_off,
    MAX(val)      FILTER (WHERE var = 'ph')         AS ph,         MAX(offset_h) FILTER (WHERE var = 'ph')         AS ph_off,
    MAX(spec)     FILTER (WHERE var = 'ph')         AS ph_spec,
    MAX(val)      FILTER (WHERE var = 'albumin')    AS alb,        MAX(offset_h) FILTER (WHERE var = 'albumin')    AS alb_off,
    MAX(val)      FILTER (WHERE var = 'paco2')      AS paco2,      MAX(offset_h) FILTER (WHERE var = 'paco2')      AS paco2_off,
    MAX(val)      FILTER (WHERE var = 'lactate')    AS lactate,    MAX(offset_h) FILTER (WHERE var = 'lactate')    AS lactate_off,
    MAX(val)      FILTER (WHERE var = 'creatinine') AS creat,      MAX(offset_h) FILTER (WHERE var = 'creatinine') AS creat_off,
    MAX(val)      FILTER (WHERE var = 'hgb')        AS hgb,        MAX(offset_h) FILTER (WHERE var = 'hgb')        AS hgb_off,
    MAX(val)      FILTER (WHERE var = 'plt')        AS plt,        MAX(offset_h) FILTER (WHERE var = 'plt')        AS plt_off,
    MAX(val)      FILTER (WHERE var = 'wbc')        AS wbc,        MAX(offset_h) FILTER (WHERE var = 'wbc')        AS wbc_off,
    MAX(val)      FILTER (WHERE var = 'glucose')    AS glucose,    MAX(offset_h) FILTER (WHERE var = 'glucose')    AS glucose_off
  FROM aecopd.m_lab_w4 GROUP BY stay_id
)
SELECT
  'MIMIC'::text                                   AS db,
  c.pop::text                                     AS pop,
  c.stay_id::text                                 AS stay_id,
  c.hadm_id::text                                 AS hadm_id,
  c.subject_id::text                              AS subject_id,
  'BIDMC'::text                                   AS hospital_id,
  NULL::text                                      AS region,
  NULL::int                                       AS teaching,
  NULL::text                                      AS bed_cat,
  CASE c.icd_version WHEN 9 THEN 'ICD-9' WHEN 10 THEN 'ICD-10' END AS coding_era,
  c.anchor_year_group::text                       AS year_grp,
  -- inclusion/exclusion
  c.incl_age::int AS incl_age, c.incl_los24::int AS incl_los24, c.incl_outcome::int AS incl_outcome,
  c.excl_covid::int AS excl_covid, c.excl_heme::int AS excl_heme, c.excl_hiv::int AS excl_hiv,
  c.excl_tx::int AS excl_tx, c.excl_cirr::int AS excl_cirr,
  (lab.alb IS NOT NULL)::int                      AS has_alb,
  (df.neut_abs IS NOT NULL AND df.lymph_abs IS NOT NULL)::int AS has_diff,
  c.aecopd_primary::int                           AS aecopd_primary,
  NULL::int                                       AS dx_by_string,
  NULL::int                                       AS dx_by_icd,
  -- outcome
  c.hospital_expire_flag::int                     AS death,
  c.discharge_location::text                      AS discharge_loc,
  -- demographics
  c.age_raw::numeric                              AS age_raw,
  c.age_gt89::int                                 AS age_gt89,
  c.gender::text                                  AS sex,
  c.race::text                                    AS race_raw,
  -- core model variables
  lab.bun AS bun,                        lab.bun_off AS bun_offset_h,
  g.gcs_min_d1::numeric AS gcs_min_d1,   g.gcs_first_d1::numeric AS gcs_first_d1,
  g.gcs_valid_first_d1::numeric AS gcs_valid_first_d1, g.gcs_motor_first_d1::numeric AS gcs_motor_first_d1,
  COALESCE(g.gcs_unable_d1, false)::int  AS gcs_unable_d1,
  CASE WHEN g.gcs_valid_first_d1 IS NOT NULL THEN 'derived.gcs' END AS gcs_source,
  vt.hr_first_d1 AS hr_first_d1,
  lab.ph AS ph,                          lab.ph_off AS ph_offset_h,
  lab.ph_spec AS ph_specimen,
  ms.af_any::int AS af_any,              ms.af_24h::int AS af_24h,
  -- immune-inflammation-nutrition markers (albumin in g/dL; converted to g/L during cleaning)
  lab.alb AS alb_gdl,                    lab.alb_off AS alb_offset_h,
  df.neut_abs::numeric AS neut_abs, df.lymph_abs::numeric AS lymph_abs, df.eos_abs::numeric AS eos_abs,
  df.neut_pct::numeric AS neut_pct, df.lymph_pct::numeric AS lymph_pct, df.eos_pct::numeric AS eos_pct,
  df.wbc_diff::numeric AS wbc_diff,      df.diff_offset_h AS diff_offset_h,
  CASE WHEN df.stay_id IS NOT NULL THEN 0 END::numeric AS wbc_lag_min,
  -- auxiliary variables
  lab.paco2 AS paco2,                    lab.paco2_off   AS paco2_offset_h,
  lab.lactate AS lactate,                lab.lactate_off AS lactate_offset_h,
  lab.creat AS creatinine,               lab.creat_off   AS creatinine_offset_h,
  lab.hgb AS hgb,                        lab.hgb_off     AS hgb_offset_h,
  lab.plt AS plt,                        lab.plt_off     AS plt_offset_h,
  lab.wbc AS wbc,                        lab.wbc_off     AS wbc_offset_h,
  vt.spo2_first_d1, vt.rr_first_d1, vt.map_first_d1, vt.temp_first_d1,
  ms.invasive_vent_d1::int               AS invasive_vent_d1,
  ms.pneumonia_any::int                  AS pneumonia_any,
  ms.steroid_before_diff::int            AS steroid_before_diff,
  ms.steroid_source::text                AS steroid_source,
  -- comparator scores
  NULL::numeric                          AS apache_iva_prob,
  ms.sapsii::numeric                     AS sapsii,
  ms.sofa::numeric                       AS sofa_d1,
  -- variables for sensitivity and secondary analyses (same names and order in both databases)
  NULL::int                              AS pneumonia_24h,      -- filled from emergency department diagnoses by 03_ed_steroid.sql
  ms.vasopressor_d1::int                 AS vasopressor_d1,
  ms.limit_tx_24h::int                   AS limit_tx_24h,
  lab.glucose::numeric                   AS glucose,
  lab.glucose_off                        AS glucose_offset_h,
  df.mono_abs::numeric                   AS mono_abs,
  1::int                                 AS med_data_avail,     -- medication records are available for the whole hospital in MIMIC-IV
  1::int                                 AS infusion_data_avail,
  1::int                                 AS careplan_data_avail
FROM aecopd.m_cohort c
LEFT JOIN lab                   ON lab.stay_id = c.stay_id
LEFT JOIN aecopd.m_diff_w4 df   ON df.stay_id  = c.stay_id
LEFT JOIN aecopd.m_gcs_d1  g    ON g.stay_id   = c.stay_id
LEFT JOIN aecopd.m_vital_d1 vt  ON vt.stay_id  = c.stay_id
LEFT JOIN aecopd.m_misc    ms   ON ms.stay_id  = c.stay_id;

CREATE UNIQUE INDEX ON aecopd.m_analytic (stay_id);

-- ---------------------------------------------------------------------
-- 7. Flow chart counts (main population; for the related population change pop = 'main' to 'related')
--    Consistency check: steps 0→1→2 + exclude COVID-19 + albumin + differential should be close to an earlier extraction (690 patients, 120 deaths)
-- ---------------------------------------------------------------------
WITH f AS (SELECT * FROM aecopd.m_analytic WHERE pop = 'main')
SELECT step, n, deaths FROM (
  SELECT 0 AS ord, '0 First ICU stay with AECOPD' AS step, COUNT(*) AS n, SUM(death) AS deaths FROM f
  UNION ALL SELECT 1, '1 Age >= 18', COUNT(*), SUM(death) FROM f WHERE incl_age = 1
  UNION ALL SELECT 2, '2 ICU>=24h', COUNT(*), SUM(death) FROM f WHERE incl_age = 1 AND incl_los24 = 1
  UNION ALL SELECT 3, '3 Exclude COVID-19', COUNT(*), SUM(death) FROM f
            WHERE incl_age = 1 AND incl_los24 = 1 AND excl_covid = 0
  UNION ALL SELECT 4, '4 Consistency check: albumin + differential (before the four disease exclusions)', COUNT(*), SUM(death) FROM f
            WHERE incl_age = 1 AND incl_los24 = 1 AND excl_covid = 0 AND has_alb = 1 AND has_diff = 1
  UNION ALL SELECT 5, '5 Exclude four conditions', COUNT(*), SUM(death) FROM f
            WHERE incl_age = 1 AND incl_los24 = 1 AND excl_covid = 0
              AND excl_heme = 0 AND excl_hiv = 0 AND excl_tx = 0 AND excl_cirr = 0
  UNION ALL SELECT 6, '6 Albumin measured', COUNT(*), SUM(death) FROM f
            WHERE incl_age = 1 AND incl_los24 = 1 AND excl_covid = 0
              AND excl_heme = 0 AND excl_hiv = 0 AND excl_tx = 0 AND excl_cirr = 0 AND has_alb = 1
  UNION ALL SELECT 7, '7 Differential measured', COUNT(*), SUM(death) FROM f
            WHERE incl_age = 1 AND incl_los24 = 1 AND excl_covid = 0
              AND excl_heme = 0 AND excl_hiv = 0 AND excl_tx = 0 AND excl_cirr = 0
              AND has_alb = 1 AND has_diff = 1
) s ORDER BY ord;
