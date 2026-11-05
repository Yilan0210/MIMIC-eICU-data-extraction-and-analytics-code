-- =====================================================================
-- MIMIC-IV 3.1 | 03_ed_steroid.sql (optional; also adds emergency department pneumonia diagnoses)
-- Runs only if the MIMIC-IV-ED module (mimiciv_ed.edstays, pyxis, diagnosis) is installed;
-- otherwise run_all.py skips this file.
-- Purpose:
--   1. Adds systemic corticosteroids given in the emergency department (pyxis dispensing time as a proxy for administration),
--      updates steroid_before_diff in aecopd.m_analytic and appends "+ED" to steroid_source
--   2. pneumonia_24h = pneumonia code among emergency department diagnoses (available within 24 h of ICU admission; used in modified DECAF)
--      0/1 for stays with a linked emergency department visit; NULL (imputed later) for direct admissions without one
-- =====================================================================

DROP TABLE IF EXISTS aecopd.m_ed_steroid;
CREATE TABLE aecopd.m_ed_steroid AS
SELECT c.stay_id, MIN(p.charttime) AS first_ed_steroid_time
FROM aecopd.m_cohort c
JOIN mimiciv_ed.edstays e ON e.hadm_id = c.hadm_id
JOIN mimiciv_ed.pyxis   p ON p.stay_id = e.stay_id
WHERE p.name ~* '(hydrocortisone|solu-?cortef|methylpred|solu-?medrol|\mmedrol|prednison|prednisolon|dexameth|decadron)'
  AND p.name !~* '(cream|oint|lotion|ophth|otic|enema|foam|supp|nasal|inhal|neb|\mdrops?\M|rectal|\meye\M|\mear\M)'
GROUP BY c.stay_id;

UPDATE aecopd.m_analytic
SET steroid_source = steroid_source || '+ED'
WHERE steroid_source NOT LIKE '%+ED';

UPDATE aecopd.m_analytic a
SET steroid_before_diff = 1
FROM aecopd.m_ed_steroid s
JOIN aecopd.m_diff_w4 d ON d.stay_id = s.stay_id
WHERE a.stay_id = s.stay_id::text
  AND s.first_ed_steroid_time < d.diff_time;

SELECT COUNT(*) AS n_with_ed_steroid_before_diff
FROM aecopd.m_ed_steroid s JOIN aecopd.m_diff_w4 d ON d.stay_id = s.stay_id
WHERE s.first_ed_steroid_time < d.diff_time;

-- Pneumonia in emergency department diagnoses (ICD-9 480–486, ICD-10 J12–J18, same codes as pneumonia_any)
DROP TABLE IF EXISTS aecopd.m_ed_pneumonia;
CREATE TABLE aecopd.m_ed_pneumonia AS
SELECT c.stay_id,
       BOOL_OR((dg.icd_version = 9  AND TRIM(dg.icd_code) ~ '^48[0-6]')
            OR (dg.icd_version = 10 AND TRIM(dg.icd_code) ~ '^J1[2-8]')) AS pneumonia_ed
FROM aecopd.m_cohort c
JOIN mimiciv_ed.edstays   e  ON e.hadm_id = c.hadm_id
LEFT JOIN mimiciv_ed.diagnosis dg ON dg.stay_id = e.stay_id
GROUP BY c.stay_id;

UPDATE aecopd.m_analytic a
SET pneumonia_24h = COALESCE(p.pneumonia_ed, false)::int
FROM aecopd.m_ed_pneumonia p
WHERE a.stay_id = p.stay_id::text;

SELECT pneumonia_24h, COUNT(*) AS n_stays
FROM aecopd.m_analytic WHERE pop = 'main' GROUP BY 1 ORDER BY 1;

-- Check: are stays without a linked emergency department record (pneumonia_24h NULL) concentrated in some years? (MIMIC-IV-ED covers a limited period)
SELECT year_grp, COUNT(*) AS n_stays, COUNT(*) FILTER (WHERE pneumonia_24h IS NULL) AS n_no_ed_link
FROM aecopd.m_analytic WHERE pop = 'main' GROUP BY 1 ORDER BY 1;

-- Check: pneumonia codes matched in emergency department diagnoses (database-wide counts; strings and counts only)
SELECT icd_version, TRIM(icd_code) AS icd_code, icd_title, COUNT(DISTINCT stay_id) AS n_ed_stays
FROM mimiciv_ed.diagnosis
WHERE (icd_version = 9 AND TRIM(icd_code) ~ '^48[0-6]') OR (icd_version = 10 AND TRIM(icd_code) ~ '^J1[2-8]')
GROUP BY 1, 2, 3 ORDER BY 4 DESC;
