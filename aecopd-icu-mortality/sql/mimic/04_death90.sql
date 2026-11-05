-- =====================================================================
-- MIMIC-IV 3.1 | 04_death90.sql: sensitivity analysis A10, all-cause death within 90 days (MIMIC-IV only)
-- Requires aecopd.m_cohort from 01_cohort.sql
-- Definition: death within 90 days of the day of ICU admission (intime) = 1, otherwise 0
--   Date of death from mimiciv_hosp.patients.dod: in-hospital deaths and out-of-hospital deaths linked from state death records
--   (out-of-hospital deaths more than 1 year after discharge are not captured; the 90-day window lies within this range, so there are no missing values)
-- Output: stay_id, death_90d; main population only (pop = 'main'); matched to the analysis population by stay_id
-- =====================================================================
SELECT c.stay_id,
       CASE WHEN p.dod IS NOT NULL AND p.dod <= (c.intime::date + 90) THEN 1 ELSE 0 END AS death_90d
FROM aecopd.m_cohort c
JOIN mimiciv_hosp.patients p ON p.subject_id = c.subject_id
WHERE c.pop = 'main'
ORDER BY c.stay_id
