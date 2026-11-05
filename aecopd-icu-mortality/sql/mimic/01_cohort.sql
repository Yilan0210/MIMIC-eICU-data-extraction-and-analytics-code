-- =====================================================================
-- MIMIC-IV 3.1 | 01_cohort.sql
-- Creates: aecopd.m_hadm_dx, aecopd.m_cohort
-- Inclusion flags, exclusion flags, and a related population for supplementary checks
-- Notes:
--   * Main population (pop = 'main'): hospital admission with ICD-9 49121 or ICD-10 J441
--   * Related population (pop = 'related'): patients never in the main population, with ICD-9 49122 or ICD-10 J440
--   * One ICU stay per patient within each population: the earliest ICU stay (intime) during a qualifying admission
--   * Independent of any earlier cohort table; the counts can be compared with an earlier extraction (690 patients, 120 deaths)
-- =====================================================================
CREATE SCHEMA IF NOT EXISTS aecopd;

DROP TABLE IF EXISTS aecopd.m_hadm_dx;
CREATE TABLE aecopd.m_hadm_dx AS
WITH d AS (
  SELECT hadm_id, subject_id, seq_num, icd_version, TRIM(icd_code) AS code
  FROM mimiciv_hosp.diagnoses_icd
)
SELECT hadm_id, subject_id,
  BOOL_OR((icd_version = 9 AND code = '49121') OR (icd_version = 10 AND code = 'J441'))                AS dx_aecopd,
  BOOL_OR((icd_version = 9 AND code = '49122') OR (icd_version = 10 AND code = 'J440'))                AS dx_related,
  BOOL_OR(seq_num = 1 AND ((icd_version = 9 AND code = '49121') OR (icd_version = 10 AND code = 'J441'))) AS aecopd_primary,
  MAX(icd_version)                                                                                    AS icd_version
FROM d
GROUP BY hadm_id, subject_id
HAVING BOOL_OR((icd_version = 9 AND code IN ('49121','49122')) OR (icd_version = 10 AND code IN ('J441','J440')));

DROP TABLE IF EXISTS aecopd.m_cohort;
CREATE TABLE aecopd.m_cohort AS
WITH stays AS (
  SELECT i.subject_id, i.hadm_id, i.stay_id, i.intime, i.outtime, i.los,
         h.dx_aecopd, h.dx_related, h.aecopd_primary, h.icd_version
  FROM mimiciv_icu.icustays i
  JOIN aecopd.m_hadm_dx h ON h.hadm_id = i.hadm_id
),
main_subj AS (SELECT DISTINCT subject_id FROM stays WHERE dx_aecopd),
cand AS (
  SELECT s.*, CASE WHEN s.dx_aecopd THEN 'main' ELSE 'related' END AS pop
  FROM stays s
  WHERE s.dx_aecopd
     OR (s.dx_related AND NOT EXISTS (SELECT 1 FROM main_subj m WHERE m.subject_id = s.subject_id))
),
ranked AS (
  SELECT c.*, ROW_NUMBER() OVER (PARTITION BY c.subject_id ORDER BY c.intime, c.stay_id) AS rn
  FROM cand c
),
idx AS (
  SELECT r.*, a.admittime, a.edregtime, a.discharge_location, a.race, a.hospital_expire_flag,
         LEAST(a.edregtime, a.admittime)                                         AS arrival,
         p.gender, p.anchor_year_group,
         p.anchor_age + (EXTRACT(YEAR FROM a.admittime)::int - p.anchor_year)    AS age_raw,
         (p.anchor_age > 89)                                                     AS age_gt89
  FROM ranked r
  JOIN mimiciv_hosp.admissions a ON a.hadm_id = r.hadm_id
  JOIN mimiciv_hosp.patients   p ON p.subject_id = r.subject_id
  WHERE r.rn = 1
),
excl AS (
  SELECT i.stay_id,
    BOOL_OR(d.icd_version = 10 AND TRIM(d.icd_code) IN ('U071','J1282'))                                   AS covid,
    BOOL_OR((d.icd_version = 9  AND TRIM(d.icd_code) ~ '^20[0-8]')
         OR (d.icd_version = 10 AND TRIM(d.icd_code) ~ '^C(8[1-9]|9[0-6])'))                               AS heme,
    BOOL_OR((d.icd_version = 9  AND (TRIM(d.icd_code) = '042' OR TRIM(d.icd_code) ~ '^V08'))
         OR (d.icd_version = 10 AND (TRIM(d.icd_code) ~ '^B20' OR TRIM(d.icd_code) = 'Z21')))              AS hiv,
    BOOL_OR((d.icd_version = 9  AND TRIM(d.icd_code) IN ('V420','V421','V426','V427','V4281','V4282','V4283','V4284'))
         OR (d.icd_version = 10 AND (TRIM(d.icd_code) ~ '^Z94[0-4]' OR TRIM(d.icd_code) ~ '^Z948[1-4]')))  AS tx,
    BOOL_OR((d.icd_version = 9  AND TRIM(d.icd_code) IN ('5712','5715','5716'))
         OR (d.icd_version = 10 AND TRIM(d.icd_code) ~ '^K(703|74[3-6])'))                                 AS cirr
  FROM idx i
  JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id = i.hadm_id
  GROUP BY i.stay_id
)
SELECT i.subject_id, i.hadm_id, i.stay_id, i.pop, i.intime, i.outtime, i.los,
       i.admittime, i.edregtime, i.arrival,
       i.icd_version, COALESCE(i.aecopd_primary, false) AS aecopd_primary,
       i.age_raw, i.age_gt89, i.gender, i.race, i.anchor_year_group,
       i.discharge_location, i.hospital_expire_flag,
       -- inclusion flags
       COALESCE(i.age_raw >= 18, false)                  AS incl_age,
       COALESCE(i.los >= 1, false)                       AS incl_los24,
       (i.hospital_expire_flag IS NOT NULL)              AS incl_outcome,
       -- exclusion flags
       COALESCE(e.covid, false) AS excl_covid, COALESCE(e.heme, false) AS excl_heme,
       COALESCE(e.hiv,   false) AS excl_hiv,   COALESCE(e.tx,   false) AS excl_tx,
       COALESCE(e.cirr,  false) AS excl_cirr
FROM idx i
LEFT JOIN excl e USING (stay_id);

CREATE UNIQUE INDEX ON aecopd.m_cohort (stay_id);
CREATE INDEX ON aecopd.m_cohort (subject_id);

SELECT pop, COUNT(*) AS n, SUM(hospital_expire_flag) AS deaths FROM aecopd.m_cohort GROUP BY pop;
