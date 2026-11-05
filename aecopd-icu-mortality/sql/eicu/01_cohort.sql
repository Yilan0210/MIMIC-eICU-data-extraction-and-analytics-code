-- =====================================================================
-- eICU-CRD 2.0 | 01_cohort.sql
-- Creates: aecopd.e_dx_hits, aecopd.e_cohort
-- Inclusion steps and exclusion flags (rows are flagged, not deleted)
-- Notes:
--   * AECOPD = diagnosisstring match OR icd9code match (union); each source is flagged separately
--     and reported in the flow chart (icd9code is often empty, so ICD codes alone would miss patients)
--   * One ICU stay per patient (uniquepid): the first ICU stay of the first qualifying hospitalisation
--   * Exclusion conditions are searched in diagnosis and pasthistory of all ICU stays of the same hospitalisation
-- =====================================================================
CREATE SCHEMA IF NOT EXISTS aecopd;

-- ---------------------------------------------------------------------
-- Step 0: ICU stays matching AECOPD
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.e_dx_hits;
CREATE TABLE aecopd.e_dx_hits AS
WITH d AS (
  SELECT patientunitstayid, diagnosispriority,
         (diagnosisstring ILIKE '%acute COPD exacerbation%')               AS by_string,
         (COALESCE(icd9code, '') ~ '(^|,\s*)(491\.21|J44\.1)\s*(,|$)')     AS by_icd
  FROM public.diagnosis
)
SELECT patientunitstayid,
       BOOL_OR(by_string)                                        AS hit_string,
       BOOL_OR(by_icd)                                           AS hit_icd,
       BOOL_OR(diagnosispriority = 'Primary')                    AS hit_primary
FROM d
WHERE by_string OR by_icd
GROUP BY patientunitstayid;

-- ---------------------------------------------------------------------
-- Steps 1–5 and exclusion flags
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS aecopd.e_cohort;
CREATE TABLE aecopd.e_cohort AS
WITH hosp AS (            -- hospitalisation level: counts if any ICU stay of the hospitalisation matches
  SELECT p.patienthealthsystemstayid,
         BOOL_OR(h.hit_string)  AS dx_by_string,
         BOOL_OR(h.hit_icd)     AS dx_by_icd,
         BOOL_OR(h.hit_primary) AS aecopd_primary
  FROM aecopd.e_dx_hits h
  JOIN public.patient p USING (patientunitstayid)
  GROUP BY p.patienthealthsystemstayid
),
ranked AS (               -- step 2: one stay per patient
  SELECT p.*, hosp.dx_by_string, hosp.dx_by_icd, hosp.aecopd_primary,
         ROW_NUMBER() OVER (PARTITION BY p.uniquepid
                            ORDER BY p.hospitaldischargeyear, p.patienthealthsystemstayid,
                                     p.unitvisitnumber, p.patientunitstayid) AS rn
  FROM public.patient p
  JOIN hosp USING (patienthealthsystemstayid)
),
idx AS (SELECT * FROM ranked WHERE rn = 1),
hosp_stays AS (           -- all ICU stays of the same hospitalisation, used to search the exclusion conditions
  SELECT i.patientunitstayid AS index_stay, p.patientunitstayid
  FROM idx i JOIN public.patient p USING (patienthealthsystemstayid)
),
dx_flags AS (
  SELECT s.index_stay,
    BOOL_OR(COALESCE(d.icd9code,'') ~ '(^|,\s*)20[0-8]\.'
         OR COALESCE(d.icd9code,'') ~ '(^|,\s*)C(8[1-9]|9[0-6])'
         OR d.diagnosisstring ~* '(leukemia|lymphoma|myeloma)')                      AS heme,
    BOOL_OR(COALESCE(d.icd9code,'') ~ '(^|,\s*)(042|V08|B20|Z21)'
         OR d.diagnosisstring ~* '\m(hiv|aids)\M')                                   AS hiv,
    BOOL_OR(d.diagnosisstring ~* 'transplant')                                       AS tx,
    BOOL_OR(COALESCE(d.icd9code,'') ~ '(^|,\s*)(571\.[256]|K70\.3|K74\.[3-6])'
         OR d.diagnosisstring ~* 'cirrhosis')                                        AS cirr
  FROM hosp_stays s JOIN public.diagnosis d ON d.patientunitstayid = s.patientunitstayid
  GROUP BY s.index_stay
),
ph_flags AS (
  SELECT s.index_stay,
    BOOL_OR(ph.pasthistorypath ~* '(leukemia|lymphoma|myeloma)') AS heme,
    BOOL_OR(ph.pasthistorypath ~* '\m(hiv|aids)\M')              AS hiv,
    BOOL_OR(ph.pasthistorypath ~* 'transplant')                  AS tx,
    BOOL_OR(ph.pasthistorypath ~* 'cirrhosis')                   AS cirr
  FROM hosp_stays s JOIN public.pasthistory ph ON ph.patientunitstayid = s.patientunitstayid
  GROUP BY s.index_stay
)
SELECT
  i.patientunitstayid, i.patienthealthsystemstayid, i.uniquepid, i.hospitalid,
  i.gender, i.ethnicity,
  i.age AS age_text,
  CASE WHEN TRIM(i.age) = '> 89'          THEN 90
       WHEN TRIM(i.age) ~ '^[0-9]+$'      THEN TRIM(i.age)::int END            AS age_raw,
  COALESCE(TRIM(i.age) = '> 89', false)                                         AS age_gt89,
  i.hospitaladmitoffset, i.unitdischargeoffset, i.hospitaldischargeyear,
  i.hospitaldischargestatus, i.hospitaldischargelocation,
  i.dx_by_string, i.dx_by_icd, COALESCE(i.aecopd_primary, false)                AS aecopd_primary,
  -- inclusion flags (steps 3–5)
  COALESCE(CASE WHEN TRIM(i.age) = '> 89' THEN true
                WHEN TRIM(i.age) ~ '^[0-9]+$' THEN TRIM(i.age)::int >= 18 END, false) AS incl_age,
  COALESCE(i.unitdischargeoffset >= 1440, false)                                AS incl_los24,
  COALESCE(i.hospitaldischargestatus IN ('Alive','Expired'), false)             AS incl_outcome,
  -- exclusion flags (eICU-CRD predates COVID-19, so excl_covid is always false)
  false                                                                         AS excl_covid,
  COALESCE(dx.heme, false) OR COALESCE(ph.heme, false)                          AS excl_heme,
  COALESCE(dx.hiv,  false) OR COALESCE(ph.hiv,  false)                          AS excl_hiv,
  COALESCE(dx.tx,   false) OR COALESCE(ph.tx,   false)                          AS excl_tx,
  COALESCE(dx.cirr, false) OR COALESCE(ph.cirr, false)                          AS excl_cirr
FROM idx i
LEFT JOIN dx_flags dx ON dx.index_stay = i.patientunitstayid
LEFT JOIN ph_flags ph ON ph.index_stay = i.patientunitstayid;

CREATE UNIQUE INDEX ON aecopd.e_cohort (patientunitstayid);

-- Quick check: source of the AECOPD match (for reporting)
SELECT dx_by_string, dx_by_icd, COUNT(*) AS n FROM aecopd.e_cohort GROUP BY 1, 2 ORDER BY 1, 2;
