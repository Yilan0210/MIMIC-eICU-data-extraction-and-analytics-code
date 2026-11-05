-- =====================================================================
-- MIMIC-IV 3.1 | 00_check_names.sql
-- Purpose: before extraction, check the codes, the charted values of heart rhythm, the names/routes of systemic corticosteroids and the emar event labels
-- The results contain only strings and counts (no patient-level data)
-- =====================================================================

-- C1. Distribution of the AECOPD codes in the main and related populations (by ICD version)
SELECT icd_version, TRIM(icd_code) AS icd_code, COUNT(DISTINCT hadm_id) AS n_hadm
FROM mimiciv_hosp.diagnoses_icd
WHERE (icd_version = 9  AND TRIM(icd_code) IN ('49121','49122','49120','4912'))
   OR (icd_version = 10 AND TRIM(icd_code) ~ '^J44')
GROUP BY 1, 2 ORDER BY 1, 2;

-- C2. Codes matched by the exclusion conditions (check and archive)
SELECT icd_version, TRIM(icd_code) AS icd_code, COUNT(DISTINCT hadm_id) AS n_hadm
FROM mimiciv_hosp.diagnoses_icd
WHERE (icd_version = 9  AND (TRIM(icd_code) ~ '^20[0-8]' OR TRIM(icd_code) = '042' OR TRIM(icd_code) ~ '^V08'
                            OR TRIM(icd_code) IN ('V420','V421','V426','V427','V4281','V4282','V4283','V4284',
                                                  '5712','5715','5716')))
   OR (icd_version = 10 AND (TRIM(icd_code) ~ '^C(8[1-9]|9[0-6])' OR TRIM(icd_code) ~ '^B20' OR TRIM(icd_code) = 'Z21'
                            OR TRIM(icd_code) ~ '^Z94[0-4]' OR TRIM(icd_code) ~ '^Z948[1-4]'
                            OR TRIM(icd_code) ~ '^K(703|74[3-6])' OR TRIM(icd_code) IN ('U071','J1282')))
GROUP BY 1, 2 ORDER BY 1, 2;

-- C3. Charted values of the heart rhythm item (itemid 220048): how atrial fibrillation is written
SELECT value, COUNT(*) AS n_rows
FROM mimiciv_icu.chartevents
WHERE itemid = 220048
GROUP BY 1 ORDER BY 2 DESC;

-- C4. Systemic corticosteroid names and event types in emar
SELECT medication, event_txt, COUNT(*) AS n_rows
FROM mimiciv_hosp.emar
WHERE medication ~* '(hydrocortisone|solu-?cortef|methylpred|solu-?medrol|\mmedrol|prednison|prednisolon|dexameth|decadron)'
GROUP BY 1, 2 ORDER BY 3 DESC;

-- C5. Systemic corticosteroid names and routes in prescriptions
SELECT drug, route, COUNT(*) AS n_rows
FROM mimiciv_hosp.prescriptions
WHERE drug ~* '(hydrocortisone|solu-?cortef|methylpred|solu-?medrol|\mmedrol|prednison|prednisolon|dexameth|decadron)'
GROUP BY 1, 2 ORDER BY 3 DESC;

-- C6. Years covered by emar (prescriptions are used when emar is incomplete)
SELECT p.anchor_year_group,
       COUNT(DISTINCT a.hadm_id)                                           AS n_hadm,
       COUNT(DISTINCT e.hadm_id)                                           AS n_hadm_with_emar
FROM mimiciv_hosp.admissions a
JOIN mimiciv_hosp.patients p USING (subject_id)
LEFT JOIN (SELECT DISTINCT hadm_id FROM mimiciv_hosp.emar) e USING (hadm_id)
GROUP BY 1 ORDER BY 1;

-- C7. Blood gas specimen types and discharge locations
SELECT specimen, COUNT(*) FROM mimiciv_derived.bg GROUP BY 1 ORDER BY 2 DESC;
SELECT discharge_location, COUNT(*) FROM mimiciv_hosp.admissions GROUP BY 1 ORDER BY 2 DESC;

-- C8. Ventilation status values: check the spelling of 'InvasiveVent' and the number of 'Tracheostomy' stays
--     Ventilated patients with a tracheostomy may be recorded as 'Tracheostomy' in the derived table
SELECT ventilation_status, COUNT(DISTINCT stay_id) AS n_stays
FROM mimiciv_derived.ventilation GROUP BY 1 ORDER BY 2 DESC;

-- C9. itemids and values of code status (treatment limitation)
SELECT di.itemid, di.label, ce.value, COUNT(*) AS n_rows
FROM mimiciv_icu.d_items di
JOIN mimiciv_icu.chartevents ce ON ce.itemid = di.itemid
WHERE di.label ILIKE '%code status%'
GROUP BY 1, 2, 3 ORDER BY 4 DESC;
-- Pneumonia codes in emergency department diagnoses are checked at the end of 03_ed_steroid.sql (skipped if MIMIC-IV-ED is not installed)
