-- =====================================================================
-- eICU-CRD 2.0 | 00_check_names.sql
-- Purpose: before extraction, check how diagnosis strings, laboratory names, nurse charting labels and drug names are actually written
-- Run each query in the eICU database; the results contain only strings and counts (no patient-level data)
-- Save each result as CSV (for example in docs/name_checks/) for the supplementary material
-- =====================================================================

-- C1. AECOPD-related diagnosis strings and ICD codes (check that no spelling of the inclusion definition is missed)
SELECT diagnosisstring, icd9code, COUNT(DISTINCT patientunitstayid) AS n_stays
FROM public.diagnosis
WHERE (diagnosisstring ILIKE '%COPD%' AND diagnosisstring ILIKE '%exacerbation%')
   OR icd9code ~ '491\.2|J44\.'
GROUP BY 1, 2 ORDER BY 3 DESC;

-- C2. Diagnosis strings matched by the exclusion conditions (check each one manually and remove false matches)
SELECT CASE
         WHEN icd9code ~ '(^|,\s*)20[0-8]\.' OR icd9code ~ '(^|,\s*)C(8[1-9]|9[0-6])'
              OR diagnosisstring ~* '(leukemia|lymphoma|myeloma)' THEN 'heme'
         WHEN icd9code ~ '(^|,\s*)(042|V08|B20|Z21)' OR diagnosisstring ~* '\m(hiv|aids)\M' THEN 'hiv'
         WHEN diagnosisstring ~* 'transplant' THEN 'tx'
         WHEN icd9code ~ '(^|,\s*)(571\.[256]|K70\.3|K74\.[3-6])' OR diagnosisstring ~* 'cirrhosis' THEN 'cirr'
       END AS category,
       diagnosisstring, icd9code, COUNT(DISTINCT patientunitstayid) AS n_stays
FROM public.diagnosis
WHERE icd9code ~ '(^|,\s*)20[0-8]\.' OR icd9code ~ '(^|,\s*)C(8[1-9]|9[0-6])'
   OR diagnosisstring ~* '(leukemia|lymphoma|myeloma|transplant|cirrhosis)'
   OR icd9code ~ '(^|,\s*)(042|V08|B20|Z21)' OR diagnosisstring ~* '\m(hiv|aids)\M'
   OR icd9code ~ '(^|,\s*)(571\.[256]|K70\.3|K74\.[3-6])'
GROUP BY 1, 2, 3 ORDER BY 1, 4 DESC;

-- C3. Past history paths matched by the exclusion conditions
SELECT pasthistorypath, COUNT(DISTINCT patientunitstayid) AS n_stays
FROM public.pasthistory
WHERE pasthistorypath ~* '(leukemia|lymphoma|myeloma|transplant|cirrhosis)'
   OR pasthistorypath ~* '\m(hiv|aids)\M'
GROUP BY 1 ORDER BY 2 DESC;

-- C4. Strings for atrial fibrillation and pneumonia (check that atrial fibrillation and flutter are separate)
SELECT 'dx' AS src, diagnosisstring AS text, icd9code, COUNT(DISTINCT patientunitstayid) AS n_stays
FROM public.diagnosis
WHERE diagnosisstring ~* '(atrial fib|atrial flutter|pneumonia)' OR icd9code ~ '427\.3|I48'
GROUP BY 1, 2, 3
UNION ALL
SELECT 'pasthistory', pasthistorypath, NULL, COUNT(DISTINCT patientunitstayid)
FROM public.pasthistory
WHERE pasthistorypath ~* '(atrial fib|atrial flutter)'
GROUP BY 1, 2, 3
ORDER BY 1, 4 DESC;

-- C5. Laboratory names (check that all names used in the scripts exist and whether similar spellings exist)
SELECT labname, COUNT(*) AS n_rows, MIN(labresult) AS min_val, MAX(labresult) AS max_val
FROM public.lab
WHERE labname IN ('albumin','BUN','pH','paCO2','lactate','creatinine','Hgb',
                  'platelets x 1000','WBC x 1000','-polys','-lymphs','-eos','-bands')
   OR labname ~* '(albumin|lymph|poly|neut|eos|wbc|platelet)'
GROUP BY 1 ORDER BY 1;

-- C6. Nurse charting labels (GCS, heart rate, respiration, temperature, SpO2)
SELECT nursingchartcelltypevallabel, nursingchartcelltypevalname, COUNT(*) AS n_rows
FROM public.nursecharting
WHERE nursingchartcelltypevallabel ~* '(glasgow|heart rate|respiratory rate|temperature|o2 sat)'
GROUP BY 1, 2 ORDER BY 1, 2;

-- C7. Non-numeric values of GCS Total (for example, how "unable to score" is written)
SELECT nursingchartvalue, COUNT(*) AS n_rows
FROM public.nursecharting
WHERE nursingchartcelltypevallabel = 'Glasgow coma score'
  AND nursingchartcelltypevalname = 'GCS Total'
  AND TRIM(nursingchartvalue) !~ '^[0-9]+(\.[0-9]+)?$'
GROUP BY 1 ORDER BY 2 DESC;

-- C8. Systemic corticosteroid names and routes (check the rule that excludes inhaled and topical preparations)
SELECT drugname, routeadmin, COUNT(*) AS n_rows
FROM public.medication
WHERE drugname ~* '(hydrocortisone|solu-?cortef|methylpred|solu-?medrol|\mmedrol|prednison|prednisolon|dexameth|decadron)'
GROUP BY 1, 2 ORDER BY 3 DESC;

-- C9. Discharge location and discharge status (is there a hospice category?)
SELECT hospitaldischargelocation, hospitaldischargestatus, COUNT(*) AS n_stays
FROM public.patient GROUP BY 1, 2 ORDER BY 3 DESC;

-- C10. APACHE versions and age format
SELECT apacheversion, COUNT(*) FROM public.apachepatientresult GROUP BY 1;
SELECT age, COUNT(*) FROM public.patient WHERE age !~ '^[0-9]+$' GROUP BY 1;

-- ===================== Additional variables =====================
-- C11. Treatment limitation (careplangeneral): group names and values
SELECT cplgroup, cplitemvalue, COUNT(DISTINCT patientunitstayid) AS n_stays
FROM public.careplangeneral
WHERE cplgroup ILIKE '%limit%'
GROUP BY 1, 2 ORDER BY 3 DESC;

-- C12. Drug names of vasopressor infusions
SELECT drugname, COUNT(*) AS n_rows
FROM public.infusiondrug
WHERE drugname ~* '(norepinephrine|levophed|epinephrine|adrenalin|dopamine|phenylephrine|neo-?synephrine|vasopressin|pitressin)'
GROUP BY 1 ORDER BY 2 DESC;

-- C13. APACHE admission diagnosis paths containing pneumonia (pneumonia_24h)
SELECT admitdxpath, COUNT(DISTINCT patientunitstayid) AS n_stays
FROM public.admissiondx
WHERE admitdxpath ~* 'pneumonia' OR admitdxname ~* 'pneumonia'
GROUP BY 1 ORDER BY 2 DESC;

-- C14. Laboratory names of glucose and monocytes
SELECT labname, COUNT(*) AS n_rows, MIN(labresult) AS min_val, MAX(labresult) AS max_val
FROM public.lab
WHERE labname ~* '(glucose|mono)'
GROUP BY 1 ORDER BY 1;

-- C15. Number of hospitals with data in each table (for hospitals without data, the variable is set to missing, not 0)
SELECT 'all_hospitals' AS tbl, COUNT(DISTINCT hospitalid) AS n_hospitals FROM public.patient
UNION ALL SELECT 'medication', COUNT(DISTINCT p.hospitalid)
  FROM public.patient p WHERE EXISTS (SELECT 1 FROM public.medication m WHERE m.patientunitstayid = p.patientunitstayid)
UNION ALL SELECT 'infusiondrug', COUNT(DISTINCT p.hospitalid)
  FROM public.patient p WHERE EXISTS (SELECT 1 FROM public.infusiondrug i WHERE i.patientunitstayid = p.patientunitstayid)
UNION ALL SELECT 'careplangeneral', COUNT(DISTINCT p.hospitalid)
  FROM public.patient p WHERE EXISTS (SELECT 1 FROM public.careplangeneral g WHERE g.patientunitstayid = p.patientunitstayid);

-- C16. Agreement between the two APACHE ventilation fields (the scripts use oobintubday1; oobventday1 for comparison)
SELECT oobintubday1, oobventday1, COUNT(*) AS n_stays
FROM public.apachepredvar GROUP BY 1, 2 ORDER BY 1, 2;
