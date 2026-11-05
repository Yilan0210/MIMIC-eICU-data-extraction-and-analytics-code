# =====================================================================
# 09_tables_calculator.R: tables, supplementary tables and the web calculator
#   Table 1  baseline characteristics (all eligible patients; by database, by outcome, by marker subset)
#   Table 2  core model: interquartile-range ORs, importance, tests of non-linearity
#   Table 3  model performance: apparent, bootstrap-corrected, IECV pooled (with 95% PI), cross-database; full and simplified models
#   Table 4  comparison with comparator scores: C, dC, calibration, net benefit, and sensitivity/specificity/PPV/NPV at predefined thresholds
#   Table 5  markers (Q3/Q4)
#   calculator/aecopd_risk_calculator.html: full model (with splines)
# =====================================================================
source("00_setup.R")
d02 <- readRDS(file.path(DIR_LOCAL, "02_dev.rds")); core <- d02$core; simp <- d02$simp; knots <- d02$knots
r03 <- readRDS(file.path(DIR_LOCAL, "03_iecv_results.rds")); stk <- readRDS(file.path(DIR_LOCAL, "03_iecv_stacked.rds"))
r04 <- readRDS(file.path(DIR_LOCAL, "04_crossdb_update.rds")); r05 <- readRDS(file.path(DIR_LOCAL, "05_markers.rds"))
a <- readRDS(file.path(DIR_LOCAL, "01_population.rds"))$a
SCORES <- c("BAP65", "DECAF_mod", "APACHE_IVa", "SAPSII", "SOFA")
NAMES <- c(full = "New model (full)", simplified = "New model (simplified)", BAP65 = "BAP-65", DECAF_mod = "Modified DECAF",
           APACHE_IVa = "APACHE IVa", SAPSII = "SAPS II", SOFA = "SOFA")
mm <- function(x) sub("^-", "−", x)
ci <- function(e, l, h, d = 2) { f <- function(z) mm(formatC(z, format = "f", digits = d))
  ifelse(is.na(e), "", ifelse(is.na(l), f(e), ifelse(l >= 0 & h >= 0, paste0(f(e), " (", f(l), "–", f(h), ")"), paste0(f(e), " (", f(l), " to ", f(h), ")")))) }
pv <- function(p) ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))

# ---- Table 1 ----
msg("Table 1")
a$female <- as.numeric(a$sex == "F"); a$nlr <- a$neut / pmax(a$lymph, ZERO_FLOOR); a$pni <- a$alb + 5 * a$lymph
T1 <- list(c("Age, years", "age", "cont"), c("Female sex", "female", "bin"),
  c("Race: White", "race4", "cat:White"), c("Race: Black", "race4", "cat:Black"), c("Race: other", "race4", "cat:Other"), c("Race: unknown", "race4", "cat:Unknown"),
  c("Core predictors", "", "head"), c("Blood urea nitrogen, mg/dL", "bun", "cont"), c("GCS motor score", "gcs", "int"),
  c("Heart rate, beats/min", "hr", "cont"), c("Respiratory rate, breaths/min", "rr", "cont"), c("Mean arterial pressure, mmHg", "map", "cont"),
  c("Temperature, °C", "temp", "cont1"), c("pH", "ph", "cont2"), c("PaCO₂, mmHg", "paco2", "cont"), c("Haemoglobin, g/dL", "hgb", "cont1"),
  c("Platelets, ×10⁹/L", "plt", "cont"), c("White cells, ×10⁹/L", "wbc", "cont1"), c("Glucose, mg/dL", "glucose", "cont"),
  c("Atrial fibrillation within 24 h", "af", "bin"), c("Invasive ventilation on day 1", "imv", "bin"),
  c("Immune-inflammation-nutrition markers", "", "head"), c("Albumin, g/L", "alb", "cont"), c("Neutrophils, ×10⁹/L", "neut", "cont2"),
  c("Lymphocytes, ×10⁹/L", "lymph", "cont2"), c("Eosinophils, ×10⁹/L", "eos", "cont2"), c("NLR", "nlr", "cont1"), c("PNI", "pni", "cont1"),
  c("Other variables", "", "head"), c("Lactate, mmol/L", "lactate", "cont2"), c("Creatinine, mg/dL", "creatinine", "cont2"), c("SpO₂, %", "spo2", "cont"),
  c("Vasopressor on day 1", "vaso", "bin"), c("Pneumonia within 24 h", "pneumonia24", "bin"), c("Systemic corticosteroid before sampling", "steroid", "bin"),
  c("Documented treatment limitation within 24 h", "limit_s18", "bin"),
  c("Severity scores", "", "head"), c("APACHE IVa predicted mortality, % (eICU-CRD)", "apache", "pct"), c("SAPS II (MIMIC-IV)", "sapsii", "cont"),
  c("SOFA, day 1 (MIMIC-IV)", "sofa", "int"), c("Outcome", "", "head"), c("In-hospital death", "death", "bin"))
cell <- function(x, type) {
  if (type == "head") return("")
  if (startsWith(type, "cat:")) { lv <- sub("cat:", "", type); ok <- !is.na(x); return(sprintf("%s (%.1f%%)", prettyNum(sum(x[ok] == lv), big.mark = ","), 100 * mean(x[ok] == lv))) }
  ok <- !is.na(x); if (!any(ok)) return("—")
  if (type == "bin") return(sprintf("%s (%.1f%%)", prettyNum(sum(x[ok] == 1), big.mark = ","), 100 * mean(x[ok] == 1)))
  x <- x[ok]; if (type == "pct") x <- 100 * x; q <- stats::quantile(x, c(.5, .25, .75))
  d <- switch(type, cont2 = 2, cont1 = 1, pct = 1, int = 0, cont = 0)
  sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), q[1], q[2], q[3])
}
groups <- list(Overall = rep(TRUE, nrow(a)), `MIMIC-IV` = a$db == "MIMIC", `eICU-CRD` = a$db == "eICU", Survived = a$death == 0, Died = a$death == 1,
               `Marker subset` = a$tested == 1, `Not in marker subset` = a$tested == 0)
t1 <- do.call(rbind, lapply(T1, function(r) { x <- a[[r[2]]]
  data.frame(Characteristic = if (r[3] == "head") r[1] else paste0("  ", r[1]), t(sapply(groups, function(g) if (r[3] == "head") "" else cell(x[g], r[3]))),
             `Missing, %` = if (r[3] == "head" || r[2] %in% c("apache", "sapsii", "sofa")) "" else sprintf("%.1f", 100 * mean(is.na(x))), check.names = FALSE) }))
t1 <- rbind(data.frame(Characteristic = "Patients, n", t(sapply(groups, function(g) prettyNum(sum(g), big.mark = ","))), `Missing, %` = "", check.names = FALSE), t1)
names(t1)[2:8] <- paste0(names(groups), " (n = ", sapply(groups, function(g) prettyNum(sum(g), big.mark = ",")), ")")

# ---- Table 2 ----
msg("Table 2")
orc <- d02$orc; imp <- d02$imp_tab
nlt <- do.call(rbind, lapply(SPLINE_VARS, function(v) { k <- paste0(v, "_nl"); z <- core$coef_ml[k] / sqrt(core$vcov_ml[k, k]); data.frame(var = v, p_nl = 2 * stats::pnorm(-abs(z))) }))
fmt_pt <- function(v, z) if (v %in% c("female", "af", "imv")) c("No", "Yes")[z + 1] else formatC(z, format = "f", digits = if (v == "ph") 2 else if (v %in% c("temp", "hgb")) 1 else 0)
t2 <- do.call(rbind, lapply(CORE_VARS, function(v) { o <- orc[orc$var == v, ]; im <- imp[imp$var == v, ]
  data.frame(Predictor = paste0(VAR_LAB[[v]], if (!is.na(VAR_UNIT[v]) && VAR_UNIT[v] != "") paste0(", ", VAR_UNIT[[v]]) else ""),
             Form = if (v %in% SPLINE_VARS) "Restricted cubic spline (3 knots)" else if (v %in% c("ln_plt", "ln_wbc", "ln_glu")) "Log, linear" else if (v %in% c("female", "af", "imv")) "Binary" else "Linear",
             Contrast = if (v %in% c("female", "af", "imv")) "Yes vs no" else paste(fmt_pt(v, o$to_raw), "vs", fmt_pt(v, o$from_raw)),
             `OR (95% CI)` = ci(o$OR, o$lo, o$hi), `Wald χ² (df)` = sprintf("%.1f (%d)", im$chi2, im$df), P = pv(im$p),
             `P non-linear` = if (v %in% SPLINE_VARS) pv(nlt$p_nl[nlt$var == v]) else "", check.names = FALSE) }))
t2 <- t2[order(-imp$chi2_minus_df[match(CORE_VARS, imp$var)]), ]

# ---- Table 3 ----
msg("Table 3")
P <- r03$pooled; pa <- P[P$clusters == "all eligible", ]
pp <- function(meas, mod, d = 2) { r <- pa[pa$measure == meas & pa$model == mod, ]; if (!nrow(r)) return(c("", ""))
  c(ci(r$est, r$lo, r$hi, d), if (is.na(r$pi_lo)) "" else paste0(mm(formatC(r$pi_lo, format = "f", digits = d)), " to ", mm(formatC(r$pi_hi, format = "f", digits = d)))) }
r02 <- readxl::read_excel(file.path(DIR_RESULT, "02_development.xlsx"), sheet = "performance")
cross <- r04$cross
cr <- function(lbl, meas, mod, d = 2) { r <- cross[cross$label == lbl & cross$measure == meas & cross$model == mod, ]; if (!nrow(r)) "" else ci(r$est, r$lo, r$hi, d) }
lblA <- unique(cross$label)[1]; lblB <- unique(cross$label)[2]
t3 <- do.call(rbind, lapply(c("full", "simplified"), function(m) { rr <- r02[r02$model == if (m == "full") "Full model" else "Simplified model", ]
  rbind(data.frame(Model = NAMES[m], Measure = "C statistic", Apparent = sprintf("%.3f", rr$C_apparent), `Bootstrap-corrected` = ifelse(is.na(rr$C_optimism_corrected), "", sprintf("%.3f", rr$C_optimism_corrected)),
                   `IECV pooled (95% CI)` = pp("logit_C", m, 3)[1], `IECV 95% PI` = pp("logit_C", m, 3)[2], `eICU → MIMIC` = cr(lblA, "logit_C", m, 3), `MIMIC → eICU` = cr(lblB, "logit_C", m, 3), check.names = FALSE),
        data.frame(Model = NAMES[m], Measure = "O/E ratio", Apparent = "1.00", `Bootstrap-corrected` = "", `IECV pooled (95% CI)` = pp("log_OE", m)[1], `IECV 95% PI` = pp("log_OE", m)[2],
                   `eICU → MIMIC` = cr(lblA, "log_OE", m), `MIMIC → eICU` = cr(lblB, "log_OE", m), check.names = FALSE),
        data.frame(Model = NAMES[m], Measure = "Calibration slope", Apparent = "", `Bootstrap-corrected` = ifelse(is.na(rr$calibration_slope_corrected), "", sprintf("%.2f", rr$calibration_slope_corrected)),
                   `IECV pooled (95% CI)` = pp("cal_slope", m)[1], `IECV 95% PI` = pp("cal_slope", m)[2], `eICU → MIMIC` = cr(lblA, "cal_slope", m), `MIMIC → eICU` = cr(lblB, "cal_slope", m), check.names = FALSE),
        data.frame(Model = NAMES[m], Measure = "Calibration-in-the-large", Apparent = "", `Bootstrap-corrected` = "", `IECV pooled (95% CI)` = pp("CITL", m)[1], `IECV 95% PI` = pp("CITL", m)[2],
                   `eICU → MIMIC` = cr(lblA, "CITL", m), `MIMIC → eICU` = cr(lblB, "CITL", m), check.names = FALSE),
        data.frame(Model = NAMES[m], Measure = "Scaled Brier score", Apparent = ifelse(is.na(rr$scaled_Brier_apparent), "", sprintf("%.3f", rr$scaled_Brier_apparent)),
                   `Bootstrap-corrected` = ifelse(is.na(rr$scaled_Brier_corrected), "", sprintf("%.3f", rr$scaled_Brier_corrected)), `IECV pooled (95% CI)` = "", `IECV 95% PI` = "",
                   `eICU → MIMIC` = "", `MIMIC → eICU` = "", check.names = FALSE),
        data.frame(Model = NAMES[m], Measure = sprintf("Net benefit at %s", paste0(100 * NB_THRESHOLDS, "%", collapse = " / ")), Apparent = "", `Bootstrap-corrected` = "",
                   `IECV pooled (95% CI)` = paste(sapply(NB_THRESHOLDS, function(pt) sprintf("%.3f", pa$est[pa$measure == sprintf("NB_%02d", 100 * pt) & pa$model == m])), collapse = " / "),
                   `IECV 95% PI` = "", `eICU → MIMIC` = "", `MIMIC → eICU` = "", check.names = FALSE)) }))

# ---- Table 4 ----
msg("Table 4")
t4a <- do.call(rbind, lapply(SCORES, function(s) { cs <- P[P$measure == "logit_C" & P$model == s, ]; if (!nrow(cs)) return(NULL)
  cf <- P[P$measure == "logit_C_same_patients" & P$model == paste0("full|", s), ]; dd <- P[P$measure == "dC_vs_score" & P$model == paste0("full-", s), ]
  oe <- P[P$measure == "log_OE" & P$model == s, ]; sl <- P[P$measure == "cal_slope" & P$model == s, ]
  nb <- function(mod) paste(sapply(NB_THRESHOLDS, function(pt) { r <- P[P$measure == sprintf("NB_%02d", 100 * pt) & P$model == mod, ]; if (nrow(r)) sprintf("%.3f", r$est) else "" }), collapse = " / ")
  data.frame(Score = NAMES[s], Clusters = sub("clusters ", "", cs$clusters), `C of score (95% CI)` = ci(cs$est, cs$lo, cs$hi, 3), `C of new model, same patients` = ci(cf$est, cf$lo, cf$hi, 3),
             `ΔC, new model − score (95% CI)` = ci(dd$est, dd$lo, dd$hi, 3), `O/E of score` = ci(oe$est, oe$lo, oe$hi), `Calibration slope of score` = ci(sl$est, sl$lo, sl$hi),
             `Net benefit 10/20/30%: score` = nb(s), `Net benefit 10/20/30%: new model` = nb(paste0("full|", s)),
             Note = if (cs$k < 3) "Two clusters: fixed-effect pooling" else if (s == "DECAF_mod") "Lacks eMRCD dyspnoea grade" else "", check.names = FALSE) }))
cls_one <- function(y, p, pt) { h <- p >= pt; tp <- sum(h & y == 1); fp <- sum(h & y == 0); fn <- sum(!h & y == 1); tn <- sum(!h & y == 0)
  lg <- function(a, b) { a <- a + 0.5; b <- b + 0.5; c(qlogis(a / (a + b)), 1 / a + 1 / b) }
  rbind(sens = lg(tp, fn), spec = lg(tn, fp), ppv = lg(tp, fp), npv = lg(tn, fn), flagged = lg(tp + fp, tn + fn)) }
SETS <- list(list(name = "All patients", rows = function(s) rep(TRUE, nrow(s)), preds = c("full", "simplified", "BAP65", "DECAF_mod")),
             list(name = "eICU-CRD with APACHE IVa", rows = function(s) !is.na(s$APACHE_IVa), preds = c("full", "APACHE_IVa")),
             list(name = "MIMIC-IV", rows = function(s) !is.na(s$SAPSII) & !is.na(s$SOFA), preds = c("full", "SAPSII", "SOFA")))
t4b <- do.call(rbind, lapply(SETS, function(st) do.call(rbind, lapply(st$preds, function(pr) do.call(rbind, lapply(NB_THRESHOLDS, function(pt) {
  per <- lapply(stk, function(s) { k <- st$rows(s) & !is.na(s[[pr]]); cls_one(s$y[k], s[[pr]][k], pt) })
  out <- sapply(rownames(per[[1]]), function(m) { r <- rubin(sapply(per, function(z) z[m, 1]), sapply(per, function(z) z[m, 2]))
    sprintf("%.1f (%.1f–%.1f)", 100 * plogis(r$est), 100 * plogis(r$lo), 100 * plogis(r$hi)) })
  k1 <- st$rows(stk[[1]]) & !is.na(stk[[1]][[pr]])
  data.frame(Population = st$name, Model = NAMES[pr], Threshold = paste0(100 * pt, "%"), n = sum(k1), Deaths = sum(stk[[1]]$y[k1]), `High risk, %` = out["flagged"],
             `Sensitivity, %` = out["sens"], `Specificity, %` = out["spec"], `PPV, %` = out["ppv"], `NPV, %` = out["npv"], check.names = FALSE, row.names = NULL) }))))))

# ---- Table 5 ----
msg("Table 5")
dl <- r05$delta; ip <- r05$iecv_pool
lbl5 <- c("E1 vs base" = "+ ln neutrophils, ln lymphocytes, albumin (primary)", "E2 vs base" = "+ eosinopenia", "E3 vs base" = "+ splines of the three markers")
t5a <- do.call(rbind, lapply(dl$comparison, function(cmp) { d <- dl[dl$comparison == cmp, ]; i <- ip[ip$comparison == cmp, ]
  data.frame(`Added to core model LP` = lbl5[[cmp]], df = d$df, `LR test P` = pv(d$LR_p), `New information fraction` = sprintf("%.3f", d$new_info_fraction),
             `ΔC bootstrap-corrected (95% CI)` = ci(d$dC_corrected, d$dC_corr_lo, d$dC_corr_hi, 3), `ΔC IECV pooled (95% CI)` = ci(i$est, i$lo, i$hi, 3),
             `Δ scaled Brier` = sprintf("%.3f", d$d_scaled_Brier_corrected), check.names = FALSE) }))
t5b <- r05$h1b; t5c <- r05$h2; t5d <- r05$s16

# ---- Calculator (full model; splines computed in the web page with the same formula) ----
msg("Calculator")
IN <- list(age = c("Age (years)", 18, 110), bun = c("Blood urea nitrogen (mg/dL)", 1, 300), hr = c("Heart rate (beats/min), first on ICU day 1", 20, 250),
           rr = c("Respiratory rate (breaths/min), first on ICU day 1", 2, 80), map = c("Mean arterial pressure (mmHg), first on ICU day 1", 20, 250),
           temp = c("Temperature (°C), first on ICU day 1", 25, 45), ph = c("pH (first blood gas)", 6.8, 7.8), paco2 = c("PaCO₂ (mmHg), same blood gas", 10, 200),
           gcs = c("GCS motor score (1–6), first on ICU day 1", 1, 6), hgb = c("Haemoglobin (g/dL)", 2, 25), plt = c("Platelets (×10⁹/L)", 1, 2000),
           wbc = c("White cells (×10⁹/L)", 0.1, 500), glucose = c("Glucose (mg/dL)", 10, 1500))
CHK <- c(female = "Female sex", af = "Atrial fibrillation within 24 h", imv = "Invasive ventilation on ICU day 1 (not via tracheostomy)")
JSX <- c(age = "Math.min(v('age'),90)", ln_bun = "Math.log(Math.max(v('bun'),1))", hr = "v('hr')", rr = "v('rr')", map = "v('map')", temp = "v('temp')", ph = "v('ph')",
         paco2 = "v('paco2')", female = "c('female')", gcs = "v('gcs')", hgb = "v('hgb')", ln_plt = "Math.log(Math.max(v('plt'),1))", ln_wbc = "Math.log(Math.max(v('wbc'),0.1))",
         ln_glu = "Math.log(Math.max(v('glucose'),10))", af = "c('af')", imv = "c('imv')")
cf <- core$coef
terms_js <- paste(vapply(CORE_VARS, function(v) { s <- sprintf("(%.8f)*X.%s", cf[[v]], v)
  if (v %in% SPLINE_VARS) s <- paste0(s, sprintf(" + (%.8f)*rcs(X.%s,[%s])", cf[[paste0(v, "_nl")]], v, paste(sprintf("%.8f", knots[[v]]), collapse = ","))); s }, character(1)), collapse = " + ")
fields <- paste(c(sprintf('<label>%s<input type="number" id="%s" min="%s" max="%s" step="any" required></label>', sapply(IN, `[`, 1), names(IN), sapply(IN, `[`, 2), sapply(IN, `[`, 3)),
                  sprintf('<label class="chk"><input type="checkbox" id="%s"> %s</label>', names(CHK), CHK)), collapse = "\n")
iecv_c <- pa[pa$measure == "logit_C" & pa$model == "full", ]
html <- paste0('<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>AECOPD ICU Mortality Risk</title>
<style>
:root{--bg:#fff;--fg:#1a1a1a;--muted:#5f6368;--line:#d0d4d9;--accent:#1f4e79;--card:#f6f8fa}
@media (prefers-color-scheme:dark){:root{--bg:#15181c;--fg:#e8eaed;--muted:#9aa0a6;--line:#3a3f45;--accent:#8ab4f8;--card:#1f2328}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,"Segoe UI",Arial,sans-serif}
main{max-width:680px;margin:0 auto;padding:24px 16px}h1{font-size:22px;margin:0 0 4px}p.sub{color:var(--muted);margin:0 0 20px}
form{display:grid;grid-template-columns:1fr 1fr;gap:12px 16px}label{display:flex;flex-direction:column;font-size:13px;color:var(--muted);gap:4px}
label.chk{flex-direction:row;align-items:center;gap:8px;color:var(--fg);grid-column:span 2}
input[type=number]{font:inherit;padding:8px;border:1px solid var(--line);border-radius:6px;background:var(--bg);color:var(--fg)}
button{grid-column:span 2;font:inherit;padding:10px;border:0;border-radius:6px;background:var(--accent);color:var(--bg);cursor:pointer}
#out{margin-top:20px;padding:16px;border-radius:8px;background:var(--card);display:none}#risk{font-size:34px;font-weight:600}
small{display:block;color:var(--muted);margin-top:20px;font-size:12px}@media (max-width:520px){form{grid-template-columns:1fr}label.chk,button{grid-column:span 1}}
</style></head><body><main>
<h1>In-hospital mortality after ICU admission for AECOPD</h1>
<p class="sub">Full model, 16 predictors available within 24 h of ICU admission. Internal-external cross-validation C statistic ',
  sprintf("%.2f (95%% CI %.2f–%.2f)", iecv_c$est, iecv_c$lo, iecv_c$hi), '.</p>
<form id="f">
', fields, '
<button type="submit">Calculate risk</button></form>
<div id="out"><div>Predicted risk of in-hospital death</div><div id="risk"></div><div id="cat"></div></div>
<small>For research use only; not a substitute for clinical judgement. Developed in US ICU databases (MIMIC-IV, eICU-CRD); applies to adults with
AECOPD who remain in the ICU for at least 24 h. Calibration varied between hospitals: update the intercept with local data before use in a new setting.</small>
</main><script>
const v=id=>parseFloat(document.getElementById(id).value), c=id=>document.getElementById(id).checked?1:0;
const p3=z=>Math.pow(Math.max(z,0),3);
const rcs=(x,k)=>(p3(x-k[0])-p3(x-k[1])*(k[2]-k[0])/(k[2]-k[1])+p3(x-k[2])*(k[1]-k[0])/(k[2]-k[1]))/Math.pow(k[2]-k[0],2);
document.getElementById("f").addEventListener("submit",e=>{e.preventDefault();
 const X={', paste(sprintf("%s:%s", names(JSX), JSX), collapse = ","), '};
 const lp=(', sprintf("%.8f", cf[["(Intercept)"]]), ') + ', terms_js, ';
 const p=1/(1+Math.exp(-lp)); document.getElementById("risk").textContent=(100*p).toFixed(1)+"%";
 document.getElementById("cat").textContent=p<0.10?"Risk group: <10%":(p<0.30?"Risk group: 10–30%":"Risk group: ≥30%");
 document.getElementById("out").style.display="block";});
</script></body></html>')
dir.create(file.path(OUT_DIR, "calculator"), showWarnings = FALSE)
writeLines(html, file.path(OUT_DIR, "calculator", "aecopd_risk_calculator.html"), useBytes = TRUE)
ex <- d02$xs[[1]][1, ]
writeLines(c(sprintf("Calculator check: predicted risk of patient 1 in imputed dataset 1 in R = %.4f", predict_dev(core, ex)),
             paste(sprintf("%s = %s", c(names(IN), names(CHK)), c(sapply(names(IN), function(k) round(ex[[k]], 4)), sapply(names(CHK), function(k) ex[[k]]))), collapse = "; ")),
           file.path(DIR_LOCAL, "09_calculator_check_example.txt"), useBytes = TRUE)   # contains individual values; local only

write_share(list(
  README = data.frame(note = c("Table1: observed values before imputation; continuous variables as median (IQR), categorical variables as n (%, denominator = non-missing)", "Table2: ORs for interquartile-range contrasts (pooled ML coefficients); full equation in 02_development.xlsx",
    "Table3: apparent, bootstrap-corrected, IECV pooled (95% CI and 95% PI) and cross-database performance of the full and simplified models", "Table4a: comparison with comparator scores; Table4b: classification at the predefined thresholds (Rubin's rules, logit scale)",
    "Table5: added value of the markers (Q3), H1b, ratio constraints and index comparison (Q4)", "This file contains summary statistics only")),
  Table1 = t1, Table2 = t2, Table3 = t3, Table4a_scores = t4a, Table4b_thresholds = t4b, Table5a_markers = t5a, Table5b_H1b = t5b,
  Table5c_ratio_constraints = t5c, Table5d_indices = t5d, eTable_model_updating = r04$upd), "09_tables.xlsx")
msg("09 done")
