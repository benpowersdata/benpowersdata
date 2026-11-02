* ══════════════════════════════════════════════════════════════════════════
* Sensitivity Analysis: IPW-Weighted Survival Analysis, 2011–2015 Cohort
* Manuscript: "The impact of neighborhood socioeconomic deprivation on
*              metastatic pancreatic cancer treatment and survival"
* Authors: Raghu, Shah, Pattnaik et al.
* ──────────────────────────────────────────────────────────────────────────
* Variable dictionary:
*   adinatquart   ADI national quartile (1=least deprived, 4=most deprived)
*   futimemonths  Follow-up time in months (continuous)
*   vitalstatus   0=alive, 1=dead
*   id            Unique patient identifier
*   age           Age at diagnosis (continuous, years)
*   race          Race (White/Black/Asian/Other/Unknown)
*   ethnicity     Hispanic ethnicity (Non-Hispanic/Hispanic/Unknown)
*   sex           Sex (Male/Female)
*   chemobin      Systemic therapy receipt (none/chemo/missing)
*   datediagnosis Date of diagnosis (Stata numeric date)
*   year_dx       Year of diagnosis (derived: gen year_dx = year(datediagnosis))
* ──────────────────────────────────────────────────────────────────────────
* NOTE: FCDS data are not publicly available due to a data use agreement
*       with the Florida Department of Health. Researchers may apply for
*       access at https://fcds.med.miami.edu
* ══════════════════════════════════════════════════════════════════════════

* ── Sensitivity Analysis: Restrict to 2011-2015 ───────────────────────────

preserve
keep if inrange(year_dx, 2011, 2015)

* ── Sample size and chemo rate ────────────────────────────────────────────

count
tab chemobin, missing
tab adinatquart chemobin, row missing

* ── Overall chemo rate among patients with known treatment status ──────────

tab chemobin if !missing(chemobin)

* ── Logistic regression: SES → chemo receipt ──────────────────────────────

logit chemobin i.adinatquart c.age i.race i.ethnicity i.sex ///
    if !missing(chemobin), vce(robust) or

margins i.adinatquart, post

* ── Re-estimate propensity score and stabilized IPW weights ───────────────
* Confounders (DAG-derived minimally sufficient adjustment set):
*   age, race, ethnicity, sex
* Insurance and chemotherapy are mediators — excluded from PS model

cap drop p1_s p2_s p3_s p4_s sw_s sw_s_trim sw_round_s

mlogit adinatquart c.age i.race i.ethnicity i.sex, ///
    baseoutcome(1) vce(robust)

predict p1_s p2_s p3_s p4_s, pr

quietly tab adinatquart, matcell(qcounts)
scalar total_s = _N

gen sw_s = .
replace sw_s = (qcounts[1,1]/total_s) / p1_s if adinatquart == 1
replace sw_s = (qcounts[2,1]/total_s) / p2_s if adinatquart == 2
replace sw_s = (qcounts[3,1]/total_s) / p3_s if adinatquart == 3
replace sw_s = (qcounts[4,1]/total_s) / p4_s if adinatquart == 4

* Trim stabilized weights at 1st/99th percentiles
winsor2 sw_s, cuts(1 99) suffix(_trim)

summarize sw_s_trim, detail

* ── IPW-weighted Cox model ─────────────────────────────────────────────────
* pweight not supported by stset; use rounded frequency weights (fw)
* for model fitting; pweight used only for point estimate check

stset futimemonths [pw=sw_s_trim], failure(vitalstatus==1) id(id)

stcox i.adinatquart, vce(robust) nolog

* ── Bootstrap CIs (500 replications, seed 42 for reproducibility) ──────────

capture program drop ipw_boot_sensitivity
program define ipw_boot_sensitivity, rclass
    preserve
    bsample

    quietly mlogit adinatquart c.age i.race i.ethnicity i.sex, ///
        baseoutcome(1)
    quietly predict bp1 bp2 bp3 bp4, pr

    quietly tab adinatquart, matcell(qc)
    scalar bt = _N
    gen double bsw = .
    replace bsw = (qc[1,1]/bt) / bp1 if adinatquart == 1
    replace bsw = (qc[2,1]/bt) / bp2 if adinatquart == 2
    replace bsw = (qc[3,1]/bt) / bp3 if adinatquart == 3
    replace bsw = (qc[4,1]/bt) / bp4 if adinatquart == 4

    quietly centile bsw, centile(1 99)
    quietly replace bsw = r(c_1) if bsw < r(c_1)
    quietly replace bsw = r(c_2) if bsw > r(c_2)

    * Round to integer frequency weights (pweight not supported in stset)
    quietly gen bsw_r = round(bsw)
    quietly replace bsw_r = 1 if bsw_r == 0

    quietly stset futimemonths [fw=bsw_r], failure(vitalstatus==1) id(id)
    quietly stcox i.adinatquart, nolog

    return scalar hr_q2 = exp(_b[2.adinatquart])
    return scalar hr_q3 = exp(_b[3.adinatquart])
    return scalar hr_q4 = exp(_b[4.adinatquart])

    restore
end

bootstrap hr_q2=r(hr_q2) hr_q3=r(hr_q3) hr_q4=r(hr_q4), ///
    reps(500) seed(42) nodots: ipw_boot_sensitivity

estat bootstrap, percentile

* ── Weighted KM curves ────────────────────────────────────────────────────

gen sw_round_s = round(sw_s_trim)
replace sw_round_s = 1 if sw_round_s == 0

stset futimemonths [fw=sw_round_s], failure(vitalstatus==1) id(id)

sts graph, by(adinatquart) ///
    title("Sensitivity Analysis: 2011–2015 Cohort") ///
    subtitle("IPW-weighted KM curves (post-FOLFIRINOX era)") ///
    xtitle("Months from diagnosis") ///
    ytitle("Overall survival probability") ///
    legend(order(1 "Q1 — Least deprived (ADI 1–35)" ///
                 2 "Q2 (ADI 36–57)" ///
                 3 "Q3 (ADI 58–79)" ///
                 4 "Q4 — Most deprived (ADI 80–100)")) ///
    risktable

* ── Median OS by quartile ─────────────────────────────────────────────────

stset futimemonths, failure(vitalstatus==1) id(id)
stci, by(adinatquart)

restore
